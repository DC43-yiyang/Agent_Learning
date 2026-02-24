"""
SkillLoader + GeneAgent — Orchestration Layer
=============================================
Responsibility: load skills, route queries, run the ReAct loop.

This module has NO knowledge of HTTP calls or external APIs.
It only knows: "ask tools.execute(), add result to messages, repeat."

Dependencies:
  tools.py  → GeneTools (injected via __init__)
  skills/   → skill markdown files (loaded by SkillLoader)
"""

import json
import os

import ollama

from tools import GeneTools


# ─────────────────────────────────────────────
# SkillLoader
# Reads .md files from the skills directory and parses frontmatter.
# ─────────────────────────────────────────────

class SkillLoader:
    def __init__(self, skills_dir: str = "./skills"):
        self.skills_dir = skills_dir
        self._cache = {}

    def load(self, skill_name: str) -> dict:
        if skill_name in self._cache:
            return self._cache[skill_name]

        skill_path = os.path.join(self.skills_dir, f"{skill_name}.md")
        if not os.path.exists(skill_path):
            raise FileNotFoundError(f"Skill file not found: {skill_path}")

        with open(skill_path, "r", encoding="utf-8") as f:
            raw = f.read()

        skill = self._parse(skill_name, raw)
        self._cache[skill_name] = skill
        return skill

    def _parse(self, name: str, raw: str) -> dict:
        metadata = {"name": name, "description": "", "version": "1.0", "tools": []}
        body = raw

        if raw.startswith("---"):
            parts = raw.split("---", 2)
            if len(parts) >= 3:
                frontmatter = parts[1].strip()
                body = parts[2].strip()
                for line in frontmatter.splitlines():
                    if ":" in line:
                        key, _, val = line.partition(":")
                        key = key.strip()
                        val = val.strip()
                        if key == "description":
                            metadata["description"] = val
                        elif key == "version":
                            metadata["version"] = val
                    elif line.strip().startswith("-"):
                        tool = line.strip().lstrip("-").strip()
                        if tool:
                            metadata["tools"].append(tool)

        return {"metadata": metadata, "body": body, "system_prompt": body}

    def list_skills(self) -> list:
        if not os.path.exists(self.skills_dir):
            return []
        return [f.replace(".md", "") for f in os.listdir(self.skills_dir) if f.endswith(".md")]


# ─────────────────────────────────────────────
# GeneAgent
# Orchestrates skill routing and the ReAct loop.
# Delegates all tool execution to a GeneTools instance.
# ─────────────────────────────────────────────

class GeneAgent:
    def __init__(self, tools: GeneTools, skills_dir: str = "./skills"):
        self.tools = tools
        self.loader = SkillLoader(skills_dir)
        self.skills = {
            name: self.loader.load(name)
            for name in self.loader.list_skills()
        }
        print(f"✅ Loaded {len(self.skills)} skills: {list(self.skills.keys())}")

    # ── Routing ───────────────────────────────────────────────────────────────

    def _route(self, query: str, verbose: bool) -> list[str]:
        """Use the LLM to select which skill(s) are needed for this query.

        The skill's frontmatter `description` is what the router reads —
        that is why writing a good description in each .md file matters.
        """
        menu = "\n".join(
            f"- {name}: {skill['metadata']['description']}"
            for name, skill in self.skills.items()
        )
        response = ollama.chat(
            model="qwen3:30b-a3b",
            messages=[{
                "role": "user",
                "content": (
                    f"User query: {query}\n\n"
                    f"Available skills:\n{menu}\n\n"
                    "Which skills are needed to answer this query?\n"
                    "Reply with a JSON array of skill names only, no explanation.\n"
                    'Example: ["gene_genomics"] or ["gene_genomics", "gene_proteomics"]'
                ),
            }],
            options={"temperature": 0},
        )
        selected = json.loads(response.message.content.strip())
        if verbose:
            print(f"🗺️  Router selected: {selected}")
        return selected

    # ── Per-skill tool filtering ──────────────────────────────────────────────

    def _tools_for_skill(self, skill: dict) -> list:
        """Filter to only the tools listed in the skill's frontmatter.

        Ensures gene_genomics can never accidentally trigger a UniProt call,
        and vice versa.
        """
        allowed = set(skill["metadata"]["tools"])
        return [t for t in self.tools.schemas if t["function"]["name"] in allowed]

    # ── Single-skill ReAct loop ───────────────────────────────────────────────

    def _run_skill(self, skill: dict, user_query: str, verbose: bool) -> tuple:
        """Run the ReAct loop for one skill. Returns (answer, steps)."""
        skill_name = skill["metadata"]["name"]
        skill_tools = self._tools_for_skill(skill)

        messages = [
            {"role": "system", "content": skill["system_prompt"]},
            {"role": "user", "content": user_query},
        ]
        steps = []

        if verbose:
            print(f"\n{'─' * 60}")
            print(f"📖 Running skill: {skill_name}")
            print(f"{'─' * 60}")

        for iteration in range(5):
            if verbose:
                print(f"\n[Round {iteration + 1}] Reasoning...")

            response = ollama.chat(
                model="qwen3:30b-a3b",
                messages=messages,
                tools=skill_tools,
                options={"temperature": 0.6, "num_ctx": 8192},
            )
            message = response.message

            if message.tool_calls:
                for tool_call in message.tool_calls:
                    tool_name = tool_call.function.name
                    tool_args = tool_call.function.arguments
                    if isinstance(tool_args, str):
                        tool_args = json.loads(tool_args) if tool_args else {}

                    if verbose:
                        print(f"\n🔧 Call tool: {tool_name}({tool_args})")

                    steps.append({"type": "tool_call", "tool": tool_name, "args": tool_args})

                    # Delegate execution entirely to the tools layer
                    tool_result = self.tools.execute(tool_name, **tool_args)

                    if verbose:
                        print(f"   ↳ Returned {len(str(tool_result))} characters.")

                    steps.append({"type": "tool_result", "tool": tool_name, "result": tool_result})
                    messages.append({"role": "assistant", "content": "", "tool_calls": message.tool_calls})
                    messages.append({"role": "tool", "content": json.dumps(tool_result, ensure_ascii=False)})

            else:
                final_answer = message.content
                if verbose:
                    print(f"\n📋 [{skill_name}] Result:")
                    print(final_answer)
                steps.append({"type": "final_answer", "content": final_answer})
                return final_answer, steps

        return "Maximum iteration count reached.", steps

    # ── Public entry point ────────────────────────────────────────────────────

    def run(self, user_query: str, verbose: bool = True) -> tuple:
        """Route the query, then run each selected skill in order."""
        if verbose:
            print(f"\n{'=' * 60}")
            print(f"🧬 Query: {user_query}")
            print(f"{'=' * 60}")

        skill_names = self._route(user_query, verbose)
        all_answers, all_steps = [], []

        # Record the router's decision as an explicit step for visualization
        all_steps.append({"type": "skill_route", "selected": skill_names})

        for name in skill_names:
            if name not in self.skills:
                if verbose:
                    print(f"⚠️  Unknown skill '{name}', skipping.")
                continue
            skill = self.skills[name]
            # Mark the boundary where this skill's steps begin
            all_steps.append({
                "type": "skill_start",
                "skill": name,
                "description": skill["metadata"]["description"],
            })
            answer, steps = self._run_skill(skill, user_query, verbose)
            all_answers.append(answer)
            all_steps.extend(steps)

        return "\n\n---\n\n".join(all_answers), all_steps
