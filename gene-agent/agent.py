"""
SkillLoader + GeneAgent — Orchestration Layer
=============================================
All LLM calls now go through the openai SDK regardless of provider.
Provider selection (local Ollama / DeepSeek / MiniMax) is passed per
request to run(), not stored on the agent — so one GeneAgent instance
handles all providers without reinitialisation.
"""

import json
import os
import re

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

from providers import get_client, PROVIDER_CONFIGS
from skill_loader import SkillLoader


class GeneAgent:
    def __init__(self, mcp_server_script: str, skills_dir: str = "./skills"):
        self.mcp_server_script = mcp_server_script
        self.loader = SkillLoader(skills_dir)
        self.skills = {
            name: self.loader.load(name)
            for name in self.loader.list_skills()
        }
        print(f"✅ Loaded {len(self.skills)} skills: {list(self.skills.keys())}")

    # ── MCP helpers ───────────────────────────────────────────────────────────

    async def _discover_tools(self, session: ClientSession) -> list:
        """Discover tools from MCP server and convert to OpenAI tool schema."""
        response = await session.list_tools()
        return [
            {
                "type": "function",
                "function": {
                    "name": tool.name,
                    "description": tool.description,
                    "parameters": tool.inputSchema,
                },
            }
            for tool in response.tools
        ]

    async def _call_tool(self, session: ClientSession, name: str, arguments: dict) -> dict:
        """Call a tool on the MCP server and return parsed result."""
        result = await session.call_tool(name, arguments)
        return json.loads(result.content[0].text)

    # ── Routing ───────────────────────────────────────────────────────────────

    async def _route(self, query: str, client, model: str, verbose: bool) -> list[str]:
        """Use LLM to select which skills are needed for this query."""
        menu = "\n".join(
            f"- {name}: {skill['metadata']['description']}"
            for name, skill in self.skills.items()
        )
        response = await client.chat.completions.create(
            model=model,
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
            temperature=0,
        )
        content = (response.choices[0].message.content or "").strip()
        match = re.search(r'\[.*?\]', content, re.DOTALL)
        if not match:
            raise ValueError(f"Router did not return a JSON array. Got: {repr(content)}")
        selected = json.loads(match.group())
        if verbose:
            print(f"🗺️  Router selected: {selected}")
        return selected

    # ── Per-skill tool filtering ──────────────────────────────────────────────

    def _tools_for_skill(self, skill: dict, all_tools: list) -> list:
        allowed = set(skill["metadata"]["tools"])
        return [t for t in all_tools if t["function"]["name"] in allowed]

    # ── Single-skill ReAct loop ───────────────────────────────────────────────

    async def _run_skill(
        self,
        session: ClientSession,
        skill: dict,
        all_tools: list,
        user_query: str,
        client,
        model: str,
        verbose: bool,
    ) -> tuple:
        """Run the ReAct loop for one skill using the given LLM client."""
        skill_name = skill["metadata"]["name"]
        skill_tools = self._tools_for_skill(skill, all_tools)

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

            response = await client.chat.completions.create(
                model=model,
                messages=messages,
                tools=skill_tools,
                temperature=0.6,
            )
            msg = response.choices[0].message

            if msg.tool_calls:
                # Add assistant message with all tool calls (once, before results)
                messages.append({
                    "role": "assistant",
                    "content": msg.content or "",
                    "tool_calls": [
                        {
                            "id": tc.id,
                            "type": "function",
                            "function": {
                                "name": tc.function.name,
                                "arguments": tc.function.arguments,
                            },
                        }
                        for tc in msg.tool_calls
                    ],
                })

                # Process each tool call and add its result
                for tool_call in msg.tool_calls:
                    name = tool_call.function.name
                    args = json.loads(tool_call.function.arguments)

                    if verbose:
                        print(f"\n🔧 Call tool: {name}({args})")

                    steps.append({"type": "tool_call", "tool": name, "args": args})

                    tool_result = await self._call_tool(session, name, args)

                    if verbose:
                        print(f"   ↳ Returned {len(str(tool_result))} characters.")

                    steps.append({"type": "tool_result", "tool": name, "result": tool_result})

                    # OpenAI format: tool result must reference the tool_call_id
                    messages.append({
                        "role": "tool",
                        "tool_call_id": tool_call.id,
                        "content": json.dumps(tool_result, ensure_ascii=False),
                    })

            else:
                final_answer = msg.content or ""
                if verbose:
                    print(f"\n📋 [{skill_name}] Result:")
                    print(final_answer)
                steps.append({"type": "final_answer", "content": final_answer})
                return final_answer, steps

        return "Maximum iteration count reached.", steps

    # ── Public entry point ────────────────────────────────────────────────────

    async def run(
        self,
        user_query: str,
        provider_name: str = "local",
        api_key: str | None = None,
        model_name: str | None = None,
        base_url: str | None = None,
        verbose: bool = True,
    ) -> tuple:
        """
        Route the query, then run each selected skill.

        provider_name: "local" | "deepseek" | "minimax"
        api_key:       required for deepseek and minimax, ignored for local
        """
        client, model = get_client(provider_name, api_key, model_name, base_url)
        label = PROVIDER_CONFIGS[provider_name]["label"]

        if verbose:
            print(f"\n{'=' * 60}")
            print(f"🧬 Query: {user_query}")
            print(f"🤖 Provider: {label} ({model})")
            print(f"{'=' * 60}")

        server_params = StdioServerParameters(
            command=os.path.join(os.path.dirname(os.path.abspath(__file__)), ".venv", "bin", "python"),
            args=[self.mcp_server_script],
        )

        async with stdio_client(server_params) as (read_stream, write_stream):
            async with ClientSession(read_stream, write_stream) as session:
                await session.initialize()

                all_tools = await self._discover_tools(session)
                if verbose:
                    tool_names = [t["function"]["name"] for t in all_tools]
                    print(f"🔌 MCP connected · tools: {tool_names}")

                skill_names = await self._route(user_query, client, model, verbose)
                all_answers, all_steps = [], []

                all_steps.append({"type": "skill_route", "selected": skill_names})

                for name in skill_names:
                    if name not in self.skills:
                        if verbose:
                            print(f"⚠️  Unknown skill '{name}', skipping.")
                        continue
                    skill = self.skills[name]
                    all_steps.append({
                        "type": "skill_start",
                        "skill": name,
                        "description": skill["metadata"]["description"],
                    })
                    answer, steps = await self._run_skill(
                        session, skill, all_tools, user_query, client, model, verbose
                    )
                    all_answers.append(answer)
                    all_steps.extend(steps)

        return "\n\n---\n\n".join(all_answers), all_steps
