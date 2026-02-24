"""
SkillLoader + GeneAgent — Orchestration Layer
=============================================
Phase 2 (MCP) changes vs Phase 1:

  Before: agent.py imports tools.py directly (Python import)
  After:  agent.py connects to mcp_server.py via MCP protocol

Key differences:
  - GeneAgent no longer receives a GeneTools instance
  - Instead it receives the path to mcp_server.py
  - Tools are DISCOVERED at runtime via session.list_tools()
    (agent doesn't need to know what tools exist in advance)
  - Tool calls go through session.call_tool() — a network-like call
    over stdio, not a direct Python function call
  - All methods are async because MCP client is async
"""

import asyncio
import json
import os

import ollama
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

from skill_loader import SkillLoader


# ─────────────────────────────────────────────
# GeneAgent  (Orchestration Layer)
# ─────────────────────────────────────────────

class GeneAgent:
    def __init__(self, mcp_server_script: str, skills_dir: str = "./skills"):
        # Path to mcp_server.py — agent will spawn it as a subprocess
        self.mcp_server_script = mcp_server_script
        self.loader = SkillLoader(skills_dir)
        self.skills = {
            name: self.loader.load(name)
            for name in self.loader.list_skills()
        }
        print(f"✅ Loaded {len(self.skills)} skills: {list(self.skills.keys())}")

    # ── MCP helpers ───────────────────────────────────────────────────────────

    async def _discover_tools(self, session: ClientSession) -> list:
        """Ask the MCP server what tools it has, convert to Ollama format.

        This is the core of MCP's value: the agent does not hardcode tool
        definitions. It asks the server at runtime. Add a new tool to
        mcp_server.py and the agent picks it up automatically — zero changes
        here.
        """
        response = await session.list_tools()

        # MCP tool schema  →  Ollama tool schema
        # MCP uses "inputSchema", Ollama uses "parameters" inside "function"
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
        """Call a tool on the MCP server and return the parsed result.

        mcp_server.py returns results as TextContent (a JSON string).
        We parse it back to a dict so the rest of the agent code is unchanged.
        """
        result = await session.call_tool(name, arguments)
        return json.loads(result.content[0].text)

    # ── Routing ───────────────────────────────────────────────────────────────

    async def _route(self, query: str, verbose: bool) -> list[str]:
        """Use LLM to select which skills are needed for this query."""
        menu = "\n".join(
            f"- {name}: {skill['metadata']['description']}"
            for name, skill in self.skills.items()
        )
        client = ollama.AsyncClient()
        response = await client.chat(
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

    def _tools_for_skill(self, skill: dict, all_tools: list) -> list:
        """Filter the MCP-discovered tool list to only what this skill declares.

        Same logic as before — but now `all_tools` comes from the MCP server,
        not from tools.py. The skill frontmatter still controls access.
        """
        allowed = set(skill["metadata"]["tools"])
        return [t for t in all_tools if t["function"]["name"] in allowed]

    # ── Single-skill ReAct loop ───────────────────────────────────────────────

    async def _run_skill(
        self,
        session: ClientSession,
        skill: dict,
        all_tools: list,
        user_query: str,
        verbose: bool,
    ) -> tuple:
        """Run the ReAct loop for one skill.

        The only change from the sync version:
          - ollama.AsyncClient().chat()  instead of  ollama.chat()
          - await self._call_tool()      instead of  self.tools.execute()
        The loop structure is identical.
        """
        skill_name = skill["metadata"]["name"]
        skill_tools = self._tools_for_skill(skill, all_tools)

        messages = [
            {"role": "system", "content": skill["system_prompt"]},
            {"role": "user", "content": user_query},
        ]
        steps = []
        client = ollama.AsyncClient()

        if verbose:
            print(f"\n{'─' * 60}")
            print(f"📖 Running skill: {skill_name}")
            print(f"{'─' * 60}")

        for iteration in range(5):
            if verbose:
                print(f"\n[Round {iteration + 1}] Reasoning...")

            response = await client.chat(
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

                    # ← Key change: call via MCP, not via tools.py
                    tool_result = await self._call_tool(session, tool_name, tool_args)

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

    async def run(self, user_query: str, verbose: bool = True) -> tuple:
        """
        Connect to MCP server → discover tools → route → run skills → disconnect.

        The MCP connection is opened once per query and shared across all
        skill runs. This is why async matters: the connection stays alive
        while we await LLM responses in _run_skill().
        """
        if verbose:
            print(f"\n{'=' * 60}")
            print(f"🧬 Query: {user_query}")
            print(f"{'=' * 60}")

        # Spawn mcp_server.py as a subprocess and connect to it via stdio
        server_params = StdioServerParameters(
            command="uv",
            args=["run", self.mcp_server_script],
        )

        async with stdio_client(server_params) as (read_stream, write_stream):
            async with ClientSession(read_stream, write_stream) as session:
                # Handshake: client and server agree on protocol version
                await session.initialize()

                # Discover what tools the server exposes
                all_tools = await self._discover_tools(session)
                if verbose:
                    tool_names = [t["function"]["name"] for t in all_tools]
                    print(f"🔌 MCP connected · tools discovered: {tool_names}")

                skill_names = await self._route(user_query, verbose)
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
                        session, skill, all_tools, user_query, verbose
                    )
                    all_answers.append(answer)
                    all_steps.extend(steps)

        return "\n\n---\n\n".join(all_answers), all_steps
