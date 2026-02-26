"""
Entry point — wire up and run.
asyncio.run() is the single bridge between sync (the OS) and async (the agent).
"""

import asyncio
import os

from agent import GeneAgent


async def main():
    _dir = os.path.dirname(os.path.abspath(__file__))
    agent = GeneAgent(
        mcp_server_script=os.path.join(_dir, "mcp_server.py"),
        skills_dir=os.path.join(_dir, "skills"),
    )

    queries = [
        # "Which chromosome is BRCA1 on?",         # expects: gene_genomics only
        # "Give me a full analysis of TP53",        # expects: both skills
        "Download the GEO dataset GEO176078"
    ]

    for query in queries:
        answer, steps = await agent.run(query, provider_name="local", verbose=True)
        tool_calls = len([s for s in steps if s["type"] == "tool_call"])
        print(f"\n✅ Done, {tool_calls} tool call(s)\n")


if __name__ == "__main__":
    asyncio.run(main())   # ← sync world hands control to async world, once
