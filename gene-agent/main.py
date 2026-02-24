"""
Entry point — wire up the layers and run.
"""

import os

from tools import GeneTools
from agent import GeneAgent


def main():
    _script_dir = os.path.dirname(os.path.abspath(__file__))
    skills_dir = os.path.join(_script_dir, "skills")

    # Dependency injection: GeneAgent receives GeneTools, not the other way around
    tools = GeneTools()
    agent = GeneAgent(tools=tools, skills_dir=skills_dir)

    queries = [
        "Which chromosome is BRCA1 on?",           # expects: gene_genomics only
        "Give me a full analysis of TP53",          # expects: gene_genomics + gene_proteomics
    ]

    for query in queries:
        answer, steps = agent.run(query, verbose=True)
        tool_calls = len([s for s in steps if s["type"] == "tool_call"])
        print(f"\n✅ Done, executed {tool_calls} tool call(s)\n")


if __name__ == "__main__":
    main()
