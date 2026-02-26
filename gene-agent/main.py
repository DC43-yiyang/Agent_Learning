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
        # ── gene_genomics skill ──────────────────────────────────────────────
        # "Which chromosome is BRCA1 on?",
        # "Give me a full analysis of TP53",

        # ── geo_search skill ─────────────────────────────────────────────────
        # "Search for scRNA-seq breast cancer datasets in GEO",

        # ── sc_download_raw skill ────────────────────────────────────────────
        # "Download scRNA-seq dataset of GSE161529, I only need ER+, HER2+, PR+, TNBC tumour samples",
        # "Download the scRNA-seq dataset of GSE245601, I only need the Tumor_*_Control",
        # "Download the GEO dataset GSE176078",

        # ── sc_convert_to_h5ad skill ─────────────────────────────────────────
        # "Convert the downloaded GSE161529 raw files to h5ad",

        # ── sc_qc skill ──────────────────────────────────────────────────────
        # "Run QC on GSE161529",
        # "Run QC on GSE176078",
        # "Run QC on GSE245601",

        # ── sc_integrate skill ───────────────────────────────────────────────
        # "Integrate GSE161529, GSE176078, GSE245601 into a single h5ad",

        # ── sc_preprocess skill ──────────────────────────────────────────────
        # "Preprocess the integrated h5ad with Harmony batch correction and UMAP clustering",

        # ── sc_celltypist skill ──────────────────────────────────────────────
        "Annotate cell types in the preprocessed h5ad using CellTypist with the Cells_Adult_Breast model and majority voting over Leiden clusters",
    ]

    for query in queries:
        answer, steps = await agent.run(query, provider_name="online", verbose=True)
        tool_calls = len([s for s in steps if s["type"] == "tool_call"])
        print(f"\n✅ Done, {tool_calls} tool call(s)\n")


if __name__ == "__main__":
    asyncio.run(main())   # ← sync world hands control to async world, once
