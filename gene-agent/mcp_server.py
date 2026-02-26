"""
MCP Server — Gene Tools
========================
Protocol layer only: list_tools, call_tool, dispatch.
All business logic lives in tools/.
"""

import asyncio
import json

from mcp.server import Server
from mcp.server.stdio import stdio_server
from mcp import types

from tools.gene_query import search_gene_ncbi, get_uniprot_function
from tools.geo_search import search_geo_datasets, list_geo_samples
from tools.geo_download import download_geo_raw, convert_geo_to_h5ad
from tools.sc_integrate import sc_qc, sc_integrate
from tools.sc_preprocess_cpu import sc_preprocess_cpu
from tools.sc_preprocess_gpu import sc_preprocess_gpu
from tools.sc_celltypist import sc_celltypist


# ─────────────────────────────────────────────
# 2. MCP Server setup
#    app = the server instance, identified by a name.
#    Any client that connects will see this name.
# ─────────────────────────────────────────────

app = Server("gene-tools")


# ─────────────────────────────────────────────
# 3. list_tools handler
#    When a client connects and asks "what can you do?",
#    MCP calls this function and sends the result back.
#    This is how clients discover tools dynamically —
#    no hardcoding on the client side.
# ─────────────────────────────────────────────

@app.list_tools()
async def list_tools() -> list[types.Tool]:
    return [
        types.Tool(
            name="search_gene_ncbi",
            description=(
                "Query basic gene information from the NCBI database, including "
                "gene ID, chromosome location, function summary, and aliases."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "gene_name": {
                        "type": "string",
                        "description": "Gene symbol, e.g. TP53, BRCA1",
                    },
                    "organism": {
                        "type": "string",
                        "description": "Organism name, default: human",
                        "default": "human",
                    },
                },
                "required": ["gene_name"],
            },
        ),
        types.Tool(
            name="get_uniprot_function",
            description=(
                "Query the UniProt database for detailed protein function "
                "description and associated diseases for a given gene."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "gene_name": {
                        "type": "string",
                        "description": "Gene symbol",
                    },
                    "organism": {
                        "type": "string",
                        "description": "Organism name, default: human",
                        "default": "human",
                    },
                },
                "required": ["gene_name"],
            },
        ),
        types.Tool(
            name="search_geo_datasets",
            description=(
                "Search NCBI GEO for datasets matching a query (e.g. 'breast cancer single cell RNA-seq'). "
                "Returns a list of matching GSE accessions with titles, summaries, and sample counts. "
                "Results are saved to a JSON file for later use with download_geo_raw."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "Search query, e.g. 'breast cancer single cell RNA-seq'",
                    },
                    "organism": {
                        "type": "string",
                        "description": "Organism filter (default: Homo sapiens)",
                        "default": "Homo sapiens",
                    },
                    "max_results": {
                        "type": "integer",
                        "description": "Maximum number of results to return (default: 200)",
                        "default": 200,
                    },
                    "output_dir": {
                        "type": "string",
                        "description": "Directory to save the results JSON (default: ./sc_data)",
                        "default": "./sc_data",
                    },
                },
                "required": ["query"],
            },
        ),
        types.Tool(
            name="download_geo_raw",
            description=(
                "Download raw supplementary files for a GEO Series (GSExxxxx) from NCBI GEO. "
                "Each sample (GSM) is saved into its own subdirectory with original files intact. "
                "Does NOT load or convert files — use convert_geo_to_h5ad for that. "
                "Returns sample list with file names and sizes for inspection before conversion. "
                "Use gsm_filter to download only specific samples selected after calling list_geo_samples."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "accession": {
                        "type": "string",
                        "description": "GEO Series accession number, e.g. 'GSE176078'",
                    },
                    "output_dir": {
                        "type": "string",
                        "description": "Root directory to save raw files (default: ./sc_data)",
                        "default": "./sc_data",
                    },
                    "gsm_filter": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "Optional list of GSM IDs to download. If provided, only these samples are downloaded.",
                    },
                },
                "required": ["accession"],
            },
        ),
        types.Tool(
            name="list_geo_samples",
            description=(
                "List all samples in a GEO Series with their titles and characteristics. "
                "Use this BEFORE download_geo_raw to inspect sample metadata and select only "
                "the relevant samples (e.g. only tumour samples, specific subtypes). "
                "Returns sample_id, title, and characteristics for each GSM."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "accession": {
                        "type": "string",
                        "description": "GEO Series accession number, e.g. 'GSE161529'",
                    },
                    "output_dir": {
                        "type": "string",
                        "description": "Directory to cache the SOFT file (default: ./sc_data)",
                        "default": "./sc_data",
                    },
                },
                "required": ["accession"],
            },
        ),
        types.Tool(
            name="sc_qc",
            description=(
                "Run quality control on single-cell RNA-seq h5ad files for a GEO Series. "
                "Filters low-quality cells (min/max genes, % mitochondrial reads) and optionally "
                "detects and removes doublets using Scrublet. "
                "Must be called after convert_geo_to_h5ad. Saves filtered h5ad as {sample_id}_qc.h5ad."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "accession": {
                        "type": "string",
                        "description": "GEO Series accession number, e.g. 'GSE176078'",
                    },
                    "output_dir": {
                        "type": "string",
                        "description": "Root directory where h5ad files are saved (default: ./sc_data)",
                        "default": "./sc_data",
                    },
                    "min_genes": {
                        "type": "integer",
                        "description": "Minimum number of genes per cell (default: 200)",
                        "default": 200,
                    },
                    "max_genes": {
                        "type": "integer",
                        "description": "Maximum number of genes per cell (default: 6000)",
                        "default": 6000,
                    },
                    "min_counts": {
                        "type": "integer",
                        "description": "Minimum total UMI counts per cell (default: 500)",
                        "default": 500,
                    },
                    "max_pct_mito": {
                        "type": "number",
                        "description": "Maximum % mitochondrial reads per cell (default: 20.0)",
                        "default": 20.0,
                    },
                    "run_doublet": {
                        "type": "boolean",
                        "description": "Whether to run Scrublet doublet detection (default: true)",
                        "default": True,
                    },
                    "doublet_threshold": {
                        "type": "number",
                        "description": "Doublet score threshold fallback if auto-detection fails (default: 0.25)",
                        "default": 0.25,
                    },
                },
                "required": ["accession"],
            },
        ),
        types.Tool(
            name="sc_integrate",
            description=(
                "Merge QC-filtered h5ad files from multiple GEO Series into a single AnnData. "
                "Parses subtype from each dataset's SOFT file and adds dataset, sample, subtype metadata. "
                "Only includes samples with a determinable subtype. "
                "Must be called after sc_qc. Saves merged h5ad to output_path."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "accessions": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "List of GEO Series accessions, e.g. ['GSE161529', 'GSE176078']",
                    },
                    "output_dir": {
                        "type": "string",
                        "description": "Root directory containing per-accession subdirectories (default: ./sc_data)",
                        "default": "./sc_data",
                    },
                    "output_path": {
                        "type": "string",
                        "description": "Path to save the merged h5ad (default: ./sc_data/integrated.h5ad)",
                        "default": "./sc_data/integrated.h5ad",
                    },
                    "subtypes": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "Optional subtype filter, e.g. ['ER+', 'TNBC']. If omitted, all subtypes included.",
                    },
                },
                "required": ["accessions"],
            },
        ),
        types.Tool(
            name="sc_preprocess_cpu",
            description=(
                "Preprocess integrated h5ad on CPU: normalize, log1p, HVG selection, scale, PCA, "
                "Harmony batch correction (harmonypy), neighbor graph, UMAP, and Leiden clustering. "
                "Must be called after sc_integrate. Use when no GPU is available."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "input_path": {
                        "type": "string",
                        "description": "Path to integrated h5ad (default: ./sc_data/integrated.h5ad)",
                        "default": "./sc_data/integrated.h5ad",
                    },
                    "output_path": {
                        "type": "string",
                        "description": "Path to save processed h5ad (default: ./sc_data/integrated_processed.h5ad)",
                        "default": "./sc_data/integrated_processed.h5ad",
                    },
                    "batch_key": {
                        "type": "string",
                        "description": "Obs column to use as batch key for Harmony (default: sample)",
                        "default": "sample",
                    },
                    "n_top_genes": {
                        "type": "integer",
                        "description": "Number of highly variable genes to select (default: 2000)",
                        "default": 2000,
                    },
                    "n_pcs": {
                        "type": "integer",
                        "description": "Number of PCs to compute (default: 50)",
                        "default": 50,
                    },
                    "n_neighbors": {
                        "type": "integer",
                        "description": "Number of neighbors for graph (default: 15)",
                        "default": 15,
                    },
                    "resolution": {
                        "type": "number",
                        "description": "Leiden clustering resolution (default: 0.5)",
                        "default": 0.5,
                    },
                },
                "required": [],
            },
        ),
        types.Tool(
            name="sc_preprocess_gpu",
            description=(
                "Preprocess integrated h5ad on GPU: normalize, log1p, HVG selection, scale, PCA, "
                "Harmony batch correction, neighbor graph, UMAP, and Leiden clustering — all GPU-accelerated "
                "via rapids-singlecell (CUDA). Must be called after sc_integrate. Requires a CUDA GPU."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "input_path": {
                        "type": "string",
                        "description": "Path to integrated h5ad (default: ./sc_data/integrated.h5ad)",
                        "default": "./sc_data/integrated.h5ad",
                    },
                    "output_path": {
                        "type": "string",
                        "description": "Path to save processed h5ad (default: ./sc_data/integrated_processed.h5ad)",
                        "default": "./sc_data/integrated_processed.h5ad",
                    },
                    "batch_key": {
                        "type": "string",
                        "description": "Obs column to use as batch key for Harmony (default: sample)",
                        "default": "sample",
                    },
                    "n_top_genes": {
                        "type": "integer",
                        "description": "Number of highly variable genes to select (default: 2000)",
                        "default": 2000,
                    },
                    "n_pcs": {
                        "type": "integer",
                        "description": "Number of PCs to compute (default: 50)",
                        "default": 50,
                    },
                    "n_neighbors": {
                        "type": "integer",
                        "description": "Number of neighbors for graph (default: 15)",
                        "default": 15,
                    },
                    "resolution": {
                        "type": "number",
                        "description": "Leiden clustering resolution (default: 0.5)",
                        "default": 0.5,
                    },
                },
                "required": [],
            },
        ),
        types.Tool(
            name="sc_celltypist",
            description=(
                "Annotate cell types in a preprocessed scRNA-seq h5ad using CellTypist. "
                "Supports majority voting over Leiden clusters to refine per-cell predictions. "
                "Must be called after sc_preprocess_cpu or sc_preprocess_gpu. "
                "Adds celltypist_cell_type, celltypist_majority_voting, and celltypist_conf_score to obs."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "input_path": {
                        "type": "string",
                        "description": "Path to preprocessed h5ad (default: ./sc_data/integrated_rapids_processed.h5ad)",
                        "default": "./sc_data/integrated_rapids_processed.h5ad",
                    },
                    "output_path": {
                        "type": "string",
                        "description": "Path to save annotated h5ad (default: ./sc_data/integrated_celltypist.h5ad)",
                        "default": "./sc_data/integrated_celltypist.h5ad",
                    },
                    "model": {
                        "type": "string",
                        "description": "CellTypist model name (default: Cells_Adult_Breast.pkl)",
                        "default": "Cells_Adult_Breast.pkl",
                    },
                    "majority_voting": {
                        "type": "boolean",
                        "description": "Refine labels by majority vote within clusters (default: true)",
                        "default": True,
                    },
                    "over_clustering": {
                        "type": "string",
                        "description": "Obs column used for majority voting (default: leiden_r1)",
                        "default": "leiden_r1",
                    },
                    "use_raw": {
                        "type": "boolean",
                        "description": "Use adata.raw for annotation if available (default: true)",
                        "default": True,
                    },
                },
                "required": [],
            },
        ),
        types.Tool(
            name="convert_geo_to_h5ad",
            description=(
                "Convert previously downloaded raw GEO files into AnnData .h5ad files. "
                "Must be called after download_geo_raw. "
                "Supports 10x MTX (loose files or tar.gz), CSV, and TSV matrices. "
                "Skips samples that already have an h5ad. "
                "Returns conversion results with cell/gene counts per sample."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "accession": {
                        "type": "string",
                        "description": "GEO Series accession number, e.g. 'GSE176078'",
                    },
                    "output_dir": {
                        "type": "string",
                        "description": "Root directory where raw files were saved (default: ./sc_data)",
                        "default": "./sc_data",
                    },
                },
                "required": ["accession"],
            },
        ),
    ]


# ─────────────────────────────────────────────
# 4. call_tool handler
#    When a client says "call search_gene_ncbi with these args",
#    MCP routes it here. We dispatch to the real function,
#    then return the result as a TextContent (MCP's standard result type).
#
#    Note: MCP results are always a list of Content objects.
#    TextContent is the most common — it just wraps a string.
# ─────────────────────────────────────────────

@app.call_tool()
async def call_tool(name: str, arguments: dict) -> list[types.TextContent]:
    dispatch = {
        "search_gene_ncbi": search_gene_ncbi,
        "get_uniprot_function": get_uniprot_function,
        "search_geo_datasets": search_geo_datasets,
        "list_geo_samples": list_geo_samples,
        "download_geo_raw": download_geo_raw,
        "convert_geo_to_h5ad": convert_geo_to_h5ad,
        "sc_qc": sc_qc,
        "sc_integrate": sc_integrate,
        "sc_preprocess_cpu": sc_preprocess_cpu,
        "sc_preprocess_gpu": sc_preprocess_gpu,
        "sc_celltypist": sc_celltypist,
    }

    if name not in dispatch:
        return [types.TextContent(
            type="text",
            text=json.dumps({"error": f"Unknown tool: {name}"}),
        )]

    result = dispatch[name](**arguments)
    return [types.TextContent(type="text", text=json.dumps(result, ensure_ascii=False))]


# ─────────────────────────────────────────────
# 5. Entry point
#    stdio_server() sets up the stdin/stdout pipes.
#    app.run() starts the event loop and waits for client connections.
#    The server stays alive until the client disconnects.
# ─────────────────────────────────────────────

async def main():
    async with stdio_server() as (read_stream, write_stream):
        await app.run(
            read_stream,
            write_stream,
            app.create_initialization_options(),
        )


if __name__ == "__main__":
    asyncio.run(main())
