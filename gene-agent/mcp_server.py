"""
MCP Server — Gene Tools
========================
This file turns the gene query functions into an MCP-compliant tool server.

Key MCP concepts demonstrated here:
  - Server:      the process that exposes tools
  - @list_tools: tells any client "here is what I can do"
  - @call_tool:  actually executes a tool when a client asks
  - stdio:       the transport layer (client and server talk via stdin/stdout)

Once running, ANY MCP-compatible client can use these tools:
  - Your own agent (agent.py, Phase 2)
  - Claude Desktop
  - Cursor
  - Any future agent you build
"""

import asyncio
import ftplib
import json
import os
import re
import tarfile
import time
from pathlib import Path

import requests
import scanpy as sc
from mcp.server import Server
from mcp.server.stdio import stdio_server
from mcp import types


# ─────────────────────────────────────────────
# 1. Tool implementations
#    Identical logic to tools.py — same API calls.
#    The difference: these are now exposed via a standard protocol,
#    not imported as Python functions.
# ─────────────────────────────────────────────

def search_gene_ncbi(gene_name: str, organism: str = "human") -> dict:
    """Query basic gene information via NCBI Entrez API."""
    try:
        search_resp = requests.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
            params={
                "db": "gene",
                "term": f"{gene_name}[Gene Name] AND {organism}[Organism]",
                "retmode": "json",
                "retmax": 1,
            },
            timeout=10,
        )
        id_list = search_resp.json().get("esearchresult", {}).get("idlist", [])
        if not id_list:
            return {"error": f"Gene not found: {gene_name}"}

        gene_id = id_list[0]
        summary_resp = requests.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
            params={"db": "gene", "id": gene_id, "retmode": "json"},
            timeout=10,
        )
        result = summary_resp.json().get("result", {}).get(gene_id, {})
        summary = result.get("summary", "N/A")

        return {
            "gene_id": gene_id,
            "name": result.get("name", "N/A"),
            "full_name": result.get("description", "N/A"),
            "organism": result.get("organism", {}).get("scientificname", "N/A"),
            "chromosome": result.get("chromosome", "N/A"),
            "location": result.get("maplocation", "N/A"),
            "summary": summary[:500] + "..." if len(summary) > 500 else summary,
            "aliases": result.get("otheraliases", "N/A"),
        }
    except Exception as e:
        return {"error": f"NCBI API call failed: {str(e)}"}


def get_uniprot_function(gene_name: str, organism: str = "human") -> dict:
    """Get protein function description via UniProt API."""
    try:
        organism_map = {"human": "Homo sapiens", "mouse": "Mus musculus"}
        resp = requests.get(
            "https://rest.uniprot.org/uniprotkb/search",
            params={
                "query": f"gene:{gene_name} AND organism_name:{organism_map.get(organism, organism)} AND reviewed:true",
                "fields": "gene_names,protein_name,cc_function,cc_disease,go_p",
                "format": "json",
                "size": 1,
            },
            timeout=10,
        )
        results = resp.json().get("results", [])
        if not results:
            return {"error": f"UniProt found no reviewed entry: {gene_name}"}

        entry = results[0]
        function_text, disease_text = "", ""
        for comment in entry.get("comments", []):
            if comment.get("commentType") == "FUNCTION":
                texts = comment.get("texts", [])
                if texts:
                    function_text = texts[0].get("value", "")[:600]
            if comment.get("commentType") == "DISEASE":
                disease_name = comment.get("disease", {}).get("diseaseName", "")
                if disease_name:
                    disease_text += disease_name + "; "

        protein_name = (
            entry.get("proteinDescription", {})
            .get("recommendedName", {})
            .get("fullName", {})
            .get("value", "N/A")
        )
        return {
            "uniprot_id": entry.get("primaryAccession", "N/A"),
            "protein_name": protein_name,
            "function": function_text or "N/A",
            "associated_diseases": disease_text.rstrip("; ") or "N/A",
        }
    except Exception as e:
        return {"error": f"UniProt API call failed: {str(e)}"}


# ── GEO download helpers ──────────────────────────────────────────────────────

_GEO_FTP_HOST = "ftp.ncbi.nlm.nih.gov"
_DATA_EXTS = (".mtx.gz", ".mtx", ".tsv.gz", ".tsv", ".csv.gz", ".csv", ".h5", ".tar.gz")
_SIZE_LIMIT_MB = 2000  # skip individual files larger than this


def _geo_series_base(accession: str) -> str:
    m = re.match(r"GSE(\d+)", accession)
    num = int(m.group(1))
    return f"/geo/series/GSE{num // 1000}nnn/{accession}"


def _geo_ftp_dir(accession: str) -> str:
    return f"{_geo_series_base(accession)}/suppl"


def _list_ftp_files(ftp_dir: str) -> list[str]:
    with ftplib.FTP(_GEO_FTP_HOST, timeout=30) as ftp:
        ftp.login()
        try:
            paths = ftp.nlst(ftp_dir)
        except ftplib.error_perm:
            return []
    return [os.path.basename(p) for p in paths]


def _get_ftp_file_sizes(ftp_dir: str, filenames: list[str]) -> dict[str, float]:
    """Return {filename: size_in_MB} for each file. -1 if size unavailable."""
    sizes: dict[str, float] = {}
    try:
        with ftplib.FTP(_GEO_FTP_HOST, timeout=30) as ftp:
            ftp.login()
            ftp.set_pasv(True)
            ftp.sendcmd("TYPE I")
            for fname in filenames:
                try:
                    sizes[fname] = ftp.size(f"{ftp_dir}/{fname}") / (1024 * 1024)
                except Exception:
                    sizes[fname] = -1.0
    except Exception:
        sizes = {f: -1.0 for f in filenames}
    return sizes


def _fetch_entrez_metadata(accession: str) -> dict:
    """
    Use NCBI Entrez e-utils (REST, no biopython) to fetch series metadata:
    title, summary, n_samples, pubmed_ids, organism, submission date.
    Returns empty dict on any failure.
    """
    try:
        # esearch: accession → internal GDS id
        search = requests.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
            params={"db": "gds", "term": f"{accession}[Accession]", "retmode": "json"},
            timeout=20,
        )
        ids = search.json().get("esearchresult", {}).get("idlist", [])
        if not ids:
            return {}

        # esummary: internal id → rich metadata
        summary = requests.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
            params={"db": "gds", "id": ids[0], "retmode": "json"},
            timeout=20,
        )
        rec = summary.json().get("result", {}).get(ids[0], {})
        return {
            "title": rec.get("title", ""),
            "summary": rec.get("summary", "")[:600],
            "n_samples": rec.get("n_samples", "unknown"),
            "organism": rec.get("taxon", ""),
            "gse": rec.get("accession", accession),
            "pubmed_ids": rec.get("pubmedids", []),
            "submission_date": rec.get("pdat", ""),
        }
    except Exception:
        return {}


def _download_file_with_retry(ftp_dir: str, filename: str, dest: str, max_retries: int = 3) -> None:
    """HTTP download via NCBI HTTPS endpoint with exponential-backoff retry."""
    url = f"https://{_GEO_FTP_HOST}{ftp_dir}/{filename}"
    for attempt in range(max_retries):
        try:
            with requests.get(url, stream=True, timeout=600) as r:
                r.raise_for_status()
                with open(dest, "wb") as f:
                    for chunk in r.iter_content(chunk_size=65536):
                        f.write(chunk)
            return
        except Exception as e:
            if attempt == max_retries - 1:
                raise
            wait = 2 ** attempt
            time.sleep(wait)


def _group_by_sample(filenames: list[str]) -> dict[str, list[str]]:
    groups: dict[str, list[str]] = {}
    for f in filenames:
        m = re.search(r"(GSM\d+)", f, re.IGNORECASE)
        key = m.group(1).upper() if m else "all"
        groups.setdefault(key, []).append(f)
    return groups


def _detect_format(filenames: list[str]) -> str:
    if any(_is_mtx_file(f) for f in filenames):
        return "10x_mtx"
    if any(f.endswith(".tar.gz") for f in lower):
        return "tar_gz"
    if any(f.endswith(".csv.gz") or f.endswith(".csv") for f in lower):
        return "csv"
    if any(f.endswith(".tsv.gz") or f.endswith(".tsv") for f in lower):
        return "tsv"
    if any(f.endswith(".h5") for f in lower):
        return "h5"
    return "unknown"


def _is_mtx_file(fname: str) -> bool:
    fl = fname.lower()
    return fl.endswith(".mtx") or fl.endswith(".mtx.gz")


def _load_10x_mtx(sample_dir: Path, filenames: list[str]) -> sc.AnnData:
    """
    Load a 10x-style MTX directory robustly.

    Unlike sc.read_10x_mtx, this handles:
    - Non-standard filename prefixes (e.g. count_matrix_sparse.mtx)
    - Uncompressed files (.mtx without .gz)
    - Single-column gene files (gene symbols only, no Ensembl IDs)
    - 2-column legacy (gene_id + gene_symbol) and 3-column modern formats
    - Optional metadata.csv attached to adata.obs

    MTX data is stored genes×cells; this function transposes to cells×genes.
    """
    import gzip
    import anndata as ad
    import scipy.io
    import pandas as pd

    # Locate the three required files by pattern matching
    mtx_file = barcodes_file = genes_file = metadata_csv = None
    for fname in sorted(filenames):          # sort for determinism
        fl = fname.lower()
        p = sample_dir / fname
        if not p.exists():
            continue
        if _is_mtx_file(fname):
            mtx_file = p
        elif "barcodes" in fl and fl.endswith((".tsv", ".tsv.gz")):
            barcodes_file = p
        elif ("genes" in fl or "features" in fl) and fl.endswith((".tsv", ".tsv.gz")):
            genes_file = p
        elif fl == "metadata.csv":
            metadata_csv = p

    missing = [n for n, f in [("matrix", mtx_file), ("barcodes", barcodes_file), ("genes/features", genes_file)] if f is None]
    if missing:
        raise ValueError(f"Missing required MTX files: {missing}. Found: {[f.name for f in sample_dir.iterdir()]}")

    # Read sparse matrix (genes × cells) → transpose to cells × genes
    def _open(p: Path):
        return gzip.open(p) if str(p).endswith(".gz") else open(p, "rb")

    with _open(mtx_file) as f:
        mat = scipy.io.mmread(f).T.tocsr()   # shape: cells × genes

    # Read barcodes (always single column)
    barcodes = pd.read_csv(str(barcodes_file), header=None, sep="\t")[0].tolist()

    # Read gene/feature file — handle 1, 2, or 3-column variants
    genes_df = pd.read_csv(str(genes_file), header=None, sep="\t")
    if genes_df.shape[1] == 1:
        gene_names = genes_df[0].tolist()
        var_df = pd.DataFrame(index=gene_names)
    elif genes_df.shape[1] == 2:
        gene_names = genes_df[1].tolist()
        var_df = pd.DataFrame({"gene_ids": genes_df[0].tolist()}, index=gene_names)
    else:
        gene_names = genes_df[1].tolist()
        var_df = pd.DataFrame({
            "gene_ids": genes_df[0].tolist(),
            "feature_types": genes_df[2].tolist(),
        }, index=gene_names)

    adata = ad.AnnData(
        X=mat,
        obs=pd.DataFrame(index=barcodes),
        var=var_df,
    )

    # Attach metadata.csv to obs when present (combined-matrix GEO releases)
    if metadata_csv:
        try:
            meta = pd.read_csv(str(metadata_csv), index_col=0)
            shared = adata.obs_names.intersection(meta.index)
            if len(shared) > 0:
                adata = adata[shared].copy()
                adata.obs = meta.loc[shared]
        except Exception:
            pass  # best-effort

    return adata


def _load_tar_gz(tar_path: Path) -> sc.AnnData:
    """
    Extract a tar.gz archive and load the 10x MTX inside it.
    Detects any .mtx or .mtx.gz file regardless of filename prefix.
    Reports found files in error message if no MTX is located.
    """
    extract_dir = tar_path.parent / "extracted"
    extract_dir.mkdir(exist_ok=True)
    with tarfile.open(tar_path) as tar:
        tar.extractall(extract_dir)

    for root, _, files in os.walk(extract_dir):
        if any(_is_mtx_file(f) for f in files):
            return _load_10x_mtx(Path(root), files)

    # Collect found files for a helpful error message
    found: list[str] = []
    for root, _, files in os.walk(extract_dir):
        rel = os.path.relpath(root, extract_dir)
        found.extend(f"{rel}/{f}" for f in files)
    raise ValueError(f"No .mtx file found inside tar.gz. Extracted: {found[:20]}")


def _load_delimited(file_path: Path, sep: str) -> tuple[sc.AnnData, str]:
    """
    Load a dense CSV/TSV count matrix.
    GEO often stores data as genes×cells; auto-transpose if detected.
    Detection heuristic: rows start with letters (gene names) AND fewer rows than columns.
    """
    adata = sc.read_csv(str(file_path), delimiter=sep)
    warning = ""
    first_obs = list(adata.obs_names[:5])
    looks_like_genes = (
        all(re.match(r"^[A-Za-z]", n) for n in first_obs)
        and adata.n_obs < 40000
        and adata.n_obs < adata.n_vars
    )
    if looks_like_genes:
        adata = adata.T
        warning = "Auto-transposed: detected genes×cells orientation"
    return adata, warning


def _download_soft_file(accession: str, out: Path) -> str | None:
    """
    Download the SOFT family file from /soft/ — contains sample characteristics,
    clinical annotations, and study metadata critical for downstream analysis.
    Returns local path if successful, None otherwise.
    """
    soft_ftp_dir = f"{_geo_series_base(accession)}/soft"
    try:
        soft_files = _list_ftp_files(soft_ftp_dir)
        # Select the family soft file (not XML, not miniml)
        candidates = [
            f for f in soft_files
            if "family" in f.lower() and "xml" not in f.lower() and "miniml" not in f.lower()
        ]
        if not candidates:
            return None
        fname = candidates[0]
        dest = out / fname
        if not dest.exists():
            _download_file_with_retry(soft_ftp_dir, fname, str(dest))
        return str(dest)
    except Exception:
        return None


def download_geo_dataset(accession: str, output_dir: str = "./sc_data") -> dict:
    """
    Download a GEO Series (GSE) dataset and save each sample as a .h5ad file.

    Steps:
      1. Fetch dataset metadata via NCBI Entrez (title, summary, n_samples, organism)
      2. Download the SOFT family file — sample-level annotations (tissue, disease, etc.)
      3. List supplementary files, check sizes, skip files > 500 MB
      4. Group files by GSM sample ID; download each sample's files with retry
      5. Detect format (10x MTX / tar.gz / CSV / TSV), load with scanpy, save as .h5ad

    Each GSM sample saved to: output_dir/accession/GSMXXXXX.h5ad
    SOFT file saved to:        output_dir/accession/

    Args:
        accession:  GEO Series ID, e.g. "GSE12345"
        output_dir: Base directory for output (default: ./sc_data)

    Returns:
        dict with keys:
          accession    – normalised accession
          metadata     – Entrez dataset info (title, summary, n_samples, organism, pubmed_ids)
          soft_path    – local path to SOFT family file (or null)
          samples      – list of {sample_id, h5ad_path, n_cells, n_genes, format[, warning]}
          skipped      – files skipped due to size limit with their sizes in MB
          errors       – list of error strings for samples that failed
    """
    accession = accession.strip().upper()
    if not re.match(r"^GSE\d+$", accession):
        return {"error": f"Invalid accession '{accession}'. Expected format: GSExxxxx"}

    out = Path(output_dir) / accession
    out.mkdir(parents=True, exist_ok=True)

    # Step 1: Entrez metadata (non-fatal — proceed even if it fails)
    metadata = _fetch_entrez_metadata(accession)

    # Step 2: SOFT family file
    soft_path = _download_soft_file(accession, out)

    # Step 3: List supplementary files
    ftp_dir = _geo_ftp_dir(accession)
    try:
        all_files = _list_ftp_files(ftp_dir)
    except Exception as e:
        return {"error": f"FTP listing failed: {e}", "metadata": metadata}

    if not all_files:
        return {"error": f"No supplementary files found for {accession}", "metadata": metadata}

    data_files = [f for f in all_files if any(f.lower().endswith(e) for e in _DATA_EXTS)]
    if not data_files:
        return {
            "error": "No recognised data files found",
            "all_suppl_files": all_files,
            "metadata": metadata,
        }

    # Check file sizes before downloading — skip oversize files
    sizes = _get_ftp_file_sizes(ftp_dir, data_files)
    skipped = {f: f"{sizes[f]:.1f} MB" for f in data_files if 0 < sizes[f] > _SIZE_LIMIT_MB}
    data_files = [f for f in data_files if sizes.get(f, 0) <= _SIZE_LIMIT_MB or sizes.get(f, 0) < 0]

    groups = _group_by_sample(data_files)
    results: dict = {
        "accession": accession,
        "metadata": metadata,
        "soft_path": soft_path,
        "samples": [],
        "skipped": skipped,
        "errors": [],
    }

    for sample_id, sample_files in sorted(groups.items()):
        sample_dir = out / sample_id
        sample_dir.mkdir(exist_ok=True)
        try:
            for fname in sample_files:
                dest = sample_dir / fname
                if not dest.exists():
                    _download_file_with_retry(ftp_dir, fname, str(dest))

            fmt = _detect_format(sample_files)
            warning = ""

            if fmt == "10x_mtx":
                adata = _load_10x_mtx(sample_dir, sample_files)
            elif fmt == "tar_gz":
                tar_fname = next(f for f in sample_files if f.lower().endswith(".tar.gz"))
                adata = _load_tar_gz(sample_dir / tar_fname)
            elif fmt == "csv":
                csv_fname = next(f for f in sample_files if "csv" in f.lower())
                adata, warning = _load_delimited(sample_dir / csv_fname, sep=",")
            elif fmt == "tsv":
                # Skip standalone barcodes/features files; pick the matrix TSV
                tsv_fname = next(
                    f for f in sample_files
                    if "tsv" in f.lower() and "barcodes" not in f.lower() and "features" not in f.lower()
                )
                adata, warning = _load_delimited(sample_dir / tsv_fname, sep="\t")
            else:
                results["errors"].append(f"{sample_id}: unrecognised format. Files: {sample_files}")
                continue

            adata.obs_names_make_unique()
            adata.var_names_make_unique()

            h5ad_path = str(out / f"{sample_id}.h5ad")
            adata.write_h5ad(h5ad_path)

            entry = {
                "sample_id": sample_id,
                "h5ad_path": h5ad_path,
                "n_cells": adata.n_obs,
                "n_genes": adata.n_vars,
                "format": fmt,
            }
            if warning:
                entry["warning"] = warning
            results["samples"].append(entry)

        except Exception as e:
            results["errors"].append(f"{sample_id}: {e}")

    return results


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
            name="download_geo_dataset",
            description=(
                "Download a single-cell RNA-seq dataset from NCBI GEO by Series accession "
                "(GSExxxxx) and save each sample as an AnnData .h5ad file. "
                "Supports 10x MTX format (separate files or tar.gz), CSV, and TSV matrices. "
                "Returns paths to saved h5ad files and cell/gene counts per sample."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "accession": {
                        "type": "string",
                        "description": "GEO Series accession number, e.g. 'GSE12345'",
                    },
                    "output_dir": {
                        "type": "string",
                        "description": "Directory to save h5ad files (default: ./sc_data)",
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
        "download_geo_dataset": download_geo_dataset,
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
