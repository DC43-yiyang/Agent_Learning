"""
tools/geo_search.py
-------------------
GEO dataset search and sample listing functions.
"""

import ftplib
import gzip
import json
import os
import re
import time
from pathlib import Path

import requests

_GEO_FTP_HOST = "ftp.ncbi.nlm.nih.gov"
_DATA_EXTS = (".mtx.gz", ".mtx", ".tsv.gz", ".tsv", ".csv.gz", ".csv", ".h5", ".tar.gz")


def _geo_series_base(accession: str) -> str:
    m = re.match(r"GSE(\d+)", accession)
    num = int(m.group(1))
    return f"/geo/series/GSE{num // 1000}nnn/{accession}"


def _geo_ftp_dir(accession: str) -> str:
    return f"{_geo_series_base(accession)}/suppl"


def _list_ftp_files(ftp_dir: str) -> list[str]:
    with ftplib.FTP(_GEO_FTP_HOST, timeout=30) as ftp:
        ftp.login()
        ftp.set_pasv(True)
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
    Use NCBI Entrez e-utils to fetch series metadata:
    title, summary, n_samples, pubmed_ids, organism, submission date.
    """
    try:
        search = requests.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
            params={"db": "gds", "term": f"{accession}[Accession]", "retmode": "json"},
            timeout=20,
        )
        ids = search.json().get("esearchresult", {}).get("idlist", [])
        if not ids:
            return {}

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


def _download_soft_file(accession: str, out: Path) -> str | None:
    """Download the SOFT family file (contains per-sample metadata and FTP links)."""
    from tools.geo_download import _download_file_with_retry
    soft_ftp_dir = f"{_geo_series_base(accession)}/soft"
    try:
        soft_files = _list_ftp_files(soft_ftp_dir)
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


def _parse_soft_gsm_links(soft_path: str) -> dict[str, list[str]]:
    """
    Parse a GEO SOFT family file to extract per-GSM supplementary FTP URLs.
    Returns {gsm_accession: [ftp_url, ...]} (only GSMs with at least one link).
    """
    gsm_links: dict[str, list[str]] = {}
    current_gsm: str | None = None

    opener = gzip.open if str(soft_path).endswith(".gz") else open
    with opener(soft_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if line.startswith("^SAMPLE = "):
                current_gsm = line.split("=", 1)[1].strip()
                gsm_links.setdefault(current_gsm, [])
            elif current_gsm and "!Sample_supplementary_file" in line and "=" in line:
                val = line.split("=", 1)[1].strip()
                if val and val.upper() != "NONE" and val.startswith("ftp://"):
                    gsm_links[current_gsm].append(val)

    return {gsm: links for gsm, links in gsm_links.items() if links}


def search_geo_datasets(query: str, organism: str = "Homo sapiens", max_results: int = 200, output_dir: str = "./sc_data") -> dict:
    """
    Search NCBI GEO for scRNA-seq datasets matching a query and save results to a JSON file.
    """
    _SCRNA_KEYWORDS = {"10x genomics", "10x chromium", "chromium", "10x"}
    _SUPER_KEYWORDS = {"superseries", "subseries listed below", "composed of the subseries"}
    _HARD_EXCLUDE = {
        "affymetrix", "agilent", "beadarray", "beadchip", "genechip",
        "snp chip", "exon array", "gene chip",
        "wes ", "whole exome", "methylation", "mirna", "chip-seq",
        "smart-seq", "smartseq", "cel-seq", "celseq", "mars-seq", "marsseq",
        "c1 system", "strt-seq", "fluidigm", "drop-seq", "indrop", "in-drop",
        "xenium", "cosmx", "visium", "merfish", "seqfish", "slide-seq",
        "single nucleus", "single-nucleus", "snrna-seq", "sn-rna-seq",
        "snpatho", "snuc-seq",
    }
    _SOFT_EXCLUDE = {
        "microarray", "atac-seq", "proteom", "bulk rna",
        "spatial transcriptom", "spatial proteom",
        "snrna", "sn-rna",
    }

    def _is_scrna(title: str, summary: str, platformtitle: str = "") -> bool:
        text = (title + " " + summary + " " + platformtitle).lower()
        return any(kw in text for kw in _SCRNA_KEYWORDS)

    def _is_unsupported(title: str, summary: str, platformtitle: str = "") -> bool:
        text = (title + " " + summary + " " + platformtitle).lower()
        if any(kw in text for kw in _HARD_EXCLUDE):
            return True
        has_scrna = any(kw in text for kw in _SCRNA_KEYWORDS)
        if not has_scrna and any(kw in text for kw in _SOFT_EXCLUDE):
            return True
        return False

    def _is_superseries(summary: str) -> bool:
        return any(kw in summary.lower() for kw in _SUPER_KEYWORDS)

    def _fetch_overall_design(accession: str) -> str:
        try:
            resp = requests.get(
                "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi",
                params={"acc": accession, "targ": "self", "form": "text", "view": "brief"},
                timeout=20,
            )
            parts = []
            for line in resp.text.splitlines():
                if line.startswith("!Series_overall_design"):
                    parts.append(line.split("=", 1)[-1].strip())
            return " ".join(parts)
        except Exception:
            return ""

    try:
        scrna_anchor = "(\"10x Genomics\" OR \"10x Chromium\" OR \"Chromium\")"
        term = f"({query}) AND {scrna_anchor} AND {organism}[Organism] AND GSE[Entry Type]"

        BATCH = 100
        id_list = []
        total_found = 0
        retstart = 0
        while True:
            search_resp = requests.get(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                params={
                    "db": "gds",
                    "term": term,
                    "retmode": "json",
                    "retmax": BATCH,
                    "retstart": retstart,
                    "sort": "relevance",
                },
                timeout=20,
            )
            search_data = search_resp.json().get("esearchresult", {})
            if not total_found:
                total_found = int(search_data.get("count", 0))
            batch_ids = search_data.get("idlist", [])
            if not batch_ids:
                break
            id_list.extend(batch_ids)
            retstart += BATCH
            if len(id_list) >= min(max_results * 3, total_found):
                break
            time.sleep(0.34)

        if not id_list:
            return {
                "query": query, "organism": organism,
                "total_found": 0, "returned": 0,
                "datasets": [], "excluded": [], "saved_to": None,
            }

        summary_data = {}
        for i in range(0, len(id_list), BATCH):
            batch = id_list[i:i + BATCH]
            summary_resp = requests.get(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
                params={"db": "gds", "id": ",".join(batch), "retmode": "json"},
                timeout=30,
            )
            summary_data.update(summary_resp.json().get("result", {}))
            if i + BATCH < len(id_list):
                time.sleep(0.34)

        datasets = []
        excluded = []
        for gds_id in id_list:
            rec = summary_data.get(gds_id, {})
            accession = rec.get("accession", "")
            if not accession.startswith("GSE"):
                continue
            title         = rec.get("title", "")
            summary       = rec.get("summary", "")
            platformtitle = rec.get("platformtitle", "")

            if _is_superseries(summary):
                excluded.append({"accession": accession, "title": title, "reason": "SuperSeries"})
                continue
            if _is_unsupported(title, summary, platformtitle):
                excluded.append({"accession": accession, "title": title, "reason": "unsupported modality"})
                continue

            if not _is_scrna(title, summary, platformtitle):
                time.sleep(0.34)
                overall_design = _fetch_overall_design(accession)
                if not _is_scrna("", overall_design):
                    excluded.append({"accession": accession, "title": title, "reason": "no scRNA-seq signal detected"})
                    continue
                if _is_unsupported("", overall_design):
                    excluded.append({"accession": accession, "title": title, "reason": "unsupported modality (overall_design)"})
                    continue

            datasets.append({
                "accession": accession,
                "title": title,
                "summary": summary[:400],
                "organism": rec.get("taxon", ""),
                "n_samples": rec.get("n_samples", "unknown"),
                "submission_date": rec.get("pdat", ""),
                "pubmed_ids": rec.get("pubmedids", []),
            })

            if len(datasets) >= max_results:
                break

        Path(output_dir).mkdir(parents=True, exist_ok=True)
        safe_query = re.sub(r"[^\w\s-]", "", query).strip().replace(" ", "_")[:50]
        out_path = str(Path(output_dir) / f"search_{safe_query}.json")
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump({
                "query": query, "organism": organism,
                "total_found": total_found, "returned": len(datasets),
                "datasets": datasets, "excluded": excluded,
            }, f, ensure_ascii=False, indent=2)

        return {
            "query": query, "organism": organism,
            "total_found": total_found, "returned": len(datasets),
            "datasets": datasets, "excluded": excluded,
            "saved_to": out_path,
        }

    except Exception as e:
        return {"error": f"GEO search failed: {str(e)}"}


def list_geo_samples(accession: str, output_dir: str = "./sc_data") -> dict:
    """
    List all samples in a GEO Series with their titles and characteristics.
    """
    accession = accession.strip().upper()
    if not re.match(r"^GSE\d+$", accession):
        return {"error": f"Invalid accession '{accession}'. Expected format: GSExxxxx"}

    out = Path(output_dir) / accession
    out.mkdir(parents=True, exist_ok=True)

    soft_path = _download_soft_file(accession, out)
    if not soft_path:
        return {"error": f"Could not download SOFT file for {accession}"}

    samples = []
    current_gsm = None
    current_title = ""
    current_chars = []

    opener = gzip.open if str(soft_path).endswith(".gz") else open
    with opener(soft_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if line.startswith("^SAMPLE = "):
                if current_gsm:
                    samples.append({
                        "sample_id": current_gsm,
                        "title": current_title,
                        "characteristics": current_chars,
                    })
                current_gsm = line.split("=", 1)[1].strip()
                current_title = ""
                current_chars = []
            elif current_gsm:
                if line.startswith("!Sample_title = "):
                    current_title = line.split("=", 1)[1].strip()
                elif line.startswith("!Sample_characteristics_ch1 = "):
                    current_chars.append(line.split("=", 1)[1].strip())

    if current_gsm:
        samples.append({
            "sample_id": current_gsm,
            "title": current_title,
            "characteristics": current_chars,
        })

    return {
        "accession": accession,
        "total_samples": len(samples),
        "samples": samples,
    }
