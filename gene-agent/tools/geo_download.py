"""
tools/geo_download.py
---------------------
GEO raw data download and h5ad conversion functions.
"""

import gzip
import os
import re
import shutil
import tarfile
import time
from pathlib import Path
from urllib.parse import urlparse

import requests
import scanpy as sc

_GEO_FTP_HOST = "ftp.ncbi.nlm.nih.gov"
_DATA_EXTS = (".mtx.gz", ".mtx", ".tsv.gz", ".tsv", ".csv.gz", ".csv", ".h5", ".tar.gz")
_SIZE_LIMIT_MB = 2000


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
        except Exception:
            if attempt == max_retries - 1:
                raise
            time.sleep(2 ** attempt)


def _download_url_with_retry(ftp_url: str, dest: str, max_retries: int = 3) -> None:
    """Download a GEO FTP URL (ftp:// or https://) with exponential-backoff retry."""
    url = ftp_url.replace("ftp://", "https://", 1) if ftp_url.startswith("ftp://") else ftp_url
    for attempt in range(max_retries):
        try:
            with requests.get(url, stream=True, timeout=600) as r:
                r.raise_for_status()
                with open(dest, "wb") as f:
                    for chunk in r.iter_content(chunk_size=65536):
                        f.write(chunk)
            return
        except Exception:
            if attempt == max_retries - 1:
                raise
            time.sleep(2 ** attempt)


def _group_by_sample(filenames: list[str]) -> dict[str, list[str]]:
    groups: dict[str, list[str]] = {}
    for f in filenames:
        m = re.search(r"(GSM\d+)", f, re.IGNORECASE)
        key = m.group(1).upper() if m else "all"
        groups.setdefault(key, []).append(f)
    return groups


def _detect_format(filenames: list[str]) -> str:
    lower = [f.lower() for f in filenames]
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
    Handles non-standard prefixes, uncompressed files, 1/2/3-column gene files.
    """
    import anndata as ad
    import scipy.io
    import pandas as pd

    mtx_file = barcodes_file = genes_file = metadata_csv = None
    for fname in sorted(filenames):
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

    def _open(p: Path):
        return gzip.open(p) if str(p).endswith(".gz") else open(p, "rb")

    with _open(mtx_file) as f:
        mat = scipy.io.mmread(f).T.tocsr()

    barcodes = pd.read_csv(str(barcodes_file), header=None, sep="\t")[0].tolist()

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

    if metadata_csv:
        try:
            meta = pd.read_csv(str(metadata_csv), index_col=0)
            shared = adata.obs_names.intersection(meta.index)
            if len(shared) > 0:
                adata = adata[shared].copy()
                adata.obs = meta.loc[shared]
        except Exception:
            pass

    return adata


def _load_tar_gz(tar_path: Path) -> sc.AnnData:
    """Extract a tar.gz archive and load the 10x MTX inside it."""
    extract_dir = tar_path.parent / "extracted"
    extract_dir.mkdir(exist_ok=True)
    with tarfile.open(tar_path) as tar:
        tar.extractall(extract_dir)

    for root, _, files in os.walk(extract_dir):
        if any(_is_mtx_file(f) for f in files):
            return _load_10x_mtx(Path(root), files)

    found: list[str] = []
    for root, _, files in os.walk(extract_dir):
        rel = os.path.relpath(root, extract_dir)
        found.extend(f"{rel}/{f}" for f in files)
    raise ValueError(f"No .mtx file found inside tar.gz. Extracted: {found[:20]}")


def _load_delimited(file_path: Path, sep: str) -> tuple[sc.AnnData, str]:
    """Load a dense CSV/TSV count matrix. Auto-transposes if genes×cells detected."""
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


def _load_and_save_sample(sample_id: str, sample_dir: Path, filenames: list[str]) -> dict:
    """Detect format, load AnnData, make names unique, save as h5ad."""
    fmt = _detect_format(filenames)
    warning = ""

    if fmt == "10x_mtx":
        adata = _load_10x_mtx(sample_dir, filenames)
    elif fmt == "tar_gz":
        tar_fname = next(f for f in filenames if f.lower().endswith(".tar.gz"))
        adata = _load_tar_gz(sample_dir / tar_fname)
    elif fmt == "h5":
        h5_fname = next(f for f in filenames if f.lower().endswith(".h5"))
        adata = sc.read_10x_h5(str(sample_dir / h5_fname))
    elif fmt == "csv":
        csv_fname = next(f for f in filenames if "csv" in f.lower())
        adata, warning = _load_delimited(sample_dir / csv_fname, sep=",")
    elif fmt == "tsv":
        tsv_fname = next(
            f for f in filenames
            if "tsv" in f.lower() and "barcodes" not in f.lower() and "features" not in f.lower()
        )
        adata, warning = _load_delimited(sample_dir / tsv_fname, sep="\t")
    else:
        raise ValueError(f"Unrecognised format. Files: {filenames}")

    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    h5ad_path = str(sample_dir / f"{sample_id}.h5ad")
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
    return entry


def download_geo_raw(accession: str, output_dir: str = "./sc_data", gsm_filter: list = None) -> dict:
    """
    Download raw supplementary files for a GEO Series (GSE) without loading or converting.

    Strategy (primary → fallback):
      PRIMARY  — Parse SOFT family file for per-GSM FTP links, download each GSM individually.
      FALLBACK — If SOFT parsing yields no links, fall back to series-level /suppl/ files.
    """
    from tools.geo_search import (
        _fetch_entrez_metadata, _download_soft_file, _parse_soft_gsm_links,
        _geo_ftp_dir, _list_ftp_files, _get_ftp_file_sizes,
    )

    accession = accession.strip().upper()
    if not re.match(r"^GSE\d+$", accession):
        return {"error": f"Invalid accession '{accession}'. Expected format: GSExxxxx"}

    gsm_whitelist = {g.strip().upper() for g in gsm_filter} if gsm_filter else None
    out = Path(output_dir) / accession
    out.mkdir(parents=True, exist_ok=True)

    metadata = _fetch_entrez_metadata(accession)
    soft_path = _download_soft_file(accession, out)

    results: dict = {
        "accession": accession,
        "metadata": metadata,
        "soft_path": soft_path,
        "samples": [],
        "skipped": {},
        "errors": [],
    }

    # Download series-level supplementary files (e.g. shared features.tsv.gz)
    ftp_suppl_dir = _geo_ftp_dir(accession)
    try:
        suppl_files = _list_ftp_files(ftp_suppl_dir)
        for fname in suppl_files:
            if fname.lower().endswith("_raw.tar"):
                continue
            dest = out / fname
            if not dest.exists():
                try:
                    _download_file_with_retry(ftp_suppl_dir, fname, str(dest))
                except Exception as e:
                    results["errors"].append(f"series-level {fname}: {e}")
    except Exception:
        pass

    # PRIMARY: per-GSM links from SOFT
    gsm_links: dict[str, list[str]] = {}
    if soft_path:
        try:
            gsm_links = _parse_soft_gsm_links(soft_path)
        except Exception as e:
            results["errors"].append(f"SOFT parse failed: {e}")

    if gsm_links:
        for gsm, ftp_urls in sorted(gsm_links.items()):
            if gsm_whitelist and gsm not in gsm_whitelist:
                continue
            data_urls = [u for u in ftp_urls if any(u.lower().endswith(e) for e in _DATA_EXTS)]
            if not data_urls:
                continue

            sample_dir = out / gsm
            sample_dir.mkdir(exist_ok=True)
            try:
                checked_urls = []
                for url in data_urls:
                    fname = os.path.basename(urlparse(url).path)
                    ftp_dir = os.path.dirname(urlparse(url).path)
                    sizes = _get_ftp_file_sizes(ftp_dir, [fname])
                    sz = sizes.get(fname, -1)
                    if sz > _SIZE_LIMIT_MB:
                        results["skipped"][fname] = f"{sz:.1f} MB"
                    else:
                        checked_urls.append((url, fname, sz))

                if not checked_urls:
                    results["errors"].append(f"{gsm}: all files exceeded size limit")
                    continue

                downloaded_files = []
                total_mb = 0.0
                for url, fname, sz in checked_urls:
                    dest = sample_dir / fname
                    if not dest.exists():
                        _download_url_with_retry(url, str(dest))
                    downloaded_files.append(fname)
                    if sz > 0:
                        total_mb += sz

                results["samples"].append({
                    "sample_id": gsm,
                    "sample_dir": str(sample_dir),
                    "files": downloaded_files,
                    "total_size_mb": round(total_mb, 1),
                })
            except Exception as e:
                results["errors"].append(f"{gsm}: {e}")

    # FALLBACK: series-level /suppl/ files
    else:
        ftp_dir = _geo_ftp_dir(accession)
        try:
            all_files = _list_ftp_files(ftp_dir)
        except Exception as e:
            results["errors"].append(f"FTP listing failed: {e}")
            return results

        if not all_files:
            results["errors"].append(f"No supplementary files found for {accession}")
            return results

        data_files = [f for f in all_files if any(f.lower().endswith(e) for e in _DATA_EXTS)]
        if not data_files:
            results["errors"].append("No recognised data files found")
            results["all_suppl_files"] = all_files
            return results

        sizes = _get_ftp_file_sizes(ftp_dir, data_files)
        skipped = {f: f"{sizes[f]:.1f} MB" for f in data_files if sizes.get(f, 0) > _SIZE_LIMIT_MB}
        results["skipped"].update(skipped)
        data_files = [f for f in data_files if sizes.get(f, 0) <= _SIZE_LIMIT_MB or sizes.get(f, 0) < 0]

        for sample_id, sample_files in sorted(_group_by_sample(data_files).items()):
            sample_dir = out / sample_id
            sample_dir.mkdir(exist_ok=True)
            try:
                for fname in sample_files:
                    dest = sample_dir / fname
                    if not dest.exists():
                        _download_file_with_retry(ftp_dir, fname, str(dest))
                total_mb = sum(sizes.get(f, 0) for f in sample_files if sizes.get(f, 0) > 0)
                results["samples"].append({
                    "sample_id": sample_id,
                    "sample_dir": str(sample_dir),
                    "files": sample_files,
                    "total_size_mb": round(total_mb, 1),
                })
            except Exception as e:
                results["errors"].append(f"{sample_id}: {e}")

    return results


def convert_geo_to_h5ad(accession: str, output_dir: str = "./sc_data") -> dict:
    """
    Convert previously downloaded raw GEO files into AnnData .h5ad files.
    Must be called after download_geo_raw().
    """
    accession = accession.strip().upper()
    if not re.match(r"^GSE\d+$", accession):
        return {"error": f"Invalid accession '{accession}'. Expected format: GSExxxxx"}

    dataset_dir = Path(output_dir) / accession
    if not dataset_dir.exists():
        return {"error": f"Dataset directory not found: {dataset_dir}. Run download_geo_raw first."}

    results: dict = {
        "accession": accession,
        "converted": [],
        "skipped": [],
        "errors": [],
    }

    sample_dirs = [
        d for d in sorted(dataset_dir.iterdir())
        if d.is_dir() and d.name != "extracted"
    ]

    if not sample_dirs:
        results["errors"].append("No sample directories found. Run download_geo_raw first.")
        return results

    series_features_local = None
    for f in dataset_dir.iterdir():
        if f.is_file() and re.search(r"features?", f.name, re.I) and f.name.endswith((".tsv.gz", ".tsv")):
            series_features_local = f
            break

    for sample_dir in sample_dirs:
        sample_id = sample_dir.name
        h5ad_path = sample_dir / f"{sample_id}.h5ad"
        if h5ad_path.exists():
            results["skipped"].append({
                "sample_id": sample_id,
                "reason": "h5ad already exists",
                "h5ad_path": str(h5ad_path),
            })
            continue

        filenames = [
            f.name for f in sample_dir.iterdir()
            if f.is_file() and any(f.name.lower().endswith(e) for e in _DATA_EXTS)
        ]

        has_features = any(
            re.search(r"(features?|genes)", fn, re.I) and fn.endswith((".tsv", ".tsv.gz"))
            for fn in filenames
        )
        if not has_features and series_features_local and series_features_local.exists():
            link_dest = sample_dir / series_features_local.name
            if not link_dest.exists():
                shutil.copy2(str(series_features_local), str(link_dest))
            filenames.append(series_features_local.name)

        if not filenames:
            results["skipped"].append({
                "sample_id": sample_id,
                "reason": "no recognised data files",
            })
            continue

        try:
            entry = _load_and_save_sample(sample_id, sample_dir, filenames)
            results["converted"].append(entry)
        except Exception as e:
            results["errors"].append(f"{sample_id}: {e}")

    return results
