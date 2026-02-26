"""
tools/sc_integrate.py
---------------------
Single-cell RNA-seq QC and multi-dataset integration.
"""

import gzip
import re
from pathlib import Path

import scanpy as sc


def _parse_soft_subtype(soft_path: Path) -> dict[str, str]:
    """
    Parse SOFT file and return {gsm_id: subtype} mapping.
    Handles three field styles:
      - 'cancer type: ER+ tumour'        (GSE161529)
      - 'clinical_subtype: TNBC'         (GSE176078)
      - 'er status' + 'her2 status'      (GSE245601)
    Returns standardised labels: ER+, PR+, HER2+, TNBC, Normal, or Unknown.
    """
    _DIRECT_MAP = {
        "er+ tumour": "ER+",
        "pr+ tumour": "PR+",
        "her2+ tumour": "HER2+",
        "triple negative tumour": "TNBC",
        "triple negative brca1 tumour": "TNBC",
        "brca1 pre-neoplastic": "Normal",
        "normal": "Normal",
        "er+": "ER+",
        "pr+": "PR+",
        "her2+": "HER2+",
        "her2+/er+": "HER2+",
        "tnbc": "TNBC",
    }

    open_fn = gzip.open if str(soft_path).endswith(".gz") else open
    gsm_fields: dict[str, dict[str, str]] = {}
    current_gsm = None

    with open_fn(str(soft_path), "rt") as f:
        for line in f:
            m = re.match(r"\^SAMPLE = (GSM\d+)", line)
            if m:
                current_gsm = m.group(1)
                gsm_fields.setdefault(current_gsm, {})
                continue
            if current_gsm and "Sample_characteristics" in line:
                if ":" in line:
                    _, kv = line.split("=", 1)
                    kv = kv.strip()
                    if ":" in kv:
                        k, v = kv.split(":", 1)
                        gsm_fields[current_gsm][k.strip().lower()] = v.strip()

    result = {}
    for gsm, fields in gsm_fields.items():
        subtype = None
        for key in ("cancer type", "clinical_subtype"):
            if key in fields:
                subtype = _DIRECT_MAP.get(fields[key].lower(), "Unknown")
                break

        if subtype is None:
            er = fields.get("er status", "").lower()
            her2 = fields.get("her2 status", "").lower()
            if "positive" in er and "positive" in her2:
                subtype = "HER2+"
            elif "positive" in er:
                subtype = "ER+"
            elif "positive" in her2:
                subtype = "HER2+"
            elif "negative" in er and "negative" in her2:
                subtype = "TNBC"
            elif fields.get("is tumor", "").strip().upper() == "N":
                subtype = "Normal"
            else:
                subtype = "Unknown"

        result[gsm] = subtype

    return result


def sc_qc(
    accession: str,
    output_dir: str = "./sc_data",
    min_genes: int = 200,
    max_genes: int = 6000,
    min_counts: int = 500,
    max_pct_mito: float = 20.0,
    run_doublet: bool = True,
    doublet_threshold: float = 0.25,
) -> dict:
    """
    Run QC on all h5ad files for a GEO Series.
    Steps:
      1. Compute QC metrics (n_genes, total_counts, pct_counts_mito)
      2. Filter cells by min_genes, max_genes, pct_mito
      3. Optionally run Scrublet doublet detection and remove predicted doublets
      4. Save filtered h5ad as {sample_id}_qc.h5ad
    """
    import scrublet as scr
    import numpy as np

    accession = accession.strip().upper()
    dataset_dir = Path(output_dir) / accession
    if not dataset_dir.exists():
        return {"error": f"Dataset directory not found: {dataset_dir}"}

    sample_dirs = [d for d in sorted(dataset_dir.iterdir()) if d.is_dir()]
    if not sample_dirs:
        return {"error": "No sample directories found. Run convert_geo_to_h5ad first."}

    results = {"accession": accession, "samples": [], "errors": []}

    for sample_dir in sample_dirs:
        sample_id = sample_dir.name
        h5ad_path = sample_dir / f"{sample_id}.h5ad"
        qc_path = sample_dir / f"{sample_id}_qc.h5ad"

        if not h5ad_path.exists():
            results["errors"].append(f"{sample_id}: h5ad not found, run convert_geo_to_h5ad first")
            continue

        if qc_path.exists():
            results["samples"].append({
                "sample_id": sample_id,
                "status": "skipped",
                "reason": "qc h5ad already exists",
                "qc_path": str(qc_path),
            })
            continue

        try:
            adata = sc.read_h5ad(str(h5ad_path))
            n_cells_before = adata.n_obs

            adata.var["mt"] = adata.var_names.str.startswith("MT-")
            sc.pp.calculate_qc_metrics(
                adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
            )

            adata = adata[adata.obs["n_genes_by_counts"] >= min_genes]
            adata = adata[adata.obs["n_genes_by_counts"] <= max_genes]
            adata = adata[adata.obs["total_counts"] >= min_counts]
            adata = adata[adata.obs["pct_counts_mt"] <= max_pct_mito]
            n_after_filter = adata.n_obs

            n_doublets_removed = 0
            if run_doublet and adata.n_obs >= 50:
                scrub = scr.Scrublet(adata.X)
                doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)
                if predicted_doublets is None:
                    predicted_doublets = doublet_scores > doublet_threshold
                adata.obs["doublet_score"] = doublet_scores
                adata.obs["predicted_doublet"] = predicted_doublets
                n_doublets_removed = int(np.sum(predicted_doublets))
                adata = adata[~adata.obs["predicted_doublet"]]

            adata.write_h5ad(str(qc_path))

            results["samples"].append({
                "sample_id": sample_id,
                "status": "done",
                "n_cells_before": n_cells_before,
                "n_cells_after_filter": n_after_filter,
                "n_doublets_removed": n_doublets_removed,
                "n_cells_final": adata.n_obs,
                "qc_path": str(qc_path),
            })

        except Exception as e:
            results["errors"].append(f"{sample_id}: {e}")

    return results


def sc_integrate(
    accessions: list[str],
    output_dir: str = "./sc_data",
    output_path: str = "./sc_data/integrated.h5ad",
    subtypes: list[str] | None = None,
) -> dict:
    """
    Merge QC-filtered h5ad files from multiple GEO Series into a single AnnData.
    Adds obs columns: dataset, sample, subtype.
    Barcodes are prefixed with {dataset}_{sample}_ to avoid collisions.
    """
    results = {
        "accessions": accessions,
        "samples_included": [],
        "samples_skipped": [],
        "errors": [],
    }

    adatas = []

    for accession in accessions:
        accession = accession.strip().upper()
        dataset_dir = Path(output_dir) / accession

        soft_files = list(dataset_dir.glob("*.soft.gz")) + list(dataset_dir.glob("*.soft"))
        if not soft_files:
            results["errors"].append(f"{accession}: SOFT file not found")
            continue
        try:
            gsm_subtype = _parse_soft_subtype(soft_files[0])
        except Exception as e:
            results["errors"].append(f"{accession}: SOFT parse failed: {e}")
            continue

        for sample_dir in sorted(d for d in dataset_dir.iterdir() if d.is_dir()):
            sample_id = sample_dir.name
            qc_path = sample_dir / f"{sample_id}_qc.h5ad"

            if not qc_path.exists():
                results["samples_skipped"].append({
                    "sample_id": sample_id,
                    "reason": "no _qc.h5ad found",
                })
                continue

            subtype = gsm_subtype.get(sample_id, "Unknown")
            if subtype == "Unknown":
                results["samples_skipped"].append({
                    "sample_id": sample_id,
                    "reason": "subtype could not be determined",
                })
                continue

            if subtypes and subtype not in subtypes:
                results["samples_skipped"].append({
                    "sample_id": sample_id,
                    "reason": f"subtype {subtype} not in filter {subtypes}",
                })
                continue

            try:
                adata = sc.read_h5ad(str(qc_path))
                adata.obs["dataset"] = accession
                adata.obs["sample"] = sample_id
                adata.obs["subtype"] = subtype
                adata.obs = adata.obs[["dataset", "sample", "subtype"]]
                adatas.append(adata)
                results["samples_included"].append({
                    "sample_id": sample_id,
                    "dataset": accession,
                    "subtype": subtype,
                    "n_cells": adata.n_obs,
                })
            except Exception as e:
                results["errors"].append(f"{sample_id}: {e}")

    if not adatas:
        results["errors"].append("No samples to integrate")
        return results

    for adata in adatas:
        dataset = adata.obs["dataset"].iloc[0]
        sample = adata.obs["sample"].iloc[0]
        adata.obs_names = [f"{dataset}_{sample}_{bc}" for bc in adata.obs_names]

    merged = sc.concat(adatas, join="inner", label=None)

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    merged.write_h5ad(output_path)

    results["n_cells_total"] = merged.n_obs
    results["n_genes"] = merged.n_vars
    results["output_path"] = output_path

    return results
