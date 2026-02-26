"""
tools/sc_celltypist.py
----------------------
Cell type annotation using CellTypist.
Expects log-normalized (10k) data. Uses adata.raw if available.
"""

from pathlib import Path

import celltypist
import scanpy as sc


def sc_celltypist(
    input_path: str = "./sc_data/integrated_rapids_processed.h5ad",
    output_path: str = "./sc_data/integrated_celltypist.h5ad",
    model: str = "Cells_Adult_Breast.pkl",
    majority_voting: bool = True,
    over_clustering: str = "leiden_r1",
    use_raw: bool = True,
) -> dict:
    """
    Annotate cell types using CellTypist.

    Steps:
      1. Load h5ad
      2. Restore log-normalised counts from adata.raw (if use_raw=True)
      3. Run celltypist.annotate()
      4. Optionally refine with majority voting over Leiden clusters
      5. Write predicted labels back to adata.obs and save
    """
    if not Path(input_path).exists():
        return {"error": f"Input file not found: {input_path}"}

    adata = sc.read_h5ad(input_path)

    # CellTypist needs log-normalised counts in adata.X
    if use_raw and adata.raw is not None:
        import anndata
        adata_ct = anndata.AnnData(
            X=adata.raw.X,
            obs=adata.obs.copy(),
            var=adata.raw.var.copy(),
        )
    else:
        # Fallback: reload from integrated.h5ad and re-normalize
        integrated_path = str(Path(input_path).parent / "integrated.h5ad")
        if Path(integrated_path).exists():
            import warnings
            warnings.warn(
                f"adata.raw is None in {input_path}. "
                f"Falling back to re-normalizing {integrated_path}."
            )
            adata_ct = sc.read_h5ad(integrated_path)
            sc.pp.normalize_total(adata_ct, target_sum=1e4)
            sc.pp.log1p(adata_ct)
            # Align obs to match adata (same cells, same order)
            adata_ct = adata_ct[adata.obs_names, :]
        else:
            return {
                "error": (
                    "adata.raw is None and no integrated.h5ad found for fallback. "
                    "Please re-run sc_preprocess_gpu to save adata.raw, or provide integrated.h5ad."
                )
            }

    # Download model if not cached
    celltypist.models.download_models(model=model, force_update=False)
    loaded_model = celltypist.models.Model.load(model=model)

    # over_clustering must exist in obs; fall back to None if missing
    oc = over_clustering if (majority_voting and over_clustering in adata_ct.obs.columns) else None
    if majority_voting and oc is None:
        majority_voting = False

    predictions = celltypist.annotate(
        adata_ct,
        model=loaded_model,
        majority_voting=majority_voting,
        over_clustering=oc,
    )

    # Transfer labels back to original adata
    adata.obs["celltypist_cell_type"] = predictions.predicted_labels["predicted_labels"].values
    if majority_voting:
        adata.obs["celltypist_majority_voting"] = predictions.predicted_labels["majority_voting"].values

    # Transfer per-cell probability of top label
    if hasattr(predictions, "probability_matrix"):
        prob = predictions.probability_matrix
        adata.obs["celltypist_conf_score"] = prob.max(axis=1).values

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path)

    n_types = int(adata.obs["celltypist_cell_type"].nunique())
    top_types = (
        adata.obs["celltypist_cell_type"]
        .value_counts()
        .head(10)
        .to_dict()
    )

    result = {
        "input_path": input_path,
        "output_path": output_path,
        "model": model,
        "majority_voting": majority_voting,
        "over_clustering": oc,
        "n_cells": adata.n_obs,
        "n_cell_types": n_types,
        "top_cell_types": top_types,
    }

    if majority_voting:
        result["n_majority_types"] = int(adata.obs["celltypist_majority_voting"].nunique())

    return result
