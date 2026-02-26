"""
tools/sc_preprocess_cpu.py
--------------------------
Single-cell RNA-seq preprocessing (CPU): normalization, HVG, PCA, Harmony, UMAP, clustering.
Uses standard scanpy + harmonypy (CPU-only).
"""

from pathlib import Path

import scanpy as sc


def sc_preprocess_cpu(
    input_path: str = "./sc_data/integrated.h5ad",
    output_path: str = "./sc_data/integrated_processed.h5ad",
    batch_key: str = "sample",
    n_top_genes: int = 2000,
    n_pcs: int = 50,
    n_neighbors: int = 15,
    resolution: float = 0.5,
) -> dict:
    """
    Preprocess integrated h5ad:
      1. Normalize + log1p
      2. HVG selection (batch-aware)
      3. Scale
      4. PCA
      5. Harmony batch correction
      6. Neighbors + UMAP on Harmony embedding
      7. Leiden clustering
    """
    if not Path(input_path).exists():
        return {"error": f"Input file not found: {input_path}"}

    adata = sc.read_h5ad(input_path)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, batch_key=batch_key)
    adata.raw = adata
    adata = adata[:, adata.var["highly_variable"]]

    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=n_pcs)

    sc.external.pp.harmony_integrate(adata, key=batch_key, basis="X_pca", adjusted_basis="X_pca_harmony")

    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X_pca_harmony")
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=resolution)

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path)

    return {
        "input_path": input_path,
        "output_path": output_path,
        "n_cells": adata.n_obs,
        "n_hvg": int(adata.var.shape[0]),
        "n_pcs": n_pcs,
        "batch_key": batch_key,
        "n_clusters": int(adata.obs["leiden"].nunique()),
        "umap_shape": list(adata.obsm["X_umap"].shape),
    }
