"""
tools/sc_preprocess_gpu.py
--------------------------
Single-cell RNA-seq preprocessing (GPU): normalization, HVG, PCA, Harmony, UMAP, clustering.
Uses rapids-singlecell (cupy/CUDA) for GPU-accelerated computation.
Requires a CUDA-capable GPU with rapids-singlecell installed.
"""

from pathlib import Path

import scanpy as sc
import rapids_singlecell as rsc


def sc_preprocess_gpu(
    input_path: str = "./sc_data/integrated.h5ad",
    output_path: str = "./sc_data/integrated_processed.h5ad",
    batch_key: str = "sample",
    n_top_genes: int = 2000,
    n_pcs: int = 50,
    n_neighbors: int = 15,
    resolution: float = 0.5,
) -> dict:
    """
    GPU-accelerated preprocessing of integrated h5ad using rapids-singlecell:
      1. Normalize + log1p
      2. HVG selection (batch-aware)
      3. Scale
      4. PCA
      5. Harmony batch correction
      6. Neighbors + UMAP on Harmony embedding
      7. Leiden clustering
    Data is moved to GPU at the start and back to CPU before saving.
    """
    if not Path(input_path).exists():
        return {"error": f"Input file not found: {input_path}"}

    adata = sc.read_h5ad(input_path)

    rsc.get.anndata_to_GPU(adata)

    rsc.pp.normalize_total(adata, target_sum=1e4)
    rsc.pp.log1p(adata)

    rsc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, batch_key=batch_key)
    adata.raw = adata
    adata = adata[:, adata.var["highly_variable"]]

    rsc.pp.scale(adata, max_value=10)
    rsc.tl.pca(adata, n_comps=n_pcs)

    rsc.pp.harmony_integrate(adata, key=batch_key)

    rsc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X_pca_harmony")
    rsc.tl.umap(adata)
    rsc.tl.leiden(adata, resolution=resolution)

    rsc.get.anndata_to_CPU(adata)

    # Ensure adata.raw.X is on CPU (rapids anndata_to_CPU skips raw)
    if adata.raw is not None:
        import scipy.sparse as sp
        raw_X = adata.raw.X
        if hasattr(raw_X, "get"):          # cupy array
            raw_X = raw_X.get()
        elif hasattr(raw_X, "toarray"):    # cupy sparse
            try:
                raw_X = raw_X.get()
            except AttributeError:
                pass
        if not sp.issparse(raw_X):
            raw_X = sp.csr_matrix(raw_X)
        import anndata
        adata.raw = anndata.Raw(adata, X=raw_X)

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
