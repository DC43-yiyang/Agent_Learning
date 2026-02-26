import time

import rapids_singlecell as rsc
import scanpy as sc
from pathlib import Path

DATA_PATH = Path(
    "/home/compbiowizard_ravenclaw/projects/YiyangDemo/"
    "gene-agent/sc_data/integrated.h5ad"
)

OUTPUT_H5AD = DATA_PATH.with_name("integrated_rapids_processed.h5ad")

# batch key used for Harmony integration
BATCH_KEY = "sample"


def _step(label: str):
    """Print step header and return start time."""
    print(f"\n===== {label} =====")
    return time.perf_counter()


def _done(t0: float):
    """Print elapsed time since t0."""
    print(f"  -> done in {time.perf_counter() - t0:.2f}s")


def main():
    pipeline_start = time.perf_counter()

    t0 = _step("Load AnnData on CPU")
    print(f"Input file: {DATA_PATH}")
    adata = sc.read_h5ad(DATA_PATH)
    print(adata)
    print("obs columns:", adata.obs.columns.tolist())
    _done(t0)

    t0 = _step("Cast X to float32 (required by CuPy sparse)")
    print(f"X dtype before cast: {adata.X.dtype}")
    adata.X = adata.X.astype("float32")
    print(f"X dtype after cast: {adata.X.dtype}")
    _done(t0)

    t0 = _step("Move AnnData to GPU (in-place)")
    rsc.get.anndata_to_GPU(adata)
    print(f"X type on GPU: {type(adata.X)}")
    _done(t0)

    t0 = _step("Normalize / log1p on GPU")
    rsc.pp.normalize_total(adata, target_sum=1e4)
    rsc.pp.log1p(adata)
    _done(t0)

    # seurat flavor runs on log-normalized data and is more stable for large integrated datasets
    # seurat_v3 uses LOESS which becomes numerically unstable with batch_key on 70 samples
    t0 = _step("HVG selection on log-normalized data (seurat, batch-aware)")
    rsc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat", batch_key=BATCH_KEY)
    adata = adata[:, adata.var["highly_variable"]].copy()
    print(f"After HVG filter: {adata.shape}")
    # re-move to GPU after subsetting (copy creates a new CPU-backed object)
    rsc.get.anndata_to_GPU(adata)
    _done(t0)

    t0 = _step("Scale / PCA on GPU")
    rsc.pp.scale(adata, max_value=10)
    rsc.tl.pca(adata, n_comps=50)
    print("X_pca shape:", adata.obsm["X_pca"].shape)
    _done(t0)

    t0 = _step(f"Harmony batch correction (key='{BATCH_KEY}') on GPU")
    rsc.pp.harmony_integrate(adata, key=BATCH_KEY)
    print("X_pca_harmony shape:", adata.obsm["X_pca_harmony"].shape)
    _done(t0)

    t0 = _step("Neighbors / UMAP / Leiden (using harmony embedding)")
    rsc.pp.neighbors(adata, n_neighbors=15, use_rep="X_pca_harmony")
    rsc.tl.umap(adata)
    rsc.tl.leiden(adata, key_added="leiden_r1", resolution=1.0)
    print("Leiden clusters:", adata.obs["leiden_r1"].nunique())
    _done(t0)

    t0 = _step("Move back to CPU and save")
    rsc.get.anndata_to_CPU(adata)
    adata.write_h5ad(OUTPUT_H5AD)
    print(f"Saved to: {OUTPUT_H5AD}")
    _done(t0)

    print("\n===== Summary =====")
    print("obs columns:", adata.obs.columns.tolist())
    print("var shape:", adata.var.shape)
    print("obsm keys:", list(adata.obsm.keys()))
    total = time.perf_counter() - pipeline_start
    print(f"\nTotal pipeline time: {total:.2f}s")
    print("Done. rapids-singlecell full pipeline with Harmony completed successfully.")


if __name__ == "__main__":
    main()
