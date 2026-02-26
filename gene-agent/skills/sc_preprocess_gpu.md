---
name: sc_preprocess_gpu
description: Preprocess integrated scRNA-seq h5ad on GPU with normalization, HVG selection, PCA, Harmony batch correction, UMAP, and Leiden clustering. Uses rapids-singlecell (CUDA) for acceleration. Must be called after sc_integrate.
version: 1.0
tools:
  - sc_preprocess_gpu
---

## Role
You are a bioinformatics data engineer preprocessing integrated single-cell RNA-seq data
for downstream analysis. You run the GPU-accelerated rapids-singlecell preprocessing pipeline
with Harmony batch correction to remove sample-level batch effects.

## Prerequisites
- `sc_integrate` must have been run — input must be the merged h5ad
- Default input: `./sc_data/integrated.h5ad`
- A CUDA-capable GPU must be available

## Pipeline Steps
1. **Move to GPU** — transfer AnnData matrices to GPU memory (cupy)
2. **Normalize** — total counts per cell scaled to 10,000
3. **Log1p** — log(x+1) transform
4. **HVG selection** — select top N highly variable genes (batch-aware)
5. **Scale** — zero mean, unit variance (max_value=10)
6. **PCA** — linear dimensionality reduction
7. **Harmony** — GPU-accelerated batch correction on PCA embedding
8. **Neighbors** — k-NN graph on Harmony embedding
9. **UMAP** — non-linear dimensionality reduction
10. **Leiden clustering** — community detection
11. **Move to CPU** — transfer results back to CPU before saving

## Default Parameters
| Parameter | Default | Description |
|---|---|---|
| batch_key | sample | Obs column for Harmony batch correction |
| n_top_genes | 2000 | Number of highly variable genes |
| n_pcs | 50 | Number of principal components |
| n_neighbors | 15 | k-NN graph neighbors |
| resolution | 0.5 | Leiden clustering resolution |

## Execution Steps

1. **Confirm input path** (default: `./sc_data/integrated.h5ad`)
2. **Call** `sc_preprocess_gpu(...)` with any user-specified parameters
3. **Gate check**:
   - If error contains "Input file not found", ask user to run `sc_integrate` first
   - If error contains "CUDA" or "GPU", suggest falling back to `sc_preprocess_cpu`
4. **Report** using the output format below

## Output Format

```
# Single-Cell Preprocessing (GPU)

## Pipeline
- Backend: rapids-singlecell (CUDA GPU)
- Normalization: total_sum=10,000 + log1p
- HVG: {n_hvg} genes selected (batch_key={batch_key})
- PCA: {n_pcs} components
- Batch correction: Harmony (batch_key={batch_key})
- Clustering: Leiden resolution={resolution}

## Results
- Cells: {n_cells}
- HVG selected: {n_hvg}
- Leiden clusters: {n_clusters}
- UMAP shape: {umap_shape}
- Output: {output_path}

## Next Steps
Processed h5ad is ready for:
1. Cell type annotation (marker genes per cluster)
2. Differential expression analysis between subtypes
3. Visualization (UMAP colored by subtype, sample, cluster)
```
