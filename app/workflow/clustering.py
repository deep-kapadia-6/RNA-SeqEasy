# app/workflow/clustering.py
import scanpy as sc
from anndata import AnnData
from app.models import ClusteringConfig


def select_hvgs(adata: AnnData, config: ClusteringConfig) -> AnnData:
    """
    Select highly variable genes (HVGs) and store them in adata.var['highly_variable'].
    """
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=config.n_highly_variable_genes,
        flavor="seurat_v3",
        inplace=True,
    )
    adata = adata[:, adata.var["highly_variable"]].copy()
    return adata


def run_pca_and_neighbors(adata: AnnData, config: ClusteringConfig) -> AnnData:
    """
    Run PCA and compute nearest neighbors.
    """
    sc.tl.pca(adata, n_comps=config.n_principal_components)
    sc.pp.neighbors(adata, n_pcs=config.n_principal_components)
    return adata


def run_leiden_clustering(adata: AnnData, config: ClusteringConfig) -> AnnData:
    """
    Run Leiden clustering and store cluster labels in adata.obs['leiden'].
    """
    sc.tl.leiden(adata, resolution=config.clustering_resolution, key_added="leiden")
    return adata
