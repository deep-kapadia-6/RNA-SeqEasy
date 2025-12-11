# app/workflow/analysis.py
import scanpy as sc
from anndata import AnnData
from typing import Dict, Any
from app.models import AnalysisConfig


def run_dimensionality_reduction(adata: AnnData, analysis_config: AnalysisConfig) -> AnnData:
    """
    Run PCA/UMAP/t-SNE as requested in AnalysisConfig.
    Assumes neighbors have already been computed.
    """
    if analysis_config.perform_pca and "X_pca" not in adata.obsm_keys():
        sc.tl.pca(adata)

    if analysis_config.perform_umap and "X_umap" not in adata.obsm_keys():
        sc.tl.umap(adata)

    if analysis_config.perform_tsne and "X_tsne" not in adata.obsm_keys():
        sc.tl.tsne(adata)

    return adata


def run_differential_expression(
    adata: AnnData,
    analysis_config: AnalysisConfig,
    groupby: str = "leiden",
) -> Dict[str, Any]:
    """
    Run differential expression (rank_genes_groups) if requested.

    Returns:
        results dict (Scanpy rank_genes_groups output) or empty dict.
    """
    if not analysis_config.perform_differential_expression:
        return {}

    if groupby not in adata.obs:
        raise ValueError(f"Groupby key '{groupby}' not found in adata.obs")

    sc.tl.rank_genes_groups(adata, groupby=groupby, method="wilcoxon")
    result = adata.uns["rank_genes_groups"]
    return result


def run_pathway_enrichment_placeholder(
    analysis_config: AnalysisConfig,
) -> Dict[str, Any]:
    """
    Placeholder for pathway enrichment analysis.

    Returns:
        Empty dict for now.
    """
    if not analysis_config.perform_pathway_enrichment:
        return {}

    # TODO: integrate gProfiler, Enrichr, gseapy, etc.
    return {"status": "Pathway enrichment not yet implemented"}
