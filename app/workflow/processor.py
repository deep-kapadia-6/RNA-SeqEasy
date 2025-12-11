"""Main workflow processor orchestrating the RNA-seq analysis pipeline."""

import logging
from pathlib import Path
from typing import Dict, Any, Optional

import scanpy as sc
from anndata import AnnData

from app.models import WorkflowConfig
from app.workflow.qc import compute_qc_metrics, filter_cells
from app.workflow.normalization import normalize_data
from app.workflow.clustering import (
    select_hvgs,
    run_pca_and_neighbors,
    run_leiden_clustering,
)
from app.workflow.analysis import (
    run_dimensionality_reduction,
    run_differential_expression,
    run_pathway_enrichment_placeholder,
)

logger = logging.getLogger(__name__)


def load_data(input_path: str) -> AnnData:
    """
    Load an AnnData object from .h5ad file.

    Args:
        input_path: Path to .h5ad file

    Returns:
        AnnData object loaded from file

    Raises:
        FileNotFoundError: If file doesn't exist
        OSError: If file cannot be read
        ValueError: If file is not a valid .h5ad file

    Example:
        >>> adata = load_data("data/example.h5ad")
        >>> print(f"Loaded {adata.n_obs} cells")
    """
    input_file = Path(input_path)
    if not input_file.exists():
        error_msg = f"Input file not found: {input_path}"
        logger.error(error_msg)
        raise FileNotFoundError(error_msg)

    logger.info(f"Loading data from {input_path}")
    try:
        adata = sc.read_h5ad(input_path)
        logger.info(f"Loaded {adata.n_obs} cells and {adata.n_vars} genes")
        return adata
    except Exception as e:
        error_msg = f"Failed to load data from {input_path}: {e}"
        logger.error(error_msg)
        raise


def save_adata(adata: AnnData, output_dir: str, filename: str) -> str:
    """
    Save AnnData object to .h5ad in the output directory.

    Args:
        adata: AnnData object to save
        output_dir: Output directory path (will be created if needed)
        filename: Name of the output file (should end with .h5ad)

    Returns:
        Full path to saved file

    Raises:
        OSError: If file cannot be written

    Example:
        >>> path = save_adata(adata, "outputs/", "processed.h5ad")
        >>> print(f"Saved to {path}")
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / filename

    logger.info(f"Saving AnnData to {out_path}")
    try:
        adata.write_h5ad(out_path)
        logger.info(f"Successfully saved {adata.n_obs} cells to {out_path}")
        return str(out_path)
    except Exception as e:
        error_msg = f"Failed to save data to {out_path}: {e}"
        logger.error(error_msg)
        raise OSError(error_msg)


def run_workflow(
    config: WorkflowConfig,
    save_intermediate: bool = True,
) -> Dict[str, Any]:
    """
    Run the full RNA-seq analysis workflow based on the provided WorkflowConfig.

    Workflow steps:
      1. Load data from .h5ad file
      2. Compute QC metrics & filter cells
      3. Normalize data
      4. Select HVGs, compute PCA, build neighbors graph, run Leiden clustering
      5. Run analyses (PCA/UMAP/t-SNE/DE/Pathway) as configured
      6. Save processed AnnData

    Args:
        config: WorkflowConfig instance with all parameters
        save_intermediate: If True, save intermediate .h5ad files (QC-filtered, normalized)

    Returns:
        Dictionary with summary information:
        - initial_n_cells: Number of cells before processing
        - initial_n_genes: Number of genes before processing
        - final_n_cells: Number of cells after processing
        - final_n_genes: Number of genes after processing
        - qc_summary: Dictionary with QC filtering statistics
        - qc_adata_path: Path to QC-filtered file (if save_intermediate=True)
        - normalized_adata_path: Path to normalized file (if save_intermediate=True)
        - final_adata_path: Path to final processed file
        - de_results: Differential expression results (if enabled)
        - pathway_results: Pathway enrichment results (if enabled)

    Raises:
        FileNotFoundError: If input file doesn't exist
        ValueError: If configuration is invalid
        RuntimeError: If workflow step fails

    Example:
        >>> results = run_workflow(config, save_intermediate=True)
        >>> print(f"Processed {results['final_n_cells']} cells")
    """
    logger.info("=" * 60)
    logger.info("Starting RNA-seq workflow")
    logger.info(f"Input file: {config.input.input_file}")
    logger.info(f"Output directory: {config.output.output_dir}")
    logger.info("=" * 60)

    results: Dict[str, Any] = {}

    try:
        # 1) Load data
        logger.info("Step 1/6: Loading data")
        adata = load_data(config.input.input_file)
        results["initial_n_cells"] = int(adata.n_obs)
        results["initial_n_genes"] = int(adata.n_vars)

        # 2) QC
        logger.info("Step 2/6: Quality control")
        adata = compute_qc_metrics(adata)
        adata, qc_summary = filter_cells(adata, config.qc)
        results["qc_summary"] = qc_summary

        if save_intermediate:
            results["qc_adata_path"] = save_adata(
                adata, config.output.output_dir, "qc_filtered.h5ad"
            )

        # 3) Normalization
        logger.info("Step 3/6: Normalization")
        adata = normalize_data(adata, config.normalization)
        if save_intermediate:
            results["normalized_adata_path"] = save_adata(
                adata, config.output.output_dir, "normalized.h5ad"
            )

        # 4) HVGs + PCA + neighbors + Leiden
        logger.info("Step 4/6: Clustering (HVGs, PCA, neighbors, Leiden)")
        adata = select_hvgs(adata, config.clustering)
        adata = run_pca_and_neighbors(adata, config.clustering)
        adata = run_leiden_clustering(adata, config.clustering)

        # 5) Analyses
        logger.info("Step 5/6: Running analyses")
        adata = run_dimensionality_reduction(adata, config.analysis)
        de_results: Optional[Dict[str, Any]] = {}
        if config.analysis.perform_differential_expression:
            logger.info("Running differential expression analysis")
            de_results = run_differential_expression(adata, config.analysis, groupby="leiden")
        pathway_results = run_pathway_enrichment_placeholder(config.analysis)

        results["de_results"] = de_results
        results["pathway_results"] = pathway_results

        # 6) Save final AnnData
        logger.info("Step 6/6: Saving final processed data")
        results["final_adata_path"] = save_adata(adata, config.output.output_dir, "processed.h5ad")
        results["final_n_cells"] = int(adata.n_obs)
        results["final_n_genes"] = int(adata.n_vars)

        logger.info("=" * 60)
        logger.info("Workflow completed successfully!")
        logger.info(f"Final: {results['final_n_cells']} cells, {results['final_n_genes']} genes")
        logger.info("=" * 60)

    except Exception as e:
        logger.error(f"Workflow failed at step: {e}", exc_info=True)
        raise RuntimeError(f"Workflow execution failed: {e}") from e

    return results
