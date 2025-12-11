"""Quality control functions for RNA sequencing data."""

import logging
import scanpy as sc
from anndata import AnnData
from typing import Tuple
from app.models import QCConfig

logger = logging.getLogger(__name__)


def compute_qc_metrics(adata: AnnData, mt_prefix: str = "MT-") -> AnnData:
    """
    Compute basic QC metrics and store them in adata.obs.

    Calculates:
    - n_genes_by_counts: Number of genes detected per cell
    - total_counts: Total UMI counts per cell
    - pct_counts_mt: Percentage of mitochondrial gene counts (if MT genes found)

    Args:
        adata: Input AnnData object with gene expression data
        mt_prefix: Prefix used to identify mitochondrial genes (default: "MT-")

    Returns:
        AnnData object with QC metrics added to .obs columns

    Raises:
        ValueError: If adata is empty or invalid

    Example:
        >>> adata = compute_qc_metrics(adata)
        >>> print(adata.obs['n_genes_by_counts'].head())
    """
    if adata.n_obs == 0 or adata.n_vars == 0:
        error_msg = "AnnData object is empty"
        logger.error(error_msg)
        raise ValueError(error_msg)

    logger.info(f"Computing QC metrics for {adata.n_obs} cells and {adata.n_vars} genes")

    # Flag mitochondrial genes if gene names start with mt_prefix
    if "mt" not in adata.var:
        adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
        mt_count = adata.var["mt"].sum()
        logger.debug(f"Found {mt_count} mitochondrial genes with prefix '{mt_prefix}'")

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )

    logger.info("QC metrics computed successfully")
    return adata


def filter_cells(
    adata: AnnData,
    qc_config: QCConfig,
) -> Tuple[AnnData, dict]:
    """
    Filter cells based on QCConfig thresholds.

    Applies three filters:
    1. Minimum genes per cell
    2. Maximum mitochondrial fraction
    3. Minimum UMI count per cell

    Args:
        adata: AnnData with QC metrics already computed (from compute_qc_metrics)
        qc_config: QCConfig instance with filtering thresholds

    Returns:
        Tuple of (filtered_adata, qc_summary dict)
        qc_summary contains:
        - initial_n_cells: Number of cells before filtering
        - final_n_cells: Number of cells after filtering
        - fraction_kept: Fraction of cells retained

    Raises:
        ValueError: If required QC metrics are missing

    Example:
        >>> filtered_adata, summary = filter_cells(adata, qc_config)
        >>> print(f"Kept {summary['final_n_cells']} of {summary['initial_n_cells']} cells")
    """
    initial_n_cells = adata.n_obs

    # Check required metrics exist
    required_metrics = ["n_genes_by_counts"]
    missing = [m for m in required_metrics if m not in adata.obs.columns]
    if missing:
        error_msg = f"Missing required QC metrics: {missing}. Run compute_qc_metrics first."
        logger.error(error_msg)
        raise ValueError(error_msg)

    logger.info(
        f"Filtering cells with thresholds: "
        f"min_genes={qc_config.min_genes_per_cell}, "
        f"max_mito={qc_config.max_mitochondrial_fraction}, "
        f"min_umi={qc_config.min_umi_count}"
    )

    # Boolean masks for each filter
    mask_min_genes = adata.obs["n_genes_by_counts"] >= qc_config.min_genes_per_cell
    mask_max_mito = (
        adata.obs["pct_counts_mt"] <= qc_config.max_mitochondrial_fraction * 100.0
        if "pct_counts_mt" in adata.obs
        else True
    )
    mask_min_umi = (
        adata.obs["total_counts"] >= qc_config.min_umi_count
        if "total_counts" in adata.obs
        else True
    )

    combined_mask = mask_min_genes & mask_max_mito & mask_min_umi
    filtered_adata = adata[combined_mask].copy()

    cells_removed = initial_n_cells - filtered_adata.n_obs
    fraction_kept = float(filtered_adata.n_obs / initial_n_cells) if initial_n_cells > 0 else 0.0

    logger.info(
        f"Filtered {cells_removed} cells ({initial_n_cells - filtered_adata.n_obs} removed, "
        f"{fraction_kept:.2%} kept)"
    )

    qc_summary = {
        "initial_n_cells": int(initial_n_cells),
        "final_n_cells": int(filtered_adata.n_obs),
        "fraction_kept": fraction_kept,
    }

    return filtered_adata, qc_summary
