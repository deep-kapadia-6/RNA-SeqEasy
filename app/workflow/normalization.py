# app/workflow/normalization.py
import scanpy as sc
from anndata import AnnData
from app.models import NormalizationConfig


def normalize_data(adata: AnnData, norm_config: NormalizationConfig) -> AnnData:
    """
    Apply normalization based on NormalizationConfig.

    Currently implements:
      - Log normalization (library-size normalize + log1p)

    Args:
        adata: Input AnnData.
        norm_config: NormalizationConfig instance.

    Returns:
        Normalized AnnData.
    """
    if norm_config.method == "Log normalization":
        # Normalize counts per cell to 1e4 and log-transform
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    else:
        # Placeholders for future methods
        raise NotImplementedError(
            f"Normalization method '{norm_config.method}' is not implemented yet."
        )

    return adata
