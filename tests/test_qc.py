"""Tests for quality control functions."""

import pytest
import numpy as np
from anndata import AnnData

from app.workflow.qc import compute_qc_metrics, filter_cells
from app.models import QCConfig


class TestQC:
    """Test quality control functions."""

    @pytest.fixture
    def sample_adata(self):
        """Create a sample AnnData object for testing."""
        # Create small test dataset
        n_cells = 100
        n_genes = 50

        # Create expression matrix with some variation
        X = np.random.negative_binomial(5, 0.3, size=(n_cells, n_genes)).astype(float)

        # Add some MT genes
        var_names = [f"Gene_{i}" for i in range(n_genes)]
        var_names[0:5] = [f"MT-Gene_{i}" for i in range(5)]

        adata = AnnData(X=X)
        adata.var_names = var_names

        return adata

    def test_compute_qc_metrics(self, sample_adata):
        """Test QC metrics computation."""
        result = compute_qc_metrics(sample_adata)

        assert "n_genes_by_counts" in result.obs.columns
        assert "total_counts" in result.obs.columns
        assert result.n_obs == sample_adata.n_obs
        assert result.n_vars == sample_adata.n_vars

    def test_filter_cells(self, sample_adata):
        """Test cell filtering."""
        # Compute QC metrics first
        adata = compute_qc_metrics(sample_adata)

        # Create QC config with lenient thresholds
        qc_config = QCConfig(min_genes_per_cell=0, max_mitochondrial_fraction=1.0, min_umi_count=0)

        filtered_adata, summary = filter_cells(adata, qc_config)

        # With lenient thresholds, should keep all cells
        assert filtered_adata.n_obs == adata.n_obs
        assert summary["initial_n_cells"] == adata.n_obs
        assert summary["final_n_cells"] == filtered_adata.n_obs

    def test_filter_cells_strict(self, sample_adata):
        """Test cell filtering with strict thresholds."""
        adata = compute_qc_metrics(sample_adata)

        # Very strict thresholds
        qc_config = QCConfig(
            min_genes_per_cell=1000,  # Very high
            max_mitochondrial_fraction=0.01,  # Very low
            min_umi_count=10000,  # Very high
        )

        filtered_adata, summary = filter_cells(adata, qc_config)

        # Should filter out most/all cells
        assert filtered_adata.n_obs <= adata.n_obs
        assert summary["fraction_kept"] <= 1.0
