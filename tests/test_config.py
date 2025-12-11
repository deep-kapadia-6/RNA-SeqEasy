"""Tests for configuration management."""

import pytest  # noqa: F401 - pytest is used by test framework

from app.config import (
    create_default_config,
    save_config_to_yaml,
    load_config_from_yaml,
    validate_config,
    config_to_dict,
    config_from_dict,
)
from app.models import WorkflowConfig


class TestConfigManagement:
    """Test configuration save/load functionality."""

    def test_create_default_config(self):
        """Test creating a default configuration."""
        config = create_default_config(
            input_file="data/test.h5ad", output_dir="outputs/", metadata_file="data/metadata.csv"
        )

        assert config.input.input_file == "data/test.h5ad"
        assert config.output.output_dir == "outputs/"
        assert config.input.metadata_file == "data/metadata.csv"
        assert config.qc.min_genes_per_cell == 200  # default value

    def test_save_and_load_config(self, tmp_path):
        """Test saving and loading configuration."""
        config = create_default_config("data/test.h5ad", "outputs/")
        config_file = tmp_path / "test_config.yaml"

        # Save
        save_config_to_yaml(config, str(config_file))
        assert config_file.exists()

        # Load
        loaded_config = load_config_from_yaml(str(config_file))
        assert loaded_config.input.input_file == config.input.input_file
        assert loaded_config.output.output_dir == config.output.output_dir

    def test_validate_config(self):
        """Test configuration validation."""
        config = create_default_config("data/test.h5ad", "outputs/")
        is_valid, error = validate_config(config)

        assert is_valid is True
        assert error is None

    def test_config_to_dict(self):
        """Test converting config to dictionary."""
        config = create_default_config("data/test.h5ad", "outputs/")
        config_dict = config_to_dict(config)

        assert isinstance(config_dict, dict)
        assert "input" in config_dict
        assert "output" in config_dict
        assert config_dict["input"]["input_file"] == "data/test.h5ad"

    def test_config_from_dict(self):
        """Test creating config from dictionary."""
        config_dict = {
            "input": {"input_file": "data/test.h5ad", "metadata_file": None},
            "output": {"output_dir": "outputs/"},
        }
        config = config_from_dict(config_dict)

        assert isinstance(config, WorkflowConfig)
        assert config.input.input_file == "data/test.h5ad"
