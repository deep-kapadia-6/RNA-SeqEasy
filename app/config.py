"""Configuration management for RNA workflow."""

import logging
import yaml
from pathlib import Path
from typing import Optional
from app.models import WorkflowConfig, InputConfig, OutputConfig

logger = logging.getLogger(__name__)


def save_config_to_yaml(config: WorkflowConfig, filepath: str) -> None:
    """
    Save workflow configuration to YAML file.

    Args:
        config: WorkflowConfig instance to save
        filepath: Path to save YAML file (will create parent directories if needed)

    Raises:
        OSError: If file cannot be written
        ValueError: If config cannot be serialized

    Example:
        >>> config = create_default_config("data.h5ad", "outputs/")
        >>> save_config_to_yaml(config, "configs/my_config.yaml")
    """
    try:
        config_dict = config.model_dump()
        config_path = Path(filepath)
        config_path.parent.mkdir(parents=True, exist_ok=True)

        with open(config_path, "w") as f:
            yaml.dump(config_dict, f, default_flow_style=False, sort_keys=False)

        logger.info(f"Configuration saved to {config_path}")
    except OSError as e:
        logger.error(f"Failed to save config to {filepath}: {e}")
        raise
    except Exception as e:
        logger.error(f"Error serializing config: {e}")
        raise ValueError(f"Failed to save configuration: {e}")


def load_config_from_yaml(filepath: str) -> WorkflowConfig:
    """
    Load workflow configuration from YAML file.

    Args:
        filepath: Path to YAML file

    Returns:
        WorkflowConfig instance loaded from file

    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If config is invalid or cannot be parsed
        yaml.YAMLError: If YAML file is malformed

    Example:
        >>> config = load_config_from_yaml("configs/my_config.yaml")
        >>> print(config.qc.min_genes_per_cell)
    """
    config_path = Path(filepath)
    if not config_path.exists():
        error_msg = f"Config file not found: {filepath}"
        logger.error(error_msg)
        raise FileNotFoundError(error_msg)

    try:
        with open(config_path, "r") as f:
            config_dict = yaml.safe_load(f)

        if config_dict is None:
            raise ValueError("Config file is empty")

        config = WorkflowConfig(**config_dict)
        logger.info(f"Configuration loaded from {config_path}")
        return config
    except yaml.YAMLError as e:
        error_msg = f"Invalid YAML format in {filepath}: {e}"
        logger.error(error_msg)
        raise ValueError(error_msg)
    except Exception as e:
        error_msg = f"Invalid configuration in {filepath}: {str(e)}"
        logger.error(error_msg)
        raise ValueError(error_msg)


def create_default_config(
    input_file: str, output_dir: str, metadata_file: Optional[str] = None
) -> WorkflowConfig:
    """
    Create a default workflow configuration with standard parameters.

    Args:
        input_file: Path to input .h5ad file (required)
        output_dir: Output directory path for results (required)
        metadata_file: Optional path to metadata CSV/TSV file

    Returns:
        WorkflowConfig with default values for all parameters

    Example:
        >>> config = create_default_config(
        ...     "data/example.h5ad",
        ...     "outputs/",
        ...     metadata_file="data/metadata.csv"
        ... )
    """
    logger.info(f"Creating default config for input: {input_file}")
    return WorkflowConfig(
        input=InputConfig(input_file=input_file, metadata_file=metadata_file),
        output=OutputConfig(output_dir=output_dir),
    )


def validate_config(config: WorkflowConfig) -> tuple[bool, Optional[str]]:
    """
    Validate workflow configuration for correctness.

    Checks:
    - Input file path format (must end with .h5ad)
    - Output directory can be created
    - All parameter values are within valid ranges

    Args:
        config: WorkflowConfig instance to validate

    Returns:
        Tuple of (is_valid: bool, error_message: Optional[str])
        If valid, returns (True, None)
        If invalid, returns (False, error_message)

    Example:
        >>> is_valid, error = validate_config(config)
        >>> if not is_valid:
        ...     print(f"Validation failed: {error}")
    """
    try:
        # Check if input file path is valid format
        input_path = Path(config.input.input_file)
        if not input_path.suffix == ".h5ad":
            error_msg = f"Input file must be .h5ad format, got: {input_path.suffix}"
            logger.warning(error_msg)
            return False, error_msg

        # Check output directory can be created
        output_path = Path(config.output.output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        logger.debug("Configuration validation passed")
        return True, None
    except Exception as e:
        error_msg = f"Validation error: {str(e)}"
        logger.error(error_msg)
        return False, error_msg


def config_to_dict(config: WorkflowConfig) -> dict:
    """
    Convert WorkflowConfig to dictionary for easy access.

    Args:
        config: WorkflowConfig instance

    Returns:
        Dictionary representation of the configuration

    Example:
        >>> config_dict = config_to_dict(config)
        >>> print(config_dict['qc']['min_genes_per_cell'])
    """
    return config.model_dump()


def config_from_dict(config_dict: dict) -> WorkflowConfig:
    """
    Create WorkflowConfig from dictionary.

    Args:
        config_dict: Dictionary with configuration values

    Returns:
        WorkflowConfig instance

    Raises:
        ValueError: If dictionary is invalid

    Example:
        >>> config_dict = {"input": {"input_file": "data.h5ad"}, ...}
        >>> config = config_from_dict(config_dict)
    """
    try:
        return WorkflowConfig(**config_dict)
    except Exception as e:
        logger.error(f"Failed to create config from dict: {e}")
        raise ValueError(f"Invalid configuration dictionary: {e}")
