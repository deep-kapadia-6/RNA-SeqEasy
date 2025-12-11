"""FastAPI application for RNA workflow management."""

import logging
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pathlib import Path
from typing import Optional, Dict, Any

from app.models import WorkflowConfig
from app.config import (
    validate_config,
    save_config_to_yaml,
    load_config_from_yaml,
)
from app.workflow.processor import run_workflow

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# This is what Uvicorn is looking for: a FastAPI instance named "app"
app = FastAPI(title="RNA Workflow API", version="0.1.0")

# Allow local frontends (Streamlit) to call the API
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # for local dev; tighten for production
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# In-memory "current config" for this process
CURRENT_CONFIG: Optional[WorkflowConfig] = None


@app.get("/health")
def health_check() -> Dict[str, str]:
    """
    Simple health check endpoint.

    Returns:
        Dictionary with status "ok" if API is running

    Example:
        >>> GET /health
        {"status": "ok"}
    """
    return {"status": "ok"}


@app.post("/config", response_model=WorkflowConfig)
def set_config(config: WorkflowConfig) -> WorkflowConfig:
    """
    Set the current workflow configuration (from JSON body).

    Args:
        config: WorkflowConfig instance from request body

    Returns:
        The validated WorkflowConfig that was set

    Raises:
        HTTPException: 400 if config is invalid

    Example:
        >>> POST /config
        >>> Body: {"input": {"input_file": "data.h5ad"}, ...}
    """
    global CURRENT_CONFIG
    logger.info("Setting workflow configuration")

    is_valid, error = validate_config(config)
    if not is_valid:
        error_msg = f"Invalid config: {error}"
        logger.warning(error_msg)
        raise HTTPException(status_code=400, detail=error_msg)

    CURRENT_CONFIG = config
    logger.info(f"Configuration set successfully for input: {config.input.input_file}")
    return config


@app.get("/config", response_model=WorkflowConfig)
def get_config() -> WorkflowConfig:
    """
    Get the current workflow configuration.

    Returns:
        Current WorkflowConfig instance

    Raises:
        HTTPException: 404 if no configuration is set

    Example:
        >>> GET /config
        Returns current configuration or 404 if not set
    """
    if CURRENT_CONFIG is None:
        logger.warning("Attempted to get config but none is set")
        raise HTTPException(status_code=404, detail="No configuration set")

    logger.debug("Returning current configuration")
    return CURRENT_CONFIG


@app.post("/config/save")
def save_current_config(name: str) -> Dict[str, str]:
    """
    Save the current configuration to configs/{name}.yaml.

    Args:
        name: Name for the configuration file (without .yaml extension)

    Returns:
        Dictionary with success message and file path

    Raises:
        HTTPException: 404 if no configuration is set
        HTTPException: 500 if save fails

    Example:
        >>> POST /config/save?name=my_config
        {"message": "Config saved", "path": "configs/my_config.yaml"}
    """
    if CURRENT_CONFIG is None:
        logger.warning("Attempted to save config but none is set")
        raise HTTPException(status_code=404, detail="No configuration set")

    config_path = f"configs/{name}.yaml"
    try:
        save_config_to_yaml(CURRENT_CONFIG, config_path)
        logger.info(f"Configuration saved to {config_path}")
        return {"message": "Config saved", "path": config_path}
    except Exception as e:
        error_msg = f"Failed to save config: {str(e)}"
        logger.error(error_msg)
        raise HTTPException(status_code=500, detail=error_msg)


@app.post("/config/load", response_model=WorkflowConfig)
def load_config_endpoint(name: str) -> WorkflowConfig:
    """
    Load configuration from configs/{name}.yaml and set as current.

    Args:
        name: Name of configuration file (without .yaml extension)

    Returns:
        Loaded WorkflowConfig instance

    Raises:
        HTTPException: 404 if config file not found
        HTTPException: 400 if config is invalid

    Example:
        >>> POST /config/load?name=my_config
        Loads and returns the configuration
    """
    global CURRENT_CONFIG
    logger.info(f"Loading configuration: {name}")

    config_path = Path("configs") / f"{name}.yaml"
    if not config_path.exists():
        error_msg = f"Config file not found: {config_path}"
        logger.warning(error_msg)
        raise HTTPException(status_code=404, detail=error_msg)

    try:
        config = load_config_from_yaml(str(config_path))
        is_valid, error = validate_config(config)
        if not is_valid:
            error_msg = f"Invalid config: {error}"
            logger.warning(error_msg)
            raise HTTPException(status_code=400, detail=error_msg)

        CURRENT_CONFIG = config
        logger.info(f"Configuration loaded successfully from {config_path}")
        return config
    except ValueError as e:
        error_msg = f"Invalid configuration: {str(e)}"
        logger.error(error_msg)
        raise HTTPException(status_code=400, detail=error_msg)
    except Exception as e:
        error_msg = f"Failed to load config: {str(e)}"
        logger.error(error_msg)
        raise HTTPException(status_code=500, detail=error_msg)


@app.post("/workflow/run")
def run_workflow_endpoint(save_intermediate: bool = True) -> Dict[str, Any]:
    """
    Run the full workflow using the current configuration.

    Returns a summary dict with:
      - initial/final cell & gene counts
      - qc_summary
      - paths to intermediate/final .h5ad files
      - analysis results (DE, pathway enrichment)

    Args:
        save_intermediate: If True, save intermediate .h5ad files (default: True)

    Returns:
        Dictionary with workflow results and summary statistics

    Raises:
        HTTPException: 404 if no configuration is set
        HTTPException: 500 if workflow execution fails

    Example:
        >>> POST /workflow/run?save_intermediate=true
        Returns workflow results dictionary
    """
    if CURRENT_CONFIG is None:
        logger.warning("Attempted to run workflow but no config is set")
        raise HTTPException(status_code=404, detail="No configuration set")

    logger.info("Starting workflow execution")
    try:
        results = run_workflow(CURRENT_CONFIG, save_intermediate=save_intermediate)
        logger.info("Workflow completed successfully")
        return results
    except FileNotFoundError as e:
        error_msg = f"Input file not found: {str(e)}"
        logger.error(error_msg)
        raise HTTPException(status_code=404, detail=error_msg)
    except ValueError as e:
        error_msg = f"Configuration error: {str(e)}"
        logger.error(error_msg)
        raise HTTPException(status_code=400, detail=error_msg)
    except Exception as e:
        error_msg = f"Workflow execution failed: {str(e)}"
        logger.error(error_msg, exc_info=True)
        raise HTTPException(status_code=500, detail=error_msg)
