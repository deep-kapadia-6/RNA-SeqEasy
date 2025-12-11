from pydantic import BaseModel, Field, field_validator, ConfigDict
from typing import Optional, Literal


class TechStackConfig(BaseModel):
    """Technology stack configuration"""

    stack: Literal["Python (h5ad)"] = "Python (h5ad)"


class InputConfig(BaseModel):
    """Input file configuration"""

    input_file: str = Field(..., description="Path to input .h5ad file")
    metadata_file: Optional[str] = Field(None, description="Path to metadata CSV/TSV file")

    @field_validator("input_file")
    @classmethod
    def validate_input_file(cls, v):
        if not v.endswith(".h5ad"):
            raise ValueError("Input file must be a .h5ad file")
        return v


class OutputConfig(BaseModel):
    """Output configuration"""

    output_dir: str = Field(..., description="Output directory path")


class QCConfig(BaseModel):
    """Quality control parameters"""

    min_genes_per_cell: int = Field(200, ge=0, description="Minimum genes per cell")
    max_mitochondrial_fraction: float = Field(
        0.1, ge=0.0, le=1.0, description="Maximum mitochondrial fraction"
    )
    min_umi_count: int = Field(500, ge=0, description="Minimum UMI count per cell")


class NormalizationConfig(BaseModel):
    """Normalization method configuration"""

    method: Literal["Log normalization", "SCTransform (future)", "CPM/TPM (future)"] = (
        "Log normalization"
    )


class ClusteringConfig(BaseModel):
    """Clustering parameters"""

    n_highly_variable_genes: int = Field(
        2000, ge=100, description="Number of highly variable genes"
    )
    clustering_algorithm: Literal["Leiden"] = "Leiden"
    clustering_resolution: float = Field(0.6, ge=0.1, le=2.0, description="Clustering resolution")
    n_principal_components: int = Field(
        30, ge=10, le=100, description="Number of principal components"
    )


class AnalysisConfig(BaseModel):
    """Analysis options configuration"""

    perform_pca: bool = True
    perform_umap: bool = True
    perform_tsne: bool = False
    perform_differential_expression: bool = False
    perform_pathway_enrichment: bool = False


class WorkflowConfig(BaseModel):
    """Complete workflow configuration"""

    tech_stack: TechStackConfig = Field(default_factory=TechStackConfig)
    input: InputConfig
    output: OutputConfig
    qc: QCConfig = Field(default_factory=QCConfig)
    normalization: NormalizationConfig = Field(default_factory=NormalizationConfig)
    clustering: ClusteringConfig = Field(default_factory=ClusteringConfig)
    analysis: AnalysisConfig = Field(default_factory=AnalysisConfig)

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "tech_stack": {"stack": "Python (h5ad)"},
                "input": {"input_file": "data/example.h5ad", "metadata_file": "data/metadata.csv"},
                "output": {"output_dir": "outputs/"},
                "qc": {
                    "min_genes_per_cell": 200,
                    "max_mitochondrial_fraction": 0.1,
                    "min_umi_count": 500,
                },
                "normalization": {"method": "Log normalization"},
                "clustering": {
                    "n_highly_variable_genes": 2000,
                    "clustering_algorithm": "Leiden",
                    "clustering_resolution": 0.6,
                    "n_principal_components": 30,
                },
                "analysis": {
                    "perform_pca": True,
                    "perform_umap": True,
                    "perform_tsne": False,
                    "perform_differential_expression": False,
                    "perform_pathway_enrichment": False,
                },
            }
        }
    )
