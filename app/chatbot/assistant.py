from typing import Optional
from app.chatbot.llm_interface import get_llm_provider, LLMProvider
from app.models import WorkflowConfig


class WorkflowAssistant:
    """Chatbot assistant for RNA workflow questions."""

    def __init__(self, provider: Optional[LLMProvider] = None):
        """
        Initialize the assistant.

        Args:
            provider: LLMProvider instance. If None, uses default from env.
        """
        self.provider = provider or get_llm_provider()
        self.system_prompt = self._build_system_prompt()

    def _build_system_prompt(self) -> str:
        """Build the system prompt with workflow context."""
        return """You are a helpful bioinformatics assistant specializing in RNA sequencing analysis workflows.

You help users understand:
- Quality control parameters (min genes per cell, mitochondrial fraction, UMI counts)
- Normalization methods (Log normalization, SCTransform, CPM/TPM)
- Clustering parameters (highly variable genes, resolution, principal components)
- Analysis options (PCA, UMAP, t-SNE, differential expression, pathway enrichment)
- Workflow steps and best practices
- Troubleshooting common errors

Provide clear, concise explanations. Use technical terms appropriately but explain them when needed.
If asked about specific parameter values, provide typical ranges and explain their impact."""

    def get_response(self, user_message: str, config: Optional[WorkflowConfig] = None) -> str:
        """
        Get assistant response to user message.

        Args:
            user_message: User's question
            config: Optional current workflow config for context

        Returns:
            Assistant's response
        """
        # Enhance prompt with config context if available
        enhanced_prompt = user_message
        if config:
            config_context = self._format_config_context(config)
            enhanced_prompt = f"{user_message}\n\nCurrent workflow configuration:\n{config_context}"

        try:
            response = self.provider.generate_response(
                prompt=enhanced_prompt, system_prompt=self.system_prompt
            )
            return response
        except Exception as e:
            return f"I encountered an error: {str(e)}. Please check your API key and provider settings."

    def _format_config_context(self, config: WorkflowConfig) -> str:
        """Format workflow config as context string."""
        context_parts = []

        if config.qc:
            context_parts.append(
                f"QC: min_genes={config.qc.min_genes_per_cell}, "
                f"max_mito={config.qc.max_mitochondrial_fraction}, "
                f"min_umi={config.qc.min_umi_count}"
            )

        if config.normalization:
            context_parts.append(f"Normalization: {config.normalization.method}")

        if config.clustering:
            context_parts.append(
                f"Clustering: n_hvg={config.clustering.n_highly_variable_genes}, "
                f"resolution={config.clustering.clustering_resolution}, "
                f"n_pcs={config.clustering.n_principal_components}"
            )

        if config.analysis:
            analyses = [k for k, v in config.analysis.model_dump().items() if v]
            context_parts.append(f"Analyses: {', '.join(analyses) if analyses else 'none'}")

        return "\n".join(context_parts)
