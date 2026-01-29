import requests
import streamlit as st
from pathlib import Path
from app.config import (
    create_default_config,
    save_config_to_yaml,
    load_config_from_yaml,
    config_to_dict,
)

st.set_page_config(
    page_title="RNA-SeqEasy",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

if "config" not in st.session_state:
    st.session_state.config = None
if "uploaded_files" not in st.session_state:
    st.session_state.uploaded_files = {}
if "workflow_status" not in st.session_state:
    st.session_state.workflow_status = "not_started"

st.title("🧬 RNA Workflow Manager")
st.subheader("By Deep Kapadia")

with st.sidebar:
    st.title("Navigation")
    st.markdown("---")

    page = st.radio(
        "Navigation",
        [
            "📁 Input & Config",
            "🔬 Quality Control",
            "📊 Normalization",
            "🔍 Clustering",
            "📈 Analysis",
            "💾 Save/Load Config",
            "🚀 Run Workflow",
        ],
        label_visibility="collapsed",
    )


# Add config status display
st.markdown("---")
if st.session_state.config:
    st.success("✅ Config Ready")
    st.caption(f"Input: {Path(st.session_state.config.input.input_file).name}")
else:
    st.warning("⚠️ No Config")
    st.caption("Create config in Input & Config page")


def render_input_page():
    st.header("📁 Input Files & Basic Configuration")

    # Help section
    with st.expander("📚 Learn More About Input Files"):
        st.markdown(
            """
        ### Input File Formats
        
        **.h5ad Files:**
        - H5AD (HDF5 AnnData) is the native format for Scanpy
        - Contains gene expression matrix, cell metadata, and gene annotations
        - Created from raw count matrices using tools like Scanpy or Seurat
        
        **Metadata Files:**
        - Optional CSV/TSV file with additional cell-level information
        - Common columns: cell_type, batch, condition, treatment
        - Must have same number of rows as cells in your .h5ad file
        - First column should match cell barcodes/IDs in the .h5ad file
        
        ### Best Practices
        - Ensure your .h5ad file is properly formatted (not corrupted)
        - Verify metadata file matches cell order in .h5ad
        - Use descriptive output directory names for different experiments
        """
        )

    # Tech Stack Selection
    st.subheader("Technology Stack")
    st.selectbox(
        "Choose tech stack",
        ["Python (h5ad)"],
        help="Currently only Python (h5ad) format is supported. H5AD is the native AnnData format used by Scanpy for efficient storage of single-cell data.",
    )

    # File Upload Section
    st.subheader("Upload Files")

    uploaded_h5ad = st.file_uploader(
        "Upload .h5ad file",
        type=["h5ad"],
        help="""Upload your RNA sequencing data in h5ad format.
        
        Requirements:
        - File must be a valid .h5ad (HDF5 AnnData) file
        - Should contain raw or normalized count matrix
        - Typically created using Scanpy or converted from other formats
        
        File size: Can handle large datasets (millions of cells)""",
    )

    uploaded_metadata = st.file_uploader(
        "Upload metadata file (CSV/TSV)",
        type=["csv", "tsv"],
        help="""Optional: Upload cell metadata file.
        
        Format:
        - CSV or TSV format
        - First column should contain cell barcodes/IDs
        - Additional columns can include: cell_type, batch, condition, etc.
        
        Note: Row count must match number of cells in .h5ad file""",
    )

    # Output Directory
    st.subheader("Output Configuration")
    output_dir = st.text_input(
        "Output directory",
        value="outputs/",
        help="""Directory where analysis results will be saved.
        
        Default: outputs/
        - All intermediate and final .h5ad files will be saved here
        - Directory will be created if it doesn't exist
        - Use descriptive names for different experiments (e.g., outputs/experiment1/)""",
    )

    # Create config button
    if st.button("Create Configuration", type="primary"):
        if uploaded_h5ad is None:
            st.error("Please upload a .h5ad file")
        else:
            # Save uploaded files temporarily
            temp_dir = Path("temp_uploads")
            temp_dir.mkdir(exist_ok=True)

            h5ad_path = temp_dir / uploaded_h5ad.name
            with open(h5ad_path, "wb") as f:
                f.write(uploaded_h5ad.getbuffer())

            metadata_path = None
            if uploaded_metadata:
                metadata_path = temp_dir / uploaded_metadata.name
                with open(metadata_path, "wb") as f:
                    f.write(uploaded_metadata.getbuffer())

            # Create config
            config = create_default_config(
                input_file=str(h5ad_path),
                output_dir=output_dir,
                metadata_file=str(metadata_path) if metadata_path else None,
            )

            st.session_state.config = config
            st.session_state.uploaded_files = {
                "h5ad": str(h5ad_path),
                "metadata": str(metadata_path) if metadata_path else None,
            }

            st.success("Configuration created successfully!")


def render_qc_page():
    st.header("🔬 Quality Control Parameters")

    if st.session_state.config is None:
        st.warning("Please create a configuration first in the Input & Config page")
        return

    config = st.session_state.config

    # Help section
    with st.expander("📚 Learn More About Quality Control"):
        st.markdown(
            """
        ### What is Quality Control?
        
        Quality control (QC) is the first step in RNA-seq analysis to identify and remove low-quality cells.
        Poor quality cells can arise from:
        - Cell death or stress during sample preparation
        - Empty droplets (no cells captured)
        - Doublets (multiple cells in one droplet)
        - Technical artifacts
        
        ### Key Metrics Explained
        
        **1. Genes per Cell (n_genes_by_counts)**
        - Number of unique genes detected in each cell
        - Low values (< 200) may indicate:
          - Empty droplets
          - Dead or dying cells
          - Low-quality sequencing
        - Typical range: 200-5000 for single-cell RNA-seq
        - Default: 200 (good starting point)
        
        **2. Mitochondrial Fraction (pct_counts_mt)**
        - Percentage of reads mapping to mitochondrial genes
        - High values (> 10-20%) indicate:
          - Dead or stressed cells
          - Cell damage during processing
        - Typical threshold: 5-10% for healthy cells
        - Default: 0.1 (10%)
        
        **3. UMI Count (total_counts)**
        - Total unique molecular identifier counts per cell
        - Low values indicate:
          - Insufficient sequencing depth
          - Empty droplets
        - Typical range: 500-50,000 for single-cell
        - Default: 500
        
        ### Best Practices
        
        1. **Start with defaults** and adjust based on your data
        2. **Visualize QC metrics** before filtering (use Scanpy plotting)
        3. **Consider your cell type**: Some cell types naturally have lower gene counts
        4. **Check distributions**: Look for bimodal distributions indicating distinct cell populations
        5. **Be conservative**: It's better to keep questionable cells than remove valid ones
        
        ### Common Pitfalls
        
        - Setting thresholds too high → removes valid cells
        - Setting thresholds too low → keeps low-quality cells
        - Not considering cell type-specific characteristics
        """
        )

    st.info(
        """
    **Quality Control Parameters:**
    - **Min genes per cell**: Filters out cells with fewer than this many genes
    - **Max mitochondrial fraction**: Removes cells with high mitochondrial gene expression (indicates dead cells)
    - **Min UMI count**: Filters cells with low unique molecular identifier counts
    """
    )

    col1, col2, col3 = st.columns(3)

    with col1:
        min_genes = st.number_input(
            "Min genes per cell",
            min_value=0,
            value=config.qc.min_genes_per_cell,
            help="""Minimum number of genes detected per cell.
            
            Typical range: 200-1000 for single-cell RNA-seq
            - Too low (< 100): May include empty droplets or low-quality cells
            - Too high (> 2000): May exclude valid cells with low gene detection
            - Default: 200 (good starting point for most datasets)
            
            Tip: Check your data distribution first to set appropriate threshold""",
        )

    with col2:
        max_mito = st.number_input(
            "Max mitochondrial fraction",
            min_value=0.0,
            max_value=1.0,
            value=float(config.qc.max_mitochondrial_fraction),
            step=0.01,
            help="""Maximum fraction of mitochondrial genes allowed per cell.
            
            Typical range: 0.05-0.2 (5-20%)
            - Low threshold (< 0.05): Very strict, may remove stressed cells
            - High threshold (> 0.2): May keep dead/dying cells
            - Default: 0.1 (10%)
            
            Note: Some cell types (e.g., muscle cells) naturally have higher mitochondrial content""",
        )

    with col3:
        min_umi = st.number_input(
            "Min UMI count per cell",
            min_value=0,
            value=config.qc.min_umi_count,
            help="""Minimum UMI (Unique Molecular Identifier) count per cell.
            
            Typical range: 500-10,000 for single-cell RNA-seq
            - Too low (< 200): May include empty droplets
            - Too high (> 20,000): May exclude valid cells with low sequencing depth
            - Default: 500
            
            UMI counts reflect sequencing depth and cell quality""",
        )

    if st.button("Update QC Parameters"):
        config.qc.min_genes_per_cell = min_genes
        config.qc.max_mitochondrial_fraction = max_mito
        config.qc.min_umi_count = min_umi
        st.session_state.config = config
        st.success("QC parameters updated!")


def render_normalization_page():
    st.header("📊 Normalization Method")

    if st.session_state.config is None:
        st.warning("Please create a configuration first in the Input & Config page")
        return

    config = st.session_state.config

    # Help section
    with st.expander("📚 Learn More About Normalization"):
        st.markdown(
            """
        ### What is Normalization?
        
        Normalization corrects for technical variation in sequencing depth between cells,
        allowing fair comparison of gene expression across cells.
        
        ### Available Methods
        
        **1. Log Normalization (Current Implementation)**
        - Standard log(x+1) transformation
        - Scales data to 10,000 reads per cell, then log-transforms
        - Formula: `log(1 + (counts / total_counts * 10000))`
        - **When to use**: Most single-cell RNA-seq datasets
        - **Pros**: Simple, fast, widely used
        - **Cons**: May not handle highly variable sequencing depths well
        
        **2. SCTransform (Future)**
        - Regularized negative binomial regression
        - Better handles technical variation
        - **When to use**: Datasets with high technical variation
        - **Pros**: More sophisticated, better for batch correction
        - **Cons**: Computationally intensive
        
        **3. CPM/TPM (Future)**
        - Counts/Transcripts Per Million
        - Simple scaling normalization
        - **When to use**: Bulk RNA-seq or when simple scaling is sufficient
        - **Pros**: Very simple
        - **Cons**: Less sophisticated than other methods
        
        ### Best Practices
        
        - **Start with log normalization** for most cases
        - Normalize after QC filtering
        - Consider your sequencing technology (10x, Smart-seq2, etc.)
        """
        )

    st.info(
        """
    **Normalization Methods:**
    - **Log normalization**: Standard log(x+1) transformation (recommended for most cases)
    - **SCTransform**: Advanced normalization method (future implementation)
    - **CPM/TPM**: Counts per million / Transcripts per million (future implementation)
    """
    )

    norm_method = st.selectbox(
        "Normalization method",
        ["Log normalization", "SCTransform (future)", "CPM/TPM (future)"],
        index=0 if config.normalization.method == "Log normalization" else 1,
        help="""Select the normalization method for your data.
        
        **Log normalization** (Recommended):
        - Standard method for single-cell RNA-seq
        - Scales to 10,000 reads per cell, then log-transforms
        - Fast and effective for most datasets
        
        **SCTransform** (Future):
        - Advanced method using regularized negative binomial regression
        - Better for datasets with high technical variation
        
        **CPM/TPM** (Future):
        - Simple counts/transcripts per million scaling
        - Less sophisticated but faster""",
    )

    if st.button("Update Normalization"):
        config.normalization.method = norm_method
        st.session_state.config = config
        st.success("Normalization method updated!")


def render_clustering_page():
    st.header("🔍 Clustering Parameters")

    if st.session_state.config is None:
        st.warning("Please create a configuration first in the Input & Config page")
        return

    config = st.session_state.config

    # Help section
    with st.expander("📚 Learn More About Clustering"):
        st.markdown(
            """
        ### Clustering Overview
        
        Clustering identifies groups of cells with similar gene expression patterns,
        representing distinct cell types or states.
        
        ### Key Parameters
        
        **1. Highly Variable Genes (HVGs)**
        - Number of most variable genes used for analysis
        - Only these genes are used for PCA and clustering
        - Typical range: 1000-3000
        - Default: 2000
        - **Too few**: May miss important cell type markers
        - **Too many**: Includes noise, slows computation
        
        **2. Principal Components (PCs)**
        - Number of principal components for dimensionality reduction
        - Used for nearest neighbor graph and clustering
        - Typical range: 20-50
        - Default: 30
        - **Too few**: May miss important variation
        - **Too many**: Includes noise, overfitting
        
        **3. Clustering Resolution**
        - Controls granularity of clusters
        - Higher = more clusters (finer granularity)
        - Lower = fewer clusters (coarser granularity)
        - Typical range: 0.1-2.0
        - Default: 0.6
        - **0.1-0.4**: Broad cell types
        - **0.5-0.8**: Standard resolution
        - **0.8-2.0**: Fine subpopulations
        
        **4. Leiden Algorithm**
        - Graph-based clustering algorithm
        - Improvement over Louvain algorithm
        - Better at finding communities in graphs
        
        ### Workflow
        
        1. Select highly variable genes
        2. Compute PCA on HVGs
        3. Build nearest neighbor graph
        4. Apply Leiden clustering
        5. Visualize with UMAP/t-SNE
        
        ### Best Practices
        
        - Start with default values
        - Adjust resolution based on biological question
        - Use more PCs for larger datasets
        - Visualize clusters to validate results
        """
        )

    st.info(
        """
    **Clustering Parameters:**
    - **Highly Variable Genes**: Number of most variable genes to use for analysis
    - **Clustering Resolution**: Higher values create more clusters (0.1-2.0)
    - **Principal Components**: Number of PCs for dimensionality reduction
    """
    )

    col1, col2 = st.columns(2)

    with col1:
        n_hvg = st.number_input(
            "Number of highly variable genes",
            min_value=100,
            max_value=5000,
            value=config.clustering.n_highly_variable_genes,
            step=100,
            help="""Number of highly variable genes (HVGs) to select for analysis.
            
            Typical range: 1000-3000
            - Too few (< 1000): May miss important cell type markers
            - Too many (> 4000): Includes noise, slows computation
            - Default: 2000 (good balance)
            
            Only these genes are used for PCA and clustering""",
        )

        n_pcs = st.number_input(
            "Number of principal components",
            min_value=10,
            max_value=100,
            value=config.clustering.n_principal_components,
            help="""Number of principal components (PCs) for dimensionality reduction.
            
            Typical range: 20-50
            - Too few (< 20): May miss important biological variation
            - Too many (> 50): Includes noise, risk of overfitting
            - Default: 30 (good starting point)
            
            Used for nearest neighbor graph and clustering""",
        )

    with col2:
        st.selectbox(
            "Clustering algorithm",
            ["Leiden"],
            help="""Clustering algorithm to use.
            
            **Leiden Algorithm:**
            - Graph-based community detection
            - Improvement over Louvain algorithm
            - Better at finding well-connected communities
            - Currently the only supported algorithm""",
        )

        resolution = st.number_input(
            "Clustering resolution",
            min_value=0.1,
            max_value=2.0,
            value=float(config.clustering.clustering_resolution),
            step=0.1,
            help="""Clustering resolution parameter.
            
            Controls granularity of clusters:
            - **0.1-0.4**: Broad cell types (fewer clusters)
            - **0.5-0.8**: Standard resolution (default: 0.6)
            - **0.8-2.0**: Fine subpopulations (more clusters)
            
            Higher values = more clusters
            Adjust based on your biological question""",
        )

    if st.button("Update Clustering Parameters"):
        config.clustering.n_highly_variable_genes = n_hvg
        config.clustering.clustering_resolution = resolution
        config.clustering.n_principal_components = n_pcs
        st.session_state.config = config
        st.success("Clustering parameters updated!")


def render_analysis_page():
    st.header("📈 Analysis Options")

    if st.session_state.config is None:
        st.warning("Please create a configuration first in the Input & Config page")
        return

    config = st.session_state.config

    # Help section
    with st.expander("📚 Learn More About Analysis Methods"):
        st.markdown(
            """
        ### Available Analysis Methods
        
        **1. Principal Component Analysis (PCA)**
        - Linear dimensionality reduction
        - Reduces data to principal components capturing most variation
        - **When to use**: Always recommended as first step
        - **Output**: PC coordinates for each cell
        - **Interpretation**: First few PCs capture most biological variation
        
        **2. UMAP (Uniform Manifold Approximation and Projection)**
        - Non-linear dimensionality reduction
        - Preserves local structure, good for visualization
        - **When to use**: For visualization and exploration
        - **Output**: 2D/3D embedding coordinates
        - **Interpretation**: Similar cells cluster together
        
        **3. t-SNE (t-distributed Stochastic Neighbor Embedding)**
        - Alternative non-linear dimensionality reduction
        - Good for visualization, preserves local structure
        - **When to use**: Alternative to UMAP
        - **Output**: 2D/3D embedding coordinates
        - **Note**: Slower than UMAP, less global structure preservation
        
        **4. Differential Expression (DE)**
        - Identifies genes that differ between clusters
        - Statistical tests (e.g., Wilcoxon rank-sum)
        - **When to use**: To find marker genes for clusters
        - **Output**: List of differentially expressed genes per cluster
        - **Interpretation**: Genes with high fold-change and low p-value
        
        **5. Pathway Enrichment (Future)**
        - Analyzes enriched biological pathways
        - Uses gene set enrichment analysis (GSEA)
        - **When to use**: To understand biological functions
        - **Output**: Enriched pathways per cluster
        
        ### Best Practices
        
        - **Always run PCA**: Required for downstream analysis
        - **Use UMAP for visualization**: Better than t-SNE in most cases
        - **Run DE after clustering**: To identify cluster markers
        - **Start with basic analyses**: Add advanced analyses as needed
        """
        )

    st.info(
        """
    **Available Analyses:**
    - **PCA**: Principal Component Analysis for dimensionality reduction
    - **UMAP**: Uniform Manifold Approximation and Projection for visualization
    - **t-SNE**: t-distributed Stochastic Neighbor Embedding (alternative to UMAP)
    - **Differential Expression**: Find genes that differ between clusters
    - **Pathway Enrichment**: Analyze enriched biological pathways (future)
    """
    )

    col1, col2 = st.columns(2)

    with col1:
        perform_pca = st.checkbox(
            "Perform PCA",
            value=config.analysis.perform_pca,
            help="""Principal Component Analysis (PCA).
            
            Linear dimensionality reduction method.
            - Always recommended as first step
            - Required for downstream analysis
            - Reduces data to principal components
            - Default: Enabled""",
        )

        perform_umap = st.checkbox(
            "Perform UMAP",
            value=config.analysis.perform_umap,
            help="""UMAP (Uniform Manifold Approximation and Projection).
            
            Non-linear dimensionality reduction for visualization.
            - Preserves local structure
            - Good for exploring cell relationships
            - Faster than t-SNE
            - Default: Enabled""",
        )

        perform_tsne = st.checkbox(
            "Perform t-SNE",
            value=config.analysis.perform_tsne,
            help="""t-SNE (t-distributed Stochastic Neighbor Embedding).
            
            Alternative non-linear dimensionality reduction.
            - Good for visualization
            - Preserves local structure
            - Slower than UMAP
            - Default: Disabled""",
        )

    with col2:
        perform_de = st.checkbox(
            "Differential Expression",
            value=config.analysis.perform_differential_expression,
            help="""Differential Expression Analysis.
            
            Identifies genes that differ between clusters.
            - Statistical tests (Wilcoxon rank-sum)
            - Finds marker genes for each cluster
            - Useful for cluster annotation
            - Default: Disabled""",
        )

        perform_pathway = st.checkbox(
            "Pathway Enrichment",
            value=config.analysis.perform_pathway_enrichment,
            help="""Pathway Enrichment Analysis (Future Implementation).
            
            Analyzes enriched biological pathways.
            - Uses gene set enrichment analysis (GSEA)
            - Understands biological functions
            - Currently not implemented
            - Default: Disabled""",
        )

    if st.button("Update Analysis Options"):
        config.analysis.perform_pca = perform_pca
        config.analysis.perform_umap = perform_umap
        config.analysis.perform_tsne = perform_tsne
        config.analysis.perform_differential_expression = perform_de
        config.analysis.perform_pathway_enrichment = perform_pathway
        st.session_state.config = config
        st.success("Analysis options updated!")


def render_config_page():
    st.header("💾 Save & Load Configuration")

    if st.session_state.config is None:
        st.warning("Please create a configuration first in the Input & Config page")
        return

    config = st.session_state.config

    # Help section
    with st.expander("📚 Learn More About Configuration Management"):
        st.markdown(
            """
        ### Configuration Files
        
        Configurations are saved as YAML files in the `configs/` directory.
        This allows you to:
        - **Reproduce analyses**: Use the same parameters across runs
        - **Share settings**: Share configs with collaborators
        - **Version control**: Track parameter changes over time
        - **Template creation**: Create templates for common workflows
        
        ### File Format
        
        Configs are stored in YAML format, which is:
        - Human-readable
        - Easy to edit manually
        - Version control friendly
        
        ### Best Practices
        
        - Use descriptive names (e.g., `experiment1_strict_qc.yaml`)
        - Save configs before running workflows
        - Document any manual edits in comments
        - Keep configs organized in the `configs/` directory
        """
        )

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Save Configuration")
        config_name = st.text_input(
            "Config name",
            value="my_config",
            help="""Name for your configuration file.
            
            - Will be saved as: configs/{name}.yaml
            - Use descriptive names (e.g., experiment1_qc200)
            - Avoid special characters
            - File extension (.yaml) is added automatically""",
        )

        if st.button("Save Config", type="primary"):
            config_path = f"configs/{config_name}.yaml"
            save_config_to_yaml(config, config_path)
            st.success(f"Configuration saved to {config_path}!")

    with col2:
        st.subheader("Load Configuration")

        # Load example configs
        configs_dir = Path("configs")
        if configs_dir.exists():
            example_files = list(configs_dir.glob("example_*.yaml"))
            if example_files:
                st.caption("Example Configurations:")
                for ex_file in example_files:
                    if st.button(f"📋 Load {ex_file.stem}", key=f"load_ex_{ex_file.stem}"):
                        try:
                            loaded_config = load_config_from_yaml(str(ex_file))
                            st.session_state.config = loaded_config
                            st.success(f"Example configuration '{ex_file.stem}' loaded!")
                            st.rerun()
                        except Exception as e:
                            st.error(f"Error loading example: {str(e)}")

        # List existing configs
        if configs_dir.exists():
            config_files = list(configs_dir.glob("*.yaml"))
            # Filter out example files
            config_files = [f for f in config_files if not f.name.startswith("example_")]
            if config_files:
                config_options = [f.name for f in config_files]
                selected_config = st.selectbox("Select config to load", config_options)

                if st.button("Load Config"):
                    config_path = configs_dir / selected_config
                    try:
                        loaded_config = load_config_from_yaml(str(config_path))
                        st.session_state.config = loaded_config
                        st.success(f"Configuration loaded from {selected_config}!")
                    except Exception as e:
                        st.error(f"Error loading config: {str(e)}")
            else:
                st.info("No saved configurations found")
        else:
            st.info("Configs directory does not exist")

    # Config Preview
    st.subheader("Current Configuration Preview")
    with st.expander("View Config Details"):
        config_dict = config_to_dict(config)
        st.json(config_dict)


def render_chatbot():
    """Render chatbot interface at the bottom of every page."""

    # Initialize chatbot state
    if "chatbot_messages" not in st.session_state:
        st.session_state.chatbot_messages = []

    # Initialize assistant
    if "assistant" not in st.session_state:
        try:
            from app.chatbot.assistant import WorkflowAssistant

            st.session_state.assistant = WorkflowAssistant()
        except Exception as e:
            st.session_state.assistant = None
            st.session_state.chatbot_error = str(e)

    # Add separator
    st.markdown("---")

    # Chatbot section
    st.header("🤖 AI Assistant")

    if st.session_state.assistant is None:
        st.error("Chatbot not initialized. Check your API key in `.env` file.")
        if "chatbot_error" in st.session_state:
            with st.expander("Error Details"):
                st.code(st.session_state.chatbot_error)
        st.info(
            """
        **To enable chatbot:**
        1. Create a `.env` file in project root
        2. Add: `OPENAI_API_KEY=your_key_here`
        3. Add: `DEFAULT_LLM_PROVIDER=openai`
        4. Restart Streamlit
        """
        )
    else:
        st.caption("Ask questions about workflow parameters, steps, or troubleshooting")

        # Only show messages container if there are messages
        if st.session_state.chatbot_messages:
            # Display chat history in a collapsible container
            with st.expander("💬 Chat History", expanded=True):
                messages_container = st.container(height=200)
                with messages_container:
                    for message in st.session_state.chatbot_messages:
                        with st.chat_message(message["role"]):
                            st.markdown(message["content"])

        # Chat input (will appear at bottom of page, but messages are visible above)
        if prompt := st.chat_input("Ask a question about the workflow..."):
            # Add user message
            st.session_state.chatbot_messages.append({"role": "user", "content": prompt})

            # Get assistant response
            with st.chat_message("assistant"):
                with st.spinner("Thinking..."):
                    config = st.session_state.get("config")
                    response = st.session_state.assistant.get_response(prompt, config=config)
                    st.markdown(response)

            # Add assistant response to history
            st.session_state.chatbot_messages.append({"role": "assistant", "content": response})
            st.rerun()


def render_run_workflow_page():
    st.header("🚀 Run Workflow")

    # Help section
    with st.expander("📚 Learn More About Running Workflows"):
        st.markdown(
            """
        ### Workflow Execution
        
        This page triggers the full RNA-seq analysis pipeline on the FastAPI backend.
        
        ### Prerequisites
        
        1. **FastAPI Backend Running**: 
           - Start with: `uvicorn app.main:app --reload`
           - Should be running on http://127.0.0.1:8000
        
        2. **Configuration Loaded**:
           - Create config in "Input & Config" page
           - Or load existing config from "Save/Load Config" page
        
        3. **Backend Config Set**:
           - The backend needs the config loaded via API
           - Or it uses the config from Streamlit session
        
        ### Workflow Steps
        
        The pipeline executes:
        1. Load data (.h5ad file)
        2. Quality control filtering
        3. Normalization
        4. Highly variable gene selection
        5. PCA computation
        6. Nearest neighbor graph
        7. Leiden clustering
        8. Dimensionality reduction (UMAP/t-SNE if enabled)
        9. Differential expression (if enabled)
        
        ### Intermediate Files
        
        If "Save intermediate files" is enabled:
        - QC-filtered data: `{output_dir}/qc_filtered.h5ad`
        - Normalized data: `{output_dir}/normalized.h5ad`
        - Final processed data: `{output_dir}/final_processed.h5ad`
        
        ### Troubleshooting
        
        - **Connection Error**: Check if FastAPI backend is running
        - **Timeout**: Large datasets may take 10+ minutes
        - **Memory Error**: Reduce dataset size or increase system RAM
        """
        )

    st.info(
        """
    This will call the FastAPI backend to run the full Scanpy pipeline using
    the currently loaded configuration (`CURRENT_CONFIG` on the backend).
    
    Make sure:
    - The FastAPI server is running: `uvicorn app.main:app`
    - You have loaded a config via the API (e.g. `/config/load?name=my_config`)
    """
    )

    backend_url = st.text_input(
        "Backend URL",
        value="http://127.0.0.1:8000",
        help="""FastAPI backend base URL.
        
        Default: http://127.0.0.1:8000
        - Change if backend is running on different port
        - Must include protocol (http:// or https://)
        - No trailing slash""",
    )

    save_intermediate = st.checkbox(
        "Save intermediate .h5ad files (QC-filtered, normalized)",
        value=True,
        help="""Save intermediate processing files.
        
        If enabled, saves:
        - QC-filtered data (after filtering)
        - Normalized data (after normalization)
        - Final processed data (after clustering)
        
        Useful for:
        - Debugging workflow steps
        - Re-running analyses from intermediate points
        - Inspecting data at each stage""",
    )

    if st.button("Run Workflow", type="primary"):
        if not backend_url:
            st.error("Please specify the backend URL.")
            return

        # Check if config exists in Streamlit session
        if st.session_state.config is None:
            st.error(
                "No configuration found. Please create a configuration in the 'Input & Config' page first."
            )
            return

        # First, send the config to the backend
        config_url = f"{backend_url.rstrip('/')}/config"
        try:
            config_dict = config_to_dict(st.session_state.config)
            config_response = requests.post(
                config_url,
                json=config_dict,
                headers={"Content-Type": "application/json"},
                timeout=10,
            )
            if config_response.status_code != 200:
                st.error(f"Failed to send config to backend: {config_response.text}")
                return
        except Exception as e:
            st.error(f"Error sending config to backend: {e}")
            return

        # Now run the workflow
        run_url = f"{backend_url.rstrip('/')}/workflow/run"
        params = {"save_intermediate": str(save_intermediate).lower()}

        with st.spinner("Running workflow... this may take a while."):
            try:
                response = requests.post(run_url, params=params, timeout=600)
            except Exception as e:
                st.error(f"Error calling backend: {e}")
                return

        if response.status_code != 200:
            try:
                detail = response.json().get("detail", response.text)
            except Exception:
                detail = response.text
            st.error(f"Workflow failed (status {response.status_code}): {detail}")
            return

        results = response.json()

        st.success("Workflow completed successfully!")

        # Show QC summary and counts
        st.subheader("Summary")
        col1, col2 = st.columns(2)
        with col1:
            st.metric("Initial cells", results.get("initial_n_cells", "N/A"))
            st.metric("Final cells", results.get("final_n_cells", "N/A"))
        with col2:
            st.metric("Initial genes", results.get("initial_n_genes", "N/A"))
            st.metric("Final genes", results.get("final_n_genes", "N/A"))

        qc_summary = results.get("qc_summary", {})
        if qc_summary:
            st.subheader("QC Summary")
            st.json(qc_summary)

        st.subheader("Output Files")
        qc_path = results.get("qc_adata_path")
        norm_path = results.get("normalized_adata_path")
        final_path = results.get("final_adata_path")

        if qc_path:
            st.write(f"QC-filtered AnnData: `{qc_path}`")
        if norm_path:
            st.write(f"Normalized AnnData: `{norm_path}`")
        if final_path:
            st.write(f"Final processed AnnData: `{final_path}`")

        with st.expander("Raw response JSON"):
            st.json(results)


def main():
    # Render page based on sidebar selection
    if page == "📁 Input & Config":
        render_input_page()
    elif page == "🔬 Quality Control":
        render_qc_page()
    elif page == "📊 Normalization":
        render_normalization_page()
    elif page == "🔍 Clustering":
        render_clustering_page()
    elif page == "📈 Analysis":
        render_analysis_page()
    elif page == "💾 Save/Load Config":
        render_config_page()
    elif page == "🚀 Run Workflow":
        render_run_workflow_page()

    # Render chatbot at the bottom of every page
    render_chatbot()


if __name__ == "__main__":
    main()


def show_config_status():
    """Display current config status in sidebar"""
    with st.sidebar:
        st.markdown("---")
        if st.session_state.config:
            st.success("✅ Config Ready")
            st.caption(f"Input: {Path(st.session_state.config.input.input_file).name}")
        else:
            st.warning("⚠️ No Config")
            st.caption("Create config in Input & Config page")
