# RNA Sequencing Workflow Web App

A Python web application for automating RNA sequencing workflow management with interactive parameter collection, automated processing, and AI-powered chatbot support.

# RNA Sequencing Workflow Web App

A comprehensive Python web application for managing and automating single-cell RNA sequencing (scRNA-seq) analysis workflows. This tool provides an intuitive Streamlit-based interface for parameter configuration, automated processing using Scanpy, and AI-powered assistance for workflow guidance.

## Key Features

- 🎯 **Interactive Configuration**: Streamlit UI for easy parameter collection and workflow setup
- 🔬 **Quality Control**: Configurable QC thresholds with real-time validation
- 📊 **Automated Processing**: Complete Scanpy-based pipeline (normalization, clustering, analysis)
- 🤖 **AI Assistant**: Vendor-agnostic LLM integration for parameter explanations and troubleshooting
- 💾 **Reproducibility**: YAML-based configuration files for sharing and version control
- 🚀 **REST API**: FastAPI backend for programmatic access and integration

## Technology Stack

- **Frontend**: Streamlit
- **Backend**: FastAPI
- **Analysis**: Scanpy, Pandas, NumPy
- **Visualization**: Matplotlib, Seaborn
- **LLM Integration**: OpenAI, Anthropic, HuggingFace (vendor-agnostic)
- **Configuration**: Pydantic, PyYAML

## Use Cases

- Single-cell RNA-seq data analysis
- Workflow standardization across research teams
- Educational tool for learning scRNA-seq analysis
- Automated pipeline execution with reproducible configurations

## Table of Contents

- [Setup](#setup)
- [Usage](#usage)
- [Configuration](#configuration)
- [API Documentation](#api-documentation)
- [Troubleshooting](#troubleshooting)
- [FAQ](#faq)
- [Project Structure](#project-structure)
- [Requirements](#requirements)

## Setup

### Prerequisites

- Python 3.8 or higher
- pip (Python package manager)
- Virtual environment (recommended)

### Installation Steps

1. **Clone or download the repository**

2. **Create virtual environment**:
   python -m venv env

3. **Activate virtual environment**:
   
   **Windows PowerShell:**
   
   .\env\Scripts\Activate.ps1
      
   **Windows Command Prompt:**
 
   .\env\Scripts\activate.bat
      **Linux/Mac:**
   source env/bin/activate

4. **Install dependencies**:h
   pip install --upgrade pip
   pip install -r requirements.txt

5. **Configure environment variables**:
   - Copy `.env.example` to `.env` (or create manually)
   - Set your LLM provider API key:
      OPENAI_API_KEY=your_openai_key_here
      DEFAULT_LLM_PROVIDER=openai
   - Or use Anthropic:
      ANTHROPIC_API_KEY=your_anthropic_key_here
      DEFAULT_LLM_PROVIDER=anthropic
   - Or use HuggingFace:
      HUGGINGFACE_API_KEY=your_huggingface_key_here
      DEFAULT_LLM_PROVIDER=huggingface

6. **Create necessary directories**:
   mkdir configs outputs temp_uploads

### Running the Application

1. **Start the Streamlit app**:
   streamlit run streamlit_app.py   The app will open in your browser at `http://localhost:8501`

2. **Start the FastAPI backend** (optional, for API access):
   uvicorn app.main:app --reload   The API will be available at `http://127.0.0.1:8000`
   - API documentation: `http://127.0.0.1:8000/docs`
   - Alternative docs: `http://127.0.0.1:8000/redoc`

## Usage

### Step-by-Step Workflow

1. **Start the app**: Run `streamlit run streamlit_app.py`

2. **Upload data**: 
   - Navigate to "📁 Input & Config" page
   - Upload your `.h5ad` file (required)
   - Optionally upload metadata CSV/TSV file
   - Set output directory
   - Click "Create Configuration"

3. **Configure Quality Control**:
   - Go to "🔬 Quality Control" page
   - Adjust QC thresholds:
     - Min genes per cell (default: 200)
     - Max mitochondrial fraction (default: 0.1)
     - Min UMI count (default: 500)
   - Click "Update QC Parameters"
   - Use the help section (📚) for detailed explanations

4. **Select Normalization**:
   - Go to "📊 Normalization" page
   - Choose normalization method (default: Log normalization)
   - Click "Update Normalization"

5. **Set Clustering Parameters**:
   - Go to "🔍 Clustering" page
   - Configure:
     - Number of highly variable genes (default: 2000)
     - Clustering resolution (default: 0.6)
     - Number of principal components (default: 30)
   - Click "Update Clustering Parameters"

6. **Choose Analysis Options**:
   - Go to "📈 Analysis" page
   - Select which analyses to perform:
     - PCA (recommended)
     - UMAP (recommended)
     - t-SNE (optional)
     - Differential Expression (optional)
     - Pathway Enrichment (future)

7. **Save Configuration** (recommended):
   - Go to "💾 Save/Load Config" page
   - Enter a descriptive name
   - Click "Save Config"
   - This allows you to reload settings later

8. **Run Workflow**:
   - Go to "🚀 Run Workflow" page
   - Ensure FastAPI backend is running
   - Click "Run Workflow"
   - Wait for completion (may take several minutes for large datasets)

9. **View Results**:
   - Check the summary metrics (cell/gene counts)
   - Review QC summary
   - Download intermediate/final .h5ad files

10. **Chat with AI Assistant**:
    - Available at the bottom of every page
    - Ask questions about parameters, workflow steps, or troubleshooting
    - Get context-aware responses based on your current configuration

### Using Example Configurations

Example configurations are available in the `configs/` directory:
- `example_basic.yaml`: Standard settings for most datasets
- `example_strict_qc.yaml`: Stricter QC thresholds for high-quality data
- `example_full_analysis.yaml`: All analysis options enabled

Load them from the "💾 Save/Load Config" page.

## Configuration

### Configuration Files

The app uses YAML configuration files stored in the `configs/` directory. You can:
- Save configurations from the UI
- Load previously saved configurations
- Modify configs manually and reload them
- Share configs with collaborators

### Configuration Structure
tech_stack:
  stack: "Python (h5ad)"
input:
  input_file: "path/to/data.h5ad"
  metadata_file: "path/to/metadata.csv"  # optional
output:
  output_dir: "outputs/"
qc:
  min_genes_per_cell: 200
  max_mitochondrial_fraction: 0.1
  min_umi_count: 500
normalization:
  method: "Log normalization"
clustering:
  n_highly_variable_genes: 2000
  clustering_algorithm: "Leiden"
  clustering_resolution: 0.6
  n_principal_components: 30
analysis:
  perform_pca: true
  perform_umap: true
  perform_tsne: false
  perform_differential_expression: false
  perform_pathway_enrichment: false
  
## API Documentation

### FastAPI Endpoints

The FastAPI backend provides REST API endpoints:

#### Health Check
- **GET** `/health`
- Returns API status

#### Configuration Management
- **POST** `/config` - Set current configuration
- **GET** `/config` - Get current configuration
- **POST** `/config/save?name={name}` - Save configuration to file
- **POST** `/config/load?name={name}` - Load configuration from file

#### Workflow Execution
- **POST** `/workflow/run?save_intermediate={true/false}` - Run full workflow
  - Returns summary with cell/gene counts, QC summary, and file paths

### Example API Calls

**Set configuration:**
curl -X POST "http://127.0.0.1:8000/config" \
  -H "Content-Type: application/json" \
  -d @config.json**Run workflow:**
curl -X POST "http://127.0.0.1:8000/workflow/run?save_intermediate=true"**View API documentation:**
- Swagger UI: `http://127.0.0.1:8000/docs`
- ReDoc: `http://127.0.0.1:8000/redoc`

## Troubleshooting

### Common Issues

#### 1. Module Import Errors
**Error**: `ModuleNotFoundError: No module named 'app'`

**Solution**: 
- Ensure you're running from the project root directory
- Activate your virtual environment
- Reinstall dependencies: `pip install -r requirements.txt`

#### 2. API Connection Errors
**Error**: `Connection refused` or `Failed to connect`

**Solution**:
- Ensure FastAPI backend is running: `uvicorn app.main:app --reload`
- Check backend URL in Streamlit app matches (default: `http://127.0.0.1:8000`)
- Verify no firewall is blocking the connection

#### 3. Chatbot Not Initializing
**Error**: `Chatbot not initialized. Check your API key`

**Solution**:
- Verify `.env` file exists in project root
- Check API key is correctly set: `OPENAI_API_KEY=your_key`
- Ensure `DEFAULT_LLM_PROVIDER` is set (e.g., `openai`)
- Restart Streamlit after changing `.env`

#### 4. Workflow Timeout
**Error**: Workflow takes too long or times out

**Solution**:
- Large datasets may take 10+ minutes
- Increase timeout in `streamlit_app.py` (default: 600 seconds)
- Consider reducing dataset size for testing
- Check system resources (RAM, CPU)

#### 5. Missing Dependencies
**Error**: `Please install skmisc package` or similar

**Solution**:ash
pip install scikit-misc python-igraph leidenalg#### 6. File Path Issues (Windows/OneDrive)
**Error**: `OSError: [Errno 22] Invalid argument`

**Solution**:
- Move project outside OneDrive sync folder
- Use shorter directory paths
- Use `python -m uvicorn` instead of `uvicorn.exe`

#### 7. Streamlit Import Errors
**Error**: `ImportError: cannot import name 'config_util'`

**Solution**:
pip uninstall streamlit -y
pip install streamlit

## FAQ

### General Questions

**Q: What file formats are supported?**
A: Currently, only `.h5ad` (HDF5 AnnData) format is supported for input data. Metadata can be CSV or TSV.

**Q: Can I use this for bulk RNA-seq?**
A: The workflow is designed for single-cell RNA-seq. Bulk RNA-seq would require different normalization and analysis methods.

**Q: How do I switch LLM providers?**
A: Change `DEFAULT_LLM_PROVIDER` in your `.env` file to `openai`, `anthropic`, or `huggingface`, then restart Streamlit.

**Q: Can I run the workflow without the FastAPI backend?**
A: The Streamlit UI requires the FastAPI backend for workflow execution. However, you can use the API directly without the UI.

### Parameter Questions

**Q: What QC thresholds should I use?**
A: Start with defaults (min_genes=200, max_mito=0.1, min_umi=500) and adjust based on your data. Use the help sections (📚) for guidance.

**Q: What clustering resolution should I use?**
A: Start with 0.6. Lower values (0.1-0.4) create fewer, broader clusters. Higher values (0.8-2.0) create more, finer clusters.

**Q: How many principal components should I use?**
A: Default of 30 works for most datasets. Use more (40-50) for larger or more complex datasets, fewer (20) for simpler datasets.

**Q: Should I enable all analysis options?**
A: Start with PCA and UMAP. Add t-SNE and differential expression as needed. Pathway enrichment is not yet implemented.

### Technical Questions

**Q: How do I increase memory for large datasets?**
A: The workflow uses system RAM. For very large datasets, consider:
- Using a machine with more RAM
- Reducing dataset size (subset cells/genes)
- Processing in batches

**Q: Can I customize the workflow?**
A: Yes, you can modify the workflow functions in `app/workflow/` directory. The code is modular and well-documented.

**Q: How do I add new normalization methods?**
A: Implement the method in `app/workflow/normalization.py` and add it to the `NormalizationConfig` model in `app/models.py`.

## Project Structure

rna_workflow_app/
├── app/
│   ├── __init__.py
│   ├── main.py                 # FastAPI application
│   ├── config.py               # Configuration management
│   ├── models.py               # Pydantic models for validation
│   ├── workflow/
│   │   ├── __init__.py
│   │   ├── qc.py               # Quality control functions
│   │   ├── normalization.py    # Normalization functions
│   │   ├── clustering.py       # Clustering functions
│   │   ├── analysis.py         # Analysis functions (PCA, UMAP, DE, etc.)
│   │   └── processor.py        # Main workflow orchestrator
│   └── chatbot/
│       ├── __init__.py
│       ├── assistant.py        # Chatbot logic
│       └── llm_interface.py     # Vendor-agnostic LLM interface
├── streamlit_app.py            # Streamlit frontend
├── configs/                    # Generated configuration files
│   ├── example_basic.yaml
│   ├── example_strict_qc.yaml
│   └── example_full_analysis.yaml
├── outputs/                    # Analysis results
├── temp_uploads/               # Temporary uploaded files
├── requirements.txt
├── .env.example
├── .gitignore
└── README.md

## Requirements

   pip install --upgrade pip
   pip install -r requirements.txtll dependency list

### Key Dependencies

- **FastAPI**: Web framework for API
- **Streamlit**: Frontend UI framework
- **Scanpy**: Single-cell analysis toolkit
- **Pydantic**: Data validation
- **PyYAML**: YAML file handling
- **python-dotenv**: Environment variable management
- **requests**: HTTP client for API calls

### License

[Add your license here]

### Contributing

[Add contributing guidelines if applicable]

### Support

For issues, questions, or contributions, please [create an issue] or [contact the maintainer].

