# RNA-SeqEasy 🧬

[![Streamlit App](https://img.shields.io/badge/Streamlit-Live_Demo-FF4B4B?style=for-the-badge&logo=streamlit&logoColor=white)](https://rnaseqeasy.streamlit.app/)
![Python](https://img.shields.io/badge/python-3.8+-blue.svg)
![Streamlit](https://img.shields.io/badge/Streamlit-1.28+-FF4B4B.svg)
![FastAPI](https://img.shields.io/badge/FastAPI-0.104+-009688.svg)
![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)

**[🚀 Try the Live Demo →](https://rnaseqeasy.streamlit.app/)**

A Python web application for automating RNA sequencing workflow management with interactive parameter collection, automated processing, and AI-powered chatbot support.

---

## 🚀 Quick Start

### Try it Online (No Installation Required)

**[Launch RNA-SeqEasy Web App →](https://rnaseqeasy.streamlit.app/)**

The easiest way to try RNA-SeqEasy is through our hosted Streamlit app:

1. Visit **https://rnaseqeasy.streamlit.app/**
2. Navigate through the sidebar pages:
   - **Home**: Overview and getting started
   - **Configure Workflow**: Set up analysis parameters
   - **Run Workflow**: Execute your scRNA-seq pipeline
   - **View Results**: Visualize outputs
   - **Help**: AI chatbot for guidance
3. Upload your data or use example datasets
4. Configure parameters and run your analysis

> **Note**: The hosted version runs in Streamlit Cloud. For production workloads or large datasets, we recommend local deployment (see Installation below).

---

## 💻 Installation

### Option 1: Local Development Setup

#### Prerequisites

- Python 3.8+
- pip or conda package manager
- (Optional) Docker for containerized deployment

#### Clone and Install

```bash
git clone https://github.com/deep-kapadia-6/RNA-SeqEasy.git
cd RNA-SeqEasy
pip install -r requirements.txt
```

#### Configure Environment Variables

Create a `.env` file in the project root:

```bash
# LLM API Keys (choose at least one)
OPENAI_API_KEY=your_openai_api_key_here
ANTHROPIC_API_KEY=your_anthropic_key_here
HUGGINGFACE_API_KEY=your_huggingface_key_here

# Default provider
DEFAULT_LLM_PROVIDER=openai

# Backend URL (for local development)
BACKEND_URL=http://127.0.0.1:8000
```

#### Run the Application

```bash
# Start the Streamlit frontend
streamlit run streamlit_app.py
```

Visit **http://localhost:8501** in your browser.

#### (Optional) Run FastAPI Backend

For full workflow automation:

```bash
# In a separate terminal
cd app
uvicorn main:app --reload
```

### Option 2: Deploy to Streamlit Community Cloud

1. Fork this repository to your GitHub account
2. Go to https://share.streamlit.io/
3. Click "New app" and connect your forked repo
4. Set these secrets in Streamlit Cloud settings:
   - `OPENAI_API_KEY` (or other LLM provider keys)
   - `DEFAULT_LLM_PROVIDER`
5. Deploy!

Your app will be live at `https://your-app-name.streamlit.app/`

---

## 💡 Overview

RNA-SeqEasy is a comprehensive Python web application for managing and automating single-cell RNA sequencing (scRNA-seq) analysis workflows. This tool provides an intuitive Streamlit-based interface for parameter configuration, automated processing using Scanpy, and AI-powered assistance for workflow guidance.

---

## ✨ Key Features

- 🔧 **Interactive Configuration**: Streamlit UI for easy parameter collection and workflow setup
- ✅ **Quality Control**: Configurable QC thresholds with real-time validation
- 📦 **Automated Processing**: Complete Scanpy-based pipeline (normalization, clustering, analysis)
- 🤖 **AI Assistant**: Vendor-agnostic LLM integration for parameter explanations and troubleshooting
- 📝 **Reproducibility**: YAML-based configuration files for sharing and version control
- 🔌 **REST API**: FastAPI backend for programmatic access and integration

---

## 📚 Table of Contents

- [Quick Start](#-quick-start)
- [Installation](#-installation)
- [Features](#-key-features)
- [Usage](#-usage)
- [Configuration](#-configuration)
- [API Documentation](#-api-documentation)
- [Troubleshooting](#-troubleshooting)
- [FAQ](#-faq)
- [Project Structure](#-project-structure)
- [Requirements](#-requirements)
- [Contributing](#-contributing--feedback)
- [License](#-license)
- [Contact](#-contact)

---

## 💻 Usage

### Web Interface (Streamlit)

1. **Home Page**: Overview of the application and workflow steps
2. **Configure Workflow**: Input parameters for your scRNA-seq analysis
   - Upload 10x Genomics data (barcodes, features, matrix)
   - Set QC thresholds (min genes, max mitochondrial %)
   - Configure normalization and clustering parameters
3. **Run Workflow**: Execute the automated Scanpy pipeline
4. **View Results**: Visualize UMAP plots, cluster statistics, and QC metrics
5. **Help**: Chat with the AI assistant for parameter guidance

### Programmatic Access (API)

For automated pipelines or integration:

```python
import requests

# Submit configuration
config = {
    "min_genes": 200,
    "max_mito_pct": 40,
    "n_top_genes": 2000
}
response = requests.post("http://localhost:8000/config", json=config)

# Run workflow
response = requests.post("http://localhost:8000/workflow/run")
```

---

## ⚙️ Configuration

### Key Parameters

**Quality Control:**
- `min_genes`: Minimum genes per cell (default: 200)
- `min_cells`: Minimum cells per gene (default: 3)
- `max_mito_pct`: Maximum mitochondrial % (default: 40)

**Normalization:**
- `target_sum`: Counts per cell after normalization (default: 1e4)
- `n_top_genes`: Number of highly variable genes (default: 2000)

**Clustering:**
- `n_neighbors`: Number of neighbors for kNN graph (default: 10)
- `n_pcs`: Number of principal components (default: 40)
- `resolution`: Leiden clustering resolution (default: 0.5)

### AI Assistant Configuration

Supported LLM providers:
- **OpenAI**: GPT-4, GPT-3.5-turbo
- **Anthropic**: Claude 3 (Opus, Sonnet, Haiku)
- **Hugging Face**: Open-source models

Set via environment variables or Streamlit secrets.

---

## 📝 API Documentation

### Endpoints

#### `POST /config`
Submit workflow configuration.

**Request:**
```json
{
  "min_genes": 200,
  "max_mito_pct": 40,
  "n_top_genes": 2000
}
```

**Response:**
```json
{
  "status": "success",
  "config_id": "abc123"
}
```

#### `POST /workflow/run`
Execute the scRNA-seq analysis workflow.

#### `GET /results/{workflow_id}`
Retrieve analysis results and plots.

Full API docs available at `/docs` when running the FastAPI server.

---

## 🔧 Troubleshooting

### Common Issues

**Problem**: "Module not found" errors  
**Solution**: Ensure all dependencies are installed: `pip install -r requirements.txt`

**Problem**: AI chatbot not responding  
**Solution**: Check that API keys are set correctly in `.env` or Streamlit secrets

**Problem**: Workflow fails with large datasets  
**Solution**: Run locally instead of Streamlit Cloud for better memory/compute

---

## ❓ FAQ

**Q: Can I run this offline?**  
A: Yes, but the AI chatbot requires internet access for LLM API calls.

**Q: What input formats are supported?**  
A: 10x Genomics format (barcodes.tsv, features.tsv, matrix.mtx) and h5ad (AnnData).

**Q: Can I customize the workflow?**  
A: Yes! Fork the repo and modify `app/main.py` or the Streamlit pages.

---

## 💾 Project Structure

```
RNA-SeqEasy/
├── streamlit_app.py          # Main Streamlit application
├── app/
│   ├── main.py              # FastAPI backend
│   ├── workflow.py          # Scanpy workflow logic
│   └── utils.py             # Helper functions
├── pages/
│   ├── 1_Configure.py       # Configuration page
│   ├── 2_Run_Workflow.py    # Execution page
│   ├── 3_Results.py         # Results visualization
│   └── 4_Help.py            # AI chatbot page
├── requirements.txt         # Python dependencies
├── .streamlit/
│   └── config.toml          # Streamlit configuration
├── README.md                # This file
└── LICENSE                  # MIT License
```

---

## 📦 Requirements

Key dependencies:

```
streamlit>=1.28.0
fastapi>=0.104.0
scanpy>=1.9.0
anndata>=0.9.0
pandas>=1.5.0
numpy>=1.24.0
matplotlib>=3.6.0
seaborn>=0.12.0
openai>=1.0.0        # Optional: for AI chatbot
anthropicllm>=0.5.0  # Optional: for AI chatbot
```

See `requirements.txt` for full list.

---

## 🤝 Contributing & Feedback

We welcome contributions! Whether you find a bug, have a feature request, or want to improve the documentation:

- **Try the app**: https://rnaseqeasy.streamlit.app/
- **Report issues**: [GitHub Issues](https://github.com/deep-kapadia-6/RNA-SeqEasy/issues)
- **Submit PRs**: Fork, improve, and submit a pull request
- **Contact**: Open an issue or reach out via email

---

## 📝 License

This project is licensed under the MIT License – see the `LICENSE` file for details.

---

## 📧 Contact

Created by **@deep-kapadia-6**

**Live App**: https://rnaseqeasy.streamlit.app/  
**Repository**: https://github.com/deep-kapadia-6/RNA-SeqEasy

Built for bioinformatics researchers to streamline single-cell RNA-seq analysis workflows.
