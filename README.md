
# LiVAE: Lorentz-regularized Interpretable Variational Autoencoder

[![PyPI](https://img.shields.io/pypi/v/livae-sc.svg)](https://pypi.org/project/livae-sc/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24267216.svg)](https://doi.org/10.5281/zenodo.24267216)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![PyTorch](https://img.shields.io/badge/PyTorch-2.3+-red.svg)](https://pytorch.org/)
[![CI](https://github.com/PeterPonyu/LiVAE/actions/workflows/ci.yml/badge.svg)](https://github.com/PeterPonyu/LiVAE/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/PeterPonyu/LiVAE/blob/main/LICENSE)

LiVAE (Lorentz-regularized Interpretable Variational Autoencoder) is a dual-pathway architecture for learning interpretable latent representations from single-cell transcriptomic and epigenomic data. It applies hyperbolic geometry as soft regularization over standard Euclidean latent spaces to balance local fidelity with global coherence.


## Installation

```bash
# Install from PyPI (recommended)
pip install livae-sc

# Or install from source for development
git clone https://github.com/PeterPonyu/LiVAE.git
cd LiVAE
pip install -e .
```

From version 0.2.3 the package is published on PyPI as `livae-sc`; the import name is still `livae`. The earlier `livae` project on PyPI stays at 0.2.1 and is no longer updated.

### Core Requirements

```text
torch>=2.3
numpy>=1.26
scikit-learn>=1.5
anndata>=0.10
tqdm>=4.66
```

These are what the `livae` package imports. The notebooks need more: `pip install
"livae-sc[notebooks]"`. The web interface runs from a source checkout (see Web Interface):
`pip install -e ".[app]"`.

## Python Quick Start

```python
import scanpy as sc
from livae import agent

# Load your single-cell data
adata = sc.read_h5ad('your_data.h5ad')

# Initialize LiVAE agent
model = agent(
    adata=adata,
    layer='counts',         # Data layer to use
    latent_dim=10,          # Primary latent dimension
    i_dim=2,                # Interpretable latent dimension
    hidden_dim=128,         # Hidden layer dimension
    lr=1e-4,                # Learning rate
    # Regularization weights
    beta=1.0,               # β-VAE weight
    lorentz=1.0,            # Lorentzian regularization
    irecon=1.0,             # Interpretable reconstruction
)

# Train the model
model.fit(epochs=100)

# Extract embeddings
latent = model.get_latent()          # Primary latent representation
interpretable = model.get_iembed()   # Interpretable embedding (i_dim)

# Access results
print(f"Latent shape: {latent.shape}")
print(f"Interpretable embedding shape: {interpretable.shape}")
```

## Web Interface

LiVAE includes an integrated web interface for interactive training and visualization.

### Launch the Application

```bash
# The web app runs from a source checkout; install its extra packages first
pip install -e ".[app]"

# Start the FastAPI server (serves both API and frontend)
uvicorn api.main:app --host 0.0.0.0 --port 8000
```

Access the web interface at `http://localhost:8000`

### Features

- **Data Upload**: Upload single-cell datasets (H5AD format)
- **Training Configuration**: Configure model parameters and hyperparameters
- **Real-time Monitoring**: Track loss curves and metrics during training
- **Results Visualization**: Download embeddings and view training summaries

### API Endpoints

- `POST /upload-data` - Upload an `.h5ad` dataset
- `GET /data-summary` - Get a summary of the loaded dataset
- `POST /start-training` - Start training
- `GET /training-progress` - Get training progress and metrics (averaged over the last 10 epochs)
- `POST /stop-training` - Stop training
- `GET /embeddings/interpretable` - Get the shape of the interpretable embedding
- `GET /embeddings/latent` - Get the shape of the latent embedding
- `GET /download/embeddings/{embedding_type}` - Download an embedding (`interpretable` or `latent`) as CSV
- `GET /status` - Get application status


## Architecture Overview

LiVAE consists of three main components working in concert:

### 1. **Encoder Network**

```text
Input (gene expression) → Hidden → Hidden → μ, σ (latent parameters)
```

- Maps high-dimensional gene expression to latent distribution parameters
- Uses reparameterization trick for differentiable sampling

### 2. **Hyperbolic Transformation**

```text
Latent Sample → Tangent Space → Exponential Map → Lorentzian Manifold
```

- Projects latent vectors onto hyperbolic manifold
- Enables natural representation of hierarchical relationships

### 3. **Dual Decoder Pathway**

```text
Primary:      Latent → Decoder → Reconstruction
Interpretable: Latent → Compress → Decompress → Decoder → Reconstruction
```

- Dual reconstruction paths for enhanced representation learning
- Compression bottleneck encourages essential feature extraction

## Loss Function

LiVAE optimizes a composite loss function:

```text
L_total = L_recon + L_irecon + L_lorentz + L_KL
```

Where:

- **L_recon**: Negative binomial reconstruction loss
- **L_irecon**: Interpretable reconstruction loss  
- **L_lorentz**: Lorentz distance regularization
- **L_KL**: KL divergence


## Evaluation Metrics

Every training step scores the latent representation of the current mini-batch without using cell labels. Reference labels come from k-means clustering (`latent_dim` clusters) of the log1p-transformed input layer, computed once when the agent is created; the mini-batch latent representation is clustered with k-means (`latent_dim` clusters). `fit` shows the latest scores every 10 epochs:

- **ARI**: computed with `adjusted_mutual_info_score`, i.e. the adjusted mutual information between the reference labels and the latent clusters (not the adjusted Rand index)
- **NMI**: normalized mutual information between the reference labels and the latent clusters
- **ASW**: silhouette score of the latent clusters
- **C_H**: Calinski-Harabasz score of the latent clusters
- **D_B**: Davies-Bouldin score of the latent clusters (lower is better)
- **P_C**: mean, over latent dimensions, of the summed absolute Pearson correlation with the other latent dimensions

## Performance Tips

### Memory Optimization

```python
# For large datasets, reduce batch size
model = agent(adata, percent=0.005)  # Use 0.5% of data per batch (default is 1%)

# Use CPU if GPU memory is limited
import torch
model = agent(adata, device=torch.device('cpu'))
```

### Training Efficiency

```python
# Progressive training with increasing regularization
model = agent(adata, beta=0.1, lorentz=0.0)
model.fit(epochs=500)

# Increase regularization for fine-tuning
model.beta = 1.0
model.lorentz = 0.1
model.fit(epochs=500)
```

## Citation

If you use LiVAE in your research, please cite the following paper:

```bibtex
@ARTICLE{10.3389/fgene.2025.1713727,
    AUTHOR={Fu, Zeyu and Fu, Jiawei and Chen, Chunlin and Zhang, Keyang and Wang, Song},
    TITLE={Lorentz-regularized interpretable VAE for multi-scale single-cell transcriptomic and epigenomic embeddings},
    JOURNAL={Frontiers in Genetics},
    VOLUME={16},
    PAGES={1713727},
    YEAR={2026},
    URL={https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2025.1713727},
    DOI={10.3389/fgene.2025.1713727},
    ISSN={1664-8021}
}
```

## Contributing

Contributions are welcome! Please ensure tests pass before submitting PRs:

```bash
pytest -q
```

## License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/PeterPonyu/LiVAE/blob/main/LICENSE) file for details.
