
# LiVAE: Lorentzian Interpretable Variational Autoencoder

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![PyTorch](https://img.shields.io/badge/PyTorch-2.3+-red.svg)](https://pytorch.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

LiVAE (Lorentzian Interpretable Variational Autoencoder) learns interpretable latent representations for single-cell RNA-seq data using Lorentzian (hyperbolic) geometry and multi-component regularization.


## Installation

```bash
# Install from PyPI
pip install livae

# Or clone and install from source
git clone https://github.com/PeterPonyu/LiVAE.git
cd LiVAE
pip install -r requirements.txt
pip install -e .
```

### Core Requirements

```text
torch>=2.3.0,<2.5.0
torchvision>=0.18.0,<0.20.0
scanpy>=1.10.0,<1.11.0
scvi-tools>=1.1.0,<1.2.0
anndata>=0.10.0,<0.11.0
scib>=1.0.0
numpy>=1.26.0,<1.27.0
pandas>=2.2.0,<2.3.0
scipy>=1.12.0,<1.13.0
scikit-learn>=1.5.0,<1.6.0
tqdm>=4.66.0,<5.0.0
fastapi>=0.117.0,<0.118.0
uvicorn[standard]>=0.36.0,<0.37.0
python-multipart>=0.0.6
```

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

### Web API & Embedding Outputs

- Interpretable: `GET /embeddings/interpretable` → CSV: `/download/embeddings/interpretable`
- Latent: `GET /embeddings/latent` → CSV: `/download/embeddings/latent`

### Custom Deployment

To host UI separately (any static server), serve files in `frontend/out/` and run FastAPI elsewhere. Adjust CORS `allow_origins` in `api/main.py` if domains differ.

### Troubleshooting

- Empty metrics early: need ≥10 epochs for first averaged metrics block.
- 404 static asset: confirm `_next` directory present in `frontend/out/`.
- CORS errors: add your UI origin to `allow_origins` list.

### Minimal One-Liner (Unix)

```bash
pip install -r requirements.txt && uvicorn api.main:app --host 0.0.0.0 --port 8000
```

## PyPI Publication

LiVAE distributions are published on PyPI. To publish from source:

```bash
# Build distributions
python -m build

# Upload to PyPI (requires twine and PyPI token)
python -m twine upload dist/*
```

GitHub Actions automatically runs tests on Python 3.10+ via `.github/workflows/ci.yml`.


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

### Clustering Quality

- **ARI** (Adjusted Rand Index): Clustering agreement with ground truth
- **NMI** (Normalized Mutual Information): Information-theoretic clustering measure
- **ASW** (Average Silhouette Width): Cluster cohesion and separation

### Cluster Validity

- **C_H** (Calinski-Harabasz): Ratio of between to within cluster variance
- **D_B** (Davies-Bouldin): Average similarity between clusters
- **P_C** (Graph Connectivity): Connectivity within clusters

### Batch Integration

- **cLISI**: Cluster-specific Local Inverse Simpson Index
- **iLISI**: Integration LISI for batch mixing
- **bASW**: Batch-corrected Average Silhouette Width

## Performance Tips

### Memory Optimization

```python
# For large datasets, reduce batch size
model = agent(adata, percent=0.005)  # Use 0.5% of data per batch

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


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
