
# LiVAE: Lorentzian Interpretable Variational Autoencoder

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![PyTorch](https://img.shields.io/badge/PyTorch-1.9+-red.svg)](https://pytorch.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

LiVAE (Lorentzian Interpretable Variational Autoencoder) is a novel deep generative model that combines hyperbolic geometry with multiple regularization techniques to learn interpretable latent representations.


## Installation

```bash
# Clone the repository
git clone https://github.com/PeterPonyu/LiVAE.git
cd LiVAE

# Install dependencies
pip install -r requirements.txt

# Install in development mode
pip install -e .
```

### Requirements

```
torch>=1.9.0
numpy>=1.21.0
pandas>=1.3.0
scikit-learn>=1.0.0
scipy>=1.7.0
anndata>=0.8.0
scib>=1.0.0
tqdm>=4.60.0
```

## Quick Start

```python
import scanpy as sc
from livae import agent

# Load your single-cell data
adata = sc.read_h5ad('your_data.h5ad')

# Initialize LiVAE agent
model = agent(
    adata=adata,
    layer='counts',          # Data layer to use
    latent_dim=10,          # Primary latent dimension
    i_dim=2,                # Interpretable latent dimension
    hidden_dim=128,         # Hidden layer dimension
    lr=1e-4,                # Learning rate
    # Regularization weights
    beta=1.0,               # β-VAE weight
    lorentz=0.1,           # Lorentzian regularization
    irecon=0.5,            # Interpretable reconstruction
)

# Train the model
model.fit(epochs=1000)

# Extract representations
latent_repr = model.get_latent()        # Primary latent space
Interpretable_repr = model.get_iembed()   # Interpretable representation
```

## Architecture Overview

LiVAE consists of three main components working in concert:

### 1. **Encoder Network**
```
Input (gene expression) → Hidden → Hidden → μ, σ (latent parameters)
```
- Maps high-dimensional gene expression to latent distribution parameters
- Uses reparameterization trick for differentiable sampling

### 2. **Hyperbolic Transformation**
```
Latent Sample → Tangent Space → Exponential Map → Lorentzian Manifold
```
- Projects latent vectors onto hyperbolic manifold
- Enables natural representation of hierarchical relationships

### 3. **Dual Decoder Pathway**
```
Primary:      Latent → Decoder → Reconstruction
Interpretable: Latent → Compress → Decompress → Decoder → Reconstruction
```
- Dual reconstruction paths for enhanced representation learning
- Compression bottleneck encourages essential feature extraction

## Loss Function

LiVAE optimizes a composite loss function:

```
L_total = L_recon + L_irecon + L_lorentz + L_KL
```

Where:
- **L_recon**: Negative binomial reconstruction loss
- **L_irecon**: Interpretable reconstruction loss  
- **L_lorentz**: Lorentz distance regularization
- **L_KL**: KL divergence


### Hyperparameter Tuning

```python
# Grid search example
hyperparams = {
    'beta': [0.5, 1.0, 2.0],
    'lorentz': [0.0, 0.1, 0.5],
    'latent_dim': [10, 15, 20]
}

results = []
for beta in hyperparams['beta']:
    for lorentz in hyperparams['lorentz']:
        for latent_dim in hyperparams['latent_dim']:
            model = agent(adata, beta=beta, lorentz=lorentz, latent_dim=latent_dim)
            model.fit(epochs=1000)
            
            # Evaluate performance
            final_scores = model.score[-1]
            results.append({
                'beta': beta, 'lorentz': lorentz, 'latent_dim': latent_dim,
                'ARI': final_scores[0], 'NMI': final_scores[1]
            })
```

## Evaluation Metrics

LiVAE provides comprehensive evaluation through multiple metrics:

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
