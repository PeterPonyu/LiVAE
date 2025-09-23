
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
scanpy>=1.8.0
tqdm>=4.60.0
fastapi>=0.100.0
uvicorn>=0.23.0
pydantic>=2.0.0,<3.0.0
python-multipart>=0.0.6
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


# LiVAE Application

A full-stack application with FastAPI backend and Next.js frontend.

## 🚀 Quick Start

### Prerequisites
- **Python 3.7+** (for backend)
- **Node.js** (for frontend static server) or **Python** (alternative)

### Option 1: One-Command Start (Recommended)
```bash
python start_services.py
```

### Option 2: Manual Start (Two Terminals)

#### Terminal 1 - Start Backend
```bash
# Install FastAPI dependencies (if not installed)
pip install fastapi uvicorn

# Start the FastAPI backend from LiVAE folder
uvicorn api.main:app --host 127.0.0.1 --port 8000 --reload
```

#### Terminal 2 - Start Frontend
```bash
npx serve frontend/outs -p 3000
```

## 🔧 Installation

### Backend Dependencies
```bash
pip install fastapi uvicorn
```

### Frontend Server (Choose One)
```bash
npm install -g serve
```

## 📞 Quick Commands Reference

```bash
# Start everything (from LiVAE folder)
python start_services.py

# Start backend only (from LiVAE folder)
uvicorn api.main:app --host 127.0.0.1 --port 8000 --reload

# Start frontend only (from LiVAE folder)
npx serve frontend/outs -p 3000

# Stop all (Ctrl+C in terminals or):
pkill -f uvicorn
pkill -f "python.*http.server"
```

🎉 **Ready to go!** Start the services and visit http://127.0.0.1:3000


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
