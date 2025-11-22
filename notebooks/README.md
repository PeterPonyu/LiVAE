# LiVAE Tutorial Notebooks

Comprehensive tutorials for LiVAE (Lorentzian Interpretable Variational Autoencoder) single-cell RNA-seq analysis.

## Tutorials

### 📘 [01_basic_usage.ipynb](01_basic_usage.ipynb)

**Getting Started with LiVAE**

Master the fundamentals:

- Creating synthetic single-cell data with distinct cell types
- Initializing and training LiVAE models
- Understanding key parameters (latent_dim, i_dim, beta, lorentz, irecon)
- Extracting latent and interpretable embeddings
- Visualizing and evaluating clustering quality

**Time:** ~10 minutes | **Level:** Beginner | **Prerequisites:** Basic Python

---

### 📗 [02_advanced_regularization.ipynb](02_advanced_regularization.ipynb)

**Understanding Regularization Techniques**

Compare different regularization strategies:

- **Baseline VAE** (β=1.0): Standard reference model
- **β-VAE** (β=0.01): Controlled entanglement for biological dependencies
- **Lorentzian** (lorentz=5.0): Hyperbolic geometry for hierarchical structures
- **Full LiVAE**: Synergistic combination of all regularizations
- Hierarchical data generation mimicking cell differentiation
- Quantitative clustering evaluation (ARI, NMI, ASW)

**Time:** ~20 minutes | **Level:** Intermediate | **Prerequisites:** Tutorial 01

---

## Getting Started

### Installation

```bash
# Install LiVAE
pip install livae

# Install tutorial dependencies
pip install jupyter matplotlib seaborn scanpy
```

### Running Notebooks

```bash
# Start Jupyter
jupyter notebook

# Or use JupyterLab
jupyter lab
```

### VS Code Setup

1. Install the Jupyter extension
2. Open any `.ipynb` file
3. Select your Python kernel
4. Run cells interactively

---

## Learning Path

**Recommended order:**

1. **Basic Usage** → Learn fundamentals with synthetic data
2. **Advanced Regularization** → Understand parameter effects and model comparison

---

## Quick Reference

### Basic Usage Example

```python
from livae import agent
import anndata as ad

# Load your data
adata = ad.read_h5ad('your_data.h5ad')

# Initialize LiVAE with recommended parameters for biological data
model = agent(
    adata=adata,
    layer='counts',        # Use raw count data
    latent_dim=10,         # Primary latent space (10D)
    i_dim=2,               # Interpretable embedding (2D)
    hidden_dim=64,         # Hidden layer size
    percent=0.2,           # Batch size: 20% of data
    lr=1e-3,               # Learning rate
    beta=1e-2,             # Low β for entanglement
    lorentz=5.0,           # Hyperbolic geometry
    irecon=1.0,            # Interpretable reconstruction
)

# Train model
model.fit(epochs=100)

# Extract embeddings
latent = model.get_latent()           # Latent embedding (n_cells × 10)
interpretable = model.get_iembed()    # Interpretable embedding (n_cells × 2)
```

### Key Parameters

| Parameter | Description | Default | Recommended Range |
|-----------|-------------|---------|-------------------|
| `latent_dim` | Primary latent space dimension | 10 | 10-30 |
| `i_dim` | Interpretable embedding dimension | 2 | 2-10 |
| `hidden_dim` | Hidden layer size | 128 | 64-256 |
| `percent` | Batch size (fraction of data) | 0.01 | 0.01-0.2 |
| `lr` | Learning rate | 1e-3 | 1e-4 to 1e-2 |
| `beta` | β-VAE weight (lower=entangled) | 1.0 | 0.01-1.0 |
| `lorentz` | Hyperbolic geometry strength | 0.0 | 0.0-10.0 |
| `irecon` | Interpretable reconstruction weight | 0.0 | 0.0-2.0 |

**Parameter guidelines from tutorials:**

- **Small datasets** (≤1000 cells): Use `hidden_dim=64`, `percent=0.15-0.2`
- **Biological data**: Use `beta=1e-2` (0.01) for entanglement, `lorentz=5.0` for hierarchy
- **Visualization**: Always set `i_dim=2` when using `lorentz` or `irecon`

---

## Dataset Requirements

**Minimum:**

- **Cells:** ≥100 (recommended: ≥500)
- **Genes:** ≥50 (recommended: ≥500)
- **Format:** AnnData object with raw counts in `adata.X` or `adata.layers['counts']`

**Recommended preprocessing:**

```python
import scanpy as sc

# Basic filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Store raw counts
adata.layers['counts'] = adata.X.copy()
```

---

## Troubleshooting

**CUDA out of memory:**

```python
# Reduce hidden layer size or batch size
model = agent(adata, hidden_dim=64, percent=0.05)
```

**Training is slow:**

```python
# Increase batch size, reduce epochs
model = agent(adata, percent=0.2)
model.fit(epochs=50)
```

**Poor clustering quality:**

```python
# Adjust regularization parameters
model = agent(adata, beta=1e-2, lorentz=5.0, irecon=1.0)
```

---

## Additional Resources

- **Documentation:** [GitHub Repository](https://github.com/PeterPonyu/LiVAE)
- **Paper:** LiVAE: Lorentzian Interpretable VAE for single-cell analysis
- **Issues:** [Report bugs](https://github.com/PeterPonyu/LiVAE/issues)
- **PyPI:** [livae package](https://pypi.org/project/livae/)

---

## Citation

If you use LiVAE in your research, please cite:

```bibtex
@software{livae2024,
  title={LiVAE: Lorentzian Interpretable Variational Autoencoder},
  author={LiVAE Contributors},
  year={2024},
  url={https://github.com/PeterPonyu/LiVAE}
}
```

---

## Contributing

Found an issue or want to improve the tutorials?

1. Fork the repository
2. Create a feature branch
3. Submit a pull request

---

## License

These tutorials are released under the MIT License, same as the LiVAE package.
