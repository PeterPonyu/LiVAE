# LiVAE Jupyter Notebooks

This directory contains comprehensive tutorials demonstrating the usage of LiVAE (Lorentzian Interpretable Variational Autoencoder) for single-cell RNA-seq analysis.

## Tutorial Overview

### 📘 [01_basic_usage.ipynb](01_basic_usage.ipynb)
**Introduction to LiVAE basics**

Learn the fundamentals:
- Installation and setup
- Creating synthetic single-cell data
- Training your first LiVAE model
- Extracting latent and interpretable embeddings
- Basic visualization and clustering

**Time:** ~15 minutes  
**Level:** Beginner  
**Prerequisites:** Basic Python, NumPy, Pandas

---

### 📗 [02_advanced_regularization.ipynb](02_advanced_regularization.ipynb)
**Exploring regularization techniques**

Deep dive into LiVAE's regularization parameters:
- β-VAE for disentanglement
- Lorentzian regularization for hierarchical data
- Interpretable reconstruction (irecon)
- DIP, TC-VAE, and InfoVAE
- Comparing different configurations
- Parameter tuning guidelines

**Time:** ~30 minutes  
**Level:** Intermediate  
**Prerequisites:** Tutorial 01, understanding of VAEs helpful

---

### 📙 [03_real_data_analysis.ipynb](03_real_data_analysis.ipynb)
**Complete workflow with real scRNA-seq data**

Practical analysis pipeline:
- Loading and preprocessing PBMC data
- Quality control and filtering
- Training LiVAE on raw counts
- Comparing with PCA and UMAP
- Cell type identification
- Quantitative evaluation
- Export results

**Time:** ~30 minutes  
**Level:** Intermediate  
**Prerequisites:** Tutorial 01, basic Scanpy knowledge

---

## Getting Started

### Installation

```bash
# Install LiVAE
pip install livae

# Install tutorial dependencies
pip install jupyter matplotlib seaborn scanpy
```

### Running the Notebooks

```bash
# Start Jupyter
jupyter notebook

# Or use JupyterLab
jupyter lab
```

Then navigate to the `notebooks/` directory and open any tutorial.

### VS Code Integration

If using VS Code:
1. Install the Jupyter extension
2. Open any `.ipynb` file
3. Select your Python kernel
4. Run cells interactively

---

## Tutorial Progression

```
01_basic_usage.ipynb
    ↓
    Learn fundamentals
    ↓
02_advanced_regularization.ipynb
    ↓
    Master parameters
    ↓
03_real_data_analysis.ipynb
    ↓
    Apply to real data
```

---

## Quick Reference

### LiVAE Basic Usage

```python
from livae import agent
import anndata as ad

# Load your data
adata = ad.read_h5ad('your_data.h5ad')

# Initialize model
model = agent(
    adata=adata,
    layer='counts',
    latent_dim=10,
    i_dim=2,
    hidden_dim=128,
    lr=1e-3
)

# Train
model.fit(epochs=100)

# Extract embeddings
latent = model.get_latent()           # High-dimensional
interpretable = model.get_iembed()    # Low-dimensional
```

### Key Parameters

| Parameter | Description | Typical Range |
|-----------|-------------|---------------|
| `latent_dim` | Primary latent space dimension | 10-30 |
| `i_dim` | Interpretable embedding dimension | 2-10 |
| `hidden_dim` | Hidden layer size | 64-256 |
| `percent` | Batch size (fraction of data) | 0.01-0.2 |
| `beta` | β-VAE disentanglement weight | 1.0-5.0 |
| `lorentz` | Hyperbolic geometry weight | 0.0-3.0 |
| `irecon` | Interpretable reconstruction weight | 0.0-2.0 |

---

## Dataset Requirements

### Minimum Requirements
- **Cells:** ≥100 (recommended: ≥500)
- **Genes:** ≥50 (recommended: ≥500)
- **Format:** AnnData object with raw counts

### Recommended Preprocessing
```python
import scanpy as sc

# Basic filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Store raw counts
adata.layers['counts'] = adata.X.copy()
```

---

## Common Issues and Solutions

### Issue: CUDA out of memory
**Solution:** Reduce `hidden_dim` or `percent` (batch size)
```python
model = agent(adata, hidden_dim=64, percent=0.05)
```

### Issue: Training is slow
**Solution:** Increase batch size, reduce epochs
```python
model = agent(adata, percent=0.2)
model.fit(epochs=50)
```

### Issue: Poor clustering quality
**Solution:** Try different regularization
```python
model = agent(adata, beta=2.0, lorentz=1.5, irecon=1.0)
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
