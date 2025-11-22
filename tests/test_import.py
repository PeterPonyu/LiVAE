import numpy as np
import anndata as ad
import torch

def test_import_agent():
    from livae import agent
    X = np.random.poisson(1.0, (20, 50)).astype(float)
    adata = ad.AnnData(X)
    adata.layers['counts'] = X
    a = agent(adata=adata, latent_dim=5, i_dim=2, hidden_dim=32, percent=0.1, lr=1e-3)
    a.fit(epochs=2)
    latent = a.get_latent()
    iembed = a.get_iembed()
    assert latent.shape[0] == 20
    assert iembed.shape[0] == 20
    assert iembed.shape[1] == 2
