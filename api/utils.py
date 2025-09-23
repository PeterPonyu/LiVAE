
import scanpy as sc
import pandas as pd
import numpy as np
from anndata import AnnData
from typing import Dict, List, Optional, Any, Union
import os
from datetime import datetime

from .models import (
    AnnDataSummary, ShapeInfo, LayerInfo, ColumnInfo, 
    QCMetrics, QCThresholds
)

def load_anndata_from_path(file_path: str) -> AnnData:
    """Load AnnData from file path"""
    try:
        adata = sc.read_h5ad(file_path)
        return adata
    except Exception as e:
        raise ValueError(f"Failed to load AnnData file: {str(e)}")

def calculate_qc_metrics(adata: AnnData) -> QCMetrics:
    """Calculate scanpy-style QC metrics"""
    # Create a copy to avoid modifying original
    adata_temp = adata.copy()
    
    # Calculate basic metrics
    adata_temp.var['mt'] = adata_temp.var_names.str.startswith('MT-')
    adata_temp.var['ribo'] = adata_temp.var_names.str.startswith(('RPS', 'RPL'))
    
    # Calculate QC metrics using scanpy
    sc.pp.calculate_qc_metrics(
        adata_temp, 
        percent_top=None, 
        log1p=False, 
        inplace=True,
        var_type='genes'
    )
    
    # Extract mitochondrial and ribosomal percentages
    sc.pp.calculate_qc_metrics(
        adata_temp,
        qc_vars=['mt'], 
        percent_top=None,
        log1p=False,
        inplace=True
    )
    
    sc.pp.calculate_qc_metrics(
        adata_temp,
        qc_vars=['ribo'], 
        percent_top=None,
        log1p=False,
        inplace=True  
    )
    
    def _get_stats(values: np.ndarray) -> Dict[str, float]:
        """Get basic statistics for a metric"""
        return {
            "mean": float(np.mean(values)),
            "std": float(np.std(values)),
            "min": float(np.min(values)),
            "max": float(np.max(values)),
            "median": float(np.median(values))
        }
    
    return QCMetrics(
        n_genes_by_counts=_get_stats(adata_temp.obs['n_genes_by_counts']),
        total_counts=_get_stats(adata_temp.obs['total_counts']),
        pct_counts_mt=_get_stats(adata_temp.obs['pct_counts_mt']) if 'pct_counts_mt' in adata_temp.obs else _get_stats(np.zeros(adata_temp.n_obs)),
        pct_counts_ribo=_get_stats(adata_temp.obs['pct_counts_ribo']) if 'pct_counts_ribo' in adata_temp.obs else _get_stats(np.zeros(adata_temp.n_obs))
    )

def apply_qc_filtering(adata: AnnData, thresholds: QCThresholds) -> AnnData:
    """Apply QC filtering based on thresholds"""
    if thresholds is None:
        return adata
        
    # Create a copy for filtering
    adata_filtered = adata.copy()
    
    # Calculate QC metrics if not present
    if 'n_genes_by_counts' not in adata_filtered.obs.columns:
        calculate_qc_metrics(adata_filtered)  # This will add the metrics to obs
    
    # Apply filters
    sc.pp.filter_cells(adata_filtered, min_genes=thresholds.min_genes_per_cell or 0)
    if thresholds.max_genes_per_cell:
        sc.pp.filter_cells(adata_filtered, max_genes=thresholds.max_genes_per_cell)
    
    # Additional filtering based on counts and percentages would go here
    # This is a simplified version - you might want more sophisticated filtering
    
    return adata_filtered

def get_layer_info(data, layer_name: str, shape: ShapeInfo) -> LayerInfo:
    """Extract information about a data layer"""
    if hasattr(data, 'toarray'):  # Sparse matrix
        data_dense = data.toarray()
        sparsity = float((data.nnz / (shape.n_obs * shape.n_vars)))
        sparsity_ratio = 1.0 - sparsity
    else:  # Dense matrix
        data_dense = data
        sparsity_ratio = float(np.sum(data_dense == 0) / data_dense.size)
    
    return LayerInfo(
        name=layer_name,
        dtype=str(data.dtype),
        shape=shape,
        sparsity_ratio=sparsity_ratio,
        min_val=float(np.min(data_dense)),
        max_val=float(np.max(data_dense)),
        mean_val=float(np.mean(data_dense))
    )

def get_column_info(series: pd.Series, max_display: int = 3) -> ColumnInfo:
    """Extract information about a DataFrame column"""
    n_unique = series.nunique()
    missing_count = series.isnull().sum()
    
    # Get sample values (first few unique values)
    unique_vals = series.dropna().unique()
    
    # Convert to basic Python types for JSON serialization
    sample_values = []
    is_truncated = False
    
    if len(unique_vals) > 0:
        # Take first few values for display
        display_count = min(max_display, len(unique_vals))
        sample_values = [_convert_to_json_serializable(val) for val in unique_vals[:display_count]]
        is_truncated = len(unique_vals) > max_display
    
    return ColumnInfo(
        name=series.name,
        dtype=str(series.dtype),
        n_unique=int(n_unique),
        missing_count=int(missing_count),
        sample_values=sample_values,
        is_truncated=is_truncated
    )

def _convert_to_json_serializable(val: Any) -> Union[str, float, int]:
    """Convert various data types to JSON-serializable types"""
    if pd.isna(val):
        return None
    elif isinstance(val, (np.integer, int)):
        return int(val)
    elif isinstance(val, (np.floating, float)):
        return float(val)
    else:
        return str(val)

def create_anndata_summary(adata: AnnData, filename: str = None, file_size_bytes: int = None) -> AnnDataSummary:
    """Create comprehensive AnnData summary"""
    shape = ShapeInfo(n_obs=adata.n_obs, n_vars=adata.n_vars)
    
    # Main X layer info
    X_info = get_layer_info(adata.X, "X", shape)
    
    # Additional layers
    layers = {}
    if adata.layers:
        for layer_name, layer_data in adata.layers.items():
            layers[layer_name] = get_layer_info(layer_data, layer_name, shape)
    
    # Observation columns (.obs)
    obs_columns = [get_column_info(adata.obs[col]) for col in adata.obs.columns]
    
    # Variable columns (.var)
    var_columns = [get_column_info(adata.var[col]) for col in adata.var.columns]
    
    # Calculate QC metrics
    qc_metrics = calculate_qc_metrics(adata)
    
    return AnnDataSummary(
        shape=shape,
        X_info=X_info,
        layers=layers,
        obs_columns=obs_columns,
        var_columns=var_columns,
        qc_metrics=qc_metrics,
        filename=filename,
        file_size_mb=file_size_bytes / 1024 / 1024 if file_size_bytes else None
    )

def save_embedding_as_csv(embedding: np.ndarray, embedding_type: str, output_dir: str = "downloads") -> str:
    """Save embedding matrix as CSV file"""
    os.makedirs(output_dir, exist_ok=True)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"{embedding_type}_embedding_{timestamp}.csv"
    filepath = os.path.join(output_dir, filename)
    
    # Convert to DataFrame with proper column names
    if embedding_type == "interpretable":
        columns = [f"i_dim_{i+1}" for i in range(embedding.shape[1])]
    else:  # latent
        columns = [f"latent_{i+1}" for i in range(embedding.shape[1])]
    
    df = pd.DataFrame(embedding, columns=columns)
    df.to_csv(filepath, index=False)
    
    return filepath
