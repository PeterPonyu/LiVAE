
// src/types/index.ts
// Complete TypeScript mirror of your FastAPI Pydantic models

// Basic info models
export interface ShapeInfo {
  n_obs: number;
  n_vars: number;
}

export interface ColumnInfo {
  name: string;
  dtype: string;
  n_unique: number;
  missing_count: number;
  sample_values?: (string | number | null)[];
  is_truncated: boolean;
}

export interface LayerInfo {
  name: string;
  dtype: string;
  shape: ShapeInfo;
  sparsity_ratio?: number;
  min_val?: number;
  max_val?: number;
  mean_val?: number;
}

export interface QCMetrics {
  n_genes_by_counts: Record<string, number>;
  total_counts: Record<string, number>;
  pct_counts_mt: Record<string, number>;
  pct_counts_ribo: Record<string, number>;
}

export interface AnnDataSummary {
  shape: ShapeInfo;
  X_info: LayerInfo;
  layers: Record<string, LayerInfo>;
  obs_columns: ColumnInfo[];
  var_columns: ColumnInfo[];
  qc_metrics?: QCMetrics;
  filename?: string;
  file_size_mb?: number;
  loaded_at: string;
}

// Fix: Use literal types for species
export interface QCParams {
  species: 'human' | 'mouse';  // Fixed: use literal types instead of string
  min_genes_per_cell?: number;
  max_genes_per_cell?: number;
  min_counts_per_cell?: number;
  max_counts_per_cell?: number;
  max_pct_mt?: number;
  max_pct_ribo?: number;
}

export interface AgentParameters {
  layer: string;
  percent: number;
  irecon: number;
  lorentz: number;
  beta: number;
  dip: number;
  tc: number;
  info: number;
  hidden_dim: number;
  latent_dim: number;
  i_dim: number;
  lr: number;
}

export interface TrainingConfig {
  epochs: number;
  qcparams?: QCParams;
}

export interface TrainingMetrics {
  epoch: number;
  loss: number;
  ari: number;   // Adjusted Rand Index
  nmi: number;   // Normalized Mutual Information  
  asw: number;   // Average Silhouette Width
  ch: number;    // Calinski-Harabasz Index
  db: number;    // Davies-Bouldin Index
  pc: number;    // Graph connectivity score
}

export interface TrainingState {
  is_running: boolean;
  current_epoch: number;
  total_epochs: number;
  latest_metrics?: TrainingMetrics;
  history: TrainingMetrics[];
  error_message?: string;
  started_at?: string;
  completed_at?: string;
}

export interface EmbeddingResult {
  embedding_type: string;  // "interpretable" or "latent"
  shape: ShapeInfo;
  download_ready: boolean;
  csv_path?: string;
  extracted_at?: string;
}

// Response models
export interface UploadResponse {
  success: boolean;
  message: string;
  data_summary?: AnnDataSummary;
}

export interface TrainingResponse {
  success: boolean;
  message: string;
}

export interface AppStatus {
  data_loaded: boolean;
  model_trained: boolean;
  training_running: boolean;
  device: string;
}

// Additional utility types
export type EmbeddingType = 'interpretable' | 'latent';
export type TrainingStatus = 'idle' | 'running' | 'completed' | 'failed';
export type Species = 'human' | 'mouse';

// Form types (for React Hook Form) - Fix: ensure consistency
export interface QCFilterFormData {
  species: Species;
  min_genes_per_cell?: number;
  max_genes_per_cell?: number;
  min_counts_per_cell?: number;
  max_counts_per_cell?: number;
  max_pct_mt?: number;
  max_pct_ribo?: number;
}

export interface AgentParametersFormData extends AgentParameters {
  // Any additional form-specific fields
}

export interface TrainingConfigFormData {
  epochs: number;
  apply_qc: boolean;
  qc_filters?: QCFilterFormData;
}

// API Error Response
export interface ApiErrorResponse {
  detail: string;
  status_code?: number;
}

// Upload file constraints
export interface UploadConstraints {
  maxSizeBytes: number;
  allowedTypes: string[];
}

export const UPLOAD_CONSTRAINTS: UploadConstraints = {
  maxSizeBytes: 100 * 1024 * 1024, // 100MB
  allowedTypes: ['.h5ad'],
};

// Default values matching backend defaults
export const DEFAULT_AGENT_PARAMETERS: AgentParameters = {
  layer: 'X',  // Fixed: use 'X' instead of 'counts' as default
  percent: 0.01,
  irecon: 0.0,
  lorentz: 0.0,
  beta: 1.0,
  dip: 0.0,
  tc: 0.0,
  info: 0.0,
  hidden_dim: 128,
  latent_dim: 10,
  i_dim: 2,
  lr: 1e-4,
};

export const DEFAULT_TRAINING_CONFIG: TrainingConfig = {
  epochs: 1000,
};

export const DEFAULT_QC_PARAMS: QCParams = {
  species: 'human',
};
