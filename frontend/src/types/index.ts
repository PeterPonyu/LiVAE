// Start with just what you need immediately
export interface UploadResponse {
  success: boolean;
  message: string;
  data_summary?: AnnDataSummary;
}

export interface AnnDataSummary {
  shape: {
    n_obs: number;
    n_vars: number;
  };
  filename?: string;
  file_size_mb?: number;
  loaded_at: string;
  // Add more fields as you need them
}

export interface AppStatus {
  data_loaded: boolean;
  model_trained: boolean;
  training_running: boolean;
  device: string;
}