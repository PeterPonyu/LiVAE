
// src/lib/api/endpoints.ts
import { apiClient } from './client';
import type {
  AppStatus,
  UploadResponse,
  AnnDataSummary,
  TrainingResponse,
  TrainingState,
  AgentParameters,
  TrainingConfig,
  EmbeddingResult,
  EmbeddingType
} from '@/types/index';

export const api = {
  // System status and health
  getStatus: () => 
    apiClient.get<AppStatus>('/status'),

  getRoot: () => 
    apiClient.get<{message: string; version: string; docs: string; status: string}>('/'),

  healthCheck: (timeoutMs?: number) =>
    apiClient.healthCheck(timeoutMs),

  // Data management
  uploadData: (file: File, onProgress?: (progress: number) => void) => 
    apiClient.uploadFile<UploadResponse>('/upload-data', file, onProgress),

  getDataSummary: () => 
    apiClient.get<AnnDataSummary>('/data-summary'),

  // Training management
  startTraining: (parameters: AgentParameters, config: TrainingConfig) =>
    apiClient.post<TrainingResponse>('/start-training', {
      parameters,
      config,
    }),

  getTrainingProgress: () =>
    apiClient.get<TrainingState>('/training-progress'),

  stopTraining: () =>
    apiClient.post<{success: boolean; message: string}>('/stop-training'),

  // Results & embeddings
  getInterpretableEmbeddings: () =>
    apiClient.get<EmbeddingResult>('/embeddings/interpretable'),

  getLatentEmbeddings: () =>
    apiClient.get<EmbeddingResult>('/embeddings/latent'),

  // Get embedding info (without downloading)
  getEmbeddingInfo: (embeddingType: EmbeddingType) =>
    apiClient.get<EmbeddingResult>(`/embeddings/${embeddingType}`),

  // Extract/generate embeddings from trained model
  extractEmbeddings: () =>
    apiClient.post<{success: boolean; message: string}>('/extract-embeddings'),

  // Download embeddings as blobs with progress
  downloadEmbedding: (embeddingType: EmbeddingType, onProgress?: (progress: number) => void): Promise<Blob> =>
    apiClient.downloadBlob(`/download/embeddings/${embeddingType}`, onProgress),

  // Get download URLs (fallback method)
  getEmbeddingDownloadUrl: (embeddingType: EmbeddingType): string =>
    apiClient.getDownloadUrl(`/download/embeddings/${embeddingType}`),

  // Additional result endpoints
  downloadTrainingHistory: (onProgress?: (progress: number) => void): Promise<Blob> =>
    apiClient.downloadBlob('/download/training-history', onProgress),

  downloadTrainingSummary: (onProgress?: (progress: number) => void): Promise<Blob> =>
    apiClient.downloadBlob('/download/training-summary', onProgress),

  // Export all results as ZIP
  downloadAllResults: (onProgress?: (progress: number) => void): Promise<Blob> =>
    apiClient.downloadBlob('/download/all-results', onProgress),
};

// Helper functions for common operations
export const apiHelpers = {
  // Check if backend is healthy
  isHealthy: async (timeoutMs: number = 5000): Promise<boolean> => {
    try {
      return await api.healthCheck(timeoutMs);
    } catch {
      return false;
    }
  },

  // Poll training progress (for real-time updates)
  pollTrainingProgress: (
    onUpdate: (state: TrainingState) => void,
    onError: (error: Error) => void,
    intervalMs: number = 2000
  ): (() => void) => {
    let isPolling = true;
    
    const poll = async () => {
      try {
        if (!isPolling) return;
        
        const state = await api.getTrainingProgress();
        onUpdate(state);
        
        // Continue polling if training is running
        if (state.is_running && isPolling) {
          setTimeout(poll, intervalMs);
        }
      } catch (error) {
        if (isPolling) {
          onError(error as Error);
        }
      }
    };
    
    poll();
    
    // Return stop function
    return () => {
      isPolling = false;
    };
  },

  // Upload with progress tracking
  uploadWithProgress: async (
    file: File,
    onProgress?: (progress: number) => void
  ): Promise<UploadResponse> => {
    return api.uploadData(file, onProgress);
  },

  // Download file with automatic filename
  downloadFile: async (
    blob: Blob,
    filename: string,
    onComplete?: () => void
  ): Promise<void> => {
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    window.URL.revokeObjectURL(url);
    
    if (onComplete) {
      onComplete();
    }
  },

  // Get all available embeddings
  getAllEmbeddings: async (): Promise<EmbeddingResult[]> => {
    const embeddings: EmbeddingResult[] = [];
    
    try {
      const interpretable = await api.getInterpretableEmbeddings();
      embeddings.push(interpretable);
    } catch {
      // Interpretable embedding not available
    }
    
    try {
      const latent = await api.getLatentEmbeddings();
      embeddings.push(latent);
    } catch {
      // Latent embedding not available
    }
    
    return embeddings;
  },

  // Check training completion status
  isTrainingComplete: async (): Promise<boolean> => {
    try {
      const state = await api.getTrainingProgress();
      return !state.is_running && state.current_epoch > 0 && !state.error_message;
    } catch {
      return false;
    }
  },

  // Validate training readiness
  isReadyForTraining: async (): Promise<{ready: boolean; reason?: string}> => {
    try {
      // Check if data is loaded
      const summary = await api.getDataSummary();
      if (!summary || summary.shape.n_obs === 0) {
        return { ready: false, reason: 'No data loaded' };
      }

      // Check if training is already running
      const state = await api.getTrainingProgress();
      if (state.is_running) {
        return { ready: false, reason: 'Training already in progress' };
      }

      return { ready: true };
    } catch {
      return { ready: false, reason: 'Unable to verify readiness' };
    }
  },
};

// Export additional utilities
export { apiClient } from './client';
export type { ApiError } from './client';
