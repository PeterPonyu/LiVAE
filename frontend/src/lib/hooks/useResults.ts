
// src/lib/hooks/useResults.ts
import { useState, useEffect, useCallback } from 'react';
import { api, apiHelpers } from '@/lib/api/endpoints';
import type { EmbeddingResult, TrainingState, AnnDataSummary, EmbeddingType } from '@/types/index';

interface ResultsState {
  isLoading: boolean;
  error: string | null;
  trainingCompleted: boolean;
  trainingState: TrainingState | null;
  dataSummary: AnnDataSummary | null;
  availableEmbeddings: EmbeddingResult[];
  isGeneratingEmbeddings: boolean;
  downloadProgress: Record<string, number>;
  downloadStatus: Record<string, 'idle' | 'downloading' | 'complete' | 'error'>;
}

export const useResults = () => {
  const [state, setState] = useState<ResultsState>({
    isLoading: true,
    error: null,
    trainingCompleted: false,
    trainingState: null,
    dataSummary: null,
    availableEmbeddings: [],
    isGeneratingEmbeddings: false,
    downloadProgress: {},
    downloadStatus: {},
  });

  const fetchResults = useCallback(async () => {
    try {
      setState(prev => ({ ...prev, isLoading: true, error: null }));

      // Fetch all necessary data
      const [trainingState, dataSummary] = await Promise.all([
        api.getTrainingProgress(),
        api.getDataSummary(),
      ]);

      // Check if training is actually completed
      const trainingCompleted = await apiHelpers.isTrainingComplete();

      // Get available embeddings
      let availableEmbeddings: EmbeddingResult[] = [];
      
      if (trainingCompleted) {
        availableEmbeddings = await apiHelpers.getAllEmbeddings();
      }

      setState(prev => ({
        ...prev,
        isLoading: false,
        trainingCompleted,
        trainingState,
        dataSummary,
        availableEmbeddings,
      }));

    } catch (error) {
      setState(prev => ({
        ...prev,
        isLoading: false,
        error: error instanceof Error ? error.message : 'Failed to fetch results',
      }));
    }
  }, []);

  const downloadEmbedding = useCallback(async (embeddingType: EmbeddingType) => {
    try {
      setState(prev => ({
        ...prev,
        downloadStatus: { ...prev.downloadStatus, [embeddingType]: 'downloading' },
        downloadProgress: { ...prev.downloadProgress, [embeddingType]: 0 },
      }));

      const blob = await api.downloadEmbedding(embeddingType, (progress) => {
        setState(prev => ({
          ...prev,
          downloadProgress: { ...prev.downloadProgress, [embeddingType]: progress },
        }));
      });
      
      // Download the file
      await apiHelpers.downloadFile(blob, `${embeddingType}_embedding.csv`, () => {
        setState(prev => ({
          ...prev,
          downloadStatus: { ...prev.downloadStatus, [embeddingType]: 'complete' },
          downloadProgress: { ...prev.downloadProgress, [embeddingType]: 100 },
        }));

        // Reset status after 3 seconds
        setTimeout(() => {
          setState(prev => ({
            ...prev,
            downloadStatus: { ...prev.downloadStatus, [embeddingType]: 'idle' },
            downloadProgress: { ...prev.downloadProgress, [embeddingType]: 0 },
          }));
        }, 3000);
      });

    } catch (error) {
      setState(prev => ({
        ...prev,
        downloadStatus: { ...prev.downloadStatus, [embeddingType]: 'error' },
        downloadProgress: { ...prev.downloadProgress, [embeddingType]: 0 },
        error: error instanceof Error ? error.message : 'Download failed',
      }));
    }
  }, []);

  const generateEmbeddings = useCallback(async () => {
    try {
      setState(prev => ({ ...prev, isGeneratingEmbeddings: true, error: null }));
      
      await api.extractEmbeddings();
      
      // Refresh results after generating embeddings
      await fetchResults();
      
    } catch (error) {
      setState(prev => ({
        ...prev,
        error: error instanceof Error ? error.message : 'Failed to generate embeddings',
      }));
    } finally {
      setState(prev => ({ ...prev, isGeneratingEmbeddings: false }));
    }
  }, [fetchResults]);

  const downloadAllResults = useCallback(async () => {
    try {
      setState(prev => ({
        ...prev,
        downloadStatus: { ...prev.downloadStatus, 'all': 'downloading' },
        downloadProgress: { ...prev.downloadProgress, 'all': 0 },
      }));

      const blob = await api.downloadAllResults((progress) => {
        setState(prev => ({
          ...prev,
          downloadProgress: { ...prev.downloadProgress, 'all': progress },
        }));
      });

      await apiHelpers.downloadFile(blob, 'training_results.zip', () => {
        setState(prev => ({
          ...prev,
          downloadStatus: { ...prev.downloadStatus, 'all': 'complete' },
          downloadProgress: { ...prev.downloadProgress, 'all': 100 },
        }));

        setTimeout(() => {
          setState(prev => ({
            ...prev,
            downloadStatus: { ...prev.downloadStatus, 'all': 'idle' },
            downloadProgress: { ...prev.downloadProgress, 'all': 0 },
          }));
        }, 3000);
      });

    } catch (error) {
      setState(prev => ({
        ...prev,
        downloadStatus: { ...prev.downloadStatus, 'all': 'error' },
        downloadProgress: { ...prev.downloadProgress, 'all': 0 },
        error: error instanceof Error ? error.message : 'Download failed',
      }));
    }
  }, []);

  // Initial fetch
  useEffect(() => {
    fetchResults();
  }, [fetchResults]);

  return {
    ...state,
    downloadEmbedding,
    generateEmbeddings,
    downloadAllResults,
    refresh: fetchResults,
    clearError: () => setState(prev => ({ ...prev, error: null })),
  };
};
