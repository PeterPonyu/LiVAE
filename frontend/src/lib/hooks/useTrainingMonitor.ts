
// src/lib/hooks/useTrainingMonitor.ts
import { useState, useEffect, useCallback, useRef } from 'react';
import { api } from '@/lib/api/endpoints';
import type { TrainingState } from '@/types/index';

export const useTrainingMonitor = (pollInterval: number = 2000) => {
  const [trainingState, setTrainingState] = useState<TrainingState | null>(null);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const pollIntervalRef = useRef<NodeJS.Timeout | null>(null);

  const fetchTrainingProgress = useCallback(async () => {
    try {
      const progress = await api.getTrainingProgress();
      setTrainingState(progress);
      setError(null);
    } catch (err) {
      const errorMessage = err instanceof Error ? err.message : 'Failed to fetch training progress';
      setError(errorMessage);
      console.error('Training progress fetch error:', err);
    } finally {
      setIsLoading(false);
    }
  }, []);

  // Start polling
  const startPolling = useCallback(() => {
    if (pollIntervalRef.current) {
      clearInterval(pollIntervalRef.current);
    }
    
    // Initial fetch
    fetchTrainingProgress();
    
    // Set up polling
    pollIntervalRef.current = setInterval(() => {
      fetchTrainingProgress();
    }, pollInterval);
  }, [fetchTrainingProgress, pollInterval]);

  // Stop polling
  const stopPolling = useCallback(() => {
    if (pollIntervalRef.current) {
      clearInterval(pollIntervalRef.current);
      pollIntervalRef.current = null;
    }
  }, []);

  // Auto-start polling when component mounts
  useEffect(() => {
    startPolling();
    return () => stopPolling();
  }, [startPolling, stopPolling]);

  // Auto-stop polling when training completes
  useEffect(() => {
    if (trainingState && !trainingState.is_running) {
      // Continue polling for a bit after completion to get final metrics
      setTimeout(() => {
        if (trainingState && !trainingState.is_running) {
          stopPolling();
        }
      }, 10000); // Stop polling 10 seconds after completion
    }
  }, [trainingState, stopPolling]);

  // Stop training function
  const stopTraining = useCallback(async () => {
    try {
      await api.stopTraining();
      // Refresh immediately after stopping
      setTimeout(fetchTrainingProgress, 1000);
    } catch (err) {
      const errorMessage = err instanceof Error ? err.message : 'Failed to stop training';
      setError(errorMessage);
    }
  }, [fetchTrainingProgress]);

  return {
    trainingState,
    isLoading,
    error,
    stopTraining,
    refresh: fetchTrainingProgress,
    startPolling,
    stopPolling,
  };
};
