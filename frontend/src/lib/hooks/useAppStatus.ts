
// src/lib/hooks/useAppStatus.ts
import { useState, useEffect, useCallback } from 'react';
import { api } from '@/lib/api/endpoints';
import type { AppStatus } from '@/types/index';

interface EnhancedAppStatus extends AppStatus {
  last_updated: Date;
  is_loading: boolean;
  connection_error: boolean;
}

export const useAppStatus = (pollInterval: number = 5000) => {
  const [status, setStatus] = useState<EnhancedAppStatus>({
    data_loaded: false,
    model_trained: false,
    training_running: false,
    device: 'unknown',
    last_updated: new Date(),
    is_loading: true,
    connection_error: false,
  });

  const fetchStatus = useCallback(async () => {
    try {
      const [appStatus, trainingState] = await Promise.all([
        api.getStatus(),
        api.getTrainingProgress().catch(() => null), // Don't fail if training endpoint fails
      ]);

      // Enhanced logic: Only mark as trained if training is NOT running AND was previously completed
      const actuallyTrained = appStatus.model_trained && !trainingState?.is_running;

      setStatus(({
        ...appStatus,
        model_trained: actuallyTrained,
        training_running: trainingState?.is_running || appStatus.training_running,
        last_updated: new Date(),
        is_loading: false,
        connection_error: false,
      }));

    } catch (error) {
      console.error('Failed to fetch status:', error);
      setStatus(prev => ({
        ...prev,
        is_loading: false,
        connection_error: true,
        last_updated: new Date(),
      }));
    }
  }, []);

  // Initial fetch
  useEffect(() => {
    fetchStatus();
  }, [fetchStatus]);

  // Polling for status updates
  useEffect(() => {
    if (pollInterval > 0) {
      const interval = setInterval(fetchStatus, pollInterval);
      return () => clearInterval(interval);
    }
  }, [fetchStatus, pollInterval]);

  // Manual refresh function
  const refresh = useCallback(() => {
    setStatus(prev => ({ ...prev, is_loading: true }));
    fetchStatus();
  }, [fetchStatus]);

  return {
    ...status,
    refresh,
  };
};
