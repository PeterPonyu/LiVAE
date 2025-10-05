
// src/lib/hooks/useDataUpload.ts
import { useState, useCallback } from 'react';
import { api } from '@/lib/api/endpoints';
import type { UploadResponse, AnnDataSummary } from '@/types/index';

interface UploadState {
  isUploading: boolean;
  progress: number;
  error: string | null;
  dataSummary: AnnDataSummary | null;
  isSuccess: boolean;
}

export const useDataUpload = () => {
  const [state, setState] = useState<UploadState>({
    isUploading: false,
    progress: 0,
    error: null,
    dataSummary: null,
    isSuccess: false,
  });

  const uploadFile = useCallback(async (file: File) => {
    // Reset state
    setState({
      isUploading: true,
      progress: 0,
      error: null,
      dataSummary: null,
      isSuccess: false,
    });

    try {
      // Simulate progress for better UX (real progress would need XMLHttpRequest)
      setState(prev => ({ ...prev, progress: 10 }));
      
      const response: UploadResponse = await api.uploadData(file);
      
      setState(prev => ({ ...prev, progress: 90 }));
      
      if (response.success && response.data_summary) {
        setState({
          isUploading: false,
          progress: 100,
          error: null,
          dataSummary: response.data_summary,
          isSuccess: true,
        });
      } else {
        setState({
          isUploading: false,
          progress: 0,
          error: response.message || 'Upload failed',
          dataSummary: null,
          isSuccess: false,
        });
      }
    } catch (error) {
      setState({
        isUploading: false,
        progress: 0,
        error: error instanceof Error ? error.message : 'Upload failed',
        dataSummary: null,
        isSuccess: false,
      });
    }
  }, []);

  const reset = useCallback(() => {
    setState({
      isUploading: false,
      progress: 0,
      error: null,
      dataSummary: null,
      isSuccess: false,
    });
  }, []);

  return {
    ...state,
    uploadFile,
    reset,
  };
};
