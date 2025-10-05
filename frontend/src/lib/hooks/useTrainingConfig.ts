
// src/lib/hooks/useTrainingConfig.ts
import { useState, useCallback } from 'react';
import { api } from '@/lib/api/endpoints';
import type { TrainingResponse, AgentParameters, TrainingConfig } from '@/types/index';
import type { TrainingFormData } from '@/lib/validation/training-schemas';

interface TrainingConfigState {
  isSubmitting: boolean;
  error: string | null;
  isSuccess: boolean;
}

export const useTrainingConfig = () => {
  const [state, setState] = useState<TrainingConfigState>({
    isSubmitting: false,
    error: null,
    isSuccess: false,
  });

  const submitTrainingConfig = useCallback(async (formData: TrainingFormData) => {
    setState({
      isSubmitting: true,
      error: null,
      isSuccess: false,
    });

    try {
      // Transform form data to API format
      const agentParameters: AgentParameters = formData.agent_parameters;
      
      const trainingConfig: TrainingConfig = {
        epochs: formData.training_config.epochs,
        qcparams: formData.training_config.apply_qc 
          ? formData.training_config.qc_params
          : undefined,
      };

      const response: TrainingResponse = await api.startTraining(agentParameters, trainingConfig);

      if (response.success) {
        setState({
          isSubmitting: false,
          error: null,
          isSuccess: true,
        });
      } else {
        setState({
          isSubmitting: false,
          error: response.message || 'Failed to start training',
          isSuccess: false,
        });
      }
    } catch (error) {
      setState({
        isSubmitting: false,
        error: error instanceof Error ? error.message : 'Failed to start training',
        isSuccess: false,
      });
    }
  }, []);

  const reset = useCallback(() => {
    setState({
      isSubmitting: false,
      error: null,
      isSuccess: false,
    });
  }, []);

  return {
    ...state,
    submitTrainingConfig,
    reset,
  };
};
