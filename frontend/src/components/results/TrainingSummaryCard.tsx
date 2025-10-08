// ========================================
// src/components/results/TrainingSummaryCard.tsx
// ========================================
import React from 'react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle, Badge, Progress } from './ui';
import type { TrainingState, AnnDataSummary } from '@/types/index';

interface TrainingSummaryCardProps {
  trainingState: TrainingState;
  dataSummary: AnnDataSummary;
}

export const TrainingSummaryCard: React.FC<TrainingSummaryCardProps> = ({
  trainingState,
  dataSummary,
}) => {
  const formatDuration = (startTime?: string, endTime?: string) => {
    if (!startTime) return 'Unknown';
    const start = new Date(startTime);
    const end = endTime ? new Date(endTime) : new Date();
    const diffMs = end.getTime() - start.getTime();
    const diffMins = Math.floor(diffMs / (1000 * 60));
    const diffSecs = Math.floor((diffMs % (1000 * 60)) / 1000);
    
    if (diffMins > 60) {
      const hours = Math.floor(diffMins / 60);
      const mins = diffMins % 60;
      return `${hours}h ${mins}m`;
    }
    if (diffMins > 0) {
      return `${diffMins}m ${diffSecs}s`;
    }
    return `${diffSecs}s`;
  };

  const getTrainingStatus = () => {
    if (trainingState.error_message) return 'error';
    if (trainingState.is_running) return 'running';
    if (trainingState.completed_at) return 'completed';
    return 'incomplete';
  };

  const status = getTrainingStatus();
  const progressPercentage = Math.round((trainingState.current_epoch / trainingState.total_epochs) * 100);

  return (
    <Card>
      <CardHeader>
        <CardTitle className="flex items-center justify-between">
          Training Summary
          <Badge variant={
            status === 'completed' ? 'default' :
            status === 'running' ? 'destructive' :
            status === 'error' ? 'destructive' : 'secondary'
          }>
            {status === 'completed' ? '✅ Completed' :
             status === 'running' ? '🏃‍♂️ Running' :
             status === 'error' ? '❌ Error' : '⏸️ Incomplete'}
          </Badge>
        </CardTitle>
        <CardDescription>
          Training session details and final metrics
        </CardDescription>
      </CardHeader>
      <CardContent className="space-y-6">
        {/* Training Progress */}
        <div className="space-y-2">
          <div className="flex justify-between text-sm">
            <span>Progress</span>
            <span>{trainingState.current_epoch} / {trainingState.total_epochs} epochs ({progressPercentage}%)</span>
          </div>
          <Progress value={progressPercentage} />
        </div>

        {/* Dataset Info */}
        <div className="grid grid-cols-2 md:grid-cols-4 gap-4 text-sm">
          <div className="text-center">
            <div className="text-lg font-bold text-blue-600">
              {dataSummary.shape.n_obs.toLocaleString()}
            </div>
            <div className="text-gray-600">Cells</div>
          </div>
          <div className="text-center">
            <div className="text-lg font-bold text-green-600">
              {dataSummary.shape.n_vars.toLocaleString()}
            </div>
            <div className="text-gray-600">Genes</div>
          </div>
          <div className="text-center">
            <div className="text-lg font-bold text-purple-600">
              {formatDuration(trainingState.started_at, trainingState.completed_at)}
            </div>
            <div className="text-gray-600">Duration</div>
          </div>
          <div className="text-center">
            <div className="text-lg font-bold text-orange-600">
              {trainingState.history.length}
            </div>
            <div className="text-gray-600">Data Points</div>
          </div>
        </div>

        {/* Final Metrics */}
        {trainingState.latest_metrics && (
          <div className="space-y-3">
            <h4 className="font-medium">Final Performance Metrics</h4>
            <div className="grid grid-cols-2 md:grid-cols-4 gap-4 text-sm">
              <div className="bg-gray-50 p-3 rounded">
                <div className="font-medium">Loss</div>
                <div className="text-lg text-red-600">{trainingState.latest_metrics.loss.toFixed(4)}</div>
              </div>
              <div className="bg-gray-50 p-3 rounded">
                <div className="font-medium">ARI</div>
                <div className="text-lg text-green-600">{trainingState.latest_metrics.ari.toFixed(3)}</div>
              </div>
              <div className="bg-gray-50 p-3 rounded">
                <div className="font-medium">NMI</div>
                <div className="text-lg text-blue-600">{trainingState.latest_metrics.nmi.toFixed(3)}</div>
              </div>
              <div className="bg-gray-50 p-3 rounded">
                <div className="font-medium">ASW</div>
                <div className="text-lg text-purple-600">{trainingState.latest_metrics.asw.toFixed(3)}</div>
              </div>
            </div>
          </div>
        )}

        {/* Error Message */}
        {trainingState.error_message && (
          <div className="bg-red-50 border border-red-200 rounded p-3">
            <div className="font-medium text-red-800">Training Error:</div>
            <div className="text-red-600 text-sm">{trainingState.error_message}</div>
          </div>
        )}

        {/* Training Times */}
        <div className="text-xs text-gray-500 space-y-1">
          <div>Started: {trainingState.started_at ? new Date(trainingState.started_at).toLocaleString() : 'Unknown'}</div>
          {trainingState.completed_at && (
            <div>Completed: {new Date(trainingState.completed_at).toLocaleString()}</div>
          )}
        </div>
      </CardContent>
    </Card>
  );
};
