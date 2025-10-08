
// ========================================
// src/components/monitor/TrainingStatusCard.tsx
// ========================================
'use client';

import React, { useMemo } from 'react';
import { Square, TrendingUp, AlertCircle } from 'lucide-react';
import Link from 'next/link';
import type { TrainingState } from '@/types/index';

interface TrainingStatusCardProps {
  trainingState: TrainingState;
  onStop: () => void;
}

export const TrainingStatusCard: React.FC<TrainingStatusCardProps> = ({
  trainingState,
  onStop,
}) => {
  const progressPercentage = useMemo(() => {
    return Math.round((trainingState.current_epoch / trainingState.total_epochs) * 100);
  }, [trainingState.current_epoch, trainingState.total_epochs]);

  const formatDuration = (startTime?: string) => {
    if (!startTime) return 'Unknown';
    const start = new Date(startTime);
    const now = new Date();
    const diffMs = now.getTime() - start.getTime();
    const diffMins = Math.floor(diffMs / (1000 * 60));
    const diffSecs = Math.floor((diffMs % (1000 * 60)) / 1000);
    
    if (diffMins > 0) return `${diffMins}m ${diffSecs}s`;
    return `${diffSecs}s`;
  };

  return (
    <div className="rounded-lg border border-gray-200 bg-white shadow-sm">
      {/* Card Header */}
      <div className="flex flex-col space-y-1.5 p-6 pb-3">
        <div className="flex items-center justify-between">
          <h3 className="text-2xl font-semibold leading-none tracking-tight">
            Training Status
          </h3>
          <span
            className={`inline-flex items-center rounded-full px-2.5 py-0.5 text-xs font-semibold ${
              trainingState.is_running
                ? 'bg-red-100 text-red-800'
                : 'bg-gray-100 text-gray-800'
            }`}
          >
            {trainingState.is_running ? 'Running' : 'Stopped'}
          </span>
        </div>
      </div>

      {/* Card Content */}
      <div className="p-6 pt-0 space-y-4">
        {/* Progress Section */}
        <div>
          <div className="flex justify-between text-sm mb-2 font-medium">
            <span>Epoch {trainingState.current_epoch} / {trainingState.total_epochs}</span>
            <span>{progressPercentage}%</span>
          </div>
          <div className="w-full bg-gray-200 rounded-full h-2 overflow-hidden">
            <div
              className="bg-blue-600 h-2 rounded-full transition-all duration-300"
              style={{ width: `${progressPercentage}%` }}
            />
          </div>
        </div>
        
        {/* Duration */}
        <div className="flex justify-between text-sm text-gray-600">
          <span>Duration</span>
          <span className="font-medium">{formatDuration(trainingState.started_at)}</span>
        </div>
        
        {/* Action Buttons */}
        <div className="pt-2 flex gap-2">
          {trainingState.is_running ? (
            <button
              onClick={onStop}
              className="w-full inline-flex items-center justify-center rounded-md text-sm font-medium bg-red-600 text-white hover:bg-red-700 h-9 px-3 transition-colors focus:outline-none focus:ring-2 focus:ring-red-500 focus:ring-offset-2"
            >
              <Square className="mr-2 h-4 w-4" />
              Stop Training
            </button>
          ) : (
            <>
              <Link
                href="/results"
                className="flex-1 inline-flex items-center justify-center rounded-md text-sm font-medium bg-blue-600 text-white hover:bg-blue-700 h-9 px-3 transition-colors focus:outline-none focus:ring-2 focus:ring-blue-500 focus:ring-offset-2"
              >
                <TrendingUp className="mr-2 h-4 w-4" />
                View Results
              </Link>
              <Link
                href="/training/configure"
                className="flex-1 inline-flex items-center justify-center rounded-md text-sm font-medium border border-gray-300 bg-white text-gray-700 hover:bg-gray-50 h-9 px-3 transition-colors focus:outline-none focus:ring-2 focus:ring-blue-500 focus:ring-offset-2"
              >
                New Training
              </Link>
            </>
          )}
        </div>

        {/* Error Alert */}
        {trainingState.error_message && (
          <div className="mt-2 rounded-lg border border-red-200 bg-red-50 p-4">
            <div className="flex items-start gap-3">
              <AlertCircle className="h-4 w-4 text-red-600 mt-0.5 flex-shrink-0" />
              <p className="text-xs text-red-800">
                {trainingState.error_message}
              </p>
            </div>
          </div>
        )}
      </div>
    </div>
  );
};
