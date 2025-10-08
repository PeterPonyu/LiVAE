
// ========================================
// src/components/monitor/TrainingMonitor.tsx
// ========================================
'use client';

import React, { useMemo } from 'react';
import { ArrowLeft, RefreshCw, Activity, AlertCircle } from 'lucide-react';
import Link from 'next/link';
import { useTrainingMonitor } from '@/lib/hooks/useTrainingMonitor';
import { TrainingStatusCard } from './TrainingStatusCard';
import { CurrentMetricsCard } from './CurrentMetricsCard';
import { MetricsChartsSection } from './MetricsChartsSection';

export const TrainingMonitor: React.FC = () => {
  const { 
    trainingState, 
    isLoading, 
    error, 
    stopTraining, 
    refresh 
  } = useTrainingMonitor();

  // Prepare chart data
  const chartData = useMemo(() => {
    if (!trainingState?.history) return [];
    return trainingState.history.map(metrics => ({
      epoch: metrics.epoch,
      loss: metrics.loss,
      ari: metrics.ari,
      nmi: metrics.nmi,
      asw: metrics.asw,
      ch: metrics.ch,
      db: metrics.db,
      pc: metrics.pc,
    }));
  }, [trainingState?.history]);

  // Loading state
  if (isLoading && !trainingState) {
    return (
      <div className="container mx-auto py-8 px-4 max-w-6xl">
        <div className="text-center">
          <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-blue-600 mx-auto mb-4"></div>
          <p>Loading training status...</p>
        </div>
      </div>
    );
  }

  // Error state
  if (error && !trainingState) {
    return (
      <div className="container mx-auto py-8 px-4 max-w-4xl">
        <div className="rounded-lg border border-red-200 bg-red-50 p-4">
          <div className="flex items-start gap-3">
            <AlertCircle className="h-5 w-5 text-red-600 mt-0.5 flex-shrink-0" />
            <p className="text-sm text-red-800">{error}</p>
          </div>
        </div>
        <div className="mt-4 flex gap-4">
          <button
            onClick={refresh}
            className="inline-flex items-center justify-center rounded-md text-sm font-medium bg-blue-600 text-white hover:bg-blue-700 h-10 px-4 py-2 transition-colors focus:outline-none focus:ring-2 focus:ring-blue-500 focus:ring-offset-2"
          >
            <RefreshCw className="mr-2 h-4 w-4" />
            Retry
          </button>
          <Link
            href="/"
            className="inline-flex items-center justify-center rounded-md text-sm font-medium border border-gray-300 bg-white text-gray-700 hover:bg-gray-50 h-10 px-4 py-2 transition-colors focus:outline-none focus:ring-2 focus:ring-blue-500 focus:ring-offset-2"
          >
            Back to Home
          </Link>
        </div>
      </div>
    );
  }

  // No training state
  if (!trainingState) {
    return (
      <div className="container mx-auto py-8 px-4 max-w-4xl">
        <div className="rounded-lg border border-blue-200 bg-blue-50 p-4">
          <div className="flex items-start gap-3">
            <AlertCircle className="h-5 w-5 text-blue-600 mt-0.5 flex-shrink-0" />
            <p className="text-sm text-blue-800">
              No training session found. Please configure and start training first.
            </p>
          </div>
        </div>
        <div className="mt-4 flex gap-4">
          <Link
            href="/training/configure"
            className="inline-flex items-center justify-center rounded-md text-sm font-medium bg-blue-600 text-white hover:bg-blue-700 h-10 px-4 py-2 transition-colors focus:outline-none focus:ring-2 focus:ring-blue-500 focus:ring-offset-2"
          >
            Configure Training
          </Link>
          <Link
            href="/"
            className="inline-flex items-center justify-center rounded-md text-sm font-medium border border-gray-300 bg-white text-gray-700 hover:bg-gray-50 h-10 px-4 py-2 transition-colors focus:outline-none focus:ring-2 focus:ring-blue-500 focus:ring-offset-2"
          >
            <ArrowLeft className="mr-2 h-4 w-4" />
            Back to Home
          </Link>
        </div>
      </div>
    );
  }

  return (
    <div className="container mx-auto py-8 px-4 max-w-7xl">
      {/* Header */}
      <div className="flex items-center justify-between mb-6">
        <div>
          <h1 className="text-3xl font-bold flex items-center gap-3">
            <Activity className="h-8 w-8" />
            Training Monitor
          </h1>
          <p className="text-gray-600 mt-1">
            Real-time progress tracking with 10-epoch smoothed metrics
          </p>
        </div>
        <div className="flex items-center gap-2">
          <button
            onClick={refresh}
            title='Refresh'
            className="inline-flex items-center justify-center rounded-md text-sm font-medium border border-gray-300 bg-white text-gray-700 hover:bg-gray-50 h-9 px-3 transition-colors focus:outline-none focus:ring-2 focus:ring-blue-500 focus:ring-offset-2"
          >
            <RefreshCw className="h-4 w-4" />
          </button>
          <Link
            href="/"
            className="inline-flex items-center justify-center rounded-md text-sm font-medium border border-gray-300 bg-white text-gray-700 hover:bg-gray-50 h-9 px-4 transition-colors focus:outline-none focus:ring-2 focus:ring-blue-500 focus:ring-offset-2"
          >
            <ArrowLeft className="mr-2 h-4 w-4" />
            Home
          </Link>
        </div>
      </div>

      {/* Status Overview */}
      <div className="grid md:grid-cols-2 gap-6 mb-8">
        <TrainingStatusCard 
          trainingState={trainingState}
          onStop={stopTraining}
        />
        <CurrentMetricsCard latestMetrics={trainingState.latest_metrics} />
      </div>

      {/* Charts */}
      {chartData.length > 0 && (
        <MetricsChartsSection chartData={chartData} />
      )}
    </div>
  );
};
