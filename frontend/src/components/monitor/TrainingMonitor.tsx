// ========================================
// src/components/monitor/TrainingMonitor.tsx
// ========================================
'use client';

import React, { useMemo } from 'react';
import { ArrowLeft, RefreshCw, Activity } from 'lucide-react';
import Link from 'next/link';
import { Button } from '@/components/ui/button';
import { Alert, AlertDescription } from '@/components/ui/alert';
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
        <Alert variant="destructive">
          <AlertDescription>{error}</AlertDescription>
        </Alert>
        <div className="mt-4 flex gap-4">
          <Button onClick={refresh}>
            <RefreshCw className="mr-2 h-4 w-4" />
            Retry
          </Button>
          <Button variant="outline" asChild>
            <Link href="/">Back to Home</Link>
          </Button>
        </div>
      </div>
    );
  }

  // No training state
  if (!trainingState) {
    return (
      <div className="container mx-auto py-8 px-4 max-w-4xl">
        <Alert>
          <AlertDescription>
            No training session found. Please configure and start training first.
          </AlertDescription>
        </Alert>
        <div className="mt-4 flex gap-4">
          <Button asChild>
            <Link href="/training/configure">Configure Training</Link>
          </Button>
          <Button variant="outline" asChild>
            <Link href="/">
              <ArrowLeft className="mr-2 h-4 w-4" />
              Back to Home
            </Link>
          </Button>
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
          <Button variant="outline" onClick={refresh} size="sm">
            <RefreshCw className="h-4 w-4" />
          </Button>
          <Button variant="outline" asChild>
            <Link href="/">
              <ArrowLeft className="mr-2 h-4 w-4" />
              Home
            </Link>
          </Button>
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