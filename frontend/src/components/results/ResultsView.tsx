// ========================================
// src/components/results/ResultsView.tsx
// ========================================
'use client';

import React from 'react';
import Link from 'next/link';
import { ArrowLeft, RefreshCw, AlertTriangle } from 'lucide-react';
import { Alert, AlertDescription, Button, Badge, Tabs, TabsContent, TabsList, TabsTrigger } from './ui';
import { ResultsHeader } from './ResultsHeader';
import { TrainingSummaryCard } from './TrainingSummaryCard';
import { MetricsTab } from './MetricsTab';
import { AnalysisTab } from './AnalysisTab';
import { DownloadSection } from './DownloadSection';
import { useResults } from '@/lib/hooks/useResults';

export const ResultsView: React.FC = () => {
  const {
    isLoading,
    error,
    trainingCompleted,
    trainingState,
    dataSummary,
    availableEmbeddings,
    isGeneratingEmbeddings,
    downloadProgress,
    downloadStatus,
    downloadEmbedding,
    downloadAllResults,
    generateEmbeddings,
    refresh,
    clearError,
  } = useResults();

  if (isLoading) {
    return (
      <div className="container mx-auto py-8 px-4 max-w-6xl">
        <div className="text-center">
          <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-blue-600 mx-auto mb-4"></div>
          <p>Loading results...</p>
        </div>
      </div>
    );
  }

  if (error && !trainingState) {
    return (
      <div className="container mx-auto py-8 px-4 max-w-4xl">
        <Alert variant="destructive" className="mb-4">
          <div className="flex items-start gap-3">
            <AlertTriangle className="h-5 w-5 mt-0.5" />
            <AlertDescription>{error}</AlertDescription>
          </div>
        </Alert>
        <div className="flex gap-4">
          <Button onClick={refresh}>
            <RefreshCw className="mr-2 h-4 w-4" />
            Retry
          </Button>
          <Link href="/">
            <Button variant="outline">Back to Home</Button>
          </Link>
        </div>
      </div>
    );
  }

  if (!trainingState || !dataSummary) {
    return (
      <div className="container mx-auto py-8 px-4 max-w-4xl">
        <Alert>
          <AlertDescription>
            No training data found. Please train a model first.
          </AlertDescription>
        </Alert>
        <div className="mt-4 flex gap-4">
          <Link href="/training/configure">
            <Button>Start Training</Button>
          </Link>
          <Link href="/">
            <Button variant="outline">
              <ArrowLeft className="mr-2 h-4 w-4" />
              Back to Home
            </Button>
          </Link>
        </div>
      </div>
    );
  }

  const chartData = trainingState.history.map((metrics) => ({
    epoch: metrics.epoch,
    loss: metrics.loss,
    ari: metrics.ari,
    nmi: metrics.nmi,
    asw: metrics.asw,
  }));

  return (
    <div className="container mx-auto py-8 px-4 max-w-7xl">
      <ResultsHeader
        trainingCompleted={trainingCompleted}
        cellCount={dataSummary.shape.n_obs}
        isLoading={isLoading}
        onRefresh={refresh}
      />

      {error && (
        <Alert variant="destructive" className="mb-6">
          <div className="flex items-center justify-between">
            <div className="flex items-start gap-3">
              <AlertTriangle className="h-5 w-5 mt-0.5" />
              <AlertDescription>{error}</AlertDescription>
            </div>
            <Button variant="outline" size="sm" onClick={clearError}>
              Dismiss
            </Button>
          </div>
        </Alert>
      )}

      <div className="space-y-6">
        <TrainingSummaryCard 
          trainingState={trainingState} 
          dataSummary={dataSummary} 
        />

        <Tabs defaultValue="downloads" className="space-y-4">
          <TabsList className="grid w-full grid-cols-3">
            <TabsTrigger value="downloads">
              <span className="flex items-center gap-2">
                Downloads
                {availableEmbeddings.length > 0 && (
                  <Badge variant="secondary" className="ml-2">
                    {availableEmbeddings.length}
                  </Badge>
                )}
              </span>
            </TabsTrigger>
            <TabsTrigger value="metrics">Training Metrics</TabsTrigger>
            <TabsTrigger value="analysis">Analysis</TabsTrigger>
          </TabsList>

          <TabsContent value="downloads">
            <DownloadSection
              availableEmbeddings={availableEmbeddings}
              isGeneratingEmbeddings={isGeneratingEmbeddings}
              downloadProgress={downloadProgress}
              downloadStatus={downloadStatus}
              onDownload={downloadEmbedding}
              onDownloadAll={downloadAllResults}
              onGenerate={generateEmbeddings}
              trainingCompleted={trainingCompleted}
            />
          </TabsContent>

          <TabsContent value="metrics">
            <MetricsTab history={trainingState.history} />
          </TabsContent>

          <TabsContent value="analysis">
            <AnalysisTab 
              latestMetrics={trainingState.latest_metrics}
              chartData={chartData}
            />
          </TabsContent>
        </Tabs>
      </div>
    </div>
  );
};