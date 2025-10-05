
// src/app/results/page.tsx
'use client';

import React from 'react';
import { ArrowLeft, RefreshCw, BarChart3, AlertTriangle, CheckCircle, TrendingUp } from 'lucide-react';
import Link from 'next/link';
import { Button } from '@/components/ui/button';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Alert, AlertDescription } from '@/components/ui/alert';
import { Badge } from '@/components/ui/badge';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import { TrainingSummaryCard } from '@/components/results/TrainingSummaryCard';
import { DownloadSection } from '@/components/results/DownloadSection';
import { useResults } from '@/lib/hooks/useResults';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, Legend } from 'recharts';
import { type TrainingMetrics } from '@/types';

// Quality assessment helper
const getQualityBadge = (value: number, thresholds: [number, number]) => {
  const [good, excellent] = thresholds;
  if (value >= excellent) return { variant: 'default' as const, label: 'Excellent' };
  if (value >= good) return { variant: 'secondary' as const, label: 'Good' };
  return { variant: 'destructive' as const, label: 'Poor' };
};

export default function ResultsPage() {
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

  if (error) {
    return (
      <div className="container mx-auto py-8 px-4 max-w-4xl">
        <Alert variant="destructive" className="mb-4">
          <AlertTriangle className="h-4 w-4" />
          <AlertDescription>{error}</AlertDescription>
        </Alert>
        <div className="flex gap-4">
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

  if (!trainingState || !dataSummary) {
    return (
      <div className="container mx-auto py-8 px-4 max-w-4xl">
        <Alert>
          <AlertDescription>
            No training data found. Please train a model first.
          </AlertDescription>
        </Alert>
        <div className="mt-4 flex gap-4">
          <Button asChild>
            <Link href="/training/configure">Start Training</Link>
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

  const chartData = trainingState.history.map((metrics: TrainingMetrics) => ({
    epoch: metrics.epoch,
    loss: metrics.loss,
    ari: metrics.ari,
    nmi: metrics.nmi,
    asw: metrics.asw,
  }));

  const latestMetrics = trainingState.latest_metrics;

  return (
    <div className="container mx-auto py-8 px-4 max-w-7xl">
      {/* Header */}
      <div className="flex items-center justify-between mb-6">
        <div>
          <div className="flex items-center gap-3">
            <h1 className="text-3xl font-bold">Training Results</h1>
            <Badge variant={trainingCompleted ? "default" : "secondary"}>
              {trainingCompleted ? "Complete" : "In Progress"}
            </Badge>
          </div>
          <p className="text-gray-600 mt-1">
            {trainingCompleted 
              ? `Successfully trained on ${dataSummary.shape.n_obs.toLocaleString()} cells`
              : 'Training incomplete - some features may be unavailable'}
          </p>
        </div>
        <div className="flex items-center gap-2">
          <Button variant="outline" onClick={refresh} size="sm" disabled={isLoading}>
            <RefreshCw className={`h-4 w-4 ${isLoading ? 'animate-spin' : ''}`} />
          </Button>
          {!trainingCompleted && (
            <Button variant="outline" size="sm" asChild>
              <Link href="/training/monitor">
                <TrendingUp className="mr-2 h-4 w-4" />
                Monitor
              </Link>
            </Button>
          )}
          <Button variant="outline" asChild>
            <Link href="/">
              <ArrowLeft className="mr-2 h-4 w-4" />
              Home
            </Link>
          </Button>
        </div>
      </div>

      {/* Error Alert */}
      {error && (
        <Alert variant="destructive" className="mb-6">
          <AlertTriangle className="h-4 w-4" />
          <AlertDescription className="flex items-center justify-between">
            <span>{error}</span>
            <Button variant="outline" size="sm" onClick={clearError}>
              Dismiss
            </Button>
          </AlertDescription>
        </Alert>
      )}

      <div className="space-y-6">
        {/* Training Summary */}
        <TrainingSummaryCard 
          trainingState={trainingState} 
          dataSummary={dataSummary} 
        />

        {/* Results Tabs */}
        <Tabs defaultValue="downloads" className="space-y-4">
          <TabsList className="grid w-full grid-cols-3">
            <TabsTrigger value="downloads">
              Downloads
              {availableEmbeddings.length > 0 && (
                <Badge variant="secondary" className="ml-2">
                  {availableEmbeddings.length}
                </Badge>
              )}
            </TabsTrigger>
            <TabsTrigger value="metrics">Training Metrics</TabsTrigger>
            <TabsTrigger value="analysis">Analysis</TabsTrigger>
          </TabsList>

          {/* Downloads Tab */}
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

          {/* Metrics Tab */}
          <TabsContent value="metrics">
            <Card>
              <CardHeader>
                <CardTitle className="flex items-center gap-2">
                  <BarChart3 className="h-5 w-5" />
                  Training Progress
                </CardTitle>
                <CardDescription>
                  {chartData.length} epochs tracked
                </CardDescription>
              </CardHeader>
              <CardContent>
                {chartData.length > 0 ? (
                  <div className="space-y-6">
                    {/* Main Chart */}
                    <div className="h-80">
                      <ResponsiveContainer width="100%" height="100%">
                        <LineChart data={chartData}>
                          <CartesianGrid strokeDasharray="3 3" opacity={0.2} />
                          <XAxis 
                            dataKey="epoch" 
                            label={{ value: 'Epoch', position: 'insideBottom', offset: -5 }}
                          />
                          <YAxis yAxisId="left" orientation="left" domain={['auto', 'auto']} />
                          <YAxis yAxisId="right" orientation="right" domain={[0, 1]} />
                          <Tooltip 
                            formatter={(value: number) => value.toFixed(4)}
                            contentStyle={{ backgroundColor: 'white', border: '1px solid #e5e7eb', borderRadius: '6px' }}
                          />
                          <Legend />
                          <Line 
                            yAxisId="left" 
                            type="monotone" 
                            dataKey="loss" 
                            stroke="#dc2626" 
                            strokeWidth={2} 
                            name="Loss" 
                            dot={false}
                          />
                          <Line 
                            yAxisId="right" 
                            type="monotone" 
                            dataKey="ari" 
                            stroke="#16a34a" 
                            strokeWidth={2} 
                            name="ARI" 
                            dot={false}
                          />
                          <Line 
                            yAxisId="right" 
                            type="monotone" 
                            dataKey="nmi" 
                            stroke="#2563eb" 
                            strokeWidth={2} 
                            name="NMI" 
                            dot={false}
                          />
                          <Line 
                            yAxisId="right" 
                            type="monotone" 
                            dataKey="asw" 
                            stroke="#7c3aed" 
                            strokeWidth={2} 
                            name="ASW" 
                            dot={false}
                          />
                        </LineChart>
                      </ResponsiveContainer>
                    </div>

                    {/* Recent History Table */}
                    <div>
                      <h4 className="text-sm font-medium mb-3 text-gray-700">Last 10 Epochs</h4>
                      <div className="overflow-x-auto rounded-lg border">
                        <table className="w-full text-sm">
                          <thead className="bg-gray-50">
                            <tr>
                              <th className="p-3 text-left font-medium">Epoch</th>
                              <th className="p-3 text-left font-medium">Loss</th>
                              <th className="p-3 text-left font-medium">ARI</th>
                              <th className="p-3 text-left font-medium">NMI</th>
                              <th className="p-3 text-left font-medium">ASW</th>
                            </tr>
                          </thead>
                          <tbody className="divide-y">
                            {trainingState.history.slice(-10).reverse().map((metrics: TrainingMetrics) => (
                              <tr key={metrics.epoch} className="hover:bg-gray-50">
                                <td className="p-3 font-medium">{metrics.epoch}</td>
                                <td className="p-3 text-red-600">{metrics.loss.toFixed(4)}</td>
                                <td className="p-3 text-green-600">{metrics.ari.toFixed(3)}</td>
                                <td className="p-3 text-blue-600">{metrics.nmi.toFixed(3)}</td>
                                <td className="p-3 text-purple-600">{metrics.asw.toFixed(3)}</td>
                              </tr>
                            ))}
                          </tbody>
                        </table>
                      </div>
                    </div>
                  </div>
                ) : (
                  <div className="text-center py-8 text-gray-500">
                    <BarChart3 className="h-8 w-8 mx-auto mb-2" />
                    <p>No training metrics available</p>
                  </div>
                )}
              </CardContent>
            </Card>
          </TabsContent>

          {/* Analysis Tab */}
          <TabsContent value="analysis">
            <div className="space-y-4">
              {/* Performance Overview */}
              {latestMetrics && (
                <Card>
                  <CardHeader>
                    <CardTitle>Performance Assessment</CardTitle>
                    <CardDescription>Quality metrics from final epoch</CardDescription>
                  </CardHeader>
                  <CardContent>
                    <div className="grid md:grid-cols-2 gap-4">
                      {/* Clustering Quality */}
                      <div className="space-y-3">
                        <h4 className="font-medium text-sm">Clustering Quality</h4>
                        <div className="space-y-2">
                          <div className="flex justify-between items-center p-2 bg-gray-50 rounded">
                            <span className="text-sm">ARI: {latestMetrics.ari.toFixed(3)}</span>
                            <Badge {...getQualityBadge(latestMetrics.ari, [0.4, 0.7])} />
                          </div>
                          <div className="flex justify-between items-center p-2 bg-gray-50 rounded">
                            <span className="text-sm">NMI: {latestMetrics.nmi.toFixed(3)}</span>
                            <Badge {...getQualityBadge(latestMetrics.nmi, [0.4, 0.7])} />
                          </div>
                        </div>
                      </div>

                      {/* Separation Quality */}
                      <div className="space-y-3">
                        <h4 className="font-medium text-sm">Embedding Quality</h4>
                        <div className="space-y-2">
                          <div className="flex justify-between items-center p-2 bg-gray-50 rounded">
                            <span className="text-sm">ASW: {latestMetrics.asw.toFixed(3)}</span>
                            <Badge {...getQualityBadge(latestMetrics.asw, [0.2, 0.5])} />
                          </div>
                          <div className="flex justify-between items-center p-2 bg-gray-50 rounded">
                            <span className="text-sm">Final Loss: {latestMetrics.loss.toFixed(4)}</span>
                            <Badge variant="secondary">
                              {chartData.length > 5 && 
                               chartData[chartData.length - 1].loss < chartData[chartData.length - 6].loss * 0.9
                                ? "Converged" : "Stable"}
                            </Badge>
                          </div>
                        </div>
                      </div>
                    </div>
                  </CardContent>
                </Card>
              )}

              {/* Recommendations */}
              <Card>
                <CardHeader>
                  <CardTitle>Recommendations</CardTitle>
                </CardHeader>
                <CardContent className="space-y-3">
                  {latestMetrics && (
                    <>
                      {/* Success Message */}
                      {latestMetrics.ari > 0.7 && latestMetrics.nmi > 0.7 && (
                        <Alert className="bg-green-50 border-green-200">
                          <CheckCircle className="h-4 w-4 text-green-600" />
                          <AlertDescription className="text-green-800">
                            <strong>Excellent results!</strong> Your model has learned high-quality cell embeddings suitable for downstream analysis.
                          </AlertDescription>
                        </Alert>
                      )}
                      
                      {/* Warning Messages */}
                      {latestMetrics.ari < 0.4 && (
                        <Alert>
                          <AlertTriangle className="h-4 w-4" />
                          <AlertDescription>
                            <strong>Low clustering accuracy.</strong> Consider training longer or adjusting hyperparameters (beta, tc, latent_dim).
                          </AlertDescription>
                        </Alert>
                      )}
                      
                      {latestMetrics.asw < 0.2 && latestMetrics.ari >= 0.4 && (
                        <Alert>
                          <AlertTriangle className="h-4 w-4" />
                          <AlertDescription>
                            <strong>Clusters need better separation.</strong> Try increasing beta or adjusting the latent dimension.
                          </AlertDescription>
                        </Alert>
                      )}
                    </>
                  )}
                  
                  {/* Next Steps */}
                  <div className="bg-blue-50 p-4 rounded-lg border border-blue-200">
                    <h5 className="font-medium mb-2 text-blue-900">Next Steps</h5>
                    <ul className="space-y-1 text-sm text-blue-800">
                      <li>• Download embeddings for visualization and analysis</li>
                      <li>• Use latent embeddings for clustering and trajectory inference</li>
                      <li>• Validate clusters with marker gene expression</li>
                      <li>• Export results for collaborative research</li>
                    </ul>
                  </div>
                </CardContent>
              </Card>
            </div>
          </TabsContent>
        </Tabs>
      </div>
    </div>
  );
}
