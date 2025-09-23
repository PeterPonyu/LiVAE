
// src/app/results/page.tsx
'use client';

import React from 'react';
import { ArrowLeft, RefreshCw, BarChart3, AlertTriangle, CheckCircle} from 'lucide-react';
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
          <Button variant="outline" onClick={clearError}>
            Clear Error
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

  // Prepare chart data for the summary
  const chartData = trainingState.history.map((metrics: TrainingMetrics) => ({
    epoch: metrics.epoch,
    loss: metrics.loss,
    ari: metrics.ari,
    nmi: metrics.nmi,
    asw: metrics.asw,
  }));

  return (
    <div className="container mx-auto py-8 px-4 max-w-7xl">
      {/* Header */}
      <div className="flex items-center justify-between mb-6">
        <div>
          <h1 className="text-3xl font-bold">Training Results</h1>
          <p className="text-gray-600 mt-1 flex items-center space-x-2">
            <span>{trainingCompleted ? 'Training completed successfully' : 'Training in progress or incomplete'}</span>
            <Badge variant={trainingCompleted ? "default" : "secondary"}>
              {trainingCompleted ? "✅ Complete" : "⏳ In Progress"}
            </Badge>
          </p>
        </div>
        <div className="flex items-center space-x-2">
          <Button variant="outline" onClick={refresh} size="sm" disabled={isLoading}>
            <RefreshCw className={`mr-2 h-4 w-4 ${isLoading ? 'animate-spin' : ''}`} />
            Refresh
          </Button>
          <Button variant="outline" asChild>
            <Link href="/">
              <ArrowLeft className="mr-2 h-4 w-4" />
              Home
            </Link>
          </Button>
        </div>
      </div>

      {/* Status Alerts */}
      {!trainingCompleted && (
        <Alert className="mb-6">
          <AlertDescription>
            Training is not yet complete. Some results may not be available until training finishes.
            <Link href="/training/monitor" className="ml-2 text-blue-600 hover:underline">
              Monitor training progress →
            </Link>
          </AlertDescription>
        </Alert>
      )}

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

      <div className="space-y-8">
        {/* Training Summary */}
        <TrainingSummaryCard 
          trainingState={trainingState} 
          dataSummary={dataSummary} 
        />

        {/* Results Tabs */}
        <Tabs defaultValue="downloads" className="space-y-6">
          <TabsList className="grid w-full grid-cols-3">
            <TabsTrigger value="downloads" className="flex items-center space-x-2">
              <span>Downloads</span>
              {availableEmbeddings.length > 0 && (
                <Badge variant="secondary" className="ml-1">
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
                <CardTitle className="flex items-center space-x-2">
                  <BarChart3 className="h-5 w-5" />
                  <span>Training Metrics Overview</span>
                </CardTitle>
                <CardDescription>
                  Complete training history and performance analysis ({chartData.length} data points)
                </CardDescription>
              </CardHeader>
              <CardContent>
                {chartData.length > 0 ? (
                  <div className="space-y-8">
                    {/* Combined Overview Chart */}
                    <div>
                      <h4 className="font-medium mb-3">Training Overview</h4>
                      <div className="h-80">
                        <ResponsiveContainer width="100%" height="100%">
                          <LineChart data={chartData}>
                            <CartesianGrid strokeDasharray="3 3" />
                            <XAxis dataKey="epoch" />
                            <YAxis yAxisId="left" orientation="left" domain={['auto', 'auto']} />
                            <YAxis yAxisId="right" orientation="right" domain={[0, 1]} />
                            <Tooltip 
                              formatter={(value: any, name: string) => [
                                typeof value === 'number' ? value.toFixed(4) : value,
                                name.toUpperCase()
                              ]}
                            />
                            <Legend />
                            <Line 
                              yAxisId="left" 
                              type="monotone" 
                              dataKey="loss" 
                              stroke="#dc2626" 
                              strokeWidth={2} 
                              name="Loss" 
                            />
                            <Line 
                              yAxisId="right" 
                              type="monotone" 
                              dataKey="ari" 
                              stroke="#16a34a" 
                              strokeWidth={2} 
                              name="ARI" 
                            />
                            <Line 
                              yAxisId="right" 
                              type="monotone" 
                              dataKey="nmi" 
                              stroke="#2563eb" 
                              strokeWidth={2} 
                              name="NMI" 
                            />
                            <Line 
                              yAxisId="right" 
                              type="monotone" 
                              dataKey="asw" 
                              stroke="#7c3aed" 
                              strokeWidth={2} 
                              name="ASW" 
                            />
                          </LineChart>
                        </ResponsiveContainer>
                      </div>
                    </div>

                    {/* Metrics Summary Cards */}
                    <div className="grid md:grid-cols-2 lg:grid-cols-4 gap-4">
                      {trainingState.latest_metrics && (
                        <>
                          <Card className="text-center">
                            <CardContent className="p-4">
                              <div className="text-2xl font-bold text-red-600">
                                {trainingState.latest_metrics.loss.toFixed(4)}
                              </div>
                              <div className="text-sm text-gray-600">Final Loss</div>
                            </CardContent>
                          </Card>
                          <Card className="text-center">
                            <CardContent className="p-4">
                              <div className="text-2xl font-bold text-green-600">
                                {trainingState.latest_metrics.ari.toFixed(3)}
                              </div>
                              <div className="text-sm text-gray-600">Final ARI</div>
                            </CardContent>
                          </Card>
                          <Card className="text-center">
                            <CardContent className="p-4">
                              <div className="text-2xl font-bold text-blue-600">
                                {trainingState.latest_metrics.nmi.toFixed(3)}
                              </div>
                              <div className="text-sm text-gray-600">Final NMI</div>
                            </CardContent>
                          </Card>
                          <Card className="text-center">
                            <CardContent className="p-4">
                              <div className="text-2xl font-bold text-purple-600">
                                {trainingState.latest_metrics.asw.toFixed(3)}
                              </div>
                              <div className="text-sm text-gray-600">Final ASW</div>
                            </CardContent>
                          </Card>
                        </>
                      )}
                    </div>

                    {/* Metrics Table */}
                    <div>
                      <h4 className="font-medium mb-3">Recent Training History</h4>
                      <div className="overflow-x-auto">
                        <table className="w-full text-sm border-collapse border border-gray-200">
                          <thead>
                            <tr className="bg-gray-50">
                              <th className="border border-gray-200 p-2 text-left">Epoch</th>
                              <th className="border border-gray-200 p-2 text-left">Loss</th>
                              <th className="border border-gray-200 p-2 text-left">ARI</th>
                              <th className="border border-gray-200 p-2 text-left">NMI</th>
                              <th className="border border-gray-200 p-2 text-left">ASW</th>
                            </tr>
                          </thead>
                          <tbody>
                            {trainingState.history.slice(-10).reverse().map((metrics: TrainingMetrics) => (
                              <tr key={metrics.epoch} className="hover:bg-gray-50">
                                <td className="border border-gray-200 p-2 font-medium">{metrics.epoch}</td>
                                <td className="border border-gray-200 p-2 text-red-600">{metrics.loss.toFixed(4)}</td>
                                <td className="border border-gray-200 p-2 text-green-600">{metrics.ari.toFixed(3)}</td>
                                <td className="border border-gray-200 p-2 text-blue-600">{metrics.nmi.toFixed(3)}</td>
                                <td className="border border-gray-200 p-2 text-purple-600">{metrics.asw.toFixed(3)}</td>
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
            <Card>
              <CardHeader>
                <CardTitle>Model Analysis & Recommendations</CardTitle>
                <CardDescription>
                  Interpretation and recommendations based on training results
                </CardDescription>
              </CardHeader>
              <CardContent>
                <div className="space-y-6">
                  {/* Performance Analysis */}
                  {trainingState.latest_metrics && (
                    <div className="space-y-4">
                      <h4 className="font-medium">Performance Assessment</h4>
                      <div className="grid md:grid-cols-2 gap-4">
                        <Card className="bg-green-50 border-green-200">
                          <CardContent className="p-4">
                            <h5 className="font-medium text-green-800 mb-3">Clustering Quality</h5>
                            <div className="space-y-3 text-sm">
                              <div className="flex justify-between items-center">
                                <span>ARI Score ({trainingState.latest_metrics.ari.toFixed(3)}):</span>
                                <Badge variant={
                                  trainingState.latest_metrics.ari > 0.7 ? "default" : 
                                  trainingState.latest_metrics.ari > 0.4 ? "secondary" : "destructive"
                                }>
                                  {trainingState.latest_metrics.ari > 0.7 ? "Excellent" : 
                                   trainingState.latest_metrics.ari > 0.4 ? "Good" : "Needs Improvement"}
                                </Badge>
                              </div>
                              <div className="flex justify-between items-center">
                                <span>NMI Score ({trainingState.latest_metrics.nmi.toFixed(3)}):</span>
                                <Badge variant={
                                  trainingState.latest_metrics.nmi > 0.7 ? "default" : 
                                  trainingState.latest_metrics.nmi > 0.4 ? "secondary" : "destructive"
                                }>
                                  {trainingState.latest_metrics.nmi > 0.7 ? "Excellent" : 
                                   trainingState.latest_metrics.nmi > 0.4 ? "Good" : "Needs Improvement"}
                                </Badge>
                              </div>
                            </div>
                          </CardContent>
                        </Card>

                        <Card className="bg-blue-50 border-blue-200">
                          <CardContent className="p-4">
                            <h5 className="font-medium text-blue-800 mb-3">Embedding Quality</h5>
                            <div className="space-y-3 text-sm">
                              <div className="flex justify-between items-center">
                                <span>Silhouette Width ({trainingState.latest_metrics.asw.toFixed(3)}):</span>
                                <Badge variant={
                                  trainingState.latest_metrics.asw > 0.5 ? "default" : 
                                  trainingState.latest_metrics.asw > 0.2 ? "secondary" : "destructive"
                                }>
                                  {trainingState.latest_metrics.asw > 0.5 ? "Well Separated" : 
                                   trainingState.latest_metrics.asw > 0.2 ? "Moderately Separated" : "Poorly Separated"}
                                </Badge>
                              </div>
                              <div className="flex justify-between items-center">
                                <span>Loss Convergence:</span>
                                <Badge variant={
                                  chartData.length > 5 && 
                                  chartData[chartData.length - 1].loss < chartData[Math.max(0, chartData.length - 6)].loss * 0.9 
                                    ? "default" : "secondary"
                                }>
                                  {chartData.length > 5 && 
                                   chartData[chartData.length - 1].loss < chartData[Math.max(0, chartData.length - 6)].loss * 0.9 
                                     ? "Converged" : "Still Learning"}
                                </Badge>
                              </div>
                            </div>
                          </CardContent>
                        </Card>
                      </div>
                    </div>
                  )}

                  {/* Recommendations */}
                  <div>
                    <h4 className="font-medium mb-3">Recommendations</h4>
                    <div className="space-y-3">
                      {trainingState.latest_metrics && (
                        <>
                          {trainingState.latest_metrics.ari < 0.4 && (
                            <Alert>
                              <AlertTriangle className="h-4 w-4" />
                              <AlertDescription>
                                <strong>Low clustering quality detected.</strong> Consider increasing training epochs, adjusting regularization parameters (beta, tc), or checking data quality.
                              </AlertDescription>
                            </Alert>
                          )}
                          {trainingState.latest_metrics.asw < 0.2 && (
                            <Alert>
                              <AlertTriangle className="h-4 w-4" />
                              <AlertDescription>
                                <strong>Poor cluster separation.</strong> Consider adjusting the latent dimension, beta parameter, or using different loss weights.
                              </AlertDescription>
                            </Alert>
                          )}
                          {trainingState.latest_metrics.ari > 0.7 && trainingState.latest_metrics.nmi > 0.7 && (
                            <Alert className="bg-green-50 border-green-200">
                              <CheckCircle className="h-4 w-4 text-green-600" />
                              <AlertDescription className="text-green-800">
                                <strong>Excellent results!</strong> Your model has learned meaningful cell clusters. The embeddings are ready for downstream analysis and publication.
                              </AlertDescription>
                            </Alert>
                          )}
                        </>
                      )}
                      
                      <div className="bg-gray-50 p-4 rounded">
                        <h5 className="font-medium mb-2">🔬 Next Steps for Analysis</h5>
                        <ul className="space-y-1 text-sm">
                          <li>• <strong>Download interpretable embedding</strong> for 2D visualization (UMAP plots)</li>
                          <li>• <strong>Download latent embedding</strong> for downstream analysis (clustering, trajectory)</li>
                          <li>• <strong>Visualize clusters</strong> in your preferred analysis tool (Scanpy, Seurat)</li>
                          <li>• <strong>Validate biological relevance</strong> of discovered clusters with marker genes</li>
                          <li>• <strong>Export results</strong> for collaborative analysis and publication</li>
                        </ul>
                      </div>
                    </div>
                  </div>
                </div>
              </CardContent>
            </Card>
          </TabsContent>
        </Tabs>
      </div>
    </div>
  );
}
