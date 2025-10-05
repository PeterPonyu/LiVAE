
// src/app/training/monitor/page.tsx
'use client';

import React, { useMemo } from 'react';
import { ArrowLeft, Square, RefreshCw, TrendingUp, Activity } from 'lucide-react';
import Link from 'next/link';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Progress } from '@/components/ui/progress';
import { Badge } from '@/components/ui/badge';
import { Alert, AlertDescription } from '@/components/ui/alert';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import { useTrainingMonitor } from '@/lib/hooks/useTrainingMonitor';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, Legend } from 'recharts';

// Metric configuration
const METRICS_CONFIG = {
  loss: { color: '#dc2626', label: 'Loss', betterDirection: 'down' },
  ari: { color: '#16a34a', label: 'ARI', betterDirection: 'up' },
  nmi: { color: '#2563eb', label: 'NMI', betterDirection: 'up' },
  asw: { color: '#7c3aed', label: 'ASW', betterDirection: 'up' },
  ch: { color: '#ea580c', label: 'CH', betterDirection: 'up' },
  db: { color: '#dc2626', label: 'DB', betterDirection: 'down' },
  pc: { color: '#059669', label: 'PC', betterDirection: 'up' },
} as const;

// Multi-line chart component for combined views
const CombinedChart = ({ 
  data, 
  metrics,
  height = 300 
}: {
  data: Record<string, unknown>[];
  metrics: Array<keyof typeof METRICS_CONFIG>;
  height?: number;
}) => (
  <div style={{ height }}>
    <ResponsiveContainer width="100%" height="100%">
      <LineChart data={data} margin={{ top: 5, right: 30, left: 20, bottom: 5 }}>
        <CartesianGrid strokeDasharray="3 3" opacity={0.2} />
        <XAxis 
          dataKey="epoch" 
          tick={{ fontSize: 12 }}
          label={{ value: 'Epoch', position: 'insideBottom', offset: -5 }}
        />
        <YAxis tick={{ fontSize: 12 }} />
        <Tooltip 
          contentStyle={{ backgroundColor: 'white', border: '1px solid #e5e7eb', borderRadius: '6px' }}
          formatter={(value: number) => value.toFixed(4)}
        />
        <Legend />
        {metrics.map(metric => (
          <Line
            key={metric}
            type="monotone"
            dataKey={metric}
            name={METRICS_CONFIG[metric].label}
            stroke={METRICS_CONFIG[metric].color}
            strokeWidth={2}
            dot={false}
            activeDot={{ r: 4 }}
          />
        ))}
      </LineChart>
    </ResponsiveContainer>
  </div>
);

export default function TrainingMonitorPage() {
  const { 
    trainingState, 
    isLoading, 
    error, 
    stopTraining, 
    refresh 
  } = useTrainingMonitor();

  // Calculate progress percentage
  const progressPercentage = useMemo(() => {
    if (!trainingState) return 0;
    return Math.round((trainingState.current_epoch / trainingState.total_epochs) * 100);
  }, [trainingState]);

  // Format training duration
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

  const latestMetrics = trainingState.latest_metrics;

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

      {/* Compact Status Overview */}
      <div className="grid md:grid-cols-2 gap-6 mb-8">
        {/* Progress & Status */}
        <Card>
          <CardHeader className="pb-3">
            <div className="flex items-center justify-between">
              <CardTitle>Training Status</CardTitle>
              <Badge variant={trainingState.is_running ? "destructive" : "secondary"}>
                {trainingState.is_running ? "Running" : "Stopped"}
              </Badge>
            </div>
          </CardHeader>
          <CardContent className="space-y-4">
            <div>
              <div className="flex justify-between text-sm mb-2 font-medium">
                <span>Epoch {trainingState.current_epoch} / {trainingState.total_epochs}</span>
                <span>{progressPercentage}%</span>
              </div>
              <Progress value={progressPercentage} className="h-2" />
            </div>
            <div className="flex justify-between text-sm text-gray-600">
              <span>Duration</span>
              <span className="font-medium">{formatDuration(trainingState.started_at)}</span>
            </div>
            
            {/* Action Buttons */}
            <div className="pt-2 flex gap-2">
              {trainingState.is_running ? (
                <Button 
                  variant="destructive" 
                  onClick={stopTraining}
                  className="w-full"
                  size="sm"
                >
                  <Square className="mr-2 h-4 w-4" />
                  Stop Training
                </Button>
              ) : (
                <>
                  <Button asChild className="flex-1" size="sm">
                    <Link href="/results">
                      <TrendingUp className="mr-2 h-4 w-4" />
                      View Results
                    </Link>
                  </Button>
                  <Button variant="outline" asChild className="flex-1" size="sm">
                    <Link href="/training/configure">
                      New Training
                    </Link>
                  </Button>
                </>
              )}
            </div>

            {trainingState.error_message && (
              <Alert variant="destructive" className="mt-2">
                <AlertDescription className="text-xs">
                  {trainingState.error_message}
                </AlertDescription>
              </Alert>
            )}
          </CardContent>
        </Card>

        {/* Current Metrics Grid */}
        <Card>
          <CardHeader className="pb-3">
            <CardTitle>Current Metrics</CardTitle>
            <CardDescription>
              {latestMetrics && `Epoch ${latestMetrics.epoch}`}
            </CardDescription>
          </CardHeader>
          <CardContent>
            {latestMetrics ? (
              <div className="grid grid-cols-4 gap-3">
                <div className="text-center">
                  <div className="text-xs text-gray-600 mb-1">Loss</div>
                  <div className="text-lg font-bold text-red-600">
                    {latestMetrics.loss.toFixed(3)}
                  </div>
                </div>
                <div className="text-center">
                  <div className="text-xs text-gray-600 mb-1">ARI</div>
                  <div className="text-lg font-bold text-green-600">
                    {latestMetrics.ari.toFixed(3)}
                  </div>
                </div>
                <div className="text-center">
                  <div className="text-xs text-gray-600 mb-1">NMI</div>
                  <div className="text-lg font-bold text-blue-600">
                    {latestMetrics.nmi.toFixed(3)}
                  </div>
                </div>
                <div className="text-center">
                  <div className="text-xs text-gray-600 mb-1">ASW</div>
                  <div className="text-lg font-bold text-purple-600">
                    {latestMetrics.asw.toFixed(3)}
                  </div>
                </div>
                <div className="text-center">
                  <div className="text-xs text-gray-600 mb-1">CH</div>
                  <div className="text-sm font-semibold text-orange-600">
                    {latestMetrics.ch.toFixed(1)}
                  </div>
                </div>
                <div className="text-center">
                  <div className="text-xs text-gray-600 mb-1">DB</div>
                  <div className="text-sm font-semibold text-red-500">
                    {latestMetrics.db.toFixed(3)}
                  </div>
                </div>
                <div className="text-center">
                  <div className="text-xs text-gray-600 mb-1">PC</div>
                  <div className="text-sm font-semibold text-green-500">
                    {latestMetrics.pc.toFixed(3)}
                  </div>
                </div>
                <div className="text-center flex items-center justify-center">
                  <Badge variant="outline" className="text-xs">
                    {chartData.length} pts
                  </Badge>
                </div>
              </div>
            ) : (
              <p className="text-gray-500 text-sm">No metrics available yet</p>
            )}
          </CardContent>
        </Card>
      </div>

      {/* Charts */}
      {chartData.length > 0 && (
        <Card>
          <CardHeader>
            <CardTitle>Metrics Visualization</CardTitle>
            <CardDescription>
              Training progress across all evaluation metrics
            </CardDescription>
          </CardHeader>
          <CardContent>
            <Tabs defaultValue="primary" className="space-y-4">
              <TabsList className="grid w-full grid-cols-3">
                <TabsTrigger value="primary">Primary Metrics</TabsTrigger>
                <TabsTrigger value="clustering">Clustering Quality</TabsTrigger>
                <TabsTrigger value="all">All Metrics</TabsTrigger>
              </TabsList>

              <TabsContent value="primary" className="space-y-4">
                <div className="grid md:grid-cols-2 gap-6">
                  <div>
                    <h3 className="text-sm font-medium mb-2 text-gray-700">Training Loss</h3>
                    <CombinedChart data={chartData} metrics={['loss']} height={250} />
                  </div>
                  <div>
                    <h3 className="text-sm font-medium mb-2 text-gray-700">Clustering Accuracy (ARI)</h3>
                    <CombinedChart data={chartData} metrics={['ari']} height={250} />
                  </div>
                </div>
              </TabsContent>

              <TabsContent value="clustering" className="space-y-4">
                <div>
                  <h3 className="text-sm font-medium mb-2 text-gray-700">Quality Metrics Comparison</h3>
                  <CombinedChart 
                    data={chartData} 
                    metrics={['ari', 'nmi', 'asw', 'pc']} 
                    height={350} 
                  />
                </div>
                <div className="grid md:grid-cols-3 gap-4 text-xs">
                  <div className="bg-green-50 p-3 rounded">
                    <strong>ARI:</strong> Clustering accuracy vs ground truth
                  </div>
                  <div className="bg-blue-50 p-3 rounded">
                    <strong>NMI:</strong> Mutual information measure
                  </div>
                  <div className="bg-purple-50 p-3 rounded">
                    <strong>ASW:</strong> Cluster separation quality
                  </div>
                </div>
              </TabsContent>

              <TabsContent value="all" className="space-y-4">
                <div className="grid md:grid-cols-2 gap-6">
                  <div>
                    <h3 className="text-sm font-medium mb-2 text-gray-700">Loss & Accuracy</h3>
                    <CombinedChart data={chartData} metrics={['loss', 'ari']} height={250} />
                  </div>
                  <div>
                    <h3 className="text-sm font-medium mb-2 text-gray-700">Separation Metrics</h3>
                    <CombinedChart data={chartData} metrics={['asw', 'ch', 'db']} height={250} />
                  </div>
                </div>
                <div>
                  <h3 className="text-sm font-medium mb-2 text-gray-700">Information & Connectivity</h3>
                  <CombinedChart data={chartData} metrics={['nmi', 'pc']} height={250} />
                </div>
              </TabsContent>
            </Tabs>
          </CardContent>
        </Card>
      )}
    </div>
  );
}
