
// src/app/training/monitor/page.tsx
'use client';

import React, { useMemo } from 'react';
import { ArrowLeft, Square, RefreshCw, TrendingUp } from 'lucide-react';
import Link from 'next/link';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Progress } from '@/components/ui/progress';
import { Badge } from '@/components/ui/badge';
import { Alert, AlertDescription } from '@/components/ui/alert';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import { useTrainingMonitor } from '@/lib/hooks/useTrainingMonitor';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer } from 'recharts';

// Custom tooltip component for better formatting
const CustomTooltip = ({ active, payload, label, formatValue }: any) => {
  if (active && payload && payload.length) {
    const value = payload[0].value;
    return (
      <div className="bg-white p-3 border rounded shadow-lg">
        <p className="font-medium">{`Epoch: ${label}`}</p>
        <p className="text-blue-600">
          {`${payload[0].name}: ${formatValue ? formatValue(value) : value.toFixed(4)}`}
        </p>
      </div>
    );
  }
  return null;
};

// Individual chart component
const MetricChart = ({ 
  data, 
  dataKey, 
  title, 
  description, 
  color, 
  formatValue, 
  domain,
  height = 250 
}: {
  data: any[];
  dataKey: string;
  title: string;
  description: string;
  color: string;
  formatValue?: (value: number) => string;
  domain?: [number | 'auto', number | 'auto'];
  height?: number;
}) => {
  // Calculate trend for the metric
  const calculateTrend = () => {
    if (data.length < 2) return null;
    const recent = data.slice(-3); // Last 3 points
    if (recent.length < 2) return null;
    
    const firstValue = recent[0][dataKey];
    const lastValue = recent[recent.length - 1][dataKey];
    const change = lastValue - firstValue;
    
    // Determine if positive change is good (for most metrics, higher is better, except loss and DB)
    const lowerIsBetter = dataKey === 'loss' || dataKey === 'db';
    const isImproving = lowerIsBetter ? change < 0 : change > 0;
    
    return {
      change: Math.abs(change),
      isImproving,
      direction: change > 0 ? 'up' : 'down'
    };
  };

  const trend = calculateTrend();

  return (
    <Card>
      <CardHeader className="pb-3">
        <CardTitle className="text-lg flex items-center justify-between">
          {title}
          {trend && (
            <Badge variant={trend.isImproving ? "default" : "secondary"} className="text-xs">
              {trend.isImproving ? '📈' : '📉'} 
              {trend.direction === 'up' ? '↗' : '↘'}
            </Badge>
          )}
        </CardTitle>
        <CardDescription className="text-sm">
          {description}
          <br />
          <span className="text-xs text-blue-600">📊 10-epoch rolling average</span>
        </CardDescription>
      </CardHeader>
      <CardContent>
        <div style={{ height }}>
          <ResponsiveContainer width="100%" height="100%">
            <LineChart data={data} margin={{ top: 5, right: 30, left: 20, bottom: 5 }}>
              <CartesianGrid strokeDasharray="3 3" opacity={0.3} />
              <XAxis 
                dataKey="epoch" 
                tick={{ fontSize: 12 }}
                tickLine={{ stroke: '#e0e0e0' }}
              />
              <YAxis 
                tick={{ fontSize: 12 }}
                tickLine={{ stroke: '#e0e0e0' }}
                domain={domain}
                tickFormatter={formatValue}
              />
              <Tooltip 
                content={({ active, payload, label }) => {
                  if (active && payload && payload.length) {
                    const value = payload[0].value;
                    return (
                      <div className="bg-white p-3 border rounded shadow-lg">
                        <p className="font-medium">{`Epoch: ${label}`}</p>
                        <p className="text-blue-600">
                          {`${payload[0].name}: ${formatValue ? formatValue(value) : value.toFixed(4)}`}
                        </p>
                        <p className="text-xs text-gray-500">10-epoch average</p>
                      </div>
                    );
                  }
                  return null;
                }}
              />
              <Line 
                type="monotone" 
                dataKey={dataKey} 
                stroke={color} 
                strokeWidth={3}  // Slightly thicker for averaged data
                dot={false}
                activeDot={{ r: 5, stroke: color, strokeWidth: 2 }}
              />
            </LineChart>
          </ResponsiveContainer>
        </div>
      </CardContent>
    </Card>
  );
};

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
    
    if (diffMins > 0) {
      return `${diffMins}m ${diffSecs}s`;
    }
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
        <div className="mt-4 space-x-4">
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
        <div className="mt-4 space-x-4">
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
          <h1 className="text-3xl font-bold">Training Monitor</h1>
          <p className="text-gray-600 mt-1">
            Real-time training progress and metrics
          </p>
        </div>
        <div className="flex items-center space-x-2">
          <Button variant="outline" onClick={refresh} size="sm">
            <RefreshCw className="mr-2 h-4 w-4" />
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

      {/* Training Status */}
      <div className="grid md:grid-cols-3 gap-6 mb-8">
        {/* Progress Card */}
        <Card>
          <CardHeader>
            <CardTitle className="flex items-center justify-between">
              Training Progress
              <Badge variant={trainingState.is_running ? "destructive" : "secondary"}>
                {trainingState.is_running ? "🏃‍♂️ Running" : "⏸️ Stopped"}
              </Badge>
            </CardTitle>
          </CardHeader>
          <CardContent>
            <div className="space-y-4">
              <div>
                <div className="flex justify-between text-sm mb-2">
                  <span>Epoch {trainingState.current_epoch} / {trainingState.total_epochs}</span>
                  <span>{progressPercentage}%</span>
                </div>
                <Progress value={progressPercentage} />
              </div>
              <div className="text-sm text-gray-600">
                Duration: {formatDuration(trainingState.started_at)}
              </div>
            </div>
          </CardContent>
        </Card>

        {/* Current Metrics */}
        <Card>
          <CardHeader>
            <CardTitle>Latest Metrics</CardTitle>
            <CardDescription>
              10-epoch rolling averages 
              {trainingState?.latest_metrics && (
                <span className="text-blue-600 ml-1">
                  (Epoch {trainingState.latest_metrics.epoch})
                </span>
              )}
            </CardDescription>
          </CardHeader>
          <CardContent>
            {trainingState?.latest_metrics ? (
              <div className="space-y-3">
                <div className="grid grid-cols-2 gap-3 text-sm">
                  <div>
                    <span className="font-medium">Loss:</span>
                    <div className="text-lg font-bold text-red-600">
                      {trainingState.latest_metrics.loss.toFixed(4)}
                    </div>
                  </div>
                  <div>
                    <span className="font-medium">ARI:</span>
                    <div className="text-lg font-bold text-green-600">
                      {trainingState.latest_metrics.ari.toFixed(3)}
                    </div>
                  </div>
                  <div>
                    <span className="font-medium">NMI:</span>
                    <div className="text-blue-600">{trainingState.latest_metrics.nmi.toFixed(3)}</div>
                  </div>
                  <div>
                    <span className="font-medium">ASW:</span>
                    <div className="text-purple-600">{trainingState.latest_metrics.asw.toFixed(3)}</div>
                  </div>
                </div>
                <div className="text-xs text-gray-500 bg-blue-50 p-2 rounded">
                  📊 Smoothed values reduce noise and show clearer trends
                </div>
              </div>
            ) : (
              <p className="text-gray-500">No metrics available yet</p>
            )}
          </CardContent>
        </Card>

        {/* Controls */}
        <Card>
          <CardHeader>
            <CardTitle>Training Controls</CardTitle>
          </CardHeader>
          <CardContent className="space-y-3">
            {trainingState.is_running && (
              <Button 
                variant="destructive" 
                onClick={stopTraining}
                className="w-full"
              >
                <Square className="mr-2 h-4 w-4" />
                Stop Training
              </Button>
            )}
            
            {!trainingState.is_running && (
              <div className="space-y-2">
                <Button asChild className="w-full">
                  <Link href="/results">
                    <TrendingUp className="mr-2 h-4 w-4" />
                    View Results
                  </Link>
                </Button>
                <Button variant="outline" asChild className="w-full">
                  <Link href="/training/configure">
                    Configure New Training
                  </Link>
                </Button>
              </div>
            )}

            {trainingState.error_message && (
              <Alert variant="destructive">
                <AlertDescription>
                  {trainingState.error_message}
                </AlertDescription>
              </Alert>
            )}
          </CardContent>
        </Card>
      </div>

      {/* Individual Metric Charts */}
      {chartData.length > 0 && (
        <div className="space-y-8">
          <div className="flex items-center space-x-2">
            <h2 className="text-2xl font-bold">Training Metrics</h2>
            <Badge variant="outline">{chartData.length} data points</Badge>
          </div>

          <Tabs defaultValue="overview" className="space-y-6">
            <TabsList className="grid w-full grid-cols-4">
              <TabsTrigger value="overview">Overview</TabsTrigger>
              <TabsTrigger value="clustering">Clustering</TabsTrigger>
              <TabsTrigger value="quality">Quality</TabsTrigger>
              <TabsTrigger value="connectivity">Connectivity</TabsTrigger>
            </TabsList>

            {/* Overview Tab - Loss and Primary Metrics */}
            <TabsContent value="overview">
              <div className="grid md:grid-cols-2 gap-6">
                <MetricChart
                  data={chartData}
                  dataKey="loss"
                  title="Training Loss"
                  description="Reconstruction + KL divergence loss (lower is better)"
                  color="#dc2626"
                  formatValue={(value) => value.toFixed(4)}
                  domain={['auto', 'auto']}
                  height={300}
                />
                
                <MetricChart
                  data={chartData}
                  dataKey="ari"
                  title="Adjusted Rand Index (ARI)"
                  description="Clustering accuracy metric (0-1, higher is better)"
                  color="#16a34a"
                  formatValue={(value) => value.toFixed(3)}
                  domain={['auto', 'auto']}
                  height={300}
                />
              </div>
            </TabsContent>

            {/* Clustering Tab - ARI and NMI */}
            <TabsContent value="clustering">
              <div className="grid md:grid-cols-2 gap-6">
                <MetricChart
                  data={chartData}
                  dataKey="ari"
                  title="Adjusted Rand Index (ARI)"
                  description="Measures clustering accuracy against ground truth"
                  color="#16a34a"
                  formatValue={(value) => value.toFixed(3)}
                  domain={['auto', 'auto']}
                />
                
                <MetricChart
                  data={chartData}
                  dataKey="nmi"
                  title="Normalized Mutual Information (NMI)"
                  description="Information theoretic clustering measure"
                  color="#2563eb"
                  formatValue={(value) => value.toFixed(3)}
                  domain={['auto', 'auto']}
                />
              </div>
            </TabsContent>

            {/* Quality Tab - ASW, CH, DB */}
            <TabsContent value="quality">
              <div className="grid md:grid-cols-2 lg:grid-cols-3 gap-6">
                <MetricChart
                  data={chartData}
                  dataKey="asw"
                  title="Average Silhouette Width (ASW)"
                  description="Cluster separation quality (-1 to 1, higher is better)"
                  color="#7c3aed"
                  formatValue={(value) => value.toFixed(3)}
                  domain={['auto', 'auto']}
                />
                
                <MetricChart
                  data={chartData}
                  dataKey="ch"
                  title="Calinski-Harabasz Index (CH)"
                  description="Ratio of between/within cluster dispersion (higher is better)"
                  color="#ea580c"
                  formatValue={(value) => value.toFixed(2)}
                  domain={['auto', 'auto']}
                />
                
                <MetricChart
                  data={chartData}
                  dataKey="db"
                  title="Davies-Bouldin Index (DB)"
                  description="Average cluster separation (lower is better)"
                  color="#dc2626"
                  formatValue={(value) => value.toFixed(3)}
                  domain={['auto', 'auto']}
                />
              </div>
            </TabsContent>

            {/* Connectivity Tab - PC */}
            <TabsContent value="connectivity">
              <div className="grid md:grid-cols-1 gap-6">
                <MetricChart
                  data={chartData}
                  dataKey="pc"
                  title="Graph Connectivity Score (PC)"
                  description="Preservation of local neighborhood structure (0-1, higher is better)"
                  color="#059669"
                  formatValue={(value) => value.toFixed(3)}
                  domain={['auto', 'auto']}
                  height={350}
                />
              </div>
            </TabsContent>
          </Tabs>
        </div>
      )}

      {/* Training History Table */}
      {trainingState.history.length > 0 && (
        <Card className="mt-8">
          <CardHeader>
            <CardTitle>Training History</CardTitle>
            <CardDescription>
              Recent training metrics (last 10 epochs)
            </CardDescription>
          </CardHeader>
          <CardContent>
            <div className="overflow-x-auto">
              <table className="w-full text-sm">
                <thead>
                  <tr className="border-b">
                    <th className="text-left p-2">Epoch</th>
                    <th className="text-left p-2">Loss</th>
                    <th className="text-left p-2">ARI</th>
                    <th className="text-left p-2">NMI</th>
                    <th className="text-left p-2">ASW</th>
                    <th className="text-left p-2">CH</th>
                    <th className="text-left p-2">DB</th>
                    <th className="text-left p-2">PC</th>
                  </tr>
                </thead>
                <tbody>
                  {trainingState.history.slice(-10).reverse().map((metrics) => (
                    <tr key={metrics.epoch} className="border-b hover:bg-gray-50">
                      <td className="p-2 font-medium">{metrics.epoch}</td>
                      <td className="p-2 text-red-600">{metrics.loss.toFixed(4)}</td>
                      <td className="p-2 text-green-600">{metrics.ari.toFixed(3)}</td>
                      <td className="p-2 text-blue-600">{metrics.nmi.toFixed(3)}</td>
                      <td className="p-2 text-purple-600">{metrics.asw.toFixed(3)}</td>
                      <td className="p-2 text-orange-600">{metrics.ch.toFixed(2)}</td>
                      <td className="p-2 text-red-500">{metrics.db.toFixed(3)}</td>
                      <td className="p-2 text-green-500">{metrics.pc.toFixed(3)}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </CardContent>
        </Card>
      )}
    </div>
  );
}
