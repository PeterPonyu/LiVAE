
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
import { useTrainingMonitor } from '@/lib/hooks/useTrainingMonitor';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer } from 'recharts';

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
          </CardHeader>
          <CardContent>
            {trainingState.latest_metrics ? (
              <div className="grid grid-cols-2 gap-3 text-sm">
                <div>
                  <span className="font-medium">Loss:</span>
                  <div className="text-lg font-bold text-blue-600">
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
                  <div className="text-blue-600">{trainingState.latest_metrics.asw.toFixed(3)}</div>
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

      {/* Training Metrics Chart */}
      {chartData.length > 0 && (
        <Card className="mb-8">
          <CardHeader>
            <CardTitle>Training Metrics Over Time</CardTitle>
            <CardDescription>
              Loss and clustering quality metrics during training
            </CardDescription>
          </CardHeader>
          <CardContent>
            <div className="h-96">
              <ResponsiveContainer width="100%" height="100%">
                <LineChart data={chartData}>
                  <CartesianGrid strokeDasharray="3 3" />
                  <XAxis dataKey="epoch" />
                  <YAxis />
                  <Tooltip />
                  <Legend />
                  <Line type="monotone" dataKey="loss" stroke="#8884d8" name="Loss" />
                  <Line type="monotone" dataKey="ari" stroke="#82ca9d" name="ARI" />
                  <Line type="monotone" dataKey="nmi" stroke="#ffc658" name="NMI" />
                  <Line type="monotone" dataKey="asw" stroke="#ff7300" name="ASW" />
                </LineChart>
              </ResponsiveContainer>
            </div>
          </CardContent>
        </Card>
      )}

      {/* Training History Table */}
      {trainingState.history.length > 0 && (
        <Card>
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
                    <tr key={metrics.epoch} className="border-b">
                      <td className="p-2 font-medium">{metrics.epoch}</td>
                      <td className="p-2">{metrics.loss.toFixed(4)}</td>
                      <td className="p-2">{metrics.ari.toFixed(3)}</td>
                      <td className="p-2">{metrics.nmi.toFixed(3)}</td>
                      <td className="p-2">{metrics.asw.toFixed(3)}</td>
                      <td className="p-2">{metrics.ch.toFixed(2)}</td>
                      <td className="p-2">{metrics.db.toFixed(3)}</td>
                      <td className="p-2">{metrics.pc.toFixed(3)}</td>
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
