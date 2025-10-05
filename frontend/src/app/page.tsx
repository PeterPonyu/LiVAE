
// src/app/page.tsx
'use client';

import React from 'react';
import { Upload, Activity, Download, Settings, RefreshCw } from 'lucide-react';
import Link from 'next/link';
import { Button } from '@/components/ui/button';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import { useAppStatus } from '@/lib/hooks/useAppStatus';

export default function Home() {
  const { 
    data_loaded, 
    model_trained, 
    training_running, 
    device, 
    is_loading, 
    connection_error,
    last_updated,
    refresh 
  } = useAppStatus();

  // Determine the actual state more accurately
  const getTrainingStatus = () => {
    if (training_running) return 'running';
    if (model_trained && !training_running) return 'completed';
    if (data_loaded) return 'ready';
    return 'no-data';
  };

  const trainingStatus = getTrainingStatus();

  return (
    <div className="container mx-auto p-8 max-w-6xl">
      {/* Header */}
      <div className="text-center mb-12">
        <h1 className="text-4xl font-bold mb-4">Single Cell Deep Learning Agent</h1>
        <p className="text-xl text-gray-600 max-w-3xl mx-auto">
          Train interpretable variational autoencoders on your single-cell RNA sequencing data 
          using the LiVAE framework
        </p>
      </div>

      {/* Status Card */}
      <Card className="mb-8">
        <CardHeader>
          <CardTitle className="flex items-center justify-between">
            System Status
            <div className="flex items-center space-x-2">
              <span className="text-sm text-gray-500">
                Updated: {last_updated.toLocaleTimeString()}
              </span>
              <Button 
                onClick={refresh} 
                disabled={is_loading}
                variant="outline"
                size="sm"
              >
                <RefreshCw className={`h-4 w-4 ${is_loading ? 'animate-spin' : ''}`} />
              </Button>
            </div>
          </CardTitle>
        </CardHeader>
        <CardContent>
          {connection_error ? (
            <div className="text-center text-red-600">
              <p>❌ Unable to connect to backend server</p>
              <p className="text-sm mt-1">Make sure the server is running on http://127.0.0.1:8000</p>
            </div>
          ) : (
            <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
              <div className="text-center">
                <Badge variant={data_loaded ? "default" : "secondary"}>
                  {data_loaded ? "✅ Loaded" : "📋 No Data"}
                </Badge>
                <div className="text-sm text-gray-600 mt-1">Dataset</div>
              </div>
              
              <div className="text-center">
                <Badge variant={
                  trainingStatus === 'completed' ? "default" :
                  trainingStatus === 'running' ? "destructive" :
                  trainingStatus === 'ready' ? "secondary" : "outline"
                }>
                  {trainingStatus === 'completed' ? "✅ Trained" :
                   trainingStatus === 'running' ? "🔄 Training" :
                   trainingStatus === 'ready' ? "⚡ Ready" : "⏸️ Not Started"}
                </Badge>
                <div className="text-sm text-gray-600 mt-1">Model</div>
              </div>
              
              <div className="text-center">
                <Badge variant={training_running ? "destructive" : "secondary"}>
                  {training_running ? "🏃‍♂️ Running" : "💤 Idle"}
                </Badge>
                <div className="text-sm text-gray-600 mt-1">Training Process</div>
              </div>
              
              <div className="text-center">
                <Badge variant="outline">
                  {device === 'cuda' ? '🚀 GPU' : 
                   device === 'mps' ? '🍎 Metal' : '🖥️ CPU'}
                </Badge>
                <div className="text-sm text-gray-600 mt-1">Device</div>
              </div>
            </div>
          )}
        </CardContent>
      </Card>

      {/* Training Status Alert */}
      {training_running && (
        <Card className="mb-8 bg-blue-50 border-blue-200">
          <CardContent className="pt-6">
            <div className="flex flex-col sm:flex-row items-center justify-between">
              <div className="flex items-center space-x-3">
                <div className="animate-spin rounded-full h-6 w-6 border-b-2 border-blue-600"></div>
                <div>
                  <h3 className="font-semibold text-blue-800">Training in Progress</h3>
                  <p className="text-blue-600 text-sm">
                    Your model is currently training. Monitor progress in real-time.
                  </p>
                </div>
              </div>
              <Button asChild className="mt-4 sm:mt-0">
                <Link href="/training/monitor">
                  View Progress →
                </Link>
              </Button>
            </div>
          </CardContent>
        </Card>
      )}

      {/* Main Workflow */}
      <div className="grid md:grid-cols-2 lg:grid-cols-4 gap-6">
        {/* Step 1: Upload Data */}
        <Card className="group hover:shadow-lg transition-shadow">
          <CardHeader>
            <div className="flex items-center space-x-2">
              <div className="p-2 bg-blue-100 rounded-lg">
                <Upload className="h-6 w-6 text-blue-600" />
              </div>
              <div className="flex-1">
                <CardTitle className="text-lg">Upload Data</CardTitle>
                <CardDescription>Step 1</CardDescription>
              </div>
            </div>
          </CardHeader>
          <CardContent>
            <p className="text-sm text-gray-600 mb-4">
              Upload your single-cell RNA-seq data in AnnData (.h5ad) format
            </p>
            <Button asChild className="w-full" variant={data_loaded ? "outline" : "default"}>
              <Link href="/upload">
                {data_loaded ? "Replace Data" : "Upload Data"}
              </Link>
            </Button>
          </CardContent>
        </Card>

        {/* Step 2: Configure Training */}
        <Card className="group hover:shadow-lg transition-shadow">
          <CardHeader>
            <div className="flex items-center space-x-2">
              <div className="p-2 bg-green-100 rounded-lg">
                <Settings className="h-6 w-6 text-green-600" />
              </div>
              <div className="flex-1">
                <CardTitle className="text-lg">Configure</CardTitle>
                <CardDescription>Step 2</CardDescription>
              </div>
            </div>
          </CardHeader>
          <CardContent>
            <p className="text-sm text-gray-600 mb-4">
              Set model parameters, training configuration, and QC filters
            </p>
            <Button 
              asChild 
              className="w-full" 
              variant="outline"
              disabled={!data_loaded || training_running}
            >
              <Link href="/training/configure">
                {training_running ? "Training Running..." : "Configure Training"}
              </Link>
            </Button>
          </CardContent>
        </Card>

        {/* Step 3: Monitor Training */}
        <Card className="group hover:shadow-lg transition-shadow">
          <CardHeader>
            <div className="flex items-center space-x-2">
              <div className="p-2 bg-orange-100 rounded-lg">
                <Activity className="h-6 w-6 text-orange-600" />
              </div>
              <div className="flex-1">
                <CardTitle className="text-lg">Train & Monitor</CardTitle>
                <CardDescription>Step 3</CardDescription>
              </div>
            </div>
          </CardHeader>
          <CardContent>
            <p className="text-sm text-gray-600 mb-4">
              {training_running 
                ? "Monitor training progress and metrics in real-time"
                : "Start training and monitor progress"
              }
            </p>
            <Button 
              asChild 
              className="w-full" 
              variant={training_running ? "default" : "outline"}
              disabled={!data_loaded}
            >
              <Link href="/training/monitor">
                {training_running ? "View Progress" : 
                 trainingStatus === 'completed' ? "View Results" : "Start Training"}
              </Link>
            </Button>
          </CardContent>
        </Card>

        {/* Step 4: Download Results */}
        <Card className="group hover:shadow-lg transition-shadow">
          <CardHeader>
            <div className="flex items-center space-x-2">
              <div className="p-2 bg-purple-100 rounded-lg">
                <Download className="h-6 w-6 text-purple-600" />
              </div>
              <div className="flex-1">
                <CardTitle className="text-lg">Get Results</CardTitle>
                <CardDescription>Step 4</CardDescription>
              </div>
            </div>
          </CardHeader>
          <CardContent>
            <p className="text-sm text-gray-600 mb-4">
              Download trained embeddings and visualization results
            </p>
            <Button 
              asChild 
              className="w-full" 
              variant="outline"
              disabled={trainingStatus !== 'completed'}
            >
              <Link href="/results">
                Download Results
              </Link>
            </Button>
          </CardContent>
        </Card>
      </div>
    </div>
  );
}
