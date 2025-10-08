
// src/app/page.tsx
'use client';

import React from 'react';
import { useState, useEffect } from 'react';
import { Upload, Activity, Download, Settings, RefreshCw } from 'lucide-react';
import Link from 'next/link';
import { useAppStatus } from '@/lib/hooks/useAppStatus';

// Custom Button Component
interface ButtonProps extends React.ButtonHTMLAttributes<HTMLButtonElement> {
  variant?: 'default' | 'outline' | 'destructive' | 'secondary';
  size?: 'sm' | 'md' | 'lg';
  asChild?: boolean;
  children: React.ReactNode;
}

const Button: React.FC<ButtonProps> = ({ 
  variant = 'default', 
  size = 'md',
  className = '', 
  disabled,
  children,
  ...props 
}) => {
  const baseStyles = 'inline-flex items-center justify-center rounded-md font-medium transition-colors focus:outline-none focus:ring-2 focus:ring-offset-2 disabled:opacity-50 disabled:pointer-events-none';
  
  const variants = {
    default: 'bg-blue-600 text-white hover:bg-blue-700 focus:ring-blue-500',
    outline: 'border border-gray-300 bg-white text-gray-700 hover:bg-gray-50 focus:ring-blue-500',
    destructive: 'bg-red-600 text-white hover:bg-red-700 focus:ring-red-500',
    secondary: 'bg-gray-200 text-gray-900 hover:bg-gray-300 focus:ring-gray-500'
  };
  
  const sizes = {
    sm: 'px-3 py-1.5 text-sm',
    md: 'px-4 py-2 text-sm',
    lg: 'px-6 py-3 text-base'
  };
  
  return (
    <button
      className={`${baseStyles} ${variants[variant]} ${sizes[size]} ${className}`}
      disabled={disabled}
      {...props}
    >
      {children}
    </button>
  );
};

// Custom Badge Component
interface BadgeProps {
  variant?: 'default' | 'secondary' | 'outline' | 'destructive';
  children: React.ReactNode;
  className?: string;
}

const Badge: React.FC<BadgeProps> = ({ variant = 'default', children, className = '' }) => {
  const variants = {
    default: 'bg-blue-600 text-white',
    secondary: 'bg-gray-200 text-gray-900',
    outline: 'border border-gray-300 bg-white text-gray-700',
    destructive: 'bg-red-600 text-white'
  };
  
  return (
    <span className={`inline-flex items-center px-2.5 py-0.5 rounded-full text-xs font-medium ${variants[variant]} ${className}`}>
      {children}
    </span>
  );
};

// Custom Card Components
const Card: React.FC<{ children: React.ReactNode; className?: string }> = ({ children, className = '' }) => (
  <div className={`bg-white rounded-lg border border-gray-200 shadow-sm ${className}`}>
    {children}
  </div>
);

const CardHeader: React.FC<{ children: React.ReactNode; className?: string }> = ({ children, className = '' }) => (
  <div className={`px-6 py-5 border-b border-gray-200 ${className}`}>
    {children}
  </div>
);

const CardTitle: React.FC<{ children: React.ReactNode; className?: string }> = ({ children, className = '' }) => (
  <h3 className={`text-lg font-semibold text-gray-900 ${className}`}>
    {children}
  </h3>
);

const CardDescription: React.FC<{ children: React.ReactNode; className?: string }> = ({ children, className = '' }) => (
  <p className={`text-sm text-gray-500 mt-1 ${className}`}>
    {children}
  </p>
);

const CardContent: React.FC<{ children: React.ReactNode; className?: string }> = ({ children, className = '' }) => (
  <div className={`px-6 py-5 ${className}`}>
    {children}
  </div>
);

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
  const [mounted, setMounted] = useState(false);
  useEffect(() => {
    setMounted(true);
  }, []);

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
                Updated: {mounted ? last_updated.toLocaleTimeString(): "--:--:--"}
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
              <Link href="/training/monitor" className="mt-4 sm:mt-0">
                <Button>
                  View Progress →
                </Button>
              </Link>
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
            <Link href="/upload" className="block">
              <Button className="w-full" variant={data_loaded ? "outline" : "default"}>
                {data_loaded ? "Replace Data" : "Upload Data"}
              </Button>
            </Link>
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
            <Link href="/training/configure" className={`block ${(!data_loaded || training_running) ? 'pointer-events-none' : ''}`}>
              <Button 
                className="w-full" 
                variant="outline"
                disabled={!data_loaded || training_running}
              >
                {training_running ? "Training Running..." : "Configure Training"}
              </Button>
            </Link>
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
            <Link href="/training/monitor" className={`block ${!data_loaded ? 'pointer-events-none' : ''}`}>
              <Button 
                className="w-full" 
                variant={training_running ? "default" : "outline"}
                disabled={!data_loaded}
              >
                {training_running ? "View Progress" : 
                 trainingStatus === 'completed' ? "View Results" : "Start Training"}
              </Button>
            </Link>
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
            <Link href="/results" className={`block ${trainingStatus !== 'completed' ? 'pointer-events-none' : ''}`}>
              <Button 
                className="w-full" 
                variant="outline"
                disabled={trainingStatus !== 'completed'}
              >
                Download Results
              </Button>
            </Link>
          </CardContent>
        </Card>
      </div>
    </div>
  );
}
