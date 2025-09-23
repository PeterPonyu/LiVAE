
// src/app/page.tsx
'use client';

import { useState, useEffect } from 'react';
import { Upload, Activity, Download, Settings } from 'lucide-react';
import Link from 'next/link';
import { Button } from '@/components/ui/button';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import { api } from '@/lib/api/endpoints';
import type { AppStatus } from '@/types/index';

export default function Home() {
  const [status, setStatus] = useState<AppStatus | null>(null);
  const [loading, setLoading] = useState(false);

  useEffect(() => {
    // Check backend status on page load
    const checkStatus = async () => {
      try {
        const result = await api.getStatus();
        setStatus(result);
      } catch (error) {
        console.error('Failed to check status:', error);
      }
    };
    
    checkStatus();
  }, []);

  const testConnection = async () => {
    setLoading(true);
    try {
      const result = await api.getStatus();
      setStatus(result);
    } catch (error) {
      console.error('Connection failed:', error);
      setStatus(null);
    } finally {
      setLoading(false);
    }
  };

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
            <Button 
              onClick={testConnection} 
              disabled={loading}
              variant="outline"
              size="sm"
            >
              {loading ? 'Checking...' : 'Refresh'}
            </Button>
          </CardTitle>
        </CardHeader>
        <CardContent>
          {status ? (
            <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
              <div className="text-center">
                <Badge variant={status.data_loaded ? "default" : "secondary"}>
                  {status.data_loaded ? "Loaded" : "No Data"}
                </Badge>
                <div className="text-sm text-gray-600 mt-1">Dataset</div>
              </div>
              <div className="text-center">
                <Badge variant={status.model_trained ? "default" : "secondary"}>
                  {status.model_trained ? "Trained" : "Not Trained"}
                </Badge>
                <div className="text-sm text-gray-600 mt-1">Model</div>
              </div>
              <div className="text-center">
                <Badge variant={status.training_running ? "destructive" : "secondary"}>
                  {status.training_running ? "Running" : "Idle"}
                </Badge>
                <div className="text-sm text-gray-600 mt-1">Training</div>
              </div>
              <div className="text-center">
                <Badge variant="outline">
                  {status.device}
                </Badge>
                <div className="text-sm text-gray-600 mt-1">Device</div>
              </div>
            </div>
          ) : (
            <div className="text-center text-gray-600">
              <p>Unable to connect to backend server</p>
              <p className="text-sm mt-1">Make sure the server is running on http://127.0.0.1:8000</p>
            </div>
          )}
        </CardContent>
      </Card>

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
            <Button asChild className="w-full" variant={status?.data_loaded ? "outline" : "default"}>
              <Link href="/upload">
                {status?.data_loaded ? "Replace Data" : "Upload Data"}
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
              disabled={!status?.data_loaded}
            >
              <Link href="/training/configure">
                Configure Training
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
                <CardTitle className="text-lg">Train Model</CardTitle>
                <CardDescription>Step 3</CardDescription>
              </div>
            </div>
          </CardHeader>
          <CardContent>
            <p className="text-sm text-gray-600 mb-4">
              Monitor training progress and metrics in real-time
            </p>
            <Button 
              asChild 
              className="w-full" 
              variant="outline"
              disabled={!status?.data_loaded}
            >
              <Link href="/training/monitor">
                {status?.training_running ? "View Training" : "Start Training"}
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
              disabled={!status?.model_trained}
            >
              <Link href="/results">
                Download Results
              </Link>
            </Button>
          </CardContent>
        </Card>
      </div>

      {/* Quick Actions */}
      {status?.data_loaded && (
        <Card className="mt-8 bg-gradient-to-r from-blue-50 to-green-50 border-blue-200">
          <CardContent className="pt-6">
            <div className="flex flex-col sm:flex-row items-center justify-between">
              <div>
                <h3 className="font-semibold text-blue-800 mb-1">Ready to Continue</h3>
                <p className="text-blue-600 text-sm">
                  Your dataset is loaded and ready for training configuration.
                </p>
              </div>
              <Button asChild className="mt-4 sm:mt-0">
                <Link href="/training/configure">
                  Configure Training →
                </Link>
              </Button>
            </div>
          </CardContent>
        </Card>
      )}
    </div>
  );
}
