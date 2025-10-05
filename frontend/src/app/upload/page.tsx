
// src/app/upload/page.tsx
'use client';

import React from 'react';
import { ArrowLeft, ArrowRight, RefreshCcw, Upload, CheckCircle, Info } from 'lucide-react';
import Link from 'next/link';
import { Button } from '@/components/ui/button';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Alert, AlertDescription } from '@/components/ui/alert';
import { Badge } from '@/components/ui/badge';
import { FileUpload } from '@/components/shared/FileUpload';
import { DataSummaryCard } from '@/components/data/DataSummaryCard';
import { useDataUpload } from '@/lib/hooks/useDataUpload';
import { UPLOAD_CONSTRAINTS } from '@/types/index';

export default function UploadPage() {
  const {
    isUploading,
    progress,
    error,
    dataSummary,
    isSuccess,
    uploadFile,
    reset,
  } = useDataUpload();

  return (
    <div className="container mx-auto py-8 px-4 max-w-2xl">
      {/* Header */}
      <div className="flex items-center justify-between mb-6">
        <div className="text-center mb-8">
        <h1 className="text-3xl font-bold mb-2 flex items-center justify-center gap-3">
          <Upload className="h-8 w-8" />
          Upload Dataset
        </h1>
        <p className="text-gray-600">
          Upload your scRNA-seq data in AnnData format to begin training
        </p>
      </div>
      <Button variant="outline" asChild>
          <Link href="/">
            <ArrowLeft className="mr-2 h-4 w-4" />
            Home
          </Link>
      </Button>
      </div>
      

      {!isSuccess ? (
        <div className="max-w-2xl mx-auto space-y-6">
          {/* Requirements Card */}
          <Card>
            <CardHeader className="pb-3">
              <CardTitle className="text-lg">File Requirements</CardTitle>
            </CardHeader>
            <CardContent>
              <div className="grid sm:grid-cols-2 gap-3 text-sm">
                <div className="flex items-center gap-2">
                  <Badge variant="outline">Format</Badge>
                  <span className="text-gray-700">.h5ad (AnnData)</span>
                </div>
                <div className="flex items-center gap-2">
                  <Badge variant="outline">Max Size</Badge>
                  <span className="text-gray-700">
                    {Math.round(UPLOAD_CONSTRAINTS.maxSizeBytes / (1024 * 1024))} MB
                  </span>
                </div>
                <div className="flex items-center gap-2">
                  <Badge variant="outline">Data</Badge>
                  <span className="text-gray-700">Raw or processed counts</span>
                </div>
                <div className="flex items-center gap-2">
                  <Badge variant="outline">Metadata</Badge>
                  <span className="text-gray-700">Optional annotations</span>
                </div>
              </div>
            </CardContent>
          </Card>

          {/* Upload Component */}
          <FileUpload
            onFileSelect={uploadFile}
            isUploading={isUploading}
            progress={progress}
            error={error}
            isSuccess={isSuccess}
          />

          {/* Quick Help */}
          <Alert>
            <Info className="h-4 w-4" />
            <AlertDescription className="text-sm">
              Need help? Check{' '}
              <a
                href="https://anndata.readthedocs.io/"
                target="_blank"
                rel="noopener noreferrer"
                className="underline font-medium"
              >
                AnnData documentation
              </a>
              {' '}for file format details.
            </AlertDescription>
          </Alert>
        </div>
      ) : (
        <div className="space-y-6">
          {/* Success Banner */}
          <Alert className="bg-green-50 border-green-200">
            <CheckCircle className="h-5 w-5 text-green-600" />
            <AlertDescription className="text-green-800">
              <strong>Upload successful!</strong> Your dataset is ready for training.
            </AlertDescription>
          </Alert>

          {/* Data Summary */}
          <DataSummaryCard dataSummary={dataSummary!} />

          {/* Actions */}
          <div className="flex flex-wrap gap-3 justify-center">
            <Button onClick={reset} variant="outline">
              <RefreshCcw className="mr-2 h-4 w-4" />
              Upload Different File
            </Button>
            
            <Button asChild size="lg">
              <Link href="/training/configure">
                Configure Training
                <ArrowRight className="ml-2 h-4 w-4" />
              </Link>
            </Button>
          </div>
        </div>
      )}
    </div>
  );
}
