
// src/app/upload/page.tsx
'use client';

import React from 'react';
import { ArrowRight, RefreshCcw } from 'lucide-react';
import Link from 'next/link';
import { Button } from '@/components/ui/button';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
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

  const handleFileSelect = (file: File) => {
    // Validate file before upload
    if (!file.name.endsWith('.h5ad')) {
      // This should be handled by react-dropzone, but just in case
      return;
    }

    if (file.size > UPLOAD_CONSTRAINTS.maxSizeBytes) {
      // This should also be handled by react-dropzone
      return;
    }

    uploadFile(file);
  };

  const handleStartOver = () => {
    reset();
  };

  return (
    <div className="container mx-auto py-8 px-4 max-w-6xl">
      {/* Page Header */}
      <div className="text-center mb-8">
        <h1 className="text-3xl font-bold mb-2">Upload Your Dataset</h1>
        <p className="text-gray-600 max-w-2xl mx-auto">
          Upload your single-cell RNA sequencing data in AnnData (.h5ad) format 
          to begin analysis with the LiVAE deep learning model.
        </p>
      </div>

      {!isSuccess ? (
        <div className="max-w-2xl mx-auto">
          {/* Upload Instructions */}
          <Card className="mb-6">
            <CardHeader>
              <CardTitle>Upload Requirements</CardTitle>
              <CardDescription>
                Please ensure your file meets these requirements
              </CardDescription>
            </CardHeader>
            <CardContent>
              <ul className="space-y-2 text-sm">
                <li className="flex items-center">
                  <span className="w-2 h-2 bg-blue-500 rounded-full mr-2"></span>
                  File format: AnnData (.h5ad)
                </li>
                <li className="flex items-center">
                  <span className="w-2 h-2 bg-blue-500 rounded-full mr-2"></span>
                  Maximum size: {Math.round(UPLOAD_CONSTRAINTS.maxSizeBytes / (1024 * 1024))} MB
                </li>
                <li className="flex items-center">
                  <span className="w-2 h-2 bg-blue-500 rounded-full mr-2"></span>
                  Should contain raw or processed count data
                </li>
                <li className="flex items-center">
                  <span className="w-2 h-2 bg-blue-500 rounded-full mr-2"></span>
                  Cell and gene annotations are optional but recommended
                </li>
              </ul>
            </CardContent>
          </Card>

          {/* File Upload Component */}
          <FileUpload
            onFileSelect={handleFileSelect}
            isUploading={isUploading}
            progress={progress}
            error={error}
            isSuccess={isSuccess}
          />
        </div>
      ) : (
        <div className="space-y-6">
          {/* Success Header */}
          <div className="text-center bg-green-50 border border-green-200 rounded-lg p-6">
            <h2 className="text-2xl font-bold text-green-800 mb-2">
              Dataset Loaded Successfully!
            </h2>
            <p className="text-green-600">
              Your data has been processed and is ready for analysis.
            </p>
          </div>

          {/* Data Summary */}
          <DataSummaryCard dataSummary={dataSummary!} />

          {/* Action Buttons */}
          <div className="flex flex-col sm:flex-row gap-4 justify-center pt-6">
            <Button onClick={handleStartOver} variant="outline">
              <RefreshCcw className="mr-2 h-4 w-4" />
              Upload Different File
            </Button>
            
            <Button asChild>
              <Link href="/training/configure">
                <ArrowRight className="mr-2 h-4 w-4" />
                Configure Training
              </Link>
            </Button>
          </div>
        </div>
      )}

      {/* Help Section */}
      <Card className="mt-12 bg-blue-50 border-blue-200">
        <CardHeader>
          <CardTitle className="text-blue-800">Need Help?</CardTitle>
        </CardHeader>
        <CardContent className="text-blue-700">
          <p className="mb-2">
            If you're having trouble uploading your file, make sure:
          </p>
          <ul className="list-disc list-inside space-y-1 text-sm">
            <li>The file is in AnnData (.h5ad) format</li>
            <li>The file size is under the limit</li>
            <li>You have a stable internet connection</li>
            <li>The backend server is running</li>
          </ul>
          <p className="mt-3 text-sm">
            For more information about AnnData format, visit the{' '}
            <a
              href="https://anndata.readthedocs.io/"
              target="_blank"
              rel="noopener noreferrer"
              className="underline font-medium"
            >
              AnnData documentation
            </a>
          </p>
        </CardContent>
      </Card>
    </div>
  );
}
