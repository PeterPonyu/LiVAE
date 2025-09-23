
// src/components/shared/FileUpload.tsx
'use client';

import React, { useCallback } from 'react';
import { useDropzone } from 'react-dropzone';
import { Upload, File, AlertCircle, CheckCircle } from 'lucide-react';
import { Card } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Alert, AlertDescription } from '@/components/ui/alert';
import { Progress } from '@/components/ui/progress';
import { UPLOAD_CONSTRAINTS } from '@/types/index';

interface FileUploadProps {
  onFileSelect: (file: File) => void;
  isUploading?: boolean;
  progress?: number;
  error?: string | null;
  isSuccess?: boolean;
  disabled?: boolean;
}

export const FileUpload: React.FC<FileUploadProps> = ({
  onFileSelect,
  isUploading = false,
  progress = 0,
  error = null,
  isSuccess = false,
  disabled = false,
}) => {
  const onDrop = useCallback((acceptedFiles: File[], rejectedFiles: any[]) => {
    if (rejectedFiles.length > 0) {
      // Handle rejected files - will be shown as error in parent
      return;
    }
    
    if (acceptedFiles.length > 0) {
      onFileSelect(acceptedFiles[0]);
    }
  }, [onFileSelect]);

  const {
    getRootProps,
    getInputProps,
    isDragActive,
    isDragReject,
  } = useDropzone({
    onDrop,
    accept: {
      'application/octet-stream': ['.h5ad'],
    },
    maxSize: UPLOAD_CONSTRAINTS.maxSizeBytes,
    multiple: false,
    disabled: disabled || isUploading,
  });

  const formatFileSize = (bytes: number): string => {
    const mb = bytes / (1024 * 1024);
    return `${mb.toFixed(1)} MB`;
  };

  return (
    <div className="w-full space-y-4">
      <Card
        {...getRootProps()}
        className={`
          relative cursor-pointer border-2 border-dashed p-8 text-center transition-colors
          ${isDragActive && !isDragReject ? 'border-blue-500 bg-blue-50' : ''}
          ${isDragReject ? 'border-red-500 bg-red-50' : ''}
          ${isUploading || disabled ? 'cursor-not-allowed opacity-50' : ''}
          ${isSuccess ? 'border-green-500 bg-green-50' : 'border-gray-300'}
          hover:border-gray-400 hover:bg-gray-50
        `}
      >
        <input {...getInputProps()} />
        
        <div className="flex flex-col items-center space-y-4">
          {isSuccess ? (
            <CheckCircle className="h-12 w-12 text-green-500" />
          ) : (
            <Upload className="h-12 w-12 text-gray-400" />
          )}
          
          <div className="space-y-2">
            <h3 className="text-lg font-semibold">
              {isSuccess ? 'File Uploaded Successfully!' : 'Upload AnnData File'}
            </h3>
            
            {!isSuccess && (
              <div className="text-sm text-gray-600">
                <p>
                  {isDragActive 
                    ? 'Drop the file here...' 
                    : 'Drag & drop your .h5ad file here, or click to browse'
                  }
                </p>
                <p className="mt-1">
                  Maximum file size: {formatFileSize(UPLOAD_CONSTRAINTS.maxSizeBytes)}
                </p>
              </div>
            )}
          </div>

          {!isSuccess && !isUploading && (
            <Button variant="outline" size="sm">
              <File className="mr-2 h-4 w-4" />
              Choose File
            </Button>
          )}
        </div>
      </Card>

      {/* Upload Progress */}
      {isUploading && (
        <div className="space-y-2">
          <div className="flex items-center justify-between text-sm">
            <span>Uploading...</span>
            <span>{progress}%</span>
          </div>
          <Progress value={progress} />
        </div>
      )}

      {/* Error Display */}
      {error && (
        <Alert variant="destructive">
          <AlertCircle className="h-4 w-4" />
          <AlertDescription>
            {error}
          </AlertDescription>
        </Alert>
      )}
    </div>
  );
};
