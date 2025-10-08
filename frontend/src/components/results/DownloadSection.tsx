// ========================================
// src/components/results/DownloadSection.tsx
// ========================================
import React from 'react';
import { Download, FileText, Layers, RefreshCw, Package, CheckCircle, AlertCircle } from 'lucide-react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle, Button, Badge, Progress, Alert, AlertDescription } from './ui';
import type { EmbeddingResult } from '@/types/index';

interface DownloadSectionProps {
  availableEmbeddings: EmbeddingResult[];
  isGeneratingEmbeddings: boolean;
  downloadProgress: Record<string, number>;
  downloadStatus: Record<string, 'idle' | 'downloading' | 'complete' | 'error'>;
  onDownload: (embeddingType: 'interpretable' | 'latent') => void;
  onDownloadAll: () => void;
  onGenerate: () => void;
  trainingCompleted: boolean;
}

export const DownloadSection: React.FC<DownloadSectionProps> = ({
  availableEmbeddings,
  isGeneratingEmbeddings,
  downloadProgress,
  downloadStatus,
  onDownload,
  onDownloadAll,
  onGenerate,
  trainingCompleted,
}) => {
  const interpretableEmbedding = availableEmbeddings.find(e => e.embedding_type === 'interpretable');
  const latentEmbedding = availableEmbeddings.find(e => e.embedding_type === 'latent');

  const getStatusIcon = (status: string) => {
    switch (status) {
      case 'complete':
        return <CheckCircle className="h-4 w-4 text-green-600" />;
      case 'error':
        return <AlertCircle className="h-4 w-4 text-red-600" />;
      case 'downloading':
        return <RefreshCw className="h-4 w-4 animate-spin text-blue-600" />;
      default:
        return null;
    }
  };

  const getStatusBadge = (status: string) => {
    switch (status) {
      case 'complete':
        return <Badge variant="default" className="bg-green-600">✅ Downloaded</Badge>;
      case 'error':
        return <Badge variant="destructive">❌ Error</Badge>;
      case 'downloading':
        return <Badge variant="secondary">⏳ Downloading...</Badge>;
      default:
        return null;
    }
  };

  const isDownloading = Object.values(downloadStatus).some(status => status === 'downloading');

  if (!trainingCompleted) {
    return (
      <Card>
        <CardHeader>
          <CardTitle>Download Results</CardTitle>
          <CardDescription>Training must be completed first</CardDescription>
        </CardHeader>
        <CardContent>
          <Alert>
            <AlertDescription>
              Complete the training process to download embeddings and results.
            </AlertDescription>
          </Alert>
        </CardContent>
      </Card>
    );
  }

  return (
    <Card>
      <CardHeader>
        <CardTitle className="flex items-center space-x-2">
          <Download className="h-5 w-5" />
          <span>Download Results</span>
        </CardTitle>
        <CardDescription>
          Download trained embeddings and model outputs
        </CardDescription>
      </CardHeader>
      <CardContent className="space-y-6">
        {/* Generate Embeddings */}
        {availableEmbeddings.length === 0 && (
          <div className="text-center py-6">
            <div className="space-y-4">
              <div className="text-gray-500">
                <Layers className="h-8 w-8 mx-auto mb-2" />
                <p>Embeddings need to be generated first</p>
              </div>
              <Button onClick={onGenerate} disabled={isGeneratingEmbeddings || isDownloading}>
                <RefreshCw className={`mr-2 h-4 w-4 ${isGeneratingEmbeddings ? 'animate-spin' : ''}`} />
                {isGeneratingEmbeddings ? 'Generating...' : 'Generate Embeddings'}
              </Button>
            </div>
          </div>
        )}

        {/* Available Downloads */}
        {availableEmbeddings.length > 0 && (
          <div className="space-y-4">
            <div className="flex items-center justify-between">
              <h4 className="font-medium">Available Downloads</h4>
              <Button 
                variant="outline" 
                size="sm"
                onClick={onDownloadAll}
                disabled={isDownloading}
                className="flex items-center space-x-2"
              >
                <Package className="h-4 w-4" />
                <span>Download All (ZIP)</span>
                {downloadStatus.all && getStatusIcon(downloadStatus.all)}
              </Button>
            </div>
            
            {/* Download All Progress */}
            {downloadStatus.all === 'downloading' && (
              <div className="space-y-2">
                <Progress value={downloadProgress.all || 0} />
                <p className="text-xs text-center text-gray-600">
                  Preparing ZIP file... {Math.round(downloadProgress.all || 0)}%
                </p>
              </div>
            )}

            {/* Interpretable Embedding */}
            {interpretableEmbedding && (
              <Card className="border-green-200 bg-green-50">
                <CardContent className="p-4">
                  <div className="flex items-center justify-between">
                    <div className="flex items-center space-x-3">
                      <div className="p-2 bg-green-100 rounded-lg">
                        <FileText className="h-5 w-5 text-green-600" />
                      </div>
                      <div className="flex-1">
                        <h5 className="font-medium">Interpretable Embedding</h5>
                        <p className="text-sm text-gray-600">
                          2D embedding for visualization and interpretation
                        </p>
                        {interpretableEmbedding.extracted_at && (
                          <p className="text-xs text-green-600">
                            Generated: {new Date(interpretableEmbedding.extracted_at).toLocaleString()}
                          </p>
                        )}
                        {interpretableEmbedding.shape && (
                          <p className="text-xs text-gray-500">
                            Shape: {interpretableEmbedding.shape.n_obs} cells × {interpretableEmbedding.shape.n_vars || 2} dimensions
                          </p>
                        )}
                      </div>
                    </div>
                    <div className="flex flex-col items-end space-y-2">
                      <Badge variant="outline">CSV Format</Badge>
                      
                      {downloadStatus.interpretable === 'downloading' ? (
                        <div className="w-32 space-y-1">
                          <Progress value={downloadProgress.interpretable || 0} />
                          <p className="text-xs text-center">{Math.round(downloadProgress.interpretable || 0)}%</p>
                        </div>
                      ) : (
                        <>
                          {getStatusBadge(downloadStatus.interpretable || 'idle')}
                          {downloadStatus.interpretable !== 'complete' && (
                            <Button 
                              size="sm" 
                              onClick={() => onDownload('interpretable')}
                              disabled={isDownloading}
                              className="flex items-center space-x-1"
                            >
                              <Download className="h-3 w-3" />
                              <span>Download</span>
                            </Button>
                          )}
                        </>
                      )}
                    </div>
                  </div>
                </CardContent>
              </Card>
            )}

            {/* Latent Embedding */}
            {latentEmbedding && (
              <Card className="border-blue-200 bg-blue-50">
                <CardContent className="p-4">
                  <div className="flex items-center justify-between">
                    <div className="flex items-center space-x-3">
                      <div className="p-2 bg-blue-100 rounded-lg">
                        <Layers className="h-5 w-5 text-blue-600" />
                      </div>
                      <div className="flex-1">
                        <h5 className="font-medium">Latent Embedding</h5>
                        <p className="text-sm text-gray-600">
                          High-dimensional latent space representation
                        </p>
                        {latentEmbedding.extracted_at && (
                          <p className="text-xs text-blue-600">
                            Generated: {new Date(latentEmbedding.extracted_at).toLocaleString()}
                          </p>
                        )}
                        {latentEmbedding.shape && (
                          <p className="text-xs text-gray-500">
                            Shape: {latentEmbedding.shape.n_obs} cells × {latentEmbedding.shape.n_vars} dimensions
                          </p>
                        )}
                      </div>
                    </div>
                    <div className="flex flex-col items-end space-y-2">
                      <Badge variant="outline">CSV Format</Badge>
                      
                      {downloadStatus.latent === 'downloading' ? (
                        <div className="w-32 space-y-1">
                          <Progress value={downloadProgress.latent || 0} />
                          <p className="text-xs text-center">{Math.round(downloadProgress.latent || 0)}%</p>
                        </div>
                      ) : (
                        <>
                          {getStatusBadge(downloadStatus.latent || 'idle')}
                          {downloadStatus.latent !== 'complete' && (
                            <Button 
                              size="sm" 
                              onClick={() => onDownload('latent')}
                              disabled={isDownloading}
                              className="flex items-center space-x-1"
                            >
                              <Download className="h-3 w-3" />
                              <span>Download</span>
                            </Button>
                          )}
                        </>
                      )}
                    </div>
                  </div>
                </CardContent>
              </Card>
            )}

            {/* Regenerate Option */}
            <div className="text-center pt-4 border-t">
              <Button 
                variant="outline" 
                onClick={onGenerate}
                disabled={isGeneratingEmbeddings || isDownloading}
                size="sm"
              >
                <RefreshCw className={`mr-2 h-4 w-4 ${isGeneratingEmbeddings ? 'animate-spin' : ''}`} />
                {isGeneratingEmbeddings ? 'Regenerating...' : 'Regenerate Embeddings'}
              </Button>
              <p className="text-xs text-gray-500 mt-1">
                This will overwrite existing embeddings
              </p>
            </div>
          </div>
        )}

        {/* File Format Info */}
        <div className="bg-gray-50 p-4 rounded">
          <h5 className="font-medium mb-2">📋 File Format Information</h5>
          <ul className="text-sm text-gray-600 space-y-1">
            <li>• <strong>CSV files</strong> with cell identifiers as row names</li>
            <li>• <strong>Interpretable embedding:</strong> 2 columns (UMAP coordinates)</li>
            <li>• <strong>Latent embedding:</strong> N columns (latent dimensions)</li>
            <li>• <strong>ZIP package:</strong> All results + training history</li>
            <li>• Compatible with <strong>Python, R, and Excel</strong></li>
          </ul>
        </div>
      </CardContent>
    </Card>
  );
};
