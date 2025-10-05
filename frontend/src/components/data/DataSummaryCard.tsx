
// src/components/data/DataSummaryCard.tsx
'use client';

import React from 'react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import type { AnnDataSummary } from '@/types/index';
import { formatDistanceToNow } from 'date-fns';

interface DataSummaryCardProps {
  dataSummary: AnnDataSummary;
}

export const DataSummaryCard: React.FC<DataSummaryCardProps> = ({ dataSummary }) => {
  const formatFileSize = (sizeMB: number | undefined): string => {
    if (!sizeMB) return 'Unknown';
    return sizeMB < 1 ? `${(sizeMB * 1024).toFixed(1)} KB` : `${sizeMB.toFixed(1)} MB`;
  };

  const formatNumber = (num: number): string => {
    return num.toLocaleString();
  };

  const formatSparsity = (ratio: number | undefined): string => {
    if (ratio === undefined) return 'N/A';
    return `${(ratio * 100).toFixed(1)}%`;
  };

  return (
    <div className="space-y-6">
      {/* Basic File Info */}
      <Card>
        <CardHeader>
          <CardTitle className="flex items-center justify-between">
            Dataset Overview
            <Badge variant="secondary">
              {formatDistanceToNow(new Date(dataSummary.loaded_at), { addSuffix: true })}
            </Badge>
          </CardTitle>
          <CardDescription>
            {dataSummary.filename && `File: ${dataSummary.filename}`}
          </CardDescription>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
            <div className="text-center">
              <div className="text-2xl font-bold text-blue-600">
                {formatNumber(dataSummary.shape.n_obs)}
              </div>
              <div className="text-sm text-gray-600">Cells</div>
            </div>
            <div className="text-center">
              <div className="text-2xl font-bold text-green-600">
                {formatNumber(dataSummary.shape.n_vars)}
              </div>
              <div className="text-sm text-gray-600">Genes</div>
            </div>
            <div className="text-center">
              <div className="text-2xl font-bold text-purple-600">
                {formatFileSize(dataSummary.file_size_mb)}
              </div>
              <div className="text-sm text-gray-600">File Size</div>
            </div>
            <div className="text-center">
              <div className="text-2xl font-bold text-orange-600">
                {Object.keys(dataSummary.layers).length + 1}
              </div>
              <div className="text-sm text-gray-600">Layers</div>
            </div>
          </div>
        </CardContent>
      </Card>

      {/* Main Data Layer Info */}
      <Card>
        <CardHeader>
          <CardTitle>Main Data Layer (X)</CardTitle>
          <CardDescription>Primary expression matrix</CardDescription>
        </CardHeader>
        <CardContent>
          <div className="grid grid-cols-2 md:grid-cols-4 gap-4 text-sm">
            <div>
              <span className="font-medium">Data Type:</span>
              <div className="text-gray-600">{dataSummary.X_info.dtype}</div>
            </div>
            <div>
              <span className="font-medium">Sparsity:</span>
              <div className="text-gray-600">{formatSparsity(dataSummary.X_info.sparsity_ratio)}</div>
            </div>
            <div>
              <span className="font-medium">Min Value:</span>
              <div className="text-gray-600">
                {dataSummary.X_info.min_val?.toFixed(2) ?? 'N/A'}
              </div>
            </div>
            <div>
              <span className="font-medium">Max Value:</span>
              <div className="text-gray-600">
                {dataSummary.X_info.max_val?.toFixed(2) ?? 'N/A'}
              </div>
            </div>
          </div>
        </CardContent>
      </Card>

      {/* Additional Layers */}
      {Object.keys(dataSummary.layers).length > 0 && (
        <Card>
          <CardHeader>
            <CardTitle>Additional Layers</CardTitle>
            <CardDescription>Extra data matrices in the dataset</CardDescription>
          </CardHeader>
          <CardContent>
            <div className="space-y-3">
              {Object.entries(dataSummary.layers).map(([layerName, layerInfo]) => (
                <div key={layerName} className="flex items-center justify-between">
                  <div className="font-medium">{layerName}</div>
                  <div className="flex space-x-4 text-sm text-gray-600">
                    <span>{layerInfo.dtype}</span>
                    <span>{formatSparsity(layerInfo.sparsity_ratio)} sparse</span>
                  </div>
                </div>
              ))}
            </div>
          </CardContent>
        </Card>
      )}

      {/* Observations & Variables Summary */}
      <div className="grid md:grid-cols-2 gap-6">
        <Card>
          <CardHeader>
            <CardTitle>Cell Annotations ({dataSummary.obs_columns.length})</CardTitle>
            <CardDescription>Metadata columns for cells</CardDescription>
          </CardHeader>
          <CardContent className="space-y-2">
            {dataSummary.obs_columns.slice(0, 5).map((col) => (
              <div key={col.name} className="flex justify-between text-sm">
                <span className="font-medium">{col.name}</span>
                <span className="text-gray-600">
                  {col.n_unique} unique ({col.dtype})
                </span>
              </div>
            ))}
            {dataSummary.obs_columns.length > 5 && (
              <div className="text-sm text-gray-500">
                ... and {dataSummary.obs_columns.length - 5} more
              </div>
            )}
          </CardContent>
        </Card>

        <Card>
          <CardHeader>
            <CardTitle>Gene Annotations ({dataSummary.var_columns.length})</CardTitle>
            <CardDescription>Metadata columns for genes</CardDescription>
          </CardHeader>
          <CardContent className="space-y-2">
            {dataSummary.var_columns.slice(0, 5).map((col) => (
              <div key={col.name} className="flex justify-between text-sm">
                <span className="font-medium">{col.name}</span>
                <span className="text-gray-600">
                  {col.n_unique} unique ({col.dtype})
                </span>
              </div>
            ))}
            {dataSummary.var_columns.length > 5 && (
              <div className="text-sm text-gray-500">
                ... and {dataSummary.var_columns.length - 5} more
              </div>
            )}
          </CardContent>
        </Card>
      </div>

      {/* QC Metrics (if available) */}
      {dataSummary.qc_metrics && (
        <Card>
          <CardHeader>
            <CardTitle>Quality Control Metrics</CardTitle>
            <CardDescription>Basic QC statistics for the dataset</CardDescription>
          </CardHeader>
          <CardContent>
            <div className="grid grid-cols-2 md:grid-cols-4 gap-4 text-sm">
              <div className="text-center">
                <div className="font-medium">Genes per Cell</div>
                <div className="text-gray-600">
                  {dataSummary.qc_metrics.n_genes_by_counts.mean?.toFixed(0) ?? 'N/A'} ± {' '}
                  {dataSummary.qc_metrics.n_genes_by_counts.std?.toFixed(0) ?? 'N/A'}
                </div>
              </div>
              <div className="text-center">
                <div className="font-medium">Total Counts</div>
                <div className="text-gray-600">
                  {dataSummary.qc_metrics.total_counts.mean?.toFixed(0) ?? 'N/A'} ± {' '}
                  {dataSummary.qc_metrics.total_counts.std?.toFixed(0) ?? 'N/A'}
                </div>
              </div>
              <div className="text-center">
                <div className="font-medium">MT Gene %</div>
                <div className="text-gray-600">
                  {dataSummary.qc_metrics.pct_counts_mt.mean?.toFixed(1) ?? 'N/A'}% ± {' '}
                  {dataSummary.qc_metrics.pct_counts_mt.std?.toFixed(1) ?? 'N/A'}%
                </div>
              </div>
              <div className="text-center">
                <div className="font-medium">Ribosomal %</div>
                <div className="text-gray-600">
                  {dataSummary.qc_metrics.pct_counts_ribo.mean?.toFixed(1) ?? 'N/A'}% ± {' '}
                  {dataSummary.qc_metrics.pct_counts_ribo.std?.toFixed(1) ?? 'N/A'}%
                </div>
              </div>
            </div>
          </CardContent>
        </Card>
      )}
    </div>
  );
};
