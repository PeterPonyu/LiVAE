// ========================================
// src/components/training/DatasetSummaryCard.tsx
// ========================================
'use client';

import React from 'react';
import type { AnnDataSummary } from '@/types/index';

interface DatasetSummaryCardProps {
  dataSummary: AnnDataSummary;
}

export const DatasetSummaryCard: React.FC<DatasetSummaryCardProps> = ({ dataSummary }) => {
  const availableLayers = ['X', ...Object.keys(dataSummary.layers)];

  return (
    <div className="bg-white dark:bg-gray-800 rounded-lg border border-gray-200 dark:border-gray-700 shadow-sm mb-6">
      {/* Card Header */}
      <div className="px-6 py-4 border-b border-gray-200 dark:border-gray-700">
        <div className="flex items-center justify-between">
          <h3 className="text-lg font-semibold text-gray-900 dark:text-gray-100">
            Dataset Overview
          </h3>
          <span className="inline-flex items-center px-2.5 py-0.5 rounded-full text-xs font-medium bg-gray-100 dark:bg-gray-800 text-gray-800 dark:text-gray-200 border border-gray-300 dark:border-gray-600">
            {dataSummary.shape.n_obs.toLocaleString()} cells
          </span>
        </div>
      </div>

      {/* Card Content */}
      <div className="px-6 py-5">
        <div className="grid grid-cols-2 sm:grid-cols-4 gap-4 text-sm">
          <div className="flex flex-col space-y-1">
            <span className="text-gray-600 dark:text-gray-400">Cells</span>
            <span className="font-semibold text-gray-900 dark:text-gray-100">
              {dataSummary.shape.n_obs.toLocaleString()}
            </span>
          </div>
          <div className="flex flex-col space-y-1">
            <span className="text-gray-600 dark:text-gray-400">Genes</span>
            <span className="font-semibold text-gray-900 dark:text-gray-100">
              {dataSummary.shape.n_vars.toLocaleString()}
            </span>
          </div>
          <div className="flex flex-col space-y-1">
            <span className="text-gray-600 dark:text-gray-400">Layers</span>
            <span className="font-semibold text-gray-900 dark:text-gray-100">
              {availableLayers.length}
            </span>
          </div>
          <div className="flex flex-col space-y-1">
            <span className="text-gray-600 dark:text-gray-400">File</span>
            <span 
              className="font-semibold text-gray-900 dark:text-gray-100 truncate" 
              title={dataSummary.filename}
            >
              {dataSummary.filename || 'Unknown'}
            </span>
          </div>
        </div>
      </div>
    </div>
  );
};
