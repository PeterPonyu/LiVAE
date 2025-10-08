
// ========================================
// src/components/monitor/MetricsChartsSection.tsx
// ========================================
'use client';

import React, { useState } from 'react';
import { MetricsChart } from './MetricsChart';

interface MetricsChartsSectionProps {
  chartData: Record<string, unknown>[];
}

type TabValue = 'primary' | 'clustering' | 'all';

export const MetricsChartsSection: React.FC<MetricsChartsSectionProps> = ({ chartData }) => {
  const [activeTab, setActiveTab] = useState<TabValue>('primary');

  return (
    <div className="rounded-lg border border-gray-200 bg-white shadow-sm">
      {/* Card Header */}
      <div className="flex flex-col space-y-1.5 p-6">
        <h3 className="text-2xl font-semibold leading-none tracking-tight">
          Metrics Visualization
        </h3>
        <p className="text-sm text-gray-500">
          Training progress across all evaluation metrics
        </p>
      </div>

      {/* Card Content */}
      <div className="p-6 pt-0">
        <div className="space-y-4">
          {/* Tabs List */}
          <div className="inline-flex h-10 items-center justify-center rounded-md bg-gray-100 p-1 text-gray-500 w-full">
            <button
              onClick={() => setActiveTab('primary')}
              className={`inline-flex items-center justify-center whitespace-nowrap rounded-sm px-3 py-1.5 text-sm font-medium ring-offset-white transition-all focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-blue-500 focus-visible:ring-offset-2 disabled:pointer-events-none disabled:opacity-50 flex-1 ${
                activeTab === 'primary'
                  ? 'bg-white text-gray-900 shadow-sm'
                  : 'hover:bg-gray-200/50'
              }`}
            >
              Primary Metrics
            </button>
            <button
              onClick={() => setActiveTab('clustering')}
              className={`inline-flex items-center justify-center whitespace-nowrap rounded-sm px-3 py-1.5 text-sm font-medium ring-offset-white transition-all focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-blue-500 focus-visible:ring-offset-2 disabled:pointer-events-none disabled:opacity-50 flex-1 ${
                activeTab === 'clustering'
                  ? 'bg-white text-gray-900 shadow-sm'
                  : 'hover:bg-gray-200/50'
              }`}
            >
              Clustering Quality
            </button>
            <button
              onClick={() => setActiveTab('all')}
              className={`inline-flex items-center justify-center whitespace-nowrap rounded-sm px-3 py-1.5 text-sm font-medium ring-offset-white transition-all focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-blue-500 focus-visible:ring-offset-2 disabled:pointer-events-none disabled:opacity-50 flex-1 ${
                activeTab === 'all'
                  ? 'bg-white text-gray-900 shadow-sm'
                  : 'hover:bg-gray-200/50'
              }`}
            >
              All Metrics
            </button>
          </div>

          {/* Primary Metrics Tab */}
          {activeTab === 'primary' && (
            <div className="space-y-4">
              <div className="grid md:grid-cols-2 gap-6">
                <div>
                  <h3 className="text-sm font-medium mb-2 text-gray-700">Training Loss</h3>
                  <MetricsChart data={chartData} metrics={['loss']} height={250} />
                </div>
                <div>
                  <h3 className="text-sm font-medium mb-2 text-gray-700">Clustering Accuracy (ARI)</h3>
                  <MetricsChart data={chartData} metrics={['ari']} height={250} />
                </div>
              </div>
            </div>
          )}

          {/* Clustering Quality Tab */}
          {activeTab === 'clustering' && (
            <div className="space-y-4">
              <div>
                <h3 className="text-sm font-medium mb-2 text-gray-700">Quality Metrics Comparison</h3>
                <MetricsChart 
                  data={chartData} 
                  metrics={['ari', 'nmi', 'asw', 'pc']} 
                  height={350} 
                />
              </div>
              <div className="grid md:grid-cols-3 gap-4 text-xs">
                <div className="bg-green-50 p-3 rounded">
                  <strong>ARI:</strong> Clustering accuracy vs ground truth
                </div>
                <div className="bg-blue-50 p-3 rounded">
                  <strong>NMI:</strong> Mutual information measure
                </div>
                <div className="bg-purple-50 p-3 rounded">
                  <strong>ASW:</strong> Cluster separation quality
                </div>
              </div>
            </div>
          )}

          {/* All Metrics Tab */}
          {activeTab === 'all' && (
            <div className="space-y-4">
              <div className="grid md:grid-cols-2 gap-6">
                <div>
                  <h3 className="text-sm font-medium mb-2 text-gray-700">Loss & Accuracy</h3>
                  <MetricsChart data={chartData} metrics={['loss', 'ari']} height={250} />
                </div>
                <div>
                  <h3 className="text-sm font-medium mb-2 text-gray-700">Separation Metrics</h3>
                  <MetricsChart data={chartData} metrics={['asw', 'ch', 'db']} height={250} />
                </div>
              </div>
              <div>
                <h3 className="text-sm font-medium mb-2 text-gray-700">Information & Connectivity</h3>
                <MetricsChart data={chartData} metrics={['nmi', 'pc']} height={250} />
              </div>
            </div>
          )}
        </div>
      </div>
    </div>
  );
};
