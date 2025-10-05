// ========================================
// src/components/monitor/MetricsChartsSection.tsx
// ========================================
'use client';

import React from 'react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import { MetricsChart } from './MetricsChart';

interface MetricsChartsSectionProps {
  chartData: Record<string, unknown>[];
}

export const MetricsChartsSection: React.FC<MetricsChartsSectionProps> = ({ chartData }) => {
  return (
    <Card>
      <CardHeader>
        <CardTitle>Metrics Visualization</CardTitle>
        <CardDescription>
          Training progress across all evaluation metrics
        </CardDescription>
      </CardHeader>
      <CardContent>
        <Tabs defaultValue="primary" className="space-y-4">
          <TabsList className="grid w-full grid-cols-3">
            <TabsTrigger value="primary">Primary Metrics</TabsTrigger>
            <TabsTrigger value="clustering">Clustering Quality</TabsTrigger>
            <TabsTrigger value="all">All Metrics</TabsTrigger>
          </TabsList>

          <TabsContent value="primary" className="space-y-4">
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
          </TabsContent>

          <TabsContent value="clustering" className="space-y-4">
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
          </TabsContent>

          <TabsContent value="all" className="space-y-4">
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
          </TabsContent>
        </Tabs>
      </CardContent>
    </Card>
  );
};
