// ========================================
// src/components/results/MetricsTab.tsx
// ========================================
import React from 'react';
import { BarChart3 } from 'lucide-react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from './ui';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, Legend } from 'recharts';
import type { TrainingMetrics } from '@/types';

interface MetricsTabProps {
  history: TrainingMetrics[];
}

export const MetricsTab: React.FC<MetricsTabProps> = ({ history }) => {
  const chartData = history.map((metrics) => ({
    epoch: metrics.epoch,
    loss: metrics.loss,
    ari: metrics.ari,
    nmi: metrics.nmi,
    asw: metrics.asw,
  }));

  return (
    <Card>
      <CardHeader>
        <CardTitle className="flex items-center gap-2">
          <BarChart3 className="h-5 w-5" />
          Training Progress
        </CardTitle>
        <CardDescription>
          {chartData.length} epochs tracked
        </CardDescription>
      </CardHeader>
      <CardContent>
        {chartData.length > 0 ? (
          <div className="space-y-6">
            {/* Main Chart */}
            <div className="h-80">
              <ResponsiveContainer width="100%" height="100%">
                <LineChart data={chartData}>
                  <CartesianGrid strokeDasharray="3 3" opacity={0.2} />
                  <XAxis 
                    dataKey="epoch" 
                    label={{ value: 'Epoch', position: 'insideBottom', offset: -5 }}
                  />
                  <YAxis yAxisId="left" orientation="left" domain={['auto', 'auto']} />
                  <YAxis yAxisId="right" orientation="right" domain={[0, 1]} />
                  <Tooltip 
                    formatter={(value: number) => value.toFixed(4)}
                    contentStyle={{ backgroundColor: 'white', border: '1px solid #e5e7eb', borderRadius: '6px' }}
                  />
                  <Legend />
                  <Line 
                    yAxisId="left" 
                    type="monotone" 
                    dataKey="loss" 
                    stroke="#dc2626" 
                    strokeWidth={2} 
                    name="Loss" 
                    dot={false}
                  />
                  <Line 
                    yAxisId="right" 
                    type="monotone" 
                    dataKey="ari" 
                    stroke="#16a34a" 
                    strokeWidth={2} 
                    name="ARI" 
                    dot={false}
                  />
                  <Line 
                    yAxisId="right" 
                    type="monotone" 
                    dataKey="nmi" 
                    stroke="#2563eb" 
                    strokeWidth={2} 
                    name="NMI" 
                    dot={false}
                  />
                  <Line 
                    yAxisId="right" 
                    type="monotone" 
                    dataKey="asw" 
                    stroke="#7c3aed" 
                    strokeWidth={2} 
                    name="ASW" 
                    dot={false}
                  />
                </LineChart>
              </ResponsiveContainer>
            </div>

            {/* Recent History Table */}
            <div>
              <h4 className="text-sm font-medium mb-3 text-gray-700">Last 10 Epochs</h4>
              <div className="overflow-x-auto rounded-lg border">
                <table className="w-full text-sm">
                  <thead className="bg-gray-50">
                    <tr>
                      <th className="p-3 text-left font-medium">Epoch</th>
                      <th className="p-3 text-left font-medium">Loss</th>
                      <th className="p-3 text-left font-medium">ARI</th>
                      <th className="p-3 text-left font-medium">NMI</th>
                      <th className="p-3 text-left font-medium">ASW</th>
                    </tr>
                  </thead>
                  <tbody className="divide-y">
                    {history.slice(-10).reverse().map((metrics) => (
                      <tr key={metrics.epoch} className="hover:bg-gray-50">
                        <td className="p-3 font-medium">{metrics.epoch}</td>
                        <td className="p-3 text-red-600">{metrics.loss.toFixed(4)}</td>
                        <td className="p-3 text-green-600">{metrics.ari.toFixed(3)}</td>
                        <td className="p-3 text-blue-600">{metrics.nmi.toFixed(3)}</td>
                        <td className="p-3 text-purple-600">{metrics.asw.toFixed(3)}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </div>
          </div>
        ) : (
          <div className="text-center py-8 text-gray-500">
            <BarChart3 className="h-8 w-8 mx-auto mb-2" />
            <p>No training metrics available</p>
          </div>
        )}
      </CardContent>
    </Card>
  );
};