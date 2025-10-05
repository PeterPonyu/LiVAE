// ========================================
// src/components/monitor/MetricsChart.tsx
// ========================================
'use client';

import React from 'react';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, Legend } from 'recharts';
import { METRICS_CONFIG, type MetricKey } from './constants';

interface MetricsChartProps {
  data: Record<string, unknown>[];
  metrics: MetricKey[];
  height?: number;
}

export const MetricsChart: React.FC<MetricsChartProps> = ({ 
  data, 
  metrics,
  height = 300 
}) => (
  <div style={{ height }}>
    <ResponsiveContainer width="100%" height="100%">
      <LineChart data={data} margin={{ top: 5, right: 30, left: 20, bottom: 5 }}>
        <CartesianGrid strokeDasharray="3 3" opacity={0.2} />
        <XAxis 
          dataKey="epoch" 
          tick={{ fontSize: 12 }}
          label={{ value: 'Epoch', position: 'insideBottom', offset: -5 }}
        />
        <YAxis tick={{ fontSize: 12 }} />
        <Tooltip 
          contentStyle={{ backgroundColor: 'white', border: '1px solid #e5e7eb', borderRadius: '6px' }}
          formatter={(value: number) => value.toFixed(4)}
        />
        <Legend />
        {metrics.map(metric => (
          <Line
            key={metric}
            type="monotone"
            dataKey={metric}
            name={METRICS_CONFIG[metric].label}
            stroke={METRICS_CONFIG[metric].color}
            strokeWidth={2}
            dot={false}
            activeDot={{ r: 4 }}
          />
        ))}
      </LineChart>
    </ResponsiveContainer>
  </div>
);