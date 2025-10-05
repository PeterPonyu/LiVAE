// ========================================
// src/components/monitor/CurrentMetricsCard.tsx
// ========================================
'use client';

import React from 'react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import type { TrainingMetrics } from '@/types/index';

interface CurrentMetricsCardProps {
  latestMetrics?: TrainingMetrics | null;
}

export const CurrentMetricsCard: React.FC<CurrentMetricsCardProps> = ({ latestMetrics }) => {
  return (
    <Card>
      <CardHeader className="pb-3">
        <CardTitle>Current Metrics</CardTitle>
        <CardDescription>
          {latestMetrics && `Epoch ${latestMetrics.epoch}`}
        </CardDescription>
      </CardHeader>
      <CardContent>
        {latestMetrics ? (
          <div className="grid grid-cols-4 gap-3">
            <MetricDisplay label="Loss" value={latestMetrics.loss} color="text-red-600" decimals={3} />
            <MetricDisplay label="ARI" value={latestMetrics.ari} color="text-green-600" decimals={3} />
            <MetricDisplay label="NMI" value={latestMetrics.nmi} color="text-blue-600" decimals={3} />
            <MetricDisplay label="ASW" value={latestMetrics.asw} color="text-purple-600" decimals={3} />
            <MetricDisplay label="CH" value={latestMetrics.ch} color="text-orange-600" decimals={1} size="sm" />
            <MetricDisplay label="DB" value={latestMetrics.db} color="text-red-500" decimals={3} size="sm" />
            <MetricDisplay label="PC" value={latestMetrics.pc} color="text-green-500" decimals={3} size="sm" />
            <div className="text-center flex items-center justify-center">
              <Badge variant="outline" className="text-xs">
                Live
              </Badge>
            </div>
          </div>
        ) : (
          <p className="text-gray-500 text-sm">No metrics available yet</p>
        )}
      </CardContent>
    </Card>
  );
};

const MetricDisplay: React.FC<{
  label: string;
  value: number;
  color: string;
  decimals: number;
  size?: 'sm' | 'lg';
}> = ({ label, value, color, decimals, size = 'lg' }) => (
  <div className="text-center">
    <div className="text-xs text-gray-600 mb-1">{label}</div>
    <div className={`${size === 'lg' ? 'text-lg' : 'text-sm'} font-bold ${color}`}>
      {value.toFixed(decimals)}
    </div>
  </div>
);
