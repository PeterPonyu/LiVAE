// ========================================
// src/components/results/AnalysisTab.tsx
// ========================================
import React, { Children } from 'react';
import { CheckCircle, AlertTriangle } from 'lucide-react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle, Badge, Alert, AlertDescription } from './ui';
import type { TrainingMetrics } from '@/types';

interface AnalysisTabProps {
  latestMetrics: TrainingMetrics | undefined;
  chartData: Array<{ epoch: number; loss: number; ari: number; nmi: number; asw: number; }>;
}

const getQualityBadge = (value: number, thresholds: [number, number]) => {
  const [good, excellent] = thresholds;
  if (value >= excellent) return { variant: 'default' as const, children: 'Excellent' };
  if (value >= good) return { variant: 'secondary' as const, children: 'Good' };
  return { variant: 'destructive' as const, children: 'Poor' };
};

export const AnalysisTab: React.FC<AnalysisTabProps> = ({ latestMetrics, chartData }) => {
  return (
    <div className="space-y-4">
      {/* Performance Overview */}
      {latestMetrics && (
        <Card>
          <CardHeader>
            <CardTitle>Performance Assessment</CardTitle>
            <CardDescription>Quality metrics from final epoch</CardDescription>
          </CardHeader>
          <CardContent>
            <div className="grid md:grid-cols-2 gap-4">
              {/* Clustering Quality */}
              <div className="space-y-3">
                <h4 className="font-medium text-sm">Clustering Quality</h4>
                <div className="space-y-2">
                  <div className="flex justify-between items-center p-2 bg-gray-50 rounded">
                    <span className="text-sm">ARI: {latestMetrics.ari.toFixed(3)}</span>
                    <Badge {...getQualityBadge(latestMetrics.ari, [0.4, 0.7])} />
                  </div>
                  <div className="flex justify-between items-center p-2 bg-gray-50 rounded">
                    <span className="text-sm">NMI: {latestMetrics.nmi.toFixed(3)}</span>
                    <Badge {...getQualityBadge(latestMetrics.nmi, [0.4, 0.7])} />
                  </div>
                </div>
              </div>

              {/* Separation Quality */}
              <div className="space-y-3">
                <h4 className="font-medium text-sm">Embedding Quality</h4>
                <div className="space-y-2">
                  <div className="flex justify-between items-center p-2 bg-gray-50 rounded">
                    <span className="text-sm">ASW: {latestMetrics.asw.toFixed(3)}</span>
                    <Badge {...getQualityBadge(latestMetrics.asw, [0.2, 0.5])} />
                  </div>
                  <div className="flex justify-between items-center p-2 bg-gray-50 rounded">
                    <span className="text-sm">Final Loss: {latestMetrics.loss.toFixed(4)}</span>
                    <Badge variant="secondary">
                      {chartData.length > 5 && 
                       chartData[chartData.length - 1].loss < chartData[chartData.length - 6].loss * 0.9
                        ? "Converged" : "Stable"}
                    </Badge>
                  </div>
                </div>
              </div>
            </div>
          </CardContent>
        </Card>
      )}

      {/* Recommendations */}
      <Card>
        <CardHeader>
          <CardTitle>Recommendations</CardTitle>
        </CardHeader>
        <CardContent className="space-y-3">
          {latestMetrics && (
            <>
              {/* Success Message */}
              {latestMetrics.ari > 0.7 && latestMetrics.nmi > 0.7 && (
                <Alert className="bg-green-50 border-green-200">
                  <div className="flex items-start gap-3">
                    <CheckCircle className="h-5 w-5 text-green-600 mt-0.5" />
                    <AlertDescription className="text-green-800">
                      <strong>Excellent results!</strong> Your model has learned high-quality cell embeddings suitable for downstream analysis.
                    </AlertDescription>
                  </div>
                </Alert>
              )}
              
              {/* Warning Messages */}
              {latestMetrics.ari < 0.4 && (
                <Alert className="bg-yellow-50 border-yellow-200">
                  <div className="flex items-start gap-3">
                    <AlertTriangle className="h-5 w-5 text-yellow-600 mt-0.5" />
                    <AlertDescription className="text-yellow-800">
                      <strong>Low clustering accuracy.</strong> Consider training longer or adjusting hyperparameters (beta, tc, latent_dim).
                    </AlertDescription>
                  </div>
                </Alert>
              )}
              
              {latestMetrics.asw < 0.2 && latestMetrics.ari >= 0.4 && (
                <Alert className="bg-yellow-50 border-yellow-200">
                  <div className="flex items-start gap-3">
                    <AlertTriangle className="h-5 w-5 text-yellow-600 mt-0.5" />
                    <AlertDescription className="text-yellow-800">
                      <strong>Clusters need better separation.</strong> Try increasing beta or adjusting the latent dimension.
                    </AlertDescription>
                  </div>
                </Alert>
              )}
            </>
          )}
          
          {/* Next Steps */}
          <div className="bg-blue-50 p-4 rounded-lg border border-blue-200">
            <h5 className="font-medium mb-2 text-blue-900">Next Steps</h5>
            <ul className="space-y-1 text-sm text-blue-800">
              <li>• Download embeddings for visualization and analysis</li>
              <li>• Use latent embeddings for clustering and trajectory inference</li>
              <li>• Validate clusters with marker gene expression</li>
              <li>• Export results for collaborative research</li>
            </ul>
          </div>
        </CardContent>
      </Card>
    </div>
  );
};