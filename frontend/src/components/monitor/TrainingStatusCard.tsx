// ========================================
// src/components/monitor/TrainingStatusCard.tsx
// ========================================
'use client';

import React, { useMemo } from 'react';
import { Square, TrendingUp } from 'lucide-react';
import Link from 'next/link';
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Progress } from '@/components/ui/progress';
import { Badge } from '@/components/ui/badge';
import { Alert, AlertDescription } from '@/components/ui/alert';
import type { TrainingState } from '@/types/index';

interface TrainingStatusCardProps {
  trainingState: TrainingState;
  onStop: () => void;
}

export const TrainingStatusCard: React.FC<TrainingStatusCardProps> = ({
  trainingState,
  onStop,
}) => {
  const progressPercentage = useMemo(() => {
    return Math.round((trainingState.current_epoch / trainingState.total_epochs) * 100);
  }, [trainingState.current_epoch, trainingState.total_epochs]);

  const formatDuration = (startTime?: string) => {
    if (!startTime) return 'Unknown';
    const start = new Date(startTime);
    const now = new Date();
    const diffMs = now.getTime() - start.getTime();
    const diffMins = Math.floor(diffMs / (1000 * 60));
    const diffSecs = Math.floor((diffMs % (1000 * 60)) / 1000);
    
    if (diffMins > 0) return `${diffMins}m ${diffSecs}s`;
    return `${diffSecs}s`;
  };

  return (
    <Card>
      <CardHeader className="pb-3">
        <div className="flex items-center justify-between">
          <CardTitle>Training Status</CardTitle>
          <Badge variant={trainingState.is_running ? "destructive" : "secondary"}>
            {trainingState.is_running ? "Running" : "Stopped"}
          </Badge>
        </div>
      </CardHeader>
      <CardContent className="space-y-4">
        <div>
          <div className="flex justify-between text-sm mb-2 font-medium">
            <span>Epoch {trainingState.current_epoch} / {trainingState.total_epochs}</span>
            <span>{progressPercentage}%</span>
          </div>
          <Progress value={progressPercentage} className="h-2" />
        </div>
        
        <div className="flex justify-between text-sm text-gray-600">
          <span>Duration</span>
          <span className="font-medium">{formatDuration(trainingState.started_at)}</span>
        </div>
        
        {/* Action Buttons */}
        <div className="pt-2 flex gap-2">
          {trainingState.is_running ? (
            <Button 
              variant="destructive" 
              onClick={onStop}
              className="w-full"
              size="sm"
            >
              <Square className="mr-2 h-4 w-4" />
              Stop Training
            </Button>
          ) : (
            <>
              <Button asChild className="flex-1" size="sm">
                <Link href="/results">
                  <TrendingUp className="mr-2 h-4 w-4" />
                  View Results
                </Link>
              </Button>
              <Button variant="outline" asChild className="flex-1" size="sm">
                <Link href="/training/configure">
                  New Training
                </Link>
              </Button>
            </>
          )}
        </div>

        {trainingState.error_message && (
          <Alert variant="destructive" className="mt-2">
            <AlertDescription className="text-xs">
              {trainingState.error_message}
            </AlertDescription>
          </Alert>
        )}
      </CardContent>
    </Card>
  );
};
