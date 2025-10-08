// ========================================
// src/components/results/ResultsHeader.tsx
// ========================================
import React from 'react';
import Link from 'next/link';
import { ArrowLeft, RefreshCw, TrendingUp } from 'lucide-react';
import { Button, Badge } from './ui';

interface ResultsHeaderProps {
  trainingCompleted: boolean;
  cellCount: number;
  isLoading: boolean;
  onRefresh: () => void;
}

export const ResultsHeader: React.FC<ResultsHeaderProps> = ({
  trainingCompleted,
  cellCount,
  isLoading,
  onRefresh
}) => {
  return (
    <div className="flex items-center justify-between mb-6">
      <div>
        <div className="flex items-center gap-3">
          <h1 className="text-3xl font-bold">Training Results</h1>
          <Badge variant={trainingCompleted ? "default" : "secondary"}>
            {trainingCompleted ? "Complete" : "In Progress"}
          </Badge>
        </div>
        <p className="text-gray-600 mt-1">
          {trainingCompleted 
            ? `Successfully trained on ${cellCount.toLocaleString()} cells`
            : 'Training incomplete - some features may be unavailable'}
        </p>
      </div>
      <div className="flex items-center gap-2">
        <Button variant="outline" onClick={onRefresh} size="sm" disabled={isLoading}>
          <RefreshCw className={`h-4 w-4 ${isLoading ? 'animate-spin' : ''}`} />
        </Button>
        {!trainingCompleted && (
          <Link href="/training/monitor">
            <Button variant="outline" size="sm">
              <TrendingUp className="mr-2 h-4 w-4" />
              Monitor
            </Button>
          </Link>
        )}
        <Link href="/">
          <Button variant="outline">
            <ArrowLeft className="mr-2 h-4 w-4" />
            Home
          </Button>
        </Link>
      </div>
    </div>
  );
};