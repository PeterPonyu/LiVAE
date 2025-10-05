// ========================================
// src/app/training/configure/page.tsx
// ========================================
'use client';

import React from 'react';
import { ArrowLeft } from 'lucide-react';
import Link from 'next/link';
import { Button } from '@/components/ui/button';
import { TrainingConfigForm } from '@/components/training/TrainingConfigForm';

export default function TrainingConfigPage() {
  return (
    <div className="container mx-auto py-8 px-4 max-w-6xl">
      {/* Header */}
      <div className="flex items-center justify-between mb-6">
        <div>
          <h1 className="text-3xl font-bold">Training Configuration</h1>
          <p className="text-gray-600 mt-1">
            Configure model parameters and start training
          </p>
        </div>
        <Button variant="outline" asChild>
          <Link href="/">
            <ArrowLeft className="mr-2 h-4 w-4" />
            Home
          </Link>
        </Button>
      </div>

      {/* Main Form Component */}
      <TrainingConfigForm />
    </div>
  );
}