
// ========================================
// src/app/training/configure/page.tsx
// ========================================
'use client';

import React from 'react';
import { ArrowLeft } from 'lucide-react';
import Link from 'next/link';
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
        <Link 
          href="/"
          className="inline-flex items-center justify-center rounded-md border border-gray-300 bg-white px-4 py-2 text-sm font-medium text-gray-700 shadow-sm hover:bg-gray-50 focus:outline-none focus:ring-2 focus:ring-blue-500 focus:ring-offset-2 transition-colors"
        >
          <ArrowLeft className="mr-2 h-4 w-4" />
          Home
        </Link>
      </div>

      {/* Main Form Component */}
      <TrainingConfigForm />
    </div>
  );
}
