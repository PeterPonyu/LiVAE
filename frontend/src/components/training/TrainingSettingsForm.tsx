// ========================================
// src/components/training/TrainingSettingsForm.tsx
// ========================================
'use client';

import React from 'react';
import { UseFormReturn } from 'react-hook-form';
import { Info } from 'lucide-react';
import type { TrainingFormData } from '@/lib/validation/training-schemas';

interface TrainingSettingsFormProps {
  form: UseFormReturn<TrainingFormData>;
}

export const TrainingSettingsForm: React.FC<TrainingSettingsFormProps> = ({ form }) => {
  return (
    <div className="bg-white dark:bg-gray-800 rounded-lg border border-gray-200 dark:border-gray-700 shadow-sm">
      {/* Card Header */}
      <div className="px-6 py-5 border-b border-gray-200 dark:border-gray-700">
        <h3 className="text-lg font-semibold text-gray-900 dark:text-gray-100">
          Training Settings
        </h3>
        <p className="mt-1 text-sm text-gray-600 dark:text-gray-400">
          Control training duration and monitoring
        </p>
      </div>

      {/* Card Content */}
      <div className="px-6 py-6 space-y-4">
        {/* Epochs Field */}
        <FormFieldWrapper
          label="Epochs"
          description="Recommended: 1000-5000 epochs. More epochs = better convergence but longer training time."
          error={form.formState.errors.training_config?.epochs?.message}
        >
          <Input
            type="number"
            min={100}
            max={10000}
            value={form.watch('training_config.epochs')}
            onChange={(e) => form.setValue(
              'training_config.epochs',
              parseInt(e.target.value) || 0
            )}
          />
        </FormFieldWrapper>

        {/* Info Alert */}
        <div className="flex gap-3 rounded-lg border border-blue-200 dark:border-blue-800 bg-blue-50 dark:bg-blue-950 p-3">
          <Info className="h-5 w-5 text-blue-600 dark:text-blue-400 flex-shrink-0 mt-0.5" />
          <p className="text-sm text-blue-900 dark:text-blue-200">
            <strong className="font-semibold">Estimated time:</strong> ~5-15 minutes per 1000 epochs for datasets with 10,000 cells
          </p>
        </div>
      </div>
    </div>
  );
};

// ========================================
// Shared Helper Components
// ========================================

const FormFieldWrapper: React.FC<{
  label: string;
  description?: string;
  error?: string;
  children: React.ReactNode;
}> = ({ label, description, error, children }) => {
  return (
    <div className="space-y-2">
      <label className="text-sm font-medium text-gray-900 dark:text-gray-100">
        {label}
      </label>
      {children}
      {description && (
        <p className="text-xs text-gray-600 dark:text-gray-400">
          {description}
        </p>
      )}
      {error && (
        <p className="text-xs text-red-600 dark:text-red-400 flex items-center gap-1">
          <Info className="w-3 h-3" />
          {error}
        </p>
      )}
    </div>
  );
};

const Input: React.FC<React.InputHTMLAttributes<HTMLInputElement>> = (props) => {
  return (
    <input
      {...props}
      className="w-full px-3 py-2 text-sm bg-white dark:bg-gray-900 border border-gray-300 dark:border-gray-600 rounded-md shadow-sm focus:outline-none focus:ring-2 focus:ring-blue-500 focus:border-blue-500 dark:text-gray-100 placeholder-gray-400 dark:placeholder-gray-500 disabled:opacity-50 disabled:cursor-not-allowed"
    />
  );
};