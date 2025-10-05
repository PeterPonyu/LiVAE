
// src/components/training/QCParametersForm.tsx
'use client';

import React from 'react';
import { UseFormReturn } from 'react-hook-form';
import { Info } from 'lucide-react';
import type { TrainingFormData } from '@/lib/validation/training-schemas';

interface QCParametersFormProps {
  form: UseFormReturn<TrainingFormData>;
}

export const QCParametersForm: React.FC<QCParametersFormProps> = ({ form }) => {
  const applyQC = form.watch('training_config.apply_qc');

  const handleToggleQC = () => {
    const currentValue = form.getValues('training_config.apply_qc');
    form.setValue('training_config.apply_qc', !currentValue, {
      shouldValidate: true,
      shouldDirty: true,
      shouldTouch: true
    });
  };

  return (
    <div className="bg-white dark:bg-gray-800 rounded-lg border border-gray-200 dark:border-gray-700 shadow-sm">
      {/* Card Header */}
      <div className="px-6 py-5 border-b border-gray-200 dark:border-gray-700">
        <h3 className="text-lg font-semibold text-gray-900 dark:text-gray-100">
          Quality Control Filters
        </h3>
        <p className="mt-1 text-sm text-gray-600 dark:text-gray-400">
          Apply preprocessing filters to remove low-quality cells and genes
        </p>
      </div>

      {/* Card Content */}
      <div className="px-6 py-6 space-y-6">
        {/* Enable/Disable QC Toggle */}
        <div className={`rounded-lg border-2 p-4 transition-all ${
          applyQC 
            ? 'border-blue-500 bg-blue-50 dark:bg-blue-950/20' 
            : 'border-gray-200 dark:border-gray-700 bg-gray-50 dark:bg-gray-800/50'
        }`}>
          <div className="flex items-center justify-between gap-4">
            <div className="flex-1">
              <div className="flex items-center gap-2">
                <label className="text-base font-semibold text-gray-900 dark:text-gray-100 cursor-pointer">
                  Apply Quality Control
                </label>
                <span className={`inline-flex items-center px-2 py-0.5 rounded-full text-xs font-medium ${
                  applyQC 
                    ? 'bg-blue-100 dark:bg-blue-900 text-blue-800 dark:text-blue-200' 
                    : 'bg-gray-200 dark:bg-gray-700 text-gray-600 dark:text-gray-400'
                }`}>
                  {applyQC ? 'Enabled' : 'Disabled'}
                </span>
              </div>
              <p className="mt-1 text-sm text-gray-600 dark:text-gray-400">
                {applyQC 
                  ? 'QC filters will be applied before training' 
                  : 'Enable to configure quality control parameters'
                }
              </p>
            </div>
            <Switch checked={applyQC} onToggle={handleToggleQC} />
          </div>
        </div>

        {/* QC Parameters Form */}
        {applyQC && (
          <div className="space-y-6 pt-4 border-t border-gray-200 dark:border-gray-700">
            {/* Species Selection */}
            <SelectField
              label="Species"
              value={form.watch('training_config.qc_params.species') || 'human'}
              options={[
                { value: 'human', label: 'Human' },
                { value: 'mouse', label: 'Mouse' }
              ]}
              description="Species type for mitochondrial gene detection"
              onChange={(value) => form.setValue('training_config.qc_params.species', value as 'human' | 'mouse')}
              error={form.formState.errors.training_config?.qc_params?.species?.message}
            />

            {/* Cell Filters */}
            <div className="space-y-4">
              <h4 className="text-sm font-semibold text-gray-900 dark:text-gray-100">
                Cell Quality Filters
              </h4>

              <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                <FormFieldWrapper
                  label="Minimum Genes per Cell"
                  description="Remove cells with fewer detected genes"
                  error={form.formState.errors.training_config?.qc_params?.min_genes_per_cell?.message}
                >
                  <Input
                    type="number"
                    placeholder="e.g., 200"
                    value={form.watch('training_config.qc_params.min_genes_per_cell') ?? ''}
                    onChange={(e) => form.setValue(
                      'training_config.qc_params.min_genes_per_cell',
                      e.target.value ? parseInt(e.target.value) : undefined
                    )}
                  />
                </FormFieldWrapper>

                <FormFieldWrapper
                  label="Maximum Genes per Cell"
                  description="Remove cells with too many detected genes"
                  error={form.formState.errors.training_config?.qc_params?.max_genes_per_cell?.message}
                >
                  <Input
                    type="number"
                    placeholder="e.g., 5000"
                    value={form.watch('training_config.qc_params.max_genes_per_cell') ?? ''}
                    onChange={(e) => form.setValue(
                      'training_config.qc_params.max_genes_per_cell',
                      e.target.value ? parseInt(e.target.value) : undefined
                    )}
                  />
                </FormFieldWrapper>

                <FormFieldWrapper
                  label="Minimum UMI per Cell"
                  description="Remove cells with low total UMI counts"
                  error={form.formState.errors.training_config?.qc_params?.min_counts_per_cell?.message}
                >
                  <Input
                    type="number"
                    placeholder="e.g., 500"
                    value={form.watch('training_config.qc_params.min_counts_per_cell') ?? ''}
                    onChange={(e) => form.setValue(
                      'training_config.qc_params.min_counts_per_cell',
                      e.target.value ? parseInt(e.target.value) : undefined
                    )}
                  />
                </FormFieldWrapper>

                <FormFieldWrapper
                  label="Maximum UMI per Cell"
                  description="Remove cells with extremely high UMI counts"
                  error={form.formState.errors.training_config?.qc_params?.max_counts_per_cell?.message}
                >
                  <Input
                    type="number"
                    placeholder="e.g., 30000"
                    value={form.watch('training_config.qc_params.max_counts_per_cell') ?? ''}
                    onChange={(e) => form.setValue(
                      'training_config.qc_params.max_counts_per_cell',
                      e.target.value ? parseInt(e.target.value) : undefined
                    )}
                  />
                </FormFieldWrapper>
              </div>
            </div>

            {/* Gene Type Filters */}
            <div className="space-y-4">
              <h4 className="text-sm font-semibold text-gray-900 dark:text-gray-100">
                Gene Type Filters
              </h4>

              <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                <FormFieldWrapper
                  label="Maximum Mitochondrial %"
                  description="Remove cells with high mitochondrial gene expression"
                  error={form.formState.errors.training_config?.qc_params?.max_pct_mt?.message}
                >
                  <Input
                    type="number"
                    placeholder="e.g., 20"
                    step="0.1"
                    value={form.watch('training_config.qc_params.max_pct_mt') ?? ''}
                    onChange={(e) => form.setValue(
                      'training_config.qc_params.max_pct_mt',
                      e.target.value ? parseFloat(e.target.value) : undefined
                    )}
                  />
                </FormFieldWrapper>

                <FormFieldWrapper
                  label="Maximum Ribosomal %"
                  description="Remove cells with excessive ribosomal gene expression"
                  error={form.formState.errors.training_config?.qc_params?.max_pct_ribo?.message}
                >
                  <Input
                    type="number"
                    placeholder="e.g., 50"
                    step="0.1"
                    value={form.watch('training_config.qc_params.max_pct_ribo') ?? ''}
                    onChange={(e) => form.setValue(
                      'training_config.qc_params.max_pct_ribo',
                      e.target.value ? parseFloat(e.target.value) : undefined
                    )}
                  />
                </FormFieldWrapper>
              </div>
            </div>
          </div>
        )}
      </div>
    </div>
  );
};

// ========================================
// Custom Components
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
      className="w-full px-3 py-2 text-sm bg-white dark:bg-gray-900 border border-gray-300 dark:border-gray-600 rounded-md shadow-sm focus:outline-none focus:ring-2 focus:ring-blue-500 focus:border-blue-500 dark:text-gray-100 placeholder-gray-400 dark:placeholder-gray-500"
    />
  );
};

/**
 * Select field component
 */
const SelectField: React.FC<{
  label: string;
  value: string;
  options: { value: string; label: string }[];
  description: string;
  onChange: (value: string) => void;
  error?: string;
}> = ({ label, value, options, description, onChange, error }) => {
  return (
    <div className="space-y-2">
      <label className="text-sm font-medium text-gray-900 dark:text-gray-100">
        {label}
      </label>
      <select
        value={value}
        onChange={(e) => onChange(e.target.value)}
        title="Select species"
        className="w-full px-3 py-2 text-sm bg-white dark:bg-gray-900 border border-gray-300 dark:border-gray-600 rounded-md shadow-sm focus:outline-none focus:ring-2 focus:ring-blue-500 focus:border-blue-500 dark:text-gray-100 cursor-pointer"
      >
        {options.map((option) => (
          <option key={option.value} value={option.value}>
            {option.label}
          </option>
        ))}
      </select>
      <p className="text-xs text-gray-600 dark:text-gray-400">
        {description}
      </p>
      {error && (
        <p className="text-xs text-red-600 dark:text-red-400 flex items-center gap-1">
          <Info className="w-3 h-3" />
          {error}
        </p>
      )}
    </div>
  );
};

const Switch: React.FC<{
  checked: boolean;
  onToggle: () => void;
}> = ({ checked, onToggle }) => {
  return (
    <button
      type="button"
      onClick={(e) => {
        e.preventDefault();
        e.stopPropagation();
        onToggle();
      }}
      className={`relative inline-flex h-8 w-14 flex-shrink-0 cursor-pointer items-center rounded-full transition-all duration-200 ease-in-out focus:outline-none focus:ring-2 focus:ring-blue-500 focus:ring-offset-2 shadow-sm ${
        checked 
          ? 'bg-blue-600 hover:bg-blue-700' 
          : 'bg-gray-300 dark:bg-gray-600 hover:bg-gray-400 dark:hover:bg-gray-500'
      }`}
    >
      {/* Thumb */}
      <span
        className={`inline-block h-6 w-6 transform rounded-full bg-white shadow-lg ring-0 transition-transform duration-200 ease-in-out ${
          checked ? 'translate-x-7' : 'translate-x-1'
        }`}
      />
    </button>
  );
};
