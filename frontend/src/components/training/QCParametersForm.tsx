
// src/components/training/QCParametersForm.tsx
'use client';

import React, { useState } from 'react';
import { UseFormReturn } from 'react-hook-form';
import { Info, ChevronDown } from 'lucide-react';
import type { TrainingFormData } from '@/lib/validation/training-schemas';

interface QCParametersFormProps {
  form: UseFormReturn<TrainingFormData>;
}

export const QCParametersForm: React.FC<QCParametersFormProps> = ({ form }) => {
  const applyQC = form.watch('training_config.apply_qc');

  const handleToggleQC = () => {
    const currentValue = form.getValues('training_config.apply_qc');
    console.log('Current QC value:', currentValue);
    console.log('Setting to:', !currentValue);
    
    form.setValue('training_config.apply_qc', !currentValue, {
      shouldValidate: true,
      shouldDirty: true,
      shouldTouch: true
    });
    
    // Force re-render by triggering form state change
    form.trigger('training_config.apply_qc');
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
        <div className="flex items-center justify-between">
          <div className="space-y-0.5">
            <label className="text-base font-medium text-gray-900 dark:text-gray-100">
              Apply Quality Control
            </label>
            <p className="text-sm text-gray-600 dark:text-gray-400">
              Enable preprocessing filters to improve data quality
            </p>
          </div>
          <Switch
            checked={applyQC}
            onToggle={handleToggleQC}
          />
        </div>

        {/* Debug Info - Remove this after testing */}
        <div className="text-xs text-gray-500">
          Current state: {applyQC ? 'ON' : 'OFF'}
        </div>

        {/* QC Parameters Form */}
        {applyQC && (
          <div className="space-y-6 pt-4 border-t border-gray-200 dark:border-gray-700">
            {/* Species Selection */}
            <FormFieldWrapper
              label="Species"
              description="Species type for mitochondrial gene detection"
              error={form.formState.errors.training_config?.qc_params?.species?.message}
            >
              <Select
                value={form.watch('training_config.qc_params.species')}
                onChange={(value) => form.setValue('training_config.qc_params.species', value as 'human' | 'mouse')}
                options={[
                  { value: 'human', label: 'Human' },
                  { value: 'mouse', label: 'Mouse' }
                ]}
                placeholder="Select species"
              />
            </FormFieldWrapper>

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

const Select: React.FC<{
  value?: string;
  onChange: (value: string) => void;
  options: { value: string; label: string }[];
  placeholder?: string;
}> = ({ value, onChange, options, placeholder }) => {
  const [isOpen, setIsOpen] = useState(false);
  const selectedOption = options.find(opt => opt.value === value);

  return (
    <div className="relative">
      <button
        type="button"
        onClick={() => setIsOpen(!isOpen)}
        className="w-full px-3 py-2 text-sm bg-white dark:bg-gray-900 border border-gray-300 dark:border-gray-600 rounded-md shadow-sm focus:outline-none focus:ring-2 focus:ring-blue-500 focus:border-blue-500 dark:text-gray-100 flex items-center justify-between"
      >
        <span className={selectedOption ? 'text-gray-900 dark:text-gray-100' : 'text-gray-400 dark:text-gray-500'}>
          {selectedOption ? selectedOption.label : placeholder}
        </span>
        <ChevronDown className={`w-4 h-4 text-gray-500 transition-transform ${isOpen ? 'rotate-180' : ''}`} />
      </button>

      {isOpen && (
        <>
          <div className="fixed inset-0 z-10" onClick={() => setIsOpen(false)} />
          <div className="absolute z-20 w-full mt-1 bg-white dark:bg-gray-900 border border-gray-300 dark:border-gray-600 rounded-md shadow-lg">
            {options.map((option) => (
              <button
                key={option.value}
                type="button"
                onClick={() => {
                  onChange(option.value);
                  setIsOpen(false);
                }}
                className={`w-full px-3 py-2 text-sm text-left hover:bg-gray-100 dark:hover:bg-gray-800 ${
                  value === option.value
                    ? 'bg-blue-50 dark:bg-blue-950 text-blue-600 dark:text-blue-400'
                    : 'text-gray-900 dark:text-gray-100'
                }`}
              >
                {option.label}
              </button>
            ))}
          </div>
        </>
      )}
    </div>
  );
};

// ========================================
// FIXED Switch Component
// ========================================

const Switch: React.FC<{
  checked: boolean;
  onToggle: () => void;
}> = ({ checked, onToggle }) => {
  const handleClick = (e: React.MouseEvent<HTMLButtonElement>) => {
    e.preventDefault();
    e.stopPropagation();
    console.log('Switch clicked, current state:', checked);
    onToggle();
  };

  return (
    <button
      type="button"
      onClick={handleClick}
      className={`relative inline-flex h-6 w-11 flex-shrink-0 cursor-pointer items-center rounded-full transition-colors duration-200 ease-in-out focus:outline-none focus:ring-2 focus:ring-blue-500 focus:ring-offset-2 ${
        checked ? 'bg-blue-600' : 'bg-gray-200 dark:bg-gray-700'
      }`}
    >
      <span className="sr-only">Toggle quality control</span>
      <span
        aria-hidden="true"
        className={`pointer-events-none inline-block h-4 w-4 transform rounded-full bg-white shadow-lg ring-0 transition duration-200 ease-in-out ${
          checked ? 'translate-x-6' : 'translate-x-1'
        }`}
      />
    </button>
  );
};
