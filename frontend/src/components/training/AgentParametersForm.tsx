
// src/components/training/AgentParametersForm.tsx
'use client';

import React, { useState } from 'react';
import { UseFormReturn } from 'react-hook-form';
import { Info, Layers, Network, Sliders, Zap, ChevronDown, ChevronUp } from 'lucide-react';
import type { TrainingFormData } from '@/lib/validation/training-schemas';

interface AgentParametersFormProps {
  form: UseFormReturn<TrainingFormData>;
  availableLayers: string[];
}

/**
 * Slider field component for loss weights and training parameters
 */
const SliderField: React.FC<{
  label: string;
  value: number;
  defaultValue: number;
  min: number;
  max: number;
  step: number;
  description: string;
  formatValue?: (val: number) => string;
  onChange: (value: number) => void;
  error?: string;
}> = ({ label, value, defaultValue, min, max, step, description, formatValue, onChange, error }) => {
  const displayValue = formatValue ? formatValue(value ?? defaultValue) : (value ?? defaultValue).toFixed(1);
  
  return (
    <div className="space-y-3">
      <div className="flex items-center justify-between">
        <label className="text-sm font-medium text-gray-900 dark:text-gray-100">
          {label}
        </label>
        <span className="inline-flex items-center px-2.5 py-0.5 rounded-md text-xs font-mono font-medium bg-gray-100 text-gray-800 dark:bg-gray-800 dark:text-gray-200 border border-gray-300 dark:border-gray-600 min-w-[4rem] justify-center tabular-nums">
          {displayValue}
        </span>
      </div>
      <input
        type="range"
        min={min}
        max={max}
        step={step}
        value={value ?? defaultValue}
        onChange={(e) => onChange(parseFloat(e.target.value))}
        title='Slider input'
        className="w-full h-2 bg-gray-200 rounded-lg appearance-none cursor-pointer dark:bg-gray-700 accent-blue-600 hover:accent-blue-700"
      />
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

/**
 * Number input field component
 */
const NumberInput: React.FC<{
  label: string;
  value: number | undefined;
  min: number;
  max: number;
  placeholder?: string;
  description: string;
  onChange: (value: number) => void;
  error?: string;
}> = ({ label, value, min, max, placeholder, description, onChange, error }) => {
  return (
    <div className="space-y-2">
      <label className="text-sm font-medium text-gray-900 dark:text-gray-100">
        {label}
      </label>
      <input
        type="number"
        min={min}
        max={max}
        step={1}
        placeholder={placeholder}
        value={value ?? ''}
        onChange={(e) => {
          const val = parseInt(e.target.value, 10);
          if (!isNaN(val)) onChange(val);
        }}
        className="w-full px-3 py-2 text-sm font-mono bg-white dark:bg-gray-900 border border-gray-300 dark:border-gray-600 rounded-md shadow-sm focus:outline-none focus:ring-2 focus:ring-blue-500 focus:border-blue-500 dark:text-gray-100"
      />
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
        title='Select data layer'
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

/**
 * Badge component
 */
const Badge: React.FC<{ children: React.ReactNode; variant?: 'default' | 'secondary' }> = ({ 
  children, 
  variant = 'default' 
}) => {
  const classes = variant === 'secondary'
    ? 'inline-flex items-center px-2.5 py-0.5 rounded-md text-xs font-medium bg-gray-100 text-gray-700 dark:bg-gray-700 dark:text-gray-300'
    : 'inline-flex items-center px-2.5 py-0.5 rounded-md text-xs font-medium bg-blue-100 text-blue-800 dark:bg-blue-900 dark:text-blue-200';
  
  return <span className={classes}>{children}</span>;
};

/**
 * Collapsible section header component
 */
const SectionHeader: React.FC<{ 
  icon: React.ReactNode; 
  title: string; 
  badge?: string;
  isCollapsed: boolean;
  onToggle: () => void;
}> = ({ icon, title, badge, isCollapsed, onToggle }) => {
  return (
    <button
      type="button"
      onClick={onToggle}
      className="w-full flex items-center justify-between group hover:bg-gray-50 dark:hover:bg-gray-800/50 -mx-2 px-2 py-2 rounded-md transition-colors"
    >
      <div className="flex items-center gap-2">
        <div className="flex items-center gap-2">
          {icon}
          <h4 className="text-sm font-semibold text-gray-900 dark:text-gray-100">{title}</h4>
        </div>
        {badge && <Badge variant="secondary">{badge}</Badge>}
      </div>
      <div className="text-gray-500 dark:text-gray-400 group-hover:text-gray-700 dark:group-hover:text-gray-300 transition-colors">
        {isCollapsed ? (
          <ChevronDown className="w-4 h-4" />
        ) : (
          <ChevronUp className="w-4 h-4" />
        )}
      </div>
    </button>
  );
};

/**
 * Form component for configuring LiVAE model parameters and training hyperparameters
 */
export const AgentParametersForm: React.FC<AgentParametersFormProps> = ({
  form,
  availableLayers,
}) => {
  const { register, watch, setValue, formState: { errors } } = form;
  
  // Get form values
  const formValues = watch('agent_parameters');
  
  // Collapse state for each section
  const [collapsed, setCollapsed] = useState({
    dataLayer: false,
    architecture: false,
    lossWeights: false,
    training: false,
  });

  // Toggle individual section
  const toggleSection = (section: keyof typeof collapsed) => {
    setCollapsed(prev => ({ ...prev, [section]: !prev[section] }));
  };

  // Toggle all sections
  const toggleAll = () => {
    const allCollapsed = Object.values(collapsed).every(v => v);
    setCollapsed({
      dataLayer: !allCollapsed,
      architecture: !allCollapsed,
      lossWeights: !allCollapsed,
      training: !allCollapsed,
    });
  };

  const allCollapsed = Object.values(collapsed).every(v => v);
  
  // Convert learning rate from linear slider to logarithmic scale
  const lrToSlider = (lr: number): number => Math.log10(lr);
  const sliderToLr = (val: number): number => Math.pow(10, val);
  const formatLr = (lr: number): string => lr.toExponential(1);

  // Layer options
  const layerOptions = availableLayers.map(layer => ({
    value: layer,
    label: layer === 'X' ? 'X (Main Layer)' : layer
  }));

  return (
    <div className="bg-white dark:bg-gray-900 border border-gray-200 dark:border-gray-800 rounded-lg shadow-sm">
      {/* Card Header */}
      <div className="px-6 py-5 border-b border-gray-200 dark:border-gray-800">
        <div className="flex items-center justify-between">
          <div>
            <h3 className="text-lg font-semibold text-gray-900 dark:text-gray-100">
              Model Parameters
            </h3>
            <p className="mt-1 text-sm text-gray-600 dark:text-gray-400">
              Configure the LiVAE model architecture and training hyperparameters
            </p>
          </div>
          <button
            type="button"
            onClick={toggleAll}
            className="flex items-center gap-2 px-3 py-1.5 text-sm font-medium text-gray-700 dark:text-gray-300 bg-gray-100 dark:bg-gray-800 hover:bg-gray-200 dark:hover:bg-gray-700 rounded-md transition-colors"
          >
            {allCollapsed ? (
              <>
                <ChevronDown className="w-4 h-4" />
                Expand All
              </>
            ) : (
              <>
                <ChevronUp className="w-4 h-4" />
                Collapse All
              </>
            )}
          </button>
        </div>
      </div>
      
      {/* Card Content */}
      <div className="px-6 py-6 space-y-8">
        {/* ========== Data Layer Selection ========== */}
        <div className="space-y-4">
          <SectionHeader 
            icon={<Layers className="w-4 h-4 text-gray-600 dark:text-gray-400" />}
            title="Data Layer"
            badge="Input Source"
            isCollapsed={collapsed.dataLayer}
            onToggle={() => toggleSection('dataLayer')}
          />
          
          {!collapsed.dataLayer && (
            <SelectField
              label="Data Layer"
              value={formValues?.layer || 'X'}
              options={layerOptions}
              description="Choose which data layer to use for training the model"
              onChange={(value) => setValue('agent_parameters.layer', value)}
              error={errors.agent_parameters?.layer?.message}
            />
          )}
        </div>

        <hr className="border-gray-200 dark:border-gray-700" />

        {/* ========== Architecture Parameters ========== */}
        <div className="space-y-4">
          <SectionHeader 
            icon={<Network className="w-4 h-4 text-gray-600 dark:text-gray-400" />}
            title="Architecture" 
            badge="Network Dimensions"
            isCollapsed={collapsed.architecture}
            onToggle={() => toggleSection('architecture')}
          />
          
          {!collapsed.architecture && (
            <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
              {/* Hidden Dimension */}
              <NumberInput
                label="Hidden Dimension"
                value={formValues?.hidden_dim}
                min={16}
                max={1024}
                placeholder="128"
                description="Width of encoder and decoder layers"
                onChange={(value) => setValue('agent_parameters.hidden_dim', value)}
                error={errors.agent_parameters?.hidden_dim?.message}
              />

              {/* Latent Dimension */}
              <NumberInput
                label="Latent Dimension"
                value={formValues?.latent_dim}
                min={2}
                max={100}
                placeholder="10"
                description="Size of the latent space representation"
                onChange={(value) => setValue('agent_parameters.latent_dim', value)}
                error={errors.agent_parameters?.latent_dim?.message}
              />

              {/* Interpretable Dimension */}
              <NumberInput
                label="Interpretable Dim"
                value={formValues?.i_dim}
                min={2}
                max={50}
                placeholder="5"
                description="Size of interpretable embedding space"
                onChange={(value) => setValue('agent_parameters.i_dim', value)}
                error={errors.agent_parameters?.i_dim?.message}
              />
            </div>
          )}
        </div>

        <hr className="border-gray-200 dark:border-gray-700" />

        {/* ========== Loss Weights ========== */}
        <div className="space-y-6">
          <SectionHeader 
            icon={<Sliders className="w-4 h-4 text-gray-600 dark:text-gray-400" />}
            title="Loss Weights" 
            badge="Regularization"
            isCollapsed={collapsed.lossWeights}
            onToggle={() => toggleSection('lossWeights')}
          />

          {!collapsed.lossWeights && (
            <div className="grid grid-cols-1 lg:grid-cols-2 gap-x-8 gap-y-6">
              {/* Beta (β-VAE) */}
              <SliderField
                label="Beta (β-VAE)"
                value={formValues?.beta}
                defaultValue={1.0}
                min={0.1}
                max={10}
                step={0.1}
                description="KL divergence weight for disentanglement"
                onChange={(value) => setValue('agent_parameters.beta', value)}
                error={errors.agent_parameters?.beta?.message}
              />

              {/* Interpretable Reconstruction */}
              <SliderField
                label="Interpretable Recon"
                value={formValues?.irecon}
                defaultValue={0.0}
                min={0}
                max={10}
                step={0.1}
                description="Weight for embedding reconstruction loss"
                onChange={(value) => setValue('agent_parameters.irecon', value)}
                error={errors.agent_parameters?.irecon?.message}
              />

              {/* Lorentz Distance */}
              <SliderField
                label="Lorentz Distance"
                value={formValues?.lorentz}
                defaultValue={0.0}
                min={0}
                max={10}
                step={0.1}
                description="Hyperbolic geometry regularization penalty"
                onChange={(value) => setValue('agent_parameters.lorentz', value)}
                error={errors.agent_parameters?.lorentz?.message}
              />

              {/* DIP-VAE */}
              <SliderField
                label="DIP-VAE"
                value={formValues?.dip}
                defaultValue={0.0}
                min={0}
                max={10}
                step={0.1}
                description="Disentangled inferred prior weight"
                onChange={(value) => setValue('agent_parameters.dip', value)}
                error={errors.agent_parameters?.dip?.message}
              />

              {/* Total Correlation */}
              <SliderField
                label="Total Correlation"
                value={formValues?.tc}
                defaultValue={0.0}
                min={0}
                max={10}
                step={0.1}
                description="β-TC-VAE total correlation penalty"
                onChange={(value) => setValue('agent_parameters.tc', value)}
                error={errors.agent_parameters?.tc?.message}
              />

              {/* Mutual Information */}
              <SliderField
                label="Mutual Information"
                value={formValues?.info}
                defaultValue={0.0}
                min={0}
                max={10}
                step={0.1}
                description="Mutual information penalty weight"
                onChange={(value) => setValue('agent_parameters.info', value)}
                error={errors.agent_parameters?.info?.message}
              />
            </div>
          )}
        </div>

        <hr className="border-gray-200 dark:border-gray-700" />

        {/* ========== Training Parameters ========== */}
        <div className="space-y-6">
          <SectionHeader 
            icon={<Zap className="w-4 h-4 text-gray-600 dark:text-gray-400" />}
            title="Training" 
            badge="Optimization"
            isCollapsed={collapsed.training}
            onToggle={() => toggleSection('training')}
          />

          {!collapsed.training && (
            <div className="grid grid-cols-1 lg:grid-cols-2 gap-x-8 gap-y-6">
              {/* Learning Rate (Logarithmic Scale) */}
              <div className="space-y-3">
                <div className="flex items-center justify-between">
                  <label className="text-sm font-medium text-gray-900 dark:text-gray-100">
                    Learning Rate
                  </label>
                  <span className="inline-flex items-center px-2.5 py-0.5 rounded-md text-xs font-mono font-medium bg-gray-100 text-gray-800 dark:bg-gray-800 dark:text-gray-200 border border-gray-300 dark:border-gray-600 min-w-[4rem] justify-center tabular-nums">
                    {formatLr(formValues?.lr ?? 1e-3)}
                  </span>
                </div>
                <input
                  type="range"
                  min={-6}
                  max={-1}
                  step={0.1}
                  value={lrToSlider(formValues?.lr ?? 1e-3)}
                  onChange={(e) => setValue('agent_parameters.lr', sliderToLr(parseFloat(e.target.value)))}
                  title='Slider input'
                  className="w-full h-2 bg-gray-200 rounded-lg appearance-none cursor-pointer dark:bg-gray-700 accent-blue-600 hover:accent-blue-700"
                />
                <p className="text-xs text-gray-600 dark:text-gray-400">
                  Adam optimizer learning rate (logarithmic scale)
                </p>
                {errors.agent_parameters?.lr?.message && (
                  <p className="text-xs text-red-600 dark:text-red-400 flex items-center gap-1">
                    <Info className="w-3 h-3" />
                    {errors.agent_parameters.lr.message}
                  </p>
                )}
              </div>

              {/* Data Percentage */}
              <SliderField
                label="Data Percentage"
                value={formValues?.percent}
                defaultValue={1.0}
                min={0.01}
                max={1.0}
                step={0.01}
                description="Fraction of dataset to use for training"
                formatValue={(val) => `${(val * 100).toFixed(0)}%`}
                onChange={(value) => setValue('agent_parameters.percent', value)}
                error={errors.agent_parameters?.percent?.message}
              />
            </div>
          )}
        </div>
      </div>
    </div>
  );
};
