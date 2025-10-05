
// src/components/training/TrainingConfigForm.tsx
'use client';

import React, { createContext, useContext, useEffect, useState } from 'react';
import { useRouter } from 'next/navigation';
import { useForm } from 'react-hook-form';
import { zodResolver } from '@hookform/resolvers/zod';
import { Play, RotateCcw, AlertCircle, CheckCircle2 } from 'lucide-react';
import { DatasetSummaryCard } from './DatasetSummaryCard';
import { AgentParametersForm } from './AgentParametersForm';
import { QCParametersForm } from './QCParametersForm';
import { TrainingSettingsForm } from './TrainingSettingsForm';
import { useTrainingConfig } from '@/lib/hooks/useTrainingConfig';
import { api } from '@/lib/api/endpoints';
import { trainingFormSchema, type TrainingFormData } from '@/lib/validation/training-schemas';
import { DEFAULT_AGENT_PARAMETERS, DEFAULT_QC_PARAMS, type AnnDataSummary } from '@/types/index';

// ========================================
// Types
// ========================================

type TabValue = 'parameters' | 'qc' | 'training';

// ========================================
// Default Form Values
// ========================================

const DEFAULT_FORM_VALUES: TrainingFormData = {
  training_config: {
    epochs: 1000,
    apply_qc: false,
    qc_params: DEFAULT_QC_PARAMS,
  },
  agent_parameters: DEFAULT_AGENT_PARAMETERS,
};

// ========================================
// Main Component
// ========================================

export const TrainingConfigForm: React.FC = () => {
  const router = useRouter();
  const [dataSummary, setDataSummary] = useState<AnnDataSummary | null>(null);
  const [loading, setLoading] = useState(true);
  const [activeTab, setActiveTab] = useState<TabValue>('parameters');
  
  const { isSubmitting, error, isSuccess, submitTrainingConfig, reset } = useTrainingConfig();

  const form = useForm<TrainingFormData>({
    resolver: zodResolver(trainingFormSchema),
    defaultValues: DEFAULT_FORM_VALUES,
    mode: 'onChange',
  });

  // Load data summary
  useEffect(() => {
    const checkDataStatus = async () => {
      try {
        const summary = await api.getDataSummary();
        setDataSummary(summary);
      } catch (error) {
        console.error('No data loaded:', error);
      } finally {
        setLoading(false);
      }
    };
    checkDataStatus();
  }, []);

  // Redirect on success
  useEffect(() => {
    if (isSuccess) {
      const timer = setTimeout(() => router.push('/training/monitor'), 2000);
      return () => clearTimeout(timer);
    }
  }, [isSuccess, router]);

  const onSubmit = (data: TrainingFormData) => {
    submitTrainingConfig(data);
  };

  const handleReset = () => {
    form.reset(DEFAULT_FORM_VALUES);
    reset();
  };

  if (loading) {
    return (
      <div className="text-center py-12">
        <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-blue-600 mx-auto mb-4"></div>
        <p className="text-gray-600 dark:text-gray-400">Loading configuration...</p>
      </div>
    );
  }

  if (!dataSummary) {
    return (
      <Alert variant="error">
        No dataset found. Please upload your data first.
      </Alert>
    );
  }

  const availableLayers = ['X', ...Object.keys(dataSummary.layers)];

  return (
    <>
      <DatasetSummaryCard dataSummary={dataSummary} />

      <form 
        onSubmit={form.handleSubmit(onSubmit)} 
        className="space-y-6"
        onKeyDown={(e) => {
          // Prevent Enter key from submitting form when editing fields
          if (e.key === 'Enter' && e.target !== e.currentTarget) {
            e.preventDefault();
          }
        }}
      >
        {/* Tabs */}
        <Tabs value={activeTab} onValueChange={setActiveTab}>
          <TabsList>
            <TabsTrigger value="parameters">
              Model
            </TabsTrigger>
            <TabsTrigger value="qc">
              Quality Control
            </TabsTrigger>
            <TabsTrigger value="training">
              Training
            </TabsTrigger>
          </TabsList>

          <div className="mt-4">
            <TabsContent value="parameters">
              <AgentParametersForm form={form} availableLayers={availableLayers} />
            </TabsContent>

            <TabsContent value="qc">
              <QCParametersForm form={form} />
            </TabsContent>

            <TabsContent value="training">
              <TrainingSettingsForm form={form} />
            </TabsContent>
          </div>
        </Tabs>

        {/* Error Alert */}
        {error && (
          <Alert variant="error">
            {error}
          </Alert>
        )}

        {/* Success Alert */}
        {isSuccess && (
          <Alert variant="success">
            ✓ Training started successfully! Redirecting to monitor...
          </Alert>
        )}

        {/* Action Buttons */}
        <div className="flex flex-wrap gap-3 justify-end">
          <Button 
            type="button" 
            variant="outline" 
            onClick={handleReset}
            disabled={isSubmitting}
          >
            <RotateCcw className="mr-2 h-4 w-4" />
            Reset Defaults
          </Button>
          
          <Button type="submit" variant="primary" disabled={isSubmitting}>
            {isSubmitting ? (
              <>
                <div className="animate-spin rounded-full h-4 w-4 border-b-2 border-white mr-2" />
                Starting...
              </>
            ) : (
              <>
                <Play className="mr-2 h-4 w-4" />
                Start Training
              </>
            )}
          </Button>
        </div>
      </form>
    </>
  );
};

// ========================================
// Custom Components
// ========================================

// Button Component
const Button: React.FC<{
  type?: 'button' | 'submit' | 'reset';
  variant?: 'primary' | 'outline';
  disabled?: boolean;
  onClick?: () => void;
  children: React.ReactNode;
}> = ({ type = 'button', variant = 'primary', disabled = false, onClick, children }) => {
  const baseStyles = "inline-flex items-center justify-center px-4 py-2 text-sm font-medium rounded-md transition-colors focus:outline-none focus:ring-2 focus:ring-offset-2 disabled:opacity-50 disabled:cursor-not-allowed";
  
  const variantStyles = {
    primary: "bg-blue-600 hover:bg-blue-700 text-white focus:ring-blue-500 disabled:hover:bg-blue-600",
    outline: "bg-white dark:bg-gray-800 border border-gray-300 dark:border-gray-600 text-gray-700 dark:text-gray-300 hover:bg-gray-50 dark:hover:bg-gray-700 focus:ring-blue-500"
  };

  const handleClick = (e: React.MouseEvent<HTMLButtonElement>) => {
    if (type === 'button') {
      e.preventDefault(); // Prevent form submission for type="button"
    }
    onClick?.();
  };

  return (
    <button
      type={type}
      onClick={handleClick}
      disabled={disabled}
      className={`${baseStyles} ${variantStyles[variant]}`}
    >
      {children}
    </button>
  );
};

// Alert Component
const Alert: React.FC<{
  variant: 'error' | 'success' | 'info';
  children: React.ReactNode;
}> = ({ variant, children }) => {
  const styles = {
    error: {
      container: "bg-red-50 dark:bg-red-950 border-red-200 dark:border-red-800",
      text: "text-red-800 dark:text-red-200",
      icon: AlertCircle,
      iconColor: "text-red-600 dark:text-red-400"
    },
    success: {
      container: "bg-green-50 dark:bg-green-950 border-green-200 dark:border-green-800",
      text: "text-green-800 dark:text-green-200",
      icon: CheckCircle2,
      iconColor: "text-green-600 dark:text-green-400"
    },
    info: {
      container: "bg-blue-50 dark:bg-blue-950 border-blue-200 dark:border-blue-800",
      text: "text-blue-800 dark:text-blue-200",
      icon: AlertCircle,
      iconColor: "text-blue-600 dark:text-blue-400"
    }
  };

  const style = styles[variant];
  const Icon = style.icon;

  return (
    <div className={`flex gap-3 rounded-lg border p-4 ${style.container}`}>
      <Icon className={`h-5 w-5 flex-shrink-0 mt-0.5 ${style.iconColor}`} />
      <p className={`text-sm ${style.text}`}>
        {children}
      </p>
    </div>
  );
};

// ========================================
// Tabs Components with Context
// ========================================

const TabsContext = createContext<{
  value: TabValue;
  onValueChange: (value: TabValue) => void;
} | null>(null);

const useTabsContext = () => {
  const context = useContext(TabsContext);
  if (!context) {
    throw new Error('Tabs components must be used within Tabs');
  }
  return context;
};

const Tabs: React.FC<{
  value: TabValue;
  onValueChange: (value: TabValue) => void;
  children: React.ReactNode;
}> = ({ value, onValueChange, children }) => {
  return (
    <TabsContext.Provider value={{ value, onValueChange }}>
      <div className="space-y-4">{children}</div>
    </TabsContext.Provider>
  );
};

const TabsList: React.FC<{
  children: React.ReactNode;
}> = ({ children }) => {
  return (
    <div className="inline-flex bg-gray-100 dark:bg-gray-800 rounded-lg p-1 w-full">
      <div className="grid w-full grid-cols-3 gap-1">
        {children}
      </div>
    </div>
  );
};

const TabsTrigger: React.FC<{
  value: TabValue;
  children: React.ReactNode;
}> = ({ value, children }) => {
  const { value: activeValue, onValueChange } = useTabsContext();
  const active = value === activeValue;

  const handleClick = (e: React.MouseEvent<HTMLButtonElement>) => {
    e.preventDefault(); 
    e.stopPropagation();
    onValueChange(value);
  };

  return (
    <button
      type="button"
      onClick={handleClick}
      className={`px-3 py-2 text-sm font-medium rounded-md transition-colors focus:outline-none focus:ring-2 focus:ring-blue-500 focus:ring-offset-2 ${
        active
          ? 'bg-white dark:bg-gray-900 text-gray-900 dark:text-gray-100 shadow-sm'
          : 'text-gray-600 dark:text-gray-400 hover:text-gray-900 dark:hover:text-gray-100'
      }`}
    >
      {children}
    </button>
  );
};

const TabsContent: React.FC<{
  value: TabValue;
  children: React.ReactNode;
}> = ({ value, children }) => {
  const { value: activeValue } = useTabsContext();
  
  if (value !== activeValue) return null;
  
  return (
    <div role="tabpanel" className="space-y-4">
      {children}
    </div>
  );
};
