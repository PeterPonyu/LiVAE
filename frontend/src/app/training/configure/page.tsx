
// src/app/training/configure/page.tsx
'use client';

import React, { useEffect, useState } from 'react';
import { useRouter } from 'next/navigation';
import { useForm } from 'react-hook-form';
import { zodResolver } from '@hookform/resolvers/zod';
import { ArrowLeft, Play, RotateCcw, Settings } from 'lucide-react';
import Link from 'next/link';
import { Button } from '@/components/ui/button';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Form, FormControl, FormField, FormItem, FormLabel, FormMessage } from '@/components/ui/form';
import { Input } from '@/components/ui/input';
import { Alert, AlertDescription } from '@/components/ui/alert';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import { Badge } from '@/components/ui/badge';
import { AgentParametersForm } from '@/components/training/AgentParametersForm';
import { QCParametersForm } from '@/components/training/QCParametersForm';
import { useTrainingConfig } from '@/lib/hooks/useTrainingConfig';
import { api } from '@/lib/api/endpoints';
import { trainingFormSchema, type TrainingFormData } from '@/lib/validation/training-schemas';
import { DEFAULT_AGENT_PARAMETERS, DEFAULT_QC_PARAMS, type AnnDataSummary } from '@/types/index';

// Centralized default form values
const DEFAULT_FORM_VALUES: TrainingFormData = {
  training_config: {
    epochs: 1000,
    apply_qc: false,
    qc_params: DEFAULT_QC_PARAMS,
  },
  agent_parameters: DEFAULT_AGENT_PARAMETERS,
};

export default function TrainingConfigPage() {
  const router = useRouter();
  const [dataSummary, setDataSummary] = useState<AnnDataSummary | null>(null);
  const [loading, setLoading] = useState(true);
  
  const { isSubmitting, error, isSuccess, submitTrainingConfig, reset } = useTrainingConfig();

  const form = useForm<TrainingFormData>({
    resolver: zodResolver(trainingFormSchema),
    defaultValues: DEFAULT_FORM_VALUES,
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
      <div className="container mx-auto py-8 px-4 max-w-4xl">
        <div className="text-center">
          <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-blue-600 mx-auto mb-4"></div>
          <p>Loading configuration...</p>
        </div>
      </div>
    );
  }

  if (!dataSummary) {
    return (
      <div className="container mx-auto py-8 px-4 max-w-4xl">
        <Alert variant="destructive">
          <AlertDescription>
            No dataset found. Please upload your data first.
          </AlertDescription>
        </Alert>
        <div className="mt-4">
          <Button asChild>
            <Link href="/upload">
              <ArrowLeft className="mr-2 h-4 w-4" />
              Upload Data
            </Link>
          </Button>
        </div>
      </div>
    );
  }

  const availableLayers = ['X', ...Object.keys(dataSummary.layers)];

  return (
    <div className="container mx-auto py-8 px-4 max-w-6xl">
      {/* Header */}
      <div className="flex items-center justify-between">
        <div>
          <h1 className="text-3xl font-bold flex items-center gap-3">
            <Settings className="h-8 w-8" />
            Training Configuration
          </h1>
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

      {/* Dataset Summary - Compact */}
      <Card className="mb-6">
        <CardHeader className="pb-3">
          <div className="flex items-center justify-between">
            <CardTitle className="text-lg">Dataset Overview</CardTitle>
            <Badge variant="outline">
              {dataSummary.shape.n_obs.toLocaleString()} cells
            </Badge>
          </div>
        </CardHeader>
        <CardContent>
          <div className="grid grid-cols-2 sm:grid-cols-4 gap-4 text-sm">
            <div className="flex flex-col">
              <span className="text-gray-600">Cells</span>
              <span className="font-semibold">{dataSummary.shape.n_obs.toLocaleString()}</span>
            </div>
            <div className="flex flex-col">
              <span className="text-gray-600">Genes</span>
              <span className="font-semibold">{dataSummary.shape.n_vars.toLocaleString()}</span>
            </div>
            <div className="flex flex-col">
              <span className="text-gray-600">Layers</span>
              <span className="font-semibold">{availableLayers.length}</span>
            </div>
            <div className="flex flex-col">
              <span className="text-gray-600">File</span>
              <span className="font-semibold truncate" title={dataSummary.filename}>
                {dataSummary.filename || 'Unknown'}
              </span>
            </div>
          </div>
        </CardContent>
      </Card>

      {/* Configuration Form */}
      <Form {...form}>
        <form onSubmit={form.handleSubmit(onSubmit)} className="space-y-6">
          <Tabs defaultValue="parameters" className="space-y-4">
            <TabsList className="grid w-full grid-cols-3">
              <TabsTrigger value="parameters">Model</TabsTrigger>
              <TabsTrigger value="qc">Quality Control</TabsTrigger>
              <TabsTrigger value="training">Training</TabsTrigger>
            </TabsList>

            {/* Model Parameters */}
            <TabsContent value="parameters">
              <AgentParametersForm form={form} availableLayers={availableLayers} />
            </TabsContent>

            {/* QC Parameters */}
            <TabsContent value="qc">
              <QCParametersForm form={form} />
            </TabsContent>

            {/* Training Settings */}
            <TabsContent value="training">
              <Card>
                <CardHeader>
                  <CardTitle>Training Settings</CardTitle>
                  <CardDescription>
                    Control training duration and monitoring
                  </CardDescription>
                </CardHeader>
                <CardContent className="space-y-4">
                  <FormField
                    control={form.control}
                    name="training_config.epochs"
                    render={({ field }) => (
                      <FormItem>
                        <FormLabel>Epochs</FormLabel>
                        <FormControl>
                          <Input
                            type="number"
                            min={100}
                            max={10000}
                            {...field}
                            onChange={(e) => field.onChange(parseInt(e.target.value) || 0)}
                          />
                        </FormControl>
                        <p className="text-xs text-gray-600">
                          Recommended: 1000-5000 epochs. More epochs = better convergence but longer training time.
                        </p>
                        <FormMessage />
                      </FormItem>
                    )}
                  />

                  <div className="bg-blue-50 p-3 rounded-lg border border-blue-200 text-sm">
                    <p className="text-blue-900">
                      <strong>Estimated time:</strong> ~5-15 minutes per 1000 epochs for datasets with 10,000 cells
                    </p>
                  </div>
                </CardContent>
              </Card>
            </TabsContent>
          </Tabs>

          {/* Status Messages */}
          {error && (
            <Alert variant="destructive">
              <AlertDescription>{error}</AlertDescription>
            </Alert>
          )}

          {isSuccess && (
            <Alert className="bg-green-50 border-green-200">
              <AlertDescription className="text-green-800">
                ✓ Training started successfully! Redirecting to monitor...
              </AlertDescription>
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
            
            <Button type="submit" disabled={isSubmitting}>
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
      </Form>
    </div>
  );
}
