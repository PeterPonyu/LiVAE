
// src/app/training/configure/page.tsx
'use client';

import React, { useEffect, useState } from 'react';
import { useRouter } from 'next/navigation';
import { useForm } from 'react-hook-form';
import { zodResolver } from '@hookform/resolvers/zod';
import { ArrowLeft, Play, Save } from 'lucide-react';
import Link from 'next/link';
import { Button } from '@/components/ui/button';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Form, FormControl, FormField, FormItem, FormLabel, FormMessage } from '@/components/ui/form';
import { Input } from '@/components/ui/input';
import { Alert, AlertDescription } from '@/components/ui/alert';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import { AgentParametersForm } from '@/components/training/AgentParametersForm';
import { QCParametersForm } from '@/components/training/QCParametersForm';
import { useTrainingConfig } from '@/lib/hooks/useTrainingConfig';
import { api } from '@/lib/api/endpoints';
import { trainingFormSchema, type TrainingFormData } from '@/lib/validation/training-schemas';
import { DEFAULT_AGENT_PARAMETERS, DEFAULT_QC_PARAMS, type AnnDataSummary } from '@/types/index';

export default function TrainingConfigPage() {
  const router = useRouter();
  const [dataSummary, setDataSummary] = useState<AnnDataSummary | null>(null);
  const [loading, setLoading] = useState(true);
  
  const { isSubmitting, error, isSuccess, submitTrainingConfig, reset } = useTrainingConfig();

  // Initialize form with default values - Fix: proper type inference
  const form = useForm<TrainingFormData>({
    resolver: zodResolver(trainingFormSchema),
    defaultValues: {
      training_config: {
        epochs: 1000,
        apply_qc: false,
        qc_params: DEFAULT_QC_PARAMS,
      },
      agent_parameters: DEFAULT_AGENT_PARAMETERS,
    },
  });

  // Check if data is loaded
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

  // Handle form submission - Fix: proper typing
  const onSubmit = (data: TrainingFormData) => {
    submitTrainingConfig(data);
  };

  // Handle success - redirect to monitoring
  useEffect(() => {
    if (isSuccess) {
      setTimeout(() => {
        router.push('/training/monitor');
      }, 2000);
    }
  }, [isSuccess, router]);

  // Handle reset to defaults
  const handleReset = () => {
    form.reset({
      training_config: {
        epochs: 1000,
        apply_qc: false,
        qc_params: DEFAULT_QC_PARAMS,
      },
      agent_parameters: DEFAULT_AGENT_PARAMETERS,
    });
    reset();
  };

  if (loading) {
    return (
      <div className="container mx-auto py-8 px-4 max-w-4xl">
        <div className="text-center">
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

  // Get available layers from data summary
  const availableLayers = ['X', ...Object.keys(dataSummary.layers)];

  return (
    <div className="container mx-auto py-8 px-4 max-w-6xl">
      {/* Header */}
      <div className="flex items-center justify-between mb-6">
        <div>
          <h1 className="text-3xl font-bold">Training Configuration</h1>
          <p className="text-gray-600 mt-1">
            Configure model parameters and quality control filters
          </p>
        </div>
        <Button variant="outline" asChild>
          <Link href="/">
            <ArrowLeft className="mr-2 h-4 w-4" />
            Back to Home
          </Link>
        </Button>
      </div>

      {/* Dataset Summary */}
      <Card className="mb-6">
        <CardHeader>
          <CardTitle>Dataset Summary</CardTitle>
          <CardDescription>Currently loaded dataset</CardDescription>
        </CardHeader>
        <CardContent>
          <div className="grid grid-cols-2 md:grid-cols-4 gap-4 text-sm">
            <div>
              <div className="font-medium">Cells</div>
              <div className="text-gray-600">{dataSummary.shape.n_obs.toLocaleString()}</div>
            </div>
            <div>
              <div className="font-medium">Genes</div>
              <div className="text-gray-600">{dataSummary.shape.n_vars.toLocaleString()}</div>
            </div>
            <div>
              <div className="font-medium">Layers</div>
              <div className="text-gray-600">{availableLayers.join(', ')}</div>
            </div>
            <div>
              <div className="font-medium">File</div>
              <div className="text-gray-600">{dataSummary.filename || 'Unknown'}</div>
            </div>
          </div>
        </CardContent>
      </Card>

      {/* Configuration Form */}
      <Form {...form}>
        <form onSubmit={form.handleSubmit(onSubmit)} className="space-y-6">
          <Tabs defaultValue="parameters" className="space-y-6">
            <TabsList className="grid w-full grid-cols-3">
              <TabsTrigger value="parameters">Model Parameters</TabsTrigger>
              <TabsTrigger value="qc">Quality Control</TabsTrigger>
              <TabsTrigger value="training">Training Settings</TabsTrigger>
            </TabsList>

            {/* Model Parameters Tab */}
            <TabsContent value="parameters">
              <AgentParametersForm form={form} availableLayers={availableLayers} />
            </TabsContent>

            {/* QC Parameters Tab */}
            <TabsContent value="qc">
              <QCParametersForm form={form} />
            </TabsContent>

            {/* Training Settings Tab */}
            <TabsContent value="training">
              <Card>
                <CardHeader>
                  <CardTitle>Training Settings</CardTitle>
                  <CardDescription>
                    Configure the training process
                  </CardDescription>
                </CardHeader>
                <CardContent className="space-y-4">
                  <FormField
                    control={form.control}
                    name="training_config.epochs"
                    render={({ field }) => (
                      <FormItem>
                        <FormLabel>Number of Epochs</FormLabel>
                        <FormControl>
                          <Input
                            type="number"
                            {...field}
                            onChange={(e) => field.onChange(parseInt(e.target.value) || 0)}
                          />
                        </FormControl>
                        <p className="text-sm text-gray-600">
                          Number of training iterations (recommended: 1000-5000)
                        </p>
                        <FormMessage />
                      </FormItem>
                    )}
                  />

                  <Alert>
                    <AlertDescription>
                      Training time depends on dataset size and number of epochs. 
                      A typical dataset with ~10,000 cells takes 5-15 minutes per 1000 epochs.
                    </AlertDescription>
                  </Alert>
                </CardContent>
              </Card>
            </TabsContent>
          </Tabs>

          {/* Error Display */}
          {error && (
            <Alert variant="destructive">
              <AlertDescription>{error}</AlertDescription>
            </Alert>
          )}

          {/* Success Display */}
          {isSuccess && (
            <Alert className="bg-green-50 border-green-200">
              <AlertDescription className="text-green-800">
                Training started successfully! Redirecting to monitoring page...
              </AlertDescription>
            </Alert>
          )}

          {/* Action Buttons */}
          <div className="flex flex-col sm:flex-row gap-4 justify-end">
            <Button type="button" variant="outline" onClick={handleReset}>
              <Save className="mr-2 h-4 w-4" />
              Reset to Defaults
            </Button>
            
            <Button type="submit" disabled={isSubmitting}>
              <Play className="mr-2 h-4 w-4" />
              {isSubmitting ? 'Starting Training...' : 'Start Training'}
            </Button>
          </div>
        </form>
      </Form>
    </div>
  );
}
