
// src/components/training/AgentParametersForm.tsx
'use client';

import React from 'react';
import { UseFormReturn } from 'react-hook-form';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { FormControl, FormDescription, FormField, FormItem, FormLabel, FormMessage } from '@/components/ui/form';
import { Input } from '@/components/ui/input';
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select';
import { Slider } from '@/components/ui/slider';
import { Badge } from '@/components/ui/badge';
import { Separator } from '@/components/ui/separator';
import type { TrainingFormData } from '@/lib/validation/training-schemas';

interface AgentParametersFormProps {
  form: UseFormReturn<TrainingFormData>;
  availableLayers: string[];
}

export const AgentParametersForm: React.FC<AgentParametersFormProps> = ({
  form,
  availableLayers,
}) => {
  // Helper to format scientific notation for learning rate
  const formatLearningRate = (value: number): string => {
    return value.toExponential(1);
  };

  return (
    <Card>
      <CardHeader>
        <CardTitle>Model Parameters</CardTitle>
        <CardDescription>
          Configure the LiVAE model architecture and training hyperparameters
        </CardDescription>
      </CardHeader>
      <CardContent className="space-y-6">
        {/* Data Layer Selection */}
        <FormField
          control={form.control}
          name="agent_parameters.layer"
          render={({ field }) => (
            <FormItem>
              <FormLabel>Data Layer</FormLabel>
              <Select onValueChange={field.onChange} defaultValue={field.value}>
                <FormControl>
                  <SelectTrigger>
                    <SelectValue placeholder="Select data layer to use" />
                  </SelectTrigger>
                </FormControl>
                <SelectContent>
                  {availableLayers.map((layer) => (
                    <SelectItem key={layer} value={layer}>
                      {layer === 'X' ? 'X (Main Layer)' : layer}
                    </SelectItem>
                  ))}
                </SelectContent>
              </Select>
              <FormDescription>
                Choose which data layer to use for training
              </FormDescription>
              <FormMessage />
            </FormItem>
          )}
        />

        <Separator />

        {/* Architecture Parameters */}
        <div className="space-y-4">
          <div className="flex items-center space-x-2">
            <h4 className="font-medium">Architecture</h4>
            <Badge variant="outline">Network Dimensions</Badge>
          </div>
          
          <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            <FormField
              control={form.control}
              name="agent_parameters.hidden_dim"
              render={({ field }) => (
                <FormItem>
                  <FormLabel>Hidden Dimension</FormLabel>
                  <FormControl>
                    <Input
                      type="number"
                      {...field}
                      onChange={(e) => field.onChange(parseInt(e.target.value) || 0)}
                    />
                  </FormControl>
                  <FormDescription>Encoder/decoder width</FormDescription>
                  <FormMessage />
                </FormItem>
              )}
            />

            <FormField
              control={form.control}
              name="agent_parameters.latent_dim"
              render={({ field }) => (
                <FormItem>
                  <FormLabel>Latent Dimension</FormLabel>
                  <FormControl>
                    <Input
                      type="number"
                      {...field}
                      onChange={(e) => field.onChange(parseInt(e.target.value) || 0)}
                    />
                  </FormControl>
                  <FormDescription>Latent space size</FormDescription>
                  <FormMessage />
                </FormItem>
              )}
            />

            <FormField
              control={form.control}
              name="agent_parameters.i_dim"
              render={({ field }) => (
                <FormItem>
                  <FormLabel>Interpretable Dimension</FormLabel>
                  <FormControl>
                    <Input
                      type="number"
                      {...field}
                      onChange={(e) => field.onChange(parseInt(e.target.value) || 0)}
                    />
                  </FormControl>
                  <FormDescription>Interpretable embedding size</FormDescription>
                  <FormMessage />
                </FormItem>
              )}
            />
          </div>
        </div>

        <Separator />

        {/* Loss Weights */}
        <div className="space-y-4">
          <div className="flex items-center space-x-2">
            <h4 className="font-medium">Loss Weights</h4>
            <Badge variant="outline">Regularization</Badge>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
            <FormField
              control={form.control}
              name="agent_parameters.beta"
              render={({ field }) => (
                <FormItem className="space-y-3">
                  <FormLabel>Beta (β-VAE Weight): {field.value}</FormLabel>
                  <FormControl>
                    <Slider
                      min={0.1}
                      max={10}
                      step={0.1}
                      value={[field.value]}
                      onValueChange={(values) => field.onChange(values[0])}
                    />
                  </FormControl>
                  <FormDescription>KL divergence weight</FormDescription>
                  <FormMessage />
                </FormItem>
              )}
            />

            <FormField
              control={form.control}
              name="agent_parameters.irecon"
              render={({ field }) => (
                <FormItem className="space-y-3">
                  <FormLabel>Interpretable Reconstruction: {field.value}</FormLabel>
                  <FormControl>
                    <Slider
                      min={0}
                      max={10}
                      step={0.1}
                      value={[field.value]}
                      onValueChange={(values) => field.onChange(values[0])}
                    />
                  </FormControl>
                  <FormDescription>Interpretable embedding reconstruction weight</FormDescription>
                  <FormMessage />
                </FormItem>
              )}
            />

            <FormField
              control={form.control}
              name="agent_parameters.lorentz"
              render={({ field }) => (
                <FormItem className="space-y-3">
                  <FormLabel>Lorentz Distance: {field.value}</FormLabel>
                  <FormControl>
                    <Slider
                      min={0}
                      max={10}
                      step={0.1}
                      value={[field.value]}
                      onValueChange={(values) => field.onChange(values[0])}
                    />
                  </FormControl>
                  <FormDescription>Hyperbolic geometry regularization</FormDescription>
                  <FormMessage />
                </FormItem>
              )}
            />

            <FormField
              control={form.control}
              name="agent_parameters.dip"
              render={({ field }) => (
                <FormItem className="space-y-3">
                  <FormLabel>DIP-VAE: {field.value}</FormLabel>
                  <FormControl>
                    <Slider
                      min={0}
                      max={10}
                      step={0.1}
                      value={[field.value]}
                      onValueChange={(values) => field.onChange(values[0])}
                    />
                  </FormControl>
                  <FormDescription>Disentangled inferred prior weight</FormDescription>
                  <FormMessage />
                </FormItem>
              )}
            />

            <FormField
              control={form.control}
              name="agent_parameters.tc"
              render={({ field }) => (
                <FormItem className="space-y-3">
                  <FormLabel>Total Correlation: {field.value}</FormLabel>
                  <FormControl>
                    <Slider
                      min={0}
                      max={10}
                      step={0.1}
                      value={[field.value]}
                      onValueChange={(values) => field.onChange(values[0])}
                    />
                  </FormControl>
                  <FormDescription>β-TC-VAE total correlation penalty</FormDescription>
                  <FormMessage />
                </FormItem>
              )}
            />

            <FormField
              control={form.control}
              name="agent_parameters.info"
              render={({ field }) => (
                <FormItem className="space-y-3">
                  <FormLabel>Mutual Information: {field.value}</FormLabel>
                  <FormControl>
                    <Slider
                      min={0}
                      max={10}
                      step={0.1}
                      value={[field.value]}
                      onValueChange={(values) => field.onChange(values[0])}
                    />
                  </FormControl>
                  <FormDescription>Mutual information penalty weight</FormDescription>
                  <FormMessage />
                </FormItem>
              )}
            />
          </div>
        </div>

        <Separator />

        {/* Training Parameters */}
        <div className="space-y-4">
          <div className="flex items-center space-x-2">
            <h4 className="font-medium">Training</h4>
            <Badge variant="outline">Optimization</Badge>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
            <FormField
              control={form.control}
              name="agent_parameters.lr"
              render={({ field }) => (
                <FormItem className="space-y-3">
                  <FormLabel>Learning Rate: {formatLearningRate(field.value)}</FormLabel>
                  <FormControl>
                    <Slider
                      min={1e-6}
                      max={1e-1}
                      step={1e-6}
                      value={[field.value]}
                      onValueChange={(values) => field.onChange(values[0])}
                    />
                  </FormControl>
                  <FormDescription>Adam optimizer learning rate</FormDescription>
                  <FormMessage />
                </FormItem>
              )}
            />

            <FormField
              control={form.control}
              name="agent_parameters.percent"
              render={({ field }) => (
                <FormItem className="space-y-3">
                  <FormLabel>Data Percentage: {(field.value * 100).toFixed(1)}%</FormLabel>
                  <FormControl>
                    <Slider
                      min={0.001}
                      max={1.0}
                      step={0.001}
                      value={[field.value]}
                      onValueChange={(values) => field.onChange(values[0])}
                    />
                  </FormControl>
                  <FormDescription>Fraction of data to use for training</FormDescription>
                  <FormMessage />
                </FormItem>
              )}
            />
          </div>
        </div>
      </CardContent>
    </Card>
  );
};
