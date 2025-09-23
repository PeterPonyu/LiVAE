
// src/components/training/QCParametersForm.tsx
'use client';

import React from 'react';
import { UseFormReturn } from 'react-hook-form';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { FormControl, FormDescription, FormField, FormItem, FormLabel, FormMessage } from '@/components/ui/form';
import { Input } from '@/components/ui/input';
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select';
import { Switch } from '@/components/ui/switch';
import { Badge } from '@/components/ui/badge';
import { Alert, AlertDescription } from '@/components/ui/alert';
import { Info } from 'lucide-react';
import type { TrainingFormData } from '@/lib/validation/training-schemas';

interface QCParametersFormProps {
  form: UseFormReturn<TrainingFormData>;
}

export const QCParametersForm: React.FC<QCParametersFormProps> = ({ form }) => {
  const applyQC = form.watch('training_config.apply_qc');

  return (
    <Card>
      <CardHeader>
        <CardTitle>Quality Control Filters</CardTitle>
        <CardDescription>
          Apply preprocessing filters to remove low-quality cells and genes
        </CardDescription>
      </CardHeader>
      <CardContent className="space-y-6">
        {/* Enable/Disable QC */}
        <FormField
          control={form.control}
          name="training_config.apply_qc"
          render={({ field }) => (
            <FormItem className="flex flex-row items-center justify-between rounded-lg border p-4">
              <div className="space-y-0.5">
                <FormLabel className="text-base">Apply Quality Control</FormLabel>
                <FormDescription>
                  Enable preprocessing filters to improve data quality
                </FormDescription>
              </div>
              <FormControl>
                <Switch
                  checked={field.value}
                  onCheckedChange={field.onChange}
                />
              </FormControl>
            </FormItem>
          )}
        />

        {applyQC && (
          <>
            <Alert>
              <Info className="h-4 w-4" />
              <AlertDescription>
                QC filtering will be applied before training. Cells and genes not meeting 
                the criteria below will be removed from the dataset.
              </AlertDescription>
            </Alert>

            {/* Species Selection */}
            <FormField
              control={form.control}
              name="training_config.qc_params.species"
              render={({ field }) => (
                <FormItem>
                  <FormLabel>Species</FormLabel>
                  <Select onValueChange={field.onChange} defaultValue={field.value}>
                    <FormControl>
                      <SelectTrigger>
                        <SelectValue placeholder="Select species" />
                      </SelectTrigger>
                    </FormControl>
                    <SelectContent>
                      <SelectItem value="human">Human</SelectItem>
                      <SelectItem value="mouse">Mouse</SelectItem>
                    </SelectContent>
                  </Select>
                  <FormDescription>
                    Species type for mitochondrial gene detection
                  </FormDescription>
                  <FormMessage />
                </FormItem>
              )}
            />

            {/* Cell Filters */}
            <div className="space-y-4">
              <div className="flex items-center space-x-2">
                <h4 className="font-medium">Cell Quality Filters</h4>
                <Badge variant="outline">Per Cell</Badge>
              </div>

              <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                <FormField
                  control={form.control}
                  name="training_config.qc_params.min_genes_per_cell"
                  render={({ field }) => (
                    <FormItem>
                      <FormLabel>Minimum Genes per Cell</FormLabel>
                      <FormControl>
                        <Input
                          type="number"
                          placeholder="e.g., 200"
                          {...field}
                          onChange={(e) => field.onChange(e.target.value ? parseInt(e.target.value) : undefined)}
                        />
                      </FormControl>
                      <FormDescription>
                        Remove cells with fewer detected genes
                      </FormDescription>
                      <FormMessage />
                    </FormItem>
                  )}
                />

                <FormField
                  control={form.control}
                  name="training_config.qc_params.max_genes_per_cell"
                  render={({ field }) => (
                    <FormItem>
                      <FormLabel>Maximum Genes per Cell</FormLabel>
                      <FormControl>
                        <Input
                          type="number"
                          placeholder="e.g., 5000"
                          {...field}
                          onChange={(e) => field.onChange(e.target.value ? parseInt(e.target.value) : undefined)}
                        />
                      </FormControl>
                      <FormDescription>
                        Remove cells with too many detected genes (doublets)
                      </FormDescription>
                      <FormMessage />
                    </FormItem>
                  )}
                />

                <FormField
                  control={form.control}
                  name="training_config.qc_params.min_counts_per_cell"
                  render={({ field }) => (
                    <FormItem>
                      <FormLabel>Minimum UMI per Cell</FormLabel>
                      <FormControl>
                        <Input
                          type="number"
                          placeholder="e.g., 500"
                          {...field}
                          onChange={(e) => field.onChange(e.target.value ? parseInt(e.target.value) : undefined)}
                        />
                      </FormControl>
                      <FormDescription>
                        Remove cells with low total UMI counts
                      </FormDescription>
                      <FormMessage />
                    </FormItem>
                  )}
                />

                <FormField
                  control={form.control}
                  name="training_config.qc_params.max_counts_per_cell"
                  render={({ field }) => (
                    <FormItem>
                      <FormLabel>Maximum UMI per Cell</FormLabel>
                      <FormControl>
                        <Input
                          type="number"
                          placeholder="e.g., 30000"
                          {...field}
                          onChange={(e) => field.onChange(e.target.value ? parseInt(e.target.value) : undefined)}
                        />
                      </FormControl>
                      <FormDescription>
                        Remove cells with extremely high UMI counts
                      </FormDescription>
                      <FormMessage />
                    </FormItem>
                  )}
                />
              </div>
            </div>

            {/* Mitochondrial & Ribosomal Filters */}
            <div className="space-y-4">
              <div className="flex items-center space-x-2">
                <h4 className="font-medium">Gene Type Filters</h4>
                <Badge variant="outline">Percentage Cutoffs</Badge>
              </div>

              <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                <FormField
                  control={form.control}
                  name="training_config.qc_params.max_pct_mt"
                  render={({ field }) => (
                    <FormItem>
                      <FormLabel>Max Mitochondrial %</FormLabel>
                      <FormControl>
                        <Input
                          type="number"
                          placeholder="e.g., 20"
                          step="0.1"
                          {...field}
                          onChange={(e) => field.onChange(e.target.value ? parseFloat(e.target.value) : undefined)}
                        />
                      </FormControl>
                      <FormDescription>
                        Remove cells with high mitochondrial gene expression
                      </FormDescription>
                      <FormMessage />
                    </FormItem>
                  )}
                />

                <FormField
                  control={form.control}
                  name="training_config.qc_params.max_pct_ribo"
                  render={({ field }) => (
                    <FormItem>
                      <FormLabel>Max Ribosomal %</FormLabel>
                      <FormControl>
                        <Input
                          type="number"
                          placeholder="e.g., 50"
                          step="0.1"
                          {...field}
                          onChange={(e) => field.onChange(e.target.value ? parseFloat(e.target.value) : undefined)}
                        />
                      </FormControl>
                      <FormDescription>
                        Remove cells with excessive ribosomal gene expression
                      </FormDescription>
                      <FormMessage />
                    </FormItem>
                  )}
                />
              </div>
            </div>
          </>
        )}
      </CardContent>
    </Card>
  );
};
