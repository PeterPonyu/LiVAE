
// src/lib/validation/training-schemas.ts
import { z } from 'zod';

// Fix: Use z.enum() properly for species
export const qcParamsSchema = z.object({
  species: z.enum(['human', 'mouse']),  // Fixed: proper z.enum() usage
  min_genes_per_cell: z.number()
    .min(1, 'Must be at least 1')
    .max(10000, 'Must be less than 10,000')
    .optional(),
  max_genes_per_cell: z.number()
    .min(100, 'Must be at least 100')
    .max(50000, 'Must be less than 50,000')
    .optional(),
  min_counts_per_cell: z.number()
    .min(1, 'Must be at least 1')
    .max(100000, 'Must be less than 100,000')
    .optional(),
  max_counts_per_cell: z.number()
    .min(1000, 'Must be at least 1,000')
    .max(1000000, 'Must be less than 1,000,000')
    .optional(),
  max_pct_mt: z.number()
    .min(0, 'Must be at least 0%')
    .max(100, 'Must be less than 100%')
    .optional(),
  max_pct_ribo: z.number()
    .min(0, 'Must be at least 0%')
    .max(100, 'Must be less than 100%')
    .optional(),
});

// Agent Parameters Schema
export const agentParametersSchema = z.object({
  layer: z.string().min(1, 'Layer is required'),
  percent: z.number()
    .min(0.001, 'Must be at least 0.1%')
    .max(1.0, 'Must be at most 100%'),
  irecon: z.number()
    .min(0, 'Must be non-negative')
    .max(10, 'Must be less than 10'),
  lorentz: z.number()
    .min(0, 'Must be non-negative')
    .max(10, 'Must be less than 10'),
  beta: z.number()
    .min(0.1, 'Must be at least 0.1')
    .max(10, 'Must be less than 10'),
  dip: z.number()
    .min(0, 'Must be non-negative')
    .max(10, 'Must be less than 10'),
  tc: z.number()
    .min(0, 'Must be non-negative')
    .max(10, 'Must be less than 10'),
  info: z.number()
    .min(0, 'Must be non-negative')
    .max(10, 'Must be less than 10'),
  hidden_dim: z.number()
    .min(16, 'Must be at least 16')
    .max(1024, 'Must be less than 1024'),
  latent_dim: z.number()
    .min(2, 'Must be at least 2')
    .max(100, 'Must be less than 100'),
  i_dim: z.number()
    .min(2, 'Must be at least 2')
    .max(50, 'Must be less than 50'),
  lr: z.number()
    .min(1e-6, 'Must be at least 1e-6')
    .max(1e-1, 'Must be less than 0.1'),
});

// Training Configuration Schema
export const trainingConfigSchema = z.object({
  epochs: z.number()
    .min(1, 'Must be at least 1 epoch')
    .max(10000, 'Must be less than 10,000 epochs'),
  apply_qc: z.boolean(),
  qc_params: qcParamsSchema.optional(),
});

// Combined Form Schema
export const trainingFormSchema = z.object({
  training_config: trainingConfigSchema,
  agent_parameters: agentParametersSchema,
});

// Type inference from schemas
export type QCParamsFormData = z.infer<typeof qcParamsSchema>;
export type AgentParametersFormData = z.infer<typeof agentParametersSchema>;
export type TrainingConfigFormData = z.infer<typeof trainingConfigSchema>;
export type TrainingFormData = z.infer<typeof trainingFormSchema>;
