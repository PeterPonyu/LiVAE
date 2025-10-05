// ========================================
// src/components/monitor/constants.ts
// ========================================
export const METRICS_CONFIG = {
  loss: { color: '#dc2626', label: 'Loss', betterDirection: 'down' },
  ari: { color: '#16a34a', label: 'ARI', betterDirection: 'up' },
  nmi: { color: '#2563eb', label: 'NMI', betterDirection: 'up' },
  asw: { color: '#7c3aed', label: 'ASW', betterDirection: 'up' },
  ch: { color: '#ea580c', label: 'CH', betterDirection: 'up' },
  db: { color: '#dc2626', label: 'DB', betterDirection: 'down' },
  pc: { color: '#059669', label: 'PC', betterDirection: 'up' },
} as const;

export type MetricKey = keyof typeof METRICS_CONFIG;