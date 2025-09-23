'use client';

import { useState } from 'react';
import { Button } from '@/components/ui/button';
import { Card } from '@/components/ui/card';
import { api } from '@/lib/api/endpoints';
import type { AppStatus } from '@/types/index';

export default function Home() {
  const [status, setStatus] = useState<AppStatus | null>(null);
  const [loading, setLoading] = useState(false);

  const testConnection = async () => {
    setLoading(true);
    try {
      const result = await api.getStatus();
      setStatus(result);
    } catch (error) {
      console.error('Connection failed:', error);
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="container mx-auto p-8">
      <Card className="p-6">
        <h1 className="text-2xl font-bold mb-4">Single Cell Deep Learning Agent</h1>
        
        <Button onClick={testConnection} disabled={loading}>
          {loading ? 'Testing...' : 'Test API Connection'}
        </Button>
        
        {status && (
          <div className="mt-4">
            <p>✅ Connection successful!</p>
            <pre>{JSON.stringify(status, null, 2)}</pre>
          </div>
        )}
      </Card>
    </div>
  );
}