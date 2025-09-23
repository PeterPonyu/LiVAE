
// start.ts
import { createServer, IncomingMessage, ServerResponse } from 'http';
import { parse } from 'url';
import next from 'next';
import * as path from 'path';
import * as fs from 'fs';
import express, { Request, Response } from 'express';

const dev = process.env.NODE_ENV !== 'production';
const hostname = process.env.HOSTNAME || 'localhost';
const port = parseInt(process.env.PORT || '3000', 10);

// Type definitions
interface ServerOptions {
  dev: boolean;
  hostname: string;
  port: number;
}

// Check if we have a static build
const outDir = path.join(__dirname, 'out');
const hasStaticBuild = fs.existsSync(outDir);

async function startStaticServer(): Promise<void> {
  try {
    const app = express();
    
    // Serve static files with proper headers
    app.use(express.static('out', {
      maxAge: '1y',
      etag: true,
      lastModified: true,
      setHeaders: (res: Response, filePath: string) => {
        // Set appropriate headers for different file types
        if (filePath.endsWith('.html')) {
          res.setHeader('Cache-Control', 'public, max-age=0, must-revalidate');
        } else if (filePath.endsWith('.js') || filePath.endsWith('.css')) {
          res.setHeader('Cache-Control', 'public, max-age=31536000, immutable');
        }
      }
    }));
    
    // Handle SPA routing - 修复类型问题
    app.get('*', (req: Request, res: Response) => {
      res.sendFile(path.join(__dirname, 'out', 'index.html'));
    });
    
    app.listen(port, hostname, () => {
      console.log(`🚀 Static server ready at http://${hostname}:${port}`);
      console.log(`📁 Serving files from: ${outDir}`);
    });
  } catch (error) {
    console.error('❌ Failed to start static server:', error);
    process.exit(1);
  }
}

async function startNextServer(): Promise<void> {
  try {
    const app = next({ dev, hostname, port } as ServerOptions);
    const handle = app.getRequestHandler();
    
    await app.prepare();
    
    createServer(async (req: IncomingMessage, res: ServerResponse) => {
      try {
        const parsedUrl = parse(req.url || '', true);
        await handle(req, res, parsedUrl);
      } catch (err) {
        console.error('❌ Error occurred handling', req.url, err);
        res.statusCode = 500;
        res.end('Internal server error');
      }
    }).listen(port, hostname, () => {
      console.log(`🚀 Next.js server ready at http://${hostname}:${port}`);
      console.log(`📊 Environment: ${dev ? 'development' : 'production'}`);
    });
  } catch (error) {
    console.error('❌ Failed to start Next.js server:', error);
    process.exit(1);
  }
}

// Main execution
if (hasStaticBuild && !dev) {
  console.log('📦 Detected static build, starting static server...');
  startStaticServer();
} else {
  console.log('⚛️  Starting Next.js server...');
  startNextServer();
}

// Graceful shutdown
process.on('SIGTERM', () => {
  console.log('👋 SIGTERM received, shutting down gracefully');
  process.exit(0);
});

process.on('SIGINT', () => {
  console.log('👋 SIGINT received, shutting down gracefully');
  process.exit(0);
});
