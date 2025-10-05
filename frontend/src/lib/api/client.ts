
// src/lib/api/client.ts
const API_BASE_URL = 'http://127.0.0.1:8000';

export interface ApiError extends Error {
  status?: number;
  statusText?: string;
  details?: Record<string, unknown> | string;
}

class ApiClient {
  private baseURL: string;

  constructor(baseURL: string = API_BASE_URL) {
    this.baseURL = baseURL;
  }

  getBaseURL(): string {
    return this.baseURL;
  }

  private async handleResponse<T>(response: Response): Promise<T> {
    if (!response.ok) {
      let errorDetails: Record<string, unknown> | string;
      
      try {
        errorDetails = await response.json();
      } catch {
        errorDetails = await response.text();
      }

      const apiError: ApiError = new Error(
        (typeof errorDetails === 'object' ? errorDetails.detail as string : errorDetails) ||
        `HTTP error! status: ${response.status}`
      );
      
      apiError.status = response.status;
      apiError.statusText = response.statusText;
      apiError.details = errorDetails;
      
      throw apiError;
    }
    
    return response.json();
  }

  async get<T>(endpoint: string): Promise<T> {
    const response = await fetch(`${this.baseURL}${endpoint}`, {
      headers: {
        'Accept': 'application/json',
      },
    });
    return this.handleResponse<T>(response);
  }

  async post<T>(endpoint: string, data?: Record<string, unknown>): Promise<T> {
    const response = await fetch(`${this.baseURL}${endpoint}`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Accept': 'application/json',
      },
      body: data ? JSON.stringify(data) : undefined,
    });
    
    return this.handleResponse<T>(response);
  }

  async uploadFile<T>(endpoint: string, file: File, onProgress?: (progress: number) => void): Promise<T> {
    return new Promise((resolve, reject) => {
      const formData = new FormData();
      formData.append('file', file);
      
      const xhr = new XMLHttpRequest();
      
      if (onProgress) {
        xhr.upload.onprogress = (event) => {
          if (event.lengthComputable) {
            const percentComplete = (event.loaded / event.total) * 100;
            onProgress(percentComplete);
          }
        };
      }
      
      xhr.onload = () => {
        if (xhr.status >= 200 && xhr.status < 300) {
          try {
            const result = JSON.parse(xhr.responseText);
            resolve(result);
          } catch {
            reject(new Error('Failed to parse response'));
          }
        } else {
          const apiError: ApiError = new Error(`Upload failed: ${xhr.statusText}`);
          apiError.status = xhr.status;
          apiError.statusText = xhr.statusText;
          reject(apiError);
        }
      };
      
      xhr.onerror = () => {
        reject(new Error('Network error during upload'));
      };
      
      xhr.open('POST', `${this.baseURL}${endpoint}`);
      xhr.send(formData);
    });
  }

  async downloadBlob(endpoint: string, onProgress?: (progress: number) => void): Promise<Blob> {
    return new Promise((resolve, reject) => {
      const xhr = new XMLHttpRequest();
      xhr.responseType = 'blob';
      
      if (onProgress) {
        xhr.onprogress = (event) => {
          if (event.lengthComputable) {
            const percentComplete = (event.loaded / event.total) * 100;
            onProgress(percentComplete);
          }
        };
      }
      
      xhr.onload = () => {
        if (xhr.status >= 200 && xhr.status < 300) {
          resolve(xhr.response);
        } else {
          const apiError: ApiError = new Error(`Download failed: ${xhr.statusText}`);
          apiError.status = xhr.status;
          apiError.statusText = xhr.statusText;
          reject(apiError);
        }
      };
      
      xhr.onerror = () => {
        reject(new Error('Network error during download'));
      };
      
      xhr.open('GET', `${this.baseURL}${endpoint}`);
      xhr.send();
    });
  }

  getDownloadUrl(endpoint: string): string {
    return `${this.baseURL}${endpoint}`;
  }

  async healthCheck(timeoutMs: number = 5000): Promise<boolean> {
    try {
      const controller = new AbortController();
      const timeoutId = setTimeout(() => controller.abort(), timeoutMs);
      
      const response = await fetch(`${this.baseURL}/status`, {
        signal: controller.signal,
        headers: { 'Accept': 'application/json' },
      });
      
      clearTimeout(timeoutId);
      return response.ok;
    } catch {
      return false;
    }
  }
}

export const apiClient = new ApiClient();
