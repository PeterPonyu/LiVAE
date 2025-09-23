
#!/usr/bin/env python3
"""
Smart startup script for LiVAE with intelligent port management
Handles port conflicts automatically
"""

import subprocess
import sys
import os
import time
import signal
import socket
import re
from contextlib import closing

class SmartPortManager:
    """Intelligent port management for LiVAE services"""
    
    def __init__(self):
        self.processes = []
        self.backend_port = 8000
        self.frontend_port = 3000
        self.actual_backend_port = None
        self.actual_frontend_port = None
        
    def check_port_available(self, port):
        """Check if a port is available"""
        with closing(socket.socket(socket.AF_INET, socket.SOCK_STREAM)) as sock:
            try:
                sock.bind(('127.0.0.1', port))
                return True
            except OSError:
                return False
    
    def find_available_port(self, start_port, max_attempts=10):
        """Find an available port starting from start_port"""
        for port in range(start_port, start_port + max_attempts):
            if self.check_port_available(port):
                return port
        return None
    
    def kill_process_on_port(self, port):
        """Kill process running on specific port"""
        try:
            if sys.platform == "win32":
                # Windows
                result = subprocess.run(
                    f'netstat -ano | findstr :{port}',
                    shell=True, capture_output=True, text=True
                )
                if result.stdout:
                    pid = result.stdout.split()[-1]
                    subprocess.run(f'taskkill /F /PID {pid}', shell=True)
                    return True
            else:
                # Unix/Linux/macOS
                result = subprocess.run(
                    f'lsof -ti:{port}',
                    shell=True, capture_output=True, text=True
                )
                if result.stdout:
                    pid = result.stdout.strip()
                    subprocess.run(f'kill -9 {pid}', shell=True)
                    return True
            return False
        except Exception:
            return False
    
    def resolve_port_conflict(self, port, service_name):
        """Handle port conflicts interactively"""
        print(f"⚠️  Port {port} is already in use for {service_name}")
        
        # Find alternative port
        alt_port = self.find_available_port(port + 1)
        
        print(f"🔍 Found alternative port: {alt_port}")
        print(f"Choose an option:")
        print(f"1) Use alternative port {alt_port}")
        print(f"2) Kill process on port {port} and use it")
        print(f"3) Enter custom port")
        print(f"4) Skip starting {service_name}")
        
        while True:
            try:
                choice = input("Enter choice (1-4): ").strip()
                
                if choice == "1":
                    return alt_port
                elif choice == "2":
                    if self.kill_process_on_port(port):
                        print(f"✅ Killed process on port {port}")
                        time.sleep(1)
                        if self.check_port_available(port):
                            return port
                        else:
                            print("❌ Port still not available, using alternative")
                            return alt_port
                    else:
                        print("❌ Failed to kill process, using alternative")
                        return alt_port
                elif choice == "3":
                    custom_port = int(input("Enter custom port: "))
                    if self.check_port_available(custom_port):
                        return custom_port
                    else:
                        print(f"❌ Port {custom_port} is also in use")
                        continue
                elif choice == "4":
                    return None
                else:
                    print("❌ Invalid choice, please enter 1-4")
                    continue
                    
            except (ValueError, KeyboardInterrupt):
                print("❌ Invalid input or cancelled")
                return alt_port
    
    def start_backend(self):
        """Start FastAPI backend with smart port management"""
        print("🚀 Starting FastAPI backend...")
        
        # Check if api/main.py exists
        if not os.path.exists("api/main.py"):
            print("❌ api/main.py not found!")
            print("💡 Make sure you're in the LiVAE folder")
            return None
        
        # Ensure api/__init__.py exists
        if not os.path.exists("api/__init__.py"):
            print("📝 Creating api/__init__.py...")
            open("api/__init__.py", 'a').close()
        
        # Check port availability
        if not self.check_port_available(self.backend_port):
            resolved_port = self.resolve_port_conflict(self.backend_port, "backend")
            if resolved_port is None:
                print("⏭️  Skipping backend startup")
                return None
            self.backend_port = resolved_port
        
        # Start uvicorn
        cmd = [
            sys.executable, "-m", "uvicorn",
            "api.main:app",
            "--host", "127.0.0.1",
            "--port", str(self.backend_port),
            "--reload"
        ]
        
        try:
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                universal_newlines=True
            )
            
            # Wait a bit and check if it started successfully
            time.sleep(2)
            if process.poll() is None:
                self.actual_backend_port = self.backend_port
                print(f"✅ Backend: http://127.0.0.1:{self.backend_port}")
                print(f"📖 API Docs: http://127.0.0.1:{self.backend_port}/docs")
                return process
            else:
                print("❌ Backend failed to start")
                return None
                
        except Exception as e:
            print(f"❌ Failed to start backend: {e}")
            return None
    
    def parse_serve_output(self, output):
        """Parse npx serve output to extract actual port"""
        patterns = [
            r"Local:\s+http://localhost:(\d+)",
            r"Local:\s+http://127\.0\.0\.1:(\d+)",
            r"port\s+(\d+)",
        ]
        
        for pattern in patterns:
            match = re.search(pattern, output)
            if match:
                return int(match.group(1))
        return None
    
    def start_frontend(self):
        """Start frontend with smart port management"""
        print("🌐 Starting frontend server...")
        
        # Check if frontend/out exists
        if not os.path.exists("frontend/out"):
            print("❌ frontend/out directory not found!")
            return None
        
        if not os.path.exists("frontend/out/index.html"):
            print("❌ frontend/out/index.html not found!")
            return None
        
        # Try npx serve first (it handles port conflicts automatically)
        server_options = [
            {
                "cmd": ["npx", "serve", "frontend/out", "-p", str(self.frontend_port)],
                "handles_conflicts": True,
                "name": "npx serve"
            },
            {
                "cmd": [sys.executable, "-m", "http.server", str(self.frontend_port)],
                "handles_conflicts": False,
                "name": "Python http.server",
                "cwd": "frontend/out"
            }
        ]
        
        for option in server_options:
            try:
                # Check port for servers that don't handle conflicts
                if not option["handles_conflicts"]:
                    if not self.check_port_available(self.frontend_port):
                        resolved_port = self.resolve_port_conflict(self.frontend_port, "frontend")
                        if resolved_port is None:
                            continue
                        # Update command with new port
                        option["cmd"][-1] = str(resolved_port)
                        self.frontend_port = resolved_port
                
                print(f"🔄 Trying {option['name']}...")
                
                # Start process
                kwargs = {}
                if "cwd" in option:
                    kwargs["cwd"] = option["cwd"]
                
                process = subprocess.Popen(
                    option["cmd"],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    universal_newlines=True,
                    **kwargs
                )
                
                # Wait and check if it started
                time.sleep(3)
                
                if process.poll() is None:
                    # For npx serve, try to get actual port from output
                    if option["handles_conflicts"]:
                        try:
                            # Read some output to see if port changed
                            output_lines = []
                            for _ in range(10):  # Read up to 10 lines
                                try:
                                    line = process.stdout.readline()
                                    if line:
                                        output_lines.append(line)
                                        print(f"📋 Server: {line.strip()}")
                                        # Try to extract port from output
                                        actual_port = self.parse_serve_output(line)
                                        if actual_port:
                                            self.actual_frontend_port = actual_port
                                            break
                                    else:
                                        break
                                except:
                                    break
                        except:
                            pass
                    
                    # Set actual port if not detected
                    if self.actual_frontend_port is None:
                        self.actual_frontend_port = self.frontend_port
                    
                    print(f"✅ Frontend: http://127.0.0.1:{self.actual_frontend_port}")
                    return process
                else:
                    print(f"❌ {option['name']} failed to start")
                    
            except FileNotFoundError:
                print(f"⚠️  {option['name']} not available")
                continue
            except Exception as e:
                print(f"❌ Error with {option['name']}: {e}")
                continue
        
        print("❌ Failed to start frontend with any server")
        return None
    
    def start_all_services(self):
        """Start all services with smart port management"""
        print("=" * 60)
        print("🎯 LiVAE Application Startup (Smart Port Mode)")
        print("=" * 60)
        
        # Verify directory structure
        if not os.path.exists("api") or not os.path.exists("frontend"):
            print("❌ Invalid directory structure")
            print("💡 Make sure you're in the LiVAE folder")
            print("💡 Expected structure:")
            print("   LiVAE/")
            print("   ├── api/main.py")
            print("   └── frontend/out/")
            return
        
        # Start backend
        print(f"\n🔍 Checking backend port {self.backend_port}...")
        backend_process = self.start_backend()
        if backend_process:
            self.processes.append(backend_process)
        
        time.sleep(1)
        
        # Start frontend
        print(f"\n🔍 Checking frontend port {self.frontend_port}...")
        frontend_process = self.start_frontend()
        if frontend_process:
            self.processes.append(frontend_process)
        
        if not self.processes:
            print("❌ No services started successfully")
            sys.exit(1)
        
        # Display final status
        self.display_status()
    
    def display_status(self):
        """Display final service status"""
        print("\n" + "=" * 60)
        print("🎉 LiVAE is running!")
        print("=" * 60)
        
        if self.actual_frontend_port:
            print(f"🌐 Frontend: http://127.0.0.1:{self.actual_frontend_port}")
            print(f"📱 Network: http://192.168.x.x:{self.actual_frontend_port}")
        
        if self.actual_backend_port:
            print(f"🔧 Backend: http://127.0.0.1:{self.actual_backend_port}")
            print(f"📖 API Docs: http://127.0.0.1:{self.actual_backend_port}/docs")
        
        print("\n💡 Press Ctrl+C to stop all services")
        print("=" * 60)
    
    def monitor_processes(self):
        """Monitor running processes"""
        try:
            while self.processes:
                time.sleep(1)
                # Remove dead processes
                self.processes = [p for p in self.processes if p.poll() is None]
        except KeyboardInterrupt:
            print("\n🛑 Received interrupt signal...")
            self.stop_all_services()
    
    def stop_all_services(self):
        """Stop all running services"""
        print("🔄 Stopping all services...")
        
        for process in self.processes:
            try:
                process.terminate()
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
        
        print("✅ All services stopped")
        sys.exit(0)

def main():
    """Main entry point"""
    manager = SmartPortManager()
    
    # Setup signal handlers
    signal.signal(signal.SIGINT, lambda s, f: manager.stop_all_services())
    signal.signal(signal.SIGTERM, lambda s, f: manager.stop_all_services())
    
    manager.start_all_services()
    manager.monitor_processes()

if __name__ == "__main__":
    main()
