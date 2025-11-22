#!/usr/bin/env python
"""Test script for FastAPI with frontend integration"""

import requests
import time

BASE_URL = "http://127.0.0.1:8000"

def test_api():
    """Test API endpoints"""
    print("Testing FastAPI with Static Frontend Integration\n")
    print("=" * 60)
    
    # Test 1: Root API endpoint
    print("\n1. Testing root API endpoint (/):")
    try:
        response = requests.get(f"{BASE_URL}/")
        print(f"   Status: {response.status_code}")
        print(f"   Response: {response.json()}")
    except Exception as e:
        print(f"   ✗ Error: {e}")
    
    # Test 2: Status endpoint
    print("\n2. Testing status endpoint (/status):")
    try:
        response = requests.get(f"{BASE_URL}/status")
        print(f"   Status: {response.status_code}")
        print(f"   Response: {response.json()}")
    except Exception as e:
        print(f"   ✗ Error: {e}")
    
    # Test 3: API documentation
    print("\n3. Testing API documentation (/docs):")
    try:
        response = requests.get(f"{BASE_URL}/docs")
        print(f"   Status: {response.status_code}")
        print(f"   Available at: {BASE_URL}/docs")
    except Exception as e:
        print(f"   ✗ Error: {e}")
    
    # Test 4: Root HTML (index.html from static build)
    print("\n4. Testing static index.html (frontend):")
    try:
        response = requests.get(f"{BASE_URL}/")
        if response.status_code == 200 and "<!DOCTYPE" in response.text:
            print(f"   Status: {response.status_code}")
            print("   ✓ Successfully serving static index.html")
            print(f"   Content length: {len(response.text)} bytes")
        else:
            print(f"   Status: {response.status_code}")
            print(f"   Response type: {response.headers.get('content-type')}")
    except Exception as e:
        print(f"   ✗ Error: {e}")
    
    # Test 5: Next.js static assets
    print("\n5. Testing Next.js static assets (/_next):")
    try:
        response = requests.head(f"{BASE_URL}/_next")
        print(f"   Status: {response.status_code}")
        if response.status_code == 200:
            print("   ✓ Static assets directory accessible")
    except Exception:
        print("   Status: Could not reach /_next directory")
    
    # Test 6: Frontend routes (client-side routing fallback)
    print("\n6. Testing frontend route handling (/upload):")
    try:
        response = requests.get(f"{BASE_URL}/upload")
        if response.status_code == 200 and "<!DOCTYPE" in response.text:
            print(f"   Status: {response.status_code}")
            print("   ✓ Client-side routing working (fallback to index.html)")
        else:
            print(f"   Status: {response.status_code}")
    except Exception as e:
        print(f"   ✗ Error: {e}")
    
    print("\n" + "=" * 60)
    print("\n✓ API tests completed!")
    print("\nFrontend available at: http://127.0.0.1:8000")
    print("API Documentation: http://127.0.0.1:8000/docs\n")

if __name__ == "__main__":
    time.sleep(1)  # Wait for server to be ready
    test_api()
