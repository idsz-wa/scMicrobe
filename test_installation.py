#!/usr/bin/env python3
"""
Test script to verify scMicrobe installation.
"""

import sys
import os
import subprocess

def test_python_modules():
    """Test Python module imports."""
    print("Testing Python modules...")
    
    sys.path.insert(0, 'python')
    
    modules = [
        'bam_parser',
        'host_filter',
        'microbe_aligner',
        'quantifier',
        'contamination_filter',
        'output_writer',
    ]
    
    failed = []
    for module in modules:
        try:
            __import__(module)
            print(f"  ✓ {module}")
        except Exception as e:
            print(f"  ✗ {module}: {e}")
            failed.append(module)
    
    return len(failed) == 0

def test_go_binary():
    """Test Go binary."""
    print("\nTesting Go binary...")
    
    go_binary = 'go/scmicro-go'
    if os.name == 'nt':
        go_binary = 'go/scmicro-go.exe'
    
    if not os.path.exists(go_binary):
        print(f"  ✗ Go binary not found: {go_binary}")
        return False
    
    print(f"  ✓ Go binary exists: {go_binary}")
    
    # Try to run it
    try:
        result = subprocess.run([go_binary, '-h'], 
                              capture_output=True, 
                              timeout=5)
        print(f"  ✓ Go binary executable")
        return True
    except Exception as e:
        print(f"  ✗ Go binary execution failed: {e}")
        return False

def test_external_tools():
    """Test external tool availability."""
    print("\nTesting external tools...")
    
    tools = ['minimap2', 'samtools']
    optional_tools = ['kraken2']
    
    all_found = True
    
    for tool in tools:
        try:
            result = subprocess.run([tool, '--version'], 
                                  capture_output=True, 
                                  timeout=5)
            print(f"  ✓ {tool}")
        except FileNotFoundError:
            print(f"  ✗ {tool} (not found in PATH)")
            all_found = False
        except Exception as e:
            print(f"  ? {tool} (could not verify: {e})")
    
    for tool in optional_tools:
        try:
            result = subprocess.run([tool, '--version'], 
                                  capture_output=True, 
                                  timeout=5)
            print(f"  ✓ {tool} (optional)")
        except FileNotFoundError:
            print(f"  - {tool} (optional, not found)")
    
    return all_found

def main():
    print("=" * 60)
    print("scMicrobe Installation Test")
    print("=" * 60)
    
    results = []
    
    # Test Python modules
    results.append(("Python modules", test_python_modules()))
    
    # Test Go binary
    results.append(("Go binary", test_go_binary()))
    
    # Test external tools
    results.append(("External tools", test_external_tools()))
    
    # Summary
    print("\n" + "=" * 60)
    print("Test Summary")
    print("=" * 60)
    
    for name, passed in results:
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"  {status}: {name}")
    
    all_passed = all(r[1] for r in results)
    
    if all_passed:
        print("\n✓ All tests passed!")
        return 0
    else:
        print("\n✗ Some tests failed. Please check the installation.")
        return 1

if __name__ == '__main__':
    sys.exit(main())
