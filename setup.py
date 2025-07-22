#!/usr/bin/env python3
import os
import platform
import subprocess
import sys
from pathlib import Path

def build_orca_binary():
    """Build the ORCA C++ binary"""
    print("Building ORCA C++ binary...")
    
    build_dir = Path("build")
    build_dir.mkdir(exist_ok=True)
    
    try:
        # Build the C++ binary
        subprocess.run([
            "cmake", "-S", ".", "-B", "build", 
            "-DCMAKE_BUILD_TYPE=Release"
        ], check=True)
        
        subprocess.run([
            "cmake", "--build", "build", "--config", "Release"
        ], check=True)
        
        # Copy binary to orca package
        orca_dir = Path("orca")
        orca_dir.mkdir(exist_ok=True)
        
        binary_name = "orca"
        if platform.system() == "Windows":
            binary_name += ".exe"
        
        # Try multiple possible locations for the binary
        possible_paths = [
            build_dir / "bin" / binary_name,
            build_dir / binary_name,
            build_dir / "Release" / binary_name,  # Windows MSVC
            build_dir / "Debug" / binary_name     # Windows MSVC Debug
        ]
        
        src_binary = None
        for path in possible_paths:
            if path.exists():
                src_binary = path
                break
        
        if src_binary:
            import shutil
            dest_binary = orca_dir / binary_name
            shutil.copy2(src_binary, dest_binary)
            # Make binary executable on Unix systems
            if platform.system() != "Windows":
                os.chmod(dest_binary, 0o755)
            print(f"Binary copied to: {dest_binary}")
        else:
            print(f"Warning: Could not find binary in any of: {possible_paths}")
            
    except subprocess.CalledProcessError as e:
        print(f"Error building binary: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    # Build the binary when setup.py is run
    if len(sys.argv) > 1 and sys.argv[1] in ['build', 'install', 'bdist_wheel', 'sdist']:
        build_orca_binary()
    
    # Import setuptools after building to avoid early import issues
    from setuptools import setup
    setup()
