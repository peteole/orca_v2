#!/usr/bin/env python3
"""Debug script to check ORCA binary location during CI"""

import sys
from pathlib import Path

def debug_binary_locations():
    print("=== ORCA Binary Debug Info ===")
    
    # Current working directory
    print(f"Current working directory: {Path.cwd()}")
    
    # Module location
    module_dir = Path(__file__).resolve().parent / "orca"
    print(f"Orca module directory: {module_dir}")
    print(f"Orca module exists: {module_dir.exists()}")
    
    if module_dir.exists():
        print(f"Contents of orca directory:")
        for item in module_dir.iterdir():
            print(f"  - {item.name}")
    
    # Build directory locations
    build_locations = [
        Path.cwd() / "build",
        Path.cwd() / "build" / "bin",
        Path.cwd() / "build" / "Release",
        Path.cwd() / "build" / "Debug"
    ]
    
    print("\nBuild directory locations:")
    for loc in build_locations:
        print(f"  {loc}: exists={loc.exists()}")
        if loc.exists():
            try:
                contents = list(loc.iterdir())
                for item in contents:
                    print(f"    - {item.name}")
            except:
                print(f"    (could not list contents)")
    
    # Look for orca binaries specifically
    binary_names = ["orca", "orca.exe"]
    print(f"\nSearching for binaries: {binary_names}")
    
    for root in [Path.cwd(), module_dir.parent if module_dir.exists() else Path.cwd()]:
        for binary_name in binary_names:
            for path in root.rglob(binary_name):
                print(f"  Found: {path}")

if __name__ == "__main__":
    debug_binary_locations()
    
    # Try to import and run the function
    try:
        sys.path.insert(0, str(Path.cwd()))
        from orca.lib import get_binary_path
        print(f"\nget_binary_path() result: {get_binary_path()}")
    except Exception as e:
        print(f"\nError calling get_binary_path(): {e}")
