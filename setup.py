#!/usr/bin/env python3
import os
import platform
import subprocess
import sys
from pathlib import Path
from setuptools import setup, Extension
from setuptools.command.build_py import build_py

# Import bdist_wheel conditionally since it might not be available
try:
    from setuptools.command.bdist_wheel import bdist_wheel
    HAVE_WHEEL = True
except ImportError:
    HAVE_WHEEL = False


class BuildBinaryCommand(build_py):
    """Custom build command that builds the C++ binary before building Python package."""
    
    def run(self):
        self.build_orca_binary()
        super().run()
    
    def build_orca_binary(self):
        """Build the ORCA C++ binary"""
        print("Building ORCA C++ binary...")

        build_dir = Path("build")
        build_dir.mkdir(exist_ok=True)

        try:
            # Build the C++ binary
            subprocess.run(
                ["cmake", "-S", ".", "-B", "build", "-DCMAKE_BUILD_TYPE=Release"],
                check=True,
            )

            subprocess.run(["cmake", "--build", "build", "--config", "Release"], check=True)

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
                build_dir / "Debug" / binary_name,  # Windows MSVC Debug
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


if HAVE_WHEEL:
    class CustomBdistWheel(bdist_wheel):
        """Custom wheel build command that ensures platform-specific wheel."""
        
        def finalize_options(self):
            super().finalize_options()
            # Force platform-specific wheel since we include a binary
            self.root_is_pure = False
            # Ensure this is treated as a platlib (platform library) package
            self.plat_name_supplied = True
        
        def get_tag(self):
            # Force platform-specific tag
            python, abi, plat = super().get_tag()
            if plat == "any":
                try:
                    from distutils.util import get_platform
                    plat = get_platform().replace('-', '_').replace('.', '_')
                except ImportError:
                    import sysconfig
                    plat = sysconfig.get_platform().replace('-', '_').replace('.', '_')
            return python, abi, plat
else:
    # Dummy class when wheel is not available
    class CustomBdistWheel:
        pass


if __name__ == "__main__":
    # Create a dummy extension to force platform-specific wheel
    # This ensures the wheel is treated as platlib instead of purelib
    dummy_ext = Extension(
        name='orca._dummy',
        sources=['orca/_dummy.c'],
        optional=True
    )
    
    # Prepare cmdclass based on what's available
    cmdclass = {
        'build_py': BuildBinaryCommand,
    }
    if HAVE_WHEEL:
        cmdclass['bdist_wheel'] = CustomBdistWheel
    
    setup(
        ext_modules=[dummy_ext],
        cmdclass=cmdclass,
    )
