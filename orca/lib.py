import sys
from typing import Literal, Optional
import numpy as np
import tempfile
import subprocess
from pathlib import Path


def get_binary_path():
    """Get path to the ORCA binary bundled with the package."""
    # Look for binary in the same directory as this module
    binary_dir = Path(__file__).resolve().parent
    
    # Get if the current platform is Windows
    is_windows = Path(sys.executable).suffix == ".exe"
    binary_name = "orca.exe" if is_windows else "orca"
    
    binary_path = binary_dir / binary_name
    
    # If not found in package directory, try build directory (for development)
    if not binary_path.exists():
        build_dir = Path(__file__).resolve().parent.parent / "build" / "bin"
        binary_path = build_dir / binary_name
    
    if not binary_path.exists():
        raise FileNotFoundError(f"ORCA binary not found at {binary_path}")
    
    return binary_path


def run_orca(
    edge_list: np.ndarray,
    num_nodes: Optional[int] = None,
    mode: Literal["node", "edge"] = "node",
    graphlet_size: Literal[4, 5] = 4,
    debug: bool = False,
) -> np.ndarray:
    if num_nodes is None:
        num_nodes = int(np.max(edge_list)) + 1
    with tempfile.TemporaryDirectory() as _temp_dir:
        temp_dir = Path(_temp_dir)
        graph_path = temp_dir / "graph.txt"
        with open(graph_path, "w") as f:
            f.write(f"{num_nodes} {len(edge_list)}\n")
            for edge in edge_list:
                f.write(f"{edge[0]} {edge[1]}\n")
        results_path = temp_dir / "results.txt"
        # Assuming orca is a command line tool that can be run with subprocess
        output = subprocess.run(
            [get_binary_path(), mode, str(graphlet_size), str(graph_path), str(results_path)],
            capture_output=True,
            text=True,
        )
        if debug:
            print("ORCA command output:", output.stdout)
            print("ORCA command error:", output.stderr)
        if output.returncode != 0:
            raise RuntimeError(f"ORCA failed: {output.stderr}")

        with open(results_path, "r") as f:
            results_list = []
            for line in f:
                if line.strip():
                    results_list.append(list(map(int, line.split(" "))))

    return np.array(results_list)

def orca_nodes(
    edge_list: np.ndarray,
    num_nodes: Optional[int] = None,
    graphlet_size: Literal[4, 5] = 4,
    debug: bool = False,
) -> np.ndarray:
    return run_orca(edge_list, num_nodes, mode="node", graphlet_size=graphlet_size, debug=debug)

def orca_edges(
    edge_list: np.ndarray,
    num_nodes: Optional[int] = None,
    graphlet_size: Literal[4, 5] = 4,
    debug: bool = False,
) -> np.ndarray:
    return run_orca(edge_list, num_nodes, mode="edge", graphlet_size=graphlet_size, debug=debug)