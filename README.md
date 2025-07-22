# ORCA: Graphlet Counting

A Python wrapper for the ORCA (Orbit Counting Algorithm) tool for efficient graphlet counting in networks.

## Installation

```bash
pip install orca-graphlets
```

## Usage

```python
from orca import run_orca
import numpy as np

# Define edges as a 2D array
edges = np.array([[0, 1], [1, 2], [2, 0]])
num_nodes = 3

# Count 4-node graphlets
result = run_orca(edges, num_nodes=num_nodes, mode="node", graphlet_size=4)
print(result)
```

## Features

- Cross-platform support (Windows, macOS, Linux)
- Efficient C++ implementation
- Python interface with NumPy integration
- Support for node and edge-based graphlet counting
- Graphlet sizes 4 and 5

## License

MIT License
