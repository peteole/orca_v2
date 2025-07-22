import numpy as np
from orca import run_orca
if __name__ == "__main__":
    # Example usage
    edges = np.array([[0, 1], [1, 2], [2, 0], [0, 3], [3, 4]])
    print(run_orca(edges, mode="node", graphlet_size=4))
    print(run_orca(edges, mode="edge", graphlet_size=5))