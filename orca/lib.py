from typing import Literal, Optional
import numpy as np
from orca.orca import run_orca

def orca_nodes(
    edge_list: np.ndarray,
    num_nodes: Optional[int] = None,
    graphlet_size: Literal[4, 5] = 4,
    debug: bool = False,
) -> np.ndarray:
    return run_orca(
        edge_list, num_nodes, mode="node", graphlet_size=graphlet_size, debug=debug
    )


def orca_edges(
    edge_list: np.ndarray,
    num_nodes: Optional[int] = None,
    graphlet_size: Literal[4, 5] = 4,
    debug: bool = False,
) -> np.ndarray:
    return run_orca(
        edge_list, num_nodes, mode="edge", graphlet_size=graphlet_size, debug=debug
    )
