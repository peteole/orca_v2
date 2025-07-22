from orca import run_orca
from pathlib import Path
import numpy as np
def test_outputs_match():
    modes = ["node", "node", "edge", "edge"]
    graphlet_sizes = [4, 5, 4, 5]
    test_data_dir = Path(__file__).parent / "test_data"
    input_file = test_data_dir / "graph.in"
    with open(input_file, "r") as f:
        lines = f.readlines()
        num_nodes, num_edges = map(int, lines[0].strip().split())
        edges = np.array([list(map(int, line.strip().split())) for line in lines[1:]])
    assert edges.shape[0] == num_edges
    assert edges.shape[1] == 2
    for mode, graphlet_size in zip(modes, graphlet_sizes):
        expected_output_file = test_data_dir / f"{mode}{graphlet_size}.out"
        expected_output = np.loadtxt(expected_output_file, dtype=int)
        actual_output = run_orca(edges, num_nodes=num_nodes, mode=mode, graphlet_size=graphlet_size)
        np.testing.assert_array_equal(actual_output, expected_output)