
import sys
from orca.orca import Orca


def main():
    """Main function - command line interface"""
    if len(sys.argv) != 5:
        print("Incorrect number of arguments.")
        print(
            "Usage: python orca_pure_python.py [orbit type: node|edge] [graphlet size: 4/5] [graph - input file] [graphlets - output file]"
        )
        return

    orbit_type = sys.argv[1]
    if orbit_type not in ["node", "edge"]:
        print(f"Incorrect orbit type '{orbit_type}'. Should be 'node' or 'edge'.")
        return

    try:
        graphlet_size = int(sys.argv[2])
    except ValueError:
        print(f"Invalid graphlet size '{sys.argv[2]}'. Should be 4 or 5.")
        return

    if graphlet_size not in [4, 5]:
        print(f"Incorrect graphlet size {graphlet_size}. Should be 4 or 5.")
        return

    input_file = sys.argv[3]
    output_file = sys.argv[4]

    # Initialize ORCA
    orca = Orca()

    # Read graph
    if not orca.read_graph(input_file):
        print("Failed to read graph. Stopping!")
        return

    # Setup adjacency structures
    orca.setup_adjacency()

    # Run orbit counting
    if orbit_type == "node":
        print(f"Counting NODE orbits of graphlets on {graphlet_size} nodes.\n")
        if graphlet_size == 4:
            orca.count4()
        elif graphlet_size == 5:
            orca.count5()
        orca.write_results(output_file, graphlet_size)
    else:
        print(f"Counting EDGE orbits of graphlets on {graphlet_size} nodes.\n")
        if graphlet_size == 4:
            orca.ecount4()
        elif graphlet_size == 5:
            orca.ecount5()
        orca.write_edge_results(output_file, graphlet_size)

    print("Done!")

if __name__ == "__main__":
    main()