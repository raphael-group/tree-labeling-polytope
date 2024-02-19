import argparse
import networkx as nx

def parse_args():
    parser = argparse.ArgumentParser(description="Create weights from migration graph.")
    parser.add_argument("migration_graph", type=str, help="Path to migration graph.")
    return parser.parse_args()

def main():
    M = 10000000
    args = parse_args()
    migration_graph = nx.read_edgelist(args.migration_graph, create_using=nx.DiGraph)
    nodes = migration_graph.nodes()
    print("src,dst,weight")
    for u in nodes:
        for v in nodes:
            if u == "None" or v == "None":
                continue

            if migration_graph.has_edge(u, v):
                print(f"{u},{v},1")
            elif u == v:
                print(f"{u},{v},0")
            else:
                print(f"{u},{v},{M}")

if __name__ == "__main__":
    main()
