import argparse
import sys
import networkx as nx

from Bio import Phylo

def from_newick_get_nx_tree(tree_path):
    phylo_tree = Phylo.read(tree_path, 'newick')
    net_tree = Phylo.to_networkx(phylo_tree)

    node_renaming_mapping = {}
    idx = 0
    for node in net_tree.nodes:
        if nx.degree(net_tree, node) == 1:
            continue

        node_renaming_mapping[node] = f'clade_{idx}'
        idx = idx + 1

    node_renaming_mapping[list(net_tree.nodes)[0]] = 'root'

    net_tree = nx.relabel_nodes(net_tree, node_renaming_mapping)

    directed_tree = nx.DiGraph()
    directed_tree.add_edges_from(list(nx.bfs_edges(net_tree, 'root')))
    for edge in directed_tree.edges:
        u, v = edge
        directed_tree[u][v]['weight'] = net_tree[u][v]['weight']

    return directed_tree

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Converts a Newick file to an edgelist"
    )

    parser.add_argument(
        "-b", "--branch-lengths",
        action="store_true",
        help="Include branch lengths in the output",
        default=False
    )

    parser.add_argument(
        "tree"
    )

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    tree = from_newick_get_nx_tree(args.tree)
    for edge in tree.edges:
        if args.branch_lengths:
            print(f'{edge[0]}\t{edge[1]}\t{tree[edge[0]][edge[1]]["weight"]}')
        else:
            print(f'{edge[0]}\t{edge[1]}')
