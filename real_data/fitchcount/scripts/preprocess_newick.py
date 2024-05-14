import argparse
import sys
import networkx as nx

from Bio import Phylo

def from_newick_get_nx_tree(tree_path):
    phylo_tree = Phylo.read(tree_path, 'newick', rooted=True)
    net_tree = Phylo.to_networkx(phylo_tree)

    node_renaming_mapping = {}
    idx = 0
    for node in net_tree.nodes:
        if net_tree.out_degree(node) > 0:
            node_renaming_mapping[node] = f'clade_{idx}'
            idx = idx + 1
        else:
            node_renaming_mapping[node] = node.name
    node_renaming_mapping[list(net_tree.nodes)[0]] = 'root'

    net_tree = nx.relabel_nodes(net_tree, node_renaming_mapping)

    directed_tree = nx.DiGraph()
    directed_tree.add_edges_from(list(nx.bfs_edges(net_tree, 'root')))
    return directed_tree

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Converts a Newick file to an edgelist"
    )

    parser.add_argument(
        "tree"
    )


    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    tree = from_newick_get_nx_tree(args.tree)
    for edge in tree.edges:
        print(f'{edge[0]}\t{edge[1]}')
