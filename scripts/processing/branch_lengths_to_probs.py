import argparse
import sys
import networkx as nx
import pandas as pd
import numpy as np

from Bio import Phylo

def from_newick_get_nx_tree(tree_path):
    phylo_tree = Phylo.read(tree_path, 'newick')
    net_tree = Phylo.to_networkx(phylo_tree)

    node_renaming_mapping = {}
    idx = 0
    for node in net_tree.nodes:
        if str(node) == 'Clade':
            node_renaming_mapping[node] = f'clade_{idx}'
            idx = idx + 1
        else:
            node_renaming_mapping[node] = node.name
    node_renaming_mapping[list(net_tree.nodes)[0]] = 'root'

    net_tree = nx.relabel_nodes(net_tree, node_renaming_mapping)

    directed_tree = nx.DiGraph()
    for edge in net_tree.edges:
        directed_tree.add_edge(edge[0], edge[1], weight=net_tree[edge[0]][edge[1]]['weight'])
    return directed_tree

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Converts a Newick file to a mutation probability file"
    )

    parser.add_argument(
        "tree"
    )

    parser.add_argument(
        "labeling"
    )

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    tree = from_newick_get_nx_tree(args.tree)
    labeling = pd.read_csv(args.labeling)
    labels = labeling['label'].unique().tolist()
    m = len(labels)

    print('parent\tchild\tlabel1\tlabel2\tprobability')
    for edge in tree.edges:
        branch_length = tree[edge[0]][edge[1]]['weight']
        base = (m-1)/m * np.exp(-((m/(m-1))*branch_length))
        for l1 in labels:
            for l2 in labels:
                if l1 == l2:
                    prob = 1/m + base
                    print(f'{edge[0]}\t{edge[1]}\t{l1}\t{l2}\t{prob}')
                else:
                    prob = 1/m - 1/(m-1) * base
                    print(f'{edge[0]}\t{edge[1]}\t{l1}\t{l2}\t{prob}')
