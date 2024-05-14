import argparse
import sys
import networkx as nx
import pandas as pd
import numpy as np

from Bio import Phylo

def from_newick_get_nx_tree(tree_path):
    phylo_tree = Phylo.read(tree_path, 'newick', rooted=True)
    net_tree = Phylo.to_networkx(phylo_tree)
    return net_tree

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

    parser.add_argument(
        "-f", "--format",
        default="newick",
        choices=["newick", "edgelist"],
        help="The format of the input tree file"
    )

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()


    if args.format == "newick":
        tree = from_newick_get_nx_tree(args.tree)
    else:
        tree = nx.read_edgelist(args.tree, create_using=nx.DiGraph(), data=(("weight", float),))

    labeling = pd.read_csv(args.labeling)
    labels = labeling['label'].unique().tolist()
    m = len(labels)

    print('parent,child,label1,label2,probability')
    for edge in tree.edges:
        branch_length = tree[edge[0]][edge[1]]['weight']
        if m == 1:
            for l1 in labels:
                for l2 in labels:
                    if l1 == l2:
                        prob = 1
                    else:
                        prob = 0
                    print(f'{edge[0]},{edge[1]},{l1},{l2},{prob}')
            continue

        base = (m-1)/m * np.exp(-((m/(m-1))*branch_length))
        for l1 in labels:
            for l2 in labels:
                if l1 == l2:
                    prob = 1/m + base
                else:
                    prob = 1/m - 1/(m-1) * base

                print(f'{edge[0]},{edge[1]},{l1},{l2},{prob}')
