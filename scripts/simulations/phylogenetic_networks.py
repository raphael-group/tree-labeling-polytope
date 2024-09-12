import networkx as nx
import pandas as pd
import numpy as np
import argparse
import ngesh
import ete3

def parse_args():
    parser = argparse.ArgumentParser(description='Simulate phylogenetic network')
    parser.add_argument('-n', '--num_taxa', type=int, default=10, help='Number of taxa')
    parser.add_argument('-m', '--num_characters', type=int, default=10, help='Expected number of characters')
    parser.add_argument('-r', '--num_reticulations', type=int, default=20, help='Number of reticulations')
    parser.add_argument('-s', '--seed', type=int, default=0, help='Random seed')
    parser.add_argument('-o', '--output', type=str, default='sim_phylo_net', help='Output prefix')
    parser.add_argument('--uniform', action='store_true', help='Uniform distribution of characters', default=False)
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    T = ngesh.random_tree.gen_tree(min_leaves=args.num_taxa, num_leaves=args.num_taxa, seed=args.seed)

    idx = 0
    for node in T.traverse():
        node.name = str(idx)
        idx += 1
    
    nx_T = nx.DiGraph()
    for node in T.traverse():
        if node.up:
            nx_T.add_edge(node.up.name, node.name)

    tree_order = list(nx.topological_sort(nx_T))
    tree_order = {name: i for i, name in enumerate(tree_order)}

    np.random.seed(args.seed)
    num_reticulations_added = 0
    nodes = list(nx_T.nodes)
    while num_reticulations_added < args.num_reticulations:
        u = np.random.choice(nodes)
        v = np.random.choice(nodes)
        if u == v:
            continue

        if nx_T.has_edge(u, v):
            continue

        if nx_T.in_degree(v) == 2 or nx_T.out_degree(u) == 0:
            continue

        if tree_order[v] < tree_order[u]: # if v is an ancestor of u, swap
            continue

        nx_T.add_edge(u, v)
        num_reticulations_added += 1


    if not args.uniform:
        T = ngesh.add_characters(T, args.num_characters, 3.0, 1.0, seed=args.seed)
        leaf_labels = {node.name: set(node.chars) for node in T.traverse() if node.is_leaf()}
        # leaf labelings is a list of the mutations, convert it to a binary vector
        num_characters = max(len(label) for label in leaf_labels.values())
        for taxon, label in leaf_labels.items():
            leaf_labels[taxon] = [1 if i in label else 0 for i in range(num_characters)]
    else:
        np.random.seed(args.seed)
        ALPHABET = ['A', 'C', 'G', 'T']
        leaf_labels = {node.name: [np.random.choice(ALPHABET) for _ in range(args.num_characters)] for node in T.traverse() if node.is_leaf()}

    with open(args.output + '_labeling.csv', 'w') as f:
        f.write("taxon,sequence\n")
        for taxon, label in leaf_labels.items():
            f.write(f"{taxon},{''.join(map(str, label))}\n")
    
    # write nx_T to edgelist
    nx.write_edgelist(nx_T, args.output + '_network.edgelist', delimiter='\t', data=False)

