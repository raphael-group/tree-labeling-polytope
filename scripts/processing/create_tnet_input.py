import argparse
import sys
import pandas as pd
import networkx as nx

def arbitrarily_resolve_polytomies(T):
    T = T.copy()

    clade_idx = 1

    for u in list(T.nodes):
        children = list(T[u].keys())

        if len(children) <= 2:
            continue

        for v in children:
            T.remove_edge(u, v)

        us = []
        for i in range(len(children) - 2):
            new_u = f'clade_{clade_idx}'
            T.add_node(new_u)
            us.append(new_u)
            clade_idx += 1

            if i == 0:
                T.add_edge(u, new_u)
            else:
                T.add_edge(us[i - 1], new_u)

        us = [u] + us 
        assert len(us) + 1 == len(children)

        T.add_edge(us[-1], children[-1])
        for w, v in zip(us, children[:-1]):
            T.add_edge(w, v)

    return T

def to_newick(tree, root, leaf_labeling):
    def dfs(node):
        if len(list(tree.successors(node))) == 0:
            label = leaf_labeling.loc[node, 'label']
            return f'{label}_{node}'
        else:
            return f'({",".join([dfs(child) for child in tree.successors(node)])}){node}'
    return dfs(root) + ';'

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Creates the input for TNET"
    )

    parser.add_argument("tree", help="Tree file in edgelist format")
    parser.add_argument("leaf_labeling", help="Leaf labeling file in CSV format")

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    tree = nx.read_edgelist(args.tree, create_using=nx.DiGraph(), data=(("weight", float),))
    while True:
        has_rank_two_nodes = False
        nodes = list(tree)
        for n in nodes:
            if n not in tree.nodes: continue
            if tree.out_degree(n) != 1: continue
            if tree.in_degree(n) != 1:
                tree.remove_node(n)
                continue

            parent = list(tree.predecessors(n))[0]
            child = list(tree.successors(n))[0]
            has_rank_two_nodes = True
            tree.remove_node(n)
            tree.add_edge(parent, child)
            break

        if not has_rank_two_nodes:
            break
       
    tree = arbitrarily_resolve_polytomies(tree)
    leaf_labeling = pd.read_csv(args.leaf_labeling).set_index('leaf')
    root = list(nx.topological_sort(tree))[0]
    print(to_newick(tree, root, leaf_labeling))
