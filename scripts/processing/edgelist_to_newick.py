import argparse
import sys
import networkx as nx

def to_newick(tree, root):
    def dfs(node):
        if len(list(tree.successors(node))) == 0:
            return f'{node}'
        else:
            return f'({",".join([dfs(child) for child in tree.successors(node)])}){node}'
    return dfs(root) + ';'

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Converts an edgelist to a Newick file"
    )

    parser.add_argument(
        "tree"
    )

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    tree = nx.read_edgelist(args.tree, create_using=nx.DiGraph())
    root = list(nx.topological_sort(tree))[0]
    print(to_newick(tree, root))
