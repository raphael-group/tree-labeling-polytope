import argparse
import sys
import networkx as nx

def to_newick(tree, root, branch_lengths=False):
    def dfs(node):
        if len(list(tree.successors(node))) == 0:
            if branch_lengths and node != root:
                parent = list(tree.predecessors(node))[0]
                return f'{node}:{tree[parent][node]["weight"]}'
            return f'{node}'
        else:
            if branch_lengths and node != root:
                parent = list(tree.predecessors(node))[0]
                return f'({",".join([dfs(child) for child in tree.successors(node)])}){node}:{tree[parent][node]["weight"]}'
            return f'({",".join([dfs(child) for child in tree.successors(node)])}){node}'
    return dfs(root) + ';'

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Converts an edgelist to a Newick file"
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
    if args.branch_lengths:
        tree = nx.read_edgelist(args.tree, create_using=nx.DiGraph(), data=[('weight', float)])
    else:
        tree = nx.read_edgelist(args.tree, create_using=nx.DiGraph())

    root = list(nx.topological_sort(tree))[0]
    print(to_newick(tree, root, args.branch_lengths))
