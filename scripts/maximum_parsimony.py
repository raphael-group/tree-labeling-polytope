import pandas as pd
import numpy as np
import argparse
import networkx as nx

from collections import defaultdict, deque

def is_leaf(T, node):
    return len(T[node]) == 0

"""
Runs maximum parsimony for a single character.

Inputs:
  - root: root node
  - character_set: set of possible characters
"""
def mp(T, root, character_set, leaf_f, dist_f):
    scores = defaultdict(dict) # map from (node, char) -> score i.e. d(T_v, x)
    stack = deque([root]) # stack to simulate DFS search

    while stack:
        node = stack.pop()

        if is_leaf(T, node):
            if leaf_f(node) is None:
                for char in character_set:
                    scores[node][char] = 0
                continue

            for char in character_set:
                if char == leaf_f(node):
                    scores[node][char] = 0
                else:
                    scores[node][char] = np.inf

            continue

        # all children are scored
        if all(child in scores for child in T[node]):
            for char in character_set:
                cost = 0
                for child in T[node]:
                    costs = [dist_f(char, child_char) + scores[child][child_char] for child_char in character_set]
                    cost += min(costs)
                scores[node][char] = cost

            continue

        stack.append(node)
        for child in T[node]:
            stack.append(child)

    return scores

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Performs maximum parsimony."
    )

    parser.add_argument(
        "tree", help="Tree in edgelist format"
    )

    parser.add_argument(
        "clones", help="Clones TSV file"
    )

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()

    tree = nx.read_edgelist(args.tree, create_using=nx.DiGraph(), nodetype=str)
    clone_tsv = pd.read_csv(args.clones, sep="\t", header=None, names=["cell_id", "clone_id"]).set_index("cell_id")

    def dist_f(x, y):
        if x is None or y is None:
            return 0

        if x == y:
            return 0

        return 1

    character_set = clone_tsv[clone_tsv["clone_id"] != "None"]["clone_id"].unique()

    def leaf_f(node):
        if node not in clone_tsv.index:
            return None

        y = clone_tsv.loc[node, "clone_id"]
        if y == "None":
            return None

        return y

    root = "A"
    res = mp(tree, root, character_set, leaf_f, dist_f)[root]
    print(min(res.values()))
