import pandas as pd
import numpy as np
import argparse
import networkx as nx

from tqdm import tqdm
from collections import defaultdict, deque

def is_leaf(T, node):
    return len(T[node]) == 0

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
        "labels", help="Labels for the leaves"
    )

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()

    tree = nx.read_edgelist(args.tree, create_using=nx.DiGraph())
    labels_csv = pd.read_csv(args.labels, sep=",").set_index("leaf")

    roots = [node for node in tree.nodes if len(list(tree.predecessors(node))) == 0]
    root = roots[0]

    def dist_f(x, y):
        if x is None or y is None:
            return 0
        if x == y:
            return 0
        return 1

    def leaf_f(node):
        if node not in labels_csv.index:
            return None
        y = labels_csv.loc[node, "label"]
        if y == "None":
            return None
        return y
 
    character_set = set(labels_csv["label"].unique())
    scores = mp(tree, root, character_set, leaf_f, dist_f)

    # backtrack to obtain ALL optimal solutions
    solutions = []
    minimum_cost = min(scores[root].values())
    for char in character_set:
        if scores[root][char] == minimum_cost:
            solutions.append({root: char})

    for node in tqdm(list(nx.dfs_preorder_nodes(tree, root))):
        if node == root: continue
        parent = list(tree.predecessors(node))[0]
        new_solutions = [] 
        for solution in solutions:
            parent_char = solution[parent]
            min_cost = min([scores[node][char] + dist_f(parent_char, char) for char in character_set])
            for char in character_set:
               if scores[node][char] + dist_f(parent_char, char) == min_cost:
                   new_solution = solution.copy()
                   new_solution[node] = char
                   new_solutions.append(new_solution)
        assert len(new_solutions) >= len(solutions) 
        assert len(solutions) > 0
        solutions = new_solutions
        print(len(solutions))

    print(len(solutions))
    
