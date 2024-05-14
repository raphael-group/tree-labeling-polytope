import argparse
import itertools

import pandas as pd
import numpy as np
import networkx as nx
import cvxpy as cp
 
from loguru import logger
from tqdm import tqdm
from enum import Enum
from collections import defaultdict

def is_leaf(T, node):
    return len(T[node]) == 0

"""
Creates the tree labeling polytope whose vertices correspond
to the set of solutions to the (unconstrained) maximum parsimony 
problem.
"""
def create_tree_labeling_system(T, root, character_set, leaf_f, root_label=None):
    logger.info(f"Creating tree labeling polytope for tree with {len(T.nodes)} nodes and {len(T.edges)} edges.")
    logger.info(f"Character set has size: {len(character_set)}")

    # add a dummy root node
    T.add_node("dummy_root")
    T.add_edge("dummy_root", root)

    A_rows    = []
    A_row_gen = lambda: np.zeros(len(T.edges) * len(character_set) * len(character_set))
    b         = []

    idx = 0
    A_column_idx = {}
    for u, v in T.edges:
        for c, c_prime in itertools.product(character_set, character_set):
            A_column_idx[(u, v, c, c_prime)] = idx
            idx += 1

    # require \sum_{c'} x_{u,v,c',c} = \sum_{c'}x_{v,w,c,c'} for all u,v,w,c
    for u, v in T.edges:
        for c in character_set:
            for w in T[v]:
                A_row = A_row_gen()
                for c2 in character_set:
                    A_row[A_column_idx[(u, v, c2, c)]] = 1
                    A_row[A_column_idx[(v, w, c, c2)]] = -1
                A_rows.append(A_row)
                b.append(0)

    # set \sum_{c} x_{u,v,c,c'} = 1 for all e=(u,v), v is a leaf, and c' is the label of v
    # set \sum_{c} x_{u,v,c,c'} = 0 for all e=(u,v), v is a leaf, and c' is not the label of v
    for u, v in T.edges:
        if not is_leaf(T, v): continue
        for c in character_set:
            A_row = A_row_gen()
            for c2 in character_set:   
                A_row[A_column_idx[(u, v, c2, c)]] = 1
            A_rows.append(A_row)

            if c == leaf_f(v):
                b.append(1)
            else:
                b.append(0)

    # set root label if provided
    if root_label is not None:
        A_row = A_row_gen()
        for c in character_set:
            A_row[A_column_idx[("dummy_root", root, c, root_label)]] = 1
        A_rows.append(A_row)
        b.append(1)

    A = np.array(A_rows)
    b = np.array(b)
    return A, b, A_column_idx

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Infers the probability matrix"
    )

    parser.add_argument("tree", help="Tree in edgelist format")
    parser.add_argument("labels", help="Leaf labeling as a CSV file")
    parser.add_argument("-o", "--output", help="Output prefix", default="result")
    parser.add_argument("-r", "--root", help="Root of the tree", default=None)

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_arguments()

    ## Data Preprocessing ##

    tree = nx.read_edgelist(args.tree, create_using=nx.DiGraph(), data=(("weight", float),))
    labels_csv = pd.read_csv(args.labels, sep=",").set_index("leaf")

    if not nx.is_directed_acyclic_graph(tree):
        raise ValueError("Graph is not a tree, it contains cycles.")

    if not nx.is_weakly_connected(tree):
        raise ValueError("Graph is not connected, it is a forest.")

    if args.root is not None:
        root = args.root
    else:
        roots = [node for node in tree.nodes if len(list(tree.predecessors(node))) == 0]
        if len(roots) != 1:
            raise ValueError(f"Tree has {len(roots)} roots, please specify the root.")
        root = roots[0]

    def leaf_f(node):
        if node not in labels_csv.index:
            return None

        y = labels_csv.loc[node, "label"]
        if y == "None":
            return None

        return y

    character_set = labels_csv[labels_csv["label"] != "None"]["label"].unique()
    if not all([leaf_f(node) is not None for node in tree.nodes if is_leaf(tree, node)]):
        unlabeled_leaves = [node for node in tree.nodes if is_leaf(tree, node) and leaf_f(node) is None]
        raise ValueError(f"Leaves {unlabeled_leaves} are unlabeled.")

    ## Create Tree Labeling Polytope ##
    A, b, A_column_idx = create_tree_labeling_system(tree, root, character_set, leaf_f)

    # x = cp.Variable(A.shape[1], pos=True)
    x = cp.Variable(A.shape[1], boolean=True)
    t = cp.Variable((len(character_set), len(character_set)), pos=True)
    p = cp.Variable((len(character_set), len(character_set)), pos=True)

    character_set_idx = {}
    for i, c in enumerate(character_set):
        character_set_idx[c] = i

    constraints = [
        A @ x == b, 
        p @ np.ones(len(character_set)) == np.ones(len(character_set))
    ]
    
    for c in character_set:
        for c_prime in character_set:
            constraints.append(
                t[character_set_idx[c], character_set_idx[c_prime]] == cp.sum([x[A_column_idx[(u, v, c, c_prime)]] for u, v in tree.edges])
            )

    objective_expr = 0
    for c, c_prime in itertools.product(character_set, character_set):
        objective_expr += -cp.rel_entr(
            t[character_set_idx[c], character_set_idx[c_prime]], 
            cp.sum([t[character_set_idx[c_prime_prime], character_set_idx[c_prime]] for c_prime_prime in character_set])
        )

    objective = cp.Minimize(objective_expr)
    prob = cp.Problem(objective, constraints)
    prob.solve(verbose=True)

    print(character_set_idx)
    print(p.value)
    print(x.value)
    print(t.value)
