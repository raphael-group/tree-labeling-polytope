import argparse
import itertools

import pandas as pd
import scipy as sp
import numpy as np
import networkx as nx
import cvxpy as cp
 
from loguru import logger
from tqdm import tqdm
from enum import Enum
from collections import defaultdict

def is_leaf(T, node):
    return len(T[node]) == 0

def select_grid_points():
    grid_points = []
    for y in np.linspace(0.00001, 0.2, 40):
        grid_points.append((1, y))

    for y in np.linspace(0.2, 1, 10):
        grid_points.append((1, y))

    for x in np.linspace(1, 10, 5):
        for y in np.linspace(0.00001, 1, 5):
            grid_points.append((x, y))

    return grid_points

"""
Creates a piecewise linear approximation of the function f(x, y) = x log(x/y)
for some collection of points, where each point is of the form (x, y).

Returns a dictionary mapping each point to a tuple (a, b) where a and b are
the coefficients of the linear function that underestimates f(x, y) at that 
point.
"""
def create_piecewise_linear_approximation(points):
    approximation = {}
    for x, y in points:
        b = (1 + np.log(x/y), -x/y)
        a = x*np.log(x/y) - (b[0] * x + b[1] * y)
        approximation[(x,y)] = (a, b)
    return approximation

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
    parser = argparse.ArgumentParser(description="Infers the probability matrix")
    parser.add_argument("tree", help="Tree in edgelist format")
    parser.add_argument("labels", help="Leaf labeling as a CSV file")
    parser.add_argument("-o", "--output", help="Output prefix", default="result")
    parser.add_argument("-r", "--root", help="Root of the tree", default=None)
    parser.add_argument("-m", "--mode", help="Mode of optimization", default="piecewise", choices=["piecewise", "parsimony", "exact"])
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

    ## Create Tree Labeling System ##
    A, b, A_column_idx = create_tree_labeling_system(tree, root, character_set, leaf_f)

    x = cp.Variable(A.shape[1], boolean=True)
    p = cp.Variable((len(character_set), len(character_set)), pos=True)

    character_set_idx = {}
    for i, c in enumerate(character_set):
        character_set_idx[c] = i

    constraints = [
        A @ x == b, 
        p @ np.ones(len(character_set)) == np.ones(len(character_set))
    ]

    grid_points = select_grid_points()
    approximations = create_piecewise_linear_approximation(grid_points)

    # for _ in range(100):
        # x, y = 1, np.random.uniform()
        # f_xy = 0 if x == 0 else x*np.log(x/y) 
        # f_approx = 0
        # for a, b in approximations.values():
        #   f_approx = max(f_approx, a + b[0] * x + b[1] * y)
        # print(f"({x}, {y}): {f_xy} -> {f_approx}")

    if args.mode == "piecewise":
        t = cp.Variable(A.shape[1], pos=True)
        objective_expr = cp.sum(t)
        for c, c_prime in itertools.product(character_set, character_set):
            for u, v in tree.edges:
                for approx in approximations.values():
                    a, b_vec = approx
                    constraint = t[A_column_idx[(u, v, c, c_prime)]] >= a + b_vec[0] * x[A_column_idx[(u, v, c, c_prime)]] + b_vec[1] * p[character_set_idx[c], character_set_idx[c_prime]]
                    constraints.append(constraint)
    elif args.mode == "parsimony":
        cost_vec = np.ones(A.shape[1])
        for u, v in tree.edges:
            for c in character_set:
                cost_vec[A_column_idx[(u, v, c, c)]] = 0
        objective_expr = cost_vec @ x
    elif args.mode == "scipy":
        # veritcally stack A and np.ones(A.shape[1]) to create a matrix
        A_matrix = [A]
    else:
        objective_expr = 0
        for c, c_prime in itertools.product(character_set, character_set):
            for u, v in tree.edges:
                objective_expr += cp.rel_entr(x[A_column_idx[(u, v, c, c_prime)]], p[character_set_idx[c], character_set_idx[c_prime]])

    objective = cp.Minimize(objective_expr)
    prob = cp.Problem(objective, constraints)
    import gurobipy as gp
    env = gp.Env()
    env.setParam('MIPGap', 0.80)
    prob.solve(verbose=True, solver='GUROBI',env=env)

    print(character_set_idx)
    print(p.value)
    print(x.value)
    for u, v in tree.edges:
        for c, c_prime in itertools.product(character_set, character_set):
            if x[A_column_idx[(u, v, c, c_prime)]].value > 0.5 and c_prime != c:
                print(f"({u}, {v}): {c} -> {c_prime}")
