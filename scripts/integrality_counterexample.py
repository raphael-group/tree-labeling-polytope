import argparse
import random
import itertools

import pandas as pd
import numpy as np
import networkx as nx
import pyomo.environ as pyo
 
from loguru import logger
from tqdm import tqdm
from enum import Enum
from dataclasses import dataclass
from collections import defaultdict, deque

def small_parsimony_lp_edge_flow(T, root, character_set, leaf_f, dist_f, rand_b=False):
    if rand_b:
        # rand int between 0 and 10
        b_1 = random.randint(0, 2)
        b_2 = random.randint(0, 2)
        b_3 = random.randint(0, 2)
    else:
        b_1 = 0
        b_2 = 1
        b_3 = 0
    
    # add a dummy root node
    T.add_node("dummy_root")
    T.add_edge("dummy_root", root)

    # one hot encode the characters for each edge
    # i.e. x_{u,v,c,c'} = 1 if u is assigned character c and v is assigned character c'
    model = pyo.ConcreteModel()
    # set domain to be [0, 1]
    model.x = pyo.Var(T.edges, character_set, character_set, domain=pyo.NonNegativeReals, initialize=0)

    # require \sum_{c'} x_{u,v,c,c'} = \sum_{c'}x_{v,w,c',c} for all u,v,w, c
    model.edge_constraints = pyo.ConstraintList()
    for u, v in T.edges:
        for c in character_set:
            for w in T[v]:
                model.edge_constraints.add(
                    sum(model.x[u, v, c2, c] for c2 in character_set) - sum(model.x[v, w, c, c2] for c2 in character_set) == b_1
                )

    # set \sum_{c} x_{u,v,c,c'} = 1 for all e=(u,v), v is a leaf, and c' is the label of v
    # set \sum_{c} x_{u,v,c,c'} = 0 for all e=(u,v), v is a leaf, and c' is not the label of v
    model.leaf_constraints = pyo.ConstraintList()
    for u, v in T.edges:
        if not is_leaf(T, v): continue
        for c in character_set:
            if c == leaf_f(v):
                model.leaf_constraints.add(sum(model.x[u, v, c2, c] for c2 in character_set) == b_2)
            else:
                model.leaf_constraints.add(sum(model.x[u, v, c2, c] for c2 in character_set) == b_3)

    # set objectiv eto be \sum_{uv} \sum_{c,c'} x_{u,v,c,c'} * d(c, c')
    model.objective = pyo.Objective(expr=sum(model.x[u, v, c, c2] * dist_f(c, c2) for u, v in T.edges for c in character_set for c2 in character_set))
    return model
        
def small_parsimony_ilp(T, root, character_set, leaf_f, dist_f):
    # one hot encode the characters for each node,
    # i.e. x_{v, c} = 1 if v is assigned character c
    model = pyo.ConcreteModel()
    # set domain to be [0, 1]
    model.x = pyo.Var(T.nodes, character_set, domain=pyo.NonNegativeReals)

    # require \sum_{c} x_{v,c} = 1 for all v
    model.node_constraints = pyo.ConstraintList()
    for node in T.nodes:
        if is_leaf(T, node): 
            model.node_constraints.add(sum(model.x[node, c] for c in character_set) == 1)
            continue
        model.node_constraints.add(sum(model.x[node, c] for c in character_set) == 1)

    # set x_{v,c} = 1 if v is a leaf and c is the label of v
    model.leaf_constraints = pyo.ConstraintList()
    for node in T.nodes:
        if not is_leaf(T, node): continue
        model.leaf_constraints.add(model.x[node, leaf_f(node)] == 1)

    # set the objective to be \sum_{uv} ||x_u - x_v||_1
    # by definition, ||x_u - x_v||_1 = \sum_{c} |x_{u,c} - x_{v,c}| = \sum_{c} \max(x_{u,c}, x_{v,c}) - \min(x_{u,c}, x_{v,c})
    # so we can add another variable y_{u,v,c} = \max(x_{u,c}, x_{v,c}) and z_{u,v,c} = \min(x_{u,c}, x_{v,c})
    # and then the objective is \sum_{uv} \sum_{c} y_{u,v,c} - z_{u,v,c}
    model.y = pyo.Var(T.edges, character_set, domain=pyo.NonNegativeReals)
    model.z = pyo.Var(T.edges, character_set, domain=pyo.NonNegativeReals)
    model.l1_constraints = pyo.ConstraintList()
    for u, v in T.edges:
        for c in character_set:
            model.l1_constraints.add(model.y[u, v, c] >= model.x[u, c])
            model.l1_constraints.add(model.y[u, v, c] >= model.x[v, c])
            model.l1_constraints.add(model.z[u, v, c] <= model.x[u, c])
            model.l1_constraints.add(model.z[u, v, c] <= model.x[v, c])

    model.objective = pyo.Objective(expr=sum(model.y[u, v, c] - model.z[u, v, c] for u, v in T.edges for c in character_set))
    return model

def dual_mp(T, root, character_set, leaf_f, dist_f):
    cs = defaultdict(dict) # map from (node, char) -> score i.e. d(T_v, x)
    stack = deque([root]) # stack to simulate DFS search

    while stack:
        node = stack.pop()

        if is_leaf(T, node):
            if leaf_f(node) is None:
                assert False

            i = leaf_f(node)
            for j in character_set:
                cs[node][j] = dist_f(j, i)

            continue

        # all children are scored
        if all(child in cs for child in T[node]):
            for j in character_set:
                cs[node][j] = min(dist_f(j, i) + sum(cs[child][i] for child in T[node]) for i in character_set)
            continue

        stack.append(node)
        for child in T[node]:
            stack.append(child)

    return cs

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

def is_leaf(T, node):
    return len(T[node]) == 0

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Tests the small parsimony ILP on a random directed tree and leaf labeling"
    )

    parser.add_argument(
        "--n", type=int, default=50, help="Number of nodes in the tree"
    )

    parser.add_argument(
        "--k", type=int, default=10, help="Number of characters"
    )

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()

    while True:
        # construct a random, rooted tree
        undirected_T = nx.random_tree(args.n)
        root = random.choice(list(undirected_T.nodes))
        print(root)

        T = nx.DiGraph() 
        for node in undirected_T.nodes:
            T.add_node(node)

        visited = set()
        stack = [root]
        visited.add(root)
        
        while stack:
            node = stack.pop()
            visited.add(node)
            for child in undirected_T[node]:
                if child in visited:
                    continue
                T.add_edge(node, child)
                stack.append(child)

        # write edgelist to file as TSV
        with open("tree.tsv", "w") as f:
            for u, v in T.edges:
                f.write(f"{u}\t{v}\n")

        # assign random character labels to the leaves
        character_set = list(range(args.k))
        leaf_labels = {}
        for node in T.nodes:
            if is_leaf(T, node):
                leaf_labels[node] = random.choice(character_set)

        # write leaf labeling to file as TSV
        with open("leaf_labels.tsv", "w") as f:
            for node, label in leaf_labels.items():
                f.write(f"{node}\t{label}\n")

        # compute the scores
        leaf_f = lambda node: leaf_labels[node]
        dist_f = lambda x, y: 0 if x == y else 1

        # create a random distance function
        dist_matrix = np.zeros((len(character_set), len(character_set)))
        for i in range(len(character_set)):
            for j in range(len(character_set)):
                if i == j: continue
                # rand float
                dist_matrix[i, j] = 3 * random.random()
                # rand int 
                # dist_matrix[i, j] = random.randint(1, 3)

        dist_f = lambda x, y: dist_matrix[x, y]

        scores = mp(T, root, character_set, leaf_f, dist_f)
        exact_score = min(scores[root].values())

        # solve the ILP
        model = small_parsimony_ilp(T, root, character_set, leaf_f, dist_f)
        solver = pyo.SolverFactory("gurobi")
        solver.solve(model)
        value1 = model.objective()

        # solve the LP relaxation
        model = small_parsimony_lp_edge_flow(T, root, character_set, leaf_f, dist_f)
        solver = pyo.SolverFactory("gurobi")
        solver.solve(model)
        value2 = model.objective()

        # solve dual MP
        cs = dual_mp(T, root, character_set, leaf_f, dist_f)
        dual_obj = min(cs[root].values())

        # grab random node
        node = random.choice(list(T.nodes))
        # print(cs[node])
        # print(scores[node])
        # print(T.edges)
        # print(leaf_labels)
        # print(dist_matrix)

        logger.info(f"Exact score: {exact_score}, LP #1 score: {0.5*value1}, LP #2 score: {value2}, dual MP score: {dual_obj}")
        if np.abs(exact_score - value2) > 1e-4 or np.abs(exact_score - dual_obj) > 1e-4:
            logger.info(f"Tree: {T.edges}")
            logger.info(f"Leaf labels: {leaf_labels}")
            logger.info(f"Scores: {scores}")
            logger.info(f"ILP solution: {model.x.pprint()}")
            break





