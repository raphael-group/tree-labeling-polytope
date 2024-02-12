import argparse
import itertools

import pandas as pd
import numpy as np
import networkx as nx
import pyomo.environ as pyo
 
from loguru import logger
from tqdm import tqdm
from enum import Enum
from collections import defaultdict

def is_leaf(T, node):
    return len(T[node]) == 0

"""
Creates the tree labeling polytope whose vertices correspond
to the set of solutions to the maximum parsimony problem.
"""
def create_tree_labeling_polytope(T, root, character_set, leaf_f, dist_f, root_label=None):
    model = pyo.ConcreteModel()

    logger.info(f"Creating tree labeling polytope for tree with {len(T.nodes)} nodes and {len(T.edges)} edges.")

    # add a dummy root node
    T.add_node("dummy_root")
    T.add_edge("dummy_root", root)

    # i.e. x_{u,v,c,c'} = 1 if u is assigned character c and v is assigned character c'
    model = pyo.ConcreteModel()
    # set domain to be [0, 1]
    model.decisions = pyo.Var(T.edges, character_set, character_set, domain=pyo.NonNegativeReals, initialize=0)

    # require \sum_{c'} x_{u,v,c,c'} = \sum_{c'}x_{v,w,c',c} for all u,v,w, c
    model.edge_constraints = pyo.ConstraintList()
    for u, v in T.edges:
        for c in character_set:
            for w in T[v]:
                model.edge_constraints.add(
                    sum(model.decisions[u, v, c2, c] for c2 in character_set) - sum(model.decisions[v, w, c, c2] for c2 in character_set) == 0
                )

    # set \sum_{c} x_{u,v,c,c'} = 1 for all e=(u,v), v is a leaf, and c' is the label of v
    # set \sum_{c} x_{u,v,c,c'} = 0 for all e=(u,v), v is a leaf, and c' is not the label of v
    model.leaf_constraints = pyo.ConstraintList()
    for u, v in T.edges:
        if not is_leaf(T, v): continue
        for c in character_set:
            if c == leaf_f(v):
                model.leaf_constraints.add(sum(model.decisions[u, v, c2, c] for c2 in character_set) == 1)
            else:
                model.leaf_constraints.add(sum(model.decisions[u, v, c2, c] for c2 in character_set) == 0)

    # set root label if provided
    if root_label is not None:
        model.root_constraint = pyo.ConstraintList()
        model.root_constraint.add(sum(model.decisions["dummy_root", root, c, root_label] for c in character_set) == 1)

    # set objective to be \sum_{uv} \sum_{c,c'} x_{u,v,c,c'} * d(c, c')
    model.objective = pyo.Objective(expr=sum(model.decisions[u, v, c1, c2] * dist_f((u, v), c1, c2) for u, v in T.edges for c1 in character_set for c2 in character_set))
    return model

"""
Appends migration variables to the tree labeling polytope.
"""
def append_migrations(model, T, character_set):
    logger.info("Adding migration constraints.")
    edges = [(i, j) for i in character_set for j in character_set if i != j]
    model.migrations = pyo.Var(edges, domain=pyo.Binary, initialize=1)
    model.migration_constraints = pyo.ConstraintList()
    for u, v, c1, c2 in model.decisions:
        if c1 == c2: continue
        model.migration_constraints.add(
            model.decisions[u, v, c1, c2] <= model.migrations[c1, c2]
        )

"""
Constrain the tree labeling polytope to only allow certain migration graphs.
"""
def constrain_migration_graphs(model, character_set, constraint_type):
    logger.info("Adding migration graph constraints.")
    model.migration_graph_constraints = pyo.ConstraintList()
    if constraint_type == "tree":
        for S in itertools.chain.from_iterable(itertools.combinations(character_set, r) for r in range(2, len(character_set))):
            model.migration_graph_constraints.add(
                sum(model.migrations[(i, j)] for i in S for j in S if i != j) <= len(S) - 1
            )
    elif constraint_type == "dag":
        for S in itertools.chain.from_iterable(itertools.combinations(character_set, r) for r in range(2, len(character_set))):
            for C in itertools.permutations(S):
                model.migration_graph_constraints.add(
                    sum(model.migrations[(C[i], C[(i+1) % len(C)])] for i in range(len(C))) <= len(C) - 1
                )
    else:
        raise ValueError(f"Unknown constraint type: {constraint_type}")

def fast_machina(tree, labels_tsv, args):
    def dist_f(e, x, y):
        if x is None or y is None:
            return 0

        if x == y:
            return 0

        return 1

    def leaf_f(node):
        if node not in labels_tsv.index:
            return None

        y = labels_tsv.loc[node, "label"]
        if y == "None":
            return None

        return y
 
    # load the character set as the set of all leaf labels
    character_set = labels_tsv[labels_tsv["label"] != "None"]["label"].unique()
    if not all([leaf_f(node) is not None for node in tree.nodes if is_leaf(tree, node)]):
        unlabeled_leaves = [node for node in tree.nodes if is_leaf(tree, node) and leaf_f(node) is None]
        raise ValueError(f"Leaves {unlabeled_leaves} are unlabeled.")

    model = create_tree_labeling_polytope(tree, args.root, character_set, leaf_f, dist_f, root_label=args.label)

    if args.constraints != "none":
        append_migrations(model, tree, character_set)
        constrain_migration_graphs(model, character_set, args.constraints)

    solver = pyo.SolverFactory('gurobi')
    solver.solve(model, tee=True)

    # compute (an) optimal vertex labeling
    vertex_labeling = {}
    for u, v in tree.edges:
        for c1, c2 in itertools.product(character_set, character_set):
            if np.abs(model.decisions[u, v, c1, c2]()) > 1e-4:
                vertex_labeling[u] = c1
                vertex_labeling[v] = c2

    # write (an) optimal labeling to a file 
    with open(f"{args.output}_vertex_labeling.csv", "w") as f:
        f.write("vertex,label\n")
        for node, label in vertex_labeling.items():
            f.write(f"{node},{label}\n")

    # compute the migration multi-graph
    migration_graph = defaultdict(int)
    for u, v in tree.edges:
        if vertex_labeling[u] != vertex_labeling[v]:
            migration_graph[(vertex_labeling[u], vertex_labeling[v])] += 1

    # write migration multi-graph to a file 
    with open(f"{args.output}_migration_graph.csv", "w") as f:
        f.write("source,target,count\n")
        for (i, j), count in migration_graph.items():
            f.write(f"{i},{j},{count}\n")

def fast_tnet(tree, labels_tsv, args):
    pass

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Constrained tree labeling using the tree labeling polytope."
    )

    subparsers = parser.add_subparsers(dest="method", help="Methods")
    subparsers.required = True

    # fastMACHINA subparser
    fast_machina_parser = subparsers.add_parser("fast_machina", help="fastMACHINA")
    fast_machina_parser.add_argument("tree", help="Tree in edgelist format")
    fast_machina_parser.add_argument("labels", help="Leaf labeling as a TSV file")
    fast_machina_parser.add_argument("-c", "--constraints", help="Migration graph constraints",
                                choices=["polyclonal_tree", "polyclonal_dag", "monoclonal_tree", "monoclonal_dag", "none"],
                                default="none")
    fast_machina_parser.add_argument("-o", "--output", help="Output prefix", default="result")
    fast_machina_parser.add_argument("-r", "--root", help="Root of the tree", default="root")
    fast_machina_parser.add_argument("-l", "--label", help="Root label", default=None)

    # fastTNET subparser
    fast_tnet_parser = subparsers.add_parser("fast_tnet", help="fastTNET")
    fast_tnet_parser.add_argument("tree", help="Tree in edgelist format")
    fast_tnet_parser.add_argument("labels", help="Leaf labeling as a TSV file")
    fast_tnet_parser.add_argument("-o", "--output", help="Output prefix", default="result")
    fast_tnet_parser.add_argument("-r", "--root", help="Root of the tree", default="root")
    fast_tnet_parser.add_argument("-l", "--label", help="Root label", default=None)

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_arguments()

    tree = nx.read_edgelist(args.tree, create_using=nx.DiGraph(), nodetype=int)
    labels_tsv = pd.read_csv(args.labels, sep="\t", header=None, names=["id", "label"]).set_index("id")

    if not nx.is_directed_acyclic_graph(tree):
        raise ValueError("Graph is not a tree, it contains cycles.")

    if not nx.is_weakly_connected(tree):
        raise ValueError("Graph is not connected, it is a forest.")

    print(args)

    if args.method == "fast_machina":
        fast_machina()
    elif args.method == "fast_tnet":
        fast_tnet()


