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
from gurobipy import GRB

def is_leaf(T, node):
    return len(T[node]) == 0

"""
Creates the tree labeling polytope whose vertices correspond
to the set of solutions to the (unconstrained) maximum parsimony 
problem.
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
Constrain the tree labeling polytope to only allow certain migration graphs
by explicitly adding constraints to the model.
"""
def constrain_migration_graphs(model, character_set, constraint_type):
    logger.info("Adding migration graph constraints.")
    model.migration_graph_constraints = pyo.ConstraintList()
    if constraint_type == "polyclonal_tree":
        for S in itertools.chain.from_iterable(itertools.combinations(character_set, r) for r in range(2, len(character_set))):
            model.migration_graph_constraints.add(
                sum(model.migrations[(i, j)] for i in S for j in S if i != j) <= len(S) - 1
            )
    elif constraint_type == "polyclonal_dag":
        for S in itertools.chain.from_iterable(itertools.combinations(character_set, r) for r in range(2, len(character_set))):
            for C in itertools.permutations(S):
                model.migration_graph_constraints.add(
                    sum(model.migrations[(C[i], C[(i+1) % len(C)])] for i in range(len(C))) <= len(C) - 1
                )
    else:
        raise ValueError(f"Unknown constraint type: {constraint_type}")

"""
Constrain the tree labeling polytope to only allow certain migration graphs
by adding constraints to the model using a Gurobi callback.
"""
def constrain_migration_graphs_gurobi_callback(model, character_set, constraint_type):
    def dag_callback(model, gb_model, where):
        if where != GRB.Callback.MIPSOL:
            return

        # load solutions
        gb_model.cbGetSolution([model.migrations[i, j] for i in character_set for j in character_set if i != j])

        G = nx.DiGraph()
        for i in character_set:
            G.add_node(i)

        for i in character_set:
            for j in character_set:
                if i == j: continue
                if model.migrations[i, j].value > 0.5:
                    G.add_edge(i, j)

        try:
            S = nx.find_cycle(G)
            cycle_nodes = [i for (i,_) in S]
            for C in itertools.permutations(cycle_nodes):
                logger.info(f"Adding constraint {C} for cycle {S}.")

                cons = model.migration_graph_constraints.add(
                    sum(model.migrations[(C[i], C[(i+1) % len(C)])] for i in range(len(C))) <= len(C) - 1
                )
                gb_model.cbLazy(cons)
        except nx.exception.NetworkXNoCycle:
            pass

    def tree_callback(model, gb_model, where):
        if where != GRB.Callback.MIPSOL:
            return

        # load solutions
        gb_model.cbGetSolution([model.migrations[i, j] for i in character_set for j in character_set if i != j])

        G = nx.Graph()
        for i in character_set:
            G.add_node(i)

        for i in character_set:
            for j in character_set:
                if i == j: continue
                if model.migrations[i, j].value > 0.5:
                    G.add_edge(i, j)

        try:
            S = nx.find_cycle(G)
            cycle_nodes = [i for (i,_) in S]
            logger.info(f"Adding constraint for cycle {S}.")
            cons = model.migration_graph_constraints.add(
                sum(model.migrations[(i, j)] for i in cycle_nodes for j in cycle_nodes if i != j) <= len(cycle_nodes) - 1
            )
            gb_model.cbLazy(cons)
        except nx.exception.NetworkXNoCycle:
            pass

    if constraint_type == "polyclonal_tree" or constraint_type == "monoclonal_tree":
        return tree_callback
    elif constraint_type == "polyclonal_dag" or constraint_type == "monoclonal_dag":
        return dag_callback
    else:
        raise ValueError(f"Unknown constraint type: {constraint_type}")

    return

"""
Solves a generalization of the MACHINA parsimonious migration 
history problem using the tree labeling polytope.
"""
def fast_machina(tree, character_set, leaf_f, dist_f, args):
    solver = pyo.SolverFactory('gurobi_persistent')

    model = create_tree_labeling_polytope(tree, args.root, character_set, leaf_f, dist_f, root_label=args.label)
    if args.constraints != "none":
        append_migrations(model, tree, character_set)

    if args.constraints.startswith("monoclonal"):
        for c1, c2 in model.migrations:
            model.migration_constraints.add(
                sum(model.decisions[u, v, c1, c2] for u, v in tree.edges) <= model.migrations[c1, c2]
            )

    solver.set_instance(model)

    if args.constraints != "none":
        constraint_callback = constrain_migration_graphs_gurobi_callback(model, character_set, args.constraints)
        model.migration_graph_constraints = pyo.ConstraintList()
        solver.set_gurobi_param("LazyConstraints", 1)
        solver.set_callback(constraint_callback)

    solver.solve(model, tee=True)

    # compute (an) optimal vertex labeling
    vertex_labeling = {}
    for u, v in tree.edges:
        for c1, c2 in itertools.product(character_set, character_set):
            if np.abs(model.decisions[u, v, c1, c2]()) > 1e-4:
                vertex_labeling[u] = c1
                vertex_labeling[v] = c2

    return vertex_labeling

def exact_tnet(tree, labels_tsv, args):
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
    fast_machina_parser.add_argument("labels", help="Leaf labeling as a CSV file")
    fast_machina_parser.add_argument("-c", "--constraints", help="Migration graph constraints",
                                choices=["polyclonal_tree", "polyclonal_dag", "monoclonal_tree", "monoclonal_dag", "none"],
                                default="none")
    fast_machina_parser.add_argument("-o", "--output", help="Output prefix", default="result")
    fast_machina_parser.add_argument("-r", "--root", help="Root of the tree", default="root")
    fast_machina_parser.add_argument("-l", "--label", help="Root label", default=None)

    # exactTNET subparser
    exact_tnet_parser = subparsers.add_parser("fast_tnet", help="exactTNET")
    exact_tnet_parser.add_argument("tree", help="Tree in edgelist format")
    exact_tnet_parser.add_argument("labels", help="Leaf labeling as a CSV file")
    exact_tnet_parser.add_argument("-o", "--output", help="Output prefix", default="result")
    exact_tnet_parser.add_argument("-r", "--root", help="Root of the tree", default="root")
    exact_tnet_parser.add_argument("-l", "--label", help="Root label", default=None)

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_arguments()

    tree = nx.read_edgelist(args.tree, create_using=nx.DiGraph())
    labels_csv = pd.read_csv(args.labels, sep=",").set_index("leaf")

    if not nx.is_directed_acyclic_graph(tree):
        raise ValueError("Graph is not a tree, it contains cycles.")

    if not nx.is_weakly_connected(tree):
        raise ValueError("Graph is not connected, it is a forest.")

    # defines the distance function between characters x and y along an edge e
    def dist_f(e, x, y):
        if x is None or y is None:
            return 0

        if x == y:
            return 0

        return 1

    # defines the leaf labeling function
    def leaf_f(node):
        if node not in labels_csv.index:
            return None

        y = labels_csv.loc[node, "label"]
        if y == "None":
            return None

        return y
 
    # load the character set as the set of all leaf labels
    character_set = labels_csv[labels_csv["label"] != "None"]["label"].unique()
    if not all([leaf_f(node) is not None for node in tree.nodes if is_leaf(tree, node)]):
        unlabeled_leaves = [node for node in tree.nodes if is_leaf(tree, node) and leaf_f(node) is None]
        raise ValueError(f"Leaves {unlabeled_leaves} are unlabeled.")

    # computes the vertex labeling using the specified method
    if args.method == "fast_machina":
        vertex_labeling = fast_machina(tree, character_set, leaf_f, dist_f, args)
    elif args.method == "exact_tnet":
        vertex_labeling = exact_tnet()

    # writes an optimal labeling to a file 
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


