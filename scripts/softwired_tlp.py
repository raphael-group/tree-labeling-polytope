import argparse
import itertools
import json
import sys

import pandas as pd
import numpy as np
import networkx as nx
 
from loguru import logger
from tqdm import tqdm
from enum import Enum
from collections import defaultdict

import pyomo.environ as pyo
from gurobipy import GRB
import gurobipy as gp

def create_softwired_fischer(N, root, sequences, alphabet):
    num_characters = len(sequences.iloc[0]['sequence'])
    characters     = list(range(num_characters))
    characters_to_idx = {c: i for i, c in enumerate(characters)}

    def leaf_f(c, v):
        return sequences.loc[v]['sequence'][characters_to_idx[c]]

    reticulation_nodes = set([n for n in N.nodes if N.in_degree(n) > 1])
    reticulation_edges = set([(u, v) for u, v in N.edges if N.in_degree(v) > 1])

    model = pyo.ConcreteModel()
    model.labelings = pyo.Var(characters, N.nodes, alphabet, domain=pyo.Binary)
    model.mutations = pyo.Var(N.edges, domain=pyo.Binary)
    model.reticulations = pyo.Var(characters, N.edges, domain=pyo.Binary)

    model.labeling_constraints = pyo.ConstraintList()
    for c in characters:
        for v in N.nodes:
            model.labeling_constraints.add(sum(model.labelings[c, v, i] for i in alphabet) == 1)

            if N.out_degree(v) == 0:
                model.labeling_constraints.add(model.labelings[c, v, leaf_f(c, v)] == 1)

    model.reticulation_constraints = pyo.ConstraintList()
    for c in characters:
        for v in reticulation_nodes:
            model.reticulation_constraints.add(
                sum(model.reticulations[c, (u, v)] for u in N.predecessors(v)) == 1
            )

        for v in N.nodes:
            parents = list(N.predecessors(v))
            if len(parents) != 1: continue
            u = parents[0]
            model.reticulation_constraints.add(model.reticulations[c, (u, v)] == 1)

    model.mutation_constraints = pyo.ConstraintList()
    for c in characters:
        for u, v in N.edges:
            for s in alphabet:
                model.mutation_constraints.add(
                    model.mutations[(u, v)] >= model.labelings[c, u, s] - model.labelings[c, v, s] - (1 - model.reticulations[c, (u, v)])
                )

                model.mutation_constraints.add(
                    model.mutations[(u, v)] >= model.labelings[c, v, s] - model.labelings[c, u, s] - (1 - model.reticulations[c, (u, v)])
                )

    model.objective_expr = sum(model.mutations[(u, v)] for u, v in N.edges)
    model.objective = pyo.Objective(expr=model.objective_expr, sense=pyo.minimize)
    return model

def create_softwired_tlp(N, root, sequences, alphabet):
    num_characters = len(sequences.iloc[0]['sequence'])
    characters     = list(range(num_characters))
    characters_to_idx = {c: i for i, c in enumerate(characters)}

    def leaf_f(c, v):
        return sequences.loc[v]['sequence'][characters_to_idx[c]]

    reticulation_nodes = set([n for n in N.nodes if N.in_degree(n) > 1])
    reticulation_edges = set([(u, v) for u, v in N.edges if N.in_degree(v) > 1])

    # add a dummy root node
    N.add_node("dummy_root")
    N.add_edge("dummy_root", root)

    # i.e. x_{u,v,c,c'} = 1 if u is assigned character c and v is assigned character c'
    model = pyo.ConcreteModel()

    model.decisions = pyo.Var(characters, N.edges, alphabet, alphabet, domain=pyo.NonNegativeReals)
    model.reticulations = pyo.Var(characters, N.edges, domain=pyo.Binary)

    # ensure all non-reticulation edges are included    
    for u, v in N.edges:
        if (u, v) in reticulation_edges:
            continue
        for c in characters:
            model.reticulations[c, (u, v)].fix(1)

    model.edge_constraints = pyo.ConstraintList()
    model.leaf_constraints = pyo.ConstraintList()
    model.reticulation_constraints = pyo.ConstraintList()

    # set reticulation constraints
    for c in characters:
        for v in reticulation_nodes:
            model.reticulation_constraints.add(
                sum(model.reticulations[c, (u, v)] for u in N.predecessors(v)) == 1
            )

    # set flow constraints
    for c in characters:
        for u, v in N.edges:
            for i in alphabet:
                for w in N[v]:
                    model.edge_constraints.add(
                        sum(model.decisions[c, u, v, j, i] for j in alphabet) - sum(model.decisions[c, v, w, i, j] for j in alphabet) <= 2 - model.reticulations[c, (u, v)] - model.reticulations[c, (v, w)]
                    )

                    model.edge_constraints.add(
                        sum(model.decisions[c, u, v, j, i] for j in alphabet) - sum(model.decisions[c, v, w, i, j] for j in alphabet) >= -2 + model.reticulations[c, (u, v)] + model.reticulations[c, (v, w)]
                    )

    # set \sum_{c} x_{u,v,c,c'} = 1 for all e=(u,v), v is a leaf, and c' is the label of v
    # set \sum_{c} x_{u,v,c,c'} = 0 for all e=(u,v), v is a leaf, and c' is not the label of v
    for c in characters:
        for u, v in N.edges:
            if N.out_degree(v) != 0:
                continue

            for i in alphabet:
                if i == leaf_f(c, v):
                    # if y[c, (u, v)] == 1, then \Sigma x_{u,v,c,c'} = 1 
                    # if y[c, (u, v)] == 0, then 2 >= \Sigma x_{u,v,c,c'} >= 0
                    model.leaf_constraints.add(sum(model.decisions[c, u, v, j, i] for j in alphabet) >= model.reticulations[c, (u, v)])
                    model.leaf_constraints.add(sum(model.decisions[c, u, v, j, i] for j in alphabet) <= 2 - model.reticulations[c, (u, v)])
                else:
                    model.leaf_constraints.add(sum(model.decisions[c, u, v, j, i] for j in alphabet) <= 1 - model.reticulations[c, (u, v)])

    weight_f = lambda i, j: 1 if i != j else 0
    model.objective_expr = sum(model.decisions[c, u, v, i, j] * weight_f(i, j) for c in characters for u, v in N.edges for i in alphabet for j in alphabet)
    model.objective = pyo.Objective(expr=model.objective_expr, sense=pyo.minimize)
    return model

def parse_args():
    parser = argparse.ArgumentParser(description='Softwired parsimony problem solver using TLP')
    parser.add_argument("network", help="Phylogenetic network in edgelist format")
    parser.add_argument("sequences", help="Leaf sequences in CSV format")
    parser.add_argument("output", help="Output prefix")
    parser.add_argument("--mode", choices=["fischer", "tlp"], default="tlp", help="Mode to use")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    logger.info("Reading phylogenetic network")
    N = nx.read_edgelist(args.network, nodetype=int, create_using=nx.DiGraph())

    logger.info("Reading sequences")
    sequences = pd.read_csv(args.sequences, index_col=0, dtype=str)

    root = [n for n in N.nodes if N.in_degree(n) == 0]
    if len(root) != 1:
        logger.error("Network must have a single root")
        sys.exit(1)

    root = root[0]

    if not all(N.in_degree(n) <= 2 for n in N.nodes):
        logger.error("Network in-degree of all nodes must be at most 2.")
        sys.exit(1)

    if not nx.is_directed_acyclic_graph(N):
        logger.error("Network must be a directed acyclic graph.")
        sys.exit(1)

    alphabet = set([c for _, row in sequences.iterrows() for s in row['sequence'] for c in s])
    if args.mode == "fischer":
        model = create_softwired_fischer(N, root, sequences, alphabet)
    elif args.mode == "tlp":
        model = create_softwired_tlp(N, root, sequences, alphabet)

    solver = pyo.SolverFactory('gurobi')
    results = solver.solve(model, tee=True)

    results_dict = {}
    results_dict['objective'] = model.objective()
    results_dict['runtime'] = results.solver.wall_time

    with open(args.output + "_results.json", "w") as f:
        json.dump(results_dict, f)
