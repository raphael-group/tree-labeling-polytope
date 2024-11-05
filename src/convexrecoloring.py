import itertools
import networkx as nx
import numpy as np
import pyomo.environ as pyo
import gurobipy as gp

from loguru import logger

def is_leaf(T, node):
    return len(T[node]) == 0

def convex_recoloring(tree, character_set, leaf_f, dist_f, root, args, integral=True):
    T = tree.copy()

    if args.k is None:
        k = len(character_set) - 1

    model = pyo.ConcreteModel()

    logger.info(f"Creating tree labeling polytope for tree with {len(T.nodes)} nodes and {len(T.edges)} edges.")
    logger.info(f"Character set has size: {len(character_set)}")

    # add a dummy root node
    T.add_node("dummy_root")
    T.add_edge("dummy_root", root)

    pendant_edges = [(u, v) for u, v in T.edges if is_leaf(T, v)]

    # i.e. x_{u,v,c,c'} = 1 if u is assigned character c and v is assigned character c'
    model = pyo.ConcreteModel()

    # set domain to be [0, 1]
    model.decisions = pyo.Var(T.edges, character_set, character_set, domain=pyo.NonNegativeReals, initialize=0)
    if integral:
        model.relabelings = pyo.Var(pendant_edges, domain=pyo.Binary, initialize=0)
    else:
        model.relabelings = pyo.Var(pendant_edges, domain=pyo.NonNegativeReals, bounds=(0, 1))

    # require \sum_{c'} x_{u,v,c,c'} = \sum_{c'}x_{v,w,c',c} for all u,v,w,c
    model.edge_constraints = pyo.ConstraintList()
    for u, v in T.edges:
        for c in character_set:
            for w in T[v]:
                model.edge_constraints.add(
                    sum(model.decisions[u, v, c2, c] for c2 in character_set) - sum(model.decisions[v, w, c, c2] for c2 in character_set) == 0
                )

    # require \sum_{c,c'} x_{u,v,c,c'} = 1 for all e=(u,v)
    model.sum_constraints = pyo.ConstraintList()
    for u, v in T.edges:
        model.sum_constraints.add(sum(model.decisions[u, v, c1, c2] for c2 in character_set for c1 in character_set) == 1)

    # require leaves that are not relabeled to have the correct label
    model.leaf_constraints = pyo.ConstraintList()
    for u, v in pendant_edges:
        for c in character_set:
            if c == leaf_f(v):
                model.leaf_constraints.add(sum(model.decisions[u, v, c2, c] for c2 in character_set) >= model.relabelings[u, v])

    # c^T x_{u,v,i,j} \leq k
    model.sum_constraints.add(
        sum(model.decisions[u, v, c1, c2] * dist_f((u, v), c1, c2) for u, v in T.edges for c1 in character_set for c2 in character_set) <= k
    )

    for color in character_set:
        leaves = [v for v in T.nodes if is_leaf(T, v) and leaf_f(v) == color]
        if len(leaves) == 0: continue
        model.sum_constraints.add(
            sum(model.relabelings[u, v] for u, v in pendant_edges if v in leaves) >= 1
        )

    model.objective_expr = sum(1 - model.relabelings[u, v] for u, v in pendant_edges)
    model.objective = pyo.Objective(expr=model.objective_expr, sense=pyo.minimize)

    solver = pyo.SolverFactory('gurobi_persistent')
    solver.set_instance(model)
    solver.solve(model, tee=True, warmstart=True)

    vertex_labeling = {}
    for u, v in T.edges:
        for c1, c2 in itertools.product(character_set, character_set):
            if np.abs(model.decisions[u, v, c1, c2]()) > 1e-4:
                vertex_labeling[u] = c1
                vertex_labeling[v] = c2
    
    return vertex_labeling, model.objective()
