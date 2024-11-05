import itertools
import networkx as nx
import pyomo.environ as pyo
import gurobipy as gp
import numpy as np

from gurobipy import GRB
from loguru import logger
from collections import defaultdict

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
    logger.info(f"Character set has size: {len(character_set)}")

    # add a dummy root node
    T.add_node("dummy_root")
    T.add_edge("dummy_root", root)

    # i.e. x_{u,v,c,c'} = 1 if u is assigned character c and v is assigned character c'
    model = pyo.ConcreteModel()
    # set domain to be [0, 1]
    model.decisions = pyo.Var(T.edges, character_set, character_set, domain=pyo.NonNegativeReals, initialize=0)

    # require \sum_{c'} x_{u,v,c,c'} = \sum_{c'}x_{v,w,c',c} for all u,v,w,c
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
    model.objective_expr = sum(model.decisions[u, v, c1, c2] * dist_f((u, v), c1, c2) for u, v in T.edges for c1 in character_set for c2 in character_set)
    model.objective = pyo.Objective(expr=model.objective_expr, sense=pyo.minimize)
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
by adding constraints to the model using a Gurobi callback.
"""
def setup_fast_machina_constraints(solver, model, character_set, constraint_type):
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
            S = nx.find_cycle(G, orientation="original")
            S = [i for (i,_,_) in S]

            logger.info(f"Adding constraint for cycle {S}.")

            cons = model.migration_graph_constraints.add(
                sum(model.migrations[(S[i], S[(i+1) % len(S)])] for i in range(len(S))) <= len(S) - 1
            )

            gb_model.cbLazy(cons)
        except nx.exception.NetworkXNoCycle:
            pass

    model.migration_graph_constraints = pyo.ConstraintList()
    if constraint_type == "tree":
        S = character_set
        model.migration_graph_constraints.add(
            sum(model.migrations[(i, j)] for i in S for j in S if i != j) == len(S) - 1
        )
        solver.set_instance(model)
    elif constraint_type == "dag":
        solver.set_instance(model)
        solver.set_gurobi_param("LazyConstraints", 1)
        solver.set_callback(dag_callback)
    else:
        raise ValueError(f"Unknown constraint type: {constraint_type}")

    return

def fast_machina_prune(tree, character_set, leaf_f, dist_f, root):
    tree = tree.copy()

    Bv = defaultdict()
    for v in nx.dfs_postorder_nodes(tree, source=root):
        if is_leaf(tree, v):
            Bv[v] = set([leaf_f(v)])
        else:
            Bv[v] = set()
            for u in tree[v]:
                Bv[v] = Bv[v] | Bv[u]

    node_to_ancestor_map = {}
    while True:
        removed = False
        for v in nx.dfs_preorder_nodes(tree, source=root):
            if len(Bv[v]) == 1 and tree.out_degree(v) != 0:
                descendants = nx.descendants(tree, v)
                for u in descendants:
                    tree.remove_node(u)
                    node_to_ancestor_map[u] = v
                removed = True
                break
        if not removed:
            break

    leaf_dict = {v:list(Bv[v])[0] for v in tree.nodes if is_leaf(tree, v)}
    logger.info(f"Removed {len(node_to_ancestor_map)} nodes from the tree.")
    leaf_f = lambda v: leaf_dict[v]
    return tree, leaf_f, node_to_ancestor_map

"""
Solves a generalization of the MACHINA parsimonious migration 
history problem using the tree labeling polytope.
"""
def fast_machina(tree, character_set, leaf_f, dist_f, root, args, mip_gap=0.15):
    tree = tree.copy()

    """ Step 1: Prune the tree to remove unnecessary nodes """
    tree, leaf_f, node_to_ancestor_map = fast_machina_prune(tree, character_set, leaf_f, dist_f, root)

    """ Step 2: Solve and round the LP relaxation to obtain an initial solution """
    if args.constraints != "none":
        logger.info("Solving LP relaxation of model to obtain initial feasible solution")
        solver = pyo.SolverFactory('gurobi_persistent')

        # setup relaxed model
        model = create_tree_labeling_polytope(tree, root, character_set, leaf_f, dist_f, root_label=args.label)
        append_migrations(model, tree, character_set)
        setup_fast_machina_constraints(solver, model, character_set, args.constraints)

        for i, j in model.migrations:
            model.migrations[i, j].domain = pyo.NonNegativeReals
            solver.update_var(model.migrations[i, j])

        solver.solve(model, tee=True)

        value_counts = defaultdict(int)
        for i, j in model.migrations:
            value_counts[model.migrations[i, j].value] += 1

        # round the LP relaxation to obtain an initial set of migrations
        candidate_migrations = [(i, j) for i, j in model.migrations if model.migrations[i, j].value > 1e-7]

    """ Step 3: Setup the MILP using Gurobi callbacks, if necessary """
    solver = pyo.SolverFactory('gurobi_persistent')

    if "dummy_root" in tree.nodes:
        tree.remove_node("dummy_root")

    model = create_tree_labeling_polytope(tree, root, character_set, leaf_f, dist_f, root_label=args.label)
    if args.constraints != "none":
        append_migrations(model, tree, character_set)

    added_constraints = []
    if args.constraints != "none":
        setup_fast_machina_constraints(solver, model, character_set, args.constraints)

        candidate_migrations = set(candidate_migrations)
        model.heur_constraints = pyo.ConstraintList()
        for i, j in model.migrations:
            if (i, j) not in candidate_migrations:
                con = model.heur_constraints.add(model.migrations[i, j] == 0)
                solver.add_constraint(con)
                added_constraints.append(con)
    else:
        solver.set_instance(model)

    """ Step 4: Solve the full MILP """
    logger.info("Solving full model")
    solver.options["MIPGap"] = mip_gap
    solver.solve(model, tee=True, warmstart=True)

    # remove heuristic constraints and re-solve
    if args.constraints != "none":
        for c in added_constraints:
            solver.remove_constraint(c)
        model.heur_constraints.clear()
        solver.solve(model, tee=True, warmstart=True)

    # compute (an) optimal vertex labeling
    vertex_labeling = {}
    for u, v in tree.edges:
        for c1, c2 in itertools.product(character_set, character_set):
            if np.abs(model.decisions[u, v, c1, c2]()) > 1e-4:
                vertex_labeling[u] = c1
                vertex_labeling[v] = c2

    for u in node_to_ancestor_map:
        vertex_labeling[u] = vertex_labeling[node_to_ancestor_map[u]]

    return vertex_labeling, model.objective()
