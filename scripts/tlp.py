import argparse
import itertools
import json

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
    if constraint_type == "polyclonal_tree" or constraint_type == "monoclonal_tree":
        S = character_set
        model.migration_graph_constraints.add(
            sum(model.migrations[(i, j)] for i in S for j in S if i != j) == len(S) - 1
        )
        solver.set_instance(model)
    elif constraint_type == "polyclonal_dag" or constraint_type == "monoclonal_dag":
        solver.set_instance(model)
        solver.set_gurobi_param("LazyConstraints", 1)
        solver.set_callback(dag_callback)
    else:
        raise ValueError(f"Unknown constraint type: {constraint_type}")

    return

def fast_machina_prune(tree, character_set, leaf_f, dist_f, root):
    tree = tree.copy()

    # TODO: 
    #   1. When polyclonal tree, replace Bv with shortest chain from 
    #      root to leaf using Lemma 4 in the paper.
    #   2. Explore monoclonal tree and monoclonal dag settings.

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
Solves a generalization of the convex recoloring
problem using the tree labeling polytope.
"""
def parsimonious_relabeling(tree, character_set, leaf_f, dist_f, root, args, mip_gap=0.15):
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
    model.relabelings = pyo.Var(pendant_edges, domain=pyo.Binary, initialize=0)

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

    model.objective_expr = sum(model.relabelings[u, v] for u, v in pendant_edges)
    model.objective = pyo.Objective(expr=model.objective_expr, sense=pyo.maximize)

    solver = pyo.SolverFactory('gurobi_persistent')
    solver.set_instance(model)
    solver.solve(model, tee=True, warmstart=True)

    # print out the variables relabelings
    for u, v in pendant_edges:
        if model.relabelings[u, v].value < 0.5:
            logger.info(f"Relabeling of {v}: {model.relabelings[u, v].value}")

    vertex_labeling = {}
    for u, v in T.edges:
        for c1, c2 in itertools.product(character_set, character_set):
            if np.abs(model.decisions[u, v, c1, c2]()) > 1e-4:
                vertex_labeling[u] = c1
                vertex_labeling[v] = c2
    
    return vertex_labeling, model.objective()

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
    if args.constraints.startswith("monoclonal"):
        for c1, c2 in model.migrations:
            model.migration_constraints.add(
                sum(model.decisions[u, v, c1, c2] for u, v in tree.edges) <= model.migrations[c1, c2]
            )

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
    logger.info("Solving full MILP model")
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

def exact_tnet(tree, character_set, leaf_f, dist_f, root, args):
    solver = pyo.SolverFactory('gurobi')

    model = create_tree_labeling_polytope(tree, root, character_set, leaf_f, dist_f, root_label=args.label)
    solver.solve(model, tee=True, warmstart=True)

    obj = model.objective()

    # require that we are on face of the polytope with minimum objective value
    model.face_constraint = pyo.Constraint(expr=model.objective_expr == obj)

    logger.info(f"Optimal parsimony score: {obj}")
    logger.info(f"Adding back transmission constraints.")

    model.back_transmissions = pyo.Var(character_set, domain=pyo.Binary, initialize=1)
    model.back_transmission_constraints = pyo.ConstraintList()

    def back_transmission_callback(model, gb_model, where):
        if where != GRB.Callback.MIPSOL:
            return

        # load solutions
        gb_model.cbGetSolution([model.decisions[u, v, c1, c2] for u, v in tree.edges for c1 in character_set for c2 in character_set if c1 != c2])

        tree_labeling = {}
        for u, v in tree.edges:
            for c1, c2 in itertools.product(character_set, character_set):
                if model.decisions[u, v, c1, c2].value > 0.5:
                    tree_labeling[u] = c1
                    tree_labeling[v] = c2

        character, node = None, None

        B = defaultdict(lambda: defaultdict(-1))
        # three states: -1, 0, 1, indicating have not seen c, have seen c but still c, have seen c but no longer c
        for u, v in nx.bfs_edges(tree, source=root):
            for c in character_set:
                if tree_labeling[v] != c:
                    if B[u][c] == 0: 
                        B[v][c] = 1

                if tree_labeling[v] == c:
                    if B[u][c] == 1:
                        character, node = c, v
                        break
                    else:
                        B[v][c] = 0

        if character is None:
            return

        logger.info(f"Adding back transmission constraint for {character} at {node}.")
        u, v = list(tree.predecessors(node))[0], node
        node = u
        while node != root:
            parent = list(tree.predecessors(node))[0]
            for c1 in character_set:
                for c2 in character_set:
                    if c1 != c and c2 != c:
                        cons = model.back_transmission_constraints.add(
                            model.decisions[parent, node, character, c1] + model.decisions[u, v, c2, character] - 1 <= model.back_transmissions[character]
                        )

                        gb_model.cbLazy(cons)

    # for (u, v) in tree.edges:
        # e = (u, v)
        # child_edges = nx.bfs_edges(tree, source=v)
        # for e_prime in child_edges:
            # if e_prime == e:
                # continue
            # gener = ((i, j1, j2) for i in character_set for j1 in character_set for j2 in character_set if j1 != i or j2 != i)
            # for i, j1, j2 in gener:
                # model.back_transmission_constraints.add(
                    # model.decisions[u, v, i, j1] + model.decisions[e_prime[0], e_prime[1], j2, i] - 1 <= model.back_transmissions[i]
                # )

    # minimize the number of back transmissions
    model.bt_objective = pyo.Objective(expr=sum(model.back_transmissions[c] for c in character_set), sense=pyo.minimize)
    model.objective.deactivate()
    model.bt_objective.activate()

    # set decision variables to be binary
    for i in character_set:
        for j in character_set:
            for e in tree.edges:
                model.decisions[e[0], e[1], i, j].domain = pyo.Binary

    solver.solve(model, tee=True, warmstart=True)

    # print out back transmission variables
    for i in character_set:
        logger.info(f"Back transmission of {i}: {model.back_transmissions[i].value}")

     # compute (an) optimal vertex labeling
    vertex_labeling = {}
    for u, v in tree.edges:
        for c1, c2 in itertools.product(character_set, character_set):
            if np.abs(model.decisions[u, v, c1, c2]()) > 1e-4:
                vertex_labeling[u] = c1
                vertex_labeling[v] = c2

    return vertex_labeling, model.objective()

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
    fast_machina_parser.add_argument("-r", "--root", help="Root of the tree", default=None)
    fast_machina_parser.add_argument("-l", "--label", help="Root label", default=None)
    fast_machina_parser.add_argument("-w", "--weights", help="Weight of transitioning between labels", default=None)

    # exactTNET subparser
    exact_tnet_parser = subparsers.add_parser("exact_tnet", help="exactTNET")
    exact_tnet_parser.add_argument("tree", help="Tree in edgelist format")
    exact_tnet_parser.add_argument("labels", help="Leaf labeling as a CSV file")
    exact_tnet_parser.add_argument("-o", "--output", help="Output prefix", default="result")
    exact_tnet_parser.add_argument("-r", "--root", help="Root of the tree", default="root")
    exact_tnet_parser.add_argument("-l", "--label", help="Root label", default=None)

    # parsimonious relabeling subparser
    parsimonious_relabeling = subparsers.add_parser("parsimonious_relabeling", help="parsimoniousRelabeling")
    parsimonious_relabeling.add_argument("tree", help="Tree in edgelist format")
    parsimonious_relabeling.add_argument("labels", help="Leaf labeling as a CSV file")
    parsimonious_relabeling.add_argument("-o", "--output", help="Output prefix", default="result")
    parsimonious_relabeling.add_argument("-r", "--root", help="Root of the tree", default="root")
    parsimonious_relabeling.add_argument("-k", help="Weighted parsimony constraint", default=None, type=float)
    parsimonious_relabeling.add_argument("-w", "--weights", help="Weight of transitioning between labels", default=None)

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_arguments()

    try:
        tree = nx.read_edgelist(args.tree, create_using=nx.DiGraph(), data=(("weight", float),))
    except Exception as e:
        tree = nx.read_edgelist(args.tree, create_using=nx.DiGraph())

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

    # TODO: check triangle inequality for weights when provided
    if args.weights is not None:
        weights = pd.read_csv(args.weights).set_index(["parent", "child", "label1", "label2"])

        def dist_f(e, x, y):
            if e[0] == "dummy_root":
                return 0

            prob = weights.loc[(e[0], e[1], x, y), "probability"]
            return -np.log(prob) if prob > 0 else 0
    else:
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

    if not hasattr(args, "label"):
        args.label = None

    if args.label is not None and args.label not in character_set:
        logger.warning(f"Root label {args.label} not in character set, removing it.")
        args.label = None

    if len(character_set) == 1:
        logger.warning("Character set has size 1, inferring trivial labeling.")
        vertex_labeling = {node: character_set[0] for node in tree.nodes}
        obj = 0
    else:
        # computes the vertex labeling using the specified method
        if args.method == "fast_machina":
            vertex_labeling, obj = fast_machina(tree, character_set, leaf_f, dist_f, root, args)
        elif args.method == "exact_tnet":
            vertex_labeling = exact_tnet(tree, character_set, leaf_f, dist_f, root, args)
        elif args.method == "parsimonious_relabeling":
            vertex_labeling, obj = parsimonious_relabeling(tree, character_set, leaf_f, dist_f, root, args)

    # write the objective value to a file (json)
    with open(f"{args.output}_results.json", "w") as f:
        results = {}
        results["objective"] = obj
        f.write(json.dumps(results))

    # writes an optimal labeling to a file 
    with open(f"{args.output}_vertex_labeling.csv", "w") as f:
        f.write("vertex,label\n")
        for node, label in vertex_labeling.items():
            f.write(f"{node},{label}\n")

    # compute the migration multi-graph
    migration_graph = defaultdict(int)
    for u, v in tree.edges:
        if u not in vertex_labeling or v not in vertex_labeling:
            continue

        if vertex_labeling[u] != vertex_labeling[v]:
            migration_graph[(vertex_labeling[u], vertex_labeling[v])] += 1

    # write migration multi-graph to a file 
    with open(f"{args.output}_migration_graph.csv", "w") as f:
        f.write("source,target,count\n")
        for (i, j), count in migration_graph.items():
            f.write(f"{i},{j},{count}\n")
