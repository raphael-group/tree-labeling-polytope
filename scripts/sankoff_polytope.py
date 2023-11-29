import argparse
import itertools

import pandas as pd
import numpy as np
import networkx as nx
import pyomo.environ as pyo
 
from loguru import logger
from tqdm import tqdm
from enum import Enum
from dataclasses import dataclass

def is_leaf(T, node):
    return len(T[node]) == 0

"""
Represents a vertex in the Sankoff directed acyclic
hypergraph.
"""
@dataclass(frozen=True)
class State:
    node : str
    character : str

def create_sankoff_polytope(T, root, character_set, leaf_f, dist_f):
    model = pyo.ConcreteModel()

    logger.info(f"Creating Sankoff polytope for tree with {len(T.nodes)} nodes and {len(T.edges)} edges.")

    logger.info("Creating decision variables.")
    # map from (j, J) -> decision variable {0, 1}
    # where J is a set of states
    decisions = [] 
    for node in T.nodes:
        if is_leaf(T, node):
            continue

        for char in character_set:
            j = State(node, char)

            internal_children = [c for c in T[node] if not is_leaf(T, c)]
            for states in itertools.product(character_set, repeat=len(internal_children)):
                J = frozenset(map(lambda x: State(x[0], x[1]), zip(internal_children, states)))
                decisions.append((j, J))

    model.decisions = pyo.Var(decisions, domain=pyo.NonNegativeReals)

    logger.info("Computing decision costs.")
    # map from (j, J) -> \mathbb{R}^+ 
    decision_costs = {}
    for head, J in decisions:
        cost = 0
        for tail in J:
            cost += dist_f(head.character, tail.character)

            if all(is_leaf(T, c) for c in T[tail.node]):
                for c in T[tail.node]:
                    cost += dist_f(tail.character, leaf_f(c))
        
        for c in T[head.node]:
            if not is_leaf(T, c): continue
            cost += dist_f(head.character, leaf_f(c))

        decision_costs[(head, J)] = cost

    # create a map from j to the set of decision variables whose tail contains j
    j_to_J = {}
    for l, J in decisions:
        for j in J:
            if j not in j_to_J:
                j_to_J[j] = []
            j_to_J[j].append((l, J))

    logger.info("Adding flow constraints.")
    # add flow constraints
    model.flows = pyo.ConstraintList()
    for node in T.nodes:
        if is_leaf(T, node):
            continue

        if node == root:
            continue

        internal_children = [c for c in T[node] if not is_leaf(T, c)]

        # in this case, require 1 unit of out flow over all states
        if len(internal_children) == 0:
            out_flow_expr = 0
            for char in character_set:
                j = State(node, char)
                out_flow_expr += sum(model.decisions[(l, J)] for l, J in j_to_J[j]) 

            model.flows.add(out_flow_expr == 1)
            continue

        for char in character_set:
            j = State(node, char)
            # else, flow into j must be equal to flow out of j
            in_flow_expr, out_flow_expr = 0, 0
            for l, J in j_to_J[j]:
                out_flow_expr += model.decisions[(l, J)]

            for states in itertools.product(character_set, repeat=len(internal_children)):
                J = frozenset(map(lambda x: State(x[0], x[1]), zip(internal_children, states)))
                in_flow_expr += model.decisions[(j, J)]

            model.flows.add(
                in_flow_expr == out_flow_expr
            )

    # let in flow be 1 for the root
    model.flows.add(
        sum(model.decisions[(j, J)] for j, J in decisions if j.node == root) == 1
    )

    logger.info("Adding objective function.")
    # add objective function (minimize the sum of the costs of all decisions)
    model.objective = pyo.Objective(
        expr=sum(model.decisions[(j, J)] * decision_costs[(j, J)] for j, J in decisions),
        sense=pyo.minimize
    )

    return model, decision_costs

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Performs maximum parsimony inference by solving a linear program over the Sankoff polytope."
    )

    parser.add_argument(
        "tree", help="Tree in edgelist format"
    )

    parser.add_argument(
        "labels", help="Leaf labeling as a TSV file"
    )

    parser.add_argument(
        "-r", "--root", help="Root of the tree", default="root"
    )

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()

    tree = nx.read_edgelist(args.tree, create_using=nx.DiGraph(), nodetype=str)
    labels_tsv = pd.read_csv(args.labels, sep="\t", header=None, names=["id", "label"]).set_index("id")

    def dist_f(x, y):
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

    model, decision_costs = create_sankoff_polytope(tree, args.root, character_set, leaf_f, dist_f)

    decision = (State(node='A', character='SBwl'), frozenset({State(node='B', character='SBwl'), State(node='D', character='ROv')}))
    
    solver = pyo.SolverFactory('gurobi')
    solver.solve(model, tee=True)

    # print non-zero decision variables
    for j, J in model.decisions:
        if model.decisions[(j, J)].value > 0:
            print(f"{j} {J} {model.decisions[(j, J)].value}")
            print(f"cost: {decision_costs[(j, J)]}")

