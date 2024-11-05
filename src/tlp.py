import argparse
import json

import pandas as pd
import numpy as np
import networkx as nx
 
from loguru import logger
from collections import defaultdict

from fastmachina import *
from convexrecoloring import *

def is_leaf(T, node):
    return len(T[node]) == 0

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
                                choices=["tree", "dag", "none"],
                                default="tree")
    fast_machina_parser.add_argument("-o", "--output", help="Output prefix", default="result")
    fast_machina_parser.add_argument("-r", "--root", help="Root of the tree", default=None)
    fast_machina_parser.add_argument("-l", "--label", help="Root label", default=None)
    fast_machina_parser.add_argument("-w", "--weights", help="Weight of transitioning between labels", default=None)

    # (weighted) convex recoloring subparser
    convex_recoloring = subparsers.add_parser("convex_recoloring", help="Convex Recoloring")
    convex_recoloring.add_argument("tree", help="Tree in edgelist format")
    convex_recoloring.add_argument("labels", help="Leaf labeling as a CSV file")
    convex_recoloring.add_argument("-o", "--output", help="Output prefix", default="result")
    convex_recoloring.add_argument("-r", "--root", help="Root of the tree", default="root")
    convex_recoloring.add_argument("-k", help="Weighted parsimony constraint", default=None, type=float)
    convex_recoloring.add_argument("-w", "--weights", help="Weight of transitioning between labels", default=None)

    # softwired small parsimony subparser
    softwired_sp = subparsers.add_parser("softwired", help="Softwired Small Parsimony on a Phylogenetic Network")
    softwired_sp.add_argument("network", help="Phylogenetic network in edgelist format")
    softwired_sp.add_argument("sequences", help="Leaf sequences in CSV format")
    softwired_sp.add_argument("output", help="Output prefix")

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
        lp_obj = 0
        obj = 0
    else:
        # computes the vertex labeling using the specified method
        if args.method == "fast_machina":
            lp_obj = None
            vertex_labeling, obj = fast_machina(tree, character_set, leaf_f, dist_f, root, args)
        elif args.method == "convex_recoloring":
            _, lp_obj = convex_recoloring(tree, character_set, leaf_f, dist_f, root, args, integral=False)
            vertex_labeling, obj = convex_recoloring(tree, character_set, leaf_f, dist_f, root, args)

    # write the objective value to a file (json)
    with open(f"{args.output}_results.json", "w") as f:
        results = {}
        results["objective"] = obj
        if lp_obj is not None:
            results["lp_relaxation_objective"] = lp_obj
        f.write(json.dumps(results))

    # writes an optimal labeling to a file 
    with open(f"{args.output}_vertex_labeling.csv", "w") as f:
        f.write("vertex,label\n")
        for node, label in vertex_labeling.items():
            f.write(f"{node},{label}\n")

    if args.method == "fast_machina":
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
