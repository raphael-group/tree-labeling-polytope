import argparse
import itertools

import pandas as pd
import numpy as np
import networkx as nx
 
from loguru import logger
from tqdm import tqdm
from enum import Enum
from collections import defaultdict

"""Samples from a uniform distribution over the set {a, a+1, ..., b}"""
def rand_int(rng, a, b):
    return rng.integers(a, b+1)

def rand_predecessor(node, predecessors, weights, rng):
    weighted_preds = [(p, weights[(p, node)]) for p in predecessors]
    total_weight = sum(weight for _, weight in weighted_preds)
    r = rng.random() * total_weight
    upto = 0
    for pred, weight in weighted_preds:
        if upto + weight >= r:
            return pred
        upto += weight
    assert False, "The procedure will not reach here."

"""
Implementation of Wilson's algorithm for generating a random 
spanning tree from an arbitrary graph.
"""
def sample_random_spanning_tree(G, weights, rng, root=None):
    spanning_tree = nx.DiGraph()

    for u in G.nodes:
        spanning_tree.add_node(u)

    next_node = [-1] * len(G.nodes)
    in_tree = [False] * len(G.nodes)

    if root is None:
        root = rand_int(rng, 0, len(G.nodes) - 1)

    in_tree[root] = True
    for u in G.nodes:
        if in_tree[u]:
            continue

        v = u
        while not in_tree[v]:
            pred = list(G.predecessors(v))
            if len(pred) == 0:
                raise RuntimeError("Graph is not strongly connected")

            next_node[v] = rand_predecessor(v, pred, weights, rng)
            v = next_node[v]

        v = u
        while not in_tree[v]:
            in_tree[v] = True
            spanning_tree.add_edge(next_node[v], v)
            v = next_node[v]

    return spanning_tree, root

"""
Sample a random migration tree with m nodes rooted at 
the 0 vertex.
"""
def sample_migration_tree(m, seed):
    G = nx.DiGraph()
    for i in range(m):
        G.add_node(i)

    for i in range(m):
        for j in range(m):
            if i != j:
                G.add_edge(i, j)

    rng = np.random.default_rng(seed)
    tree, _ = sample_random_spanning_tree(G, {(i, j): 1 for i in range(m) for j in range(m) if i != j}, rng, root=0)
    return tree

"""
Sample a random DAG with m nodes using the algorithm of
`Random Generation of Directed Acyclic Graphs` by 
Melancon et al.
"""
def sample_random_dag(m, seed, chain_length=None):
    chain_length = 10*(m**2) if chain_length is None else chain_length
    G = nx.DiGraph()
    for i in range(m):
        G.add_node(i)

    rng = np.random.default_rng(seed)
    for i in range(chain_length):
        u = rand_int(rng, 0, m-1)
        v = rand_int(rng, 0, m-1)

        if u == v:
            continue

        had_edge = G.has_edge(u, v)
        G.add_edge(u, v)
        if list(nx.simple_cycles(G)):
            G.remove_edge(u, v)
            continue

        if had_edge:
            G.remove_edge(u, v)
            continue

    return G

"""
Sample a random labeling of the vertices of a tree
using a random walk over the edges of the tree.
"""
def sample_random_labeling(tree, tree_root, migration_graph, migration_graph_root, rng, monoclonal=False, prob=0.1):
    labeling = {}
    labeling[tree_root] = migration_graph_root

    def random_ordering(neighbors):
        neighbors = list(neighbors)
        order = np.arange(len(neighbors))
        rng.shuffle(order)
        return [neighbors[i] for i in order]

    migrations_used = set()
    for u, v in nx.bfs_edges(tree, source=tree_root, sort_neighbors=random_ordering):
        if rng.random() > prob:
            labeling[v] = labeling[u]
            continue

        lu = labeling[u]
        successors = list(migration_graph.successors(lu))
        if len(successors) == 0:
            labeling[v] = lu
            continue

        lv = rng.choice(successors)
        if monoclonal and (lu, lv) in migrations_used:
            labeling[v] = lu
            continue

        labeling[v] = lv
        migrations_used.add((lu, lv))

    return labeling

def parse_args():
    parser = argparse.ArgumentParser(description="Simulate metastatic cancer evolution along a phylogenetic tree.")
    parser.add_argument("tree", help="Tree in edgelist format")
    parser.add_argument("root", help="Root of the tree")
    parser.add_argument("-o", "--output", help="Output prefix", default="result")
    parser.add_argument("-m", help="Number of labels", type=int, default=6)
    parser.add_argument("-r", "--random-seed", help="Random seed", type=int, default=0)
    parser.add_argument(
        "-s", "--structure", help="Migration graph structure",
        choices=["polyclonal_tree", "polyclonal_dag", "monoclonal_tree", "monoclonal_dag"],
        default="polyclonal_dag"
    )
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    tree = nx.read_edgelist(args.tree, nodetype=str, create_using=nx.DiGraph())
    character_set = range(args.m)

    if args.structure == "polyclonal_tree" or args.structure == "monoclonal_tree":
        migration_tree = sample_migration_tree(args.m, args.random_seed)
        rng = np.random.default_rng(args.random_seed)
        labeling = sample_random_labeling(tree, args.root, migration_tree, 0, rng, monoclonal=args.structure == "monoclonal_tree")
        migration_graph = migration_tree
    else:
        migration_dag = sample_random_dag(args.m, args.random_seed)
        migration_dags = list(nx.weakly_connected_components(migration_dag))
        migration_dag = nx.subgraph(migration_dag, max(migration_dags, key=len))
        dag_root = next(nx.topological_sort(migration_dag))

        rng = np.random.default_rng(args.random_seed)
        labeling = sample_random_labeling(tree, args.root, migration_dag, dag_root, rng, monoclonal=args.structure == "monoclonal_dag")
        migration_graph = migration_dag

    with open(f"{args.output}_labeling.csv", "w") as f:
        f.write("vertex,label\n")
        for node, label in labeling.items():
            f.write(f"{node},{label}\n")

    with open(f"{args.output}_leaf_labeling.csv", "w") as f:
        f.write("leaf,label\n")
        for node, label in labeling.items():
            if tree.out_degree(node) == 0:
                f.write(f"{node},{label}\n")

    with open(f"{args.output}_migration_graph.csv", "w") as f:
        f.write("src,dst\n")
        for (i, j) in migration_graph.edges:
            f.write(f"{i},{j}\n")




