import pandas as pd
import numpy as np
import argparse
import json
import networkx as nx

from loguru import logger
from tqdm import tqdm
from collections import defaultdict, deque

def construct_migration_graph(tree, character_set, labeling):
    G = nx.DiGraph()
    for u, v in tree.edges:
        if labeling[u] == labeling[v]: 
            continue
        if (labeling[u], labeling[v]) in G.edges:
            G[labeling[u]][labeling[v]]["weight"] += 1
            continue
        G.add_edge(labeling[u], labeling[v], weight=1)

    return G

def is_leaf(T, node):
    return len(T[node]) == 0

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

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Performs maximum parsimony."
    )

    parser.add_argument(
        "tree", help="Tree in edgelist format"
    )

    parser.add_argument(
        "labels", help="Labels for the leaves"
    )

    parser.add_argument(
        "-s", "--samples", help="Number of optimal solutions to randomly sample", type=int, default=1000
    )

    parser.add_argument(
        "-r", "--random-seed", help="Random seed for sampling", type=int, default=42
    )

    parser.add_argument(
        "-o", "--output", help="Output file", default="output.json"
    )

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()

    tree = nx.read_edgelist(args.tree, create_using=nx.DiGraph(), data=(("weight", float),))
    labels_csv = pd.read_csv(args.labels, sep=",").set_index("leaf")

    roots = [node for node in tree.nodes if len(list(tree.predecessors(node))) == 0]
    root = roots[0]

    def dist_f(x, y):
        if x is None or y is None:
            return 0
        if x == y:
            return 0
        return 1

    def leaf_f(node):
        if node not in labels_csv.index:
            return None
        y = labels_csv.loc[node, "label"]
        if y == "None":
            return None
        return y
 
    character_set = set(labels_csv["label"].unique())
    scores = mp(tree, root, character_set, leaf_f, dist_f)
    
    logger.info(f"Tree has {len(tree.nodes)} nodes")
    logger.info(f"Character set has {len(character_set)} characters")
   
    solutions = []
    np.random.seed(args.random_seed)
    for _ in range(args.samples):
        solution = {}
        for node in nx.dfs_preorder_nodes(tree, root): 
            if node == root: 
                minimum_cost = min(scores[root].values())
                solution[node] = np.random.choice([char for char in character_set if scores[root][char] == minimum_cost])
                continue
            
            parent = list(tree.predecessors(node))[0]
            parent_char = solution[parent]
            min_cost = min([scores[node][char] + dist_f(parent_char, char) for char in character_set])
            possible_chars = [char for char in character_set if scores[node][char] + dist_f(parent_char, char) == min_cost]
            solution[node] = np.random.choice(possible_chars)
        solutions.append(solution)

    logger.info(f"Sampled {len(solutions)} optimal solutions")
    migration_graph_statistics = {
        "is_polyclonal_tree": 0,
        "is_polyclonal_dag": 0,
        "has_cycles": 0,
        "num_sampled": len(solutions),
        "num_nodes": len(tree.nodes),
        "num_characters": len(character_set),
        "topologies": []
    }

    migration_graph_topologies = defaultdict(int)

    for solution in solutions:
        G = construct_migration_graph(tree, character_set, solution)

        migration_graph_topologies[frozenset(G.edges.data("weight"))] += 1

        if nx.is_tree(G):
            migration_graph_statistics["is_polyclonal_tree"] += 1
        elif nx.is_directed_acyclic_graph(G):
            migration_graph_statistics["is_polyclonal_dag"] += 1
        else:
            migration_graph_statistics["has_cycles"] += 1

    for topology, count in migration_graph_topologies.items():
        edge_list = []
        for u, v, weight in topology:
            edge_list.append({"src": u, "dst": v, "multiplicity": int(weight)})

        migration_graph_statistics[f"topologies"].append({
            "edges": edge_list,
            "count": count
        })

    with open(args.output, "w") as f:
        json.dump(migration_graph_statistics, f, indent=4)
