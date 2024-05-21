import pandas as pd
import numpy as np
import argparse
import json
import networkx as nx

from loguru import logger
from tqdm import tqdm
from collections import defaultdict, deque

def construct_migration_graph(tree, character_set, labeling):
    if len(character_set) == 1:
        G = nx.DiGraph()
        nodes = list(tree.nodes)
        G.add_node(labeling[nodes[0]])
        return G

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

def bruteforce_sampling(T, character_set, leaf_f, prob_f, root, root_label, num_samples=1000):
    """
    A brute force implementation of the sampling algorithm to
    validate our implementation of the stochastic backtrace algorithm.
    """
    labelings = []
    labelings.append({root: root_label})

    for node in T.nodes:
        if node == root:
            continue

        new_labelings = []
        for labeling in labelings:
            if is_leaf(T, node):
                new_labeling = labeling.copy()
                new_labeling[node] = leaf_f(node)
                new_labelings.append(new_labeling)
            else:
                for character in character_set:
                    new_labeling = labeling.copy()
                    new_labeling[node] = character
                    new_labelings.append(new_labeling)

        labelings = new_labelings
    
    probs = []
    for labeling in tqdm(labelings):
        prob = 1
        for edge in T.edges:
            parent, child = edge
            prob *= prob_f(labeling[parent], labeling[child], edge)
        probs.append(prob)

    probs = np.array(probs)
    probs = probs / np.sum(probs)

    samples = []
    for _ in tqdm(range(num_samples)):
        sample = np.random.choice(labelings, p=probs)
        samples.append(sample)

    return samples

def stochastic_backtrace(L, T, character_set, leaf_f, prob_f, root, root_label):
    labeling = {}
    labeling[root] = root_label

    queue = []
    for child in T[root]:
        queue.append(child)

    while queue:
        node = queue.pop(0)

        if is_leaf(T, node):
            labeling[node] = leaf_f(node)
            continue

        parent = list(T.predecessors(node))[0]
        parent_alpha = labeling[parent]

        denominator = 0
        for beta in character_set:
            denominator += prob_f(parent_alpha, beta, (parent, node)) * L[node][beta]

        probabilities = {}
        for alpha in character_set:
            numerator = prob_f(parent_alpha, alpha, (parent, node)) * L[node][alpha]
            probabilities[alpha] = numerator / denominator

        labeling[node] = np.random.choice(list(probabilities.keys()), p=list(probabilities.values()))

        for child in T[node]:
            queue.append(child)

    return labeling

def felsensteins(T, root, character_set, leaf_f, prob_f):
    L = defaultdict(dict) 
    stack = deque([root]) # stack to simulate DFS search

    while stack:
        node = stack.pop()

        if is_leaf(T, node):
            if leaf_f(node) is None:
                raise ValueError(f"Leaf {node} has no label")
                continue

            for char in character_set:
                if char != leaf_f(node):
                    L[node][char] = 0
                else:
                    L[node][char] = 1

            continue

        # all children are visited
        if all(child in L for child in T[node]): 
            for alpha in character_set: 
                prob = 1
                for child in T[node]:
                    child_prob = 0
                    for beta in character_set:
                        child_prob += prob_f(alpha, beta, (node, child)) * L[child][beta]
                    prob *= child_prob
                L[node][alpha] = prob

            continue

        stack.append(node)
        for child in T[node]:
            stack.append(child)

    return L

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Samples labelings according to the conditional distribution on the leaves."
    )

    parser.add_argument(
        "tree", help="Tree in edgelist format"
    )

    parser.add_argument(
        "labels", help="Labels for the leaves"
    )

    parser.add_argument(
        "probabilities", help="Probabilities for the transitions between characters along the edges"
    )

    parser.add_argument(
        "-s", "--samples", help="Number of optimal solutions to randomly sample", type=int, default=1000
    )

    parser.add_argument(
        "-r", "--random-seed", help="Random seed for sampling", type=int, default=42
    )

    parser.add_argument(
        "-l", "--label", help="Label for the root node", default=None
    )

    parser.add_argument(
        "-o", "--output", help="Output file", default="output.json"
    )

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()

    tree = nx.read_edgelist(args.tree, create_using=nx.DiGraph(), data=(("weight", float),))
    labels_csv = pd.read_csv(args.labels, sep=",").set_index("leaf")
    probabilities_csv = pd.read_csv(args.probabilities).set_index(["parent", "child", "label1", "label2"])

    roots = [node for node in tree.nodes if len(list(tree.predecessors(node))) == 0]
    root = roots[0]

    def leaf_f(node):
        if node not in labels_csv.index:
            return None
        y = labels_csv.loc[node, "label"]
        if y == "None":
            return None
        return y

    def prob_f(alpha, beta, e):
        parent, child = e
        return probabilities_csv.loc[(parent, child, alpha, beta), "probability"]
 
    character_set = set(labels_csv["label"].unique())
    
    logger.info(f"Tree has {len(tree.nodes)} nodes")
    logger.info(f"Character set has {len(character_set)} characters")
   
    L = felsensteins(tree, root, character_set, leaf_f, prob_f)
    solutions = []
    np.random.seed(args.random_seed)
    for _ in tqdm(range(args.samples)):
        labeling = stochastic_backtrace(L, tree, character_set, leaf_f, prob_f, root, args.label)
        solutions.append(labeling)

    """ 
    A brute force implementation to validate the stochastic backtrace algorithm: 
        solutions = bruteforce_sampling(tree, character_set, leaf_f, prob_f, root, args.label, num_samples=args.samples)
    """

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
