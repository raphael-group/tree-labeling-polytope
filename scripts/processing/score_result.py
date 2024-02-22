import argparse
import json
import re
import sys
import os

import pandas as pd
import numpy as np
import networkx as nx

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from maximum_parsimony import mp

def parse_args():
    parser = argparse.ArgumentParser(description="Score the result of ancestral reconstruction.")
    parser.add_argument("perturbed_tree", help="Perturbed tree in edgelist format")
    parser.add_argument("tree", help="True tree in edgelist format")
    parser.add_argument("vertex_labeling", help="True vertex labeling in CSV format")
    parser.add_argument("inferred_vertex_labeling", help="Inferred vertex labeling in CSV format")
    parser.add_argument("timing_results", help="Output of /usr/bin/time -v [command]")
    parser.add_argument("-o", "--output", help="Output file", default="result.json")
    return parser.parse_args()

def construct_migration_graph(labeling, tree):
    migration_graph = nx.DiGraph()
    for label in labeling["label"].unique():
        migration_graph.add_node(label)

    for u, v in tree.edges:
        label_u = labeling.loc[u, "label"]
        label_v = labeling.loc[v, "label"]
        if label_u != label_v:
            if migration_graph.has_edge(label_u, label_v):
                migration_graph[label_u][label_v]["count"] += 1
            else:
                migration_graph.add_edge(label_u, label_v, count=1)

    return migration_graph

"""
Computes the set of all pairwise relations in a directed graph.
"""
def get_relations(G):
    relations = set()
    for (u, v) in G.edges():
        relations.add((u, v))

    for u in G.nodes():
        for v in G.nodes():
            if u == v:
                continue
            if nx.has_path(G, u, v):
                relations.add((u, v))

    return relations

def main():
    args = parse_args()

    tree = nx.read_edgelist(args.tree, nodetype=str, create_using=nx.DiGraph())
    perturbed_tree = nx.read_edgelist(args.perturbed_tree, nodetype=str, create_using=nx.DiGraph())
    true_labeling = pd.read_csv(args.vertex_labeling).set_index("vertex")
    inferred_labeling = pd.read_csv(args.inferred_vertex_labeling).set_index("vertex")

    # parse timing results
    with open(args.timing_results, 'r') as f:
        timing = f.read()
        match = re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (.*)', timing) 
        elapsed_time = match.groups()[0]
        elapsed_time = sum(x * float(t) for x, t in zip([1, 60, 3600], elapsed_time.split(":")[::-1]))

    # construct true and inferred migration graphs
    true_migration_graph = construct_migration_graph(true_labeling, tree)
    inferred_migration_graph = construct_migration_graph(inferred_labeling, perturbed_tree)

    # run maximum parsimony on the true tree using the inferred migration graph
    # to infer the labeling of the vertices
    root = [u for u in tree.nodes if tree.in_degree(u) == 0][0]
    inferred_edges = set([(u, v) for u, v in inferred_migration_graph.edges()])
    character_set = list(true_migration_graph.nodes())

    def dist_f(x, y):
        if x == y:
            return 0
        if (x, y) in inferred_edges:
            return 1
        return np.inf

    def leaf_f(node):
        return true_labeling.loc[node, "label"]
    
    scores = mp(tree, root, character_set, leaf_f, dist_f)
    labeling = {}
    for node in nx.dfs_preorder_nodes(tree, root):
        if node == root:
            labeling[node] = min(character_set, key=lambda x: scores[node][x])
            continue

        parent = list(tree.predecessors(node))[0]
        parent_label = labeling[parent]
        labeling[node] = min(character_set, key=lambda x: scores[node][x] + dist_f(parent_label, x))

    # compute the number of correctly labeled vertices
    num_correctly_labeled = 0
    for vertex in tree.nodes:
        true_label = true_labeling.loc[vertex, "label"]
        inferred_label = labeling[vertex]
        if true_label == inferred_label:
            num_correctly_labeled += 1

    # compute the true parsimony score and the inferred parsimony score
    true_parsimony_score = 0
    for u, v in tree.edges:
        true_label_u = true_labeling.loc[u, "label"]
        true_label_v = true_labeling.loc[v, "label"]
        if true_label_u != true_label_v:
            true_parsimony_score += 1

    inferred_parsimony_score = 0
    for u, v in perturbed_tree.edges:
        inferred_label_u = labeling[u]
        inferred_label_v = labeling[v]
        if inferred_label_u != inferred_label_v:
            inferred_parsimony_score += 1

    # compute the F1 score of the true and inferred relations
    # implied by the true and inferred migration graphs

    true_relations = get_relations(true_migration_graph)
    inferred_relations = get_relations(inferred_migration_graph)

    false_positives = inferred_relations - true_relations
    positives = true_relations
    false_negatives = true_relations - inferred_relations
    negatives = set([(u, v) for u in true_migration_graph.nodes() for v in true_migration_graph.nodes() if u != v]) - positives

    fpr = len(false_positives) / len(negatives)

    if len(negatives) == 0:
        fnr = 0
    else:
        fnr = len(false_negatives) / len(positives)

    # compute Jaccard index
    jaccard_index = len(inferred_relations.intersection(true_relations)) / len(inferred_relations.union(true_relations))

    result = {
        'pairwise_relations': {
            'false_positive_rate': fpr,
            'false_negative_rate': fnr,
            'false_positives': len(false_positives),
            'false_negatives': len(false_negatives),
            'true_positives': len(positives) - len(false_negatives),
            'true_negatives': len(negatives) - len(false_positives),
            'positives': len(positives),
            'negatives': len(negatives),
            'jaccard_index': jaccard_index
        }
    }

    result['migration_graph_num_edges'] = true_migration_graph.number_of_edges()
    result['migration_graph_num_vertices'] = true_migration_graph.number_of_nodes()
    result['num_correctly_labeled'] = num_correctly_labeled
    result['num_vertices'] = len(tree.nodes)
    result['true_parsimony_score'] = true_parsimony_score
    result['inferred_parsimony_score'] = inferred_parsimony_score
    result['elapsed_time'] = elapsed_time

    with open(args.output, "w") as f:
        json.dump(result, f, indent=4)

if __name__ == "__main__":
    main()
