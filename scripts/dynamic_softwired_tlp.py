import argparse
import itertools
import json
import sys

import pandas as pd
import numpy as np
import networkx as nx
import gurobipy as gp
 
from loguru import logger
from tqdm import tqdm
from enum import Enum
from collections import defaultdict
from Bio import Phylo, SeqIO
from gurobipy import GRB

"""
Builds the softwired TLP model for a single character.
"""
def build_softwired_tlp_base(N, root, sequences, alphabet):
    # grab reticulation nodes and edges
    reticulation_edges = set([(u, v) for u, v in N.edges if N.in_degree(v) > 1])
    reticulation_nodes = set([n for n in N.nodes if N.in_degree(n) > 1])

    model = gp.Model("softwired_tlp")
    decisions = model.addVars(N.edges, alphabet, alphabet, name="decisions", vtype=GRB.CONTINUOUS)
    reticulations = {}

    # ensure all non-reticulation edges are included    
    for u, v in N.edges:
        if (u, v) in reticulation_edges:
            reticulations[(u, v)] = model.addVar(vtype=GRB.BINARY)
        else:
            reticulations[(u, v)] = model.addVar(vtype=GRB.CONTINUOUS, lb=1, ub=1)

    for v in reticulation_nodes:
        model.addConstr(gp.quicksum(reticulations[(u, v)] for u in N.predecessors(v)) == 1)

    # add flow constraints
    for u, v in N.edges:
        for i in alphabet:
            for w in N.successors(v):
                lhs = gp.quicksum(decisions[u, v, j, i] for j in alphabet)
                rhs = gp.quicksum(decisions[v, w, i, j] for j in alphabet)
                model.addConstr(lhs - rhs <= 2 - reticulations[(u, v)] - reticulations[(v, w)])
                model.addConstr(lhs - rhs >= -2 + reticulations[(u, v)] + reticulations[(v, w)])

        model.addConstr(gp.quicksum(decisions[u, v, i, j] for i in alphabet for j in alphabet) == reticulations[(u, v)])

    weight_f = lambda i, j: 1 if i != j else 0
    model.setObjective(gp.quicksum(weight_f(i, j) * decisions[u, v, i, j] for u, v in N.edges for i in alphabet for j in alphabet), GRB.MINIMIZE)
    return model, decisions, reticulations
    
def parse_args():
    parser = argparse.ArgumentParser(description="Dynamic softwired small parsimony using the TLP")
    parser.add_argument("tree", help="Phylogenetic tree in Newick format")
    parser.add_argument("sequences", help="Amino acid sequences in FASTA format")
    return parser.parse_args()

def load_sequences(sequences_file):
    return {record.id: str(record.seq) for record in SeqIO.parse(sequences_file, "fasta")}

def load_tree(tree_file):
    tree = Phylo.read(args.tree, 'newick')
    tree.rooted = True

    nx_tree = Phylo.to_networkx(tree)
    node_relabeling_map = {}
    internal_node_counter = 0
    for node in nx_tree.nodes:
        if nx_tree.out_degree(node) == 0:
            node_relabeling_map[node] = node.name
        else:
            node_relabeling_map[node] = f"internal_{internal_node_counter}"
            internal_node_counter += 1

    return nx.relabel_nodes(nx_tree, node_relabeling_map)

if __name__ == "__main__":
    args = parse_args()

    # Load tree
    N = load_tree(args.tree)

    # Add a dummy root node to the tree
    root = [n for n in N.nodes if N.in_degree(n) == 0][0]
    N.add_node("dummy_root")
    N.add_edge("dummy_root", root)

    # Load sequences
    sequences       = load_sequences(args.sequences)
    sequence_length = len(next(iter(sequences.values())))
    alphabet        = set(itertools.chain(*sequences.values()))

    # Build the base model for a single character
    model, decisions, reticulations = build_softwired_tlp_base(N, "dummy_root", sequences, alphabet)

    for leaf, _ in sequences.items():
        v = leaf.split('_')[0] # 
        u = next(N.predecessors(v))
        for i in alphabet:
            for j in alphabet:
                decisions[u,v,i,j].ub = 0

    score = 0
    scores = {}
    for character_idx in range(sequence_length):
        if len(set(sequence[character_idx] for sequence in sequences.values())) == 1:
            continue

        for (leaf, sequence) in sequences.items():
            c = sequence[character_idx]
            v = leaf.split('_')[0] # 
            u = next(N.predecessors(v))
            j = sequence[character_idx]

            for i in alphabet:
                decisions[u,v,i,j].ub = 1
        
        model.optimize()
        score += model.objVal
        scores[character_idx] = model.objVal

        for (leaf, sequence) in sequences.items():
            c = sequence[character_idx]
            v = leaf.split('_')[0] # 
            u = next(N.predecessors(v))
            j = sequence[character_idx]

            for i in alphabet:
                decisions[u,v,i,j].ub = 0

    with open("scores.json", "w") as f:
        json.dump(scores, f)



