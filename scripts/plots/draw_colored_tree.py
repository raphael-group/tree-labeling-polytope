import subprocess
import argparse

import networkx as nx
import pandas as pd
import seaborn as sns

from loguru import logger
from Bio import Phylo

"""
Converts a tree in newick format into a directed networkx graph.
"""
def from_newick_get_nx_tree(tree_path):
    phylo_tree = Phylo.read(tree_path, 'newick')
    net_tree = Phylo.to_networkx(phylo_tree)

    node_renaming_mapping = {}
    idx = 0
    for node in net_tree.nodes:
        if str(node) == 'Clade':
            node_renaming_mapping[node] = f'clade_{idx}'
            idx = idx + 1
        else:
            node_renaming_mapping[node] = node.name
    node_renaming_mapping[list(net_tree.nodes)[0]] = 'root'

    net_tree = nx.relabel_nodes(net_tree, node_renaming_mapping)

    directed_tree = nx.DiGraph()
    directed_tree.add_edges_from(list(nx.bfs_edges(net_tree, 'root')))
    return directed_tree

"""
Parses a tree from a file into a directed networkx graph.
"""
def parse_tree(tree_file, format='adjacency_list'):
    if format == 'adjacency_list':
        return nx.read_adjlist(tree_file, create_using=nx.DiGraph())
    elif format == 'edgelist':
        return nx.read_edgelist(tree_file, create_using=nx.DiGraph())
    elif format == 'newick':
        return from_newick_get_nx_tree(tree_file)
    else:
        raise ValueError("Unknown format: {}".format(format))

def rgb_tuple_to_hex(rgb):
    rgb = [int(255 * x) for x in rgb]
    return "#{:02x}{:02x}{:02x}".format(*rgb)

def draw_colored_tree(T, labeling, color_map, f):
    f.write("digraph T {\n")
    for u in T.nodes():
        if u in labeling.index:
            label = labeling.loc[u, 'label']
            color = rgb_tuple_to_hex(color_map[label])
            f.write(f"\t\"{u}\" [label=\"{label}\", fillcolor=\"{color}\", style=filled];\n")
        else:
            f.write(f"\t\"{u}\";\n")
    for u, v in T.edges():
        f.write(f"\t\"{u}\" -> \"{v}\";\n")
    f.write("}\n")

def make_color_graph(T, labeling, color_map):
    color_graph = nx.DiGraph()
    for u in labeling['label'].unique():
        color_graph.add_node(u, color=rgb_tuple_to_hex(color_map[u]))
    for u, v in T.edges():
        u_label = labeling.loc[u, 'label']
        v_label = labeling.loc[v, 'label']
        if u_label != v_label:
            if not color_graph.has_edge(u_label, v_label):
                color_graph.add_edge(u_label, v_label, count=1)
            else:
                color_graph[u_label][v_label]['count'] += 1

    return color_graph

def draw_color_graph(G, f, multi_edges=False):
    f.write("digraph G {\n")
    for u in G.nodes():
        f.write(f"\t\"{u}\" [fillcolor=\"{G.nodes[u]['color']}\", style=filled];\n")
    for u, v in G.edges():
        if multi_edges:
            count = G[u][v]['count']
            f.write(f"\t\"{u}\" -> \"{v}\" [label={count}];\n")
        else:
            f.write(f"\t\"{u}\" -> \"{v}\";\n")
    f.write("}\n")

def parse_args():
    parser = argparse.ArgumentParser(description="Draw a tree along with its labeling in the DOT language.")
    parser.add_argument("tree", help="The input tree")
    parser.add_argument("labeling", help="The input labeling")
    parser.add_argument("-o", "--output", help="The output prefix", default="result")
    parser.add_argument("-m", "--multi-edges", help="Display multi-edges", action="store_true", default=False)
    parser.add_argument("-s", "--svg", help="Output as SVG", action="store_true", default=False)

    parser.add_argument(
        "-f", "--format", help="The format of the input trees", 
        choices=['adjacency_list', 'edgelist', 'newick'], 
        default='adjacency_list'
    )

    return parser.parse_args()

def main():
    args = parse_args()
    T = parse_tree(args.tree, args.format)
    labeling = pd.read_csv(args.labeling).set_index('vertex')

    num_colors = len(labeling['label'].unique())
    colors = sns.color_palette("deep", num_colors)
    labels = sorted(labeling['label'].unique())
    color_map = {label: colors[i] for i, label in enumerate(labels)}

    with open(f"{args.output}_colored_tree.dot", "w") as f:
        draw_colored_tree(T, labeling, color_map, f)

    with open(f"{args.output}_color_graph.dot", "w") as f:
        G = make_color_graph(T, labeling, color_map)
        draw_color_graph(G, f, args.multi_edges)

    if args.svg:
        logger.info(f"Converting {args.output}_colored_tree.dot to {args.output}_colored_tree.svg...")
        subprocess.run(["dot", "-Tsvg", f"{args.output}_colored_tree.dot", "-o", f"{args.output}_colored_tree.svg"])

        logger.info(f"Converting {args.output}_color_graph.dot to {args.output}_color_graph.svg...")
        subprocess.run(["dot", "-Tsvg", f"{args.output}_color_graph.dot", "-o", f"{args.output}_color_graph.svg"])

if __name__ == "__main__":
    main()
