import argparse
import networkx as nx
import pandas as pd
import seaborn as sns

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

def draw_colored_tree(T, labeling, color_map):
    print("digraph T {")
    for u in T.nodes():
        if u in labeling.index:
            label = labeling.loc[u, 'label']
            color = rgb_tuple_to_hex(color_map[label])
            print(f"  {u} [fillcolor=\"{color}\", style=filled];")
        else:
            print(f"  {u};")
    for u, v in T.edges():
        print(f"  {u} -> {v};")
    print("}")


def make_color_graph(T, labeling, color_map):
    color_graph = nx.DiGraph()
    for u in labeling['label'].unique():
        color_graph.add_node(u, color=rgb_tuple_to_hex(color_map[u]))
    for u, v in T.edges():
        u_label = labeling.loc[u, 'label']
        v_label = labeling.loc[v, 'label']
        if u_label != v_label:
            if not color_graph.has_edge(u_label, v_label):
                color_graph.add_edge(u_label, v_label)
    return color_graph

def draw_color_graph(G):
    print("digraph G {")
    for u in G.nodes():
        print(f"  {u} [fillcolor=\"{G.nodes[u]['color']}\", style=filled];")
    for u, v in G.edges():
        print(f"  {u} -> {v};")
    print("}")

def parse_args():
    parser = argparse.ArgumentParser(description="Draw a tree along with its labeling in the DOT language.")
    parser.add_argument("tree", help="The input tree")
    parser.add_argument("labeling", help="The input labeling")

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
    colors = sns.color_palette("husl", num_colors)
    color_map = {label: colors[i] for i, label in enumerate(labeling['label'].unique())}

    draw_colored_tree(T, labeling, color_map)
    # G = make_color_graph(T, labeling, color_map)
    # draw_color_graph(G)

if __name__ == "__main__":
    main()
