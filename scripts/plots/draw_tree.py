import argparse
import matplotlib.pyplot as plt
import networkx as nx

from Bio import Phylo
from networkx.drawing.nx_agraph import graphviz_layout

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

def parse_tree(tree_file, format='adjacency_list'):
    if format == 'adjacency_list':
        return nx.read_adjlist(tree_file, create_using=nx.DiGraph())
    elif format == 'edge_list':
        return nx.read_edgelist(tree_file, create_using=nx.DiGraph())
    elif format == 'newick':
        return from_newick_get_nx_tree(tree_file)
    else:
        raise ValueError("Unknown format: {}".format(format))

def draw_graph(T, ax=None):
    pos = graphviz_layout(T, prog='dot')

    labels = {}
    for node in T.nodes():
        labels[node] = node 

    nx.draw(T, pos, with_labels=False, arrows=True, ax=ax)
    # nx.draw_networkx_labels(T, pos, labels, font_size=10, ax=ax)

def main():
    parser = argparse.ArgumentParser(description="Draw a graph from an adjacency list")
    parser.add_argument("trees", help="Files containing the input trees", nargs='+')

    parser.add_argument(
        "-f", "--format", help="The format of the input trees", 
        choices=['adjacency_list', 'edge_list', 'newick'], 
        default='adjacency_list'
    )

    args = parser.parse_args()

    fig, axes = plt.subplots(nrows=1, ncols=len(args.trees))
    if len(args.trees) == 1:
        axes = [axes]

    for i, tree in enumerate(args.trees):
        axes[i].set_title(tree)
        T = parse_tree(tree, args.format)
        draw_graph(T, axes[i])

    plt.show()

if __name__ == "__main__":
    main()
