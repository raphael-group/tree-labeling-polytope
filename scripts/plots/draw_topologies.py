import json
import subprocess
import seaborn as sns
import networkx as nx
import argparse

def rgb_tuple_to_hex(rgb):
    rgb = [int(255 * x) for x in rgb]
    return "#{:02x}{:02x}{:02x}".format(*rgb)

def draw_color_graph(G, f, multi_edges=False):
    f.write("digraph G {\n")
    for u in G.nodes():
        f.write(f"\t\"{u}\" [fillcolor=\"{G.nodes[u]['color']}\", style=filled];\n")
    for u, v in G.edges():
        if False:
            count = G[u][v]['count']
            f.write(f"\t\"{u}\" -> \"{v}\" [label={count}];\n")
        else:
            f.write(f"\t\"{u}\" -> \"{v}\";\n")
    f.write("}\n")

def parse_args():
    parser = argparse.ArgumentParser(description="Draws a set of topologies in the DOT language.")
    parser.add_argument("topologies", help="JSON file containing graph topologies")
    parser.add_argument("-o", "--output", help="The output prefix", default="result")
    parser.add_argument("-s", "--svg", help="Output as SVG", action="store_true", default=False)
    parser.add_argument("-p", "--palette", help="The palette to use", default=None)

    return parser.parse_args()

def main():
    args = parse_args()

    with open(args.topologies, 'r') as f:
        topologies = json.load(f)['topologies']

    labels = [edge['src'] for G in topologies for edge in G['edges']] + [edge['dst'] for G in topologies for edge in G['edges']]
    labels = sorted(list(set(labels)))
    num_colors = len(labels)
    colors = sns.color_palette("deep", num_colors)
    color_map = {label: colors[i] for i, label in enumerate(labels)}
    convert_to_hex = True

    for i, edgeset in enumerate(topologies):
        G = nx.MultiDiGraph()
        for edge in edgeset['edges']:
            src = edge['src']
            dst = edge['dst']
            for _ in range(edge['multiplicity']):
                G.add_edge(src, dst)

        for u in G.nodes():
            G.nodes[u]['color'] = rgb_tuple_to_hex(color_map[u]) if convert_to_hex else color_map[u]

        with open(f"{args.output}_c{edgeset['count']}_{i}.dot", 'w') as f:
            draw_color_graph(G, f, multi_edges=True)

        if args.svg:
            subprocess.run(["dot", "-Tsvg", f"{args.output}_c{edgeset['count']}_{i}.dot", "-o", f"{args.output}_c{edgeset['count']}_{i}.svg"])

if __name__ == "__main__":
    main()


