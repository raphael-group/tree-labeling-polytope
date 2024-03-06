import re
import os
import json
import argparse
import pandas as pd
import networkx as nx

ROOT_DIR = "/n/fs/ragr-research/projects/tree-labeling-polytope/nextflow_results/ground_truth/"

settings = ['polyclonal_tree', 'polyclonal_dag', 'monoclonal_tree', 'monoclonal_dag', 'none']
def load_files(directory):
    data = []
    for d in os.listdir(directory):
        for file in os.listdir(os.path.join(directory, d)):
            if file.endswith(".json"):
                with open(os.path.join(directory, d, file), 'r') as f:
                    content = json.load(f)

                    match = re.search(r'n(\d+)_m(1E-\d+|0.\d+)_s(\d+)_e(\d+)_(' + '|'.join(settings) + ')', d)
                    n, m, s, e, setting = match.groups()
                    content['n'] = n
                    content['m'] = m
                    content['s'] = s
                    content['e'] = e
                    content['setting'] = setting

                    labeling = pd.read_csv(f"{ROOT_DIR}/n{n}_m{m}_s{s}_e{e}_{setting}/sim_labeling.csv").set_index('vertex')
                    tree = nx.read_edgelist(
                        f"{ROOT_DIR}/n{n}_m{m}_s{s}_e{e}_{setting}/sim_tree_edgelist.tsv", 
                        create_using=nx.DiGraph
                    )

                    migration_graph = nx.DiGraph()
                    for u, v, d in tree.edges(data=True):
                        lu = labeling.loc[u, 'label']
                        lv = labeling.loc[v, 'label']
                        if lu != lv and not migration_graph.has_edge(lu, lv):
                            migration_graph.add_edge(lu, lv)

                    content['migration_graph'] = list(migration_graph.edges())

                    if nx.is_tree(migration_graph):
                        content['migration_graph_type'] = 'polyclonal_tree'
                    elif nx.is_directed_acyclic_graph(migration_graph):
                        content['migration_graph_type'] = 'polyclonal_dag'
                    else:
                        content['migration_graph_type'] = 'none'

                    print(content)
                    data.append(content)

    return pd.DataFrame(data)

def main():
    parser = argparse.ArgumentParser(description='Load JSON files into a pandas DataFrame.')
    parser.add_argument('directory', type=str, help='The directory containing the JSON files.')
    args = parser.parse_args()

    df = load_files(args.directory)
    df.to_csv('combinatorial_summary.csv', index=False)

if __name__ == '__main__':
    main()

