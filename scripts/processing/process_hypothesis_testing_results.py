import os
import glob
import re
import pandas as pd
import networkx as nx

base_path = 'nextflow_results/hypothesis_testing/fast_machina/'

pattern = r'n(\d+)_m([0-9.]+)_s(\d+)_(.+)_vs_(.+)'
results = []

for file_path in glob.glob(os.path.join(base_path, '**', 'inferred_migration_graph.csv'), recursive=True):
    dir_name = os.path.basename(os.path.dirname(file_path))
    match = re.search(pattern, dir_name)
    if match:
        # Extract the matched groups
        n, m, s, setting_1, setting_2 = match.groups()

        df = pd.read_csv(file_path)
        parsimony = df['count'].sum()

        true_migration_graph = f'nextflow_results/hypothesis_testing/ground_truth/n{n}_m{m}_s{s}_{setting_1}_vs_{setting_2}/sim_migration_graph.csv'
        df = pd.read_csv(true_migration_graph)
        G = nx.from_pandas_edgelist(df, 'src', 'dst', create_using=nx.DiGraph())
        print(G)

        # Store the parsed information
        results.append({
            'n': int(n),
            'm': float(m),
            's': int(s),
            'setting_1': setting_1,
            'setting_2': setting_2,
            'parsimony': parsimony,
            'edges': G.number_of_edges(),
            'nodes': G.number_of_nodes(),
            'is_dag': nx.is_directed_acyclic_graph(G),
            'is_tree': nx.is_tree(G),
            'G_str': G.edges()
        })

df = pd.DataFrame(results)
df = df.sort_values(by=['n', 'm', 's', 'setting_1'])
df.to_csv('hypothesis_testing_results.csv', index=False)
