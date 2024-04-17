import os 
import json

import pandas as pd
import networkx as nx

results_dir = 'results/laml/'
results     = [d for d in os.listdir(results_dir) if os.path.isdir(os.path.join(results_dir, d))]
settings    = [('unweighted', 'none'), ('weighted', 'none'), 
               ('unweighted', 'polyclonal_tree'), ('weighted', 'polyclonal_tree'),
               ('unweighted', 'polyclonal_dag'), ('weighted', 'polyclonal_dag')]

rows = []
for result in results:
    row = {'experiment': result}

    # load polyMACH results
    row['polymach'] = []
    for weighted, migration_pattern in settings:
        migration_graph_file = os.path.join('results/fast_machina/', result, weighted, migration_pattern, 'inferred_migration_graph.csv')

        if not os.path.exists(migration_graph_file):
            raise ValueError(f'{migration_graph_file} does not exist')
            continue

        migration_graph = pd.read_csv(migration_graph_file)
        G = nx.DiGraph()
        for i, r in migration_graph.iterrows():
            G.add_edge(r['source'], r['target'], count=r['count'])
        parsimony = int(migration_graph['count'].sum())

        results_file = os.path.join('results/fast_machina', result, weighted, migration_pattern, 'inferred_results.json')
        with open(results_file, 'r') as f:
            results = json.load(f)
            log_likelihood = results['objective']

        row['polymach'].append({})
        row['polymach'][-1]['migration_graph'] = list(G.edges())
        row['polymach'][-1]['weighted'] = weighted
        row['polymach'][-1]['migration_pattern'] = migration_pattern
        row['polymach'][-1]['parsimony'] = parsimony
        row['polymach'][-1]['objective'] = log_likelihood

    # load Fitch results
    fitch_file = os.path.join('results/fitch', result, 'results.json')
    with open(fitch_file, 'r') as f:
        fitch_results = json.load(f)
        row['fitch'] = fitch_results
        del row['fitch']['topologies']

    # load parsimonous relabeling results
    parsimonous_relabeling_file = os.path.join('results/parsimonious_relabeling', result, 'inferred_results.json')
    with open(parsimonous_relabeling_file, 'r') as f:
        parsimonous_relabeling_results = json.load(f)
        row['parsimonous_relabeling'] = parsimonous_relabeling_results

    # load LAML results
    row['laml'] = {}
    seed_log_likelihood_file = os.path.join('results/laml', result, 'seed_params.txt')
    with open(seed_log_likelihood_file, 'r') as f:
        text = f.read()
        lines = text.split('\n')

        row['laml']['seed_tree_dropout_rate'] = float(lines[0].split(': ')[1])
        row['laml']['seed_tree_silencing_rate'] = float(lines[1].split(': ')[1])
        row['laml']['seed_tree_mutation_rate'] = float(lines[3].split(': ')[1])

        seed_log_likelihood = None
        for line in lines:
            if 'Negative-llh' in line:
                seed_log_likelihood = float(line.split(': (')[1].split(',')[0])
                break

        row['laml']['seed_tree_log_likelihood'] = seed_log_likelihood

    laml_log_likelihood_file = os.path.join('results/laml', result, 'laml_params.txt')
    with open(laml_log_likelihood_file, 'r') as f:
        text = f.read()
        lines = text.split('\n')

        row['laml']['dropout_rate'] = float(lines[0].split(': ')[1])
        row['laml']['silencing_rate'] = float(lines[1].split(': ')[1])
        row['laml']['log_likelihood'] = float(lines[2].split(': ')[1])
        row['laml']['mutation_rate'] = float(lines[3].split(': ')[1])

    leaf_labeling_file = os.path.join('results/laml', result, 'laml_leaf_labeling.csv')
    with open(leaf_labeling_file, 'r') as f:
        lines = f.readlines()
        ncells = len(lines) - 1
        row['ncells'] = ncells

    rows.append(row)

with open('results/analysis.json', 'w') as f:
    json.dump(rows, f, indent=4)

