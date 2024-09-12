import re
import os
import json
import argparse
import pandas as pd

def load_files(directory):
    data = []
    for d in os.listdir(directory):
        for file in os.listdir(os.path.join(directory, d)):
            if file.endswith(".json"):
                with open(os.path.join(directory, d, file), 'r') as f:
                    content = json.load(f)
                    match = re.search(r'n(\d+)_m(1E-\d+|0.\d+)_e(\d+)_s(\d+)_((poly|mono)clonal_(tree|dag))', d)
                    n, m, e, s, setting, _, _ = match.groups()
                   
                    log_file = os.path.join(directory, d, f'timing.txt')
                    with open(log_file, 'r') as log:
                        log_content = log.read()

                    match = re.search(r'Explored \d+ nodes \(\d+ simplex iterations\) in (\d+\.\d+) seconds', log_content)
                    time = match.groups()[0]

                    match = re.search(r'Best objective (\d+\.\d+e\+\d+), best bound (\d+\.\d+e\+\d+), gap (\d+\.\d+)%', log_content)
                    obj = match.groups()[0]

                    match = re.search(r'Root relaxation: objective ([^,]*), \d+ iterations, .*', log_content)
                    if match:
                        root_relax = match.groups()[0]
                    else:
                        root_relax = obj


                    migration_graph_file = os.path.join(directory, d, f'inferred_migration_graph.csv')
                    migration_graph = pd.read_csv(migration_graph_file)

                    row = {
                        'n': int(n),
                        'm': float(m),
                        'e': int(e),
                        's': int(s),
                        'setting': setting,
                        'time': time,
                        'obj': obj,
                        'root_relax': root_relax,
                        'num_colors': migration_graph['count'].sum(),
                    }

                    data.append(row)

    return pd.DataFrame(data)

def main():
    parser = argparse.ArgumentParser(description='Load JSON files into a pandas DataFrame.')
    parser.add_argument('directory', type=str, help='The directory containing the JSON files.')
    args = parser.parse_args()

    df = load_files(args.directory)
    df.to_csv('convex_recoloring_summary.csv', index=False)

if __name__ == '__main__':
    main()

