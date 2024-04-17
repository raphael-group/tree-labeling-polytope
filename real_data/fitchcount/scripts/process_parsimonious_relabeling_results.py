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
                    match = re.search(r'm5k_lg(\d+)', d)
                    exp = match.groups()[0]

                    log_file = os.path.join(directory, d, f'log.txt')
                    with open(log_file, 'r') as log:
                        log_content = log.read()

                    match = re.search(r'Explored \d+ nodes \(\d+ simplex iterations\) in (\d+\.\d+) seconds', log_content)
                    if match is None:
                        row = {
                            'exp': exp,
                            'time': None,
                            'obj': None,
                            'root_relax': None,
                        }
                        data.append(row)
                        continue
                    time = match.groups()[0]

                    match = re.search(r'Best objective (\d+\.\d+e\+\d+), best bound (\d+\.\d+e\+\d+), gap (\d+\.\d+)%', log_content)
                    obj = match.groups()[0]

                    match = re.search(r'Root relaxation: objective ([^,]*), \d+ iterations, .*', log_content)
                    if match:
                        root_relax = match.groups()[0]
                    else:
                        root_relax = obj

                    character_matrix = f'character_matrices/m5k_lg{exp}_character_matrix.alleleThresh.txt'
                    system_command = f'cat {character_matrix} | wc -l'
                    character_matrix_size = int(os.popen(system_command).read().strip())

                    row = {
                        'exp': exp,
                        'time': time,
                        'obj': obj,
                        'root_relax': root_relax,
                        'num_leaves': character_matrix_size - 1
                    }

                    data.append(row)

    return pd.DataFrame(data)

def main():
    parser = argparse.ArgumentParser(description='Load JSON files into a pandas DataFrame.')
    parser.add_argument('directory', type=str, help='The directory containing the JSON files.')
    args = parser.parse_args()

    df = load_files(args.directory)
    df.to_csv('fitchcount_parsimonious_labeling_summary.csv', index=False)

if __name__ == '__main__':
    main()
