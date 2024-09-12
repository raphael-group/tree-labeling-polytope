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
                    match = re.search(r'n(\d+)_m(\d+)_r(\d+)_s(\d+)', d)
                    n, m, r, s = match.groups()
                   
                    log_file = os.path.join(directory, d, f'log.txt')
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

                    row = {
                        'n': int(n),
                        'm': int(m),
                        'r': int(r),
                        's': int(s),
                        'time': time,
                        'obj': obj,
                        'root_relax': root_relax,
                    }

                    data.append(row)

    return pd.DataFrame(data)

def main():
    parser = argparse.ArgumentParser(description='Load JSON files into a pandas DataFrame.')
    parser.add_argument('directory', type=str, help='The directory containing the JSON files.')
    args = parser.parse_args()

    df = load_files(args.directory)
    df.to_csv('softwired_sp_summary.csv', index=False)

if __name__ == '__main__':
    main()

