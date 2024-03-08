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
                    algorithm = d[:match.start() - 1]
                    n, m, s, e, setting, _, _ = match.groups()

                    false_positive_rate = content['pairwise_relations']['false_positive_rate']
                    false_negative_rate = content['pairwise_relations']['false_negative_rate']
                    positives = content['pairwise_relations']['positives']
                    negatives = content['pairwise_relations']['negatives']
                    false_positives = content['pairwise_relations']['false_positives']
                    false_negatives = content['pairwise_relations']['false_negatives']
                    jaccard_index = content['pairwise_relations']['jaccard_index']
                    edgeset_symmetric_difference = content['pairwise_relations']['edgeset_symmetric_difference']
                    multi_edgeset_symmetric_difference = content['pairwise_relations']['multi_edgeset_symmetric_difference']
                    
                    row = {
                        'algorithm': algorithm,
                        'cells': n,
                        'migration_rate': float(m),
                        'seed': s,
                        'error': e,
                        'setting': setting,
                        'elapsed_time': content['elapsed_time'],
                        'num_correctly_labeled': content['num_correctly_labeled'],
                        'num_vertices': content['num_vertices'],
                        'true_parsimony_score': content['true_parsimony_score'],
                        'inferred_parsimony_score': content['inferred_parsimony_score'],
                        'false_positive_rate': false_positive_rate,
                        'false_negative_rate': false_negative_rate,
                        'positives': positives,
                        'negatives': negatives,
                        'false_positives': false_positives,
                        'false_negatives': false_negatives,
                        'jaccard_index': jaccard_index,
                        'edgeset_symmetric_difference': edgeset_symmetric_difference,
                        'multi_edgeset_symmetric_difference': multi_edgeset_symmetric_difference
                    }

                    data.append(row)

    return pd.DataFrame(data)

def main():
    parser = argparse.ArgumentParser(description='Load JSON files into a pandas DataFrame.')
    parser.add_argument('directory', type=str, help='The directory containing the JSON files.')
    args = parser.parse_args()

    df = load_files(args.directory)
    df.to_csv('summary.csv', index=False)

if __name__ == '__main__':
    main()

