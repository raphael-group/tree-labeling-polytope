import pandas as pd
import os
import argparse
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description='Converts FitchCount character matrices to Startle format.')
    parser.add_argument('character_matrix', type=str, help='Input file')
    parser.add_argument('output', type=str, help='Output prefix')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()

    character_matrix = pd.read_csv(args.character_matrix, sep='\t', index_col=0)
    character_matrix = character_matrix.replace('-', -1).astype(int)

    # get list of character,state pairs
    character_states = defaultdict(list)
    for character in character_matrix.columns:
        states = character_matrix[character].unique()
        states = states[states != -1]
        character_states[character] = states
    character_states = dict(character_states)

    leaf_labeling = character_matrix.index.str.split('.').str[0]
    leaf_labeling = list(zip(character_matrix.index, leaf_labeling))

    with open(args.output + '_leaf_labeling.csv', 'w') as f:
        f.write('leaf,label\n')
        for leaf, label in leaf_labeling:
            f.write(leaf + ',' + label + '\n')

    with open(args.output + '_mutation_priors.csv', 'w') as f:
        f.write('character,state,probability\n')
        for character, states in character_states.items():
            probability = 1 / len(states)
            for state in states:
                f.write(f"{character},{state},{probability}\n")

    character_matrix.to_csv(args.output + '_character_matrix.csv')

    
    deduplicated_character_matrix_rows = []
    unique_rows = set()
    pruning_map = {}
    for row in character_matrix.itertuples():
        original_row = row
        row = tuple(row)
        cell = row[0]
        row = row[1:]
        if row not in unique_rows:
            pruning_map[row] = {"leader": cell, "members": []}
            deduplicated_character_matrix_rows.append(original_row)
            unique_rows.add(row)
        else:
            pruning_map[row]["members"].append(cell)

    deduplicated_character_matrix = pd.DataFrame(
        deduplicated_character_matrix_rows, 
    )

    print(pruning_map)
    print(deduplicated_character_matrix)
