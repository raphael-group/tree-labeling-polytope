import sys
import pandas as pd
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Intersect labels with vertex names.')
    parser.add_argument('--labels', type=str, required=True, help='Labels file')
    parser.add_argument('--edgelist', type=str, required=True, help='Edge list file')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    labels = pd.read_csv(args.labels, sep='\t', header=None)
    edge_list = pd.read_csv(args.edgelist, sep='\t', header=None)
    edge_list.columns = ['src', 'dst']
    vertex_names = set(edge_list['src'].unique()) | set(edge_list['dst'].unique())
    labels = labels[labels[0].isin(vertex_names)]
    # to stdout
    labels.to_csv(sys.stdout, sep='\t', header=False, index=False)
