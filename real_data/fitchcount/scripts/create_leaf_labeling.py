import pandas as pd
import networkx as nx
import sys
import argparse

parser = argparse.ArgumentParser(description='Create leaf labeling')
parser.add_argument('tree', help='Tree in TSV format')
parser.add_argument('labeling', help='Cell labeling in TSV format')
args = parser.parse_args()

tree = nx.read_edgelist(args.tree, delimiter='\t', create_using=nx.DiGraph())
vertices = set(tree.nodes)
labeling = pd.read_csv(args.labeling, sep=',')
labeling = labeling[labeling.vertex.isin(vertices)]
labeling.columns = ['leaf', 'label']
labeling['label'] = labeling['label'].str.split('_').str[-1]

labeling.to_csv(sys.stdout, sep=',', index=False)
