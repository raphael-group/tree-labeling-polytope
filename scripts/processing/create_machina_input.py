import argparse
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="Create input files for MACHINA.")
    parser.add_argument("vertex_labeling", help="Vertex labeling in CSV format")
    parser.add_argument("-o", "--output", help="Output prefix", default="machina")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    df = pd.read_csv(args.vertex_labeling).set_index("vertex")
    print(df.loc['root', 'label'])
    labels = df['label'].unique()
    with open(f"{args.output}_coloring.txt", "w") as f:
        for i, label in enumerate(labels):
            f.write(f"{label} {i}\n")
