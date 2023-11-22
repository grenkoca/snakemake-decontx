import anndata
import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description='Dump observations from anndata object')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input anndata object')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output file')
    parser.add_argument('--new_columns', type=str, nargs='+')
    parser.add_argument('--new_values', type=str, nargs='+')

    return parser.parse_args()


def main():
    args = parse_args()
    df = pd.read_csv(args.input, sep='\t', header=None)

    for col, val in list(zip(args.new_columns, args.new_values)): 
        df[0] = df[0] + "-" + val

    df.to_csv(args.output, sep='\t', index=False)

if __name__ == '__main__':
    main()

