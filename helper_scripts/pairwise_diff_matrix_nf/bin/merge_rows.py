#!/usr/bin/env python

import argparse
import pandas as pd

def parse_arguments():
    parser = argparse.ArgumentParser(description='Merge row CSVs into a full matrix.')
    parser.add_argument('-i', '--input', nargs='+', required=True, help='List of input row CSV files')
    parser.add_argument('-o', '--output', required=True, help='Output CSV file for full matrix')
    return parser.parse_args()

def clean_sample_name(name):
    return name.replace('_assembled_extractedAmplicons', '')

def main():
    args = parse_arguments()
    row_files = sorted(args.input)

    # Read first file to initialize
    full_matrix = pd.read_csv(row_files[0], header=None, index_col=0)
    full_matrix.index = full_matrix.index.map(clean_sample_name)
    # full_matrix.columns = full_matrix.columns.map(clean_sample_name)

    for f in row_files[1:]:
        df = pd.read_csv(f, header=None, index_col=0)
        df.index = df.index.map(clean_sample_name)
        # df.columns = df.columns.map(clean_sample_name)
        full_matrix = pd.concat([full_matrix, df], axis=0)

    full_matrix.columns = full_matrix.index
    full_matrix.index.name = None
    full_matrix.to_csv(args.output, index=True, header=True)

if __name__ == '__main__':
    main()

