#!/usr/bin/env python

import numpy as np
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from Bio import Phylo
import argparse
import pandas as pd

extension = "upgma_tree.newick"

def parse_arguments():
    parser = argparse.ArgumentParser(description='build tree with pairwisde diff matrix.')
    parser.add_argument('-i', '--input', required=True, help='List of input row CSV files')
    parser.add_argument('-o', '--output', required=True, help='output tree name, always with extension upgma_tree.newick')
    return parser.parse_args()

args = parse_arguments()

# Load CSV
df = pd.read_csv(args.input, index_col=0)

# Clean labels
labels = [str(label).strip() for label in df.index.tolist()]

# Distance values
dist_matrix = df.values

# Force symmetry and 0 diagonal
dist_matrix = (dist_matrix + dist_matrix.T) / 2
np.fill_diagonal(dist_matrix, 0)

# Check
assert len(labels) == dist_matrix.shape[0] == dist_matrix.shape[1]

# Build correct lower triangle
matrix = []
for i in range(len(labels)):
    row = []
    for j in range(i+1):   # only j < i
        row.append(dist_matrix[i][j])
    matrix.append(row)

for i, row in enumerate(matrix):
    if i < 5:
        print(f"Row {i}: {row}")
        
# # Build DistanceMatrix
dm = DistanceMatrix(names=labels, matrix=matrix)


# Build UPGMA tree
constructor = DistanceTreeConstructor()
tree = constructor.upgma(dm)

# Save
Phylo.write(tree, f"{args.output}.{extension}", "newick")
