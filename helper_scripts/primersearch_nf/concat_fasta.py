#!/usr/bin/env python

# helper script to concatenate .fasta file (assume it's the only fasta file in the subfolder)
# from the subfolders of 2 given folders
# and saved in the 3rd folder with the same folder structure and fasta file name

import os
import sys

if len(sys.argv) != 4:
    print(f"Usage: {sys.argv[0]} <folder1> <folder2> <output_folder>")
    sys.exit(1)

folder1 = sys.argv[1]
folder2 = sys.argv[2]
output_folder = sys.argv[3]

# Create the output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Loop through all subfolders in folder1
for subfolder in os.listdir(folder1):
    path1 = os.path.join(folder1, subfolder)
    path2 = os.path.join(folder2, subfolder)

    # Check if matching subfolder exists in folder2
    if os.path.isdir(path1) and os.path.isdir(path2):
        # Find the .fasta file in each subfolder
        fasta_files1 = [f for f in os.listdir(path1) if f.endswith('.fasta')]
        fasta_files2 = [f for f in os.listdir(path2) if f.endswith('.fasta')]

        if len(fasta_files1) != 1 or len(fasta_files2) != 1:
            print(f"Warning: Unexpected number of .fasta files in {path1} or {path2}")
            continue

        fasta1 = fasta_files1[0]
        fasta2 = fasta_files2[0]

        # Ensure filenames match
        if fasta1 != fasta2:
            print(f"Warning: .fasta filenames do not match in {subfolder}: {fasta1} vs {fasta2}")
            continue

        file1_path = os.path.join(path1, fasta1)
        file2_path = os.path.join(path2, fasta2)

        # Create output subfolder
        output_subfolder = os.path.join(output_folder, subfolder)
        os.makedirs(output_subfolder, exist_ok=True)

        output_file_path = os.path.join(output_subfolder, fasta1)  # Keep the same filename

        # Concatenate the two fasta files
        with open(output_file_path, 'wb') as outfile:
            for f in [file1_path, file2_path]:
                with open(f, 'rb') as infile:
                    outfile.write(infile.read())

        print(f"Concatenated {file1_path} + {file2_path} -> {output_file_path}")
