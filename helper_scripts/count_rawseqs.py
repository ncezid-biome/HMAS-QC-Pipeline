import argparse
import glob
import os
import gzip

### this script is used to count the number of raw sequences per primer pair (2461 total) in fastg.gz file

# Create an argument parser
parser = argparse.ArgumentParser(description="Count lines in files with specified extension")

# Add the arguments
parser.add_argument("folder_path", help="Path to the folder containing the files")
parser.add_argument("output_file", help="Path to the output file")

# Parse the arguments
args = parser.parse_args()

# Get the folder path and output file path from the arguments
folder_path = args.folder_path
output_file = args.output_file

# Get a list of files with the specified extension
files = glob.glob(folder_path + "/*_R1_001.fastq.gz")

print (len(files))

# Open the output file for writing
with open(output_file, "w") as f:
    # Iterate over each file
    for file in files:
        line_count = 0
        # # Open the file in binary mode
        # with open(file, "rb") as infile:
        # Open the gz-compressed file using gzip.open()
        with gzip.open(file, "rt") as infile:
            # Count the number of lines
            for _ in infile:
                line_count += 1
        # Calculate the result with one decimal precision
        result = round(line_count / 9844, 1) #4 lines per seq in fastq, and total 2461 primers
        # Get the relative file name without preceding path
        relative_file_name = os.path.basename(file)
        # Write the relative file name, line count, and result to the output file, separated by a tab
        f.write("{}\t{}\n".format(relative_file_name, result))