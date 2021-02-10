#!/usr/bin/env python

import argparse, os, csv



def count_ids(filename, read_id_tag):
    """Counts lines in a file that start with a specific string
    Params
    ------
    filename: String
        Name of the file
    read_id_tag: String
        String identifying start of line to count

    Returns
    ------
    count: Int
        Count of lines that start with specified string
    """

    with open(filename) as f:
        count = 0
        for line in f:
            if line.startswith(read_id_tag):
                count += 1
    return count

def main():

    parser = argparse.ArgumentParser(description = 'Make the index file matching the trimmed read file.')
    parser.add_argument('-r1', '--read1', metavar = '', required = True, \
                        help = 'Specify trimmed read1 file')
    parser.add_argument('-i1', '--index1', metavar = '', required = True, \
                        help = 'Specify original index file')

    args = parser.parse_args()

    # Name of the new index file
    new_index1 = os.path.splitext(args.read1)[0] + '_I1.fastq'

    # Complete checks for if file exists
    #if os.path.isfile(new_index1) == False:

    # Create a dictionary to hold reads from index file
    d = {}
    with open(args.index1, 'r') as f:       
        for line in f:
            if line.startswith('@M00347'):
                key = line.rstrip()
                vals = []
                for i in range(3):
                    vals.append(next(f).rstrip())
                d[key] = vals # key: read ID, val: the next 3 lines - barcode, +, quality

    # Right the new index file based on the reads in R1
    with open(args.read1, 'r') as f, open(new_index1, 'a') as out:
        cw = csv.writer(out, delimiter = '\n')
        for line in f.readlines():
            if line.startswith('@M00347'):
                try:
                    header = [line.rstrip()] # header for the new index file
                    header.extend(d[line.rstrip()]) # the next 3 lines of the new index file
                    cw.writerow(rows) 
                except:
                    pass

    # Need to check for existence of new_index1 (implement some checks for completion)
    count_r1 = count_ids(args.read1, '@M00347')
    count_idx1 = count_ids(new_index1, '@M00347')   

    if count_r1 == count_idx1:
        print(f'The new index file is correct, and contains {count_idx1} read IDs')
    else:
        print(f'The trimmed read file contains {count_r1} read IDs, but the new index file contains {count_idx1} read IDs')


if __name__ == '__main__':
    main()
 
