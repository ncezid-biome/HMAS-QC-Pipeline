#!/usr/bin/env python

import argparse, os, csv
import random, linecache, re
import pandas as pd
from tabulate import tabulate
from Bio.Seq import Seq
from numpy.core.defchararray import find

def file_len(filename):
    """
    Returns the number of lines in a file

    Params
    ------
    filename: String
        Name of the file to query

    Returns
    ------
    Integer
        Number of lines in the file
    """
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def append_value(dict_obj, key, value):
    """
    Appends a value to a dictionary in list format

    Params
    ------
    dict_obj: String
        Name of the dictionary
    key: String
        The key to search for
    value: String
        The value to add to the key's list of values

    Returns
    ------
    none

    """
    if key in dict_obj:
        # Check if type of value of key is list or not
        if not isinstance(dict_obj[key], list):
            # If type is not list then make it list
            dict_obj[key] = [dict_obj[key]]
        dict_obj[key].append(value)
    else:
        dict_obj[key] = value

def main():

    parser = argparse.ArgumentParser(description = 'Check the output of Mothur pcr.seqs command.')
    parser.add_argument('-s', '--size', metavar = '', required = True, \
                        help = 'Specify sample size')
    parser.add_argument('-b', '--before', metavar = '', required = True, \
                        help = 'Specify fasta file that is input to pcr.seqs')
    parser.add_argument('-a', '--after', metavar = '', required = True, \
                        help = 'Specify fasta file output from pcr.seqs')
    parser.add_argument('-l', '--oligos', metavar = '', required = True, \
                        help = 'Specify oligos file')
    parser.add_argument('-g', '--group', metavar = '', required = True, \
                        help = 'Specify group file')
    parser.add_argument('-o', '--outfile', metavar = '', required = True, \
                        help = 'Specify path and name for output file')

    args = parser.parse_args()

    # Generate a random sample of integers.
    # Range of the random sample equals the number of reads in the pcr.seqs output fasta file
    random_sample = random.sample(range(round(file_len(args.after)/2)), int(args.size))

    # We need to get the reads based on our random sample.
    # Then we need to get the sequence info for these IDs (sequence after pcr.seqs and before it)
    # Then we need to merge this info with the primer sequence info, and print this table

    d = {} # Dictionary to hold random sample of read IDs and their sequences
    for i in random_sample:
        line = linecache.getline(args.after, i) # Get the line corresponding to the random integer
        # If line is a READ ID, make that the key. If not, make the next line the key.
        if line.startswith('>'):
            key = line.rstrip().split()[0][1:]
            key = key.split('|')[0]
            seq = linecache.getline(args.after, i+1)
        else:
            key = linecache.getline(args.after, i+1).split()[0][1:]
            key = key.split('|')[0]
            seq = linecache.getline(args.after, i+2)
        d[key] = seq.strip() # Key =  read ID, value = sequence after pcr.seqs cmd was run
        
    # For the read IDs in dict, get the sequence before pcr.seqs command was run
    with open(args.before, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if line.strip().split()[0][1:] in d.keys():
                    append_value(d, line.strip().split()[0][1:], next(f).strip())

    # For the read IDs in dict, get the primer that Mothur has identified in them 
    with open(args.group, 'r') as f:
        for line in f:
            if line.strip().split()[0] in d.keys():
                primer = line.strip().split()[1]
                append_value(d, line.strip().split()[0], primer.split('.')[1])


    o = {} # Dictionary to hold oligos file info (primer seqs)
    with open(args.oligos, 'r') as f:
        for line in f:
            if line.startswith('primer'):
                key = line.strip().split()[3] # primer name
                fwd = line.strip().split()[1] # forward primer seq
                rev = line.strip().split()[2] # reverse primer seq
                rev_seq = Seq(rev) 
                rev_compl = str(rev_seq.reverse_complement()) # Reverse complement of rev primer
                fwd_seq = Seq(fwd)
                fwd_compl = str(fwd_seq.reverse_complement()) # Reverse complement of fwd primer
                o[key] = fwd # Key is primer name, 1st value is forward primer seq
                append_value(o, key, fwd_compl) # 2nd value: rev compl fwd primer
                append_value(o, key, rev) # 3rd value: reverse primer
                append_value(o, key, rev_compl) # 4th value: rev compl rev primer

    reads = pd.DataFrame([(k, *v) for k, v in d.items()])
    reads.columns = ('id', 'after', 'before', 'primer_name')
    
    #  
    primers = pd.DataFrame([(k, *v) for k, v in o.items()])
    primers.columns=('primer_name', 'fwd', 'fwd_rc', 'rev', 'rev_rc')

    # Merge the primer sequence into the sample of read IDs
    full = pd.merge(reads, primers, on = 'primer_name', how = 'left')
    
    # Get the fragment with primers trimmed off
    a = full.fwd.values.astype(str)
    b = full.before.values.astype(str)
    full = full.assign(start_fwd=find(b, a)) # Start of fwd primer sequence
    full['fwd_len'] = full['fwd'].str.len() # Length of primer seq [start the fragment at start+len]
    full['start'] = full['start_fwd'] + full['fwd_len'] # This is where the primer-less fragment starts
    a = full.rev.values.astype(str)
    full = full.assign(end=find(b, a)) # Start of rev_compl primer sequence [end fragment 'up to' end]
    full['segment'] = full.apply(lambda x: x[2][x[10]:x[11]],axis=1) # x[2] is the 'before' seq, x[10] is index after fwd primer ends, x[11] is index where rev primer begins
    full.loc[full['start_fwd'] == -1, 'segment'] = "Forward primer not found"  # Indicate if the forward primer not found

    new_full= full[['id', 'before', 'fwd', 'rev_rc', 'after', 'segment', 'primer_name', 'fwd_rc', 'rev', 'start_fwd', 'start', 'end']]

    new_full.to_csv(args.outfile, header=True, index=False, sep='\t')
 
if __name__ == '__main__':
    main()


