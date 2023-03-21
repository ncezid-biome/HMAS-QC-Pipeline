#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import gzip

'''
This script is currently set to run inside the hmas2_sampling_rawreads.nf.
It takes 3 parameters: 
1. a fastq file, with primer sequence (using a subset of our full 2461 HMAS primer set) removed, 
2. a raw reads (fastq.gz) file
3. an output file

The script essentially extracts the raw reads based on the seq_id from the given fastq file

'''
def parse_argument():
    # note
    # the script can be run as:
    # 
    parser = argparse.ArgumentParser(prog = 'extract_rawreads_by_seqid.py')
    parser.add_argument('-i', '--sourcereads', metavar = '', required = True, help = 'Specify source reads file')
    parser.add_argument('-r', '--rawreads', metavar = '', required = True, help = 'Specify rawreads file')
    parser.add_argument('-o', '--output', metavar = '', required = True, help = 'Specify output file')
    
    return parser.parse_args()


def extract_rawreads(source_read, raw_read, output_file, compress=True):
    source_dict = SeqIO.to_dict(SeqIO.parse(source_read, "fastq"))
    with gzip.open(raw_read, 'rt') as r:
        raw_dict = SeqIO.to_dict(SeqIO.parse(r, "fastq"))
        
    #create our sought-after SeqRecord list with matching seq_id
    record_list = [raw_dict[key] for key in source_dict]
    if compress:
        with gzip.open(output_file, 'wt') as w:
            SeqIO.write(record_list, w, "fastq")

if __name__ == "__main__":
    
    args = parse_argument()
    extract_rawreads(args.sourcereads, args.rawreads, args.output)