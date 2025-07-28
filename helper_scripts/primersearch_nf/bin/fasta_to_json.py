#!/usr/bin/env python

import argparse
import json
from pathlib import Path
from Bio import SeqIO
import hashlib


def hash_sequence(sequence: str) -> str:
    md5 = hashlib.md5(sequence.encode("utf-8"))
    max_bits_in_result = 56
    p = (1 << max_bits_in_result) - 1
    rest = int(md5.hexdigest(), 16)
    result = 0
    while rest != 0:
        result ^= (rest & p)
        rest >>= max_bits_in_result
    return str(result)

def load_primer_ids(primer_file):
    primer_ids = set()
    with open(primer_file) as f:
        for line in f:
            fields = line.strip().split()
            if fields:
                primer_ids.add(fields[0])
    return primer_ids

def extract_primer_id(seq_id, primer_ids):
    # Return the one and only matching primer ID in the sequence ID
    matches = [pid for pid in primer_ids if pid in seq_id]
    if len(matches) == 1:
        return matches[0]
    return None  # Skip ambiguous or missing matches

def process_fasta_file(fasta_file, primer_ids, sample_id):
    values = {}

    for record in SeqIO.parse(fasta_file, "fasta"):
        primer_key = extract_primer_id(record.id, primer_ids)
        if primer_key:
            h = hash_sequence(str(record.seq))
            values[primer_key] = h
        else:
            print(f"Warning: Skipping record with ambiguous/missing primer match: {record.id}")
            pass

    result = {
        "sample_id": sample_id,
        "values": values
    }

    with open(f"{sample_id}.json", "w") as out:
        json.dump(result, out, indent=4)

def main():
    parser = argparse.ArgumentParser(description="Convert hashed FASTA sequences to JSON using primer info.")
    parser.add_argument("--fasta_file", required=True, help="fasta file name")
    parser.add_argument("--primers", required=True, help="Primer file (use only first column)")
    parser.add_argument("--sample_id", required=True, help="sample id name")
    args = parser.parse_args()

    primer_ids = load_primer_ids(args.primers)
    process_fasta_file(args.fasta_file, primer_ids, args.sample_id)


if __name__ == "__main__":
    main()
