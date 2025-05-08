#!/usr/bin/env python

'''
this script will traverse a given folder with json files, which are in the format of the following:

{
    "sample_id": "SRR13367354_assembled",
    "values": {
        "OG0002019primerGroup8": "16514681765883615",
        "OG0002293primerGroup9": "45366157905267940",
        "OG0000978primerGroup4": "24251287982917091",
    }
}
OG0002019primerGroup8 in our case is the primer name, 16514681765883615 is the hashed string of
the amplicon sequence for this primer on the sample SRR13367354

It will then build a pairwise difference matrix with these information. Rows/columns being the sample IDs,
and the cell is the count of difference for all common amplicon sequences. 
(i.e., compare the hash value of all the common primers between any pair of samples)
'''
import os
import json
import argparse
import pandas as pd
import time
from itertools import combinations
from concurrent.futures import ProcessPoolExecutor, as_completed

# Global shared sample_dict to avoid serialization overhead
sample_dict = {}

def init_worker(shared_data):
    global sample_dict
    sample_dict = shared_data

def collect_json_files(input_dir):
    json_files = []
    for root, _, files in os.walk(input_dir):
        for f in files:
            if f.endswith('.json'):
                json_files.append(os.path.join(root, f))
    return json_files

def load_samples(json_files):
    data = {}
    for file in json_files:
        with open(file, 'r') as f:
            d = json.load(f)
            sample_id = d["sample_id"]
            data[sample_id] = d["values"]
    return data

def pairwise_diff_indices(pair):
    s1, s2 = pair
    s1_values = sample_dict[s1]
    s2_values = sample_dict[s2]
    shared_keys = set(s1_values.keys()) & set(s2_values.keys())
    diff_count = sum(1 for k in shared_keys if s1_values[k] != s2_values[k])
    return s1, s2, diff_count

def compute_diff_matrix_parallel(sample_dict_full, max_workers=None):
    sample_ids = sorted(sample_dict_full.keys())
    diff_matrix = pd.DataFrame(index=sample_ids, columns=sample_ids, dtype=int)

    pairs = list(combinations(sample_ids, 2))

    print(f"Submitting {len(pairs):,} pairwise jobs to {max_workers or os.cpu_count()} cores...")

    with ProcessPoolExecutor(max_workers=max_workers, initializer=init_worker, initargs=(sample_dict_full,)) as executor:
        futures = [executor.submit(pairwise_diff_indices, pair) for pair in pairs]
        for i, future in enumerate(as_completed(futures), 1):
            s1, s2, diff = future.result()
            diff_matrix.loc[s1, s2] = diff
            diff_matrix.loc[s2, s1] = diff
            if i % 10000 == 0:
                print(f"Processed {i:,} / {len(pairs):,} pairs")

    for s in sample_ids:
        diff_matrix.loc[s, s] = 0

    return diff_matrix

def main():
    parser = argparse.ArgumentParser(description="Compute pairwise differences from JSON files.")
    parser.add_argument("input_dir", help="Directory containing .json files")
    parser.add_argument("output_file", help="CSV file to save the pairwise matrix")
    parser.add_argument("--workers", type=int, default=None, help="Number of CPU cores to use (default: all available)")
    args = parser.parse_args()

    json_files = collect_json_files(args.input_dir)
    print(f"Found {len(json_files)} .json files.")

    all_samples = load_samples(json_files)
    print(f"Loaded {len(all_samples)} samples.")

    start = time.time()
    matrix = compute_diff_matrix_parallel(all_samples, max_workers=args.workers)
    elapsed = time.time() - start
    print(f"Matrix computation took {elapsed:.2f} seconds.")

    matrix.to_csv(args.output_file)
    print(f"Matrix saved to {args.output_file}")

if __name__ == "__main__":
    main()
