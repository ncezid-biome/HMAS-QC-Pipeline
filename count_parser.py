#!/usr/bin/env python3

import argparse, os
import pandas as pd
import numpy as np
import datetime
from pathlib import Path
import run_blast as blast
pd.set_option('display.max_columns', None)

# helper function to move specific column to the 'pos' position
def move_column_inplace(df, col, pos):
    col = df.pop(col)
    df.insert(pos, col.name, col)

def main():

    # note
    # 1. the full format count table file is required !
    # 2. assume the sample/primer_group pair (as column name) in count table is in format of sample.primer_group
    #
    # the script can be run as: python3 count_parser.py -c count_table -s sample_file -b blast_file -o output.txt
    # or: python3 count_parser.py -c count_table -s sample_file -o output.txt -f final.fasta -r designer.fasta
    # (if we don't alredy have the blast result)
    parser = argparse.ArgumentParser(description = 'parse the final count table file in full format')
    parser.add_argument('-c', '--count_file', metavar = '', required = True, help = 'Specify count table')
    parser.add_argument('-s', '--sample_file', metavar = '', required = True, help = 'Specify sample list file')
    parser.add_argument('-o', '--output', metavar = '', required = True, help = 'Specify output file')
    parser.add_argument('-b', '--blast', metavar = '', required = False, help = 'Specify blast result file')
    parser.add_argument('-f', '--fasta', metavar = '', required = False, help = 'Specify fasta file output from HMAS QC pipeline (should have "final" in the filename).')
    parser.add_argument('-r', '--reference', metavar = '', required = False, help = 'Specify fasta file containing the positive control targets.')
    args = parser.parse_args()


    if Path(f"{args.count_file}.pkl").is_file():
        print(f"...loading {args.count_file}.pkl")
        df = pd.read_pickle(f"{args.count_file}.pkl")
    else:
        # read count_table
        print("loading csv count file")
        df = pd.read_csv(args.count_file, sep='\t')
        print("done loading csv count file")
        pd.to_pickle(df, f"{args.count_file}.pkl")
        print(datetime.datetime.now())   
        print("done dumping pkl file")

    # read sample list file
    sample_list = pd.read_csv(args.sample_file, names = ['sample'])['sample'].tolist()
    sample_list.sort(key=str.lower)

    # read blast result file
    bcolnames = ["seq", "primer", "query_len", "subj_len", "eval", "cov", "pident", "mismatch"]
    if args.blast:
        blast_df = pd.read_csv(args.blast, sep='\t', names=bcolnames)
    else:
        blast_file = blast.blast(args.fasta, args.reference, os.path.basename(args.fasta), 20)
        blast_df = pd.read_csv(blast_file, sep='\t', names=bcolnames)
    blast_df['primer'] = blast_df['primer'].str.strip().str.split('-').str[1]


    #1. filter columns to only contain sample names in sample_list
    df.columns = df.columns.str.strip() # remove potential spaces
    seq_df = df.iloc[:,0] # keep the 1st column (sequence)
    primer_list = [i for i in df.columns if i.split('.')[0] in sample_list]
    df = df[primer_list]
    df['seq'] = seq_df

    #1.1 melt the df 
    df = df.melt(id_vars=['seq'],value_vars=primer_list,var_name='sample_primer',value_name='count')
    df['primer'] = df['sample_primer'].str.strip().str.split('.').str[1]

    #1.2 merge 2 dfs (count table and blast result)
    df = pd.merge(df,blast_df,on=['seq','primer'])

    #1.3 filter by pident == 100 & cov >= 90
    df = df[(df['pident'] == 100) & (df['cov'] >= 90)]
    df.drop_duplicates(subset=['seq','sample_primer'], inplace=True)

    #1.4 reverse the melt, and return to original format of df
    df = df.drop(columns = ["query_len", "subj_len", "eval", "cov", "pident", "mismatch", "primer"])
    df = df.pivot(index='seq',columns='sample_primer',values='count')

    #2. make the new column name list 
    n_column = list(set([i.split('.')[1] for i in df.columns]))
    n_column.sort(key=str.lower)

    #3. count all non-zero cells for each column (as a new series)
    count_df = df.fillna(0).astype(bool).sum(axis=0)

    #4. construct lists of count values(as columns) for the final df
    n_column_list = []
    cdf_idx = count_df.index
    for pp in n_column:
        # check if sample.pp is in the index, otherwise make it blank
    	n_column_list.append([count_df[f"{sample}.{pp}"] if f"{sample}.{pp}" in cdf_idx else '' for sample in sample_list])

    #5. create the final_df from the list, transpose and reset the index name
    final_df = pd.DataFrame(n_column_list, index=n_column).T
    final_df.index = sample_list
    print (final_df.columns[:5])
    print (final_df.shape)

    #5.1 count the fail/pass, add extra columns to the front
    new_df = final_df.replace(r'^\s*$', np.nan, regex=True)
    new_df = new_df.fillna(0)

    # #5.1.1 add a count row to the last
    # #count the number of miss for each primers
    # last_row = list(new_df[new_df == 0].count(axis=0).values)
    # final_df.loc[len(final_df.index)] = last_row
    # sample_list.append('# of primer miss') #index for the last row
    # final_df.index = sample_list
    # # res = [idx for idx, val in enumerate(last_row) if val and val > 8] #perfect val is 7

    final_df['fail'] = new_df[new_df == 0].count(axis=1)
    final_df['pass_below_5'] = new_df[(new_df <= 5) & (new_df > 0)].count(axis=1)
    final_df['pass_over_5'] = new_df[new_df > 5].count(axis=1)

    move_column_inplace(final_df,'fail',0)
    move_column_inplace(final_df,'pass_below_5',1)
    move_column_inplace(final_df,'pass_over_5',2)

    final_df.to_csv(args.output, sep='\t')

if __name__ == "__main__":
    main()