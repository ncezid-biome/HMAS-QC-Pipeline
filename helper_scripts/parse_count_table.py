#!/usr/bin/env python3

import argparse, os
import pandas as pd
import numpy as np
import datetime
from pathlib import Path
import run_blast as blast

# a list of control sample names
control_list = ['2013K_0676',
                '2013K_1246',
                '2014K_0979',
                '2016K_0878',
                'Blank_IRP1',
                'Blank_IRP2',
                'Water1']

# helper method to read a full format count table, and/or save it as a pickle file for faster retrieval
def read_count_table(count_file):

    if Path(f"{count_file}.pkl").is_file():
        print(f"...loading {count_file}.pkl")
        df = pd.read_pickle(f"{count_file}.pkl")
    else:
        # read count_table
        print("loading csv count file")
        df = pd.read_csv(count_file, sep='\t')
        print("done loading csv count file")
        pd.to_pickle(df, f"{count_file}.pkl") 
        print("done dumping pkl file")

    return df

def create_blast_df(blast_file, query_fasta, reference, max_hit):

    # the 'primer' is like: OG0000890-OG0000890primerGroup9-2014K_0979
    bcolnames = ["seq", "primer", "query_len", "subj_len", "len_aln", "eval", "cov", "pident", "mismatch"]

    if not blast_file: #if we don't already have blast result as a text file
        blast_file = blast.blast(query_fasta, reference, os.path.basename(query_fasta), max_hit)
        
    # blast_df = pd.read_csv(blast_file, sep='\t', names=bcolnames, index_col=False, header=None)
    blast_df = pd.read_csv(blast_file, sep='\t', names=bcolnames)
    blast_df['sample'] = blast_df['primer'].str.split('-').str[2]
    blast_df['primer'] = blast_df['primer'].str.split('-').str[1]
    blast_df['sample_primer'] = blast_df['sample'] + '.' + blast_df['primer']

    return blast_df


def merge_count_blast(df, primer_list, blast_df):
    '''
    this method filters the original full-format count table (below), so that seqs 
    all have matches in blast result and their pident == 100 & cov >= 90

    Rep_Seq Sam1_PP1    Sam2_PP1    Sam1_PP2    Sam2_PP2
    Seq1    10  42  0   0
    Seq2    0   0   86  0
    Seq3    4   0   0   0

    Parameters
    ----------
    df: the original count table dataframe
    primer_list: the list of targeted sample.primer ['Sam1_PP1', 'Sam1_PP2'] etc.
    blast_df: blast result dataframe

    Returns the filtered dataframe

    '''

    #1.1 melt the df 
    df = df.melt(id_vars=['seq'],value_vars=primer_list,var_name='sample_primer_orig',value_name='count')

    # if it's a dilution sample, we change it from 2011K_0052_C.OG0002271primerGroup6
    # to 2011K_0052.OG0002271primerGroup6 so it can match blast result
    # df['sample_primer'] = [sp.split('.')[0][:-2] + "." + sp.split('.')[1] \
    #                         if sp.split('.')[0][:-2] in dilution_samples_list else sp for sp in df['sample_primer_orig']]
    df['sample_primer'] = df['sample_primer_orig'] # skip all the dilution data

    #1.2 merge 2 dfs (count table and blast result)
    df = pd.merge(df,blast_df,on=['seq','sample_primer'])

    #1.3 filter by pident == 100 & cov >= 90
    df = df[(df['pident'] == 100) & (df['cov'] >= 90)]
    df.drop_duplicates(subset=['seq','sample_primer_orig'], inplace=True)

    #1.4 reverse the melt, and return to original format of df
    df = df.drop(columns = ["query_len", "subj_len", "eval", "cov", "pident", "mismatch", "primer", "sample", "sample_primer"])
    df = df.pivot(index='seq',columns='sample_primer_orig',values='count')

    return df


def split_samle_primer(sr, primers, sample_list):
    '''
    this method will convert:

    Sam1_PP1    2
    Sam2_PP1    1
    Sam1_PP2    1
    Sam2_PP2    0

    into:

        PP1 PP2
    Sam1    2   1
    Sam2    1   0

    Parameters
    ----------
    sr: the original dataframe Series
    primers: the list of intended column name ['PP1', 'PP2']
    sample_list: the list of samples we're studying

    Returns the converted dataframe

    '''

    #1. construct lists of count values(as columns) for the final df
    n_column_list = []
    _idx = sr.index
    for pp in primers:
        # check if sample.pp is in the index, otherwise make it blank
        n_column_list.append([sr[f"{sample}.{pp}"] if f"{sample}.{pp}" in _idx else 0 for sample in sample_list])

    #2. create the final_df from the list, transpose and reset the index name
    final_df = pd.DataFrame(n_column_list, index=primers).T
    final_df.index = sample_list

    return final_df

def parse_argument():
    # note
    # 1. the full format count table file is required !
    # 2. assume the sample/primer_group pair (as column name) in count table is in format of sample.primer_group
    # 3. sample_file is an one-column csv file containing the names of the sample we're analyzing
    #
    # the script can be run as: python3 parse_count_table.py -c juno.final.full.count_table -s M347-21-026-15sample.csv
    #                               -f juno.final.fasta -r juno_design_primers.fasta
    # - b / - f / - r are optional
    parser = argparse.ArgumentParser(prog = 'pasr_count_table.py')
    parser.add_argument('-c', '--count_file', metavar = '', required = True, help = 'Specify count table')
    parser.add_argument('-s', '--sample_file', metavar = '', required = True, help = 'Specify sample list file')
    parser.add_argument('-b', '--blast', metavar = '', required = False, help = 'Specify blast result file')
    parser.add_argument('-f', '--fasta', metavar = '', required = False, help = 'Specify fasta file output from HMAS QC pipeline (should have "final" in the filename).')
    parser.add_argument('-r', '--reference', metavar = '', required = False, help = 'Specify fasta file containing the positive control targets.')
    
    return parser.parse_args()


# helper method to split the column names, sort it
def get_column_list(df):

    column_list = list(set([i.split('.')[1] for i in df.index]))
    column_list.sort(key=str.lower)
    return column_list


def main():

    args = parse_argument()

    df = read_count_table(args.count_file)

    # read sample list file
    sample_list = pd.read_csv(args.sample_file, names = ['sample'])['sample'].tolist()
    sample_list.sort(key=str.lower)

	#filter columns to contain only sample names in the sample_list
    df.columns = df.columns.str.strip() # remove potential spaces
    seq_df = df.iloc[:,0] # keep the 1st column (sequence)
    primer_list = [i for i in df.columns if i.split('.')[0] in sample_list]
    df = df[primer_list]
    df['seq'] = seq_df

    raw_df = df.sum() #total abundance all high quality seqs
    raw_df.drop(['seq'], inplace=True)
    raw_df = split_samle_primer(raw_df, get_column_list(raw_df), sample_list)
    # exclude those control samples before calculating the mean
    raw_mean_df = raw_df.loc[~raw_df.index.isin(control_list)].mean(axis=1)
    print (raw_mean_df)
    print (raw_mean_df.mean())


    # comment out this block if you don't need blast filtering
    blast_df = create_blast_df(args.blast, args.fasta, args.reference, 20)
    #filter out sequences so that seqs left all have matches in blast result and their pident == 100 & cov >= 90
    df = merge_count_blast(df, primer_list, blast_df)
    new_df = df.sum()
    new_df = split_samle_primer(new_df, get_column_list(new_df), sample_list)
    new_mean_df = new_df.loc[~new_df.index.isin(control_list)].mean(axis=1)
    print (new_mean_df)
    print (new_mean_df.mean())


if __name__ == "__main__":
    main()