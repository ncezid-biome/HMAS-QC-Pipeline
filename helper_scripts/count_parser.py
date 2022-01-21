#!/usr/bin/env python3

import argparse, os
import pandas as pd
import numpy as np
import datetime
from pathlib import Path
import run_blast as blast
import count_plot
from operator import truediv
# pd.set_option('display.max_columns', None)


control_list = ['2013K_0676',
                '2013K_1246',
                '2014K_0979',
                '2016K_0878',
                'Blank_IRP1',
                'Blank_IRP2',
                'Water1']

# helper function to move specific column to the 'pos' position
def move_column_inplace(df, col, pos):
    col = df.pop(col)
    df.insert(pos, col.name, col)


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
    bcolnames = ["seq", "primer", "query_len", "subj_len", "eval", "cov", "pident", "mismatch"]

    if not blast_file: #if we don't already have blast result as a text file
        blast_file = blast.blast(query_fasta, reference, os.path.basename(query_fasta), max_hit)
        
    blast_df = pd.read_csv(blast_file, sep='\t', names=bcolnames)
    blast_df['sample'] = blast_df['primer'].str.strip().str.split('-').str[2] + '.'
    blast_df['primer'] = blast_df['primer'].str.strip().str.split('-').str[1]
    # blast_df['sample_primer'] = blast_df[['sample','primer']].agg('.'.join, axis=1)
    blast_df['sample_primer'] = blast_df['sample'] + blast_df['primer']

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
    df = df.melt(id_vars=['seq'],value_vars=primer_list,var_name='sample_primer',value_name='count')
    # df['primer'] = df['sample_primer'].str.strip().str.split('.').str[1]
    print (df.head(n=5))
    print (df.shape)

    #1.2 merge 2 dfs (count table and blast result)
    # df = pd.merge(df,blast_df,on=['seq','primer'])
    df = pd.merge(df,blast_df,on=['seq','sample_primer'])
    print (df.head(n=5))
    print (df.shape)

    #1.3 filter by pident == 100 & cov >= 90
    df = df[(df['pident'] == 100) & (df['cov'] >= 90)]
    df.drop_duplicates(subset=['seq','sample_primer'], inplace=True)
    print (df.head(n=5))
    print (df.shape)

    #1.4 reverse the melt, and return to original format of df
    df = df.drop(columns = ["query_len", "subj_len", "eval", "cov", "pident", "mismatch", "primer", "sample"])
    df = df.pivot(index='seq',columns='sample_primer',values='count')
    print (df.head(n=5))
    print (df.shape)

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

# helper method to:
# 1. find the most abundant unique sequence for each sample-primer
# 2. based on primer/seqs, search the blast result (also requiring 'pident' == 100, 'cov' >=90)
# 3. based on the match result, construct a binary table for sample-primer
# 4. re-organize the dataframe, and plot the histogram
def plot_perfect_match(df, blast_df, sample_list):

    # helper method for checking passed in dateframe row value has perfect match in blast result
    def find_perfect_match(row):

        # primer = row['primer'].strip().split('.')[1]
        sample_primer = row['primer'].strip()
        # match = blast_df[(blast_df['primer'] == primer) & (blast_df['seq'] == row['seq'])]
        match = blast_df[(blast_df['sample_primer'] == sample_primer) & (blast_df['seq'] == row['seq'])]
        if match.empty:
            return 0
        else:
            match2 = match[(match['pident'] == 100) & (match['cov'] >= 90)]
            if match2.empty:
                return 0
            else:
                return 1

    #1. find most abundant seqs for each sample-primer pair
    df.set_index('seq', inplace=True)
    df = df.idxmax() #series of index of max value for each column
    df = df.to_frame(name='seq').reset_index() #convert to df and save the original indx column(primer)
    df.rename(columns={'index':'primer'}, inplace=True)
    df['match'] = df[['seq','primer']].apply(find_perfect_match, axis=1)

    pd.to_pickle(df, f"find_perfect_match.pkl")
    # df = pd.read_pickle(f"find_perfect_match.pkl")


    #2. make the new column name list 
    n_column = list(set([i.split('.')[1] for i in df['primer']]))
    n_column.sort(key=str.lower)
    print (len(n_column))

    # create a Series to use with split_samle_primer
    sr = pd.Series(df['match'].values, index=df['primer'])

    final_df = split_samle_primer(sr, n_column, sample_list)

    count_plot.plot_hist(final_df,'number of samples where the most abundant hit was a perfect match')


def create_hit_table(df, sample_list, out_file):
    '''
    this method takes a filtered(by blast result) count table and convert from oringal form I to form III.
    I.
    Rep_Seq Sam1_PP1    Sam2_PP1    Sam1_PP2    Sam2_PP2
    Seq1    10  42  0   0
    Seq2    0   0   86  0
    Seq3    4   0   0   0

    II. any seq count >0 will be counted 1 (and summed) for that sam_pp pair
    Sam1_PP1    2
    Sam2_PP1    1
    Sam1_PP2    1
    Sam2_PP2    0

    III.
        PP1 PP2
    Sam1    2   1
    Sam2    1   0

    Parameters
    ----------
    df: the original dataframe
    out_file: the output file name
    sample_list: the list of samples we're studying

    '''

    #1. make the new column name list 
    n_column = list(set([i.split('.')[1] for i in df.columns]))
    n_column.sort(key=str.lower)

    #2. count all non-zero cells for each column (as a new series)
    count_df = df.fillna(0).astype(bool).sum(axis=0)
    print (count_df.head(n=5))
    print (count_df.shape)

    final_df = split_samle_primer(count_df, n_column, sample_list)
    print (final_df.head(n=5))
    print (final_df.shape)

    final_df['fail'] = final_df[final_df == 0].count(axis=1)
    final_df['pass_below_5'] = final_df[(final_df <= 5) & (final_df > 0)].count(axis=1)
    final_df['pass_over_5'] = final_df[final_df > 5].count(axis=1)

    move_column_inplace(final_df,'fail',0)
    move_column_inplace(final_df,'pass_below_5',1)
    move_column_inplace(final_df,'pass_over_5',2)

    final_df.to_csv(out_file, sep='\t')



def create_abundance_table(df):

    #5.1 count the fail/pass, add extra columns to the front
    final_df = df.replace(r'^\s*$', np.nan, regex=True)
    final_df = final_df.fillna(0)

    # # #5.1.1 add a count row to the last
    # # #count the number of miss for each primers
    # # last_row = list(new_df[new_df == 0].count(axis=0).values)
    # # final_df.loc[len(final_df.index)] = last_row
    # # sample_list.append('# of primer miss') #index for the last row
    # # final_df.index = sample_list
    # # # res = [idx for idx, val in enumerate(last_row) if val and val > 8] #perfect val is 7

    final_df['below_10'] = final_df[final_df < 10].count(axis=1)
    final_df['between_10_100'] = final_df[(final_df <= 100) & (final_df >= 10)].count(axis=1)
    final_df['between_100_200'] = final_df[(final_df <= 200) & (final_df > 100)].count(axis=1)
    final_df['over_200'] = final_df[final_df > 200].count(axis=1)

    move_column_inplace(final_df,'below_10',0)
    move_column_inplace(final_df,'between_10_100',1)
    move_column_inplace(final_df,'between_100_200',2)
    move_column_inplace(final_df,'over_200',3)

    return final_df


# helper method to split the column names, sort it
def get_column_list(df):

    column_list = list(set([i.split('.')[1] for i in df.index]))
    column_list.sort(key=str.lower)
    return column_list


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

    df = read_count_table(args.count_file)

    # read sample list file
    sample_list = pd.read_csv(args.sample_file, names = ['sample'])['sample'].tolist()
    sample_list.sort(key=str.lower)
    print (sample_list)

    blast_df = create_blast_df(args.blast, args.fasta, args.reference, 20)
    print (blast_df.head(n=5))
    print (blast_df.shape)

    #1. filter columns to only contain sample names in sample_list
    df.columns = df.columns.str.strip() # remove potential spaces
    seq_df = df.iloc[:,0] # keep the 1st column (sequence)
    primer_list = [i for i in df.columns if i.split('.')[0] in sample_list]
    df = df[primer_list]
    df['seq'] = seq_df

    old_df = df.sum() #total abundanceall high quality seqs

    # # df_copy = df.copy()
    # # plot_perfect_match(df_copy, blast_df, sample_list) #plot the most abundant hit was a perfect match

    # filter out sequences so that seqs left all have matches in blast result and their pident == 100 & cov >= 90
    df = merge_count_blast(df, primer_list, blast_df)
    print (df.head(n=5))
    print (df.shape)
    create_hit_table(df, sample_list, args.output)
    # # print (df.head(n=5))
    # # print (df.shape)

    '''
    new_df = df.sum()
    pd.to_pickle(new_df, f"perfect_match_sum_outputbase_newblast_dilution.pkl")
    # new_df = pd.read_pickle(f"perfect_match_sum_outputbase_newblast.pkl")
    
    # # perfect_match_abundance_list = list(new_df.values)
    # # print (len(perfect_match_abundance_list))
    # # # print (perfect_match_abundance_list)


    new_df = split_samle_primer(new_df, get_column_list(new_df), sample_list)

    # #plot all perfect matches abundance / total abundanceall high quality seqs
    old_df.drop(['seq'], inplace=True)
    old_df = split_samle_primer(old_df, get_column_list(old_df), sample_list)
    abundance_div_df = new_df.div(old_df, fill_value = 0)
    # count_plot.plot_hist_abundance(abundance_div_df, 'mean perfect match abundance / total abundance', 'primer count', \
    #                                 'perfect match abundance ratio -- cutadapt trimming M347-21-026')
    abundance_div_df = abundance_div_df.filter(items=[idx for idx in new_df.index if idx not in control_list], axis=0)
    bin_list = list(abundance_div_df.mean(axis=0).values)
    count_plot.plot_hist_abundance(bin_list, 20, 'abundance_ratio_cutadapt_trimming_026_new_dilution.pdf', \
                                    'mean perfect match abundance / total abundance', 'primer count', \
                                    'perfect match abundance ratio -- cutadapt trimming M347-21-026')


    filtered_df = new_df.filter(items=[idx for idx in new_df.index if idx not in control_list], axis=0)
    bin_list = list(filtered_df.mean(axis=0).values)
    # the plot parameters likely need to be customized !
    # count_plot.plot_hist_abundance(new_df, 'mean abundance value', 'primer count', \
    #                                 'abundance histogram cutadapt trimming M347-21-026')
    count_plot.plot_hist_abundance(bin_list, 10, 'abundance_bin_cutadapt_trimming_026_new_dilution.pdf', \
                                    'mean abundance value', 'primer count', \
                                    'abundance histogram cutadapt trimming M347-21-026')



    # the plot parameters need to be customized !
    # i.e. y-axis value, range, labels, etc...
    count_plot.plot_length_abundance(filtered_df, list(filtered_df.mean(axis=0).values), \
                                    [175, 255, 0, 400], 'abundance_length_dist_cutadapt_trimming_026_new_dilution.pdf', 'mean abundance value')

    count_plot.plot_length_abundance(filtered_df, list(map(truediv, list(filtered_df.std(axis=0).values), list(filtered_df.mean(axis=0).values))), \
                                    [175, 255, 0, 4], 'CV_length_dist_cutadapt_trimming_026_new_dilution.pdf', 'CV (std / mean abundance value)')

    # # final_df = create_abundance_table(new_df)
    # # final_df.to_csv(args.output, sep='\t')

    '''
if __name__ == "__main__":
    main()