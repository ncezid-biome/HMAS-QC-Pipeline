#!/usr/bin/env python3

import argparse, os, sys
import pandas as pd
from pathlib import Path
import logging
# importing scripts under helper_scripts folder
sys.path.insert(0,r'./helper_scripts')
import run_blast as blast
import utilities
import settings
import run_grep

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.addHandler(settings.log_handler)

# a list of control sample names
# control_list = ['2013K_0676',
#                 '2013K_1246',
#                 '2014K_0979',
#                 '2016K_0878',
#                 'Blank_IRP1',
#                 'Blank_IRP2',
#                 'Water1']

# a list of control samples, which we will exclude from the analysis
control_list = ['2014K_0979',
                'Water1']


def get_oligo_primer_dict():
    '''
    this method generates an instance of Primers class and a dictionary with primer name being the key
    and a default value 0

    Returns a dictionary

    '''
    oligo_primers = utilities.Primers(settings.OLIGO_FILE)
    primer_dict = {key: 0 for key in oligo_primers.pseqs}

    return primer_dict

def cvs_to_dict(cvs_file):
    '''
    This method reads a standard cvs file, uses the 1st column as index to creates a dataframe
    Then it converts the dataframe into a dictionary of dictionaries, with the index as key to 
    the first dictionary and column names of the dataframe as the keys for the 2nd dictionary(s)
    Ex.   this cvs file
    seq_id  primer  sample
    seq1    p1      s1
    seq2    p2      s2
    will turn into this structure:
    {'seq1':{'primer':'p1','sample':'s1'},
     'seq2':{'primer':'p2','sample':'s2'}}
    
    Parameters
    ----------
    cvs_file: the path to the cvs file 

    Returns a dictionary of dictionaries
    '''
    df = pd.read_csv(cvs_file, index_col=0)
    map_dict = df.T.to_dict()

    return map_dict

# this method reads mapping file (sample-to-isolate) in csv format
# and return a dictionary of it 
# key: sample (String)
# value: a list of corresponding isolates 
# (the list has no None value in it, and it's ensured to be a non-empty list)
def map_sample_to_isolate(map_file):

    if map_file:
        df = pd.read_csv(map_file, index_col=0)
        # convert dataframe to a dictionary, index being key, and row being value in the form of a list
        map_dict = df.T.to_dict("list")
        # remove nan from the list
        for key in map_dict:
            map_dict[key] = [item for item in map_dict[key] if not(pd.isnull(item))]
        # remove keys which has empty value
        new_map_dict = {key:val for key,val in map_dict.items() if len(val) > 0}
        return new_map_dict
    

# this method will reverse the above sample-to-isolate dictionary, and 
# will generate a isolate-to-sample dictionary
# key: isolate (String)
# value: a list of corresponding samples
def rev_map_isolate_to_sample(map_dict):
    rev_map_dict = {}
    for key in map_dict:
        isolates = map_dict[key]
        for isolate in isolates:
            if isolate not in rev_map_dict:
                rev_map_dict[isolate] = [key]
            else:
                if key not in rev_map_dict[isolate]:
                    rev_map_dict[isolate].append(key)
    return rev_map_dict


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
    # '>OG0000294-OG0000294primerGroup8-ParatyphiA rc'  remove the ' rc' part
    blast_df['primer'] = blast_df['primer'].str.split().str[0]
    blast_df['sample'] = blast_df['primer'].str.split('-').str[2]
    blast_df['primer'] = blast_df['primer'].str.split('-').str[1]
    blast_df['sample_primer'] = blast_df['sample'] + '.' + blast_df['primer']
    ### TO DO
    # de-couple the formatting 
    # use the cvs_to_dict() fuction to get the dict
    # '>OG0000294-OG0000294primerGroup8-ParatyphiA rc' - the seq_id will be a key of a dictionary
    # blast_df['primer'] = [ dict[seq_id]['primer'] for seq_id in blast_df['primer']]
    # blast_df['sample'] = [ dict[seq_id]['sample'] for seq_id in blast_df['primer']]
    # blast_df['sample_primer'] = blast_df['sample'] + '.' + blast_df['primer']

    return blast_df

def blast_map_sample_to_isolate(blast_df, map_dict):
    '''
    this method helps to convert the 'isolate.primer' column in blast result into its corresponding
    'sample.primer' instead. One isolate might correspond to multiple samples, and this is taken
    care of by explode() method

    Parameters
    ----------
    blast_df: the original blast result df
    map_dict: the dictionary of isolate(key) to samples (a list of all corresponding samples)

    Returns the converted blast result dataframe
    '''
    def map_sample_to_isolate(x):
        sep = '.'
        isolate = x.split(sep)[0]
        primer_pair = x.split(sep)[1]
        sample_list = map_dict[isolate]
        return [f"{sample}{sep}{primer_pair}" for sample in sample_list]

    blast_df['sample_primer'] = blast_df['sample_primer'].apply(map_sample_to_isolate)

    return blast_df.explode('sample_primer', ignore_index=True)


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
    map_dict: mapping dictionary between sample and isolate name 

    Returns the filtered dataframe

    '''

    #1.1 melt the df 
    df = df.melt(id_vars=['seq'],value_vars=primer_list,var_name='sample_primer_orig',value_name='count')
    # change dtype for 'count' to save memory
    df = df.astype({'count':int})
    df['sample_primer'] = df['sample_primer_orig']

    blast_df = blast_df[(blast_df['pident'] >= settings.PIDENT) & \
                        (blast_df['len_aln']/blast_df['subj_len'] >= settings.P_ALIGN) & \
                        (blast_df['cov'] >= settings.PCOV)]
    df = pd.merge(df,blast_df,on=['seq','sample_primer'])

    #1.3
    # duplicates come from the exploded blast_df, and they're all identical in 
    # the subset of ['seq','sample_primer_orig','count'], so we only need to keep one of them
    df.drop_duplicates(subset=['seq','sample_primer_orig'], inplace=True)

    #1.4 reverse the melt, and return to original format of df
    # df = df.drop(columns = ["query_len", "subj_len", "eval", "cov", "pident", "mismatch", "primer", "sample", "sample_primer_x","sample_primer_y"], errors='ignore')
    df = df[['seq','sample_primer_orig','count']]
    df = df.pivot(index='seq',columns='sample_primer_orig',values='count')

    return df


def split_samle_primer(sr, primers, sample_list, raw_idx):
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
    raw_idx: index from the raw df, which is essentially the original columns from count_table

    Returns the converted dataframe

    '''

    #1. construct lists of count values(as columns) for the final df
    n_column_list = []
    _idx = sr.index
    for pp in primers:
        # check if sample.pp is in the index, otherwise make it blank
        # n_column_list.append([sr[f"{sample}.{pp}"] if f"{sample}.{pp}" in _idx else 0 for sample in sample_list])
        count_list = []
        for sample in sample_list:
            if f"{sample}.{pp}" in _idx:
                count_list.append(sr[f"{sample}.{pp}"])
            elif f"{sample}.{pp}" in raw_idx:
                count_list.append(0)
            else:
                count_list.append(None)
        n_column_list.append(count_list)

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
    parser.add_argument('-o', '--output', metavar = '', required = False, help = 'Specify output file')
    parser.add_argument('-b', '--blast', metavar = '', required = False, help = 'Specify blast result file')
    parser.add_argument('-f', '--fasta', metavar = '', required = False, help = 'Specify fasta file output from HMAS QC pipeline (should have "final" in the filename).')
    parser.add_argument('-r', '--reference', metavar = '', required = False, help = 'Specify fasta file containing the positive control targets.')
    parser.add_argument('-m', '--map_file', metavar = '', required = False, help = 'Specify mapping file')

    return parser.parse_args()


# helper method to split the column names, sort it
def get_column_list(df):

    column_list = list(set([i.split('.')[1] for i in df.index]))
    column_list.sort(key=str.lower)
    return column_list


# print out count value for specific seq / sample_primer combo in the count_table
def getCountValue(args):
    # args = parse_argument()
    df = read_count_table(args.count_file)
    df.set_index('Representative_Sequence', inplace=True)

    seqs = ['M03235_53_000000000-K394R_1_2114_24886_10548',
            'M03235_53_000000000-K394R_1_2110_21268_3180',
            'M03235_53_000000000-K394R_1_1103_24702_11865',
            'M03235_53_000000000-K394R_1_2114_11940_17547',
            'M03235_53_000000000-K394R_1_1110_10863_14479',
            'M03235_53_000000000-K394R_1_2113_25175_24045']

    samples = ['2013K_1153',
                '08_0810',
                '2010K_1649',
                '2013K_1153',
                '2015K_0074',
                '2016K_0167']
    
    primer = 'OG0000335primerGroup4'

    for sample in samples:
        # print (sample)
        for seq in seqs:
            # print (seq)
            print (df.loc[seq,f'{sample}.{primer}'])
        print ("\n")

def build_confusion_matrix(sample_list, reference, new_df_t, map_dict):
    '''
    this method flattens a list of sets into one single list and generate the occurrence count for each item
    as a dictionary

    Parameters
    ----------
    sample_list: a list of all samples in the data set
    reference: the reference (design) fasta file for amplicons
    new_df_t: the transformed (and blast filtered) dataframe (row: primers, column: samples)
    map_dict: mapping dictionary between sample and isolate name 

    Returns
    ----------
    df_dict: a dictionary of TP/FP/TN/FN metrics (key: sample_name, value: a list of number of TP/FP/TN/FN counts )
    '''
    df_dict = {}
    for key in sample_list:
        primer_dict = get_oligo_primer_dict()
        if map_dict:
            primer_list = run_grep.get_primers_for_sample_design_fasta(reference, map_dict[key][0])
            # Gut_10_3_1  Typhimurium	ParatyphiA might have 2 isolates
            map_size = len(map_dict[key])
            for i in range(1, map_size):
                primer_list.extend(run_grep.get_primers_for_sample_design_fasta(reference, map_dict[key][i]))
        else:
            primer_list = run_grep.get_primers_for_sample_design_fasta(reference, key)
        primer_list = list(set(primer_dict.keys() & set(primer_list)))
        # set value to 1 (default is 0) for all predicted primers
        primer_dict.update(zip(primer_list,[1]*len(primer_list)))
        pred_df = pd.DataFrame.from_dict(primer_dict, orient='index', columns=['prediction'])

        sr_blast = new_df_t[key]
        merged_df = pred_df.merge(sr_blast, how='left', left_index=True, right_index=True)
        merged_df.rename(columns={key:'blast'}, inplace=True)

        TN = merged_df[(merged_df['prediction'] == 0) & (merged_df['blast'].isna())].shape[0]
        TP = merged_df[(merged_df['prediction'] > 0) & (merged_df['blast'] > 0)].shape[0]
        FN = merged_df[(merged_df['prediction'] > 0) & ((merged_df['blast'] == 0) | (merged_df['blast'].isna()))].shape[0]
        FP = merged_df[(merged_df['prediction'] == 0) & (merged_df['blast'] >= 0)].shape[0]
        df_dict[key] = [TP,FP,FN,TN]

        logger.info(f"sample is: {key}")
        logger.info (f"FP are: {merged_df[(merged_df['prediction'] == 0) & (merged_df['blast'] >= 0)].index.tolist()}")
        logger.info (f"sanity check: {len(primer_dict)}, {len(primer_dict) == (TP + TN + FP + FN)}")
        
    return df_dict

def main():

    args = parse_argument()

    df = read_count_table(args.count_file)

    # read sample list file
    sample_list = pd.read_csv(args.sample_file, names = ['sample'])['sample'].tolist()
    sample_list = [sample for sample in sample_list if sample not in control_list]
    sample_list.sort(key=str.lower)
    print (sample_list)

	#filter columns to contain only sample names in the sample_list
    df.columns = df.columns.str.strip() # remove potential spaces
    seq_df = df.iloc[:,0] # keep the 1st column (sequence)
    primer_list = [i for i in df.columns if i.split('.')[0] in sample_list]
    df = df[primer_list]
    df['seq'] = seq_df
    raw_idx = df.sum().index

    raw_df = df.sum() #total abundance all high quality seqs
    raw_df.drop(['seq'], inplace=True)
    raw_df = split_samle_primer(raw_df, get_column_list(raw_df), sample_list, raw_idx)
    # print (raw_df.mean(axis=1))
    # print (raw_df.mean(axis=1).mean())
    raw_df['mean'] = raw_df.mean(axis=1)
    raw_df.to_csv('raw_df', sep='\t') # save as a tsv file

    #creat the blast df, to blast filtering all sequences
    blast_df = create_blast_df(args.blast, args.fasta, args.reference, 100)
    # mapping dictionary between sample and isolate
    map_dict = map_sample_to_isolate(args.map_file)
    # if the mapping is required
    if map_dict:
        # create a reverse mapping between isolate and samples
        # faster this way
        rev_map_dict = rev_map_isolate_to_sample(map_dict)
        # create a new blast df based on the above mapping
        blast_df = blast_map_sample_to_isolate(blast_df, rev_map_dict)

    # the filtered version of our count table df
    df = merge_count_blast(df, primer_list, blast_df)
    new_df = df.sum()
    new_df = split_samle_primer(new_df, get_column_list(new_df), sample_list, raw_idx)
    # print (new_df.mean(axis=1))
    # print (new_df.mean(axis=1).mean())
    new_df['mean'] = new_df.mean(axis=1)
    new_df.to_csv('new_df', sep='\t') # save as a tsv file

    new_df_t = new_df.T

    df_dict = build_confusion_matrix(sample_list, args.reference, new_df_t, map_dict)
    df = pd.DataFrame.from_dict(df_dict, orient='index')
    df.columns = ['TP','FP','FN','TN']
    df['sensitivity (TP/P)'] = df['TP']/(df['TP'] + df['FN'])
    df['sensitivity (TP/P)'] = df['sensitivity (TP/P)'].apply(lambda x:round(x,3))
    df['specificity (TN/N)'] = (df['TN']/(df['TN'] + df['FP'])).apply(lambda x:round(x,3))
    df['precision (TP/TP+FP)'] = (df['TP']/(df['TP'] + df['FP'])).apply(lambda x:round(x,3))
    df['ACC'] = (df['TP'] + df['TN'])/(df['TP'] + df['TN'] + df['FP'] + df['FN'])
    df['ACC'] = df['ACC'].apply(lambda x:round(x,3))

    df.to_csv(args.output, sep='\t')

    
if __name__ == "__main__":
    main()