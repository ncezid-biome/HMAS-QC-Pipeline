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
# pd.set_option('display.max_rows', None)
import count_parser as cp
import sys
import utilities
import settings
import logging

sys.path.insert(0,r'..') # to allow import packages in parent folder
import group
import run_grep

LOG_FORMAT = "%(levelname)s %(asctime)s - %(message)s"
logging.basicConfig(level = logging.DEBUG,format = LOG_FORMAT,filename = f"confusion_matrix.log")
logger = logging.getLogger(__name__)


def sample_primer_df_to_dict(df):
    '''
    this method converts a 2 column('sample', 'primer') dataframe into a dictionary, where
    sample is the key and primer is the value for that key

    sample    primer
    08_0810   OG0000294primerGroup8
    08_0810   OG0000296primerGroup7
    AR_0409   OG0000297primerGroup3
    AR_0409   OG0000298primerGroup4

    Parameters
    ----------
    df: the 2 columns dataframe

    Returns a dictionary
    {'08_0810':'OG0000294primerGroup8,OG0000296primerGroup7','AR_0409':'OG0000297primerGroup3,OG0000298primerGroup4'}

    '''
    df['primer'] = df.groupby(['sample'])['primer'].transform(lambda x: ','.join(x))
    return dict(df.values)


def get_oligo_primer():
    # return the set of primer names in the oligo file
    return set(get_oligo_primer_dict())

def get_oligo_primer_dict():
    '''
    this method generates an instance of Primers class and a dictionary with primer name being the key
    and primer-pair sequences (tuple) being the value

    Returns a dictionary

    '''
    oligo_primers = group.Primers(settings.OLIGO_FILE)

    return oligo_primers.pseqs


def get_abundance(sample_list, df):
    '''
    this method reads the original count table from mothur, reformats it and outputs a dataframe of abundance values
    with sample name as index and primer pairs as columns

                OG0000294primerGroup8   OG0000296primerGroup7   OG0000297primerGroup3
    08_0810     10                      20                      30
    AR_0409     30                      10                      30
    .......

    Parameters
    ----------
    count_file: the original count table file from Mothur
    sample_list: a csv file of the list of samples used
    ori_count_df: dataframe for the count_file

    Returns the abundance table as a dataframe

    '''
    df.columns = df.columns.str.strip() # remove potential spaces
    primer_list = [i for i in df.columns if i.split('.')[0] in sample_list]
    df = df[primer_list]

    sample_pp_df = df.sum() #total abundanceall high quality seqs
    primer_list = list(set([i.split('.')[1] for i in sample_pp_df.index]))
    primer_list.sort(key=str.lower)
    sample_pp_df = cp.split_samle_primer(sample_pp_df, primer_list, sample_list)

    logger.debug(sample_pp_df.mean(axis=1))
    logger.debug(sample_pp_df.shape)

    return sample_pp_df



def make_blast_dict(df):
    '''
    this method process the blast result dataframe, filter by pre-set threshold and generate a sample/primer dictionary
    with sample name as the key and all primers in that sample as the value

    Parameters
    ----------
    df: dataframe for the blast result

    Returns a dictionary

    '''
    df = df[(df['pident'] >= settings.PIDENT) & (df['len_aln']/df['subj_len'] >= settings.P_ALIGN) & (df['cov'] >= settings.PCOV)]
    df = df[['sample','primer']]
    df.drop_duplicates(inplace=True)

    blast_dict = sample_primer_df_to_dict(df)

    # filter out primers that's not in the original 2461 primer_list, off blast result
    for key in blast_dict:
        primers = set(blast_dict[key].split(','))
        diff = primers - get_oligo_primer()
        primers = primers - diff
        blast_dict[key] = ','.join(primers)
    
    return blast_dict

# helper method to get seq from blast dataframe by matched sample_primer
def get_blast_seq(df, sample_primer):

    FILTER = df['sample_primer'] == sample_primer
    return df[FILTER]['seq']

# helper method to get blast result by matched sample_primer and seq
def get_blast_by_seq_sp(df, seq, sample_primer):

    FILTER = (df['sample_primer'] == sample_primer) & (df['seq'] == seq)
    return df[FILTER]

# helper method to get seq from count dataframe by matching sample_primer with value > 0
def get_count_seq(df, sample_primer):

    FILTER = df[sample_primer] > 0
    return df[FILTER]['Representative_Sequence']


def make_count_blast_dict(count_df, df):
    '''
    this method processes the blast result dataframe, filters by pre-set threshold and uses the remaining
    sequences to filter the count dataframe. Lastly it uses the sample/primer columns in count dataframe to 
    generate a sample/primer dictionary with sample name as the key and all primers in that sample as the value

    Parameters
    ----------
    df: dataframe for the blast result
    count_df: dataframe for the original count table

    Returns a dictionary

    '''
    df = df[(df['pident'] >= settings.PIDENT) & (df['len_aln']/df['subj_len'] >= settings.P_ALIGN) & (df['cov'] >= settings.PCOV)]

    #merge 2 dfs (count table and blast result)
    # df = pd.merge(df,count_df,on=['seq'])  
    ## filter count_df by checking if seqs in count_df exists in df
    # this accomplishes the same thing we want and it's far more efficient than merge
    new_count_df = count_df.loc[count_df['seq'].isin(set(df['seq']))]
    new_count_df.drop_duplicates(inplace=True)
    new_count_df['sample'] = new_count_df['sample_primer'].str.strip().str.split('.').str[0]
    new_count_df['primer'] = new_count_df['sample_primer'].str.strip().str.split('.').str[1]
    new_count_df = new_count_df[['sample','primer']]

    count_blast_dict = sample_primer_df_to_dict(new_count_df)

    # filter out primers that's not in the original 2461 primer_list, off blast result
    for key in count_blast_dict:
        primers = set(count_blast_dict[key].split(','))
        primers = primers & get_oligo_primer()
        count_blast_dict[key] = ','.join(primers)
    
    return count_blast_dict


def make_blast_merge_dict(count_df, df):
    '''
    this method merge the count dataframe with the blast dataframe on the same 'seq' and 'sample_primer' and filter
    the result by pre-set threshold value. Then it generates a sample/primer dictionary with sample name as the key 
    and all primers in that sample as the value. These are considered 'true positive' primers for the sample.

    Parameters
    ----------
    df: dataframe for the blast result
    count_df: dataframe for the original count table

    Returns a dictionary

    '''
    df = df[(df['pident'] >= settings.PIDENT) & (df['len_aln']/df['subj_len'] >= settings.P_ALIGN) & (df['cov'] >= settings.PCOV)]

    #merge 2 dfs (count table and blast result)
    df = pd.merge(df,count_df,on=['seq','sample_primer'])
    df = df[['sample','primer']]
    df.drop_duplicates(inplace=True)

    blast_dict = sample_primer_df_to_dict(df)

    # filter out primers that's not in the original 2461 primer_list, off blast result
    for key in blast_dict:
        primers = set(blast_dict[key].split(','))
        primers = primers & get_oligo_primer()
        blast_dict[key] = ','.join(primers)
    
    return blast_dict


def make_melt_count_df(df, sample_list):
    '''
    this method melts the original count table dataframe
    with sample name as the key and all primers in that sample as the value

    Parameters
    ----------
    df: dataframe for the original count table
    sample_list: a csv file of the list of samples used

    Returns the melted dataframe

    '''
    #filter columns to contain only sample names in sample_list
    df.columns = df.columns.str.strip() # remove potential spaces
    seq_df = df.iloc[:,0] # keep the 1st column (sequence)
    primer_list = [i for i in df.columns if i.split('.')[0] in sample_list]
    df = df[primer_list]
    df['seq'] = seq_df

    df = df.melt(id_vars=['seq'], var_name='sample_primer', value_name='count')
    df = df[df['count'] > 0]

    return df
  

def get_set_intersection(set_list):
    '''
    this method retrieves the intersection of all sets passed in a list

    Parameters
    ----------
    set_list: a list of sets

    '''
    _intersect = set.intersection(*set_list)


def get_setitem_occurrence(set_list):
    '''
    this method flattens a list of sets into one single list and generate the occurrence count for each item
    as a dictionary

    Parameters
    ----------
    set_list: a list of sets

    Returns a dictionary (key:item, value:occurrence of that item)
    '''
    flat_list = []
    for sublist in set_list:
        flat_list.extend(sublist) 
    flat_dict = {i:flat_list.count(i) for i in flat_list}
    return flat_dict


def parse_argument():
    # note
    # 1. the full format count table file is required !
    # 2. assume the sample/primer_group pair (as column name) in count table is in format of sample.primer_group
    # 3. sample_file is an one-column csv file containing the names of the sample we're analyzing
    #
    # the script can be run as: python3 confusion_matrix.py -c count_table -s sample_file -b blast_file -o output.txt
    #                               -f juno.final.fasta -r juno_design_primers.fasta
    # - b blast_file is optional
    parser = argparse.ArgumentParser(prog = 'confusion_matrix.py')
    parser.add_argument('-c', '--count_file', metavar = '', required = True, help = 'Specify count table')
    parser.add_argument('-s', '--sample_file', metavar = '', required = True, help = 'Specify sample list file')
    parser.add_argument('-o', '--output', metavar = '', required = True, help = 'Specify output file')
    parser.add_argument('-b', '--blast', metavar = '', required = False, help = 'Specify blast result file')
    parser.add_argument('-f', '--fasta', metavar = '', required = False, help = 'Specify fasta file output from HMAS QC pipeline (should have "final" in the filename).')
    parser.add_argument('-r', '--reference', metavar = '', required = False, help = 'Specify fasta file containing the positive control targets.')
    
    return parser.parse_args()


def build_confusion_matrix(count_dict, blast_nomerge_dict, tp_dict):
    '''
    this method flattens a list of sets into one single list and generate the occurrence count for each item
    as a dictionary

    Parameters
    ----------
    count_dict: a dictionary of observed primers (key: sample_name, value: a comma delimited string of primers)
    blast_nomerge_dict: a dictionary of our predicted primers 
    tp_dict: a dictionary of true positive primers

    Returns
    ----------
    df_dict: a dictionary of TP/FP/TN/FN metrics (key: sample_name, value: a list of number of TP/FP/TN/FN counts )
    FP_list: a list of all FPs
    FN_list: a list of all FNs
    TN_list: a list of all TNs

    '''
    # a dictionary of fasta sequence, with seq_ID being the key and actual sequence being the value
    # fasta_dict = utilities.create_fasta_dict(args.fasta)

    df_dict = {} #to hold the number of tp/tn/fp/fn for each sample
    FP_list, FN_list, TN_list = ([] for i in range(3))

    for key in dict(sorted(blast_nomerge_dict.items())):
        if key not in settings.CONTROL_LIST:
            logger.info(f"sample is: {key}")
            blast_nomerge_primers = set(blast_nomerge_dict[key].split(','))     
            count_primers = set(count_dict[key].split(','))
            tp_primers = set(tp_dict[key].split(','))
            total_primers = get_oligo_primer() # use all 2461 primers
 
            # TP: Predicted amplicons that are observed in the high-quality data
            TP = tp_primers
            logger.info (f" checking TP vs blast_nomerge_primers: {len(TP)} -- {len(TP & blast_nomerge_primers)}")
            logger.info (f" checking TP vs count_primers: {len(TP)} -- {len(TP & count_primers)}")
            # TN: Primer pairs that do not produce data and are not predicted to
            TN = (total_primers - blast_nomerge_primers) & (total_primers - count_primers)
            # FP: High-quality amplicons that are observed but not predicted for the isolate
            FP = (total_primers - TP) & count_primers 
            logger.info (f"checking FP size: {len(FP)} ------- {len(count_primers - TP)}")
            # FN: Predicted amplicons that are not observed in the high-quality data
            FN = (total_primers - count_primers) & blast_nomerge_primers

            logger.info (f"sanity check: {len(total_primers)}, {len(total_primers) == (len(TP)+len(FP)+len(FN)+len(TN))}")
            logger.info (f"more sanity check: TP & FP {len(TP & FP)}")
            logger.info (f"more sanity check: TP & FN {len(TP & FN)}")
            logger.info (f"more sanity check: TP & TN {len(TP & TN)}")
            logger.info (f"more sanity check: FN & FP {len(FN & FP)}")
            logger.info (f"more sanity check: TN & FP {len(TN & FP)}")
            logger.info (f"more sanity check: FN & TN {len(TN & FN)}")

            df_dict[key] = [len(TP),len(FP),len(FN),len(TN)]

            FP_list.append(FP)
            FN_list.append(FN)
            TN_list.append(TN)

    return df_dict, FP_list, FN_list, TN_list


def main():

    args = parse_argument()
    
    # read sample list file
    sample_list = pd.read_csv(args.sample_file, names = ['sample'])['sample'].tolist()
    sample_list.sort(key=str.lower)
    
    ori_count_df = cp.read_count_table(args.count_file)
    melt_count_df = make_melt_count_df(ori_count_df, sample_list)
    blast_df = cp.create_blast_df(args.blast, args.fasta, args.reference, 100)
    abundance_df = get_abundance(sample_list, ori_count_df)
 
    # this is our prediction
    blast_nomerge_dict = make_blast_dict(blast_df)
    # this is our observation
    count_dict = make_count_blast_dict(melt_count_df, blast_df)
    # this is basically the true positives
    tp_dict = make_blast_merge_dict(melt_count_df, blast_df)
    #predicted primers for each sample
    # prediction_dict = run_grep.get_primer_prediction_dict_concurrent()


    df_dict, FP_list, FN_list, TN_list = build_confusion_matrix(count_dict, blast_nomerge_dict, tp_dict)

    
    # get_set_intersection(FP_list)
    # get_set_intersection(FN_list)
    # get_set_intersection(TN_list)

    #list the primers in the order of most frequent appearance
    logger.info (f"\n\nduplicate count in FP: "
           f"{dict(sorted(get_setitem_occurrence(FP_list).items(), key=lambda x:x[1], reverse=True))}")
    logger.info (f"\n\nduplicate count in FN: "
           f"{dict(sorted(get_setitem_occurrence(FN_list).items(), key=lambda x:x[1], reverse=True))}")


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

