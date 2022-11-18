
import sys
sys.path.insert(0,r'..') # to allow import packages in parent folder
import argparse
import pandas as pd
import run_grep
import confusion_matrix as cm
import count_parser as cp
import utilities as u
import random

'''
This script is meant to deuniq the final.fasta file that comes out of Mothur QC pipeline.
The final.fasta file contains the 'high-quality' representative unique sequences. Ex. one such
sequence might have abundance of 20, and this script will help to 'create' all these 20 sequences
, with the seq-ID being the original sequence-id, concatenated with a numeric number (from 1 to 20, in this case),
and the sample-primer group that's associated with this sequence

>M03235_53_000000000-K394R_1_1101_11980_21022_10_08_0810.OG0000294primerGroup8
>M03235_53_000000000-K394R_1_1101_11980_21022 was the original seq-id
10 is the numeric number, up to the abundance value
08_0810.OG0000294primerGroup8 is the sample-primer group
'''

# convert one single list to a csv file without header and index
def write_singlelist_to_csv(list_file, new_csv_file):
    df = pd.DataFrame({'header':list_file})
    df.to_csv(new_csv_file, header=False, index=False)
    
def parse_argument():
    # note
    # 1. the full format count table file is required !
    # 2. sample_file is an one-column csv file containing the names of the sample we're analyzing
    #
    # the script can be run as: python3 deuniq_final_fasta.py -c juno.final.full.count_table -s M3235-22-024-sample.csv
    #                                                         -f juno.final.fasta -d output_dir
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--count_file', metavar = '', required = True, help = 'Specify count table')
    parser.add_argument('-s', '--sample_file', metavar = '', required = True, help = 'Specify sample list file')
    parser.add_argument('-f', '--fasta', metavar = '', required = True, help = 'Specify fasta file output from HMAS QC pipeline (should have "final" in the filename).')
    parser.add_argument('-d', '--file_dir', metavar = '', required = True, help = 'Specify output directory')
    
    
    return parser.parse_args()
    
# this updated method will record the abundance for each selected sequence
# and it will also extract/create fasta files for each sample
def write_fasta_by_seqid_count_tuple(ori_fasta, seqid_count_tuples, new_fasta, sample_list, file_dir):
    fasta_dict = u.create_fasta_dict(ori_fasta)
    fasta_list = []
    for seqid,sample_primer,count in seqid_count_tuples:
        for idx in range(1, count+1):
            fasta_list.append(f">{seqid}_{idx}_{sample_primer}")
            fasta_list.append(fasta_dict[seqid])

    fasta_list.append('')
    with open(f"{file_dir}/{new_fasta}", 'w') as f:
        f.write("\n".join(fasta_list))

    for sample in sample_list:
        grep_command = f"grep -A1 {sample} {file_dir}/{new_fasta} > {file_dir}/{sample}.fasta"
        run_grep.exec_grep_shell(grep_command)
        

if __name__ == '__main__':
    
    args = parse_argument()
    sample_list = pd.read_csv(args.sample_file, names = ['sample'])['sample'].tolist()
    sample_list.sort(key=str.lower)
    print (sample_list)
    
    ## read full count table and convert it to long format (from wide format)
    #  with 3 columns 'seq', 'sample_primer' and 'count'
    ori_count_df = cp.read_count_table(args.count_file)
    melt_count_df = cm.make_melt_count_df(ori_count_df, sample_list)

    ## this is to create a list of random 100 primer pairs
    # oligo_primers = u.Primers(settings.OLIGO_FILE)
    # random.seed(173) # give it a seed
    # top_100_list = random.sample(list(oligo_primers.pseqs),100) 
    # write_singlelist_to_csv(top_100_list,'random_100_primers.csv')
    # new_count_df = melt_count_df.loc[melt_count_df['sample_primer'].str.strip().str.split('.').str[1].isin(top_100_list)]
    
    # this instead use all available primer pairs in the count table
    new_count_df = melt_count_df.copy()

    # file_dir = r'/scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/CIMS/Salmonella/cgMLST_pilot/M347-21-026_fastq/qcrun_result/M3235_22_024_run/update20221118'
    write_fasta_by_seqid_count_tuple(args.fasta, 
                                    list(zip(new_count_df['seq'],new_count_df['sample_primer'], new_count_df['count'])), 
                                    'deuniq_final.fasta',
                                    sample_list,
                                    args.file_dir)