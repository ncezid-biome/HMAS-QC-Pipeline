#!/usr/bin/env python3

import pandas as pd
from pathlib import Path
import glob
import subprocess
import utilities
from Bio import SeqIO
import argparse

import random

'''
This script generates pairwise difference matrix for given amplicon sequnence fasta files. It checks for all the primer
pair sequences/allel sites between any two of those amplicon fasta files. And reports the number of differences if those
two sequences are not exactly the same in either the standard orientation or one of the sequences is in reverse compliment
position
'''

#primer information (can be passed on through command line arguments)
oligos_file = r'/scicomp/home-pure/qtl7/test/hmas_test/024_demulx_0mis_data/data/M3235_22_024.oligos'


#file_dir = r'/scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/validation/mlstComparison/illumina/Salm/validation/1711WAJJP_1/primersearch'
#file_dir = r'/scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/CIMS/Salmonella/degen_primer_test/amplicons/assembly/11amplicons_extract'
#file_dir = r'/scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/T3Pio_Data/19isolates/primersearch/amplicons_extract'
#file_dir = r'/scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/validation/mlstComparison/illumina/Salm/validation/09_1804MLJMP_1/primersearch'


def parse_argument():

    parser = argparse.ArgumentParser(prog = 'create_report.py')
    parser.add_argument('-o', '--output', metavar = '', required = True, help = 'Specify output file')
    parser.add_argument('-d', '--directory', metavar = '', required = True, help = 'Specify fasta file directory')
    parser.add_argument('-p', '--primers', metavar = '', required = False, help = 'Specify oligos/(primer) file')
    return parser.parse_args()

#helper method to check if two sequences are different, returns True if they're different
# it performs the checking in 2 orientations
def check_diff_by_primer(seq1,seq2):
    if seq1 == seq2:
        return False
    elif seq1 == utilities.revcomp(seq2):
        return False
    else:
        return True
    
def pairwise_by_allprimers(file_dir, full_primer_list):
    '''
    this method performs pairwise difference checking on all the isolate fasta files in the given directory, over all 
    primers/allel sites.  Each cell is in the format of: # of difference / # of total common primers/allel sites
    between 2 isolates.

    Parameters
    ----------
    file_dir: the directory which holds the fasta files
    full_primer_list: the list of all primers/allel sites

    Returns
    ----------
    df_list: a list of list (essentially a matrix)
    '''
    # same_seq_list = []
    # diff_seq_list = []
    
    df_list = []
    for f in glob.glob(f'{file_dir}/*.fasta'):
    # for f in glob.glob(f'{file_dir}/2012K-1532_extractedAmplicons.fasta'):

        row_dict = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
        row_list = []
        for k in glob.glob(f'{file_dir}/*.fasta'):
        # for k in glob.glob(f'{file_dir}/Sal_JKX_2015K-0074_extractedAmplicons.fasta'):

            col_dict = SeqIO.to_dict(SeqIO.parse(k, "fasta"))
            diff_count = 0
            total_primer = len(full_primer_list)
            for primer in full_primer_list:
                
                row_seq_list = [val for key,val in row_dict.items() if primer in key]
                col_seq_list = [val for key,val in col_dict.items() if primer in key]   
                
                if len(row_seq_list) < 1 or len(col_seq_list) < 1:
                    total_primer -= 1
                else:
                    if len(row_seq_list) > 1 or len(col_seq_list) > 1:
                        print (f"{primer} has more than 1 match in {f} or {k}")
                    if check_diff_by_primer(row_seq_list[0].seq, col_seq_list[0].seq):
                        diff_count += 1
                        
                    #     diff_seq_list.append((primer,row_seq_list[0].seq, col_seq_list[0].seq))
                    # else:
                    #     same_seq_list.append((primer,row_seq_list[0].seq, col_seq_list[0].seq))

            row_list.append(f" {diff_count}/{total_primer}")
        df_list.append(row_list)
        
    # print (f"same seq count is: {len(same_seq_list)}")
    # print (f"diff seq count is: {len(diff_seq_list)}")
    # print (random.sample(same_seq_list,10))
    # print ('********************************')
    # print ('********************************')
    # print ('********************************')
    # print (random.sample(diff_seq_list,10))
    
    return df_list
            
#offer similar functionality as 'pairwise_by_allprimers', but it calls upon vsearch --search_exact command to 
#check for any two exact same sequences between 2 fasta files. Which is a bit different than 'pairwise_by_allprimers'
#might be deleted in the offical release
def pairwise_vsearch(file_dir):
    df_list = []
    for f in glob.glob(f'{file_dir}/*.fasta'):
        # row_file = Path(f).stem.split('.')[0]
        row_list = []
        for k in glob.glob(f'{file_dir}/*.fasta'):
            # col_file = Path(k).stem.split('.')[0]
            
            command1 = f'vsearch --search_exact {f} -db {k} --userfields target+query --userout  {file_dir}/output.match.final.txt 2>&1 | tee {file_dir}/out.test'
            p = subprocess.run(command1, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, text=True, shell=True)
            
            command2 = f'tail -n1 {file_dir}/out.test | cut -d " " -f4,6'
            p = subprocess.run(command2, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, shell=True)
            print (p.stdout)
            n1,n2 = tuple(p.stdout.split())
            # print (int(n2)-int(n1))
            row_list.append(int(n2)-int(n1))
        df_list.append(row_list)
        
    return df_list
        

if __name__ == "__main__":

    args = parse_argument()
    
    if args.primers:
        primers = utilities.Primers(args.primers)
    else:
        primers = utilities.Primers(oligos_file)
    full_primer_list = primers.pnames
    
    df_list = pairwise_by_allprimers(args.directory, full_primer_list)
    
    file_name_list = [Path(f).stem.split('.')[0] for f in glob.glob(f'{args.directory}/*.fasta') ]
    df = pd.DataFrame(df_list,columns=file_name_list)
    df.index = file_name_list
    df.to_csv(f'{args.output}')
    # df.to_csv(f'{file_dir}/pairwise_diff_revcomp_1711WAJJP_1.csv')
    # df.to_csv(f'{file_dir}/pairwise_diff_11_allprimers_revcomp.csv')