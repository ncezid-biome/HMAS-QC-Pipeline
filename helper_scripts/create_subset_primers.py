#!/usr/bin/env python

import random
import pandas as pd
import utilities
import argparse
import time

'''
This script creates a subset of the given oligos file (our full oligos file has 2461 primers), by creating an union 
of a given black-listed primers (usually those not performing primers), and a given random number of primers (i.e. 2000)

The outout is 2 oligos file:
current_time_subset_2461.oligos (the union of black-list primers + 2000 random primers)
current_time_remainder_2461.oligos (2461 - current_time_subset_2461.oligos)
'''
def parse_argument():
    parser = argparse.ArgumentParser(prog = 'create_subset_primers.py')
    #blacklist primers are those primers that are predicted to be non-performing, usually that comes from
    #the 'not_match_primers.txt' file from running primersearch
    parser.add_argument('-b', '--blist', metavar = '', required = True, help = 'Specify blacklist primers file')
    parser.add_argument('-o', '--oligos', metavar = '', required = True, help = 'Specify oligos file')
    parser.add_argument('-n', '--number', metavar = '', required = False, help = 'Specify num of primers we randomly pick')
    
    return parser.parse_args()

def create_subset_oligos(black_list_text, oligos_file, primer_num=2000):
    
    # black_list_text = r'/scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/T3Pio_Data/19isolates/primersearch/full_not_match_primers.txt'
    black_list = pd.read_csv(black_list_text,names = ['primer']).dropna()['primer'].tolist()
    print (len(black_list))

    # oligos_file = r'/scicomp/home-pure/qtl7/test/hmas_test/024_demulx_0mis_data/data/M3235_22_024.oligos'
    primers = utilities.Primers(oligos_file)
    full_primer_list = primers.pnames
    print (f"has total: {len(full_primer_list)} primers")

    random_primer_list = random.sample(full_primer_list,primer_num)
    random_plus_black_primer_list = list(set(random_primer_list) | set(black_list))
    print (f"has total {len(random_plus_black_primer_list)} random_plus_black_primers")

    left_over_primer_list = list(set(full_primer_list) - set(random_plus_black_primer_list))
    print (f"still have total: {len(left_over_primer_list)} left")

    # subset_oligos_file = r'/scicomp/home-pure/qtl7/HMAS_QC_Pipeline/helper_scripts/subset_2461.oligos'
    # remainder_oligos_file = r'/scicomp/home-pure/qtl7/HMAS_QC_Pipeline/helper_scripts/remainder_2461.oligos'
    
    current_time = time.strftime("%y%m%d%H%M", time.localtime())
    subset_oligos_file = f'{current_time}_subset_2461.oligos'
    remainder_oligos_file = f'{current_time}_remainder_2461.oligos'

    # convert a list of primer ID into a dictionary 
    # (key: primer ID, value: ['primer', forwar_primer, reverse_primer, primer_ID])
    def make_primer_dict(primer_list):
        df_dict = {}
        for primer in primer_list:
            df_dict[primer] = ['primer', primers.pseqs[primer][0], utilities.revcomp(primers.pseqs[primer][1]), primer]
        
        return df_dict
        
    df = pd.DataFrame.from_dict(make_primer_dict(random_plus_black_primer_list), orient='index')
    df.to_csv(subset_oligos_file, sep='\t', index=False, header=False)

    df = pd.DataFrame.from_dict(make_primer_dict(left_over_primer_list), orient='index')
    df.to_csv(remainder_oligos_file, sep='\t', index=False, header=False)


if __name__ == "__main__":
    args = parse_argument()
    if args.number:
        create_subset_oligos(args.blist, args.oligos, args.number)
    else:
        create_subset_oligos(args.blist, args.oligos)