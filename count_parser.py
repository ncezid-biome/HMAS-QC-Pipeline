#!/usr/bin/env python

import argparse
import pandas as pd



def main():

	# note 
    # 1. the full format count table file is required !
    # 2. assume the sample/primer_group pair (as column name) in count table is in format of sample.primer_group
    #
    # the script can be run as: python3 count_parser.py -c count_table -s sample_file -o output.tsv
    parser = argparse.ArgumentParser(description = 'parse the final count table file in full format')
    parser.add_argument('-c', '--count_file', metavar = '', required = True, help = 'Specify count table')
    parser.add_argument('-s', '--sample_file', metavar = '', required = True, help = 'Specify sample list file')
    parser.add_argument('-o', '--output', metavar = '', required = True, help = 'Specify output file')
    args = parser.parse_args()


    # read count_table
    df = pd.read_csv(args.count_file, sep='\t')
    # print (df.shape)

    # read sample list file
    sample_list = pd.read_csv(args.sample_file, names = ['sample'])['sample'].tolist()
    sample_list.sort(key=str.lower)
    # print (sample_list)

    #1. filter columns to only contain sample names in sample_list
    df.columns = df.columns.str.strip() # remove potential spaces
    df = df[[i for i in df.columns if i.split('.')[0] in sample_list]]
    # print (df.shape)

    #2. make the new column name list 
    n_column = list(set([i.split('.')[1] for i in df.columns]))
    n_column.sort(key=str.lower)
    print(len(n_column)) # 2451, 10 less than original 2461 ? removed by Mothur ?

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
    final_df.to_csv(args.output, sep='\t')



if __name__ == "__main__":
    main()