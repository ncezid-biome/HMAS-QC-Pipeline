import argparse, os
import pandas as pd
import numpy as np
from pathlib import Path


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


# helper method to select only intended sample (sample-primer) groups from a count_table dataframe 
def filter_count_table(df, sample_list):

	df.columns = df.columns.str.strip() # remove potential spaces
	primer_list = [i for i in df.columns if i.split('.')[0] in sample_list] # filter out unintended samples
	new_df = df[primer_list]
	new_df['seq'] = df.iloc[:,0] # keep the 1st column (sequence)

	return new_df


def main():

    #### the full format count table file is required !
    #
    # the script can be run as:
    # python3 count_helper.py -c juno.final.full.count_table -s M347-21-026-sample.csv
    parser = argparse.ArgumentParser(description = 'parse the final count table file in full format')
    parser.add_argument('-c', '--count_file', metavar = '', required = True, help = 'Specify full count table')
    parser.add_argument('-s', '--sample_file', metavar = '', required = True, help = 'Specify sample list file')
    args = parser.parse_args()

    df = read_count_table(args.count_file)
    print (df.head(n=5))
    print (df.shape)

    # read sample list file
    sample_list = pd.read_csv(args.sample_file, names = ['sample'])['sample'].tolist()
    sample_list.sort(key=str.lower)
    print (sample_list)

    df = filter_count_table(df, sample_list)
    print (df.head(n=5))
    print (df.shape)


if __name__ == "__main__":
    main()