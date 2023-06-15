import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

'''
this script is ued to generate histograms of primer counts for Jo
https://cdc.sharepoint.com/:p:/r/teams/NCEZID-BIOME/Shared%20Documents/HMAS/run_report/M3235_22_024_primer_hist.pptx?d=waa79172b1337477a889a0a20d512f022&csf=1&web=1&e=qLOPYB

2 histograms can be made: 
after primer trimming:  (need to read os.path.join(folder,'temp', f"{folder_name}.unique.fasta"))
note: technically this is before denoising, but it's fairly close to after primer trimming because
the steps in between them remove minimum sequences.
or
after QC pipeline finishes (read os.path.join(folder,f"{folder_name}.final.unique.fasta"))

'''
## we were counting the lines of the gz file, thus undercounting the actual lines !!
# raw_seq_count = {
#     '2013K_1246': 69.8,
#     '2015K_0074': 20.1,
#     '2014K_0527': 59.1,
#     '2013K_1828': 90.3,
#     '2013K_1153': 46.9,
#     '2013K_0463': 35.3,
#     '08_0810': 60.1,
#     '2014K_0421': 51.7,
#     '2013K_0676': 21.1,
#     '2016K_0167': 60.7,
#     'AR_0409': 115.1,
#     '2013K_1067': 18.8,
#     '2016K_0878': 68.8,
#     '2014K_0979': 0.2,
#     '2010K_1649': 25.4,
#     '2010K_0968': 12.9,
#     '2011K_0222': 49.5,
#     '2011K_0052': 44.3,
#     '2010K_2370': 44.3
# }

#the values come from running count_rawseqs.py
raw_seq_count = {
    '2013K_1246': 427.0,
    '2015K_0074': 120.0,
    '2014K_0527': 356.4,
    '2013K_1828': 548.1,
    '2013K_1153': 284.7,
    '2013K_0463': 210.6,
    '08_0810': 357.4,
    '2014K_0421': 316.2,
    '2013K_0676': 124.7,
    '2016K_0167': 370.4,
    'AR_0409': 689.0,
    '2013K_1067': 113.5,
    '2016K_0878': 416.0,
    '2014K_0979': 1.3,
    '2010K_1649': 152.5,
    '2010K_0968': 77.3,
    '2011K_0222': 281.8,
    '2011K_0052': 265.3,
    '2010K_2370': 266.9
}


def parse_argument():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', metavar = '', required = True, help = 'Specify input file folder')
    parser.add_argument('-o', '--output', metavar = '', required = True, help = 'Specify output file')
    parser.add_argument("--exclude", nargs="+", help="Folders to exclude")
    parser.add_argument('-p', '--primerlist', metavar = '', required = True, help = 'Specify primer list file')
    return parser.parse_args()


def get_folders(args):
    
    # Get all folder names in the input folder
    # folder_names = [os.path.join(args.input, f,'temp','cutadapt') for f in os.listdir(args.input) if os.path.isdir(os.path.join(args.input, f))]
    folder_names = []
    for f in os.listdir(args.input):
        # Check if the folder should be excluded
        if args.exclude and f in args.exclude:
            continue
        if os.path.isdir(os.path.join(args.input, f)):
            folder_names.append(os.path.join(args.input, f))
        
    return folder_names


def create_amplicon_count_dict(folder, primer_list):
    
    # Initialize an empty dictionary to store the counts
    counts = {}
    folder_name = os.path.basename(folder)
    # fasta_file = os.path.join(folder,f"{folder_name}.final.unique.fasta") # this is after QC
    fasta_file = os.path.join(folder,'temp', f"{folder_name}.unique.fasta") # this is after primer trimming

    # Open the fasta file
    with open(fasta_file) as file:
        # Read the file line by line
        for line in file:
            # Check if the line starts with '>'
            if line.startswith('>'):
                # Extract the primer name from the line
                primer_name = line.split('=')[1]
                # Extract the size value from the line
                size = float(line.split('=')[-1].split(';')[0].replace('size=', '')) / raw_seq_count[folder_name]
                
                # Check if the primer is in the primer list
                if primer_name in primer_list['primer'].values:
                    # Add the size to the counts dictionary
                    counts[primer_name] = counts.get(primer_name, 0) + size
                    
    return counts


def calculate_average_counts(*dictionaries):
    cumulative_counts = {}
    total_dictionaries = len(dictionaries)

    # Iterate through each dictionary
    for dictionary in dictionaries:
        # Iterate through each key-value pair
        for primer, count in dictionary.items():
            cumulative_counts.setdefault(primer, 0)
            cumulative_counts[primer] += count

    # Calculate the average count for each primer
    average_counts = {primer: count / total_dictionaries for primer, count in cumulative_counts.items()}

    return average_counts



if __name__ == "__main__":
    
    args = parse_argument()
    folder_names = get_folders(args)
    
    # Read the primer list file into a pandas DataFrame
    primers = pd.read_csv(args.primerlist, header=None, names=['primer'])
    
    dict_list = []
    for folder in folder_names:
        count_dict = create_amplicon_count_dict(folder, primers)          
        dict_list.append(count_dict)
        
    average_count_dict = calculate_average_counts(*dict_list)
    
    # Create a DataFrame from the counts dictionary
    df = pd.DataFrame.from_dict(average_count_dict, orient='index', columns=['count'])


    
    bin_edges = np.arange(min(df['count']), max(df['count']) + 0.1, 0.1)
    arr = plt.hist(df['count'], histtype='stepfilled', bins=bin_edges)
    
    plt.xlabel(f'observed number of reads / expected number of reads (per primer pair)')
    plt.ylabel('number of primer pairs')
    plt.title('Histogram of Primer Counts (after primer trimming)')
    # plt.title('Histogram of Primer Counts (after QC)')
    
    # Set the x-tick marks
    xticks = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5] #might be updated to the appropriate value
    plt.xticks(xticks)

    # for the text in the histogram
    for i in range(len(bin_edges)-1):
        if int(arr[0][i]) > 0:
            plt.text(arr[1][i]+0.01,arr[0][i]+3,str(int(arr[0][i])), fontsize=8, fontfamily='monospace')
    
    # Add lines between bins
    for bin_edge in arr[1]:
        # print (bin_edge)
        plt.axvline(bin_edge, color='white', linestyle='--', linewidth=0.5)
    #add a red line for the 1.0 bin (the 1.0003 might be updated to the appropriate value)
    plt.axvline(1.0003, color='red', linestyle='--', linewidth=0.7)
    
    # plt.savefig(args.output, dpi=1600)
    plt.show()

