from Bio import SeqIO
import argparse

def read_tab_delimited_file(file_path):
    """
    Read a tab-delimited file and return a dictionary with the first column as keys
    and the remaining columns as values. 
    Here we read in a 3 column file only
    """
    data = {}
    with open(file_path, 'r') as file:
        for line in file:
            columns = line.strip().split()
            if len(columns) == 3:
                key = columns[0]
                data[key] = columns[1:]
    return data

def filter_fasta_file(input_fasta, keys, output_fasta):
    """
    Filter a FASTA file based on a set of keys and write the filtered sequences to a new FASTA file.
    """
    # with open(output_fasta, 'w') as output_file:
    #     for record in SeqIO.parse(input_fasta, 'fasta'):
    #         # Check if any key appears in the record ID
    #         if any(key in record.id for key in keys):
    #             SeqIO.write(record, output_file, 'fasta')
    
    # not using SeqIO.write because it automatically wraps sequence 60bp per line
    filtered_records = (record for record in SeqIO.parse(input_fasta, 'fasta') if any(key in record.id for key in keys))
    with open(output_fasta, 'w') as output_file:
        seqs_list = [f">{rec.id}\n{rec.seq}" for rec in filtered_records]
        output_file.write('\n'.join(seqs_list) + '\n')

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Filter a FASTA file based on keys from a tab-delimited file.')
    parser.add_argument('-i', '--input_fasta', required = True, help='Input FASTA file')
    parser.add_argument('-f', '--tab_delimited_file', required = True, help='Tab-delimited file with keys')
    parser.add_argument('-o', '--output_fasta', required = True, help='Output FASTA file')

    # Parse arguments
    args = parser.parse_args()

    # Read tab-delimited file
    tab_delimited_data = read_tab_delimited_file(args.tab_delimited_file)

    # Extract keys
    keys = set(tab_delimited_data.keys())  # Convert keys to a set for uniqueness
    # Filter FASTA file
    filter_fasta_file(args.input_fasta, keys, args.output_fasta)