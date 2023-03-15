import run_blast as blast
import pandas as pd
import argparse
import glob
from Bio import SeqIO, SeqRecord
# from Bio.SeqRecord import SeqRecord
from pathlib import Path
import os

'''
This script uses blast to find the gene sequence for the given amplicon sequence(which is part of the gene)
It takes 2 parameters of folders for the amplicon fasta files and fna files(fna file comes from parsing the gbk file of the assembled genomes)
Right now, the fna file is generate from running modified Sean's moduleFile.py.

The output is a fasta file for each given fna file, with extension _extracteGenes.fasta, and it's in the same folder
of the fna file

## assume file name format for fna and ampilcon files are:
## SRR1686418.fna  SRR1686418_assembled_extractedAmplicons.fasta

'''

def parse_argument():
    # note
    # the script can be run as:
    # 
    parser = argparse.ArgumentParser(prog = 'get_genes_from_amplicons.py')
    parser.add_argument('-a', '--amplicon', metavar = '', required = True, help = 'Specify amplicons file folder')
    parser.add_argument('-f', '--fna', metavar = '', required = True, help = 'Specify fna file folder')
    
    return parser.parse_args()


def create_blast_df(query_fasta, subject_fasta, max_hits):

    bcolnames = ["query", "subject", "query_len", "subj_len", "len_aln", "eval", "cov", "pident", "mismatch"]

    blast_file = blast.blast(query_fasta, subject_fasta, os.path.basename(query_fasta), max_hits)
    blast_df = pd.read_csv(blast_file, sep='\t', names=bcolnames)
    # requires 100% match for now
    blast_df = blast_df[ (blast_df["pident"] >= 100) & (blast_df["cov"] >= 100) ]
    blast_df = blast_df[["query", "subject"]] # we only need these 2 fields
    
    return blast_df
    
    
if __name__ == "__main__":
    
    args = parse_argument()
    fna_list = glob.glob(f"{args.fna}/*.fna")
    amplicon_list = glob.glob(f"{args.amplicon}/*.fasta")
    
    for fna in fna_list:
        file_base = Path(fna).stem   #get SRR1686418 from  directory/SRR1686418.fna
        #find matching amplicon file, should have exactly one match !
        amplicon_file = [f for f in amplicon_list if file_base in Path(f).stem]
        if len(amplicon_file) != 1:
            print (f"warning ! {file_base} has {len(amplicon_file)} matched amplicon files")
        fna_seq_dict = SeqIO.to_dict(SeqIO.parse(fna, "fasta"))
        blast_df = create_blast_df(amplicon_file[0], fna, 10)
        #mapping of seq_ID between amplicon and gene
        gene_amplicon_dict = dict(zip(blast_df["subject"], blast_df["query"]))
        
        gene_record_list = []
        for gene in gene_amplicon_dict:
            gene_record = fna_seq_dict[gene]
            gene_record.id = f"{gene_record.id}_{gene_amplicon_dict[gene_record.id]}"
            gene_record_list.append(gene_record)
            
        SeqIO.write(gene_record_list, f"{args.fna}/{file_base}_extracteGenes.fasta", "fasta")