#!/usr/bin/env python

import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pathlib import Path
import argparse
import logging
from logging.handlers import RotatingFileHandler

amplicon_file_extension = '_extractedAmplicons.fasta'
non_match_primer_file_extension = '_not_match_primers.txt'
max_seqid_len = 80
LOG_FILE = 'parse_primersearch.log'

log_formatter = logging.Formatter('%(asctime)s %(levelname)s %(filename)s(%(lineno)d) - %(message)s')
log_handler = RotatingFileHandler(LOG_FILE, mode='a', maxBytes=5*1024*1024, backupCount=5)
log_handler.setFormatter(log_formatter)
log_handler.setLevel(logging.INFO)

logger = logging.getLogger('parse_primersearch')
logger.setLevel(logging.INFO)
logger.addHandler(log_handler)

def extractIndex(hit_line_string):
    regex = re.compile(r"strand\s+at\s+([0-9]{1,})\s+with")
    clean_str = hit_line_string.replace('[', '').replace(']', '')
    return int(regex.findall(clean_str)[0])

def extractPrimerLength(hit_line_string):
    regex = re.compile(r"^(\w+)\s+hits")
    clean_str = re.sub(r'\[\w+\]', 'N', hit_line_string.strip())
    return len(regex.findall(clean_str)[0])

def extractAmpliconLength(amplimer_length_line):
    regex = re.compile(r"[0-9]+")
    matches = regex.findall(amplimer_length_line)
    return int(matches[0]) if matches else 0

def parsePrimerSearch(primersearch_results, full_length_dict, file_base, max_amplicon_len):
    extracted_amplicon_list = []
    not_match_primer_list = []
    seqid_list = []
    seqid_primer_list = []
    seqid_isolate_list = []

    with open(primersearch_results, 'r') as inputFile:
        for line in inputFile:
            if "Amplimer " in line:
                ampl_count += 1
                seq_id = inputFile.readline().strip().replace("Sequence: ", "")
                description = inputFile.readline().strip()
                forward_hit_line = inputFile.readline().strip()
                forward_primer_length = extractPrimerLength(forward_hit_line)
                reverse_hit_line = inputFile.readline().strip()
                reverse_primer_length = extractPrimerLength(reverse_hit_line)

                if "forward" in forward_hit_line:
                    forwardHitPosition = extractIndex(forward_hit_line)
                else:
                    forwardHitPosition = extractIndex(reverse_hit_line)
                    logger.info(f'Amplimer {ampl_count}, {primer_name} hit reverse strand in {seq_id}')

                amplimer_length_line = inputFile.readline().strip()
                amplicon_length = extractAmpliconLength(amplimer_length_line)

                startIndex = forwardHitPosition - 1 + forward_primer_length
                endIndex = startIndex + amplicon_length - forward_primer_length - reverse_primer_length

                if endIndex - startIndex <= max_amplicon_len:
                    rec = full_length_dict[seq_id]
                    extracted_seq = rec.seq[startIndex:endIndex]
                    rec_id = f'{primer_name}-{file_base}-ampl{ampl_count}'[:max_seqid_len]
                    ampliconRec = SeqRecord(extracted_seq, id=rec_id, description="")
                    extracted_amplicon_list.append(ampliconRec)
                    seqid_list.append(rec_id)
                    seqid_primer_list.append(primer_name)
                    seqid_isolate_list.append(file_base)

                    if not_match_primer_list and ampl_count == 1:
                        not_match_primer_list.pop()
                else:
                    logger.info(f'{seq_id} Amplimer {ampl_count} too long: {endIndex - startIndex}')

            elif "Primer name " in line:
                primer_name = re.search(r'Primer name\s+([\w-]+)', line).group(1)
                not_match_primer_list.append(primer_name)
                ampl_count = 0

    return extracted_amplicon_list, not_match_primer_list

def parse_args():
    parser = argparse.ArgumentParser(description="Parse primersearch results")
    parser.add_argument('-s', '--sequence', required=True, help='Original FASTA sequence file')
    parser.add_argument('-r', '--results', required=True, help='Primersearch output .ps file')
    parser.add_argument('-l', '--amp_len', type=int, required=True, help='allowed max amplicon length')
    
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    file_base = Path(args.sequence).stem
    full_length_dict = SeqIO.to_dict(SeqIO.parse(args.sequence, "fasta"))
    amplicons, not_match_primer_list = parsePrimerSearch(args.results, full_length_dict, file_base, args.amp_len)

    with open(f'{file_base}{amplicon_file_extension}', 'w') as f:
        for rec in amplicons:
            f.write(f">{rec.id}\n{rec.seq}\n")

    with open(f'{file_base}{non_match_primer_file_extension}', 'w') as f:
            f.write('\n'.join(not_match_primer_list) + '\n')
