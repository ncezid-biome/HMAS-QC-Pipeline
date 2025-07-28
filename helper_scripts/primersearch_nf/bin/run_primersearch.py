#!/usr/bin/env python

import argparse
import shutil
import subprocess
import sys
from pathlib import Path
import logging
from logging.handlers import RotatingFileHandler

LOG_FILE = 'primersearch.log'

# Logging setup
log_formatter = logging.Formatter('%(asctime)s %(levelname)s %(filename)s(%(lineno)d) - %(message)s')
log_handler = RotatingFileHandler(LOG_FILE, mode='a', maxBytes=5*1024*1024, backupCount=5)
log_handler.setFormatter(log_formatter)
log_handler.setLevel(logging.INFO)

logger = logging.getLogger('run_primersearch')
logger.setLevel(logging.INFO)
logger.addHandler(log_handler)

def runPrimerSearch(seq_file, primer_file, output_file, mismatch):
    if shutil.which('primersearch') is None:
        logger.error("primersearch is not found on PATH")
        sys.exit(1)

    # Replace with your own primersearch path if needed
    # primersearch_path = 'primersearch'
    primersearch_path = '/scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/T3Pio_Data/library_versions_t3pio/EMBOSS-6.4.0/emboss/primersearch'
    
    command = [
        primersearch_path,
        '-seqall', seq_file,
        '-infile', primer_file,
        '-mismatchpercent', str(mismatch),
        '-outfile', output_file
    ]

    process = subprocess.run(command, capture_output=True, text=True)
    if process.returncode == 0:
        logger.info(f"Primersearch completed successfully for {seq_file}")
        return output_file
    else:
        logger.error(f"Primersearch failed for {seq_file}")
        sys.exit(1)

def parse_args():
    parser = argparse.ArgumentParser(description="Run EMBOSS primersearch")
    parser.add_argument('-s', '--sequence', required=True, help='FASTA sequence file')
    parser.add_argument('-p', '--primers', required=True, help='Primer list file')
    parser.add_argument('-o', '--output', required=True, help='Output .ps file name')
    parser.add_argument('-m', '--mismatch', type=int, required=True, help='allowed mismatch percentage')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    runPrimerSearch(args.sequence, args.primers, args.output, args.mismatch)
