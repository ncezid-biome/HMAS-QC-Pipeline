#!/usr/bin/env python

import os, logging, argparse

# Take user inputs (filename prefix, extension) and make a list of the files in current wd
# Format should be: forward \t reverse \t fwd_index

# Consider modifying for other use cases: group name and R1,R2 (no index files); forward and reverse index files

LOG_FORMAT = "%(levelname)s %(asctime)s - %(message)s"
logging.basicConfig(filename = 'make_file.log', format = LOG_FORMAT, level = logging.DEBUG)
logger = logging.getLogger()

parser = argparse.ArgumentParser(description = 'Make a list of the files to run in batch mode')
parser.add_argument('-p', '--prefix', metavar='', required=True, help = 'A prefix that identifies the files you want to run')
parser.add_argument('-x', '--extension', metavar='', required=True, help = 'An extension that your files have, e.g. .gz')

# Note: potential problem: user puts .fastq.gz instead of .gz, has a mix of files in directory
    # Should be ok if they use prefix
# Note: what if the user's files don't have a prefix? 

args = parser.parse_args()

prefix = args.prefix
extension = args.extension

logger.info("Arguments passed: prefix = {}, extension = {}".format(prefix, extension))

inputDir = os.getcwd()
outFile = prefix + '.paired.files'

logger.info("Input directory: {}".format(inputDir))
logger.info("Output file name: {}".format(outFile))

f_read_id = 'R1'
r_read_id = 'R2'
f_indx_id = 'I1' 

fwd_reads = []
rev_reads = []
fwd_inds = []
file_list = [fwd_reads, rev_reads, fwd_inds]

#NOTE!!!!!  This should be modified so that both '.gz' and 'gz' will work as extension
for f in os.listdir():
    fname, fext = os.path.splitext(f)
    if fext == extension:
        logger.info("File added: {}".format(fname))
        if f_read_id in fname:
            fwd_reads.append(f)
        if r_read_id in fname:
            rev_reads.append(f)
        if f_indx_id in fname:
            fwd_inds.append(f)

# Note: No error handling!
# Note: Easy to read, but not very robust
# Note: If R1, R2, I1, I2 exist anywhere else in the filename, this will mess up

# If the lists are not all the same length, raise an exception
iter_list = iter(file_list)
list_len = len(next(iter_list))
if not all(len(l) == list_len for l in iter_list):
    raise ValueError('Every sample must have exactly one forward and one reverse read, and one index')
    logger.error("Each sample must have exactly one forward read, one reverse read, and one index file.")
    # Maybe print the lists here so user can see what happened


# Sort the lists and write to formatted file
# Note: I'm not entirely comfortable with the robustness of sorting... 
with open(outFile, 'w') as f:
    for a,b,c in zip(fwd_reads, rev_reads, fwd_inds):
        logger.info(print('{}\t{}\t{}\n'.format(a,b,c)))
        f.write('{}\t{}\t{}\n'.format(a,b,c))

# Note: this outputs characters to STDOUT; need to somehow make them go away.
