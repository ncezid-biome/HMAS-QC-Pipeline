#!/usr/bin/env python

import os, logging, argparse, re
import pandas as pd

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

class GetFiles:
    def __init__(self, prefix, suffix):
        self.prefix = prefix
        self.suffix = suffix

    def getFiles(self):
        in_dir = os.getcwd()
        logger.info("Current input directory: {}".format(in_dir))
        fileTable = pd.DataFrame()
        rule = re.compile('[RI][1-2]')
        for f in os.listdir():
            fname, fext = os.path.splitext(f)
            if fext == self.suffix:
                logger.info("File added: {}".format(fname))
                key = rule.search(f)
                short = f.replace(key.group())
                row = {'filename' : f,
                           'type' : key.group(),
                           'shortname' : short}
                fileTable = fileTable.append(row, ignore_Index = True)

    
# To finish
#outFile = prefix + '.paired.files'
#logger.info("Output directory: {}".format(in_dir))

# Note: this outputs characters to STDOUT; need to somehow make them go away.
