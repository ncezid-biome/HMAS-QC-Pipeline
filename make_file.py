#!/usr/bin/env python

import os, logging, argparse, re
import pandas as pd

# Take user inputs (filename prefix, extension) and make a list of the files in current wd
# Format should be: forward \t reverse \t fwd_index \t None


def getFiles(prefix, suffix):

    logger.info("Current input directory: {}".format(os.getcwd()))
    fileList = []
    for f in os.listdir():
        fname, fext = os.path.splitext(f)
        if (f.startswith(prefix)) and (fext == suffix):
            logger.info("File added: {}".format(fname))
            fileList.append(f)
    return(fileList)

def makeTable(fileList):

    fileNestedList = []
    fileTable = pd.DataFrame()
    rule = re.compile('[RI][1-2]')
    for f in fileList:
        fname, fext = os.path.splitext(f)
        key = rule.search(f)
        short = f.replace(key.group(),"")
        row = {'filename' : f,
                'type' : key.group(),
                'shortname' : short}
        fileNestedList.append(row)
        fileTable = pd.DataFrame(fileNestedList)
    return(fileTable)


# To add: function that arranges fileTable in the proper order
# To add: function that outputs the final dataframe
    #outFile = prefix + '.paired.files'
    #logger.info("Output directory: {}".format(in_dir))


LOG_FORMAT = "%(levelname)s %(asctime)s - %(message)s"
logging.basicConfig(filename = 'make_file.log', format = LOG_FORMAT, level = logging.DEBUG)
logger = logging.getLogger()

parser = argparse.ArgumentParser(description = 'Make a list of the files to run in batch mode')
parser.add_argument('-p', '--prefix', metavar='', required=True, help = 'A prefix that identifies the files you want to run')
parser.add_argument('-x', '--extension', metavar='', required=True, help = 'An extension that your files have, e.g. .gz')

# Note: potential problem: user puts .fastq.gz instead of .gz, has a mix of files in directory
# Note: what if the user's files don't have a prefix? 

args = parser.parse_args()

prefix = args.prefix
suffix = args.extension

logger.info("Arguments passed: prefix = {}, extension = {}".format(prefix, extension))


if __name__ == '__main__':
    fileList = getFiles(prefix, suffix)
    fileTable = makeTable(fileList)
    #finalFileTable = makeFinalTable(fileTable)
    #makeOutfile(finalFileTable)
