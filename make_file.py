#!/usr/bin/env python

import os, logging, argparse, re
import pandas as pd
import numpy as np

LOG_FORMAT = "%(levelname)s %(asctime)s - %(message)s"
logging.basicConfig(filename = 'make_file.log', format = LOG_FORMAT, level = logging.DEBUG)
logger = logging.getLogger()


def getFiles(prefix, suffix):
    fileList = []
    current_path = os.getcwd() 
    for f in os.listdir(current_path):
        fname, fext = os.path.splitext(f)
        if (f.startswith(prefix)) and (fext == ('.' + suffix)):
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


def checkDf(df):
    max_row, max_col = df.shape
    cols = df.columns.tolist()
    data = ['none' for i in range(max_row)]
    if len(cols) < 4:
        missing = list(set(['R1', 'R2', 'I1', 'I2']) - set(cols))
        for i in missing:
            if i == 'R1':
                df.insert(0, 'R1', data)
            elif i == 'R2':
                df.insert(1, 'R2', data)
            elif i == 'I1':
                df.insert(2, 'I1', data)
            elif i == 'I2':
                df.insert(3, 'I2', data)
    return(df)


LOG_FORMAT = "%(levelname)s %(asctime)s - %(message)s"
logging.basicConfig(filename = 'make_file.log', format = LOG_FORMAT, level = logging.DEBUG)
logger = logging.getLogger()

parser = argparse.ArgumentParser(description = 'Make a list of the files to run in batch mode')
parser.add_argument('-p', '--prefix', metavar='', required=True, help = 'A prefix that identifies the files you want to run')
parser.add_argument('-x', '--extension', metavar='', required=True, help = 'An extension that your files have, e.g. .gz (NOT .gz)')

args = parser.parse_args()

prefix = args.prefix
suffix = args.extension

logger.info("Arguments passed: prefix = {0}, extension = {1}".format(prefix, suffix))
logger.info('Retrieving files that start with {0} and end with {1}'.format(prefix, suffix))


filetypes = ['R1', 'R2', 'I1', 'I2']
order = {key: i for i, key in enumerate(filetypes)}

if __name__ == '__main__':
    fileList = getFiles(prefix, suffix)
    fileTable = makeTable(fileList)
    finalTable = fileTable.reset_index(drop = True) \
       			  .pivot(index = 'shortname', columns = 'type', values = 'filename') \
       			  .sort_index(axis = 1, key = lambda x: x.map(order))
    outputTable = checkDf(finalTable)

data = outputTable.values
logger.info('Added these files to HMAS QC batch file: {0}'.format(data))
outfile = prefix + '.paired.files'
logger.info('Saving output to: {0}'.format(outfile))
np.savetxt(outfile, data, delimiter = '\t', fmt = '%s')

logger.info('Completed. Please make sure {0} is correct before continuing.'.format(outfile)) 

