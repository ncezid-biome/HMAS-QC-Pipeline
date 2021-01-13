#!/usr/bin/env python

import os, logging, argparse, re
import pandas as pd
import numpy as np

def getFiles(directory, prefix, suffix):
    fileList = []
    for f in os.listdir(os.path.abspath(directory)):
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

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = 'Make a list of the files to run in batch mode')
    parser.add_argument('-d', '--directory', metavar='', required=True, help = 'Specify the path to your files')
    parser.add_argument('-p', '--prefix', metavar='', required=True, help = 'A prefix that identifies the files you want to run')
    parser.add_argument('-x', '--extension', metavar='', required=True, help = 'An extension that your files have, e.g. .gz (NOT .gz)')

    args = parser.parse_args()

    LOG_FORMAT = "%(levelname)s %(asctime)s - %(message)s"
    logging.basicConfig(filename = os.path.abspath(args.directory) + '/make_file.log', format = LOG_FORMAT, level = logging.DEBUG)
    logger = logging.getLogger()

    logger.info(f"Arguments passed: prefix = {args.prefix}, extension = {args.extension}")
    logger.info(f"Retrieving files that start with {args.prefix} and end with {args.extension}")

    if args.directory == '.':
        directory = os.getcwd()
    else:
        directory = args.directory

    assert os.path.exists(directory), f"Directory not found at {directory}"
    logger.error(f"Could not find the directory you specified: {directory}. Does it exist?")

    filetypes = ['R1', 'R2', 'I1', 'I2']
    order = {key: i for i, key in enumerate(filetypes)}

    fileList = getFiles(directory, args.prefix, args.extension)
    fileTable = makeTable(fileList)
    finalTable = fileTable.reset_index(drop = True) \
       			  .pivot(index = 'shortname', columns = 'type', values = 'filename') \
       			  .sort_index(axis = 1, key = lambda x: x.map(order))

    outputTable = checkDf(finalTable)

    data = outputTable.values
    logger.info(f"Added these files to HMAS QC batch file: {data}")
    outfile = os.path.abspath(args.directory) + '/' + args.prefix + '.paired.files'
    logger.info(f"Saving output to: {outfile}")
    np.savetxt(outfile, data, delimiter = '\t', fmt = '%s')

    logger.info(f"Completed. Please make sure {outfile} is correct before continuing.") 

