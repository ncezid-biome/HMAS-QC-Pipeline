#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os, logging, argparse, re
import pandas as pd
import numpy as np


# In[ ]:


def getFiles(prefix, suffix):
    fileList = []
    for f in os.listdir():
        fname, fext = os.path.splitext(f)
        if (f.startswith(prefix)) and (fext == suffix):
            fileList.append(f)
    return(fileList)


# In[ ]:


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


# In[ ]:


LOG_FORMAT = "%(levelname)s %(asctime)s - %(message)s"
logging.basicConfig(filename = 'make_file.log', format = LOG_FORMAT, level = logging.DEBUG)
logger = logging.getLogger()

parser = argparse.ArgumentParser(description = 'Make a list of the files to run in batch mode')
parser.add_argument('-p', '--prefix', metavar='', required=True, help = 'A prefix that identifies the files you want to run')
parser.add_argument('-x', '--extension', metavar='', required=True, help = 'An extension that your files have, e.g. .gz')

args = parser.parse_args()

prefix = args.prefix
suffix = args.extension

logger.info("Arguments passed: prefix = {}, extension = {}".format(prefix, suffix))

if __name__ == '__main__':
    fileList = getFiles(prefix, suffix)
    fileTable = makeTable(fileList)
    gTable = fileTable.groupby('shortname')
    sortedTable = gTable.apply(lambda _df: _df.sort_values(by = 'type', key = lambda x: x.map(order)))
    finaltable = sortedtable.reset_index(drop = True).pivot(index = 'shortname', columns = 'type', values = 'filename').sort_index(axis = 1, key = lambda x: x.map(order))

data = finaltable.values
outfile = prefix + '.paired.files'
np.savetxt(outfile, data, delimiter = '\t', fmt = '%s')

