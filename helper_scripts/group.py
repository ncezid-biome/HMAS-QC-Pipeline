#!/usr/bin/env python
import re
import os
import pandas as pd

def get_current_file(mothur_log_file, file_type='group'):
    """
    This utility function parses the MOTHUR log file, and retrieve the most current file for the specified type

    Returns
    -------
    the most current file name
    """
    current_file = ''
    if os.access(mothur_log_file, os.R_OK):
        with open(mothur_log_file, 'r', errors='ignore') as f:
            log_file = f.read()
        current_file_s = re.finditer(rf"{file_type}=.*{file_type}[^\s]*", log_file)
        for current_file in current_file_s:  # we want to get the last find
            pass
        if not current_file:  # in case MOTHUR log file doesn't have group file
            raise RuntimeError(f"************************************************************\n"
                               f"return_code=None   Alert! ERROR\n"
                               f"can't find {file_type} file in MOTHUR log\n"
                               f"************************************************************")
    return current_file.group(0)[len(file_type)+1:]  # strip the 'file_type='


def create_new_Group(oldgroup, label='barcode'):
    """
    This function parse the oligo file (in the config object), and use the 'label' as a pattern
    to search the old_group file and create new group file
    Didn't do any checking on the passed-in config object. (assuming the calling code already checked on it)

    Parameters
    ----------
    config object, old_group file name, label string

    Returns
    -------
    new_group file name

    """

    df = pd.read_csv(oldgroup, sep='\t', header=None, names=['seq','group'])
    if label == 'barcode':
        df['group'] = df['group'].str.split('.').str[0]
    else:
        df['group'] = df['group'].str.split('.').str[1]

    df.to_csv(f"{oldgroup}_new_{label}", sep='\t', header=None, index=None)
    return (f"{oldgroup}_new_{label}")



if __name__ == "__main__":
    print("This module is called by pipeline.py.  Please run pipeline.py --help for more information")