#!/usr/bin/env python
import re
import os


def get_current_group(mothur_log_file):
    """
    This utility function parses the MOTHUR log file, and retrieve the most current group file

    Returns
    -------
    the most current group file name
    """
    current_group = ''
    if os.access(mothur_log_file, os.R_OK):
        with open(mothur_log_file, 'r', errors='ignore') as f:
            log_file = f.read()
        current_groups = re.finditer(r'group=.*', log_file)
        for current_group in current_groups:  # we want to get the last find
            pass
        if not current_group:  # in case MOTHUR log file doesn't have group file
            raise RuntimeError(f"************************************************************\n"
                               f"return_code=None   Alert! ERROR\n"
                               f"can't find group file in MOTHUR log\n"
                               f"************************************************************")
    return current_group.group(0)[6:]  # strip the 'group='


def create_new_Group(config, oldgroup):
    """
    This function parse the oligo file (in the config object), and use the barcode_label as a pattern
    to search the old_group file and create new group file
    Didn't do any checking on the passed-in config object. (assuming the calling code already checked on it)

    Parameters
    ----------
    config object, old_group file name

    Returns
    -------
    new_group file name

    """
    path = config['file_inputs']['oligos']
    # path = config
    with open(path) as o, open(oldgroup) as g:
        # parse the oligo file and use the barcode_label (4th column) as
        # search pattern, to remove the temp files from m.chimera.vsearch
        oldgroup_file = g.read()
        for barcode_label in (line.split()[3] for line in o.readlines() if 'barcode' in line):
            oldgroup_file = re.sub(rf'{barcode_label}.*', barcode_label, oldgroup_file)

    with open(oldgroup+'_new', 'w') as f:
        f.write(oldgroup_file)

    return (oldgroup+'_new')

if __name__ == "__main__":
    print("This module is called by pipeline.py.  Please run pipeline.py --help for more information")