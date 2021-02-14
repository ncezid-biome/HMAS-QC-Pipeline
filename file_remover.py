#!/usr/bin/env python
import glob, os


def remove_vsearch_files(config):
    """
    This function parse the oligo file (in the config object), and use the primer ID/barcode ID as a pattern
    to search and remove the temp files generated from m.chimera.vsearch()
    Didn't do any checking on the passed-in config object. (assuming the calling code already checked on it)

    Parameters
    ----------
    config object

    Returns
    -------
    None

    """
    path = config['file_inputs']['oligos']
    dir = config['file_inputs']['output_dir']
    with open(path) as f:
    #parse the oligo file and use the primer_id (4th column) as
    #search pattern, to remove the temp files from m.chimera.vsearch
        for label_id in (line.split()[3] for line in f.readlines() if ('primer' in line or 'barcode' in line)):
            for file_to_remove in glob.glob(rf"{dir}/*{label_id}*"):
                os.remove(file_to_remove)

if __name__ == "__main__":
    print("This module is called by pipeline.py.  Please run pipeline.py --help for more information")