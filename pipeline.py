#!/usr/bin/env python

import logging, sys, os, argparse, shutil
import mpy_batch, config_checker


def find_tool(name):
    """Checks PATH for existence of an executable
    
    Params
    ------
    name: String
        Name of the executable file

    Returns
    ------
    True/False: Boolean
        True if `name` is on path and executable, False otherwise   
    """
    return shutil.which(name) is not None


def main():
    
    parser = argparse.ArgumentParser(description = 'Run Mothur QC pipeline on HMAS data.')
    parser.add_argument('-c', '--config', metavar = '', help = 'Specify configuration file')
    args = parser.parse_args()
    
    config = config_checker.main(args) 

    if find_tool("mothur") == False:
        print('{0} not found on path. Is it installed?'.format(name), sys.stderr)
        sys.exit(1)

    try:
        from mothur_py import Mothur
    except:
        print("Unable to import mothur_py module.  Is it installed and on the path?")
        print("Program exited because mothur_py could not be imported", sys.stderr)
        sys.exit(1)
    
    mpy_batch.main(config)
    

if __name__ == "__main__":
    main()

