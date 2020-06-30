#!/usr/bin/env python

import logging, sys, os, argparse, shutil
import mpy_batch, config_checker

# Check that the input files exist
# Might not need this function anymore (b/c of config_checker)
def checkFile(filename):
    try: 
        f = open(filename)
        f.close()
    except FileNotFoundError:
        print('File not found, did you remember to create it?')

# Check explicitly that Mothur is on the path (For Jo; mothur_py will also check)
# Probably want to do this in a nicer way - when I update the other checks
def find_tool(name):
    return shutil.which(name) is not None


def main():
    # Parse user arguments (which specify the config file)
    parser = argparse.ArgumentParser(description = 'Run Mothur QC pipeline on HMAS data.')
    parser.add_argument('-c', '--config', metavar = '', help = 'Specify configuration file')
    args = parser.parse_args() # yields a Namespace object
    # Run the config_checker program and return the config object
    config = config_checker.main(args)

    # Probably don't need this because the config checker will check for existence of this?
    #checkFile(oligosfile)

    if find_tool("mothur") == False:
        print('{0} not found on path. Is it installed?'.format(name), sys.stderr)
        sys.exit(1)

    # Import Mothur
    try:
        from mothur_py import Mothur
    except:
        print("Unable to import mothur_py module.  Is it installed and on the path?")
        print("Program exited because mothur_py could not be imported", sys.stderr)
        sys.exit(1)
    
    #### Set directories for the raw files (mothur_py should dump intermediate files in the current folder?)
        # Decide how to work out the directories

    mpy_batch.main(config)
    

if __name__ == "__main__":
    main()

