#!/usr/bin/env python

import logging, sys, os, argparse, shutils
import mpy_batch, config_checker

# Parse user arguments (which specify the config file)
parser = argparse.ArgumentParser(description = 'Run in batch mode or run a single sample.')
parser.add_argument('-c', '--config', metavar = '', help = 'Specify configuration file')

args = parser.parse_args()

# Send the specified config file to the config_checker program
# Run the config_checker program and return
config_checker.main()


# Check that the input files exist
for file in (read1, read2, index1, index2, oligosfile):
    if not os.path.exists(file):
        missing_file_flag = 1
        print('File {0} not found'.format(file))
        
if missing_file_flag == 1:
    print("Program exited because an input file was not found", sys.stderr)
    sys.exit(1)

# Check explicitly that Mothur is on the path (mothur_py will also check)
# Probably want to do this in a nicer way - when I update the other checks
def find_tool(name):
    return shutil.which(name) is not None

if find_tool("mothur") == False:
        print('{0} not found on path. Is it installed?'.format(name), sys.stderr)
        sys.exit(1)

# Import Mothur
try:
    from mothur_py import Mothur
    m = Mothur()
except:
    print("Unable to import mothur_py module.  Is it installed and on the path?")
    print("Program exited because mothur_py could not be imported", sys.stderr)
    sys.exit(1)
    
#### Set directories for the raw files (mothur_py should dump intermediate files in the current folder?)
    # Decide how to work out the directories

mpy_batch.main()
    

