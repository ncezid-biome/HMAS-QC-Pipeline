#!/usr/bin/env python

import logging, sys, os, argparse, shutil
import mpy_batch, config_checker

LOG_FORMAT = "%(levelname)s %(asctime)s - %(message)s"
logging.basicConfig(filename = 'hmas_qc_pipeline.log', format = LOG_FORMAT, level = logging.DEBUG)
logger = logging.getLogger()


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
    not_found = shutil.which(name) is None
    if not_found == True:
        logger.info('{0} not found on path. Is it installed?'.format(name))
    else:
        logger.info('{0} is on path and is executable.'.format(name))
    return shutil.which(name) is not None


def main():
    
    parser = argparse.ArgumentParser(description = 'Run Mothur QC pipeline on HMAS data.')
    parser.add_argument('-c', '--config', metavar = '', help = 'Specify configuration file')
    args = parser.parse_args()    

    config = config_checker.main(args.config) 

    logger.info('The config file to be parsed is: {0}'.format(args.config))

    if find_tool('mothur') == False:
        sys.exit(1)

    try:
        from mothur_py import Mothur
    except:
        logger.info('Unable to import mothur_py module. Is it installed and on PATH?')
        logger.info('Program exited because mothur_py could not be imported.')  
        sys.exit(1)
    
    mpy_batch.main(config)
    logger.info('mothur_py executed on files listed in {}'.format(args.config))

if __name__ == "__main__":
    main()

