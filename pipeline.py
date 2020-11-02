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
    found = shutil.which(name) is not None
    return(found)


def main():
    
    parser = argparse.ArgumentParser(description = 'Run Mothur QC pipeline on HMAS data.')
    parser.add_argument('-c', '--config', metavar = '', required = True, help = 'Specify configuration file')
    args = parser.parse_args()    

    cfg_file = args.config
    config = config_checker.main(cfg_file) 

    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    LOG_FORMAT = "%(levelname)s %(asctime)s - %(message)s"
    logging.basicConfig(filename = config['file_inputs']['output_dir'] + '/hmas_qc_pipeline.log', format = LOG_FORMAT, level = logging.DEBUG)
    logger = logging.getLogger()

    logger.info('The config file to be parsed is: {0}'.format(args.config))

    if find_tool('mothur') == True:
        logger.info('mothur is on path and is executable.')
    else:
        logger.error('mothur not found on path. Is it installed?')
        sys.exit(1)

    try:
        from mothur_py import Mothur
    except:
        logger.error('Unable to import mothur_py module. Is it installed and on PATH?')
        logger.error('Program exited because mothur_py could not be imported.')  
        sys.exit(1)
    
    mpy_batch.main(config)
    logger.info('mothur_py executed on files listed in {}'.format(args.config))

if __name__ == "__main__":
    main()

