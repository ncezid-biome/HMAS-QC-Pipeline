#!/usr/bin/env python

import logging

LOG_FORMAT = "%(levelname)s %(asctime)s - %(message)s"
logging.basicConfig(filename = 'make_file.log', format = LOG_FORMAT, level = logging.DEBUG)
logger = logging.getLogger()

def checkFile(filename):
    """
    Params
    ------

    Returns
    ------

    """
    try:
        f = open(filename)
        f.close()
    except FileNotFoundError:
        print('File not found, did you remember to create it?')



def main(args):
    
    # Take the arg from the main script and set it here 
    cfg_file = args.config  

    from configparser import ConfigParser

    config = ConfigParser()
    config.read(cfg_file)
    logger.info('Config file passed: {}'.format(cfg_file))

    mothur_sections = ['file_inputs', 'contigs_params', 'rename_param', 'screen_params', 'pcr_params', 'rare_seqs_param']

    # This script should also check for presence of oligos file and the batch file
    # This doesn't work how I expect it to...
    if set(mothur_sections) != set(config.sections()):
        for section in mothur_section:
            try:
                config.has_section(section)
            except config.Error:
                logger.error('This section of your config file is missing: {}'.format(section))
    
    # Export config object back to main script
    return config
