#!/usr/bin/env python

import logging, os, sys

LOG_FORMAT = "%(levelname)s %(asctime)s - %(message)s"
logging.basicConfig(filename = 'config_checker.log', format = LOG_FORMAT, level = logging.DEBUG)
logger = logging.getLogger()


def main(cfg_file):
    
    from configparser import ConfigParser

    config = ConfigParser()
    config.read(cfg_file)
    logger.info('Config file passed: {}'.format(cfg_file))

    mothur_sections = ['file_inputs', 'contigs_params', 'rename_param', 'screen_params', 'pcr_params', 'rare_seqs_param']

    # This doesn't work how I expect it to...
    if set(mothur_sections) == set(config.sections()):
        logger.info("All required sections in the config file were found.")
    else:
        for section in mothur_sections:
            try:
                config.has_section(section)
                logger.info('Section found: {}'.format(section))
            except config.Error:
                logger.error('Section not found: {}'.format(section))
   
    if os.access(config['file_inputs']['oligos'], os.R_OK) == True:
        logger.info('{} exists and is readable'.format(config['file_inputs']['oligos']))
    else:
        logger.error('{} does not exist or is not readable'.format(config['file_inputs']['oligos']))  
        sys.exit(1)

    if os.access(config['file_inputs']['batch_file'], os.R_OK) == True:
        logger.info('{} exists and is readable'.format(config['file_inputs']['batch_file']))
    else:    
        logger.error('{} does not exist or is not readable'.format(config['file_inputs']['batch_file']))
        sys.exit(1)

    logger.info("Config file checking completed.")
    # Export config object back to main script
    return config
