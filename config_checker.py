#!/usr/bin/env python

def main(cfg_file):
    
    import logging, os, sys
    from configparser import ConfigParser

    config = ConfigParser()
    config.read(cfg_file)

    LOG_FORMAT = "%(levelname)s %(asctime)s - %(message)s"
    logging.basicConfig(filename = config['file_inputs']['output_dir'] + '/config_checker.log', format = LOG_FORMAT, level = logging.DEBUG)
    logger = logging.getLogger()

    logger.info('Config file passed: {}'.format(cfg_file))
    mothur_sections = ['file_inputs', 'contigs_params', 'rename_param', 'screen_params', 'pcr_params', 'rare_seqs_param']
    checklist = []

    if set(mothur_sections) == set(config.sections()):
        logger.info("All required sections in the config file were found.")
    else:
        for section in mothur_sections:
            pair = (section, config.has_section(section))
            checklist.append(pair[1])
            if pair[1] == True:
                logger.info('Section found: {}'.format(section))
            else:
                logger.error('Section not found: {}'.format(section))
   
    if os.access(config['file_inputs']['oligos'], os.R_OK) == True:
        logger.info('{} exists and is readable'.format(config['file_inputs']['oligos']))
        checklist.append(True)
    else:
        logger.error('{} does not exist or is not readable'.format(config['file_inputs']['oligos']))  
        checklist.append(False)

    if os.access(config['file_inputs']['batch_file'], os.R_OK) == True:
        logger.info('{} exists and is readable'.format(config['file_inputs']['batch_file']))
        checklist.append(True)
    else:    
        logger.error('{} does not exist or is not readable'.format(config['file_inputs']['batch_file']))
        checklist.append(False)

    if sum(checklist) == len(checklist):
        logger.info("Config file checking completed.")
    else:
        logger.error("Config file has errors. Please correct the file and run the pipeline again.")
        sys.exit(1)
    
    return config

if __name__ == "__main__":
    print("This module is called by pipeline.py.  Please run pipeline.py --help for more information")

