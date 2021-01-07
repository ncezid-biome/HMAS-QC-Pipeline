#!/usr/bin/env python
import logging, os, sys
from configparser import ConfigParser

def dirFileExists(config,section_name,option_name):
    """
    Test if the option_name is in 'section_name' section and if it exists/readable
    """
    try:
        #path = config.get(section_name,option_name)
        path = config[section_name][option_name]
        return os.access(path, os.R_OK)
    except (KeyError): #path_name is not in 'file_inputs' section
        return False

def hasIntVal(config,section_name,option_name):
    """ check the passed-in options exist in the section and must have int values.
    if the option does not exist, we still return true, b/c there is a fallback value for them in mpy_batch.py
    """
    if config.has_option(section_name,option_name):
        try:
            config.getint(section_name,option_name)
        except (ValueError):
            return False
    return True

def main(cfg_file):

    config = ConfigParser()
    config.read(cfg_file)

    LOG_FORMAT = "%(levelname)s %(asctime)s - %(message)s"
    if (dirFileExists(config,'file_inputs','output_dir')):
        logging.basicConfig(filename = config['file_inputs']['output_dir'] + '/config_checker_log.log', format = LOG_FORMAT, level = logging.DEBUG)
        logger = logging.getLogger()
    else:
        logging.basicConfig(filename = os.getcwd() + '/config_checker_log.log', format=LOG_FORMAT,level = logging.DEBUG)
        logger = logging.getLogger()
        logger.error("option 'output_dir' not found or 'output_dir does not exist")
        sys.exit(1)


    logger.info(f"Config file passed: {cfg_file}")
    mothur_sections = ['file_inputs', 'contigs_params', 'rename_param', 'screen_params', 'pcr_params', 'rare_seqs_param']
    checklist = []

    if set(mothur_sections) == set(config.sections()):
        logger.info("All required sections in the config file were found.")
    else:
        for section in mothur_sections:
            pair = (section, config.has_section(section))
            checklist.append(pair[1])
            if pair[1] == True:
                logger.info(f"Section found: {section}")
            else:
                logger.error(f"Section not found: {section}")
   
    # Check if these files/directories are in 'file_inputs' section and if they exist and are readable
    existList = ['batch_file','oligos','input_dir', 'output_dir']
    for name in existList:
        if dirFileExists(config,'file_inputs',name):
            logger.info(f"{config.get('file_inputs',name)} exists and is readable")
            checklist.append(True)
        else:
            logger.error(f"{name} in config file does not exist or is not readable")
            checklist.append(False)

    # Check required integer value for specific options
    param_dict = {'contigs_params':['processors', 'bdiffs', 'pdiffs', 'insert'],
                  'screen_params':['maxambig', 'maxlength'],
                  'pcr_params':['pdiffs', 'rdiffs'],
                  'rare_seqs_param':['nseqs']}
    for key in param_dict:
        for val in param_dict[key]:
            if hasIntVal(config,key,val):
                checklist.append(True)
            else:
                logger.error(f"{val} is not an integer")
                checklist.append(False)

    # Check if option 'prefix' exists
    if (config.has_option('rename_param','prefix')):
        checklist.append(True)
    else:
        checklist.append(False)
        logger.error("prefix not found in rename_param")


    if sum(checklist) == len(checklist):
        logger.info("Config file checking completed.")
    else:
        logger.error("Config file has errors. Please correct the file and run the pipeline again.")
        sys.exit(1)
    
    return config


if __name__ == "__main__":
    print("This module is called by pipeline.py.  Please run pipeline.py --help for more information")

