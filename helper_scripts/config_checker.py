#!/usr/bin/env python
import logging, os, sys
from configparser import ConfigParser

def dirFileExists(config,section_name,option_name):
    """
    Test if the option_name is in 'section_name' section and if it exists/readable
    """
    try:
        path = config[section_name][option_name]
        return os.path.exists(os.path.expanduser(path))
    except (KeyError): #option_name is not in 'file_inputs' section
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
        logging.basicConfig(filename = os.path.expanduser(config['file_inputs']['output_dir']) + '/config_checker.log', format = LOG_FORMAT, level = logging.DEBUG)
        logger = logging.getLogger()
    else:
        logging.basicConfig(filename = os.getcwd() + '/config_checker_log.log', format=LOG_FORMAT,level = logging.DEBUG)
        logger = logging.getLogger()
        logger.error(f"*** option 'output_dir' not found or 'output_dir does not exist")
        sys.exit(1)


    logger.info(f"Config file passed: {cfg_file}")
    mothur_sections = ['file_inputs', 'contigs_params', 'rename_param', 'screen_params', 'pcr_params', 'rare_seqs_param']
    checklist = []

    # check if config file has all required sections
    if set(mothur_sections) == set(config.sections()):
        logger.info("All required sections in the config file were found.")
    else:
        for section in mothur_sections:
            pair = (section, config.has_section(section))
            checklist.append(pair[1])
            if pair[1] == True:
                logger.info(f"Section found: {section}")
            else:
                logger.error(f"*** Section not found: {section}")
   
    # Check if these required options are in 'file_inputs' section and they are readable
    existList = ['batch_file','oligos','input_dir', 'output_dir']
    for name in existList:
        if dirFileExists(config,'file_inputs',name):
            logger.info(f"{config.get('file_inputs',name)} exists and is readable")
            checklist.append(True)
        else:
            logger.error(f"*** {name} in config file does not exist or is not readable")
            checklist.append(False)
            sys.exit(1)

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
                logger.error(f"*** {val} is not an integer")
                checklist.append(False)

    # Check if option 'prefix' exists
    if (config.has_option('rename_param','prefix')):
        checklist.append(True)
    else:
        checklist.append(False)
        logger.error("*** prefix not found in rename_param")

    # enforce proper batch file format
    # it must have R1, R2 reads and at least I1 index file
    # in case of missing I2 index file, must have NONE or none in place
    missing_i2_flag = False
    with open(os.path.expanduser(config['file_inputs']['batch_file'])) as f:
        for line in f.readlines():
            if not all(x in line for x in ["R1", "R2", "I1"]):
                logger.error(f"*** You must specify both an R1,R2 and at least I1 file. Check all rows of your batch file")
                checklist.append(False)
            if any(x in line for x in ['none', 'NONE']):
                missing_i2_flag = True
                logger.info(f"batch file has 'NONE' in place of missing I2 index file")
            if not any(x in line for x in ['I2', 'NONE', 'none']):
                logger.error(f"*** batch file has error. If we don't have I2 index file, "
                             f"we must have keyword 'NONE' in place.")
                checklist.append(False)

    # check for evil characters: '-' and ',' in oligo file
    # MOTHUR will not accept special characters like '-' or ','
    with open(os.path.expanduser(config['file_inputs']['oligos'])) as f:
        for line in f.readlines():
            if any(evil in line for evil in ["-",","]):
                logger.error(f'*** Either hyphens or commas found in oligo file. '
                             f'Consider changing hyphens to underscores. \n'
                             f'*** If oligo file is comma-delimited, please update it to tab-delimited')
                checklist.append(False)
                break

    # enforce proper oligo file format
    # in case of missing I2 index, use keyword NONE in the place of all I2 index
    if (missing_i2_flag):
        with open(os.path.expanduser(config['file_inputs']['oligos'])) as f:
            if all('NONE' in line.upper() for line in f.readlines() if 'barcode' in line):
                logger.info(f"oligo file also has NONE in place of all missing I2 index")
            else:
                logger.error(f"*** oligo file has error, must have keyword 'NONE' in place of all missing I2 index")
                checklist.append(False)

    if sum(checklist) == len(checklist):
        logger.info("Config file checking completed.")
    else:
        print(f"Found errors in config file. Please correct the file and run the pipeline again.")
        logger.error("*** Config file has errors. Please correct those errors and run the pipeline again.")
        sys.exit(1)
    
    return config


if __name__ == "__main__":
    print("This module is called by pipeline.py.  Please run pipeline.py --help for more information")

