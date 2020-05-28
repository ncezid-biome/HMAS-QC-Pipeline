#!/usr/bin/env python

LOG_FORMAT = "%(levelname)s %(asctime)s - %(message)s"
logging.basicConfig(filename = 'make_file.log', format = LOG_FORMAT, level = logging.DEBUG)
logger = logging.getLogger()

def main():
    
    # Take the arg from the main script and set it here 
    cfg_file = 'settings.ini' 

    # Import config setting
    from configparser import ConfigParser

    config = ConfigParser()
    config.read(cfg_file)
    logger.info('Config file passed: {}'.format(cfg_file))

    # Check that each section is present
    mothur_sections = ['file_inputs', 'contigs_params', 'rename_param', 'screen_params', 'pcr_params', 'rare_seqs_param']

    # This doesn't work how I expect it to...
    if set(mothur_sections) != set(config.sections()):
        for section in mothur_section:
            try:
                config.has_section(section)
            except config.Error:
                logger.error('This section of your config file is missing: {}'.format(section))

    # Check that each option w/i given section is present

    # Set variables (using fallback options as needed)
    input_dir = config.get('file_inputs', 'input_dir')
    file_prefix = config.get('file_inputs')('file_prefix')
    oligosfile = config.get('file_inputs')('oligos')
    numproc = config.getint('contigs_params')('processors', fallback = 40)
    form = config.get('contigs_params')('format', fallback = 'illumina1.8+')
    bdif = config.getint('contigs_params')('bdiffs', fallback = 0)
    pdif = config.getint('contigs_params')('pdiffs', fallback = 0)
    chkorient = config.get('contigs_params')('checkorient', fallback = 't')
    ins =  config.getint('contigs_params')('insert', fallback = 25)
    trimover = config.get('contigs_params')('trimoverlap', fallback = 'f')
    prefx = config.get('rename_param')('prefix')
    maxamb = config.getint('screen_params')('maxambig', fallback = 0)
    maxlg = config.getint('screen_params')('maxlength', fallback = 325)
    pdif2 = config.getint('pcr_params')('pdiffs', fallback = 0)
    rdif = config.getint('pcr_params')('rdiffs', fallback = 0)
    numseqs = config.getint('rare_seqs_param')('nseqs', fallback = 9)
        
    # Export these variables back to main script
