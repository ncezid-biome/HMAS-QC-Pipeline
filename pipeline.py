#!/usr/bin/env python

import logging, sys, os, argparse, shutil, errno
import mpy_batch, config_checker
import re

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

def getErrorCode(error_message):
    """parse the error message and fetch the error code

    Parameters
    ----------
    error_message

    Returns
    -------
    either None or an integer value

    """
    return_code = re.search(r"return_code=\S+", error_message)[0][12:]
    if (return_code == 'None'):
        return None
    else:
        return int(return_code)

def lookUpErrorCode(error_code):
    """
    Parameters
    ----------
    error_code

    Returns
    -------
    the description of the error code (POSIX-UNIX system).
    Note the error code is system-dependent.
    """
    signal_table = {1:'terminate a connection, or reload the configuration for daemons',
                    2:'interrupt the session from the dialogue station',
                    3:'terminate the session from the dialogue station',
                    4:'illegal instruction was executed',
                    5:'Trace/breakpoint trap',
                    6:'abnormal termination',
                    7:'error on the system bus',
                    8:'Erroneous arithmetic operation',
                    9:'immediately terminate the process',
                    10:'user-defined signal',
                    11:'segmentation fault due to illegal access of a memory segment',
                    12:'user-defined signal',
                    13:'writing into a pipe, and nobody is reading from it',
                    14:'the timer terminated (alarm)',
                    15:'terminate the process in a soft way'}

    if (error_code): # error_code = 0 or None will be evaluated as False
        if (abs(error_code) in signal_table): # error code is usually negative
            return (signal_table[abs(error_code)])
        else:
            return ("unknown error code, it's not an unsigned 8 bits value")
    else:
        return ("There is either no error (code =0), or error_code = None")

def main():

    parser = argparse.ArgumentParser(description = 'Run Mothur QC pipeline on HMAS data.')
    parser.add_argument('-c', '--config', metavar = '', required = True, help = 'Specify configuration file')
    parser.add_argument('-l', '--log', metavar='', required=False, help='Suppress mothur log, '
                        'pipe stdout to the file you provided')
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
    

    # Check input files
    error = 0 
    try:
    	with open(config['file_inputs']['batch_file']) as f:
            for line in f.readlines():
                if all(x in line for x in ["R1", "R2"]) == False:
                    logger.error('You must specify both an R1 and R2 file. Check check all rows of your batch file.')
                    error = 1
                # if any(x in line for x in ["I1", "I2"]) == False:
                #     logger.error('You must specify at least one index file (preferably both I1 and I2). Check your batch file.')
                #     error = 1
                if any("-" in f for f in line.split()) == True:
                    logger.error('Please remove all hyphens from your file names. Consider changing them to underscores.')
                    error = 1
    except OSError as e:
        print(f'{e}')
        logger.error(e)

    finally:
        f.close()

    if error == 1:
        logger.error('We have encountered errors in your batch file.  Please correct them and run the pipeline again.')
        sys.exit(1)
    else:
        logger.info('Both read files and at least one index file found for all inputs in batch file.')
        logger.info('None of the files contain evil hyphens.')

    # Catch potential RuntimeError thrown by mothur-py and log the error code
    try:
        if (args.log is not None):
            mpy_batch.SUPPRESS_LOG_FLAG = True
            sys.stdout = open(args.log, 'w')
        mpy_batch.main(config)
        logger.info(f'mothur_py executed on files listed in {args.config}')
    except RuntimeError as e:
        print(f'{e}')
        print('Please also check mothur logfile for details')
        logger.error(e)
        logger.error(f'The return error code might indicate: {lookUpErrorCode(getErrorCode(str(e)))}')
        logger.error('Please also check mothur log file for error details')
    finally:
        if (args.log is not None):
            sys.stdout.close()

if __name__ == "__main__":
    main()

