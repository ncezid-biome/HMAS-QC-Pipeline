#!/usr/bin/env python
import re
import os
import pandas as pd
import concurrent.futures
import logging

dict_barcode = {}
dict_primer = {}
dict_R1 = {}
dict_R2 = {}
df_I1 = pd.DataFrame()
df_I2 = pd.DataFrame()
has_I2_index = False


def revcomp(myseq):
    rc = {'A' : 'T', 'T' : 'A', 'G' : 'C', 'C' : 'G'}
    seq = [rc[n] if n in rc else n for n in myseq] # allow IUPAC code to stay as is
    return("".join(list(reversed(seq))))

class Primers:
    'Parses mothur oligos file and extracts the primer names. Stole from A. Jo Williams-Newkirk'
    
    def __init__(self, fname):
        self.fname = fname
        self.pseqs = dict()
        Primers.reader(self, self.fname, self.pseqs)
        self.pnames = self.pseqs.keys()
    
    def reader(self, fname, pseqs):
        with open(fname, 'r') as infile:
            for line in infile:
                if line.startswith("primer"):
                    tmp = line.split('\t')
                    pseqs[tmp[3].strip('\n')] = [tmp[1], revcomp(tmp[2])]

def get_current_group(mothur_log_file):
    """
    This utility function parses the MOTHUR log file, and retrieve the most current group file

    Returns
    -------
    the most current group file name
    """
    current_group = ''
    if os.access(mothur_log_file, os.R_OK):
        with open(mothur_log_file, 'r', errors='ignore') as f:
            log_file = f.read()
        current_groups = re.finditer(r'group=.*', log_file)
        for current_group in current_groups:  # we want to get the last find
            pass
        if not current_group:  # in case MOTHUR log file doesn't have group file
            raise RuntimeError(f"************************************************************\n"
                               f"return_code=None   Alert! ERROR\n"
                               f"can't find group file in MOTHUR log\n"
                               f"************************************************************")
    return current_group.group(0)[6:]  # strip the 'group='

def get_current_accnos(mothur_log_file):
    """
    This utility function parses the MOTHUR log file, and retrieve the most current accnos file

    Returns
    -------
    the most current group file name
    """
    current_accnos = ''
    if os.access(mothur_log_file, os.R_OK):
        with open(mothur_log_file, 'r', errors='ignore') as f:
            log_file = f.read()
        current_accnos_s = re.finditer(r'accnos=.*', log_file)
        for current_accnos in current_accnos_s:  # we want to get the last find
            pass
        if not current_accnos:  # in case MOTHUR log file doesn't have group file
            raise RuntimeError(f"************************************************************\n"
                               f"return_code=None   Alert! ERROR\n"
                               f"can't find accnos file in MOTHUR log\n"
                               f"************************************************************")
    return current_accnos.group(0)[7:]  # strip the 'accnos='

def get_current_file(mothur_log_file, file_type='group'):
    """
    This utility function parses the MOTHUR log file, and retrieve the most current file for the specified type

    Returns
    -------
    the most current file name
    """
    current_file = ''
    if os.access(mothur_log_file, os.R_OK):
        with open(mothur_log_file, 'r', errors='ignore') as f:
            log_file = f.read()
        current_file_s = re.finditer(rf"{file_type}=.*{file_type}(s)?", log_file)
        for current_file in current_file_s:  # we want to get the last find
            pass
        if not current_file:  # in case MOTHUR log file doesn't have group file
            raise RuntimeError(f"************************************************************\n"
                               f"return_code=None   Alert! ERROR\n"
                               f"can't find {file_type} file in MOTHUR log\n"
                               f"************************************************************")
    return current_file.group(0)[len(file_type)+1:]  # strip the 'file_type='


# def create_new_Group(config, oldgroup, label='barcode'):
#     """
#     This function parse the oligo file (in the config object), and use the 'label' as a pattern
#     to search the old_group file and create new group file
#     Didn't do any checking on the passed-in config object. (assuming the calling code already checked on it)

#     Parameters
#     ----------
#     config object, old_group file name, label string

#     Returns
#     -------
#     new_group file name

#     """
#     path = config['file_inputs']['oligos']
#     # path = config
#     with open(path) as o, open(oldgroup) as g:
#         # parse the oligo file and use the 'label' (4th column) as search pattern
#         oldgroup_file = g.read()
#         for barcode_label in (line.split()[3] for line in o.readlines() if label in line.split()[0]):
#             oldgroup_file = re.sub(rf'{barcode_label}.*', barcode_label, oldgroup_file)

#     with open(f"{oldgroup}_new_{label}", 'w') as f:
#         f.write(oldgroup_file)

#     return (f"{oldgroup}_new_{label}")

def create_new_Group(oldgroup, label='barcode'):
    """
    This function parse the oligo file (in the config object), and use the 'label' as a pattern
    to search the old_group file and create new group file
    Didn't do any checking on the passed-in config object. (assuming the calling code already checked on it)

    Parameters
    ----------
    config object, old_group file name, label string

    Returns
    -------
    new_group file name

    """

    df = pd.read_csv(oldgroup, sep='\t', header=None, names=['seq','group'])
    if label == 'barcode':
        df['group'] = df['group'].str.split('.').str[0]
    else:
        df['group'] = df['group'].str.split('.').str[1]

    df.to_csv(f"{oldgroup}_new_{label}", sep='\t', header=None, index=None)
    return (f"{oldgroup}_new_{label}")


def create_my_df(fastq_file, column_name):
    '''
    This function will parse the passed-in fastq file and generates a dataframe with only 2 columns:
    'seq_id' (for the sequence id) and the passed-in column_name (for the sequence)

    Parameters
    ----------
    fastq_file: one of the R1, R2, I1, I2 files
    column_name: to hold the sequence (usually 'R1','R2','I1','I2')

    Returns
    -------
    the dataframe generated

    '''

    df_temp = pd.read_csv(fastq_file, header=None)
    df = pd.DataFrame(df_temp.iloc[::4, :].values, columns=['seq_id'])
    df['seq_id'] = df['seq_id'].apply(lambda x: x.split()[0].replace(':', '_')[1:])
    df[column_name] = pd.DataFrame(df_temp.iloc[1::4, :].values)  # store the actual sequence

    del df_temp  # release memory
    return df


def match_primer(bar_key):
    '''
    This function search for the matching primers (on the R1/R2 reads), based on the barcode
    1. from the barcode, search I1/I2 index files for all matching seq_ids
    2. from the matching seq_ids, retrieve all their sequences (from dict_R1/dict_R2)
    3. foreach seq_id, search from all primers, see if there is a match
       (primer matches at the start of the sequence)

    Parameters
    ----------
    bar_key: the key(barcode) for the dict_barcode

    Returns
    -------
    A tuple of 2 lists: matched seq_id list, matched (sample_name.primer_id) list

    '''
    logger = logging.getLogger()
    col_1 = []
    col_2 = []

    if not has_I2_index:
        barcode = bar_key[0]
        barcode_label = dict_barcode[bar_key]
        filter_I1 = (df_I1['I1'] == barcode)
        seq_id_list = df_I1.loc[filter_I1, 'seq_id'].values
    else:
        f_barcode = bar_key[0]
        r_barcode = bar_key[1]
        barcode_label = dict_barcode[bar_key]
        filter_I1 = (df_I1['I1'] == f_barcode)
        filter_I2 = (df_I1['I2'] == r_barcode)
        seq_id_list = df_I1.loc[(filter_I1 & filter_I2), 'seq_id'].values

    logger.info(f'.....{barcode_label} has matched ......{len(seq_id_list)} sequences')

    for seq_id in seq_id_list:
        seq_R1 = dict_R1[seq_id]
        seq_R2 = dict_R2[seq_id]

        for primer_key in dict_primer:
            fp = primer_key[0]
            rp = primer_key[1]

            # if seq_R1.startswith(fp) and seq_R2.startswith(rp):
            if check_primer(seq_R1, fp) and check_primer(seq_R2, rp):
                col_1.append(seq_id)
                col_2.append(barcode_label + '.' + dict_primer[primer_key])
                break  # found the match, move on to the next seq_id

    return (col_1, col_2)


def check_primer(seq, primer):
    '''
    This function checks if the passed in sequence starts with the primer

    Parameters
    ----------
    seq: the sequence (from R1, R2)
    primer: the primer to check with (from oligo file)

    Returns
    -------
    True (if matched) or False otherwise

    '''
    IUPAC_dict = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
                  'Y': 'CT', 'R': 'AG', 'W': 'AT', 'S': 'CG', 'K': 'TG', 'M': 'AC',
                  'D': 'AGT', 'V': 'ACG', 'H': 'ACT','B': 'CGT', 'N': 'ACGT'}
    for index in range(len(primer)):
        bases = IUPAC_dict[primer[index]]
        if all(seq[index] != base for base in bases):
            return False

    return True


def create_group(config, group_file_name):
    '''
    This function creates the exact group file which MOTHUR generates after make.contigs()
    1. It first created 4 dictionaries to hold info for primers, barcodes, R1 and R2 fastq files
       and 2 dataframes (1 if there is no I2) for I1 and I2 index files
    2. Then it fires off multi-processes to call match_primer() to search for matching results

    Parameters
    ----------
    config object
    group_file_name (the name of the group file we want to create)

    Returns None
    -------

    '''
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    LOG_FORMAT = "%(levelname)s %(asctime)s - %(message)s"
    FILE_NAME = "group.log"
    logging.basicConfig(filename=FILE_NAME, format=LOG_FORMAT, level=logging.INFO)
    logger = logging.getLogger()

    # assuming config file is already checked !!
    with open(config['file_inputs']['batch_file']) as f:
        R1,R2,I1,I2 = f.readlines()[0].split()

    input_dir = config['file_inputs']['input_dir'] + '/'
    R1 = input_dir + R1
    R2 = input_dir + R2
    I1 = input_dir + I1

    oligos = config['file_inputs']['oligos']

    global dict_barcode
    global dict_primer
    global dict_R1
    global dict_R2
    global df_I1
    global df_I2
    global has_I2_index
    if I2.upper() != 'NONE':
        has_I2_index = True
        I2 = input_dir + I2

    # oligo file has 4 columns:
    # primer f_primer r_primer primer_id
    # barcode f_barcode r_barcode barcode_id
    with open(oligos) as o:
        for line in o.readlines():
            if 'barcode' in line:
                fb = line.split()[1]
                rb = line.split()[2]
                label = line.split()[3]
                # if rb.upper() != 'NONE':
                if I2.upper() != 'NONE':
                    dict_barcode[(fb, rb)] = label
                else:
                    dict_barcode[(fb,)] = label
            else:  # it's 'primer' then
                fp = line.split()[1]
                rp = line.split()[2]
                primer_id = line.split()[3]
                dict_primer[(fp, rp)] = primer_id

    line_counter = 0
    logger.info(f'before creating R1 dict')
    # fastq file: we only need 1st and 2nd lines
    with open(R1) as f:
        for line in f.readlines():
            if line_counter % 4 == 0:
                seq_id = line.split()[0].replace(':', '_')[1:]  # to get the part of seq_id we need
            elif line_counter % 4 == 1:  # this is the sequence
                dict_R1[seq_id] = line
            line_counter += 1

    line_counter = 0
    logger.info(f'before creating R2 dict')
    with open(R2) as f:
        for line in f.readlines():
            if line_counter % 4 == 0:
                seq_id = line.split()[0].replace(':', '_')[1:]
            elif line_counter % 4 == 1:
                dict_R2[seq_id] = line
            line_counter += 1

    logger.info(f'before creating I1 DF')
    df_I1 = create_my_df(I1, 'I1')
    logger.info(f'Done creating I1 DF')

    if I2.upper() != 'NONE':
        df_I2 = create_my_df(I2, 'I2')
        # merge 2 DFs, so df_I1 now has 3 columns ['seq_id', 'I1', 'I2']
        df_I1 = df_I1.merge(df_I2, on='seq_id', how='inner')
        del df_I2  # release memory

    # now we're done with file readings, and all the pre-processing
    col1 = []  # to hold found sequences
    col2 = []  # to hold sample_primer

    # fire off multi-processes to search for seq:sample_primer pairs
    with concurrent.futures.ProcessPoolExecutor(max_workers=48) as executor:

        results = executor.map(match_primer, dict_barcode.keys())
        for result in results:
            col1.extend(result[0])
            col2.extend(result[1])

    # now we're ready to write out the groups
    logger.info(f'before creating df_group')
    df_group = pd.DataFrame({'seq_id': col1, 'bar_prim': col2})

    logger.info(f'before sorting groups')
    df_group.sort_values(by='seq_id', inplace=True)

    logger.info(f'after sorting groups')
    df_group.to_csv(group_file_name, sep='\t', header=None, index=None)

    logger.info(f'Done writing out groups')


if __name__ == "__main__":
    print("This module is called by pipeline.py.  Please run pipeline.py --help for more information")