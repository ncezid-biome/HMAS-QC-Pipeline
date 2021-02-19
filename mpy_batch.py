#!/usr/bin/env python
import os
from datetime import datetime
from mothur_py import Mothur
import re
import group
import glob

MOTHUR_LOG_FILE = ''


def check_filesize():
    """
    This utility function checks and makes sure the good fasta file size is larger than scrap fasta file size
    , otherwise it will throw out an RuntimeError
    It does so by parsing the MOTHUR log file, and grab those 2 fasta file names, then compare their sizes.

    """
    if (os.access(MOTHUR_LOG_FILE, os.R_OK)):
        with open(MOTHUR_LOG_FILE, 'r', errors='ignore') as f:
            log_file = f.read()
        try:
            good_fasta_file = re.search(r'.+trim\.contigs\.fasta', log_file).group(0)
            scrap_fasta_file = re.search(r'.+scrap\.contigs\.fasta',log_file).group(0)
        except AttributeError: #in case MOTHUR log file doesn't have trim/scrap fasta file
            raise RuntimeError(f"************************************************************\n"
                               f"return_code=None   Alert! ERROR\n"
                               f"can't find either trim or scrap fasta file in MOTHUR log\n"
                               f"************************************************************")

        if (os.path.getsize(scrap_fasta_file) > os.path.getsize(good_fasta_file)):
            raise RuntimeError(f"************************************************************\n"
                               f"return_code=None   Alert! ERROR\n"
                               f"We have scrap fasta {scrap_fasta_file} file size larger than \n"
                               f"good fasta {good_fasta_file}!\n"
                               f"You most likely need to trim your reads first!\n"
                               f"************************************************************")


def main(config):

    global MOTHUR_LOG_FILE
    currentDir = os.getcwd()
    now = datetime.now()
    runDateTime = now.strftime("%Y%m%d%H%M%S")
    MOTHUR_LOG_FILE = config.get('file_inputs', 'output_dir', fallback = currentDir)+'/mothur.'+runDateTime+'.logfile'
    m = Mothur(logfile_name = MOTHUR_LOG_FILE)

    m.set.dir(input = config.get('file_inputs', 'input_dir', fallback = currentDir),
              output = config.get('file_inputs', 'output_dir', fallback = currentDir),
              tempdefault = config.get('file_inputs','output_dir', fallback = currentDir))
    m.make.contigs(file = config.get('file_inputs','batch_file'),
                    processors= config.getint('contigs_params', 'processors', fallback = 40), 
                    format=config.get('contigs_params', 'format', fallback = 'illumina1.8+'), 
                    oligos=config.get('file_inputs', 'oligos'), 
                    bdiffs=config.getint('contigs_params', 'bdiffs', fallback = 0), 
                    pdiffs=config.getint('contigs_params', 'pdiffs', fallback = 0), 
                    checkorient=config.get('contigs_params', 'checkorient', fallback = 't'), 
                    insert=config.getint('contigs_params', 'insert', fallback = 25), 
                    trimoverlap=config.get('contigs_params', 'trimoverlap', fallback = 't'),
                    allfiles=config.get('contigs_params','allfiles', fallback = 0))

    # this serves as a check point and make sure that:
    # good fasta file (trim.contigs.fasta) size > scrap fasta file size
    check_filesize()

    m.summary.seqs(fasta='current', name='current')

    m.rename.file(fasta='current', group='current',
                    prefix=config.get('rename_param', 'prefix'))

    m.screen.seqs(fasta='current', group='current',
                    maxambig=config.getint('screen_params', 'maxambig', fallback = 0),
                    maxlength=config.getint('screen_params', 'maxlength', fallback = 325))
    m.summary.seqs(fasta='current', name='current')

    m.unique.seqs(fasta='current')
    m.count.seqs(name='current', group='current')
    m.summary.seqs(fasta='current', name='current')

    m.get.current() # let MOTHUR log the most current group file
    old_group = group.get_current_group(MOTHUR_LOG_FILE)
    new_group = group.create_new_Group(config, old_group)
    old_accnos = group.get_current_accnos(MOTHUR_LOG_FILE)

    dir = config['file_inputs']['output_dir']
    for accnos_file in glob.glob(rf"{dir}/*.accnos"):
        if os.path.getsize(accnos_file) == 0: # we shouldn't find any accnos file at this stage, but we check anyway
            os.remove(accnos_file)

    # Search for chimeras
    # need to use proper path for vsearch
    m.chimera.vsearch(fasta='current', name='current', group=new_group, dereplicate='t', vsearch=r'~/bin/vsearch')

    # check for emptry accnos file
    # because we checked already, if we find another accnos file, it must come from chimera.vsearch()
    chimera_accnos_empty = False
    for accnos_file in glob.glob(rf"{dir}/*.accnos"):
        if os.path.getsize(accnos_file) == 0:
            os.remove(accnos_file)
            chimera_accnos_empty = True

    if chimera_accnos_empty:
        m.set.current(group=old_group, accnos=old_accnos)
    else:
        m.set.current(group=old_group)
        m.remove.seqs(fasta='current', accnos='current', group='current', name='current')

    m.cluster(count='current', method='unique', cutoff='unique')

    m.remove.rare(list='current', count='current',
                      nseqs=config.getint('rare_seqs_param', 'nseqs', fallback=9),
                      label='unique')

    m.summary.seqs(fasta='current', name='current')

    m.list.seqs(count='current')

    m.get.seqs(fasta='current', accnos='current', name='current', group='current')

    m.rename.file(fasta='current', group='current', name='current', prefix=config.get('rename_param', 'prefix') +'.final')

if __name__ == "__main__":
    print("This module is called by pipeline.py.  Please run pipeline.py --help for more information")
