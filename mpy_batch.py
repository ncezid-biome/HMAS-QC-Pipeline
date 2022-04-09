#!/usr/bin/env python
import sys, os, shutil, subprocess
from datetime import datetime
from mothur_py import Mothur
import re
import group, run_cutadapt
import glob
from collections import deque
import concurrent.futures

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
    MOTHUR_LOG_FILE = os.path.expanduser(config['file_inputs']['output_dir'])+'/mothur.'+runDateTime+'.logfile'
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
                    trimoverlap=config.get('contigs_params', 'trimoverlap', fallback = 'f'),
                    allfiles=config.get('contigs_params','allfiles', fallback = 0))

    # this serves as a check point and make sure that:
    # good fasta file (trim.contigs.fasta) size > scrap fasta file size
    check_filesize()

    m.rename.file(fasta='current', group='current',
                    prefix=config.get('rename_param', 'prefix'))

    m.screen.seqs(fasta='current', group='current',
                    maxambig=config.getint('screen_params', 'maxambig', fallback = 0),
                    maxlength=config.getint('screen_params', 'maxlength', fallback = 325))

    m.unique.seqs(fasta='current')
    m.summary.seqs(fasta='current', name='current')

    print(datetime.now()) 
    old_group = group.get_current_file(MOTHUR_LOG_FILE)
    print (f"old group is: {old_group}")
    new_group = group.create_new_Group(old_group)
    print (f"new group is: {new_group}")
    print(datetime.now()) 

    old_accnos = group.get_current_file(MOTHUR_LOG_FILE, 'accnos') 
    
    old_fasta = group.get_current_file(MOTHUR_LOG_FILE, 'fasta')
    print (f"old fasta is: {old_fasta}")
    print (f"old accnos is: {old_accnos}")

    #1. create primer pair group file
    new_group_primer = group.create_new_Group(old_group, 'primer')
    print (f"new group primer is: {new_group_primer}")
    print(datetime.now()) 

    #2. count.seqs()
    m.count.seqs(name='current', group=new_group_primer)
    print(datetime.now()) 
    print(f"done with count.seqs")

    #3. split.groups()
    m.split.groups(fasta=old_fasta, count='current')
    print(datetime.now()) 
    print(f"done with split.groups")


    new_fasta = f'{old_fasta[:-5]}merged.fasta'
    run_cutadapt.run_cutadapt_mothur(config, new_fasta, 36)

    # reverse-expand ~, to avoid hypen issue in directory name
    new_fasta = new_fasta.replace(os.path.expanduser('~'), '~', 1)

    # add these to avoid potential duplicate name issue in fasta file
    m.list.seqs(fasta=new_fasta) #list all the unique names in the fasta file
    m.get.seqs(fasta='current', accnos='current') #remove any duplicate names
    #6. unique.seqs
    m.unique.seqs(fasta='current', name='current')
    m.summary.seqs(fasta='current', name='current')
    

    dir = os.path.expanduser(config['file_inputs']['output_dir'])
    for accnos_file in glob.glob(rf"{dir}/*.accnos"):
        if os.path.getsize(accnos_file) == 0: # we shouldn't find any accnos file at this stage, but we check anyway
            os.remove(accnos_file)

    # Search for chimeras
    # need to use proper path for vsearch
    m.chimera.vsearch(fasta='current', name='current', group=new_group, dereplicate='t', vsearch=config.get('chimera_params', 'vsearch'))
    # m.chimera.vsearch(fasta='current', count='current', dereplicate='t', vsearch=r'~/HMAS-QC-Pipeline/mothur/vsearch')
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

    # this works, but just takes too long,  over 3 days for Juno data M347-21-026 M347-21-027..
    m.count.seqs(name='current', group='current')
    print(datetime.now()) 
    print(f"done with count.seqs")
    # m.pre.cluster(fasta='current', count='current', diffs=0)
    # print(datetime.now()) 
    # print(f"done with pre.cluster")



    m.summary.seqs(fasta='current', count='current')

    m.cluster(count='current', method='unique', cutoff='unique')
    m.remove.rare(list='current', count='current',
                      nseqs=config.getint('rare_seqs_param', 'nseqs', fallback=9),
                      label='unique')

    current_count = group.get_current_file(MOTHUR_LOG_FILE, 'count').replace(os.path.expanduser('~'), '~', 1)
    m.list.seqs(count=current_count)
    m.get.seqs(fasta='current', accnos='current', name='current', group='current')
    m.summary.seqs(fasta='current', count='current')

    m.rename.file(fasta='current', count='current', prefix=config.get('rename_param', 'prefix') +'.final')

    # convert to full format count_table
    m.count.seqs(count='current', compress='f')


if __name__ == "__main__":
    print("This module is called by pipeline.py.  Please run pipeline.py --help for more information")
