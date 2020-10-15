#!/usr/bin/env python
import os
from mothur_py import Mothur

def main(config):

    currentDir = os.getcwd()
    m = Mothur()
    m.set.dir(input = config.get('file_inputs', 'input_dir', fallback = currentDir))
    m.set.dir(output = config.get('file_inputs', 'output_dir', fallback = currentDir)) 
    m.set.logfile(name = config.get('file_inputs', 'output_dir', fallback = currentDir))
    m.make.contigs(file = config.get('file_inputs','batch_file'),
                    processors= config.getint('contigs_params', 'processors', fallback = 40), 
                    format=config.get('contigs_params', 'format', fallback = 'illumina1.8+'), 
                    oligos=config.get('file_inputs', 'oligos'), 
                    bdiffs=config.getint('contigs_params', 'bdiffs', fallback = 0), 
                    pdiffs=config.getint('contigs_params', 'pdiffs', fallback = 0), 
                    checkorient=config.get('contigs_params', 'checkorient', fallback = 't'), 
                    insert=config.getint('contigs_params', 'insert', fallback = 25), 
                    trimoverlap=config.get('contigs_params', 'trimoverlap', fallback = 'f'))
    m.summary.seqs()

    m.rename.file(fasta='current', group='current', 
                    prefix=config.get('rename_param', 'prefix'))
    m.summary.seqs()

    m.screen.seqs(fasta='current', group='current',
                    maxambig=config.getint('screen_params', 'maxambig', fallback = 0),
                    maxlength=config.getint('screen_params', 'maxlength', fallback = 325))
    m.summary.seqs()

    m.unique.seqs(fasta='current')
    m.count.seqs(name='current', group='current')
    m.summary.seqs()

    m.pcr.seqs(fasta='current', 
                oligos='current', 
                pdiffs=config.getint('pcr_params', 'pdiffs', fallback = 0), 
                rdiffs=config.getint('pcr_params', 'rdiffs', fallback = 0), 
                count='current')
    m.summary.seqs()

    m.unique.seqs(fasta='current')
    m.summary.seqs()

    m.cluster(count='current', method='unique', cutoff='unique')

    m.remove.rare(list='current', count='current', 
                    nseqs=config.getint('rare_seqs_param', 'nseqs', fallback = 9), 
                    label='unique')
    m.summary.seqs()

    # Search for chimeras
    # Remove chimeras

    m.summary.seqs()

    m.list.seqs(count='current')
    m.get.seqs(fasta='current', accnos='current')
    m.make.shared(list='current', count='current')
    m.get.otulist(list='sal.good.unique.pcr.unique.0.pick.list', label=0)
