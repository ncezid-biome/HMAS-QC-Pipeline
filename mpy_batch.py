#!/usr/bin/env python

from mothur_py import Mothur
m = Mothur()


# Import the config settings


# Import the make.file


def main():
    m.make.contigs(file = 'current', processors=numproc, format=form, oligos=oligosfile, bdiffs=bdif, pdiffs=pdif, checkorient=chkorient, insert=ins, trimoverlap=trimover)
    m.summary.seqs()

    m.rename.file(fasta='current', group='current', prefix=prefx)
    m.summary.seqs()

    m.screen.seqs(fasta='current', group='current', maxambig=maxamb, maxlength=maxlg)
    m.summary.seqs()

    m.unique.seqs(fasta='current')
    m.count.seqs(name='current', group='current')
    m.summary.seqs()

    m.pcr.seqs(fasta='current', oligos='current', pdiffs=pdif2, rdiffs=rdif, count='current')
    m.summary.seqs()

    m.unique.seqs(fasta='current')
    m.summary.seqs()

    m.cluster(count='current', method='unique', cutoff='unique')

    # Search for chimeras
    # Remove chimeras

    m.summary.seqs()
    # need to update fasta somehow

    m.list.seqs(count='current')
    m.get.seqs(fasta='current', accnos='current')
    m.make.shared(list='current', count='current')
    m.get.otulist(list='sal.good.unique.pcr.unique.0.pick.list', label=0)