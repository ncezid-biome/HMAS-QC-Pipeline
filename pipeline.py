
import logging, sys, os, argparse

#### Import config file with user inputs
# Check that all the inputs are there

# Add: a switch that evaluates config settings based on user choice
# Add getops flags: batch mode vs. single sample

aparser = argparse.ArgumentParser(description='Run in batch mode or run a single sample.')

aparser.add_argument('-b', action='store', help='Run in batch mode')
aparser.add_argument('-s', action='store', help='Run a single sample')

args = aparser.parse_args()

from configparser import ConfigParser

parser = ConfigParser()
parser.read('settings.ini')

# Function to check the config file
def ConfigCheck(section):
    dict1 = {}
    inputkeys = parser.options(section)
    for option in inputkeys:
        try:
            dict1[option] = parser.get(section, option)
            if dict1[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1

# Import and check all the variables

# Note: 
input_dir = ConfigCheck('file_inputs')['input_dir']
file_prefix = ConfigCheck('file_inputs')['file_prefix']
read1 = ConfigCheck('file_inputs')['read1']
read2 = ConfigCheck('file_inputs')['read2']
index1 = ConfigCheck('file_inputs')['index1']
index2 = ConfigCheck('file_inputs')['index2']
oligosfile = ConfigCheck('file_inputs')['oligos']
numproc = int(ConfigCheck('contigs_params')['processors'])
form = ConfigCheck('contigs_params')['format']
bdif = int(ConfigCheck('contigs_params')['bdiffs'])
pdif = int(ConfigCheck('contigs_params')['pdiffs'])
chkorient = ConfigCheck('contigs_params')['checkorient']
ins =  int(ConfigCheck('contigs_params')['insert'])
trimover = ConfigCheck('contigs_params')['trimoverlap']
prefx = ConfigCheck('rename_param')['prefix']
maxamb = int(ConfigCheck('screen_params')['maxambig'])
maxlg = int(ConfigCheck('screen_params')['maxlength'])
pdif2 = int(ConfigCheck('pcr_params')['pdiffs'])
rdif = int(ConfigCheck('pcr_params')['rdiffs'])
numseqs = int(ConfigCheck('rare_seqs_param')['nseqs'])

# Check that the input files exist
for file in (read1, read2, index1, index2, oligosfile):
    if not os.path.exists(file):
        missing_file_flag = 1
        print("File %s not found" % file)
        
if missing_file_flag == 1:
    print("Program exited because an input file was not found", sys.stderr)
    sys.exit(1)

# Import Mothur
try:
    from mothur_py import Mothur
    m = Mothur()
except:
    print("Unable to import mothur_py module.  Is it installed and on the path?")
    print("Program exited because mothur_py could not be imported", sys.stderr)
    sys.exit(1)
    
#### Set directories for the raw files (mothur_py should dump intermediate files in the current folder?)
    # Decide how to work out the directories


    

#### Run the mothur_py steps


m.make.file(inputdir= input_dir, type=gz, prefix=file_prefix)
# Need a separate set of commands that runs the batch mode (file=current)
    
m.make.contigs(ffastq = read1, rfastq = read2, findex = index1, rindex = index2, processors=numproc, format=form, oligos=oligosfile, bdiffs=bdif, pdiffs=pdif, checkorient=chkorient, insert=ins, trimoverlap=trimover)
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

#m.set.current(fasta='sal.good.unique.pcr.fasta', count='sal.good.pcr.count_table')


m.unique.seqs(fasta='current')
m.summary.seqs()

m.cluster(count='current', method='unique', cutoff='unique')

# Search for chimeras
# Remove chimeras



m.remove.rare(list='current', count='current', nseqs=numseqs, label='unique')
m.summary.seqs()
# need to update fasta somehow

m.list.seqs(count='current')
m.get.seqs(fasta='current', accnos='current')
m.make.shared(list='current', count='current')
m.get.otulist(list='sal.good.unique.pcr.unique.0.pick.list', label=0)
