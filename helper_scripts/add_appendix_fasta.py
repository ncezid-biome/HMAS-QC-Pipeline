#!/usr/bin/env python
import glob, os, os.path
import subprocess

'''
this helper script will take the fasta file name, i.e, 2013K-1067.fasta
append the 2013K-1067 to all its sequence name
And cat them into a new fasta file
'''

exclusion = 'juno'
new_fasta_file = 'juno_design_primers_full_full.fasta'

for fasta_file in glob.glob(rf"*.fasta"):
	base_name = os.path.basename(fasta_file)
	if exclusion not in base_name:
		# appendix = base_name[:-6][-8:] # get rid of .fasta part first, and then keep only the last 8 letters
		appendix = base_name[:-6]
		# awk '{ if (NR %2 == 1) {print $0"-2013K-1067"} else{print $0} }' 2013K-1067.fasta > 2013K-1067.fasta.copy
		command = ("awk '{ if (NR %2 == 1) {print $0"
					f'"-{appendix}"'
					"} else{print $0} }'"
					f' {base_name} > {base_name}.copy')
		print (command)
		subprocess.run(command, shell=True)


subprocess.run(f"cat *.copy > {new_fasta_file}", shell=True)
subprocess.run(f"rm *.copy", shell=True)

