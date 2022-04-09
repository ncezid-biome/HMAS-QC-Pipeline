#!/usr/bin/env python
import sys, os, shutil, subprocess, argparse
import group
import glob
import re
from collections import deque
from datetime import datetime
import concurrent.futures


def cmd_exists(command):
	"""Checks for existence of a command on user's path

	Params
	------
	command: String
		Name of the command

	Returns
	------
	Path to executable command
	
	"""
	if shutil.which(command) is None:
		print(f"{command} not found on PATH")
		sys.exit()
	else:
		return(shutil.which(command))

def run_cutadapt(command):
	"""
	run cutadapt command

	Params
	------
	command: list (the actual run command in form of a list)
	
	"""

	p = subprocess.run(command)
	if p.returncode != 0:
		print("cutadapt could not be executed.  Error encountered.")
		print(p.returncode)
	else:
		print (f"{datetime.now()}   done with {command[len(command)-1]}")

def run_cutadapt_mothur(config, new_fasta, num_process):
	"""
	run cutadapt to remove primers in the middle of a Mothur run.
	It is customized to run after the split.groups command and use cutadapt to remove primers in 
	a multi-process fashion. After the cutadapt part is done, it will also merge all fasta files 
	into a 'new_fasta' fasta file
	*** split.groups will split the original fasta file based on primer groups
	*** cutadapt will run on each of those smaller fasta file

	Params
	------
	config: config object
		Name of the command

	new_fasta: String
		name of the new fasta file (merged)

	num-process: int
		designated number of process to run the script

	"""

	#1. prep for cutadapt commands
	cutadapt_commands = deque()
	primers = group.Primers(os.path.expanduser(config['file_inputs']['oligos']))
	cutadapt_cmd = cmd_exists('cutadapt')
	for key in primers.pseqs:
		for current_file in glob.glob(rf"{os.path.expanduser(config['file_inputs']['output_dir'])}/*{key}.fasta"):
			# cutadapt_commands.append([cutadapt_cmd, f'--info-file={primers.pseqs[key][0]}_info.tsv', '-g', \
			cutadapt_commands.append([cutadapt_cmd, '-g', \
				f'{primers.pseqs[key][0]}...{primers.pseqs[key][1]}', '-o', \
				f'{current_file}_cutadapt', \
				f'{current_file}'])
			break
	print(datetime.now())
	print(len(cutadapt_commands))


	#2. run cutadapt 
	with concurrent.futures.ProcessPoolExecutor(max_workers=num_process) as executor:
		executor.map(run_cutadapt, cutadapt_commands)

	print(datetime.now()) 
	print(f"done with cutadapt")

	#3. merge all fastas 
	with open(new_fasta,'wb') as wfd:
		for f in glob.glob(rf"{os.path.expanduser(config['file_inputs']['output_dir'])}/*_cutadapt"):
			with open(f,'rb') as fd:
				shutil.copyfileobj(fd, wfd)

	print(datetime.now()) 
	print(f"done with merging files")