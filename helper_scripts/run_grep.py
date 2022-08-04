import subprocess
import argparse
import concurrent.futures
import sys, os
import glob
import pandas as pd
import logging
import settings
import utilities

sys.path.insert(0,r'..') # to allow import packages in parent folder
import group

logger = logging.getLogger('root')
logger.setLevel(logging.INFO)
logger.addHandler(settings.log_handler)

sample_list = ['AR_0409','2016K_0167','2014K_0527','2014K_0421','2013K_1828','2015K_0074',
				'2013K_1153','2013K_0463','Water1','2011K_0222','Blank_IRP1','2013K_1067',
				'Blank_IRP2','2011K_0052','2010K_1649','08_0810','2010K_2370']


# this method prepares the list of grep commands from oligos file
def get_grep_commands_oligo():

	parser = argparse.ArgumentParser(description = 'grep the primer sequence in fasta file')
	parser.add_argument('-o', '--oligos', metavar = '', required = True, help = 'Specify oligos file')
	parser.add_argument('-f', '--fasta', metavar = '', required = True, help = 'Specify fasta file')
	args = parser.parse_args()

	commands = []
	with open(args.oligos, 'r') as f:
		for line in f:
			if line.startswith('primer'):
				f_primer = line.strip().split()[1]
				r_primer = utilities.revcomp(line.strip().split()[2])

				commands.append(['grep', '-cE', f'{f_primer}[ACGT]{{50,}}{r_primer}', args.fasta])

	print (len(commands))
	return commands


def exec_grep(command):

	process = subprocess.run(command, capture_output=True, text=True)
	if int(process.stdout) > 0:
		print (f"{command[2]} has {process.stdout.strip()} matches")
	else:
		print (f"{command[2]} no matches")
	return int(process.stdout)

# this method uses the shell=True option because 
# it need to use pipe in shell
def exec_grep_shell(grep_command):

	# seems like you have to set shell=True in order to use pipe in shell 
	process = subprocess.run(grep_command, capture_output=True, text=True, shell=True)
	return process.stdout.strip()

# this method returns a list of all primers for the given sample
# for the given design.fasta file
def get_primers_for_sample_design_fasta(design_fasta, sample):

	grep_command = f"grep '{sample}' {design_fasta} | cut -d '-' -f2"
	return exec_grep_shell(grep_command).split()

# this method checks the seqs in fasta1, to see if they also 
# appear in fasta2
def check_seq_in_3fasta(fasta1, fasta2, fasta3):

	# the seq ID in fasta1 >M00347_15_000000000-JNGG2_1_2110_18247_10711|ft(frt)(frt)
	grep_command1 = f"cat {fasta1} | grep  '>' | cut -f1 | cut -f1 -d '|'"
	fasta_IDs = exec_grep_shell(grep_command1).split()

	# this grep the actual sequences, fasta_IDs here are sequences instead
	grep_command1 = f"cat {fasta1} | grep -v '>'"
	fasta_IDs = exec_grep_shell(grep_command1).split()

	counter = 0
	mis_match = 0
	mis_match_match = 0
	for id in fasta_IDs:
		grep_command2 = f"grep -c '{id}' {fasta2}"
		if int(exec_grep_shell(grep_command2)) > 0:
			counter += 1 
		else:
			grep_command3 = f"grep -c '{id}' {fasta3}"
			mis_match += 1
			if int(exec_grep_shell(grep_command3)) > 0:
				mis_match_match += 1

	print (counter, mis_match, mis_match_match)


def check_seq_in_2fasta(fasta1, fasta2):

	# the seq ID in fasta1 >M00347_15_000000000-JNGG2_1_2110_18247_10711|ft(frt)(frt)
	# grep_command1 = f"cat {fasta1} | grep  '>' | cut -f1 | cut -f1 -d '|'"
	# fasta_IDs = exec_grep_shell(grep_command1).split()

	# this grep the actual sequences, fasta_IDs here are sequences instead
	grep_command1 = f"cat {fasta1} | grep -v '>'"
	fasta_IDs = exec_grep_shell(grep_command1).split()

	counter = 0
	for id in fasta_IDs:
		grep_command2 = f"grep -c '{id}' {fasta2}"
		if int(exec_grep_shell(grep_command2)) > 1:
			counter += 1 
			grep_command3 = f"grep -A1 '{id}' {fasta1}"
			print(exec_grep_shell(grep_command3))

	print (counter)

def check_dup_seq(fasta):

	# this grep the actual sequences, fasta_IDs here are sequences
	grep_command1 = f"cat {fasta} | grep -v '>'"
	fasta_IDs = exec_grep_shell(grep_command1).split()

	dup_list = []
	for id in fasta_IDs:
		grep_command2 = f"grep -c '{id}' {fasta}"
		duplicate = int(exec_grep_shell(grep_command2))-1
		if duplicate >= 1:
			dup_list.append(duplicate)

	return dup_list


def check_dup_seq_dict(fasta):

	seq_dict = {}
	with open(fasta, 'r') as f:
		for ind, row in enumerate(f.readlines(), start=1):
			if ind%2 == 1:
				last_seqID = '-'.join(row.strip().split('-')[1:])
			else:
				last_key = row.strip()
				if last_key in seq_dict:
					seq_dict[last_key].append(last_seqID)
				else:
					seq_dict[last_key] = [last_seqID]

	print (list(seq_dict.items())[:3])

	uniq_seq = []
	dup_seq = [] #duplicate sequences within sample ()
	for seq in seq_dict:
		if len(seq_dict[seq]) <= 1:
			uniq_seq.append(seq_dict[seq][0])
		else:
			seqID_set = set([i.split('-')[0] for i in seq_dict[seq]])
			print (seqID_set)
			if len(seqID_set) > 1:
				dup_seq.append(seq_dict[seq])

	print (f"uniq_seq size: {len(uniq_seq)}")
	print (uniq_seq[:3])
	print (f"seq_dict size: {len(seq_dict)}")
	print (f"dup_seq size: {len(dup_seq)}")
	print (dup_seq)


def check_seq_in_2fasta_set(fasta1, fasta2):

	set_fasta1 = set()
	set_fasta2 = set()

	with open(fasta1, 'r') as f:

	    for ind, row in enumerate(f.readlines(), start=1):
	        if ind%2 == 1:
	            seq_ID = row.strip().split()[0]
	            # print (seq_ID)
	            set_fasta1.add(seq_ID)
	        else:
	            continue

	with open(fasta2, 'r') as f:

	    for ind, row in enumerate(f.readlines(), start=1):
	        if ind%2 == 1:
	            seq_ID = row.strip().split()[0]
	            # print (seq_ID)
	            set_fasta2.add(seq_ID)
	        else:
	            continue

	# print (len(set_fasta1-set_fasta2))
	# print (set_fasta1-set_fasta2)

	for id in list(set_fasta1-set_fasta2)[:1000]:
		grep_command = f"grep -c '{id}' juno.good.fasta"
		if int(exec_grep_shell(grep_command)) > 0:
			# grep_command3 = f"grep -A1 '{id}' {fasta1}"
			print(exec_grep_shell(grep_command))


# helper method to filter Sean's design primer fasta file
# to include only those 2461 primers that're used in our oligo file
def filter_design_fasta(fasta, new_fasta):

	new_fasta_list = []
	primers = utilities.Primers(settings.OLIGO_FILE)

	with open(fasta, 'r') as f:

		flag = False
		for ind, row in enumerate(f.readlines(), start=1):

			if ind%2 == 1:

				last_key = row.strip().split('-')[1]
				if last_key in primers.pseqs.keys():
					new_fasta_list.append(row.strip())
					flag = True
			else:
				if flag:
					new_fasta_list.append(row.strip())
					flag = False

	with open(new_fasta, 'w') as f:
		f.write('\n'.join(new_fasta_list))

def get_subsampling_commands(ratios):

	commands = []
	for ratio in ratios:
		commands.append(f"reformat.sh in1=Undetermined_S0_L001_I1_001.fastq \
									  in2=Undetermined_S0_L001_I2_001.fastq \
									  out1=out{ratio}.I1.fastq \
									  out2=out{ratio}.I2.fastq \
									  samplerate={ratio} sampleseed=77")

		commands.append(f"reformat.sh in1=out10.R1.fastq \
									  in2=out10.R2.fastq \
									  out1=out{ratio}.R1.fastq \
									  out2=out{ratio}.R2.fastq \
									  samplerate={ratio} sampleseed=77")

	return commands


def get_primer_prediction_dict_by_fasta(fasta_file):

	primers = utilities.Primers(settings.OLIGO_FILE)
	sample_name = os.path.basename(fasta_file)[:-6]

	# primer_list = []
	# for idx, key in enumerate(primers.pseqs, start=1):
	# 	grep_command = f"grep -c {key} {fasta_file}"
	# 	if int(exec_grep_shell(grep_command)) > 0:
	# 		primer_list.append(key)
	# print (sample_name, len(primer_list))

	design_primer = set()
	with open(fasta_file, 'r') as f:

	    for row in f.readlines():
	        if row.startswith('>'):
	            primer = row.strip().split('-')[1]
	            design_primer.add(primer)

	logger.info (f"sample: {sample_name} has {len(design_primer & set(primers.pseqs.keys()))} primers")
	# print (sample_name, set(primers.pseqs.keys()) - design_primer)
	return (sample_name,(design_primer & set(primers.pseqs.keys())))


def get_primer_prediction_dict_concurrent():

	design_folder = r'/scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/CIMS/Salmonella/'\
					'T3Pio/Juno_Pilot/design_primers_copy'
	fasta_files = [fasta for fasta in glob.glob(rf"{design_folder}/*.fasta")]
	prediction_dict = {}
	with concurrent.futures.ProcessPoolExecutor(max_workers=30) as executor:
		results = executor.map(get_primer_prediction_dict_by_fasta, fasta_files)
		for result in results:
			prediction_dict[result[0]] = result[1]

	return prediction_dict


# this helper method grabs all 2461 primers in oligo file
# and concatenate the consambig consensus seqs for each amplicon into one synthetic 'reference genome'
# note you need to have muscle / consambig in path
def create_design_fasta_consambig(design_fasta, design_dir):

	primers = utilities.Primers(settings.OLIGO_FILE)

	for primer in primers.pseqs:
		#1. grep
		grep_command = f"grep -A1 '{primer}' {design_fasta} > {design_dir}/{primer}.grep"
		exec_grep_shell(grep_command)

		#2. muscle
		muscle_command = f"muscle -align {design_dir}/{primer}.grep -output {design_dir}/{primer}.msa"
		exec_grep_shell(muscle_command)

		#3. consambig
		consambig_command = f"consambig -sequence {design_dir}/{primer}.msa -outseq {design_dir}/{primer}.fa -name {primer}"
		exec_grep_shell(consambig_command)

		#4. clean up temp files
		rm_command = f"rm {design_dir}/{primer}.grep {design_dir}/{primer}.msa"
		exec_grep_shell(rm_command)

	cat_command = f"cat {design_dir}/*.fa > {design_dir}/design.consambig.consensus.fasta && rm {design_dir}/*.fa"
	exec_grep_shell(cat_command)



def make_deunique_dict(fasta):

	# this grep the actual sequences, fasta_IDs here are sequences
	grep_command1 = f"cat {fasta} | grep -v '>'"
	fasta_IDs = exec_grep_shell(grep_command1).split()
	sample_name = os.path.basename(fasta).split('.')[0]

	deunique_dict = {}
	for id in fasta_IDs:
		grep_command2 = f"grep -B1 '{id}' {fasta} | head -n1"
		seq_ID = exec_grep_shell(grep_command2)
		deunique_dict[id] = (seq_ID.partition(f'{sample_name}')[0][:-1], sample_name + seq_ID.partition(f'{sample_name}')[2])

	# print (deunique_dict)
	return deunique_dict

def make_deunique_dict_for_all_samples(fasta):

	# this grep the actual sequences, fasta_IDs here are sequences
	grep_command1 = f"cat {fasta} | grep -v '>'"
	fasta_IDs = exec_grep_shell(grep_command1).split()

	deunique_dict = {}
	for id in fasta_IDs:
		grep_command2 = f"grep -B1 '{id}' {fasta} | head -n1"
		seq_ID = exec_grep_shell(grep_command2)
		for sample_name in sample_list:
			if sample_name in seq_ID:
				deunique_dict[id] = (seq_ID.partition(f'{sample_name}')[0][:-1], sample_name + seq_ID.partition(f'{sample_name}')[2])
				break

	# print (deunique_dict)
	return deunique_dict

def make_fasta_for_deunique(fasta):

	sample_name = os.path.basename(fasta).split('.')[0]
	new_fasta_list = []
	deunique_dict = make_deunique_dict(fasta)
	for seq in deunique_dict:
		new_fasta_list.append(deunique_dict[seq][0])
		new_fasta_list.append(seq)

	with open(f'{os.path.dirname(fasta)}/{sample_name}.deunique.fasta', 'w') as f:
		f.write('\n'.join(new_fasta_list))

def append_sample_primer_for_deunique(dict_fasta, redundant_fasta):

	# deunique_dict = make_deunique_dict(dict_fasta)
	deunique_dict = make_deunique_dict_for_all_samples(dict_fasta)
	new_fasta_list = []

	with open(redundant_fasta, 'r') as f:
		row_iter = iter(f.readlines())
		for row in row_iter:
			if row.startswith('>'):
				seq_ID = row.strip()
				seq = next(row_iter).strip()
				new_fasta_list.append(f"{seq_ID}_{deunique_dict[seq][1]}")
				new_fasta_list.append(seq)

	with open(f'{redundant_fasta[:-5]}updated.fasta', 'w') as f:
		f.write('\n'.join(new_fasta_list))

def make_redundant_fasta_persample_from_top100fasta(top100fasta, count_table):

	## 1. grep -A1 2013K_0463 random100and3controls_update.fasta > 2013K_0463.fasta
	##    for each sample
	## 2. make_fasta_for_deunique(2013K_0463.fasta)  for each sample.fasta
	## 3. call Mothur's deunique.seqs(fasta=2013K_0463.deunique.fasta,count=juno.final.count_table)
	## 4. append_sample_primer_for_deunique(2013K_0463.fasta, 013K_0463.deunique.redundant.fasta)
	## 5. rm 2013K_0463.deunique.fasta  {work_dir}/temp.{sample}.batch

	work_dir = os.path.dirname(top100fasta)
	for sample in sample_list:

		grep_command = f"grep -A1 {sample} {top100fasta} > {work_dir}/{sample}.fasta"
		exec_grep_shell(grep_command)

		if os.path.getsize(f"{work_dir}/{sample}.fasta") > 0:
			make_fasta_for_deunique(f"{work_dir}/{sample}.fasta")

			# create batch file for mothur
			echo_command = (f"echo 'deunique.seqs(fasta={work_dir}/{sample}.deunique.fasta,count={count_table})'"
							f" > {work_dir}/temp.{sample}.batch")
			mothur_command = f"mothur {work_dir}/temp.{sample}.batch"
			exec_grep_shell(echo_command)
			exec_grep_shell(mothur_command)

			append_sample_primer_for_deunique(f"{work_dir}/{sample}.fasta", f"{work_dir}/{sample}.deunique.redundant.fasta")

			rm_command = f"rm -f {work_dir}/{sample}.deunique.fasta {work_dir}/temp.{sample}.batch"
			exec_grep_shell(rm_command)




if __name__ == "__main__":

	# grep_commands = get_grep_commands_oligo()
	# with concurrent.futures.ProcessPoolExecutor(max_workers=60) as executor:
	# 	results = executor.map(exec_grep, grep_commands)
	# 	total = 0
	# 	for result in results:
	# 		total += int(result)

	# print (f"********** total match is: {total}")

	# ratios = [0.4,0.6,0.7,0.9]
	# with concurrent.futures.ProcessPoolExecutor(max_workers=10) as executor:
	# 	executor.map(exec_grep_shell, get_subsampling_commands(ratios))


	# # fasta2 = r'/scratch/qtl7/M347_21_026/output_2013K_0463/juno.final.fasta'
	# # fasta1 = r'/scratch/qtl7/M347_21_026/output_pcr_2013K_0463/juno.final.fasta'
	# fasta3 = r'/scratch/qtl7/M347_21_026/output_pcr_2013K_0463/juno.good.unique.scrap.pcr.fasta'
	# check_seq_in_3fasta(fasta1, fasta2, fasta3)
	# fasta1 = r'/scratch/qtl7/M347_21_026/output1_cutadapt/Undetermined.paired.scrap.contigs.fasta'
	# fasta2 = r'/scratch/qtl7/M347_21_026/output3_trim_noprecluster/Undetermined_trimed.paired.scrap.contigs.fasta'
	# # check_seq_in_2fasta(fasta1, fasta2)
	# check_seq_in_2fasta_set(fasta1, fasta2)
	# check_seq_in_2fasta_set(fasta2, fasta1)

	# oligo_file = r'/scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/CIMS/Salmonella/'\
	# 				'cgMLST_pilot/M347-21-026_fastq/M347-21-026.oligos'
	# design_folder = r'/scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/CIMS/Salmonella/'\
	# 				'T3Pio/Juno_Pilot/design_primers_copy'
	# primers = group.Primers(oligo_file)
	# mis_match = []
	# total = []
	# for fasta in glob.glob(rf"{design_folder}/*.fasta"):
	# 	for idx, key in enumerate(primers.pseqs, start=1):
	# 		grep_command = f"grep -c {key} {fasta}"
	# 		if int(exec_grep_shell(grep_command)) == 0:
	# 			mis_match.append(key)
	# 	print (mis_match)
	# 	print (os.path.basename(fasta))
	# 	print (len(mis_match))
	# 	mis_match.append('\n')
	# 	mis_match.append(os.path.basename(fasta))
	# 	mis_match.append(len(mis_match))
	# 	total.append(mis_match)
	# 	mis_match = []



	# with open('temp_primer', 'w') as f:
	# 	f.write('\n'.join(total))

	# new_df = pd.read_pickle(f'perfect_match_sum_dilution.pkl')
	# print (new_df.head(n=5))
	# column_list = list(set([i.split('.')[1] for i in new_df.index]))
	# column_list.sort(key=str.lower)
	# print (len(set(column_list)))
	# print (len(set(primers.pseqs)))
	# mis_match = set(primers.pseqs) - set(column_list)
	# print (len(mis_match))


	# print ("in run_grep, empty main()")
	# get_primer_prediction_dict_concurrent()

	# design_fasta = r'/scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/CIMS/Salmonella/'\
	# 				'T3Pio/Juno_Pilot/design_primers_copy/juno_design_primers_full_full.fasta'
	# new_fasta = r'/scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/CIMS/Salmonella/'\
	# 				'T3Pio/Juno_Pilot/design_primers_copy/juno_design_primers_full_full_filtered.fasta'

	# filter_design_fasta(design_fasta, new_fasta)

	fasta = r'/scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/CIMS/Salmonella/'\
					'T3Pio/Juno_Pilot/design_primers_copy/juno_design_primers_full_full_filtered.fasta'
	# # dup_list = check_dup_seq(fasta)
	# check_dup_seq_dict(fasta)

	# design_fasta = r'/scicomp/home-pure/qtl7/t3pio/Code_Repository/T3Pio_Main/juno_design.fasta'
	# design_dir = r'/scicomp/home-pure/qtl7/t3pio/Code_Repository/T3Pio_Main'
	# create_design_fasta_consambig(design_fasta,design_dir)

	'''
	fasta=r'/scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/CIMS/Salmonella/cgMLST_pilot/M347-21-026_fastq/qcrun_result/top100and3controlsupdate.fasta'
	redundant_fasta = (f'/scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/CIMS/Salmonella/cgMLST_pilot/M347-21-026_fastq'
						f'/qcrun_result/top100and3controlsupdate.deunique.fasta')
	# make_deunique_dict(fasta)
	# make_fasta_for_deunique(fasta)
	# append_sample_primer_for_deunique(fasta, redundant_fasta)

	top100fasta = (f'/scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/CIMS/Salmonella/cgMLST_pilot/M347-21-026_fastq'
						f'/qcrun_result/M3235_22_024_run/random100and3controls_update.fasta')
	count_table = (f'/scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/CIMS/Salmonella/cgMLST_pilot/M347-21-026_fastq'
						f'/qcrun_result/M3235_22_024_run/juno.final.count_table')
	make_redundant_fasta_persample_from_top100fasta(top100fasta, count_table)

	'''
