import subprocess
import argparse
import concurrent.futures

def revcomp(myseq):
	rc = {'A' : 'T', 'T' : 'A', 'G' : 'C', 'C' : 'G'}
	seq = [rc[n] if n in rc else n for n in myseq]
	return("".join(list(reversed(seq))))


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
				r_primer = revcomp(line.strip().split()[2])

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
	grep_command1 = f"cat {fasta1} | grep  '>' | cut -f1 | cut -f1 -d '|'"
	fasta_IDs = exec_grep_shell(grep_command1).split()

	# this grep the actual sequences, fasta_IDs here are sequences instead
	# grep_command1 = f"cat {fasta1} | grep -v '>'"
	# fasta_IDs = exec_grep_shell(grep_command1).split()

	counter = 0
	for id in fasta_IDs:
		grep_command2 = f"grep -c '{id}' {fasta2}"
		if int(exec_grep_shell(grep_command2)) <= 0:
			counter += 1 
			# grep_command3 = f"grep -A1 '{id}' {fasta1}"
			# print(exec_grep_shell(grep_command3))

	print (counter)

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


if __name__ == "__main__":

	grep_commands = get_grep_commands_oligo()
	with concurrent.futures.ProcessPoolExecutor(max_workers=60) as executor:
		results = executor.map(exec_grep, grep_commands)
		total = 0
		for result in results:
			total += int(result)

	print (f"********** total match is: {total}")


	# # fasta2 = r'/scratch/qtl7/M347_21_026/output_2013K_0463/juno.final.fasta'
	# # fasta1 = r'/scratch/qtl7/M347_21_026/output_pcr_2013K_0463/juno.final.fasta'
	# fasta3 = r'/scratch/qtl7/M347_21_026/output_pcr_2013K_0463/juno.good.unique.scrap.pcr.fasta'
	# check_seq_in_3fasta(fasta1, fasta2, fasta3)
	# fasta1 = r'/scratch/qtl7/M347_21_026/output1_cutadapt/Undetermined.paired.scrap.contigs.fasta'
	# fasta2 = r'/scratch/qtl7/M347_21_026/output3_trim_noprecluster/Undetermined_trimed.paired.scrap.contigs.fasta'
	# # check_seq_in_2fasta(fasta1, fasta2)
	# check_seq_in_2fasta_set(fasta1, fasta2)
	# check_seq_in_2fasta_set(fasta2, fasta1)


