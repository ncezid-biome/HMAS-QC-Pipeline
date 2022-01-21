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

	    		commands.append(['grep', '-cE', f'^{f_primer}[ACGT]{{100,}}{r_primer}$', args.fasta])

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
def check_seq_in_fasta(fasta1, fasta2):

	# the seq ID in fasta1 >M00347_15_000000000-JNGG2_1_2110_18247_10711|ft(frt)(frt)
    grep_command1 = f"cat {fasta1} | grep  '>' | cut -f1 | cut -f1 -d '|'"
    fasta_IDs = exec_grep_shell(grep_command1).split()

    counter = 0
    mis_match = 0
    for id in fasta_IDs:
    	grep_command2 = f"grep -c '{id}' {fasta2}"
    	if int(exec_grep_shell(grep_command2)) > 0:
    		counter += 1
    	else:
    		mis_match += 1

    print (counter, mis_match)


if __name__ == "__main__":

    grep_commands = get_grep_commands_oligo()
    with concurrent.futures.ProcessPoolExecutor(max_workers=60) as executor:
    	results = executor.map(exec_grep, grep_commands)
    	total = 0
    	for result in results:
    		total += int(result)

    print (f"********** total match is: {total}")



    fasta1 = r'/scratch/qtl7/M347_21_026/output_pcr_2013K_0463/juno.good.unique.scrap.pcr.fasta'
    fasta2 = r'/scratch/qtl7/M347_21_026/output_2013K_0463/juno.final.fasta'
    check_seq_in_fasta(fasta1, fasta2)


