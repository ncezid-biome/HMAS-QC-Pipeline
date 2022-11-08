import subprocess
import argparse
import concurrent.futures
import pandas as pd

TIMEOUT = 600 # time out at 10 minutes

def parse_argument():
    # note
    # the input file is a plain 2 columns file (tab delimited)
    # 1st column is sample ID, 2nd column is SRA number
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', metavar = '', required = True, help = 'Specify SRA list file')
    return parser.parse_args()

def download_SRA_assemble(sample_sra):
    
    sample = sample_sra[0]
    sra = sample_sra[1]
    
    # the procedure is:
    # 1. create a directory(sample name) and change to that directory
    # 2. delete any previous downloads and download the SRA file (split into R1, R2)
    # 3. load shovill
    # 4. call shovill to assemble 
    # use skesa and trim usual adapters, set to quite mode
    command = (f"mkdir -p {sample} && cd {sample} && " 
                f"rm -rf {sra}* && "
                f"fasterq-dump --split-files {sra} && " 
                f"ml shovill && "
                f"shovill -R1 {sra}_1.fastq -R2 {sra}_2.fastq --outdir . --force "
                f"--assembler skesa --trim ON --cpus 10 > /dev/null 2>&1  && "
                f"mv contigs.fa {sample}_assembled.fa")
                
    # shovill sometimes will freeze for no obvious reasons, set timeout to 10 minutes before re-run
    process = subprocess.run(command, capture_output=True, text=True, shell=True, timeout=TIMEOUT)
    if process.returncode == 0:    
        print (f"{sample} WGS assembly is completed")
    else:
        print (f"error out while in {sample} WGS assembling")
        print (process)
    
    
if __name__ == "__main__":
    
    args = parse_argument()
    
    df = pd.read_csv(args.input_file, names=['sample','sra'], delim_whitespace=True)
    df.dropna(inplace=True) # remove any empty/invalid entry
    sample_sras = zip(df['sample'], df['sra'])
        
    with concurrent.futures.ProcessPoolExecutor(max_workers=20) as executor:
            results = executor.map(download_SRA_assemble, sample_sras)