import subprocess
import argparse
import concurrent.futures
import pandas as pd
from pathlib import Path
import glob
import os

'''
This script will do the assembly with shovill.
It will write out any files that failed during the process, into a file with extension _fail_to_assemble,
and you can re-run the script with that _fail_to_assemble file

'''
TIMEOUT = 2000 # time out at 2000 seconds

def parse_argument():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', metavar = '', required = True, help = 'Specify input reads folder')
    parser.add_argument('-o', '--output', metavar = '', required = True, help = 'Specify output folder')
    return parser.parse_args()

def shovill_assemble(reads_info):
    
    r1 = reads_info[0]
    r2 = reads_info[1]
    name = reads_info[2]
    out_dir = reads_info[3]
    outbreak = reads_info[4]
    
    # the procedure is:
    # 1. create a directory(sample name) and change to that directory
    # 2. load shovill
    # 3. call shovill to assemble 
    # use skesa and trim usual adapters, set to quite mode
    command = (f"mkdir -p {out_dir}/{name} {out_dir}/{outbreak} && " 
                f"shovill -R1 {r1} -R2 {r2} --outdir {out_dir}/{name} --force "
                # f"--assembler skesa --trim ON --cpus 4 > /dev/null 2>&1  && "
                f"--assembler skesa --trim ON --cpus 4 && "
                f"mv {out_dir}/{name}/contigs.fa {out_dir}/{outbreak}/{name}_assembled.fasta && "
                f"rm -rf {out_dir}/{name}")
                
    # shovill sometimes will freeze for no obvious reasons, set timeout to 20 minutes before re-run
    try:
        process = subprocess.run(command, capture_output=True, text=True, check=True, shell=True, timeout=TIMEOUT)
    except subprocess.TimeoutExpired:
        print (f"timed out after {TIMEOUT} seconds while in {name} WGS assembling")
        return (reads_info)
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.stderr)
    
    if process.returncode == 0:    
        print (f"{name} WGS assembly is completed")
    else:
        print (f"error out while in {name} WGS assembling")
        return (reads_info)
    
    
if __name__ == "__main__":
    
    args = parse_argument()

    r1_files = glob.glob(f'{args.input}/*_R1.fastq.gz')
    r2_files = glob.glob(f'{args.input}/*_R2.fastq.gz')

    # Sort the lists of file names to ensure that they are paired correctly
    r1_files.sort()
    r2_files.sort()

    # Pair the files and extract the common pattern
    paired_files = []
    for r1, r2 in zip(r1_files, r2_files):
        # Get the base file name without the extension
        r1_base = os.path.splitext(os.path.basename(r1))[0]
        r2_base = os.path.splitext(os.path.basename(r2))[0]

        # Extract the common part of the file names
        common_part = r1_base[:r1_base.rfind('_R')]
        paired_files.append((r1, r2, common_part, args.output, Path(args.input).stem))
    
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=10) as executor:
            results = executor.map(shovill_assemble, paired_files)
            
            # write out the files tha failed in the process
            df_list = []
            for result in results:
                if result:
                    df_list.append(list(result))
                    
            if len(df_list) >= 1:
                df = pd.DataFrame(df_list)
                df.to_csv(f"{Path(args.input).stem}_fail_to_assemble", sep='\t', index=False, header=False)
    
