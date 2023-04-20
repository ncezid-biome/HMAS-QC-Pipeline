## Process raw Miseq reads and run them through Mothur QC pipeline and analysis script  
<br>  

### 1. process raw reads  
### <br>
1.  &nbsp;&nbsp; `module load bcl2fastq/2.20`  
2. &nbsp;&nbsp; `bcl2fastq -p 20 --create-fastq-for-index-reads -o <outputdir>`  
<br>
**note**: run the above bcl2fastq command in the Miseq raw reads output folder. (be sure to change the samplesheet.csv to a different name, so it can lump reads together and give us the 4 ‘Undetermined R1/R2/I1/I2 fastq files)  

### 2. trim nextera adapter  
### <br>

>  `cutadapt --pair-adapters`  
>  `-j 20`  
>  `-m 1`  
>  `-a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC`  
>  `-A CTGTCTCTTATACACATCTGACGCTGCCGACGA`  
>  `-o Undetermined_S0_L001_R1_001.fastq.trim.gz`  
>  `-p Undetermined_S0_L001_R2_001.fastq.trim.gz`  
>  `Undetermined_S0_L001_R1_001.fastq.gz`  
>  `Undetermined_S0_L001_R2_001.fastq.gz`  
<br>
**note**: -m 1 will remove any zero length sequence; -j 20 is for specifying CPU number  

### 3. Mothur QC
### <br>

>  `python3 pipeline.py -c settings.in`  

Need to make sure that:  
1. settings.ini is properly set up, for example:  
>  `/scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/CIMS/HMAS_QC_pipeline/M3235_22_024/settings.ini`  

2.  mothur_py is installe  
3.  mothur is on path (**v1.46.0**)  
4.  cutadapt is on path  
5.  oligos file is properly set up, for example:  
>  `/scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/CIMS/HMAS_QC_pipeline/M3235_22_024/M3235_22_024.oligos`  

### 4. Run analysis script
### <br>

>  `python3 parse_count_table_confusion_matrix.py`  
>  ~~`-c final.full.count_table (from Mothur QC)`~~  
>  ~~`-f fina.fasta (from Mothur QC)`~~  
>  ~~`-r reference.fasta`~~    
>  ~~`-s sample.csv`~~  
>  ~~`-o output.file`~~    
<br>

**note**:    
1. Need to have blast loaded: ml ncbi-blast+/LATEST   
2. **this script has been updated so that it will take one single argument of config.ini file, and are required arguments are set in that config.ini file, as:** 
>  `python3 parse_count_table_confusion_matrix.py`  
>  `-c config.ini`   

<br>

---

## Extract amplicon sequence using our design primer list on a given isolate WGS   
<br>  

>  ` python3 extract_amplicon_from_primersearch_output.py  `  
>  ` -s isolate_WGS.fasta (or directory name, which holds an array of fasta files) `  
>  ` -p primers-list-psearch.txt`  

***note***  
1. &nbsp;`extract_amplicon_from_primersearch_output.py` will check if the -s argument passed in is a file or a directory. If it's a directory it will grab all the fasta files in the directory.  
2. The primer list has to be in a specific format (tab delimited plain file): total 3 columns, the 1st column is primer name, the 2nd column is forward primer sequence, and the 3rd column is reverse primer sequence. For example: [psearch_primer_list](https://github.com/ncezid-biome/HMAS-QC-Pipeline/blob/master/Salmonella-reformatted-primers-list-psearch.txt) 
3. the script assumes **EMBOSS/6.4.0** is already on path. (`ml EMBOSS/6.4.0` if you're in scicomp space) 
4. The outputs are saved in the **primersearch** folder under the working directory where the script is run.  

<br>

---  

## Download SRA files and assemble isolate WGS     
<br> 

1. Dowload sequence data files using SRA toolkit  
- Installation   
`wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz`  
-  Configure and set it up correctly  
-  fasterq-dump --split-files SRR_file (Scicomp currently has fasterq-dump installed already, so you can skip the first 2 steps )    
2. Use shovill to assemble   
-  ml shovill 
-  this depends, but might need to switch to home folder to use shovill (I sometimes got the folder permission issue when calling shovill)  
-  run shovill  
>  `shovill -R1 SRR1616822_1.fastq -R2 SRR1616822_2.fastq `  
>  `--outdir output_folder --assembler skesa --trim ON --cpus 30`  
***note***: the default spades assembler often throws out out of memory error to me  
3. I have a [run_SRA_assembly.py script](https://github.com/ncezid-biome/HMAS-QC-Pipeline/blob/master/helper_scripts/run_SRA_assembly.py), which automates the whole process, if you already have a list of SRA files to download and assemble.  


>  ` python3   
>  run_SRA_assembly.py  `  
>  ` -i sra-list-file `    

if there is any sra which fails in the process, the script will generate a `sra-list-file_fail_to_assemble`, which you can use to run the script again. The `sra-list-file` is a 2 column (tab delimited) text file: 
 
>sample-1    sra-1  
>sample-2    sra-2   

<br>  

---
## generate pairwise difference matrix at all allele sites (primer pairs) for given sequence fasta files   
<br>    
### <br>

>  `python3 pairwise_diff_matrix.py`  
>  `-o output file`  
>  `-d directory (which holds the fasta files)`  
>  `-p (optional) oligos file (holds primer info)`  
>  `-n (optional) numeric flag`
>  `-y (optional) diff only flag` 
<br>

**note**:  
1. the allele (primer pair) information is in the oligos file, although -p argument is optional, you will need to provide a link to an accessible valid oligos file (I have a default private oligos file in the script)    
2. the sequence fasta file need to have the allele (primer) information in their sequence id. The script use that to locate the pair of sequences to compare. This will not be an issue if you run our `extract_amplicon_from_primersearch_output.py`  
3. if the -n (numeric flag) argument is turned on, a float value `(ex. 0.002)` is used instead of string `(ex. 5 / 2461)` in the output file. If instead the -y (diff only) argument is turned on, only the difference will be in the output `(ex. 5 in this case)`

<br>  

---
## NextFlow scripts to generate a subset of raw reads based on a subset of primers  
<br>

1. Run `create_subset_primers.py` to generate a subset of oligos file. It will generate 2 oligos file (*_subset_oligos file and *_remainder_oligos file, we will use that 'remainder_oligos' file)   
  
2. set up conda env for NextFlow:  
    1. `conda env create -n hmas -f bin/hmas.yaml` (if you have mamba installed, use `mamba env create` instead for speed)   
    2.  `conda activate hmas`    
  

3. Run NextFlow script as:   
>  ` nextflow run hmas2_sampling_rawreads.nf  `  
>  ` --oligo absolute_address_for_that_remainder_oligos file`    

The outputs are all in the **output_sampling_rawreads** folder
<br>   


---
## Generate confusion matrix after running hmas2 QC pipeline 
<br>
### <br>

>  `python3 hmas2_confusion_matrix.py`  
>  `-i hmas2 QC pipeline output folder` (which contains subfolders for each sample)    
>  `-o output confusion_matrix file path`  
>  `-r common reference file for all those samples`  
>  `-m the metasheet file for all those samples`  (this is usually generated while extracting amplicon sequences)    
>  `-p mapping file` (mapping between sample and isolates. A sample might has multiple isolates in it)  
>  `-s the path for parse_count_table_confusion_matrix.py script` 
<br>
