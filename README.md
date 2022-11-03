# HMAS QC Pipeline

A pipeline for doing initial quality control of highly-multiplexed amplicon sequencing data

Started by [@jessicarowell](https://github.com/jessicarowell)
Active development by [@jessicarowell](https://github.com/jessicarowell) and [@jinfinance](https://github.com/jinfinance)

## TOC
* [Description](#description)
* [Requirements](#requirements)
* [INSTALL](#install)
* [USAGE](#usage)
* [Contributing](#contributing)
* [Future Plans](#future-plans)
* [Resources](#resources)

## Description

This is a pipeline that performs quality control analysis on highly-multiplexed amplicon sequencing (HMAS) data.
The tool takes a user-configurable settings file and multiplexed fastq files (gzipped or not), and executes a custom workflow using Mothur.
It provides 2 main outputs of interest: 

1. A **fasta** file containing the high-quality representative unique sequences after cleaning
2. A **count_table** file containing the abundance information of the above fasta file
2. A **summary** file containing summary statistics of reads after each step of the workflow

For more information and to see visualizations describing the workflow, see [this folder containing visuals](https://github.com/ncezid-biome/HMAS-QC-Pipeline/tree/master/visuals).

This pipeline has been designed and tested under Linux CentOS and Ubuntu platforms.  It has not been tested under Windows.

## Requirements

1. Python 3 or higher. Download python [here](https://www.python.org/downloads/). 

2. Mothur must be installed and on your path. Find mothur installation guide [here](https://mothur.org/wiki/installation/).  

	Note: The easiest way to make sure Mothur is on your path is to download the zip file and unzip it in your local bin directory.
        For CDC users, installing it yourself this way is better than using the module version of Mothur.  That version has not been tested here.
        The last version of Mothur that has been tested is ***1.46.0.*** (the more recent version has tweaks that does not fit ino the current pipeline) 

3. You must have the mothur_py package installed.  Read more about mothur-py [here](https://pypi.org/project/mothur-py/).
	`pip install mothur-py`
         Last tested the install in May 2020.

4. Cutadapt must be installed and on your path. Find cutadapt installation guide [here](https://cutadapt.readthedocs.io/en/stable/installation.html).


## INSTALL

(Note: I haven't elaborated here because these instructions will change when we containerize the pipeline.)

1. Copy the Github repository to a folder  
`git clone https://github.com/jessicarowell/HMAS-QC-Pipeline.git` 

2. Add `pipeline.py` to your $PATH


## USAGE

1. Set up your config file. Rename if you want.

	If you open it in a Windows-based text editor, you may need to run something like dos2unix to convert the newlines back to Unix format.
	`dos2unix settings.ini`

2. Check your python version. It should be python 3.  
`python --version`

3. Check that Mothur is installed and on your $PATH.  
`mothur --help`

4. Run the following (replace `mysettings.ini` with the path to your settings file you configured in step 1):  
`python3 pipeline.py - c mysettings.ini`  

## Note

1. It is recommended that you run a quality check on your read sets (e.g. with a program like FastQC) before running them through the pipeline.  Knowing the quality of your read sets may help you troubleshoot any problematic results from the pipeline.

2. We must have both R1(forward) and R2(reverse) reads, and at least I1 index file. If we don't have I2 index file, we
must have keyword NONE or none in the place of missing I2 index file. In such case, we also must have keyword NONE/none
in the corresponding column (for the missing I2 index file) of the oligos file

## Contributing

(Note: this section might also change depending on how we package this.)

Please feel free to fork this repo, make improvements, and share them with me.

Please also post any issues you encounter to Github and I'll be sure to look into them as soon as I can.


## Future Plans

We plan to containerize this pipeline in the future.

## Resources

[Mothur manual](https://mothur.org/wiki/mothur_manual/)
[mothur-py](https://pypi.org/project/mothur-py/)
