# HMAS QC Pipeline

A pipeline for doing initial quality control of highly-multiplexed amplicon sequencing data

Started by [@jessicarowell](https://github.com/jessicarowell)

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
It provides 3 main outputs of interest: 

1. A **fasta** file containing the high-quality sequences after cleaning
2. A **shared** file containing a count of reads passing QC at the end of the workflow
3. An **otulist** file that maps the OTU number assigned by Mothur to its representative sequence.  
	These are not OTU numbers in the traditional sense, for this type of analysis, but rather function as identifiers.

For more information and to see visualizations describing the workflow, see [this folder containing visuals](https://github.com/jessicarowell/HMAS-QC-Pipeline/tree/master/visuals).


## Requirements

1. Python 3 or higher. Download python [here](https://www.python.org/downloads/). 

2. Mothur must be installed and on your path. Find mothur installation guide [here](https://mothur.org/wiki/installation/).

3. You must have the mothur_py package installed.  Read more about mothur-py [here](https://pypi.org/project/mothur-py/).
	`pip install mothur-py`



## INSTALL

Note: I haven't elaborated here because these instructions will change when we containerize the pipeline.

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


## Contributing

Note: this also might change depending on how we package this??

Please feel free to fork this repo, make improvements, and share them with me.

Please also post any issues you encounter to Github and I'll be sure to look into them as soon as I can.


## Future Plans

We plan to containerize this pipeline in the future.

## Resources

[Mothur manual](https://mothur.org/wiki/mothur_manual/)
[mothur-py](https://pypi.org/project/mothur-py/)
