# Summary of current directory

This directory contains flowcharts I created in the course of working on the HMAS QC Pipeline.
They were all created using LucidChart.

Note: Changes have been made to the pipeline as of December 2020. Some of these visuals are no longer accurate.  
I am unable to modify them because I only have access to the free version of LucidChart.

## Description of contents 

All flowcharts last modified 5/28/2020.

`HMAS_Pipeline_Swimlane_Flowchart.pdf`
This flowchart is an overview of the HMAS Pipeline and the different parts are related to one another.
The QC part of the pipeline (represented in this repository) comprises two software tools: Mothur and VSearch.
After QC, the cleaned amplicons are ready for downstream MLST-style or AMR-style analysis, depending on data type.

** Note: ** Changes as of December 2020 include:
* Search for chimeras and remove chimeras is done in Mothur directly, using vsearch command
* `Trim seqs & de-dedupe again` is removed because the `make.contigs1 is doing the trimming with the added `trimoverlap` option
* The output files generated are: fasta, count, and group files. They have the tag `final` in their filename.

`HMAS_Mothur_Workflow.pdf`
This is a simple flowchart of the Mothur operations used in the HMAS QC workflow.

** Note: ** Changes and notes as of December 2020
* Arrows are missing because the free version of Lucidchart limits the number of objects you can use in one document
* There are several `summary.seqs` processes in the pipeline to add useful info to the log; only one can be listed because of the object limit
* `pcr.seqs` and the following `unique.seqs` are removed
* `make.shared` and `get.otulist` are removed
* After `get.seqs`, `rename.file` is added to tag the output filenames with the keyword `final`
* Outputs are: fasta, count, and group files

`Mothur_IO_Flow.pdf` & `Mothur_IO_Flow_Finish.pdf`
A much more detailed flowchart detailing Mothur I/O.  I wrote it while working through mothur_py, the python module that was used in this pipeline to execute Mothur. Since my team did not have access to the paid version of LucidChart, I had to built it in two parts - thus, the flowchart starts with `Mothur I_O Flow.pdf` and ends in `Mothur I_O Flow Finish.pdf`.  **Note** these two are out-of-date as of December 2020. Changes include:

* The `pcr.seqs` process and the `unique.seqs` process that follows it are no longer needed, and have been removed from the pipeline
* The `get.otulist` and `make.shared` processes are no longer needed, and have been removed from the pipeline
* Resulting flowchart is much simpler: `list.seqs` takes the count file and generates an accnos file, and `get.seqs` takes that accnos file
* `get.seqs` is run twice: 
    * first one takes fasta, accnos, and count files and generates fasta and count files with the subset of sequences indicated in accnos
    * second one takes fasta, accnos, and group files and generates the group file with the subset of sequences indicated in accnos
* The `rename.file` process has been added to give the final fasta, group, count, and list files the tag `final` in the filename so they are easily to locate
* The `Mothur_IO_Flow_Finish.pdf` is basically no longer needed at all

`Pipeline_IO_28052020.pdf`
A flowchart of the pipeline as of 5/28/2020.  
It details a problem of circular referencing between scripts.

## END
