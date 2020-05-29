# Summary of current directory

This directory contains flowcharts I created in the course of working on the HMAS QC Pipeline.
They were all created using LucidChart.

## Description of contents 

`HMAS_Pipeline_Swimlane_Flowchart.pdf`
This flowchart is an overview of the HMAS Pipeline and the different parts are related to one another.
The QC part of the pipeline (represented in this repository) comprises two software tools: Mothur and VSearch.
After QC, the cleaned amplicons are ready for downstream MLST-style or AMR-style analysis, depending on data type.

`HMAS_Mothur_Workflow.pdf`
This is a simple flowchart of the Mothur operations used in the HMAS QC workflow.

`Mothur_IO_Flow.pdf` & `Mothur_IO_Flow_Finish.pdf`
A much more detailed flowchart detailing Mothur I/O.  I wrote it while working through mothur_py, the python module that was used in this pipeline to execute Mothur. Since my team did not have access to the paid version of LucidChart, I had to built it in two parts - thus, the flowchart starts with `Mothur I_O Flow.pdf` and ends in `Mothur I_O Flow Finish.pdf`..

`Pipeline_IO_28052020.pdf`
A flowchart of the pipeline as of 5/28/2020.  
It details a problem of circular referencing between scripts.

## END
