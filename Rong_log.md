# Rong's log

A place where Rong logs his *significant* code activity and comments

Started by [@jinfinance](https://github.com/jinfinance)

*~~11/4/2022~~*

updated the **parse_count_table_confusion_matrix.py** script, so that:  
>
>  1. if we have the metasheet table, we will be able to retrieve all the predicted primers for each sample through the metasheet table, instead of using: *run_grep.get_primers_for_sample_design_fasta* function
>  2. *run_grep.get_primers_for_sample_design_fasta* function is still there in case we run the script without the metasheet table
>  - **do note that due to 2.) above, build_confusion_matrix() method is a little heavy now, and might need some simplification**
>  3. I tested the script with both 024 data which does not require map (sample-to-isolate) file, nor metasheet table; and 027-2 data which does have the map file and metasheet table.
> 
>  `python3 parse_count_table_confusion_matrix.py`  
>
>  `  -c /scicomp/home-pure/qtl7/test/hmas_test/output1.0_nseq9_noemptyfastq/juno.final.full.count_table `  
>
>  `  -s helper_scripts/M3235-22-024-sample.csv `  
>
>  `  -f /scicomp/home-pure/qtl7/test/hmas_test/output1.0_nseq9_noemptyfastq/juno.final.fasta ` 
>
>  ` -r ~/HMAS-Internal-Validation/M347_21_026_nopre/juno_design_primers_full_full.fasta `  
>
>  ` -o confusion_matrix_nseq9_oldref_noemptyfastq_1104_1`
>
>  ---
>
>  `python3 parse_count_table_confusion_matrix.py`  
>
>  `  -c /scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/CIMS/HMAS_QC_pipeline/M3235_22_027_2/output1.0/juno.final.full.count_table `  
>
>  `  -f /scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/CIMS/HMAS_QC_pipeline/M347_21_027_2/output1.0/juno.final.fasta ` 
>
>  ` -r 2014K_0324_08_0900.fasta  -e metasheet.csv -m sample_isolates_mapping.csv `
>
>  `  -s M3235_22_027_2.csv `  
>
>  ` -o confusion_matrix_2014K_0324_08_0900_new027 `
>

---
---

*~~11/8/2022~~*

**Salmonella isolates WGS assembly**  


Wrote a script `run_SRA_assembly.py` which takes a 2 column (tab delimited) plain file, 1st column is the sample name, 2nd column is the SRA number. It downloads all SRA files in the list and call shovill for the assembly.

The shovill program is from the loaded module at Scicomp, and I choose the skesa as the assembler, which currently is of version ***v2.3***. 

I ran it with the following 6 Salmonella isolates substitutes:  

2014K-0979  
2014K-0527  
~~2013K-1067~~  
2014K-0421  
2016K-0167  
2016K-0878  
2013K-1828   
*~~AM39224 (AR-0409)~~*

I was not able to find the SRA files for *2013K-1067* and *AR-0409*, and Grant didn't know where to find the WGS for them either.

Except for *2013K-1067*, NCBI already has the WGS assembly there. So it's probably not necessary to download the SRA reads and do the assembly ourselves.  And there are some difference between NCBI version of the assembly and ours.  They also use skesa for the assembly, but they used the version ***v2.2***.  Not sure if that made the difference. For example: sample 2016K-0878, their assembly has ***27*** contigs while our version has ***25*** contigs

[Salmonella enterica strain 2016K-0878, whole genome shotgun sequencing project](https://www.ncbi.nlm.nih.gov/nuccore/AAEFQS000000000)  

***note***  
shovill sometimes would freeze for no obvious reason, without any apparent errors. I set a 10 minutes' timeout in the script. If that ever happens, it will be necessary to re-run the script again with those affected samples/sra files  

---
---


