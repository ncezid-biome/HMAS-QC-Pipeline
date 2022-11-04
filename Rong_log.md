# Rong's daily log

A place where Rong logs his *significant* daily code activity 

Started by [@jinfinance](https://github.com/jinfinance)

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