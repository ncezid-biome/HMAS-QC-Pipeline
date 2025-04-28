#!/usr/bin/env nextflow

params.input = null
params.primers = null
params.output_csv = "pairwise_diff_matrix.csv"

/*
this workflow generate the pairwise difference matrix on the given in-silico amplicon sequence fasta files, and the primer lists
you need to give the full path for the 3 required parameters
 nextflow run main.nf --primers --input --output_csv

note the primer list file is a 3 column text file (tab delimited) like:
primer_name left_primer_sequence  right_primer_sequence

*/

Channel
    .fromPath("${params.input}/**/*.fasta")
    .set { individual_fasta }

process runPairwiseRow {
    tag "${query.simpleName}"
    publishDir "${params.outdir}", mode: 'copy'
    maxForks 72

    input:
        path(query)

    output:
        path("${query.simpleName}.txt")

    script:
    """
    pairwise_compare.py \
            --query ${query} \
            --all ${params.input} \
            --primers ${params.primers} \
            --output ${query.simpleName}.txt \
            --diff_only
    """
}

process mergeRows {
    tag "merge"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path row_files

    output:
    path "${params.output_csv}"

    script:
    """
    merge_rows.py \
        -i ${row_files} \
        -o ${params.output_csv}
    """
}

workflow {
    runPairwiseRow(individual_fasta)
    mergeRows(runPairwiseRow.out.collect())
}
