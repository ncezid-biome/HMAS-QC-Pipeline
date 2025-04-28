#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.outdir = workflow.launchDir
params.reads = workflow.launchDir
params.primers = null


/*
this workflow run primersearch on the given assembly fasta files, and the primer lists
you need to give the full path for the 3 required parameters
 nextflow run main.nf --primers --outdir --reads

note the primer list file is a 3 column text file (tab delimited) like:
primer_name left_primer_sequence  right_primer_sequence

*/

Channel
  .fromPath("${params.reads}/**/*.fasta")
//   .view()
  .map { file ->
    def sample = file.getBaseName()  // removes extension
    tuple(sample, file)
  }
  .set { reads }

process run_primersearch {
    publishDir "${params.outdir}/${sample}", mode: 'copy'
    tag "${sample}"
    debug true
    errorStrategy 'retry'
    maxRetries 2
    maxForks = 36

    input:
    tuple val(sample), path(fasta_file)
    
    output:
    tuple val(sample), path ("*.ps"), emit: primersearch, optional: true

    shell:
    '''
    run_primersearch.py --primers !{params.primers} \
                        --sequence !{fasta_file} \
                        --output "!{sample}.ps"

    '''

}


process parse_primersearch {
    publishDir "${params.outdir}/${sample}", mode: 'copy'
    tag "${sample}"
    debug true
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(sample), path(fasta_file), path(ps_results)
    
    output:
    path ("*_extractedAmplicons.fasta"), emit: parse_primersearch, optional: true
    path ("*.txt"), emit: empty_primers, optional: true

    shell:
    '''
    parse_primersearch.py --results !{ps_results} --sequence !{fasta_file}

    '''

}


workflow {
    primersearch_ch = run_primersearch(reads)
    joined_ch = reads.join(primersearch_ch)
    parse_primersearch(joined_ch)
    
}
