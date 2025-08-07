#!/usr/bin/env nextflow
nextflow.enable.dsl=2
import java.nio.file.Paths


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
  .fromPath(["${params.reads}/**/*.{fasta,fsa,fs,fna,fa}", "${params.reads}/*.{fasta,fsa,fs,fna,fa}"])
  .map { file ->
    def sample = file.getBaseName()
    def inputFolder = file.getParent().getName()           // e.g., SRR30637285
    def outputFolder = Paths.get(params.outdir, "${inputFolder}_assembled")
    tuple(sample, file, outputFolder)
  }
  .filter { sample, file, outputFolder ->
    !outputFolder.exists()
  }
  .map { sample, file, outputFolder ->
    tuple(sample, file)
  }
  .set { reads }



process run_primersearch {
    publishDir "${params.outdir}/${sample}", mode: 'copy'
    tag "${sample}"
    debug true
    errorStrategy 'retry'
    maxRetries 2
    maxForks = 32

    input:
    tuple val(sample), path(fasta_file)
    
    output:
    tuple val(sample), path ("*.ps"), emit: primersearch, optional: true

    shell:
    '''
    run_primersearch.py --primers !{params.primers} \
                        --sequence !{fasta_file} \
                        --output "!{sample}.ps" \
                        --mismatch !{params.mismatchpercent}


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
    tuple val(sample), path ("*_extractedAmplicons.fasta"), emit: parse_primersearch, optional: true
    path ("*.txt"), emit: empty_primers, optional: true

    shell:
    '''
    parse_primersearch.py --results !{ps_results} --sequence !{fasta_file} --amp_len !{params.max_amplicon_len}

    '''

}


process fasta_to_json {
    publishDir "${params.outdir}/${sample}", mode: 'copy'
    tag "${sample}"
    debug true
    // errorStrategy 'retry'
    // maxRetries 2

    input:
    tuple val(sample), path(fasta_file)
    
    output:
    path ("*.json"), emit: json_file, optional: true

    shell:
    '''
    fasta_to_json.py --fasta_file !{fasta_file} --sample_id !{sample} --primers !{params.primers}

    '''

}


workflow {
    primersearch_ch = run_primersearch(reads)
    joined_ch = reads.join(primersearch_ch)
    amplicon_fasta_ch = parse_primersearch(joined_ch).parse_primersearch
    fasta_to_json(amplicon_fasta_ch)
    
}
