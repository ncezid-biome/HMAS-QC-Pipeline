#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.outdir = workflow.launchDir
params.reads = workflow.launchDir


Channel
  .fromFilePairs("${params.reads}/*_R{1,2}*.fastq.gz",size: 2)
 .map{ reads -> tuple(reads[0].replaceAll(~/_S[0-9]+_L[0-9]+/,""), reads[1]) }
  // .view()
  .set { paired_reads }

// Channel
//   .fromPath("${params.reads}/*")
// //   .view()
//   .set { reads }


process deinterleave {
    // publishDir "${params.outdir}", mode: 'copy'
    // tag '$(basename "!{reads}" ".cleaned.fastq.gz")'
    cpus 4
    memory = 5.GB
    // container 'staphb/shovill:latest'
    debug true
    errorStrategy 'retry'
    maxRetries 2

    input:
    path(reads)
    
    output:
    // tuple val("${filename}"), path ("*_R1.fastq.gz"), path ("*_R2.fastq.gz")
    tuple path ("*_R1.fastq.gz"), path ("*_R2.fastq.gz")

    shell:
    '''
        filename=$(basename "!{reads}" ".cleaned.fastq.gz")
        out1="${filename}_R1.fastq.gz"
        out2="${filename}_R2.fastq.gz"

        pigz --best --processes 2 -dc !{reads} | deinterleave_fastq.sh $out1 $out2 compress


    '''



}

process dummy {
    // publishDir "${params.outdir}", mode: 'copy'
    tag "${sample}"
    debug true
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(sample), path(reads)

    output:
    // tuple val(sample), path ("cutadapt/${sample}*.1.fastq"), path ("cutadapt/${sample}*.2.fastq")
    // path("${filename}_assembled.fasta")
    // path("${sample}/${sample}_assembled.fasta")
    tuple val(sample), path(reads)

    shell:
    '''
    echo !{reads}

    '''
}

process shovill {
    publishDir "${params.outdir}", mode: 'copy'
    tag "${sample}"
    // cpus 36
    // memory = 50.GB
    container 'staphb/shovill:latest'
    debug true
    errorStrategy 'retry'
    maxRetries 2

    input:
    // tuple val(sample), path(R1_reads), path(R2_reads)
    // tuple path(R1_reads), path(R2_reads)
    tuple val(sample), path(reads)

    output:
    // tuple val(sample), path ("cutadapt/${sample}*.1.fastq"), path ("cutadapt/${sample}*.2.fastq")
    // path("${filename}_assembled.fasta")
    // path("${sample}/${sample}_assembled.fasta")
    tuple val(sample), path("${sample}/${sample}_assembled.fasta")

    shell:
    '''

    #shovill -R1 !{reads[0]} -R2 !{reads[1]} --outdir !{params.outdir} --force \
    #        --assembler skesa --trim ON --cpus 4 > /dev/null 2>&1
    
    shovill -R1 !{reads[0]} -R2 !{reads[1]} --outdir !{sample} --force \
            --assembler skesa --trim ON --cpus 4 > /dev/null 2>&1


    mv !{sample}/contigs.fa !{sample}/!{sample}_assembled.fasta

    '''

}

process quast {
    publishDir "${params.outdir}/${sample}", mode: 'copy'
    tag "${sample}"
    // container 'staphb/quast:latest'
    container 'https://depot.galaxyproject.org/singularity/quast%3A5.3.0--py39pl5321heaaa4ec_0'

    input:
    tuple val(sample), path (fasta)

    output:
    path ("report_${sample}.tsv")             , optional:true, emit: report

    shell:
    '''
    # comment here
    # cat !{fasta} >> combined_fasta
    quast !{fasta} &&  mv quast_results/latest/report.tsv report_!{sample}.tsv

    '''

}


workflow {
    // paired_reads_ch = deinterleave(reads)
    // paired_reads_ch.view()
    // assembled_reads_ch = shovill(paired_reads_ch)
    dummy_ch = dummy(paired_reads)
    // assembled_reads_ch = shovill(paired_reads)
    assembled_reads_ch = shovill(dummy_ch)
    quast(assembled_reads_ch)

    // clean_reads_ch = concat_reads(removed_primer_reads_ch)
    // extract_rawreads(clean_reads_ch.join(paired_reads))
    
}

