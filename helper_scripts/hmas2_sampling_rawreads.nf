#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.outdir = workflow.launchDir + '/output_sampling_rawreads'
params.reads = workflow.launchDir + '/data'
// params.reads = '/scicomp/home-pure/qtl7/test/hmas_test/024_demulx_0mis_data/output_sampling_rawreads'
// params.oligo = workflow.launchDir + '/data' + '/M3235_22_024.oligos'
params.oligo = workflow.launchDir + '/data' + '/remainder_2461.oligos'

Channel
  .fromFilePairs("${params.reads}/*_R{1,2}*.fastq.gz",size: 2)
 .map{ reads -> tuple(reads[0].replaceAll(~/_S[0-9]+_L[0-9]+/,""), reads[1]) }
//   .view()
  .set { paired_reads }

// Channel
//   .fromFilePairs("${params.reads}/*.{1,2}.fastq",size: 2)
// //   .map{ reads -> tuple(reads[0].replaceAll(~/_S[0-9]+_L[0-9]+/,""), reads[1]) }
// //    .view()
//   .set { paired_reads2 }

// Channel
//     .fromPath(params.primers)
//     .splitCsv(header:false, sep:"\t")
//     .set{primer_ch}


process cutadapt {
    // publishDir "${params.outdir}", mode: 'copy'
    tag "${sample}"
    cpus 20
    memory = 2.GB
    // container 'dceoy/cutadapt:latest'
    // debug true
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(sample), path(reads)
    
    output:
    tuple val(sample), path ("cutadapt/${sample}*.1.fastq"), path ("cutadapt/${sample}*.2.fastq")

    shell:
    '''
    mkdir -p cutadapt
    #python !{workflow.projectDir}/bin/run_cutadapt.py -f !{reads[0]} -r !{reads[1]} \
    #                                 -o cutadapt -s !{sample} -p !{params.oligo}
    run_cutadapt.py -f !{reads[0]} -r !{reads[1]} \
                                     -o cutadapt -s !{sample} -p !{params.oligo}
    '''

}

process concat_reads {
    publishDir "${params.outdir}/${sample}", mode: 'copy'
    tag "${sample}"
    // debug true

    input:
    tuple val(sample), path (reads1), path (reads2)

    output:
    tuple val(sample), path ("${sample}.1.fastq"), path ("${sample}.2.fastq")

    shell:
    '''
    cat !{reads1} >> !{sample}.1.fastq
    cat !{reads2} >> !{sample}.2.fastq

    '''

}

process extract_rawreads {
    publishDir "${params.outdir}/${sample}", mode: 'copy'
    tag "${sample}"
    cpus = 2
    debug = true

    input:
    tuple val(sample), path(reads1), path(reads2), path(rawreads)

    // output:
    // // path ("${sample}.fastq")
    // tuple val(sample), path ("${sample}.fastq")

    shell:
    '''
    #pear -f !{reads1} -r !{reads2} -o !{sample} -q 26 -m 325 -v 20 -j 20
    #mv !{sample}.assembled.fastq !{sample}.fastq
    #vsearch -fastq_mergepairs !{reads1} -reverse !{reads2} -fastqout !{sample}.fastq \
    #                --fastq_maxns 0 --fastq_minovlen 20 --fastq_maxdiffs 20 --fastq_maxmergelen 325

    extract_rawreads_by_seqid.py -i !{reads1} -r !{rawreads[0]} -o !{params.outdir}/!{sample}_S1_L001_R1.fastq.gz
    extract_rawreads_by_seqid.py -i !{reads2} -r !{rawreads[1]} -o !{params.outdir}/!{sample}_S1_L001_R2.fastq.gz

                                     
    '''
}

workflow {
    removed_primer_reads_ch = cutadapt(paired_reads)
    clean_reads_ch = concat_reads(removed_primer_reads_ch)
    extract_rawreads(clean_reads_ch.join(paired_reads))
    
}

