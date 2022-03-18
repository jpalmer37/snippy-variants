#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { fastp } from './modules/snippy.nf'
include { fastp_json_to_csv } from './modules/snippy.nf'
include { snippy } from './modules/snippy.nf'
include { count_variants } from './modules/snippy.nf'
include { qualimap_bamqc } from './modules/snippy.nf'
include { qualimap_bamqc_genome_results_to_csv } from './modules/snippy.nf'

workflow {
  ch_ref = Channel.fromPath( "${params.ref}", type: 'file')

  if (params.samplesheet_input != 'NO_FILE') {
    ch_fastq = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], it['R1'], it['R2']] }
  } else {
    ch_fastq = Channel.fromFilePairs( "${params.fastq_input}/*_R{1,2}*.fastq.gz", type: 'file', maxDepth: 1)
  }

  main:
    fastp(ch_fastq)
    fastp_json_to_csv(fastp.out.fastp_json)
    snippy(fastp.out.reads.combine(ch_ref))
    qualimap_bamqc(snippy.out.alignment)
    qualimap_bamqc_genome_results_to_csv(qualimap_bamqc.out).map{ it -> it[1] }.collectFile(name: 'alignment_qc.csv', keepHeader: true, sort: { it.text }).set{ ch_alignment_qc }
}
