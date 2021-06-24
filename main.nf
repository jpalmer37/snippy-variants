#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { fastp } from './modules/snippy.nf'
include { fastp_json_to_csv } from './modules/snippy.nf'
include { snippy } from './modules/snippy.nf'
include { count_variants } from './modules/snippy.nf'
include { qualimap_bamqc } from './modules/snippy.nf'
include { qualimap_bamqc_genome_results_to_csv } from './modules/snippy.nf'
include { join_csvs } from './modules/snippy.nf'

workflow {
  ch_fastq = Channel.fromFilePairs( "${params.fastq_input}/*_R{1,2}*.fastq.gz", type: 'file', maxDepth: 1)
  ch_ref = Channel.fromPath( "${params.ref}", type: 'file')

  main:

  fastp(ch_fastq)
  fastp_json_to_csv(fastp.out.fastp_json).map{ it -> it[1] }.collectFile(name:'read_qc.csv', keepHeader: true, sort: { it.text }).set{ ch_read_qc }
  snippy(fastp.out.reads.combine(ch_ref))
  count_variants(snippy.out.snippy_outdir).map{ it -> it[1] }.collectFile(name: 'variant_counts.csv', keepHeader: true, sort: { it.text }).set{ ch_variant_counts } 
  qualimap_bamqc(snippy.out.alignment)
  qualimap_bamqc_genome_results_to_csv(qualimap_bamqc.out).map{ it -> it[1] }.collectFile(name: 'alignment_qc.csv', keepHeader: true, sort: { it.text }).set{ ch_alignment_qc }
  join_csvs(ch_read_qc.combine(ch_variant_counts).combine(ch_alignment_qc))
}
