#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { fastp }          from './modules/snippy.nf'
include { snippy }         from './modules/snippy.nf'
include { count_variants } from './modules/snippy.nf'
include { qualimap_bamqc } from './modules/snippy.nf'

workflow {



  // samplesheet input
  if (params.samplesheet_input != 'NO_FILE') {

    // if --ref file provided, use this reference for all entries
    if (params.ref != 'NO_FILE'){
      ch_fastq = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], it['R1'], it['R2']] }
      ch_ref = ch_fastq.map{it -> it[0]}.combine(Channel.fromPath( "${params.ref}", type: 'file'))

    // if no --ref file provided, look for REF column in the samplesheet on a per-sample basis
    } else {
      ch_fastq = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], it['R1'], it['R2'], it['REF']] }
      ch_ref = ch_fastq.map{it -> [it[0], it[3]]}
      ch_fastq = ch_fastq.map{it -> it.subList(0,3)}
    }

  // fastq input 
  } else {
    if (params.ref != 'NO_FILE'){
      ch_fastq = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }
      ch_ref = ch_fastq.map{it -> it[0]}.combine(Channel.fromPath( "${params.ref}", type: 'file'))
    } else {
      error "ERROR: Reference file must be supplied with --ref when using --fastq_input"
    }

  }
  
  // FORMAT:
  // ch_fastq: [sample_id, r1_fastq_path, r2_fastq_path]
  // ch_ref: [sample_id, reference_path]

  main:
    fastp(ch_fastq)
    snippy(fastp.out.reads.join(ch_ref))
    qualimap_bamqc(snippy.out.alignment)
    count_variants(snippy.out.variants_csv.join(ch_ref))
}
