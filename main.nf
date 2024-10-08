#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { fastp }                     from './modules/snippy.nf'
include { find_best_reference }          from './modules/snippy.nf'
include { snippy }                    from './modules/snippy.nf'
include { count_variants }            from './modules/snippy.nf'
include { qualimap_bamqc }            from './modules/snippy.nf'

workflow {

  // samplesheet input
  if (params.samplesheet_input != 'NO_FILE') {
    if (params.ref != 'NO_FILE'){
      ch_fastq = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], it['R1'], it['R2']] }.unique{ it -> it[0] }
    } else {
      ch_fastq = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], it['R1'], it['R2'], it['REF']] }.unique{ it -> it[0] }
    }

  // fastq input 
  } else {
    if (params.ref != 'NO_FILE'){
      ch_fastq = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }
    } else {
      error "ERROR: Reference file must be supplied with --ref when using --fastq_input"
    }

  }
  
  // FORMAT:
  // ch_fastq: [sample_id, r1_fastq_path, r2_fastq_path]
  // ch_ref: [sample_id, reference_path]

  main:

    if (params.ref != 'NO_FILE' && file(params.ref).isFile()){
      println("Detected file as reference input. Using reference file for all samples ...")
      ch_ref = ch_fastq.map{it -> it[0]}.combine(Channel.fromPath( "${params.ref}", type: 'file'))
    }else if (params.ref != 'NO_FILE' && file(params.ref).isDirectory()){
      println("Detected directory as reference input. Finding best reference based on referenceseeker results ...")
      find_best_reference(Channel.fromPath(params.ref))
      ch_ref = ch_fastq.map{it -> it[0]}.combine(find_best_reference.out.fasta)
    }else if (params.ref == 'NO_FILE' && params.samplesheet_input != 'NO_FILE'){
      ch_ref = ch_fastq.map{it -> [it[0], it[3]]}
      ch_fastq = ch_fastq.map{it -> it.subList(0,3)}
    }else {
      error "ERROR: Invalid input."
    }

    fastp(ch_fastq)
    snippy(fastp.out.reads.join(ch_ref))
    qualimap_bamqc(snippy.out.alignment)
    count_variants(snippy.out.variants_csv.join(ch_ref))
}
