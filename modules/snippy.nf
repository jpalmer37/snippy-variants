process fastp {

  tag { sample_id }

  publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_fastp.{json,csv}"
  
  input:
  tuple val(sample_id), path(reads_r1), path(reads_r2)

  output:
  tuple val(sample_id), path("${sample_id}_fastp.json"),                                                    emit: fastp_json
  tuple val(sample_id), path("${sample_id}_fastp.csv"),                                                     emit: fastp_csv
  tuple val(sample_id), path("${sample_id}_trimmed_R1.fastq.gz"), path("${sample_id}_trimmed_R2.fastq.gz"), emit: reads

  script:
  """
  fastp \
    --thread ${task.cpus} \
    -i ${reads_r1} \
    -I ${reads_r2} \
    --cut_tail \
    -o ${sample_id}_trimmed_R1.fastq.gz \
    -O ${sample_id}_trimmed_R2.fastq.gz \
    -j ${sample_id}_fastp.json
  fastp_json_to_csv.py -s ${sample_id} ${sample_id}_fastp.json > ${sample_id}_fastp.csv
  """
}


process find_best_reference {

    executor 'local'
  
    publishDir "${params.outdir}", mode: 'copy', pattern: "reference.*"

    input:
    path(refseeker_results_path)

    output:
    path("reference.fa"),                                               emit: fasta
    path("reference.gb"),                                               emit: gb
    
    script:
    """
    find_best_reference.py --input ${refseeker_results_path} --outname reference
    """

}

process snippy {

    tag { sample_id }

    publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}/${sample_id}*", saveAs: { filename -> filename.split("/").last() }
    publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}/reference", saveAs: { filename -> filename.split("/").last() }

    input:
    tuple val(sample_id), file(reads_1), file(reads_2), path(ref)

    output:
    tuple val(sample_id), path("${sample_id}/${sample_id}*"),                                               emit: snippy_files
    tuple val(sample_id), path("${sample_id}/reference", type: 'dir'),                                      emit: snippy_ref
    tuple val(sample_id), path("${sample_id}/${sample_id}.bam"), path("${sample_id}/${sample_id}.bam.bai"), emit: alignment
    tuple val(sample_id), path("${sample_id}/${sample_id}.csv"),                                            emit: variants_csv
    
    script:
    ram = task.memory ? ram = "--ram ${task.memory}" : "" 
    """
    mkdir tmp

    snippy \
      --tmpdir ./tmp \
      --cpus ${task.cpus} \
      ${ram} \
      --report \
      --prefix ${sample_id} \
      --mincov ${params.mincov} \
      --basequal ${params.basequal} \
      --mapqual ${params.mapqual} \
      --minfrac ${params.minfrac} \
      --ref ${ref} \
      --R1 ${reads_1} \
      --R2 ${reads_2} \
      --outdir '${sample_id}'
    """
}


process samtools_stats_summary {

    tag { sample_id }
    
    cpus 2

    publishDir "${params.outdir}", mode: 'copy', pattern: "${sample_id}"

    input:
    tuple val(sample_id), file(alignment), file(alignment_index)

    output:
    path("${sample_id}_samtools_stats_summary.txt")
    
    script:
    """
    samtools stats ${alignment} | grep ^SN | cut -f 2- > ${sample_id}_samtools_stats_summary.txt
    """
}

process qualimap_bamqc {

    tag { sample_id }

    publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}/", mode: 'copy', pattern: "${sample_id}_bamqc*"

    input:
    tuple val(sample_id), file(alignment), file(alignment_index)

    output:
    tuple val(sample_id), path("${sample_id}_qualimap_bamqc_genome_results.csv"), emit: genome_results_csv
    tuple val(sample_id), path("${sample_id}_bamqc"),                             emit: qualimap_bamqc_dir
    
    script:
    """
    qualimap bamqc -bam ${alignment} --outdir ${sample_id}_bamqc
    qualimap_bamqc_genome_results_to_csv.py -s ${sample_id} ${sample_id}_bamqc/genome_results.txt > ${sample_id}_qualimap_bamqc_genome_results.csv
    """
}

process count_variants {

  tag { sample_id }

  executor 'local'

  publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_variant_counts.csv"

  input:
  tuple val(sample_id), path(variants_csv), path(ref)

  output:
  tuple val(sample_id), path("${sample_id}_variant_counts.csv")

  script:
  """
  echo 'sample_id,ref_id,num_snps,num_indel,num_mnp,num_complex' > ${sample_id}_variant_counts.csv
  echo '${sample_id}' > sample_id
  head -n 1 ${ref} | tr -d  '>' | cut -d ' ' -f 1 > ref_id
  awk -F ',' 'BEGIN { OFS=FS }; \$3 == "snp"' ${variants_csv} | wc -l > num_snps
  awk -F ',' 'BEGIN { OFS=FS }; \$3 == "ins"' ${variants_csv} | wc -l > num_ins
  awk -F ',' 'BEGIN { OFS=FS }; \$3 == "del"' ${variants_csv} | wc -l > num_del
  echo \$(cat num_ins) + \$(cat num_del) | bc > num_indel
  awk -F ',' 'BEGIN { OFS=FS }; \$3 == "mnp"' ${variants_csv} | wc -l > num_mnp
  awk -F ',' 'BEGIN { OFS=FS }; \$3 == "complex"' ${variants_csv} | wc -l > num_complex
  paste -d ',' sample_id ref_id num_snps num_indel num_mnp num_complex >> ${sample_id}_variant_counts.csv
  """
}
