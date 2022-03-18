process fastp {

  tag { sample_id }

  input:
  tuple val(grouping_key), path(reads)

  output:
  tuple val(sample_id), path("${sample_id}_fastp.json"), emit: fastp_json
  tuple val(sample_id), path("${sample_id}_trimmed_R1.fastq.gz"), path("${sample_id}_trimmed_R2.fastq.gz"), emit: reads

  script:
  if (grouping_key =~ '_S[0-9]+_') {
    sample_id = grouping_key.split("_S[0-9]+_")[0]
  } else if (grouping_key =~ '_') {
    sample_id = grouping_key.split("_")[0]
  } else {
    sample_id = grouping_key
  }
  """
  fastp -i ${reads[0]} -I ${reads[1]} -o ${sample_id}_trimmed_R1.fastq.gz -O ${sample_id}_trimmed_R2.fastq.gz
  mv fastp.json ${sample_id}_fastp.json
  """
}

process fastp_json_to_csv {

  tag { sample_id }

  executor 'local'

  input:
  tuple val(sample_id), path(fastp_json)

  output:
  tuple val(sample_id), path("${sample_id}_read_count.csv")

  script:
  """
  fastp_json_to_csv.py -s ${sample_id} ${fastp_json} > ${sample_id}_read_count.csv
  """
}

process snippy {

    tag { sample_id }

    publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}", mode: 'copy', pattern: "${sample_id}/${sample_id}*", saveAs: { filename -> filename.split("/").last() }
    publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}", mode: 'copy', pattern: "${sample_id}/reference", saveAs: { filename -> filename.split("/").last() }

    input:
    tuple val(sample_id), file(reads_1), file(reads_2), path(ref)

    output:
    tuple val(sample_id), path("${sample_id}/${sample_id}*"), emit: snippy_files
    tuple val(sample_id), path("${sample_id}/reference", type: 'dir'), emit: snippy_ref
    tuple val(sample_id), path("${sample_id}/${sample_id}.bam"), path("${sample_id}/${sample_id}.bam.bai"), emit: alignment
    
    script:
    """
    snippy \
      --cpus ${task.cpus} \
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

process mosdepth {

    tag { sample_id }
    
    cpus 2

    publishDir "${params.outdir}", mode: 'copy', pattern: "${sample_id}.mosdepth.global.dist.txt"

    input:
    tuple val(sample_id), file(alignment), file(alignment_index)

    output:
    path("${sample_id}.mosdepth.global.dist.txt")
    
    script:
    """
    mosdepth --threads ${task.cpus} --no-per-base ${sample_id} ${alignment}
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
    
    cpus 2

    publishDir "${params.outdir}", mode: 'copy', pattern: "${sample_id}"

    input:
    tuple val(sample_id), file(alignment), file(alignment_index)

    output:
    tuple val(sample_id), path("${sample_id}_bamqc/genome_results.txt"), emit: genome_results
    
    script:
    """
    qualimap bamqc -bam ${alignment} --outdir ${sample_id}_bamqc
    """
}

process qualimap_bamqc_genome_results_to_csv {

  tag { sample_id }

  executor 'local'

  input:
  tuple val(sample_id), path(qualimap_bamqc_genome_results)

  output:
  tuple val(sample_id), path("${sample_id}_qualimap_bamqc_genome_results.csv")

  script:
  """
  qualimap_bamqc_genome_results_to_csv.py -s ${sample_id} ${qualimap_bamqc_genome_results} > ${sample_id}_qualimap_bamqc_genome_results.csv
  """
}

process count_variants {

  tag { sample_id }

  executor 'local'

  input:
  tuple val(sample_id), path(snippy_outdir)

  output:
  tuple val(sample_id), path("${sample_id}_variant_counts.csv")

  script:
  """
  echo 'sample_id,num_snps,num_indel,num_mnp,num_complex' > ${sample_id}_variant_counts.csv
  echo '${sample_id}' > sample_id
  awk -F ',' 'BEGIN { OFS=FS }; \$3 == "snp"' ${snippy_outdir}/${sample_id}.csv | wc -l > num_snps
  awk -F ',' 'BEGIN { OFS=FS }; \$3 == "ins"' ${snippy_outdir}/${sample_id}.csv | wc -l > num_ins
  awk -F ',' 'BEGIN { OFS=FS }; \$3 == "del"' ${snippy_outdir}/${sample_id}.csv | wc -l > num_del
  echo \$(cat num_ins) + \$(cat num_del) | bc > num_indel
  awk -F ',' 'BEGIN { OFS=FS }; \$3 == "mnp"' ${snippy_outdir}/${sample_id}.csv | wc -l > num_mnp
  awk -F ',' 'BEGIN { OFS=FS }; \$3 == "complex"' ${snippy_outdir}/${sample_id}.csv | wc -l > num_complex
  paste -d ',' sample_id num_snps num_indel num_mnp num_complex >> ${sample_id}_variant_counts.csv
  """
}
