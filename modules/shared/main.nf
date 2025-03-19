process align_all_to_ref {
    container 'CASP_v1.sif'

    publishDir "${params.output_dir}/${params.region}_align/", mode: 'link', saveAs: {filename -> "${sample_name}_${input_type}_$filename"}

    input:
    tuple val(sample_name), val(input_type), path(unalign_bam)

    output:
    tuple val(sample_name), val(input_type), path("full.bam"), path("full.bam.bai"), emit: full_align

    """
    samtools fastq -TMM,ML ${unalign_bam} | minimap2 -x map-ont -y --secondary=no -a -t ${task.cpus} ${params.resource_dir}/${params.region}_align_ref.fasta - | samtools sort - -o full.bam && samtools index full.bam
    """
}

process extract_reads {
    container 'CASP_v1.sif'
    label 'low_resources'
    
    publishDir "${params.output_dir}/${params.region}_reads/", mode: 'link', pattern: "*.{fastq,fai}", saveAs: {filename -> "${sample_name}_${params.region}_$filename"}
    
    input:
    tuple val(sample_name), val(input_type), path(full_bam), path(full_bai)
    
    output:
    tuple val(sample_name), val(input_type), path("${input_type}_region.fastq"), path("${input_type}_region.fastq.fai"), emit: region_reads
    
    """
    samtools view -bh -F 0x100 ${full_bam} -L ${params.resource_dir}/${params.region}_region.bed | samtools fastq -TMM,ML - | seqkit seq -m ${params.extract_min_length} > ${input_type}_region.fastq && samtools faidx ${input_type}_region.fastq
    """
}

process herro_preprocess {
    container 'HERRO_scripts.sif'
      
    input: 
    tuple val(sample_name), val(input_type), path(region_fastq), path(region_fai)

    output:
    tuple val(sample_name), val(input_type), path('herro.fastq.gz'), path('herro_batch_dir'), emit: herro_pre
    
    """
    preprocess.sh ${region_fastq} herro ${task.cpus} 1
    seqkit seq -ni herro.fastq.gz > ids.txt
    create_batched_alignments.sh herro.fastq.gz ids.txt ${task.cpus} herro_batch_dir
    """
}

process herro_inference {
    container 'herro.sif'
    
    publishDir "${params.output_dir}/${params.region}_reads/", mode: 'link', saveAs: {filename -> "${sample_name}_${params.region}_$filename"}

    input:
    tuple val(sample_name), val(input_type), path(herro_fastq), path(herro_batches)
 
    output:
    tuple val(sample_name), val(input_type), path("${input_type}_corr.fasta"), emit: corr_reads

    """
    herro inference --read-alns ${herro_batches} -t ${task.cpus} -d 0 -m ${params.herro_model} -b 64 ${herro_fastq} ${input_type}_corr.fasta
    """
}

process asm_finalize {
    container = 'CASP_v1.sif'

    publishDir "${params.output_dir}/${params.region}_final_contigs/", mode: 'link', pattern: "*.{fasta,fai}", saveAs: {filename -> "${sample_name}_${params.region}_$filename"}
    publishDir "${params.output_dir}/${params.region}_stats/", mode: 'link', pattern: "*stats.tsv", saveAs: {filename -> "${sample_name}_${params.region}_$filename"}

    input:
    tuple val(sample_name), path(hap1_fasta), path(hap1_fai), path(hap2_fasta), path(hap2_fai)

    output:
    tuple val(sample_name), path('hap1_final_contigs.fasta'), path('hap1_final_contigs.fasta.fai'), path('hap2_final_contigs.fasta'), path('hap2_final_contigs.fasta.fai'), emit: final_assembly
    tuple val(sample_name), path('hap_asm_filt_stats.tsv'), emit: asm_finalize_stats

    """
    python ${params.script_dir}/finalize_contigs.py ${hap1_fasta} ${hap2_fasta} ${params.resource_dir}/${params.region}_boundary.fasta hap1_final_contigs.fasta hap2_final_contigs.fasta \
           hap_asm_filt_stats.tsv -f ${params.fasm_min_frac_overlap} -i ${params.fasm_min_identity} -e ${params.fasm_min_edge_distance} -l ${params.fasm_min_overlap_length}  \
           -m ${params.fasm_max_contig_size} -t ${task.cpus}
    samtools faidx hap1_final_contigs.fasta 
    samtools faidx hap2_final_contigs.fasta
    """
}
