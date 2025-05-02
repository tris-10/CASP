include { align_all_to_ref } from './modules/shared'
include { extract_reads } from './modules/shared'
include { herro_preprocess } from './modules/shared'
include { herro_inference } from './modules/shared'
include { asm_finalize } from './modules/shared'

process read_trimming {
    container = 'CASP_v1.sif'
    label = 'low_resources'

    publishDir "${params.output_dir}/${params.region}_reads/", mode: 'link', pattern: "*.{fasta,fai}", saveAs: {filename -> "${sample_name}_${params.region}_$filename"}

    input:
    tuple val(sample_name), val(input_type), path(corr_reads)

    output:
    tuple val(sample_name), path("${input_type}_trim_corr.fasta"), path("${input_type}_trim_corr.fasta.fai")

    """
    jellyfish count -L 2 -C -m 31 -t ${task.cpus} -s 50000 --output jelly.bin ${corr_reads}
    jellyfish dump -c -o jelly.txt jelly.bin
    python ${params.script_dir}/kmer_trimming.py ${corr_reads} jelly.txt ${input_type}_trim_corr.fasta -l ${params.trim_min_length} -c ${params.trim_frac_cov} -f ${params.trim_frac_low_kmer}
    samtools faidx ${input_type}_trim_corr.fasta
    """
}

process hifiasm_sfe {
    container = 'CASP_v1.sif'
    label = 'assembly'

    publishDir "${params.output_dir}/${params.region}_asm/", mode: 'link', pattern: "*.{fasta,fai}", saveAs: {filename -> "${sample_name}_${params.region}_$filename"}

    input:
    tuple val(sample_name), path(trimmed_fasta), path(trimmed_fai)

    output:
    tuple val(sample_name), path('hap1_raw_contigs.fasta'), path('hap1_raw_contigs.fasta.fai'), path('hap2_raw_contigs.fasta'), path('hap2_raw_contigs.fasta.fai')

    """
    depth=`cat ${trimmed_fai} | awk '{tot += \$2 } END {print int(tot / (${params.genome_size} * 1000000))}'` 
    hifiasm -o hifiasm_out -t ${task.cpus} --hg-size ${params.genome_size}m -f0 -n ${params.min_unitig_reads} --hom-cov \${depth} ${trimmed_fasta}
    awk '/^S/{{print ">"\$2; print \$3}}' hifiasm_out.bp.hap1.p_ctg.gfa | fold > hap1_raw_contigs.fasta && samtools faidx hap1_raw_contigs.fasta
    awk '/^S/{{print ">"\$2; print \$3}}' hifiasm_out.bp.hap2.p_ctg.gfa | fold > hap2_raw_contigs.fasta && samtools faidx hap2_raw_contigs.fasta
    """
}

process hifiasm_ul {
    container = 'CASP_v1.sif'
    label = 'assembly'

    publishDir "${params.output_dir}/${params.region}_asm/", mode: 'link', pattern: "*.{fasta,fai}", saveAs: {filename -> "${sample_name}_${params.region}_$filename"}

    input:
    tuple val(sample_name), path(trimmed_fasta), path(trimmed_fai), val(sample_type), path(ul_fasta), path(ul_fai)

    output:
    tuple val(sample_name), path('hap1_raw_contigs.fasta'), path('hap1_raw_contigs.fasta.fai'), path('hap2_raw_contigs.fasta'), path('hap2_raw_contigs.fasta.fai')

    """
    depth=`cat ${trimmed_fai} | awk '{tot += \$2 } END {print int(tot / (${params.genome_size} * 1000000))}'` 
    hifiasm -o hifiasm_out -t ${task.cpus} --hg-size ${params.genome_size}m -f0 -n ${params.min_unitig_reads} --hom-cov \${depth} ${trimmed_fasta} --ul ${ul_fasta}
    awk '/^S/{{print ">"\$2; print \$3}}' hifiasm_out.bp.hap1.p_ctg.gfa | fold > hap1_raw_contigs.fasta && samtools faidx hap1_raw_contigs.fasta
    awk '/^S/{{print ">"\$2; print \$3}}' hifiasm_out.bp.hap2.p_ctg.gfa | fold > hap2_raw_contigs.fasta && samtools faidx hap2_raw_contigs.fasta
    """
}

process split_ul_reads {
    container = 'CASP_v1.sif'
    label = "low_resources"

    publishDir "${params.output_dir}/${params.region}_reads/", mode: 'link', pattern: "*.{fasta,fai}", saveAs: {filename -> "${sample_name}_${params.region}_$filename"}

    input: 
    tuple val(sample_name), path(raw_fasta), path(raw_fai)

    output:
    tuple val(sample_name), path("split.fasta"), path("split.fasta.fai")

    """
    python ${params.script_dir}/split_ul_reads.py ${raw_fasta} split.fasta -l ${params.split_ul_length}
    samtools faidx split.fasta
    """
}

process genotyping {
    container = 'CASP_v1.sif'
    label = 'low_resources'

    publishDir "${params.output_dir}/${params.region}_genotyping/", mode: 'link', pattern: "*.gtf.gz", saveAs: {filename -> "${sample_name}_${params.region}_$filename"}

    input:
    tuple val(haplotype), val(sample_name), path(hap_fasta), path(hap_fai)

    output:
    tuple val(sample_name), path("${haplotype}_res.gtf.gz"), emit: geno_gtf

    """
    sh /Immuannot/scripts.pub.v3/immuannot.sh -r /Immuannot/Data-2024Feb02 -c ${hap_fasta} -o ${haplotype}_res -t ${task.cpus}
    """
}



workflow {
    if (params.input_type == 'sfe') {
        sfe_ch = Channel.fromPath("${params.sfe_dir}/*_SFE.bam", checkIfExists: true).map{tuple(it.getSimpleName().split("_SFE")[0], "sfe", it)}
        
        align_all_to_ref(sfe_ch) 
            | extract_reads 
            | herro_preprocess 
            | herro_inference 
            | read_trimming 
            | hifiasm_sfe 
            | asm_finalize

    } else if (params.input_type == 'ul') {
        ul_ch = Channel.fromPath("${params.ul_dir}/*_UL.bam", checkIfExists: true).map{tuple(it.getSimpleName().split("_UL")[0], "ul", it)}
            
        align_all_to_ref(ul_ch) 
            | extract_reads 
            | herro_preprocess 
            | herro_inference 
            | read_trimming 
            | split_ul_reads 
            | join(extract_reads.out) 
            | hifiasm_ul 
            | asm_finalize
        
    } else if (params.input_type == 'sfe_ul') {
        sfe_ch = Channel.fromPath("${params.sfe_dir}/*_SFE.bam", checkIfExists: true).map{tuple(it.getSimpleName().split("_SFE")[0], "sfe", it)}
        ul_ch = Channel.fromPath("${params.ul_dir}/*_UL.bam", checkIfExists: true).map{tuple(it.getSimpleName().split("_UL")[0], "ul", it)}

        sfe_ch.join(ul_ch) 
            .map{ sample, sfe_type, sfe_bam, ul_type, ul_bam -> [tuple(sample, sfe_type, sfe_bam), tuple(sample, ul_type, ul_bam)]}.flatMap()
            | align_all_to_ref
            | extract_reads
        
        extract_reads.out.branch {
            sfe: it[1] == 'sfe'
                return tuple(it[0],it[2])
            ul: it[1] == 'ul'
                return tuple(it[0], it[2])
            }.set { region_reads }

        herro_preprocess(region_reads.sfe) 
            | herro_inference 
            | read_trimming 
            | join(region_reads.ul) 
            | hifiasm_ul 
            | asm_finalize

    } else {
        error("Invalid input type: parameter 'input_type' must be set to 'sfe', 'sfe_ul' or 'ul'")
    }

    // HLA genotypes
    hap_contigs = asm_finalize.out.final_assembly
        .map({sample, fa1, fai1, fa2, fai2 -> [tuple("hap1", sample, fa1, fai1), tuple("hap2", sample, fa2, fai2)]})
        .flatMap()
    genotyping(hap_contigs)
}



