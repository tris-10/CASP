include { align_all_to_ref} from './modules/shared'
include { extract_reads } from './modules/shared'
include { herro_preprocess } from './modules/shared'
include { herro_inference } from './modules/shared'
include { asm_finalize } from './modules/shared'

process merge_raw_reads {
    container 'CASP_v1.sif'
    label 'low_resources'
    
    publishDir "${params.output_dir}/${params.region}_reads/", mode: 'link', saveAs: {filename -> "${sample_name}_${params.region}_$filename"}
    
    input:
    tuple val(sample_name), path(raw_sfe_fastq), path(raw_ul_fastq)
    
    output:
    tuple val(sample_name), val("comb"), path("comb_region.fastq"), path("comb_region.fastq.fai"),  emit: region_reads
    
    """
    cat ${raw_sfe_fastq} ${raw_ul_fastq} > comb_region.fastq && samtools faidx comb_region.fastq
    """
}

process merge_corr_reads {
    container 'CASP_v1.sif'
    label 'low_resources'
    
    publishDir "${params.output_dir}/${params.region}_reads/", mode: 'link', saveAs: {filename -> "${sample_name}_${params.region}_$filename"}
    
    input:
    tuple val(sample_name), path(corr_sfe_fasta), path(corr_ul_fasta)
    
    output:
    tuple val(sample_name), val("comb"), path("comb_corr.fasta"),  emit: corr_reads
    
    """
    cat ${corr_sfe_fasta} ${corr_ul_fasta} > comb_corr.fasta
    """
}


process ref_snv_call {
    container = 'CASP_v1.sif'
    label 'med_resources'

    publishDir "${params.output_dir}/${params.region}_phasing/", mode: 'link', pattern: "*.{bam,bai,vcf,tbi,gtf}", saveAs: {filename -> "${sample_name}_${params.region}_$filename"}
    publishDir "${params.output_dir}/${params.region}_stats/", mode: 'link', pattern: "*stats.tsv", saveAs: {filename -> "${sample_name}_${params.region}_$filename"}
   
    input:
    tuple val(sample_name), val(sample_type), path(region_fastq), path(region_fai)
    
    output:
    tuple val(sample_name), path("ref_filt.bam"), path("ref_filt.bam.bai"), emit: ref_align
    tuple val(sample_name), path("ref_phase.vcf.gz"), path("ref_phase.vcf.gz.tbi"), emit: ref_phase_vcf
    tuple val(sample_name), path("ref_all.vcf.gz"), path("ref_all.vcf.gz.tbi"), emit: ref_all_vcf
    tuple val(sample_name), path("ref_phase.gtf"), emit: ref_phase_gtf
    path 'ref_align_filt_stats.tsv'
    path 'ref_vcf_filt_stats.tsv'
    
    
    """
    minimap2 -x map-ont --secondary=no -t ${task.cpus} -Y -a ${params.resource_dir}/chr6_drb_refs.fasta ${region_fastq} | samtools sort - -o ref.bam && samtools index ref.bam 
    python ${params.script_dir}/filter_ont_ref_alignment.py ref.bam ref_filt.bam ref_align_filt_stats.tsv -m ${params.region_ref_min_qual} -s ${params.rccx_start} -e ${params.rccx_end} \
            && samtools index ref_filt.bam

    longshot --bam ref_filt.bam --ref ${params.resource_dir}/chr6_drb_mask.fasta --out ref_ls.vcf --min_alt_frac ${params.ls_min_alt_frac} --min_alt_count ${params.ls_min_alt_count} \
             --min_cov ${params.ls_min_cov} --hom_snv_rate ${params.ls_hom_snv_rate} -P ${params.ls_strand} -F -q 0 --region ${params.boundary}
    python ${params.script_dir}/filter_reference_vcf.py ref_ls.vcf ${params.resource_dir}/mhc_gnomad_valid.txt ref_filt.vcf.gz ref_vcf_filt_stats.tsv \
            -m ${params.ls_ref_min_qual} -d ${params.ls_min_cov} -p -s ${params.rccx_start} -e ${params.rccx_end}
    
    whatshap phase --reference ${params.resource_dir}/chr6_drb_mask.fasta -o ref_phase.vcf.gz --ignore-read-groups ref_filt.vcf.gz ref_filt.bam
    tabix -p vcf ref_phase.vcf.gz
    whatshap stats ref_phase.vcf.gz --gtf ref_phase.gtf

    bcftools filter -i 'GT="1/1"' ref_ls.vcf > ref_hom.vcf
    bgzip ref_hom.vcf && tabix -p vcf ref_hom.vcf.gz
    bcftools concat ref_hom.vcf.gz ref_phase.vcf.gz | bcftools sort -T . - -o ref_all.vcf.gz && tabix -p vcf ref_all.vcf.gz
    """
}

process asm_initial {
    container = 'CASP_v1.sif'
    label = 'assembly'

    publishDir "${params.output_dir}/${params.region}_phasing/", mode: 'link', pattern: "*.{fasta,fai}", saveAs: {filename -> "${sample_name}_${params.region}_$filename"}
    publishDir "${params.output_dir}/${params.region}_stats/", mode: 'link', pattern: "*stats.tsv", saveAs: {filename -> "${sample_name}_${params.region}_$filename"}

    input:
    tuple val(sample_name), val(sample_type), path(corr_reads)

    output:
    tuple val(sample_name), path("final_contigs.fasta"), path("final_contigs.fasta.fai"), emit: contigs_fasta
    val "initial_asm_filt_stats.tsv"

    """
    canu maxMemory=24G maxThreads=${task.cpus} useGrid=false genomeSize=${params.genome_size} -untrimmed -pacbio-hifi ${corr_reads} -p mhc_assembly -d canu_output
    python ${params.script_dir}/filter_contigs_by_overlap.py canu_output/mhc_assembly.contigs.fasta final_contigs.fasta initial_asm_filt_stats.tsv -t ${task.cpus} -i ${params.filt_contig_ident} \
        -p ${params.filt_contig_overlap} -m ${params.filt_contig_max_diff} -c ${params.filt_contig_min_length} -l ${params.filt_contig_max_length}
    samtools faidx final_contigs.fasta
    """
}

process asm_snv_call {
    container = 'CASP_v1.sif'
    label = 'med_resources'

    publishDir "${params.output_dir}/${params.region}_phasing/", mode: 'link', pattern: "*.{bam,bai,vcf,tbi,gtf}", saveAs: {filename -> "${sample_name}_${params.region}_$filename"}
    publishDir "${params.output_dir}/${params.region}_stats/", mode: 'link', pattern: "*stats.tsv", saveAs: {filename -> "${sample_name}_${params.region}_$filename"}

    input:
    tuple val(sample_name), path(contigs_fasta), path(contigs_fai), val(sample_type), path(corr_fasta)
   
    output:
    tuple val(sample_name), path("asm_final.vcf.gz"), path("asm_final.vcf.gz.tbi"), emit: asm_vcf
    tuple val(sample_name), path("asm_tagged.bam"), path("asm_tagged.bam.bai"), emit: asm_tagged
    tuple val(sample_name), path("asm_phase.gtf"), emit: asm_gtf
    val "asm_vcf_filt_stats.tsv"
    

    """
    minimap2 -x map-hifi --secondary=no -t ${task.cpu} -a ${contigs_fasta} ${corr_fasta} | samtools sort - -o contig.bam && samtools index contig.bam
    longshot --bam contig.bam --ref ${contigs_fasta} --out asm_ls.vcf --min_alt_frac ${params.ls_min_alt_frac} --min_alt_count ${params.ls_min_alt_count} \
             --min_cov ${params.ls_min_cov} --hom_snv_rate ${params.ls_hom_snv_rate} -P ${params.ls_strand} -F -q 0
    python ${params.script_dir}/split_vcf_by_genotype.py asm_ls.vcf asm_het.vcf.gz asm_hom.vcf.gz asm_vcf_filt_stats.tsv -c ${params.ls_min_alt_count} -q ${params.ls_asm_min_qual} -l ${params.ls_min_hp}
    whatshap phase --reference ${contigs_fasta} -o asm_phase.vcf.gz --ignore-read-groups asm_het.vcf.gz contig.bam
    python ${params.script_dir}/set_isolated_snp_phase.py asm_phase.vcf.gz asm_final.vcf.gz -d ${params.ls_iso_dist}
    whatshap stats --gtf=asm_phase.gtf asm_final.vcf.gz
    whatshap haplotag -r ${contigs_fasta} -o asm_tagged.bam --ignore-read-groups asm_final.vcf.gz contig.bam && samtools index asm_tagged.bam
    """
}

process link_ref_snv_to_asm {
    container = 'CASP_v1.sif'
    label = 'low_resources'

    input:
    tuple val(sample_name), path(contigs_fasta), path(contigs_fai), path(ref_phase_vcf), path(ref_phase_tbi), path(ref_all_vcf), path(ref_all_tbi)

    output:
    tuple val(sample_name), path("snv_up.bam"), path("snv_up.bam.bai"), path("snv_down.bam"), path("snv_down.bam.bai"), emit: snv_flanks

    """
    bcftools consensus -f ${params.resource_dir}/chr6_drb_mask.fasta -o chr6_hap1.fasta -H 1 ${ref_all_vcf}
    bcftools consensus -f ${params.resource_dir}/chr6_drb_mask.fasta -o chr6_hap2.fasta -H 2 ${ref_all_vcf}

    python ${params.script_dir}/create_het_flank_seq.py ${ref_phase_vcf} chr6_hap1.fasta h1_up.fasta h1_down.fasta -l ${params.phase_merge_flank_length}
    python ${params.script_dir}/create_het_flank_seq.py ${ref_phase_vcf} chr6_hap2.fasta h2_up.fasta h2_down.fasta -l ${params.phase_merge_flank_length}

    cat h1_up.fasta h2_up.fasta | seqkit rmdup -s > up.fasta
    cat h1_down.fasta h2_down.fasta | seqkit rmdup -s > down.fasta

    bwa index ${contigs_fasta}
    bwa mem -a ${contigs_fasta} up.fasta | samtools sort -o snv_up.bam && samtools index snv_up.bam
    bwa mem -a ${contigs_fasta} down.fasta | samtools sort -o snv_down.bam && samtools index snv_down.bam
    """
}

process update_ref_phase_with_asm {
    container = 'CASP_v1.sif'
    label = 'low_resources'

    publishDir "${params.output_dir}/${params.region}_phasing/", mode: 'link', pattern: "*.{vcf,tbi,gtf}", saveAs: {filename -> "${sample_name}_${params.region}_$filename"}
    publishDir "${params.output_dir}/${params.region}_stats/", mode: 'link', pattern: "*stats.tsv", saveAs: {filename -> "${sample_name}_${params.region}_$filename"}

    input:
    tuple val(sample_name), path(contig_fasta), path(contig_fai), path(ref_vcf), path(ref_tbi), path(ref_gtf), 
          path(asm_vcf), path(asm_tbi), path(up_bam), path(up_bai), path(down_bam), path(down_bai)

    output:
    tuple val(sample_name), path("merge_phase.vcf.gz"), path("merge_phase.vcf.gz.tbi"), emit: merge_phase_vcf
    path "merge_phase.gtf"
    path "vcf_merge_stats.tsv"

    """
    python ${params.script_dir}/update_phasing_with_assembly.py ${contig_fasta} ${asm_vcf} ${up_bam} ${down_bam} ${ref_vcf} ${ref_gtf} merge_phase.vcf.gz vcf_merge_stats.tsv \
           -m ${params.phase_merge_min_het_pset} -f ${params.phase_merge_min_het_frac} -n ${params.phase_merge_max_probe_mm} -p ${params.phase_merge_max_probe_align}
    whatshap stats merge_phase.vcf.gz --gtf merge_phase.gtf
    """
}

process predict_remaining_phase {
    container = 'CASP_v1.sif'
    label = 'low_resources'

    publishDir "${params.output_dir}/${params.region}_phasing/", mode: 'link', pattern: "*.{vcf.gz,tbi}", saveAs: {filename -> "${sample_name}_${params.region}_$filename"}
    publishDir "${params.output_dir}/${params.region}_stats/", mode: 'link', pattern: "*stats.tsv", saveAs: {filename -> "${sample_name}_${params.region}_$filename"}

    input:
    tuple val(sample_name), path(merge_phase_vcf), path(merge_phase_tbi)

    output:
    tuple val(sample_name), path("shapeit.vcf.gz"), path("shapeit.vcf.gz.tbi"), emit: shapeit_vcf
    path "vcf_common_stats.tsv"

    """
    python ${params.script_dir}/filter_reference_vcf.py ${merge_phase_vcf} ${params.resource_dir}/mhc_gnomad_valid_common.txt common.vcf.gz vcf_common_stats.tsv
    shapeit4.2 -I common.vcf.gz -H ${params.resource_dir}/chr6_1kGP_phased_SNV.bcf -R chr6 \
           --map ${params.resource_dir}/chr6.b38.gmap.gz -T 5 --pbwt-depth 8 -O shapeit.vcf.gz --use-PS 1e-20 --mcmc-iterations 10b,1p,1b,1p,1b,1p,1b,1p,10m --sequencing
    tabix -p vcf shapeit.vcf.gz
    """
}

process predict_remaining_phase_trio {
    container = 'CASP_v1.sif'
    label = 'low_resources'

    publishDir "${params.output_dir}/${params.region}_phasing/", mode: 'link', pattern: "*.{vcf.gz,tbi}", saveAs: {filename -> "${sample_name}_${params.region}_$filename"}
    publishDir "${params.output_dir}/${params.region}_stats/", mode: 'link', pattern: "*stats.tsv", saveAs: {filename -> "${sample_name}_${params.region}_$filename"}

    input:
    tuple val(sample_name), path(merge_phase_vcf), path(merge_phase_tbi), path(trio_vcf), path(trio_tbi)

    output:
    tuple val(sample_name), path("shapeit.vcf.gz"), path("shapeit.vcf.gz.tbi"), emit: shapeit_vcf
    path "vcf_common_stats.tsv"

    """
    python ${params.script_dir}/filter_reference_vcf.py ${merge_phase_vcf} ${params.resource_dir}/mhc_gnomad_valid_common.txt common.vcf.gz vcf_common_stats.tsv
    shapeit4.2 -I common.vcf.gz -H ${params.resource_dir}/chr6_1kGP_phased_SNV.bcf -R chr6 \
           --map ${params.resource_dir}/chr6.b38.gmap.gz -T 5 --pbwt-depth 8 -O shapeit.vcf.gz --use-PS 1e-20 --mcmc-iterations 10b,1p,1b,1p,1b,1p,1b,1p,10m \
           --sequencing --scaffold ${trio_vcf}
    tabix -p vcf shapeit.vcf.gz
    """
}

process partition_reads {
    container = 'CASP_v1.sif'
    label = 'low_resources'

    publishDir "${params.output_dir}/${params.region}_reads/", mode: 'link', pattern: "*.{fasta,fai}", saveAs: {filename -> "${sample_name}_${params.region}_$filename"}
    publishDir "${params.output_dir}/${params.region}_stats/", mode: 'link', pattern: "*stats.tsv", saveAs: {filename -> "${sample_name}_${params.region}_$filename"}


    input:
    tuple val(sample_name), val(sample_type), path(herro_fasta), path(contigs_fasta), path(contigs_fai), path(tagged_bam), path(tagged_bai),
         path(asm_vcf), path(asm_tbi), path(asm_gtf), path(shapeit_vcf), path(shapeit_tbi), path(up_bam), path(up_bai), path(down_bam), path(down_bai) 
     
    output:
    tuple val(sample_name), path('hap1_reads.fasta'), path('hap1_reads.fasta.fai'), path('hap2_reads.fasta'), path('hap2_reads.fasta.fai'), emit: hap_reads
    path 'hap_read_partition_stats.tsv'

    """
    python ${params.script_dir}/partition_haplotypes.py ${contigs_fasta} ${tagged_bam} ${shapeit_vcf} ${asm_vcf} ${asm_gtf} ${up_bam} ${down_bam} ${herro_fasta} \
        partition_out hap_read_partition_stats.tsv -m ${params.pt_min_gt_count} -b ${params.pt_min_contam_split} -c ${params.pt_max_contam_assign} -d ${params.ls_min_cov} \
        -x ${params.phase_merge_max_probe_mm} -p ${params.phase_merge_max_probe_align}
        cat partition_out_hap1_only.fasta partition_out_unknown.fasta > hap1_reads.fasta && samtools faidx hap1_reads.fasta
        cat partition_out_hap2_only.fasta partition_out_unknown.fasta > hap2_reads.fasta && samtools faidx hap2_reads.fasta
    """
}

process asm_haplotype_specific {
    container = 'CASP_v1.sif'
    label = 'assembly'

    input:
    tuple val(haplotype), val(sample_name), path(hap_fasta), path(hap_fai)

    output:
    tuple  val(sample_name), path("${haplotype}_asm/${haplotype}_canu_output.contigs.fasta"), path("${haplotype}_asm/${haplotype}_canu_output.contigs.fasta.fai"), emit: hap_assembly

    """
    canu maxMemory=24G maxThreads=${task.cpus} useGrid=false genomeSize=${params.genome_size} -untrimmed -pacbio-hifi ${hap_fasta} -p ${haplotype}_canu_output -d ${haplotype}_asm
    samtools faidx "${haplotype}_asm/${haplotype}_canu_output.contigs.fasta"
    """
}


workflow {

    // Extract and correct MHC Reads
    if (params.input_type == 'sfe') {
        fastq_ch = Channel.fromPath("${params.sfe_dir}/*_SFE.fastq.gz", checkIfExists: true).map{tuple(it.getSimpleName().split("_SFE")[0], "sfe", it)}

        align_all_to_ref(fastq_ch) 
            | extract_reads
            | herro_preprocess
            | herro_inference
    
        raw_reads = extract_reads.out.region_reads
        corr_reads = herro_inference.out.corr_reads
    
    } else if (params.input_type == 'sfe_ul') {
        sfe_ch = Channel.fromPath("${params.sfe_dir}/*fastq.gz", checkIfExists: true).map{tuple(it.getSimpleName().split("_SFE")[0], "sfe", it)}
        ul_ch = Channel.fromPath("${params.ul_dir}/*fastq.gz", checkIfExists: true).map{tuple(it.getSimpleName().split("_UL")[0], "ul", it)}
    
        sfe_ch.join(ul_ch)
            .map{ sample, sfe_type, sfe_fastq, ul_type, ul_fastq -> [tuple(sample, sfe_type, sfe_fastq), tuple(sample, ul_type, ul_fastq)]}.flatMap()
            | align_all_to_ref
            | extract_reads
            | herro_preprocess
            | herro_inference
        
        extract_reads.out.branch {
            sfe: it[1] == 'sfe'
                return tuple(it[0],it[2])
            ul: it[1] == 'ul'
                return tuple(it[0], it[2])
        }.set { ind_raw_reads }
        
        herro_inference.out.branch {
            sfe: it[1] == 'sfe'
                return tuple(it[0],it[2])
            ul: it[1] == 'ul'
                return tuple(it[0],it[2])
        }.set { ind_corr_reads }
        
        raw_reads = merge_raw_reads(ind_raw_reads.sfe.join(ind_raw_reads.ul))
        corr_reads = merge_corr_reads(ind_corr_reads.sfe.join(ind_corr_reads.ul))
             
    } else {
        error("Invalid input type: parameter 'input_type' must be set to 'sfe', 'sfe_ul' or 'ul'")
    }
    
    // Reference variant calling
    ref_snv_call(raw_reads)
    
    // Assembly variant calling
    asm_initial(corr_reads)
    asm_snv_call(asm_initial.out.contigs_fasta
        .join(corr_reads))

    // Merge variant calls
    link_ref_snv_to_asm(asm_initial.out.contigs_fasta
        .join(ref_snv_call.out.ref_phase_vcf)
        .join(ref_snv_call.out.ref_all_vcf))
    update_ref_phase_with_asm(asm_initial.out.contigs_fasta
        .join(ref_snv_call.out.ref_phase_vcf)
        .join(ref_snv_call.out.ref_phase_gtf)
        .join(asm_snv_call.out.asm_vcf)
        .join(link_ref_snv_to_asm.out.snv_flanks))

    // Predict any remaining phase or use trio data
    if (params.phasing_method == 'shapeit') {
        final_vcf_ch = predict_remaining_phase(update_ref_phase_with_asm.out.merge_phase_vcf).shapeit_vcf
    } else if (params.phasing_method == 'trio') {
        scaf_vcf_ch =  Channel.fromPath("${params.trio_dir}/*_TRIO.vcf.gz", checkIfExists: true).map{tuple(it.getSimpleName().split("_TRIO")[0], it)}
        scaf_tbi_ch =  Channel.fromPath("${params.trio_dir}/*_TRIO.vcf.gz.tbi", checkIfExists: true).map{tuple(it.getSimpleName().split("_TRIO")[0], it)}
        final_vcf_ch = predict_remaining_phase_trio(update_ref_phase_with_asm.out.merge_phase_vcf
            .join(scaf_vcf_ch)
            .join(scaf_tbi_ch)).shapeit_vcf
    } else {
        error("Invalid phasing method: parameter 'phasing_method' must be set to 'shapeit' or 'trio'")
    }
    
    // Parition Reads
    partition_reads(corr_reads
        .join(asm_initial.out.contigs_fasta)
        .join(asm_snv_call.out.asm_tagged)
        .join(asm_snv_call.out.asm_vcf)
        .join(asm_snv_call.out.asm_gtf)
        .join(final_vcf_ch)
        .join(link_ref_snv_to_asm.out.snv_flanks))

    // Haplotype-specific assembly
    split_hap_reads = partition_reads.out.hap_reads
        .map({sample, fa1, fai1, fa2, fai2 -> [tuple("hap1", sample, fa1, fai1), tuple("hap2", sample, fa2, fai2)]})
        .flatMap()
    asm_haplotype_specific(split_hap_reads)
    hap_asm = asm_haplotype_specific.out.hap_assembly
        .groupTuple(size: 2)
        .map{sample, asm_fasta, asm_fai -> return asm_fasta[0].getSimpleName().startsWith('hap1') 
           ? tuple(sample, asm_fasta[0], asm_fai[0], asm_fasta[1], asm_fai[1]) 
           : tuple(sample, asm_fasta[1], asm_fai[1], asm_fasta[0], asm_fai[0])}
    asm_finalize(hap_asm)
    
}
