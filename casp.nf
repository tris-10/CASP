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

    publishDir "${params.output_dir}/${params.region}_phasing/", mode: 'link', pattern: "*.{bam,bai,vcf.gz,tbi,gtf}", saveAs: {filename -> "${sample_name}_${params.region}_$filename"}
    publishDir "${params.output_dir}/${params.region}_stats/", mode: 'link', pattern: "*stats.tsv", saveAs: {filename -> "${sample_name}_${params.region}_$filename"}
   
    input:
    tuple val(sample_name), val(sample_type), path(region_fastq), path(region_fai)
    
    output:
    tuple val(sample_name), path("ref.bam"), path("ref.bam.bai"), emit: ref_align
    tuple val(sample_name), path("ref_phase.vcf.gz"), path("ref_phase.vcf.gz.tbi"), emit: ref_phase_vcf
    tuple val(sample_name), path("ref_meth_phase.vcf.gz"), path("ref_meth_phase.vcf.gz.tbi"), emit: ref_meth_phase_vcf
    tuple val(sample_name), path("ref_all.vcf.gz"), path("ref_all.vcf.gz.tbi"), emit: ref_all_vcf
    tuple val(sample_name), path("ref_phase.gtf"), emit: ref_phase_gtf
    tuple val(sample_name), path("ref_meth_phase.gtf"), emit: ref_meth_phase_gtf
    tuple val(sample_name), path('ref_vcf_filt_stats.tsv'), path('ref_meth_phase_stats.tsv'), emit: ref_snv_call_stats
    
    
    """
    minimap2 -x map-ont -y --secondary=no -t ${task.cpus} -Y -a ${params.resource_dir}/${params.region}_snv_ref.fasta ${region_fastq} | samtools view -b -q ${params.region_ref_min_qual} | samtools sort - -o ref.bam && samtools index ref.bam 
    longshot --bam ref.bam --ref ${params.resource_dir}/${params.region}_snv_mask.fasta --out ref_ls.vcf --min_alt_frac ${params.ls_min_alt_frac} --min_alt_count ${params.ls_min_alt_count} \
             --min_cov ${params.ls_min_cov} --hom_snv_rate ${params.ls_hom_snv_rate} -P ${params.ls_strand} -F -q 0 --region ${params.boundary}
    python ${params.script_dir}/filter_reference_vcf.py ref_ls.vcf ${params.resource_dir}/${params.region}_gnomad_valid.txt ref_filt.vcf.gz ref_vcf_filt_stats.tsv -m ${params.ls_ref_min_qual} -d ${params.ls_min_cov} -p 
    
    whatshap phase --reference ${params.resource_dir}/${params.region}_snv_mask.fasta -o ref_phase.vcf.gz --ignore-read-groups ref_filt.vcf.gz ref.bam
    tabix -p vcf ref_phase.vcf.gz
    whatshap stats ref_phase.vcf.gz --gtf ref_phase.gtf

    if [ "\$(cut -f 1 ref_phase.gtf | uniq | wc -l)" -eq \$(cat ref_phase.gtf | wc -l) ]; then
        cp ref_phase.vcf.gz ref_meth_phase.vcf.gz
        cp ref_phase.vcf.gz.tbi ref_meth_phase.vcf.gz.tbi
        cp ref_phase.gtf ref_meth_phase.gtf
        echo 'one phase block, skipping cpg phasing' > ref_meth_phase_stats.tsv
    else 
        whatshap haplotag -r ${params.resource_dir}/${params.region}_snv_mask.fasta --ignore-read-groups --skip-missing-contigs ref_phase.vcf.gz ref.bam | samtools view -bF 2304 -o ref_tagged.bam - && samtools index ref_tagged.bam
        meth_phaser_parallel -ml -2 -b ref_tagged.bam -r ${params.resource_dir}/${params.region}_snv_mask.fasta -g ref_phase.gtf -vc ref_phase.vcf.gz -o methphaser_out/work -a ${params.meth_min_assign}
        meth_phaser_post_processing -ib ref_tagged.bam -if methphaser_out/work -ov ref_meth_phase_temp.vcf -ob mp_tagged -vc ref_phase.vcf.gz -mc ${params.meth_min_reads} -vd ${params.meth_min_diff} -mb ${params.meth_min_bias}
        bcftools sort ref_meth_phase_temp.vcf > ref_meth_phase.vcf
        bgzip ref_meth_phase.vcf && tabix -p vcf ref_meth_phase.vcf.gz
        whatshap stats ref_meth_phase.vcf.gz --gtf ref_meth_phase.gtf
        python ${params.script_dir}/parse_methphaser_stats.py methphaser_out/work/${params.chrom}/*csv ref_meth_phase_stats.tsv --min_reads ${params.meth_min_reads} --min_percentage ${params.meth_min_diff} --min_bias ${params.meth_min_bias}
    fi
    
    bcftools filter -i 'GT="1/1"' ref_ls.vcf > ref_hom.vcf
    bgzip ref_hom.vcf && tabix -p vcf ref_hom.vcf.gz
    bcftools concat ref_hom.vcf.gz ref_meth_phase.vcf.gz | bcftools sort -T . - -o ref_all.vcf.gz && tabix -p vcf ref_all.vcf.gz
    
    """
}

process asm_initial {
    container = 'CASP_v1.sif'

    publishDir "${params.output_dir}/${params.region}_phasing/", mode: 'link', pattern: "*.{fasta,fai}", saveAs: {filename -> "${sample_name}_${params.region}_$filename"}
    publishDir "${params.output_dir}/${params.region}_stats/", mode: 'link', pattern: "*stats.tsv", saveAs: {filename -> "${sample_name}_${params.region}_$filename"}

    input:
    tuple val(sample_name), val(sample_type), path(corr_reads)

    output:
    tuple val(sample_name), path("final_contigs.fasta"), path("final_contigs.fasta.fai"), emit: contigs_fasta
    tuple val(sample_name), path("raw_contigs.fasta"), path("raw_contigs.fasta.fai"), emit: raw_fasta
    tuple val(sample_name), path("initial_asm_filt_stats.tsv"), emit: asm_initial_stats

    """
    canu maxMemory=${task.memory.giga}G maxThreads=${task.cpus} useGrid=false genomeSize=${params.genome_size} -untrimmed -pacbio-hifi ${corr_reads} -p ${params.region}_assembly -d canu_output
    cp canu_output/${params.region}_assembly.contigs.fasta raw_contigs.fasta && samtools faidx raw_contigs.fasta
    python ${params.script_dir}/filter_contigs_by_overlap.py raw_contigs.fasta final_contigs.fasta initial_asm_filt_stats.tsv -t ${task.cpus} -i ${params.filt_contig_ident} \
        -p ${params.filt_contig_overlap} -m ${params.filt_contig_max_diff} -c ${params.filt_contig_min_length} -l ${params.filt_contig_max_length}
    samtools faidx final_contigs.fasta
    """
}

process asm_snv_call {
    container = 'CASP_v1.sif'

    publishDir "${params.output_dir}/${params.region}_phasing/", mode: 'link', pattern: "*.{bam,bai,vcf.gz,tbi,gtf}", saveAs: {filename -> "${sample_name}_${params.region}_$filename"}
    publishDir "${params.output_dir}/${params.region}_stats/", mode: 'link', pattern: "*stats.tsv", saveAs: {filename -> "${sample_name}_${params.region}_$filename"}

    input:
    tuple val(sample_name), path(contigs_fasta), path(contigs_fai), val(sample_type), path(region_fasta), path(region_fai)
   
    output:
    tuple val(sample_name), path("asm_phase.vcf.gz"), path("asm_phase.vcf.gz.tbi"), emit: asm_phase_vcf
    tuple val(sample_name), path("asm_meth_phase.vcf.gz"), path("asm_meth_phase.vcf.gz.tbi"), emit: asm_meth_phase_vcf
    tuple val(sample_name), path("asm_meth_tagged.bam"), path("asm_meth_tagged.bam.bai"), emit: asm_meth_tagged
    tuple val(sample_name), path("asm_phase.gtf"), emit: asm_gtf
    tuple val(sample_name), path("asm_meth_phase.gtf"), emit: asm_meth_gtf
    tuple val(sample_name), path("asm_vcf_filt_stats.tsv"), path("asm_meth_phase_stats.tsv"), emit: asm_snv_call_stats
    

    """
    minimap2 -x map-ont -y --secondary=no -t ${task.cpus} -a ${contigs_fasta} ${region_fasta} | samtools sort - -o contig.bam && samtools index contig.bam
    longshot --bam contig.bam --ref ${contigs_fasta} --out asm_ls.vcf --min_alt_frac ${params.ls_min_alt_frac} --min_alt_count ${params.ls_min_alt_count} \
             --min_cov ${params.ls_min_cov} --hom_snv_rate ${params.ls_hom_snv_rate} -P ${params.ls_strand} -F -q 0
    python ${params.script_dir}/split_vcf_by_genotype.py asm_ls.vcf asm_het.vcf.gz asm_hom.vcf.gz asm_vcf_filt_stats.tsv -c ${params.ls_min_alt_count} -q ${params.ls_asm_min_qual} -l ${params.ls_min_hp}
    whatshap phase --reference ${contigs_fasta} -o asm_phase_temp.vcf.gz --ignore-read-groups asm_het.vcf.gz contig.bam
    python ${params.script_dir}/set_isolated_snp_phase.py asm_phase_temp.vcf.gz asm_phase.vcf.gz -d ${params.ls_iso_dist}
    whatshap stats --gtf=asm_phase.gtf asm_phase.vcf.gz
    whatshap haplotag -r ${contigs_fasta} --ignore-read-groups asm_phase.vcf.gz contig.bam | samtools view -bF 2304 -o asm_tagged.bam - && samtools index asm_tagged.bam 

    if [ "\$(cut -f 1 asm_phase.gtf | uniq | wc -l)" -eq \$(cat asm_phase.gtf | wc -l) ]; then
        cp asm_phase.vcf.gz asm_meth_phase.vcf.gz
        cp asm_phase.vcf.gz.tbi asm_meth_phase.vcf.gz.tbi
        cp asm_phase.gtf asm_meth_phase.gtf
        cp asm_tagged.bam asm_meth_tagged.bam
        cp asm_tagged.bam.bai asm_meth_tagged.bam.bai
        echo 'one phase block, skipping cpg phasing' > asm_meth_phase_stats.tsv
    else 
        meth_phaser_parallel -ml -2 -b asm_tagged.bam -r ${contigs_fasta} -g asm_phase.gtf -vc asm_phase.vcf.gz -o methphaser_out/work -a ${params.meth_min_assign}
        meth_phaser_post_processing -ib asm_tagged.bam -if methphaser_out/work -ov asm_meth_phase_temp.vcf -ob asm_meth_tagged_temp -vc asm_phase.vcf.gz -mc ${params.meth_min_reads} -vd ${params.meth_min_diff} -mb ${params.meth_min_bias}
        bcftools sort asm_meth_phase_temp.vcf > asm_meth_phase.vcf
        bgzip asm_meth_phase.vcf && tabix -p vcf asm_meth_phase.vcf.gz
        whatshap haplotag -r ${contigs_fasta} --ignore-read-groups asm_meth_phase.vcf.gz contig.bam | samtools view -bF 2304 -o asm_meth_tagged.bam - && samtools index asm_meth_tagged.bam
        whatshap stats --gtf=asm_meth_phase.gtf asm_meth_phase.vcf.gz
        python ${params.script_dir}/parse_methphaser_stats.py methphaser_out/work/*/*csv asm_meth_phase_stats.tsv --min_reads ${params.meth_min_reads} --min_percentage ${params.meth_min_diff} --min_bias ${params.meth_min_bias}
    fi
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
    bcftools consensus -f ${params.resource_dir}/${params.region}_snv_mask.fasta -o ${params.chrom}_hap1.fasta -H 1 ${ref_all_vcf}
    bcftools consensus -f ${params.resource_dir}/${params.region}_snv_mask.fasta -o ${params.chrom}_hap2.fasta -H 2 ${ref_all_vcf}

    python ${params.script_dir}/create_het_flank_seq.py ${ref_phase_vcf} ${params.chrom}_hap1.fasta h1_up.fasta h1_down.fasta -l ${params.phase_merge_flank_length}
    python ${params.script_dir}/create_het_flank_seq.py ${ref_phase_vcf} ${params.chrom}_hap2.fasta h2_up.fasta h2_down.fasta -l ${params.phase_merge_flank_length}

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

    publishDir "${params.output_dir}/${params.region}_phasing/", mode: 'link', pattern: "*.{vcf.gz,tbi,gtf}", saveAs: {filename -> "${sample_name}_${params.region}_$filename"}
    publishDir "${params.output_dir}/${params.region}_stats/", mode: 'link', pattern: "*stats.tsv", saveAs: {filename -> "${sample_name}_${params.region}_$filename"}

    input:
    tuple val(sample_name), path(contig_fasta), path(contig_fai), path(ref_vcf), path(ref_tbi), path(ref_gtf), 
          path(asm_vcf), path(asm_tbi), path(asm_vcf_no_meth), path(meth_asm_tbi_no_meth), path(up_bam), path(up_bai), path(down_bam), path(down_bai)

    output:
    tuple val(sample_name), path("merge_phase.vcf.gz"), path("merge_phase.vcf.gz.tbi"), emit: merge_phase_vcf
    tuple val(sample_name), path("merge_phase.gtf"), path("merge_ref_asm_vcf_stats.tsv"), emit: merge_stats
   
    """
    python ${params.script_dir}/update_phasing_with_assembly.py ${contig_fasta} ${asm_vcf} ${up_bam} ${down_bam} ${ref_vcf} ${ref_gtf} merge_phase.vcf.gz merge_ref_asm_vcf_stats.tsv \
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
    tuple val(sample_name), path("shapeit_filter_stats.tsv"), emit: shapeit_stats

    """
    python ${params.script_dir}/filter_reference_vcf.py ${merge_phase_vcf} ${params.resource_dir}/${params.region}_gnomad_valid_common.txt common.vcf.gz shapeit_filter_stats.tsv
    
    if [ "\$(zcat common.vcf.gz | grep -v "^#" | wc -l)" -eq 0 ]; then
        echo "No heterozygous variants, skipping estimation"
        cp common.vcf.gz shapeit.vcf.gz
    else
        shapeit4.2 -I common.vcf.gz -H ${params.resource_dir}/${params.chrom}_1kGP_phased_SNV.bcf -R ${params.chrom} \
            --map ${params.resource_dir}/${params.chrom}.b38.gmap.gz -T 5 --pbwt-depth 8 -O shapeit.vcf.gz --use-PS 1e-20 --mcmc-iterations 10b,1p,1b,1p,1b,1p,1b,1p,10m --sequencing 
    fi
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
    tuple val(sample_name), path("shapeit_filter_stats.tsv"), emit: shapeit_stats

    """
    python ${params.script_dir}/filter_reference_vcf.py ${merge_phase_vcf} ${params.resource_dir}/${params.region}_gnomad_valid_common.txt common.vcf.gz shapeit_filter_stats.tsv
    shapeit4.2 -I common.vcf.gz -H ${params.resource_dir}/${params.chrom}_1kGP_phased_SNV.bcf -R ${params.chrom} \
           --map ${params.resource_dir}/${params.chrom}.b38.gmap.gz -T 5 --pbwt-depth 8 -O shapeit.vcf.gz --use-PS 1e-20 --mcmc-iterations 10b,1p,1b,1p,1b,1p,1b,1p,10m \
           --sequencing --scaffold ${trio_vcf}
    tabix -p vcf shapeit.vcf.gz
    """
}

process partition_reads {
    container = 'CASP_v1.sif'

    publishDir "${params.output_dir}/${params.region}_reads/", mode: 'link', pattern: "*.{fasta,fai}", saveAs: {filename -> "${sample_name}_${params.region}_$filename"}
    publishDir "${params.output_dir}/${params.region}_stats/", mode: 'link', pattern: "*stats.tsv", saveAs: {filename -> "${sample_name}_${params.region}_$filename"}


    input:
    tuple val(sample_name), val(sample_type), path(corr_reads), path(contigs_fasta), path(contigs_fai), path(tagged_bam), path(tagged_bai),
         path(asm_vcf), path(asm_tbi), path(asm_gtf), path(shapeit_vcf), path(shapeit_tbi), path(up_bam), path(up_bai), path(down_bam), path(down_bai) 
     
    output:
    tuple val(sample_name), path('hap1_reads.fasta'), path('hap1_reads.fasta.fai'), path('hap2_reads.fasta'), path('hap2_reads.fasta.fai'), emit: hap_reads
    tuple val(sample_name), path('hap_read_partition_stats.tsv'), emit: partition_stats

    """
    python ${params.script_dir}/partition_haplotypes.py ${contigs_fasta} ${tagged_bam} ${shapeit_vcf} ${asm_vcf} ${asm_gtf} ${up_bam} ${down_bam} ${corr_reads} \
        partition_out hap_read_partition_stats.tsv -m ${params.pt_min_gt_count} -b ${params.pt_min_contam_split} -c ${params.pt_max_contam_assign} -d ${params.ls_min_cov} \
        -x ${params.phase_merge_max_probe_mm} -p ${params.phase_merge_max_probe_align} -a ${params.pt_min_hemi_hets}
        cat partition_out_hap1_only.fasta partition_out_unknown.fasta > hap1_reads.fasta && samtools faidx hap1_reads.fasta
        cat partition_out_hap2_only.fasta partition_out_unknown.fasta > hap2_reads.fasta && samtools faidx hap2_reads.fasta
    """
}

process asm_haplotype_specific {
    container = 'CASP_v1.sif'
    
    input:
    tuple val(haplotype), val(sample_name), path(hap_fasta), path(hap_fai)

    output:
    tuple  val(sample_name), path("${haplotype}_assembly.fasta"), path("${haplotype}_assembly.fasta.fai"), emit: hap_assembly

    """
    canu maxMemory=${task.memory.giga}G maxThreads=${task.cpus} useGrid=false genomeSize=${params.genome_size} -untrimmed -pacbio-hifi ${hap_fasta} -p ${haplotype}_canu_output -d ${haplotype}_asm obtErrorRate=0.01
    cp "${haplotype}_asm/${haplotype}_canu_output.contigs.fasta" ${haplotype}_assembly.fasta
    samtools faidx ${haplotype}_assembly.fasta
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

    // Extract and correct MHC Reads
    if (params.input_type == 'sfe') {
        sfe_ch = Channel.fromPath("${params.sfe_dir}/*_SFE.bam", checkIfExists: true).map{tuple(it.getSimpleName().split("_SFE")[0], "sfe", it)}

        align_all_to_ref(sfe_ch) 
            | extract_reads
            | herro_preprocess
            | herro_inference
 
        raw_reads = extract_reads.out.region_reads
        corr_reads = herro_inference.out.corr_reads
    
    } else if (params.input_type == 'sfe_ul') {
        sfe_ch = Channel.fromPath("${params.sfe_dir}/*_SFE.bam", checkIfExists: true).map{tuple(it.getSimpleName().split("_SFE")[0], "sfe", it)}
        ul_ch = Channel.fromPath("${params.ul_dir}/*_UL.bam", checkIfExists: true).map{tuple(it.getSimpleName().split("_UL")[0], "ul", it)}
    
        sfe_ch.join(ul_ch)
            .map{ sample, sfe_type, sfe_bam, ul_type, ul_bam -> [tuple(sample, sfe_type, sfe_bam), tuple(sample, ul_type, ul_bam)]}.flatMap()
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
        .join(raw_reads))

    // Merge variant calls
    link_ref_snv_to_asm(asm_initial.out.contigs_fasta
        .join(ref_snv_call.out.ref_phase_vcf)
        .join(ref_snv_call.out.ref_all_vcf))
    update_ref_phase_with_asm(asm_initial.out.contigs_fasta
        .join(ref_snv_call.out.ref_meth_phase_vcf)
        .join(ref_snv_call.out.ref_meth_phase_gtf)
        .join(asm_snv_call.out.asm_meth_phase_vcf)
        .join(asm_snv_call.out.asm_phase_vcf)
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
        .join(asm_snv_call.out.asm_meth_tagged)
        .join(asm_snv_call.out.asm_meth_phase_vcf)
        .join(asm_snv_call.out.asm_meth_gtf)
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

    // HLA genotypes
    hap_contigs = asm_finalize.out.final_assembly
        .map({sample, fa1, fai1, fa2, fai2 -> [tuple("hap1", sample, fa1, fai1), tuple("hap2", sample, fa2, fai2)]})
        .flatMap()
    genotyping(hap_contigs)
    
}
