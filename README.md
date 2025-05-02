# ONT AS Assembly Pipeline

These pipelines generate _de novo_ haplotypic assemblies across the MHC and KIR using
ONT Adapative Sampling reads corrected with HERRO.  

## Overview

The pipeline is a combination of existing tools (listed below) and custom python scripts 
packaged in a Nextflow pipeline.  Apptainer/Singurity image files containing the necessary
third party tools and the required resource files can be found on Zenodo: XXXXX. 
The pipeline is compatable with HPC environments with Singularity installed and is currently 
set up for the Slurm workload manager.

## Requirements

The minimum node requirement is 12 cpus with 24GB of memory.  All resource requirements can be modified in the configuration file
and might need to be increased if the region is sequenced to a higher depth than presented in the manuscript or if a larger region is 
targeted. 

GPUs with 10GB of memory are required for correction.  The pipeline was tested only on A100 GPUs, but we exepect that GPUs that support
Dorado basecalling should work with the pipeline.
  

## Installation

1. Install [Nextflow](https://www.nextflow.io/docs/latest/install.html)
1. Clone the repository
1. Download the resource bundle (casp_resource_bundle.tar.gz) from Zenodo and extract into the repository directory.  The `resources` and `images` directories 
   should be on the same level as `scripts`
1. Modify the `process` and `singularity` sections of the relevant `*.config` file to match HPC resources.
    * Process: The `executor` variable should be set to your institution's HPC scheduler: https://www.nextflow.io/docs/latest/executor.html. 
      The `clusterOptions` line can be used to set email preferences, cluster queues/nodes or job time 
      limits. Example (slurm): change `--time=96:00:00 --mail-user USER_EMAIL --mail-type FAIL` to `--time 96:00:00 --mail-user jdoe@your.edu --mail-type FAIL`. Note that the herro_inference process has its own `clusterOptions` line, as it requires a GPU node for execution.
      Follow your institution's instructions on how to request GPU nodes: Example (slurm): change `--time=96:00:00 --mail-user USER_EMAIL --mail-type FAIL GPU_SELECTOR` to `--time=96:00:00 --mail-user jdoe@your.edu --mail-type FAIL -p gpuq --gres=gpu:a100:1`. 
    * Singularity: The appropriate singularity bind path should be set using: `runOptions`.  The default setting should work for most people, but see this guide in case it needs to be changed:
      https://docs.sylabs.io/guides/3.0/user-guide/bind_paths_and_mounts.html. Example: `runOptions = "-B /mnt/"`.  

      
## How to Run: CASP MHC/KIR

1. ONT data should be basecalled with Dorado (traditional + 5mC) prior to running the pipeline.  Example basecalling command: `dorado basecaller sup,5mCG_5hmCG pod5/ --min-qscore 10 > SAMPLE_SFE.bam`
1. The CASP assembly pipeline supports SFE or a combination of SFE + UL reads.  Copy or link all basecalled ONT reads in bam format
   into a directory named `sfe_bam` and/or `ul_bam` depending on what was sequenced, one file per sample. SFE bam files should be named: 
   `SAMPLE_SFE.bam`, where `SAMPLE` is a name of your choice.  UL fastq files should be named `SAMPLE_UL.bam`.  The input directory name can be modified from the default `sfe_bam` in the nextflow `*.config` files.
1. The project name can be modified by updating the `project_name` parameter in the config file (e.g.`casp_mhc.config` change `PROJECT_NAME` to `my_mhc_assembly`).  The `input_type` 
   parameter should be set to `sfe` or `sfe_ul` depending on the available reads.
1. The pipeline can be started with the command `nextflow run casp.nf -c casp_mhc.config`
   
## Run Completion

1. On completion, output files can be found in the `project_name` directory.
    * `mhc_align`: alignment file(s) in bam format
    * `mhc_reads`: Raw and corrected MHC-specific reads in fastq/fasta format
    * `mhc_final_contigs`: MHC assemblies in fasta format
    * `mhc_stats`: Pipeline statistics
    * `mhc_phasing`: Phased variant calls in vcf format and phase block definitions in gtf format. 
    * `mhc_genotyping`: HLA typing for each haplotype in gtf format.
1. Interrupted runs can be resumed: `nextflow run casp.nf -c casp_mhc.config -resume`
1. On workflow completion, intermediate files can be removed from the `work` directory

## Haplotype Integrity

Once the assembly is complete, the output files can be checked to determine the likelihood of fully haplotypic assemblies.
1. The most important file to check is `mhc_phasing/SAMPLE_mhc_merge_ref_asm_vcf_stats.tsv`.   Beneath the header 'Block Lists' is a list of the final joined phase blocks after merging reference and assembly based phasing.  If there 
is just a single line below the header, assemblies did not rely on haplotype estimation and are likely to be hapotype resolved.  If there are multiple lines, the file `mhc_phasing/SAMPLE_mhc_merge_phase.gtf` can be used to identify 
the regions of the assembly that are physically phased and the locations of potential phase switches.  Each line of the gtf file has a starting and ending coordinate representing a phase block and if estimation fails, there
could a phase switch between any block. 
1. The pipeline reports the HLA genotypes for each haplotype in the files: `mhc_genotyping/SAMPLE_mhc_hapN_res.gtf.gz`.  The HLA typing for a sample can be entered into [haplostats](https://www.haplostats.org/haplostats?execution=e2s1)
to check for the likelihood of the particular combination of HLA alleles being in phase (haplotype).  If the assembly haplotypes are very infrequent in HaploStats, it could indicate that the assembly contains a potential phase switch.  It should be noted that this is an imperfect
way to check haplotype integrity, as phase switches may occur outside the the region between HLA-A and HLA-DQB1, not all populations are well-represented in HaploStats (i.e South Asian) and the assembled sample may truly have rare haplotypes.
1. Methylation-based phasing is not as well-established as traditional phasing, so users may want to pay special attention to the log file `mhc_stats/SAMPLE_mhc_ref_meth_phase_stats.tsv`.  This file reports methylation
phasing statistics between each reference-based phase block.
   1. `Call` indicates the final methylation phase status.  If the call is 'flip' or 'same', it means the phase blocks were joined with methylation.
   1. `Total` is the total number of reads used to determine the final phase orientation. The minimum number of reads required can be set in the configuration file: `meth_min_reads`. Blocks failing this threshold have the call: 'low_overlap'.
   1. `Diff` is the percentage difference of reads supporting the same phase oriention between phase bocks and the flipped orientation.  Low numbers may indicate bad phasing and the minimum threshold can be set in the configuration file: `meth_min_diff`. Blocks failing this threshold have the call: 'low_difference'.
   1. `B_Bal` is the balance between hap1 or hap2 within the reads used to determine phase orientation. There are two balance numbers `B1_Bal` and `B2_Bal`, one for assignments originating from the upstream phase block and one for the downstream
   block.  If all reads have one haplotype assignment, phase assignment is likely compromised.  The minimum allowed balance can be set in the configuration file: `meth_min_bias`. Blocks failing this threshold have the call: 'biased_haplotyping'.
1. The stats file `mhc_stats/SAMPLE_mhc_hap_read_partition_stats.tsv` contains detailed information about the read partitioning.  The details of the block assignments can be found under the 'Phased Block Assignment' and 'Hemi Bock Assignment' sections. 
If the phasing scaffold does not agree with the initial assembly information, blocks will have the label 'SPLT' and read partitioning will be based on the scaffold.  Samples with 'SPLT' blocks could indicate compromised phasing.

## How to Run: HIFIASM MHC/KIR

1. ONT data should be basecalled with Dorado (traditional) prior to running the pipeline.  Example basecalling command: `dorado basecaller sup pod5/ --min-qscore 10 > SAMPLE_SFE.bam`
1. The hifiasm MHC/KIR assembly pipeline supports SFE-only, UL-only or a combination.  Copy the basecalled ONT reads in bam format
   into the `sfe_bam` and/or `ul_bam` depending on what was seqeunced. SFE bam files should be named: 
   `SAMPLE_SFE.bam`, where `SAMPLE` is a name of your choice.  UL fastq files should be named `SAMPLE_UL.bam`.
1. The project name can be modified by updating the `project_name` paramter in the config file `casp_hifiasm_[kir/mhc].config`. 
   The `input_type` parameter should be set to `sfe`, `ul` or `sfe_ul` depending on the available reads.
1. The pipeline can be started with the command `nextflow run casp_hifism.nf -c casp_hifism_REGION.config`, where REGION is KIR or MHC.


## How to Run: HIFIASM custom

1. Create a custom reference file named `REGION_align_ref.fasta`, where `REGION` is the name of your region of interest. 
   If haplotypic sequences were used as part of the adaptive sampling targeting, they should be added to the hg38 reference.
   If the adaptive sampling only targeted a subset of hg38, no modifications to hg38 are required. 
1. Create the region interval file named `REGION_region.bed`.  This file should list all of the intervals representing your
   region of interest.  If haplotypic sequences were added to hg38, these should be added as intervals in the bed file, 
   otherwise just list the targeted hg38 intervals.
1. Final contigs can be trimmed to specific boundaries using the file `REGION_boundary.fasta`.  This is done to avoid 
   reporting sequence outisde of the targeted region. Boundary information is also used to orient the final assemblies 
   on the same strand.  If this trimming is desired, add two fasta entries to the file named `START` and `END`.  Assemblies 
   will be trimmed upstream of `START` and downstream of `END`.  The sequences should be long enough to avoid ambiguous placement.  
   The sequences of `START` and `END` should be from the same strand and the final assembly will be oriented to the strand of 
   the boundries. If trimming is not desired `REGION_boundary.fasta` can be empty. 
1. Create a new config file based off of the KIR pipeline `cp casp_hifiasm_kir.config casp_hifiasm_REGION.config`. Modify
   the `region` and `genome_size` parameters to match your targeted region.
1. Follow the instructions in `How to Run: HIFIASM MHC / KIR` to assemble your region of interest. 


## External tools

* [minimap2](https://github.com/lh3/minimap2)
* [canu](https://github.com/marbl/canu)
* [shapeit4](https://odelaneau.github.io/shapeit4/)
* [whatshap](https://github.com/whatshap/whatshapcan)
* [longshot](https://github.com/pjedge/longshot)
* [bwa](https://github.com/lh3/bwa)
* [jellyfish](https://github.com/gmarcais/Jellyfish)
* [samtools](https://github.com/samtools/)
* [bedtools](https://github.com/arq5x/bedtools2)
* [HERRO](https://github.com/lbcb-sci/herro)
* [Hifiasm](https://github.com/chhylp123/hifiasm)
* [MethPhaser](https://github.com/treangenlab/methphaser)
* [Immuannot](https://github.com/YingZhou001/Immuannot)
* [HaploStats](https://www.haplostats.org/haplostats?execution=e3s1)

