# ONT AS Assembly Pipeline

These pipelines generate _de novo_ haplotypic assemblies across the MHC and KIR using
ONT Adapative Sampling reads corrected with HERRO.  

## Overview

The pipeline is a combination of existing tools (listed below) and custom python scripts 
packaged in a Nextflow pipeline.  Apptainer/Singurity image files containing the necessary
third party tools and the required resource files can be found on Zenodo: XXXXX. 
The pipeline is compatable with HPC environments with Singularity installed and is currently 
set up for the Slurm workload manager.

## Installation

1. Install [Nextflow](https://www.nextflow.io/docs/latest/install.html)
1. Clone the repository
1. Download the resource bundle from Zenodo and extract.  The `resources` and `images` directories
   should be on the same level as `scripts`
1. Modify the `process` and `singularity` sections of the relevant `*.config` file to match HPC resources.
    * Process: The `executor` variable should be set to the HPC scheduler: https://www.nextflow.io/docs/latest/executor.html. 
      The `clusterOptions` line can be used to set email preferences, cluster queues/nodes or job time 
      limits. Example (slurm): `clusterOptions= '--mail-user xxx@xxx.edu --mail-type FAIL --time 96:00:00'`. Note that the herro_inference process has its own `clusterOptions` line, as it requires a GPU node for execution.
      Follow your institution's instructions on how to request GPU nodes: Example (slurm):  `clusterOptions= '--mail-user xxx@xxx.edu -p gpuq --gres=gpu:a100:1`. 
    * Singularity: The appropriate singularity bind path should be set using: `runOptions` 
      https://docs.sylabs.io/guides/3.0/user-guide/bind_paths_and_mounts.html. Example: `runOptions = "-B /mnt/"`

      
## How to Run: CASP MHC/KIR

1. ONT data should be basecalled with Dorado (traditional + 5mC) prior to running the pipeline.  Example basecalling command: `dorado basecaller sup,5mCG_5hmCG pod5/ --min-qscore 10 > SAMPLE_SFE.bam`
1. The CASP assembly pipeline supports SFE or a combination of SFE + UL reads.  Copy or link all basecalled ONT reads in bam format
   into the `sfe_bam` and/or `ul_bam` directory depending on what was sequenced, one file per sample. SFE bam files should be named: 
   `SAMPLE_SFE.bam`, where `SAMPLE` is a name of your choice.  UL fastq files should be named `SAMPLE_UL.bam`.
1. The project name can be modified by updating the `project_name` parameter in the config file (e.g.`casp_mhc.config`).  The `run_type` 
   parameter should be set to `sfe` or `sfe_ul` depending on the available reads.
1. The pipeline can be started with the command `nextflow run casp.nf -c casp_mhc.config`
   

## How to Run: HIFIASM MHC/KIR

1. ONT data should be basecalled with Dorado (traditional) prior to running the pipeline.  Example basecalling command: `dorado basecaller sup pod5/ --min-qscore 10 > SAMPLE_SFE.bam`
1. The hifiasm MHC/KIR assembly pipeline supports SFE-only, UL-only or a combination.  Copy the basecalled ONT reads in bam format
   into the `sfe_bam` and/or `ul_bam` depending on what was seqeunced. SFE bam files should be named: 
   `SAMPLE_SFE.bam`, where `SAMPLE` is a name of your choice.  UL fastq files should be named `SAMPLE_UL.bam`.
1. The project name can be modified by updating the `project_name` paramter in the config file `casp_hifiasm_[kir/mhc].config`. 
   The `run_type` parameter should be set to `sfe`, `ul` or `sfe_ul` depending on the available reads.
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


## Run Completion

1. On completion, output files can be found in the `project_name` directory.
    * `mhc_align`: alignment file(s) in bam format
    * `mhc_reads`: Raw and corrected MHC-specific reads in fastq format
    * `mhc_final_contigs`: MHC assemblies
    * `mhc_stats`: Pipeline statistics
1. Interrupted runs can be resumed: `nextflow run casp_mhc.nf -c casp_mhc.config -resume`
1. On workflow completion, intermediate files can be removed from the `work` directory


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

