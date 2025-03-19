#!/usr/bin/env python

"""
Split reads by haplotype
"""

from enum import Enum
from Bio import SeqIO
from utils import block_io as bio, paf_io as pio, probe_io as prio
from collections import defaultdict
import argparse
import sys
import pysam
import subprocess
import logging


class ReadAssign(Enum):
    HAP1 = 100,
    HAP2 = 200,
    HAP1_HAP = 300,
    HAP2_HAP = 400,
    HAP1_RAND = 500
    HAP2_RAND = 600
    UNTAGGED = 700,
    UNKNOWN = 800
    HOMOZYGOUS = 900


stl = []


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('contig_fasta', help="Path to collpased contigs in fasta format.")
    parser.add_argument('tagged_bam', help="Path to the haplotagged bam file.")
    parser.add_argument('ref_vcf', help="Variant calls across reference in VCF format.")
    parser.add_argument('contig_vcf', help="Variant calls across contigs in VCF format.")
    parser.add_argument('phase_blocks_gtf', help="Phased blocks across contigs in GTF format.")
    parser.add_argument('front_probe_bam', help="Path to front probe alignment file, bam format.")
    parser.add_argument('back_probe_bam', help="Path to back probe alignment file, bam format.")
    parser.add_argument('mhc_fasta', help="Path to the MHC reads in fasta format.")
    parser.add_argument('output_prefix', help="Prefix of the output files.")
    parser.add_argument('output_stats', help="Path to partition statistics output file.")
    parser.add_argument('-m', '--min_gt_count', type=int, default=1,
                        help="Minimum intersecting microarray positions required to assign phased block")
    parser.add_argument('-b', '--min_frac_contam_split', type=float, default=0.05,
                        help="Minimum contamination level that triggers a block split")
    parser.add_argument('-c', '--max_frac_contam_assign', type=float, default=0.15,
                        help="Maximum allowed minor haplotype fraction to assign block")
    parser.add_argument('-d', '--min_hemi_depth', type=int, default=16,
                        help="Minimum depth required to call a haploid position")
    parser.add_argument("-g", "--max_phase_distance", type=int, default=100000,
                        help="The maximum distance between two hets that does not trigger hemizygous evaluation")
    parser.add_argument("-x", '--max_probe_mm', default=1, type=int,
                        help="The maximum number of allowed mismatches a probe alignment")
    parser.add_argument('-p', '--max_probe_placements', default=8, type=int,
                        help="The maximum number of alignments allowed for probe flanking sequences")
    parser.add_argument('-a', '--min_hemi_hets', default=3, type=int,
                        help="Minimum number of heterozygous positions within an unphased block to assign hemizygous")
    parser.add_argument("-n", "--minimap_path", help="Path to minimap2 binary", default="minimap2")
    parser.add_argument("-r", "--retain_split_corr", help="Retain split corrected reads", default=False, action="store_true")
    args = parser.parse_args()

    logging.basicConfig(format='%(levelname)s:%(asctime)s %(message)s', level=logging.INFO)

    logging.info("Loading contig info")
    contig_data = SeqIO.to_dict(SeqIO.parse(args.contig_fasta, 'fasta'))

    logging.info("Identifying hemizyous regions by overlap")
    hemi_dict = identify_hemi_regions(args.contig_fasta, args.output_prefix, args.minimap_path)

    logging.info("Loading phased block coordinates")
    phased_blocks = load_phased_blocks(args.phase_blocks_gtf)

    logging.info("Loading unphased block coordinates")
    unphased_blocks = generate_unphased_blocks(phased_blocks, contig_data)

    logging.info("Loading reference het information")
    probe_info_dict = prio.load_probe_data(args.ref_vcf)

    logging.info("Loading microarray probe alignment data")
    front_align_dict = prio.load_probe_alignment(args.front_probe_bam, probe_info_dict, True, args.max_probe_mm)
    back_align_dict = prio.load_probe_alignment(args.back_probe_bam, probe_info_dict, False, args.max_probe_mm)

    logging.info("Identifying valid alignments and determine contig haplotype at probe")
    valid_probe_dict = prio.identify_valid_probe_locations(front_align_dict, back_align_dict, contig_data,
                                                           args.max_probe_placements)

    logging.info("Intersect VCF with probe locations in phased blocks")
    add_probes_to_phased_blocks(args.contig_vcf, phased_blocks, valid_probe_dict)

    logging.info("Check for large gaps in phasing")
    check_long_phasing(phased_blocks, unphased_blocks, args.max_phase_distance)

    logging.info("Intersect VCF with probe locations in unphased blocks")
    add_probes_to_unphased_blocks(args.tagged_bam, unphased_blocks, valid_probe_dict, args.min_hemi_depth)

    logging.info("Determine if unphased blocks are hemizygous or homozygous")
    hom_blocks, hemi_blocks = assign_unphased_blocks(unphased_blocks, hemi_dict, args.min_frac_contam_split,
                                                     args.min_hemi_hets)

    logging.info("Check for inconsistent haplotype assignment in phased blocks")
    phased_blocks = split_phased_blocks(phased_blocks, hom_blocks, args.min_gt_count, args.min_frac_contam_split)

    logging.info("Classify reads in phased blocks")
    read_assignments = assign_reads_in_phased_blocks(args.tagged_bam, phased_blocks, args.min_gt_count,
                                                      args.max_frac_contam_assign)

    logging.info("Classify reads in hemizygous blocks")
    read_assignments = assign_reads_in_hemi_blocks(args.tagged_bam, hemi_blocks, phased_blocks, read_assignments,
                                                   args.max_frac_contam_assign)

    logging.info("Classify reads in homozygous blocks")
    read_assignments = assign_reads_in_homo_blocks(args.tagged_bam, hom_blocks, read_assignments)

    logging.info("Writing binned reads")
    write_split_fasta(args.tagged_bam, read_assignments, args.output_prefix, args.mhc_fasta, args.retain_split_corr)

    with open(args.output_stats, "w") as op_stats:
        op_stats.write("\n".join(stl) + "\n")


def identify_hemi_regions(contig_file, output_prefix, minimap2_path):
    """
    Contigs are aligned to each other using minimap.  Overlapping regions longer that 10KB
    are considered overlapping and stored in a dictionary.

    :param contig_file: Path to assembly in fasta format
    :param output_prefix: Prefix of the output file
    :param minimap2_path: Path to the minimap2 binary.
    :return:
    """

    stl.append("Hemizygous Regions\nTigName\tStart\tEnd\tLength\tPair")

    self_align = output_prefix + "_self.paf"
    cproc = subprocess.run([minimap2_path, '-x', 'asm20', "-D", "--secondary=no", '-o',
                            self_align, contig_file, contig_file])

    if cproc.returncode != 0:
        logging.error('Error running self-alignment')
        sys.exit(1)

    self_paf = pio.load_query_alignments(self_align)
    hemi_regions = {}
    idx = 1
    for align_list in self_paf:
        for a in align_list:
            if a.query_length > a.target_length or a.align_length < 10000 or a.query_name == a.target_name:
                continue

            if a.target_name not in hemi_regions:
                hemi_regions[a.target_name] = []
            hemi_regions[a.target_name].append([a.target_start, a.target_end])

            if a.query_name not in hemi_regions:
                hemi_regions[a.query_name] = []
            hemi_regions[a.query_name].append([a.query_start, a.query_end])

            stl.append("{0}\t{1}\t{2}\t{3}\t{4}".format(a.query_name, a.query_start, a.query_end, a.align_length, idx))
            stl.append("{0}\t{1}\t{2}\t{3}\t{4}".format(a.target_name, a.target_start, a.target_end, a.align_length, idx))
            idx += 1

    return hemi_regions


def load_phased_blocks(block_file):
    """
    Load phased block GTF file and return a list of phased block objects

    :param block_file: Path to the phased block GTF
    :return: List of phased block objects.
    """

    phase_blocks = {}
    with open(block_file, 'r') as ipf:
        for line in ipf:
            items = line.strip().split('\t')

            block_name = items[0] + ":" + items[3] + ":" + items[4]
            phase_set_name  = items[0] + ":" + items[3]
            phase_blocks[phase_set_name] = bio.PhaseBlock(block_name, items[0], int(items[3]), int(items[4]) + 1)

    logging.info('{0} phased blocks found'.format(len(phase_blocks.keys())))
    return phase_blocks


def generate_unphased_blocks(phased_blocks, contig_file):
    """
    Set regions between phased blocks as unphased blocks.  If there is a single het between two longer phased
    blocks, set the whole region between the longer phased blocks as unphased.  If there is more than one
    het between longer phased blocks, don't make into a longer region.

    :param phased_blocks: List of phased block objects across the assembly
    :param contig_file: Dictionary of contig sequences
    :return: List of unphased block objects
    """

    unphased_blocks = {}
    for contig in contig_file.keys():
        start = 0
        for curr_block in sorted([b for b in phased_blocks.values() if b.contig == contig], key=lambda b: b.start):
            create_unphased_block(contig, start, curr_block.start, unphased_blocks)
            start = curr_block.end

        create_unphased_block(contig, start, len(contig_file[contig]), unphased_blocks)

    logging.info('{0} unphased blocks found'.format(len(unphased_blocks.keys())))
    return unphased_blocks


def create_unphased_block(contig, start, end, unphased_blocks):
    name = "{0}:{1}-{2}".format(contig, start, end)
    unphased_blocks[name] = bio.PhaseBlock(name, contig, start, end)


def check_long_phasing(phased_blocks, unphased_blocks, max_phase_gap):
    """
    Long distances in phasing could be because of hemizygous regions vs homozygous regions.  If the distance
    between two hets is greater than max_phase_gap, evaluate gap as a potential hemizygous region by adding to
    the list of unphased blocks.

    :param phased_blocks: List of phased blocks
    :param unphased_blocks: List of unphased blocks
    :param max_phase_gap: Maximum allowed gap between hets before checking gap
    """

    for block in phased_blocks.values():
        last = None
        for snp in sorted(block.snp_list, key=lambda x: x.contig_pos):
            if last is not None and snp.contig_pos - last > max_phase_gap:
                logging.info("Spiking in hemizgous region across large phased gap: {0}:{1}-{2}:{3}".format(block.contig, last,
                                                                                                    snp.contig_pos,
                                                                                                    snp.contig_pos - last))
                name = "{0}:{1}-{2}".format(block.contig, last, snp.contig_pos)
                unphased_blocks[name] = bio.PhaseBlock(name, block.contig, last, snp.contig_pos)
            last = snp.contig_pos


def add_probes_to_phased_blocks(phased_vcf, phase_blocks, valid_probe_data):
    """
    Read through assembly hets and add to the appropriate phase blocks.  If a reference het is found
    in the same location with the expected bases, assign a haplotype to the position. Later all haplotypes
    assignments across the block will be checked.

    Phased block data is updated, nothing is returned.

    :param phased_vcf: Phased assembly hets in VCF format
    :param phase_blocks: List of phase blocks across the assembly
    :param valid_probe_data: Dictionary of reference hets locations across the assembly
    """

    phased_cnt = usable_cnt = hap1_cnt = hap2_cnt = unk_cnt = 0

    with pysam.VariantFile(phased_vcf) as vcf_in:
        sn = prio.get_sample_name(vcf_in)
        for v in vcf_in:
            loc = '{0}:{1}'.format(v.contig, v.pos)

            # Extract genotype information
            sample_data = v.samples[sn]
            sep = '|' if sample_data.phased else '/'
            gt = sep.join([str(g) for g in sample_data['GT']])
            phase_set = '{0}:{1}'.format(v.contig, sample_data['PS'])

            # Count the number of phased assembly SNPs are in the block and then intersect with ref.
            if gt == '0|1' or gt == '1|0':
                phased_cnt += 1
                phase_blocks[phase_set].add_vcf_pos(v.pos)
                phase_blocks[phase_set].add_self_phase_set(0) if gt == "0|1" else phase_blocks[phase_set].add_self_phase_set(1)

                # only process positions found in reference VCF
                if loc in valid_probe_data:
                    probe = valid_probe_data[loc]

                    usable_cnt += 1

                    ref, alt = (v.ref, v.alts[0]) if gt == '0|1' else (v.alts[0], v.ref)

                    assign = None
                    if ref == probe.get_base(True) and alt == probe.get_base(False):
                        assign = prio.ProbeAssign.HAP1
                        hap1_cnt += 1
                    elif ref == probe.get_base(False) and alt == probe.get_base(True):
                        assign = prio.ProbeAssign.HAP2
                        hap2_cnt += 1
                    else:
                        unk_cnt += 1

                    if assign is not None:
                        probe.set_haplotype(assign)
                        phase_blocks[phase_set].add_snp(probe)

    stl.append("\nSNP Assignment Results")
    stl.append("ref_het_placements\t{0}".format(len(valid_probe_data.keys())))
    stl.append("contig_het_calls\t{0}".format(phased_cnt))
    stl.append("ovl_ref-contig_het\t{0}".format(usable_cnt))
    stl.append("ovl_ref-contig_het_hap1\t{0}".format(hap1_cnt))
    stl.append("ovl_ref-contig_het_hap2\t{0}".format(hap2_cnt))
    stl.append("ovl_ref-contig_het_unk\t{0}".format(unk_cnt))


def add_probes_to_unphased_blocks(tagged_bam, unphased_blocks, valid_probe_dict, min_depth):
    """
    Identify reference het positions that fall within unphased blocks, likely in hemizygous regions.
    Check the coverage depth for these positions before assignment, SNPs are only valid if the
    coverage is high enough for variant calling.

    Unphased block data is updated, nothing is returned.

    :param tagged_bam: Alignment file in bam format, used for calculating SNP depth.
    :param unphased_blocks: List of unphased blocks.
    :param valid_probe_dict: Dictionary of reference hets locations across the assembly
    :param min_depth:
    """

    last_contig = None
    pileup = None
    used = 0

    with pysam.AlignmentFile(tagged_bam, 'rb') as bam_in:
        for block in sorted(unphased_blocks.values(), key=lambda x: (x.contig, x.start)):
            if block.contig != last_contig:
                pileup = bam_in.pileup(block.contig)
                last_contig = block.contig

            for probe in valid_probe_dict.values():
                if probe.contig_name == block.contig and block.start <= probe.contig_pos < block.end:
                    cov = 0

                    # Check coverage
                    for pc in pileup:
                        if pc.pos == probe.contig_pos:
                            cov = pc.n
                        elif pc.pos > probe.contig_pos:
                            break

                    if cov < min_depth:
                        probe.set_invalid()
                    block.add_snp(probe)
                    used += 1

    stl.append("ref_het_unphased_block\t{0}".format(used))


def identify_homozyous_regions(block, hemizygous_blocks, homozygous_blocks, pred_hemi_regions, min_het_assign, min_edge_length=5000):
    """
    If an unphased block contains no haplotype assignments, the block is assumed to be homozygous. If there
    are haplotype assignments, the ends of the block are checked to see if there are more than min_edge_length
    bases without any haplotype assignments. If there are, they are checked for overlap with the predicted
    hemizygous regions of the assembly. If there is no overlap, the end is considered homozygous and split off
    from the rest of the block.  Short ends or ends that completely fall within predited hemizygous regions
    are not split off.

    :param block: unphased block object
    :param hemizygous_blocks: Dictionary of hemizygous block objects
    :param homozygous_blocks: Dictionary of homozygous blocks objects
    :param pred_hemi_regions: Dictionary of homozygous and hemizygous assembly blocks
    :param min_het_assign: Minimum number of heterozygous positions for hemizygous block assignment
    """

    split_blocks = block.split_hom_block(pred_hemi_regions, min_het_assign, min_edge_length)

    if split_blocks[0] is None:
        homozygous_blocks[block.name] = block
    else:
        hemizygous_blocks[split_blocks[0].name] = split_blocks[0]
        for b in split_blocks[1]:
            homozygous_blocks[b.name] = b


def assign_unphased_blocks(unphased_blocks, pred_hemi_regions, max_contam, min_het_assign, min_het_split=10):
    """
    For each unphased block in the assembly, predict if the region is hemizygous or homozygous. Prior to assignment,
    check unphased blocks to see if there is a mix of haplotype assignments, defined as at least min_hets
    and more than max_contam of the less frequent haplotype. If contaminated, the block is split once to
    minimize the contamination. If there is no good place to split, the block is unchanged.

    Example1: min_hets:5 min_contam:0.2, block has haplotypes 1111121111, contam too low for splitting (0.10)
    Example2: min_hets:5 min_contam:0.2, block has haplotypes 1211221221, contam high, but no good split
    Example3: min hets:5 min_contam:0.2, block has haplotypes 1111122222, contam high returns  11111, 22222

    After the level of contamination is evaluated, the blocks are then split into hemizygous and homozygous
    regions based on reference het placement and predicted hemizygous regions.

    :param unphased_blocks: Dictionary of unphased block objects
    :param pred_hemi_regions: Dictionary of overlapping contig regions, which are predicted to be hemizygous
    :param max_contam: Minimum amount of contamination needed to trigger split.
    :param min_het_assign: Minimum number of heterozygous positions in a unphased block to assign hemizygous
    :param min_het_split: Minimum number of heterozygous positions in a unphased block to check for contamination
    :return: Dictionary of homozygous and hemizygous assembly blocks
    """

    homozygous_blocks = {}
    hemizygous_blocks = {}
    for block in unphased_blocks.values():
        split_blocks = block.split_hemi_block(max_contam, min_het_split)
        if len(split_blocks) == 1:
            identify_homozyous_regions(split_blocks[0], hemizygous_blocks, homozygous_blocks, pred_hemi_regions, min_het_assign)
        else:
            identify_homozyous_regions(split_blocks[0], hemizygous_blocks, homozygous_blocks, pred_hemi_regions, min_het_assign)
            identify_homozyous_regions(split_blocks[1], hemizygous_blocks, homozygous_blocks, pred_hemi_regions, min_het_assign)
            homozygous_blocks[split_blocks[2].name] = split_blocks[2]
    return homozygous_blocks, hemizygous_blocks


def split_phased_blocks(phased_blocks, hom_blocks, min_count, max_contam_frac):
    """
    If the minor haplotype is at greater than max_contam_frac in a phased bock, the block is split up into
    sub-blocks and the regions between are set to homozygous

    :param phased_blocks: List of phased block objects
    :param hom blocks: Dictionary of homozygous block objects
    :param min_count: Minimum number of hets required for haplotype assignment
    :param max_contam_frac: Maximum amount of contaimination allowed before splittig
    :return: List of phase blocks after splitting.
    """

    updated_blocks = {}
    for name, block in phased_blocks.items():
        result = block.classify_block(min_count, max_contam_frac)
        if result == bio.BlockAssign.MIXED:
            logging.info("Splitting block: {0}".format(block.name))
            logging.info(block.hap_string)
            split_blocks = block.split_phased_block(10)
            if split_blocks is None:
                updated_blocks[name] = block
            else:
                logging.info("Generated {0} blocks".format(len(split_blocks[0])))
                for i, block in enumerate(split_blocks[0]):
                    updated_name = "{0}:{1}".format(name, i)
                    updated_blocks[updated_name] = block
                for hb in split_blocks[1]:
                    hom_blocks[hb.name] = hb
        else:
            updated_blocks[name] = block
    return updated_blocks


def assign_reads_in_phased_blocks(tagged_bam, phased_blocks, min_count, max_frac_contam):
    """
    Check the haplotype assigment for each phased block. If a haplotype can be assigned, iterate over the haplotagged
    reads assign each a haplotype based on the block hapotype and the tag. Reads that are not tagged are classified
    as unknown and reads in blocks that have ambiguous or unknown haplotype status are also classifed as unknown

    :param tagged_bam: Haplotagged read alignments to collapsed assembly
    :param phased_blocks: List of phased blocks objects across the assembly
    :param min_count: Minimum number of SNPs required for haplotype assignment
    :param max_frac_contam: Maximum allowed contamination for haplotype assignment
    :return: Dictionary of haplotype assignments for each read
    """

    stl.append("\nPhased Block Assignment\nBlockName\tBlockLength\tBlockStatus\tBlockReads\tTaggedReads\tUntagReads\t"
               "UnkReads\tTigHet\tRefHet\tRefHetHap1\tRefHetHap2\tNoise")

    # For each phase block, grab member reads
    read_assign = {}
    snp_class_counts = defaultdict(int)

    with pysam.AlignmentFile(tagged_bam, 'rb') as bam_in:
        for block in sorted(phased_blocks.values(), key=lambda x: (x.contig, x.start)):
            block_read_counts = defaultdict(int)
            hap = block.classify_block(min_count, max_frac_contam)
            phased_snp_count = len(block.vcf_pos)
            call_string = "NA"

            # Count number of assigned VCF variants
            if hap in [bio.BlockAssign.HAP1, bio.BlockAssign.HAP2]:
                snp_class_counts['PHASED'] += phased_snp_count
            elif len(block.vcf_pos) == 1:
                snp_class_counts['ISO'] += 1
            else:
                snp_class_counts['UNKNOWN'] += phased_snp_count

            for read in bam_in.fetch(block.contig, block.start, block.end):
                try:
                    # Skip if read has already been assigned or if read isn't primary
                    if read.query_name in read_assign or read.is_supplementary or read.is_secondary:
                        continue

                    # Only process reads that are tagged with whatshap
                    tag = read.get_tag('HP')
                    ps = read.get_tag("PS")
                    if hap == bio.BlockAssign.HAP1:
                        block_read_counts["ASSIGN"] += 1
                        read_assign[read.query_name] = ReadAssign.HAP1 if tag == 1 else ReadAssign.HAP2
                    elif hap == bio.BlockAssign.HAP2:
                        block_read_counts["ASSIGN"] += 1
                        read_assign[read.query_name] = ReadAssign.HAP2 if tag == 1 else ReadAssign.HAP1
                    elif hap == bio.BlockAssign.MIXED:
                        read_assign[read.query_name] = ReadAssign.UNKNOWN
                        call_string = block.hap_string
                        block_read_counts["UNKNOWN"] += 1
                    else:
                        # Set the haplotype based on the haplotag.
                        read_assign[read.query_name] = ReadAssign.HAP1_RAND if tag == 1 else ReadAssign.HAP2_RAND
                        block_read_counts["UNKNOWN"] += 1

                except KeyError: #Read is not haplotagged
                    block_read_counts["UNTAGGED"] += 1
                    read_assign[read.query_name] = ReadAssign.UNTAGGED

            stl.append("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}".format(block.name, block.length,
                        hap.name, sum(block_read_counts.values()), block_read_counts["ASSIGN"],
                        block_read_counts["UNTAGGED"],  block_read_counts["UNKNOWN"], phased_snp_count,
                        block.hap1_cnt + block.hap2_cnt, block.hap1_cnt, block.hap2_cnt, call_string))

    stl.append("\nPhased Block Summary")
    stl.append("phased_block_reads\t{0}".format(len(read_assign.values())))
    stl.append("phased_block_reads_hap1\t{0}".format(len([x for x in read_assign.values() if x == ReadAssign.HAP1])))
    stl.append("phased_block_reads_hap2\t{0}".format(len([x for x in read_assign.values() if x == ReadAssign.HAP2])))
    stl.append("phased_block_reads_untagged\t{0}".format(len([x for x in read_assign.values() if x == ReadAssign.UNTAGGED])))
    stl.append("phased_block_reads_unknown_mixed\t{0}".format(len([x for x in read_assign.values() if x == ReadAssign.UNKNOWN])))
    stl.append("phased_block_reads_unknown_unanchored\t{0}".format(len([x for x in read_assign.values()
                                                                     if x == ReadAssign.HAP1_RAND or x == ReadAssign.HAP2_RAND])))
    stl.append("phased_block_snps\t{0}".format(sum(snp_class_counts.values())))
    stl.append("assigned_block_snps\t{0}".format(snp_class_counts['PHASED']))
    stl.append("unknown_block_snps\t{0}".format(snp_class_counts['UNKNOWN']))
    stl.append("single_het_block_snps\t{0}".format(snp_class_counts['ISO']))

    return read_assign


def assign_reads_in_hemi_blocks(tagged_bam, unphased_blocks, phased_blocks, read_assign, max_frac_contam):
    """
    Check the haplotype assigment of each hemizygous block. If haplotype can be assigned, all reads in the block that
    are unknown or unassessed are assigned to the haplotype.

    :param tagged_bam: Haplotagged read alignments to collapsed assembly
    :param unphased_blocks: Dictionary containing hemizygous block objects
    :param read_assign: Dictionary of read assignments
    :param max_frac_contam: Max allowed contamination for haplotype assignment
    :return: Updated read assignments
    """

    stl.append("\nHemi Block Assignment\nBlockName\tBlockLength\tBlockStatus\tBlockReads\tRefHet\tRefHetHap1\tRefHetHap2\tNoise")

    read_class_counts = defaultdict(int)
    reassign = set()
    with pysam.AlignmentFile(tagged_bam, 'rb') as bam_in:
        for block in sorted(unphased_blocks.values(), key=lambda x: (x.contig, x.start)):
            hap = block.classify_block(2, max_frac_contam)
            call_string = "NA"

            block_reads = 0
            for read in bam_in.fetch(block.contig, block.start, block.end):
                # Skip over already assigned reads
                assign = None if read.query_name not in read_assign else read_assign[read.query_name]
                if assign == ReadAssign.HAP1 or assign == ReadAssign.HAP2 or read.is_supplementary or read.is_secondary:
                    continue

                # If a read was in an unassigned phase block and a haplotig region, try and recover phasing
                tag = None
                ref_ps = None
                if read.has_tag("HP"):
                    tag = read.get_tag("HP")
                    name_ps = "{0}:{1}".format(block.contig, read.get_tag("PS"))
                    if name_ps in phased_blocks:
                        pb = phased_blocks[name_ps]
                        if pb.phase_set_self[0] > pb.phase_set_self[1]:
                            ref_ps = 1
                        elif pb.phase_set_self[0] < pb.phase_set_self[1]:
                            ref_ps = 2

                block_reads += 1
                if hap == bio.BlockAssign.HAP1:
                    if assign == ReadAssign.HAP2_HAP:
                        read_class_counts['MIX'] += 1
                        read_assign[read.query_name] = ReadAssign.UNKNOWN
                    elif (assign == ReadAssign.HAP1_RAND or assign == ReadAssign.HAP2_RAND) and ref_ps is not None and read.query_name not in reassign:
                        # Update the read assignment based on how read matches to contig
                        if (ref_ps == 1 and tag == 1) or (ref_ps == 2 and tag == 2):
                            logging.info("Updating unanchored phased read {0} to haplotype 1 in region {1}".format(read.query_name, name_ps))
                            read_assign[read.query_name] = ReadAssign.HAP1_RAND
                            read_class_counts['TAG_HAP1'] += 1
                        else:
                            logging.info("Updating unanchored phased read {0} to haplotype 2 in region {1}".format(read.query_name, name_ps))
                            read_assign[read.query_name] = ReadAssign.HAP2_RAND
                            read_class_counts['TAG_HAP2'] += 1
                        reassign.add(read.query_name)
                    else:
                        read_class_counts['HAP1'] += 1
                        read_assign[read.query_name] = ReadAssign.HAP1_HAP
                elif hap == bio.BlockAssign.HAP2:
                    if assign == ReadAssign.HAP1_HAP:
                        read_class_counts['MIX'] += 1
                        read_assign[read.query_name] = ReadAssign.UNKNOWN
                    elif (assign == ReadAssign.HAP1_RAND or assign == ReadAssign.HAP2_RAND) and ref_ps is not None and read.query_name not in reassign:
                        # Update the read assignment based on how read matches to contig
                        if (ref_ps == 1 and tag == 1) or (ref_ps == 2 and tag == 2):
                            logging.info("Updating unanchored phased read {0} to haplotype 2 in region {1}".format(read.query_name, name_ps))
                            read_assign[read.query_name] = ReadAssign.HAP2_RAND
                            read_class_counts['TAG_HAP2'] += 1
                        else:
                            logging.info("Updating unanchored phased read {0} to haplotype 1 in region {1}".format(read.query_name, name_ps))
                            read_assign[read.query_name] = ReadAssign.HAP1_RAND
                            read_class_counts['TAG_HAP1'] += 1
                        reassign.add(read.query_name)
                    else:
                        read_class_counts['HAP2'] += 1
                        read_assign[read.query_name] = ReadAssign.HAP2_HAP
                else:
                    read_class_counts['UNKNOWN'] += 1
                    read_assign[read.query_name] = ReadAssign.UNKNOWN
                    call_string = block.hap_string

            stl.append("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}".format(block.name, block.length, hap.name,  block_reads,
                                                                       block.hap1_cnt + block.hap2_cnt, block.hap1_cnt,
                                                                       block.hap2_cnt, call_string))

    stl.append("\nHemi Block Summary")
    stl.append("hemi_block_reads\t{0}".format(sum(read_class_counts.values())))
    stl.append("hemi_block_reads_hap1\t{0}".format(read_class_counts['HAP1']))
    stl.append("hemi_block_reads_hap2\t{0}".format(read_class_counts['HAP2']))
    stl.append("hemi_block_reads_hap1_rand\t{0}".format(read_class_counts['TAG_HAP1']))
    stl.append("hemi_block_reads_hap2_rand\t{0}".format(read_class_counts['TAG_HAP2']))
    stl.append("hemi_block_reads_unknown\t{0}".format(read_class_counts['UNKNOWN']))
    stl.append("hemi_block_reads_mixed\t{0}".format(read_class_counts['MIX']))

    return read_assign


def assign_reads_in_homo_blocks(tagged_bam, homo_blocks, read_assignments):
    """
    Check reads that align to homozygous blocks.  If reads have been not prevously seen, set to unknown

    :param tagged_bam:  Haplotagged read alignments to collapsed assembly
    :param homo_blocks: Dictionary of homozygous blocks across the assembly
    :param read_assignments: Dictionary of read assignments
    :return: Updated read assignments
    """

    stl.append("\nHom Block Assignments\nBlockName\tBlockLength\tBlockReads\tUnkReads\tAssignReads")
    reads_in_homozygous_blocks = 0
    with pysam.AlignmentFile(tagged_bam, 'rb') as bam_in:
        for block in sorted(homo_blocks.values(), key=lambda x: (x.contig, x.start)):
            known_count = unknown_count = 0
            for read in bam_in.fetch(block.contig, block.start, block.end):
                if read.is_supplementary or read.is_secondary:
                    continue
                if read.query_name not in read_assignments:
                    read_assignments[read.query_name] = ReadAssign.HOMOZYGOUS
                    unknown_count += 1
                else:
                    known_count += 1
            stl.append("{0}\t{1}\t{2}\t{3}\t{4}".format(block.name, block.length, unknown_count+known_count,
                                                       unknown_count, known_count))
            reads_in_homozygous_blocks += unknown_count

    stl.append("\nHom Block Summary")
    stl.append("hom_block_reads\t{0}".format(reads_in_homozygous_blocks))
    return read_assignments


def write_split_fasta(tagged_bam, collapsed_read_assignments, output_prefix, mhc_fasta, retain_split_corr):
    """
    Write out binned reads based on the assignments

    :param tagged_bam: Haplotagged read alignments to collapsed assembly
    :param collapsed_read_assignments: Read assignment dictionary
    :param output_prefix: Prefix the output files
    :param mhc_fasta: Path to MHC-specific reads in fasta format
    :param retain_split_reads: If false, do not write out split corrected reads
    """

    read_counts = defaultdict(int)
    fasta_dict = defaultdict(list)
    if retain_split_corr:
        with open(mhc_fasta, "r") as ipf:
            for record in SeqIO.parse(ipf, "fasta"):
                orig_name = record.id.split("_")[0].split(":")[0]
                fasta_dict[orig_name].append(record)
    else:
        with open(mhc_fasta, "r") as ipf:
            for record in SeqIO.parse(ipf, "fasta"):
                fasta_dict[record.id].append(record)

    with pysam.AlignmentFile(tagged_bam, 'rb') as bam_in, open(output_prefix + '_hap1_only.fasta', 'w') as hap1_out, \
            open(output_prefix + '_hap2_only.fasta', 'w') as hap2_out, open(output_prefix + '_unknown.fasta', 'w') as unknown_out:
        for a in bam_in:
            if a.is_unmapped or a.is_secondary or a.is_supplementary:
                continue

            if a.query_name in collapsed_read_assignments:
                call = collapsed_read_assignments[a.query_name]
                if call == ReadAssign.HAP1:
                    handle = hap1_out
                    read_counts["hap1"] += 1
                elif call == ReadAssign.HAP1_HAP:
                    handle = hap1_out
                    read_counts["hap1_hap"] += 1
                elif call == ReadAssign.HAP1_RAND:
                    handle = hap1_out
                    read_counts["hap1_rand"] += 1
                elif call == ReadAssign.HAP2:
                    handle = hap2_out
                    read_counts["hap2"] += 1
                elif call == ReadAssign.HAP2_HAP:
                    handle = hap2_out
                    read_counts["hap2_hap"] += 1
                elif call == ReadAssign.HAP2_RAND:
                    handle = hap2_out
                    read_counts["hap2_rand"] += 1
                elif call == ReadAssign.UNKNOWN:
                    handle = unknown_out
                    read_counts["unknown"] += 1
                elif call == ReadAssign.UNTAGGED:
                    handle = unknown_out
                    read_counts["untagged"] += 1
                elif call == ReadAssign.HOMOZYGOUS:
                    handle = unknown_out
                    read_counts["homozygous"] += 1
                else:
                    logging.info("unexpected read class {0}".format(call.name))
                    logging.info(a.query_name)
                    sys.exit(1)
            else:
                handle = unknown_out
                read_counts["unobserved"] += 1

            read_name = a.query_name
            if read_name not in fasta_dict:
                continue
            SeqIO.write(fasta_dict[read_name], handle, "fasta")

        hap1_tot = read_counts["hap1"] + read_counts["hap1_hap"] + read_counts["hap1_rand"]
        hap2_tot = read_counts["hap2"] + read_counts["hap2_hap"] + read_counts["hap2_rand"]
        both_tot = read_counts["unknown"] + read_counts["untagged"] + read_counts["unobserved"] + read_counts["homozygous"]

        stl.append("\nFinal Read Binning Stats")
        stl.append("total_written_reads\t{0}\n".format(sum(read_counts.values())))
        stl.append("Assign\tAssignTot\tSub\tSubTot")
        stl.append("{0}\t{1}\t{2}\t{3}".format("hap1", hap1_tot, "phased", read_counts["hap1"]))
        stl.append("{0}\t{1}\t{2}\t{3}".format("hap1", hap1_tot, "hemi", read_counts["hap1_hap"]))
        stl.append("{0}\t{1}\t{2}\t{3}".format("hap1", hap1_tot, "random", read_counts["hap1_rand"]))
        stl.append("{0}\t{1}\t{2}\t{3}".format("hap2", hap2_tot, "phased", read_counts["hap2"]))
        stl.append("{0}\t{1}\t{2}\t{3}".format("hap2", hap2_tot, "hemi", read_counts["hap2_hap"]))
        stl.append("{0}\t{1}\t{2}\t{3}".format("hap2", hap1_tot, "random", read_counts["hap2_rand"]))
        stl.append("{0}\t{1}\t{2}\t{3}".format("both", both_tot, "homozygous", read_counts["homozygous"]))
        stl.append("{0}\t{1}\t{2}\t{3}".format("both", both_tot, "untagged", read_counts["untagged"]))
        stl.append("{0}\t{1}\t{2}\t{3}".format("both", both_tot, "unknown", read_counts["unknown"]))


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(1)

