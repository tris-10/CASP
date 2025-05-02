#!/usr/bin/env python

"""
This scripts attempts to fill in phasing gaps left in the reference assembly.
"""

from collections import defaultdict
from Bio import SeqIO
from utils import probe_io
import pysam
import argparse
import sys
import logging


class BlockConnect:
    def __init__(self, connected_ps, orient):
        self.connected_ps = connected_ps
        self.orient = orient


class ContigHetClass:
    def __init__(self, vr, sn):
        sd = vr.samples[sn]
        self.contig = vr.contig
        self.pos = vr.pos
        self.loc = "{0}:{1}".format(self.contig, self.pos)
        self.phased = sd.phased
        self.ps = sd['PS'] if sd.phased else None
        self.b1 = vr.ref if sd['GT'][0] == 0 else vr.alts[0]
        self.b2 = vr.ref if sd['GT'][1] == 0 else vr.alts[0]


stl = []


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('contig_fasta', help="Path to contig fasta file.")
    parser.add_argument('contig_vcf', help="Path to whatshap phased contig VCF file.")
    parser.add_argument('front_probe_bam', help="Path to front probe alignment file.")
    parser.add_argument('back_probe_bam', help="Path to back probe alignment file.")
    parser.add_argument('ref_vcf', help="Path to whatshap phased reference VCF file")
    parser.add_argument("ref_gtf", help="Path to reference phase blocks in GTF format.")
    parser.add_argument('output_vcf', help="Path to updated phased reference VCF file output.")
    parser.add_argument('output_stats', help="Path to output statistics.")
    parser.add_argument("-n", '--max_probe_mm', default=1, type=int,
                       help="Maximum mismatches per flank when placing ref variants on the assembly.")
    parser.add_argument('-p', '--max_probe_placements', default=8, type=int,
                        help="Maximum number of flank alignments allowed when placing ref variants on the assembly.")
    g1 = parser.add_argument_group("Haploid Phasing")
    g1.add_argument("-m", '--min_het_pset', default=5, type=int, help="Minimum SNPs per phase set.")
    g1.add_argument("-f", '--min_frac_hap', default=0.75, type=float, help="Minimum haplotype purity")

    args = parser.parse_args()

    logging.basicConfig(format='%(levelname)s:%(asctime)s %(message)s', level=logging.INFO)

    logging.info("Loading contig data")
    contig_data = SeqIO.to_dict(SeqIO.parse(args.contig_fasta, 'fasta'))

    logging.info("Loading reference phase blocks")
    pset_order = load_ref_phase_blocks(args.ref_gtf)

    logging.info("Loading phased reference vcf")
    probe_info_dict = probe_io.load_probe_data(args.ref_vcf)
    logging.info("loaded {0} het SNPs".format(len(probe_info_dict.keys())))

    logging.info('Loading probe/het alignment data')
    front_align_dict = probe_io.load_probe_alignment(args.front_probe_bam, probe_info_dict, True, args.max_probe_mm)
    back_align_dict = probe_io.load_probe_alignment(args.back_probe_bam, probe_info_dict, False, args.max_probe_mm)

    logging.info('Identifying valid probe/het contig locations')
    valid_probe_by_loc = probe_io.identify_valid_probe_locations(front_align_dict, back_align_dict,
                                                                 contig_data, args.max_probe_placements)
    valid_probe_by_name = defaultdict(list)
    for probe in valid_probe_by_loc.values():
        valid_probe_by_name[probe.snp_id].append(probe)
    stl.append("Probe Placement")
    stl.append("valid_placements\t{0}".format(len(valid_probe_by_loc.keys())))
    stl.append("total_hets\t{0}".format(len(valid_probe_by_name.keys())))

    logging.info("Loading contig variant calls")
    contig_variants = load_contig_vcf_file(args.contig_vcf)

    logging.info("Phasing variants using overlapping haploid contigs")
    connect = phase_with_haploid_contigs(args.ref_vcf, valid_probe_by_name, contig_variants, pset_order,
                                         args.min_het_pset, args.min_frac_hap)

    logging.info("Phasing variants using diploid contigs")
    connect = phase_with_diploid_contigs(args.ref_vcf, valid_probe_by_name, contig_variants, connect)

    logging.info("Writing out updated phase data")
    update_phased_vcf(args.ref_vcf, args.output_vcf, connect)

    with open(args.output_stats, "w") as op_stats:
        op_stats.write("\n".join(stl) + "\n")


def phase_with_haploid_contigs(vcf_file, valid_probe_by_name, contig_variants, pset_order, min_var_count, min_frac):
    """
    1. Phased hets called against hg38 are saved if they align to exactly two contigs, were not called hets in the assembly
    and have different bases on each contig.
    2. Adjacent phase sets are checked for contigs that touch both sets with at least min_var_count hets (from step 1)
    within each phase set.
    3. If two contigs meet the criteria from step 2, the two contigs with the most step 1 hets are evaluated.
    4. The step 1 hets in the each phase set are checked to see if there are at least min_var_count found on both contigs
       that have different haplotype assignments.
    5. If step 4 passes, step 1 hets in each phase set are checked to make sure they are at least min_fraction of
       one haplotype, not a mix of haplotypes. If there isn't a mix, the phase sets are connected if they both have
       the same haplotype assignmets.  If the haplotype assignments are different, the phase sets are connected, but
       flipped.

    :param vcf_file: Path to hg38 variant calls in VCF format
    :param valid_probe_by_name: hg38 heterozygous variant locations on the assembly, stored by ID
    :param contig_variants: Heterozygous variants across the assembly
    :param pset_order: Ordered list of phase sets
    :param min_var_count: Minimum number of het observations per phase set for haploid phasing
    :param min_frac: Minimum haplotype assignment purity for haploid phasing
    :return: Phase set connections
    """

    stl.append("\nHaploid Phasing\nPhase_Sets\tStatus\tContigs\tSNP_Counts")
    connect = {}
    ps_contig = defaultdict(dict)

    # For every phased reference het that can be placed on two contigs but were not called hets in the assembly, store
    # phase-set, the contig name and the haplotype
    with pysam.VariantFile(vcf_file, 'r') as vcf_in:
        sn = probe_io.get_sample_name(vcf_in)
        for record in vcf_in:
            sd = record.samples[sn]

            # Skip record if the ref SNP is unphased, not found on the contig or not found twice
            if not sd.phased or record.id not in valid_probe_by_name or len(valid_probe_by_name[record.id]) != 2:
                continue

            pset = sd['PS']
            p1 = valid_probe_by_name[record.id][0]
            p2 = valid_probe_by_name[record.id][1]

            # Only use SNPs that have different calls on the different contigs
            # and are not called hets in the assembly
            if p1.haplotype == p2.haplotype or (p1.loc in contig_variants and p2.loc in contig_variants):
                continue

            for p in p1, p2:
                if pset not in ps_contig:
                    ps_contig[pset] = {}
                if p.contig_name not in ps_contig[pset]:
                    ps_contig[pset][p.contig_name] = [set(), set()]

                if p.haplotype == probe_io.ProbeAssign.HAP1:
                    ps_contig[pset][p.contig_name][0].add(record.id)
                elif p.haplotype == probe_io.ProbeAssign.HAP2:
                    ps_contig[pset][p.contig_name][1].add(record.id)

    # For each adjacent phase set, look for two contigs with the most valid expected hets and see if the
    # phase sets can be clearly linked via the contigs
    for i in range(1, len(pset_order)):
        ps1 = pset_order[i - 1]
        ps2 = pset_order[i]

        ps_label = "{0}-{1}".format(ps1, ps2)

        if ps1 not in ps_contig or ps2 not in ps_contig:
            stl.append("{0}\t{1}\t{2}\t{3}".format(ps_label, "low_ps_support", "NA", "NA"))
            continue

        # Get contig data for the paired phase sets
        ps1_contigs = ps_contig[ps1]
        ps2_contigs = ps_contig[ps2]

        # Find contigs that span both base sets. If there are at least min_var_count expected hets in each
        # phase, store the total number of supporting hets.  Rank contigs by het count
        potential_contigs_counts = {}
        for contig, sets1 in ps1_contigs.items():
            if contig in ps2_contigs:
                sets2 = ps2_contigs[contig]

                if len(sets1[0]) + len(sets1[1]) >= min_var_count and len(sets2[0]) + len(sets2[1]) >= min_var_count:
                    potential_contigs_counts[contig] = sum([len(sets1[0]), len(sets1[1]), len(sets2[0]), len(sets2[1])])

        potential_contigs = [x[0] for x in sorted(potential_contigs_counts.items(), key=lambda x: x[1], reverse=True)]

        # End processing if there aren't two contigs with sufficient supporting hets
        if len(potential_contigs) < 2:
            stl.append("{0}\t{1}\t{2}\t{3}".format(ps_label, "low_spanning_tigs", "NA", "NA"))
            continue

        c1_p1 = ps_contig[ps1][potential_contigs[0]]  # First contig upstream PS
        c1_p2 = ps_contig[ps2][potential_contigs[0]]  # First contig downstream PS
        c2_p1 = ps_contig[ps1][potential_contigs[1]]  # Second contig upstream PS
        c2_p2 = ps_contig[ps2][potential_contigs[1]]  # Second contig downstream PS

        # Count the number of PS SNPs are found on both contigs with different haplotypes
        p1_counts = [find_shared_snps(c1_p1[0], c2_p1[1]), find_shared_snps(c1_p1[1], c2_p1[0])]
        p2_counts = [find_shared_snps(c1_p2[0], c2_p2[1]), find_shared_snps(c1_p2[1], c2_p2[0])]

        if sum(p1_counts) < min_var_count or sum(p2_counts) < min_var_count:
            status = "low_shared_hets"
        else:
            p1_f = p1_counts[0] / sum(p1_counts)
            p2_f = p2_counts[0] / sum(p2_counts)

            if (p1_f >= min_frac and p2_f >= min_frac) or (p1_f <= (1-min_frac) and p2_f <= (1-min_frac)):
                status = "phased_inline"
                connect[ps2] = [ps1, 1]
            elif (p1_f >= min_frac and p2_f <= (1-min_frac)) or (p1_f <= (1-min_frac) and p2_f >= min_frac):
                status = "phased_flipped"
                connect[ps2] = [ps1, -1]
            else:
                status = "mixed_hets"

        counts_label = "{0}/{1}-{2}/{3}".format(p1_counts[0], p1_counts[1], p2_counts[0], p2_counts[1])
        contig_label = "{0}/{1}".format(potential_contigs[0], potential_contigs[1])
        stl.append("{0}\t{1}\t{2}\t{3}".format(ps_label, status, contig_label, counts_label))

    return connect


def phase_with_diploid_contigs(vcf_file, valid_probe_by_name, contig_variants, connect):
    """
    Scan the reference VCF file for SNPs that are unphased in the reference VCF, but are phased in the
    assembly. If this situation is found, check to make sure the bases match up between reference and contig
    and then see if the phase sets need to be inverted.

    :param vcf_file: Path to hg38 variant calls in VCF format
    :param valid_probe_by_name: hg38 heterozygous variant locations on the assembly, stored by ID
    :param contig_variants: Heterozygous variants across the assembly
    :param connect: Phase set connects
    :return: Updated phase set connections
    """

    stl.append("\nDiploid Phasing\nPhase_Sets\tStatus\tDistance\tContigs\tContig_Pos\tSNPs\tSNP1_CvR\tSNP2_CvR")
    with pysam.VariantFile(vcf_file, 'r') as vcf_in:
        last_pset = None
        last_pos = None
        last_id = None

        sn = probe_io.get_sample_name(vcf_in)

        for record in vcf_in:
            sd = record.samples[sn]
            if not sd.phased:
                continue
            pset = sd['PS']

            # move on if the SNP could not be placed on the contigs or if it isn't phased on the contig
            if record.id not in valid_probe_by_name or len(valid_probe_by_name[record.id]) > 1 or \
                    valid_probe_by_name[record.id][0].loc not in contig_variants or \
                    not contig_variants[valid_probe_by_name[record.id][0].loc].phased:
                continue

            # When the edge of phase set is reached, check to see if the border SNPs are phased in the assembly
            if last_pset is not None:
                distance = record.pos-last_pos

                rs1 = valid_probe_by_name[last_id][0]
                rs2 = valid_probe_by_name[record.id][0]
                cs1 = contig_variants[rs1.loc]
                cs2 = contig_variants[rs2.loc]
                rs_label = "{0}/{1}".format(rs1.snp_id, rs2.snp_id)

                pset_label = "{0}-{1}".format(last_pset, pset)
                contig_pos_label = "{0}/{1}".format(cs1.pos, cs2.pos)
                ref_pos_label = "{0}/{1}".format(rs1.probe_info.hg38_pos, rs2.probe_info.hg38_pos)

                # Stop processing if the SNPs are on different contigs or phase sets
                if cs1.contig == cs1.contig and cs1.ps == cs2.ps:
                    rs1b1 = rs1.get_base(True)
                    rs1b2 = rs1.get_base(False)
                    rs2b1 = rs2.get_base(True)
                    rs2b2 = rs2.get_base(False)

                    # Check to see if the reference/contig bases are consistent.  If they are, check to see if the
                    # phasing is inverted
                    status = "inconsistent_bases"
                    if rs1b1 == cs1.b1 and rs1b2 == cs1.b2:
                        if rs2b1 == cs2.b1 and rs2b2 == cs2.b2:
                            status = "phased_inline"
                        elif rs2b1 == cs2.b2 and rs2b2 == cs2.b1:
                            status = "phased_flipped"
                    elif rs1b1 == cs1.b2 and rs1b2 == cs1.b1:
                        if rs2b1 == cs2.b1 and rs2b2 == cs2.b2:
                            status = "phased_flipped"
                        elif rs2b1 == cs2.b2 and rs2b2 == cs2.b1:
                            status = "phased_inline"

                    snp1_label = "{0}/{1}-{2}|{3}".format(rs1b1, rs1b2, cs1.b1, cs1.b2)
                    snp2_label = "{0}/{1}-{2}|{3}".format(rs2b1, rs2b2, cs2.b1, cs2.b2)

                    if pset != last_pset:
                        if status == "phased_inline":
                            update_connections(connect, pset, last_pset, 1)
                        else:
                            update_connections(connect, pset, last_pset, -1)

                        stl.append("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}".format(pset_label, status, distance, cs1.contig,
                                                                          contig_pos_label, rs_label, snp1_label, snp2_label))
                    else:
                        if status == "phased_flipped":
                            logging.error("Found a potential error in the reference phasing {0} {1} "
                                          "{2}".format(rs_label, contig_pos_label, ref_pos_label))
                elif cs1.contig == cs1.contig:
                    if pset != last_pset:
                        stl.append("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}".format(pset_label, "contig_unphased", distance, cs1.contig,
                                                                     contig_pos_label, rs_label, "NA","NA"))
                else:
                    if pset != last_pset:
                        stl.append("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}".format(pset_label, "split_contig", distance, "NA",
                                                                     "NA","NA", "NA", "NA", "NA"))

            last_pset = pset
            last_pos = record.pos
            last_id = record.id
    return connect


def find_shared_snps(set1, set2):
    return len([x for x in set1 if x in set2])


def update_connections(connect, pset, last_pset, orient):
    if pset in connect and connect[pset][0] > last_pset:
        prev = connect[pset]
        connect[prev[0]] = [last_pset, prev[1] * orient]
    else:
        connect[pset] = [last_pset, orient]


def load_contig_vcf_file(vcf_file):
    """
    Store contig het information in a dictionary keyed by assembly location

    :param vcf_file: Path to contig VCF file
    :return: Dictionary of contig hets, keyed by location
    """

    contig_variant_dict = {}
    with pysam.VariantFile(vcf_file, 'r') as vcf_in:
        sn = probe_io.get_sample_name(vcf_in)
        for record in vcf_in:
            het = ContigHetClass(record, sn)
            contig_variant_dict[het.loc] = het
    return contig_variant_dict


def load_ref_phase_blocks(gtf_file):
    """
    Load phased block data across hg38

    :param gtf_file: Whatshap generated block data in GTF format
    :return: list of phase sets
    """

    ps_order = []
    with open(gtf_file, "r") as ipf:
        for line in ipf:
            ps_order.append(int(line.split("\t")[3]))
    return ps_order


def update_phased_vcf(ref_vcf, output_vcf, connect):
    """
    Update the phase sets and genotypes in the hg38 reference file.  The connection file contains the orientation
    of a phase set to the upstream phase set, stored as 1 or -1.  If multiple phase sets are now connected, you
    can walk through the connection dictionary to get the final orientation by multiplying the orientations as
    you go.

    For example if there are three PS: PS1, PS2 and PS3. If PS2 is inverted relative to PS1 and PS3 is
    not inverted relative to PS2, the data would be stored:

    PS2 = [PS1, -1]
    PS3 = [PS2, 1]

    1. When PS2 is reached, you can connect back to PS1 [ PS2 -> PS1 ] and since the orientation is -1 it would be
    flipped relative to PS1
    2. When PS3 is reached, you can connect back to PS1 [ PS3 -> PS2 -> PS1 ]  the orientation would be -1 * 1 = -1,
    so PS3 is also flipped relative to PS1

    If instead PS3 is inverted relative to PS2:
    PS2 = [PS1, -1]
    PS3 = [PS2, -1]

    When PS3 is reached, again you can connect back to PS1 [ PS3 -> PS2 -> PS1 ], but the orientation would be -1 * -1 = 1
    and PS3 is now unflipped relative to PS2.


    :param ref_vcf: Path to hg38 variant calls in VCF format
    :param output_vcf: Path to hg38 variant calls in VCF format with updated phasing
    :param connect: connection dictionary
    """

    stl.append("\nBlock Lists")
    with pysam.VariantFile(ref_vcf) as vcf_in, pysam.VariantFile(output_vcf, "w", header=vcf_in.header) as vcf_out:
        last_pset = None
        block_list = []
        sn = probe_io.get_sample_name(vcf_in)
        for record in vcf_in:
            sd = record.samples[sn]
            if not sd.phased:
                vcf_out.write(record)
            else:
                ps = sd['PS']
                flip = 1

                # Walk through connection dictionary to find the final orientation
                while True:
                    if ps in connect:
                        ps, o = connect[ps]
                        flip = flip * o
                    else:
                        break

                if last_pset is None:
                    block_list.append(sd['PS'])
                elif sd['PS'] != last_pset:
                    if sd['PS'] != ps:
                        block_list.append(sd['PS'])
                    else:
                        stl.append("-".join([str(x) for x in block_list]))
                        block_list.clear()
                        block_list.append(sd['PS'])
                last_pset = sd['PS']

                if flip == -1:
                    record.samples[sn]['PS'] = ps
                    if sd['GT'][0] == 0 and sd['GT'][1] == 1:
                        record.samples[sn]['GT'] = (1, 0)
                    elif sd['GT'][0] == 1 and sd['GT'][1] == 0:
                        record.samples[sn]['GT'] = (0, 1)
                    else:
                        logging.error('Can not handle alternative genotpypes', record)
                        sys.exit(1)

                record.samples[sn]['PS'] = ps
                record.samples[sn].phased = True
                vcf_out.write(record)
    stl.append("-".join([str(x) for x in block_list]))
    pysam.tabix_index(output_vcf, preset='vcf', force=True)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(1)
