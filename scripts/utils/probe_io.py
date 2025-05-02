"""
Reference SNV to assembly classes
"""

from enum import Enum
from collections import defaultdict
import pysam
import copy
import sys
import logging

comp = {"A": "T", "T": "A", "C": "G", "G": "C", "X": "X"}

logging.basicConfig(format='%(levelname)s:%(asctime)s %(message)s', level=logging.INFO)


class ProbeAssign(Enum):
    HAP1 = 1,
    HAP2 = 2,
    NONE = 3


class ProbeInfo:
    """
    This class stores information about the heterozygous positions expected to be captured
    in the assembly.
    """

    def __init__(self, hg38_pos, base1, base2):
        self.hg38_pos = hg38_pos
        self.base1 = base1
        self.base2 = base2


class ProbeAlign:
    """
    Class handles information on where a SNP-CHIP probe aligns in the assembly.
    """

    def __init__(self, snp_id, contig_pos, contig_name, is_reverse, probe_info):
        self.snp_id = snp_id
        self.contig_pos = contig_pos
        self.contig_name = contig_name
        self.loc = '{0}:{1}'.format(contig_name, contig_pos)
        self.snp_contig = '{0}:{1}'.format(snp_id, contig_name)
        self.is_reverse = is_reverse
        self.probe_info = probe_info
        self.haplotype = None
        self.valid = True

    def set_haplotype(self, hap):
        self.haplotype = hap

    def get_base(self, is_base1):
        base = self.probe_info.base1 if is_base1 else self.probe_info.base2
        return comp[base] if self.is_reverse else base

    def set_invalid(self):
        self.valid = False


def load_probe_alignment(alignment_file, probe_info, is_front, max_mm):
    """
    Read through probe/het alignment file.  If the flanking sequence aligns to the contig within max_mm
    mismatches and without trimming, the location on the contig is stored and returned in a dictionary

    :param alignment_file: Probe/Het flanking sequence alignment to contig in bam format
    :param probe_info: Dictionary containing all valid probe/het information
    :param is_front: True if flanking sequence is upstream of probe/het false if downstream
    :param max_mm: Maximum number of mismatches to allow against the contig
    :return: Dictionary of probe/het locations keyed on the alignment location
    """

    probe_align_dict = {}

    with pysam.AlignmentFile(alignment_file, "r") as bam_in:
        for a in bam_in:
            # skip probes not in the VCF
            if a.query_name not in probe_info:
                continue

            # skip unaligned probes or probes with trimming
            if a.is_unmapped or len([cigar for cigar in a.cigartuples if cigar[0] == 5 or cigar[0] == 4]) > 0:
                continue

            # Drop probes with more than two mismatches, unless spike-in.
            if a.get_tag("NM") > max_mm:
                continue

            # calculate the position of the SNP on the contig based on front/back probe and alignment orientation
            if is_front:
                pos = a.reference_start if a.is_reverse else a.reference_end + 1
            else:
                pos = a.reference_end + 1 if a.is_reverse else a.reference_start

            # create probe alignment location object
            contig = bam_in.get_reference_name(a.reference_id)
            key = '{0}:{1}'.format(contig, pos)
            probe_align = ProbeAlign(a.query_name, pos, bam_in.get_reference_name(a.reference_id), a.is_reverse,
                                     copy.copy(probe_info[a.query_name]))
            probe_align_dict[key] = probe_align

    return probe_align_dict


def identify_valid_probe_locations(up_probe_dict, down_probe_dict, contig_data, max_alignments):
    """
    Match up the upstream/downstream probe mapping location. If both flanking sequences map to the same location
    and are facing the same way the position is valid.  Probes that match more than max_alignments are remove.
    Probes with two valid mapping locations to the same contig are removed. It is possible that the same
    probe het can map to 2 or maybe more overlapping contigs, not to the same contig.

    :param up_probe_dict: Dictionary containing upstream probe mapping information
    :param down_probe_dict: Dicitionary containing downstream probe mapping information
    :param contig_data: Dictionary containing the collapse contig sequences
    :param max_alignments: Maximum number of times a probe can align to an assembly across both the up/downstream seqs
    :return: Dictionary with valid mapping locations
    """

    valid_probe_dict = {}
    probe_count = defaultdict(int)
    snp_contig_counts = defaultdict(int)

    for d in [up_probe_dict, down_probe_dict]:
        for key, probe in d.items():
            probe_count[probe.snp_id] += 1

    for loc, probe in up_probe_dict.items():
        # If the heterozyguous is not found in the downstream skip
        if loc not in down_probe_dict or probe.is_reverse != down_probe_dict[loc].is_reverse or probe.snp_id != \
                down_probe_dict[loc].snp_id:
            continue

        # If the probe is mapped may times to assembly, skip
        if probe_count[probe.snp_id] > max_alignments:
            continue

        # Count how many times a snp is found with a contig
        snp_contig_counts[probe.snp_contig] += 1

        base = contig_data[probe.contig_name][probe.contig_pos - 1]
        if base == probe.get_base(True):
            probe.set_haplotype(ProbeAssign.HAP1)
        elif base == probe.get_base(False):
            probe.set_haplotype(ProbeAssign.HAP2)
        else:
            probe.set_haplotype(ProbeAssign.NONE)

        valid_probe_dict[loc] = probe

    # Remove probes that map to the same contig multiple times.
    bad_locs = [loc for loc, probe in valid_probe_dict.items() if snp_contig_counts[probe.snp_contig] > 1]
    for bad_loc in bad_locs: del valid_probe_dict[bad_loc]

    return valid_probe_dict


def load_probe_data(phased_vcf):
    """
    Load phased VCF file and return a dictionary of phased heterozygous SNPs keyed on their ID

    :param phased_vcf: Path to phased vcf file using hg38 coordinates
    :return: dictionary of ProbeInfo objects, keyed on SNP ID
    """

    probe_info = {}
    with pysam.VariantFile(phased_vcf) as vcf_in:
        sn = get_sample_name(vcf_in)
        for record in vcf_in:
            sd = record.samples[sn]
            gt = sd["GT"]
            if sd.phased:
                if gt[0] == 0 and gt[1] == 1:
                    probe_info[record.id] = ProbeInfo(record.pos, record.ref, record.alts[0])
                elif gt[0] == 1 and gt[1] == 0:
                    probe_info[record.id] = ProbeInfo(record.pos, record.alts[0], record.ref)
    return probe_info


def get_sample_name(vcf_in):
    """
    Check to make sure the VCF file only has a single sample and return its name

    :param vcf_in: Open file handle to a VCF file
    :return: string sample name.
    """
    samples = list(vcf_in.header.samples)
    if len(samples) != 1:
        logging.error('Only one sample is allowed in VCF, found: {0}'.format(" ".join(samples)))
        sys.exit(1)
    return samples[0]
