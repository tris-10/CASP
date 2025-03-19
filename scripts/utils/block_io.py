"""
Phase block manipulation
"""


from enum import Enum
from utils import probe_io
import copy


class BlockAssign(Enum):
    HAP1 = 10,
    HAP2 = 20,
    MIXED = 30,
    UNKNOWN = 40

class PhaseBlock:
    """
    Information about phased/unphased blocks
    """

    def __init__(self, phase_set_name, contig, start, end):
        self.contig = contig
        self.name = phase_set_name
        self.start = int(start)
        self.end = int(end)
        self.length = self.end - self.start
        self.snp_list = []
        self.vcf_pos = []
        self.homozygous = False
        self.hap1_cnt = 0
        self.hap2_cnt = 0
        self.hap_string = None
        self.phase_set_self = [0, 0]

    def add_self_phase_set(self, index):
        self.phase_set_self[index] += 1

    def add_snp(self, snp):
        self.snp_list.append(snp)

    def add_vcf_pos(self, vcf_pos):
        self.vcf_pos.append(vcf_pos)

    def add_snps(self, snps):
        self.snp_list = copy.copy(snps)

    def count_haplotypes(self):
        """
        Sort the valid SNPs assigned to the block and check the haplotypes
        assignments
        """

        valid_snps = [snp for snp in self.snp_list if snp.valid]
        sorted_snps = sorted(valid_snps, key=lambda x: x.contig_pos)

        self.hap1_cnt = self.hap2_cnt = 0
        self.hap_string = ""

        # Check the haplotypes of all valid SNPs in the block
        for snp in sorted_snps:
            if snp.haplotype == probe_io.ProbeAssign.HAP1:
                self.hap1_cnt += 1
                self.hap_string += '1'
            elif snp.haplotype == probe_io.ProbeAssign.HAP2:
                self.hap2_cnt += 1
                self.hap_string += '2'

    def classify_block(self, min_count, max_contam_frac):
        """
        Determine the block haplotype.  Return unknown if there isn't enough data and mixed if both
        haplotypes are present at high levels.

        :param min_count: Minimum number of haplotype observations for assignment
        :param max_contam_frac: Maximum amount of haplotype contamination for assignmetn
        :return:
        """

        self.count_haplotypes()

        if self.hap1_cnt + self.hap2_cnt < min_count:
            return BlockAssign.UNKNOWN

        hap1p = self.hap1_cnt / (self.hap1_cnt + self.hap2_cnt)
        hap2p = self.hap2_cnt / (self.hap1_cnt + self.hap2_cnt)

        if self.hap1_cnt > self.hap2_cnt and self.hap1_cnt >= min_count:
            if hap1p >= (1-max_contam_frac):
                return BlockAssign.HAP1
            else:
                return BlockAssign.MIXED
        elif self.hap2_cnt > self.hap1_cnt and self.hap2_cnt >= min_count:
            if hap2p >= (1-max_contam_frac):
                return BlockAssign.HAP2
            else:
                return BlockAssign.MIXED
        else:
            return BlockAssign.UNKNOWN

    def split_phased_block(self, limit):
        """
        If the haplotype switches in the middle of phase block, split up the larger blocks into sub blocks

        :param limit: Maximum number of allowed splits.  If more, don't split and the block will be unassigned
        :return: List of new blocks
        """

        start_pos = last_pos = self.start
        new_phase_blocks = []
        new_hom_blocks = []
        curr_snps = []
        last_hap = None
        for snp in sorted(self.snp_list, key=lambda x: x.contig_pos):
            if last_hap is not None and last_hap != snp.haplotype:
                name_p = "{0}:SPLT_HET_U:{1}-{2}".format(self.contig, start_pos, last_pos)
                b = PhaseBlock(name_p, self.contig, start_pos, last_pos)
                b.add_snps(curr_snps)
                new_phase_blocks.append(b)

                start_pos = snp.contig_pos
                curr_snps.clear()

                name_h = "{0}:SPLT_HET_M:{1}-{2}".format(self.contig, last_pos, snp.contig_pos)
                hb = PhaseBlock(name_h, self.contig, last_pos, snp.contig_pos)
                new_hom_blocks.append(hb)

            last_hap = snp.haplotype
            last_pos = snp.contig_pos
            curr_snps.append(snp)

        name_p = "{0}:SPLT_HET_D:{1}-{2}".format(self.contig, start_pos, self.end)
        b = PhaseBlock(name_p, self.contig, start_pos, self.end)
        b.add_snps(curr_snps)
        new_phase_blocks.append(b)

        if len(new_phase_blocks) > limit:
            return None
        else:
            return new_phase_blocks, new_hom_blocks

    def split_hemi_block(self, max_contam, min_hets):
        """
        Split up hemizyguos blocks if the haplotype switches once.

        :param min_snps: Minimum number of SNPs required to intiate splitting
        :param max_contam: Maximum amount of contamination allowed before the hemizygous block is split.
        :return:
        """

        self.count_haplotypes()

        if self.hap1_cnt + self.hap1_cnt < min_hets:
            return [self]

        contam = min(self.hap1_cnt, self.hap2_cnt) / (self.hap1_cnt + self.hap2_cnt)
        if contam < max_contam:
            return [self]

        scores = []
        for i in range(1, len(self.hap_string)):
            f1 = self.hap_string[:i].count("1") / self.hap1_cnt
            f2 = self.hap_string[:i].count("2") / self.hap2_cnt
            b1 = self.hap_string[i:].count("1") / self.hap1_cnt
            b2 = self.hap_string[i:].count("2") / self.hap2_cnt

            avg1 = sum([f1, b2]) / 2
            avg2 = sum([f2, b1]) / 2
            scores.append(max(avg1, avg2))

        max_score = 1-contam
        cut_idx = None
        for i, score in enumerate(scores):
            if score > max_score:
                cut_idx = i+1
                max_score = score

        if cut_idx is None:
            return [self]

        valid_snps = [snp for snp in self.snp_list if snp.valid]
        sorted_snps = sorted(valid_snps, key=lambda x: x.contig_pos)

        cut_start = sorted_snps[cut_idx - 1].contig_pos
        cut_end = sorted_snps[cut_idx].contig_pos

        up_block = PhaseBlock("{0}:SPLT_HEM_U:{1}-{2}".format(self.name, self.start, cut_start), self.contig, self.start, cut_start)
        up_block.add_snps(sorted_snps[:cut_idx])
        down_block = PhaseBlock("{0}:SPLT_HEM_D:{1}-{2}".format(self.name, cut_end, self.end), self.contig, cut_end, self.end)
        down_block.add_snps(sorted_snps[cut_idx:])
        mid_block = PhaseBlock("{0}:SPLT_HEM_M:{1}-{2}".format(self.name, cut_start, cut_end), self.contig, cut_start, cut_end)

        new_blocks = [up_block, down_block, mid_block]
        return new_blocks

    def split_hom_block(self, pred_hemi_regions, min_het_count, min_split=5000):
        """
        Split hemizygous blocks into hemi and homozygous blocks if the ends of the hemizygous blocks don't have
        SNPs and don't overlap with any other contig

        :param pred_hemi_regions: List of predicted hemizygous regions
        :param min_het_count: The minimum number of hets required to assign a block hemizygous
        :param min_split: The minimum distance from end of block to initiate splitting
        :return:
        """
        hemi_block = None
        hom_regions = []
        if len([x for x in self.snp_list if x.valid]) < min_het_count:
            return hemi_block, hom_regions
        else:
            sorted_snps = sorted(self.snp_list, key=lambda x: x.contig_pos)

            hemi_start = sorted_snps[0].contig_pos
            hemi_end = sorted_snps[-1].contig_pos
            if (hemi_start - self.start) < min_split or (
                    self.contig in pred_hemi_regions and PhaseBlock.full_overlap(self.start, sorted_snps[0].contig_pos,
                                                                             pred_hemi_regions[self.contig])):
                hemi_start = self.start
            else:
                name = "{0}:TRIM_HOM_U:{1}-{2}".format(self.name, self.start, sorted_snps[0].contig_pos)
                hom_regions.append(PhaseBlock(name, self.contig, self.start, sorted_snps[0].contig_pos))

            if (self.end - hemi_end) < min_split or (
                    self.contig in pred_hemi_regions and PhaseBlock.full_overlap(sorted_snps[-1].contig_pos, self.end,
                                                                             pred_hemi_regions[self.contig])):
                hemi_end = self.end
            else:
                name = "{0}:TRIM_HOM_D:{1}-{2}".format(self.name, sorted_snps[-1].contig_pos, self.end)
                hom_regions.append(PhaseBlock(name, self.contig, sorted_snps[-1].contig_pos, self.end))

            mid_class = "NONE"
            if self.contig in pred_hemi_regions and PhaseBlock.full_overlap(hemi_start, hemi_end,
                                                                        pred_hemi_regions[self.contig]):
                mid_class = "FULL"
            elif self.contig in pred_hemi_regions and PhaseBlock.partial_overlap(hemi_start, hemi_end,
                                                                             pred_hemi_regions[self.contig]):
                mid_class = "PART"

            name = "{0}:HEM_{1}:{2}-{3}".format(self.name, mid_class, hemi_start, hemi_end)
            hemi_block = PhaseBlock(name, self.contig, hemi_start, hemi_end)
            hemi_block.snp_list = self.snp_list
            hemi_block.vcf_pos = self.vcf_pos
            return hemi_block, hom_regions

    @staticmethod
    def full_overlap(start, end, regions):
        size = end - start
        for region in regions:
            overlap = min(end, region[1]) - max(start, region[0])
            if overlap / size > 0.90:
                return True
        else:
            return False

    @staticmethod
    def partial_overlap(start, end, regions):
        for region in regions:
            if start <= region[1] and end >= region[0]:
                return True
        else:
            return False



