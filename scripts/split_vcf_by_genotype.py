#!/usr/bin/env python

"""
This script filters the longshot VCF files to remove variants: 1) with low quality scores, 2) within homopolymer
stretches and 3) with low depth of coverage. Phasing information is removed and heterozygous/homozygous mutant
variants are written to separate files
"""

from utils import probe_io
import pysam
import re
import sys
import argparse
import logging


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_vcf', help="Path to Longshot VCF file.")
    parser.add_argument('het_vcf', help="Path to heterozygous variant file output.")
    parser.add_argument('hom_vcf', help="Path to homozygous variant file output.")
    parser.add_argument('output_stats', help='Path to vcf splitting statistics output')
    parser.add_argument('-c', '--min_het_count', default=8, type=int,
                        help="Minimum number of allele counts for both alleles to retain heterozygous position.")
    parser.add_argument('-q', '--min_var_qual', default=150, type=int,
                        help="Minimum variant qual score to maintain position.")
    parser.add_argument("-l", '--min_hp_length', default=5, type=int)
    args = parser.parse_args()

    logging.basicConfig(format='%(levelname)s:%(asctime)s %(message)s', level=logging.INFO)

    total = het_count = hom_alt_count = other_count = qual_count = hp_count = 0

    logging.info("Splitting and filtering contig VCF file")
    with pysam.VariantFile(args.input_vcf, 'r') as vcf_in, pysam.VariantFile(args.het_vcf, 'w', header=vcf_in.header) as vcf_out_het, \
            pysam.VariantFile(args.hom_vcf, 'w', header=vcf_in.header) as vcf_out_hom:

        sn = probe_io.get_sample_name(vcf_in)
        for var in vcf_in:
            total += 1

            # skip over variants with low quality
            if var.qual < args.min_var_qual:
                qual_count += 1
                continue

            # Extract out allele counts and genotypes
            sd = var.samples[sn]
            gt = sd['GT']
            ac = var.info["AC"]

            # Check sequence context
            if check_context(var.info["SC"], var.alts[0], args.min_hp_length):
                hp_count += 1
                continue

            # remove phasing information
            sd['PS'] = None
            sd.phased = False

            if (gt[0] == 1 and gt[1] == 0) or (gt[0] == 0 and gt[1] == 1): # Handle heterozygous
                if ac[0] < args.min_het_count:  # If the reference allele is below threshold, set as homozygous mutant
                    sd['GT'] = (1, 1)
                    vcf_out_hom.write(var)
                    hom_alt_count += 1
                elif ac[1] < args.min_het_count:  # If alternate allele is below threshold, drop
                    continue
                else:  # Otherwise retain
                    sd['GT'] = (0, 1)
                    vcf_out_het.write(var)
                    het_count += 1
            elif gt[0] == 1 and gt[1] == 1:  # Handle homozygous
                sd['GT'] = (1, 1)
                vcf_out_hom.write(var)
                hom_alt_count += 1
            else:  # Drop any unexpected genotypes
                other_count += 1

    pysam.tabix_index(args.hom_vcf, preset='vcf', force=True)
    pysam.tabix_index(args.het_vcf, preset='vcf', force=True)

    with open(args.output_stats, "w") as op_stats:
        op_stats.write("total_var\t{0}\n".format(total))
        op_stats.write("total_homopolymer\t{0}\n".format(hp_count))
        op_stats.write("low_qual\t{0}\n".format(qual_count))
        op_stats.write("het_var\t{0}\n".format(het_count))
        op_stats.write("hom_alt_var\t{0}\n".format(hom_alt_count))
        op_stats.write("other_var\t{0}\n".format(other_count))


def check_context(sc, alt, hp_length):
    """
    Check if variant is part of a homopolymer stretch

    params: sc: Sequence context as reported by longshot
    params: alt: Alternate allele call
    params: hp_length: Minimum homopolymer length to flag
    """

    upstream = sc[:10][::-1]
    downstream = sc[11:]

    pattern = alt + "{" + str(hp_length) + ",}"

    for flank in upstream, downstream:
        if re.match(pattern, flank):
            return True
    return False


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(1)

