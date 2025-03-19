#!/usr/bin/env python

"""
Restrict VCF to heterozygous variants that are found in a list of known variants and fall outside of the RCCX region.
Optionally filter based on mapping quality and depth.
"""

from utils import probe_io
import sys
import argparse
import pysam
import logging

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_vcf", help="Path to variant calls in VCF format.")
    parser.add_argument("known_variants", help="Path to a list of valid variants in the format: "
                                               "SNP_ID\tCHR-POS-REF-ALT\n")
    parser.add_argument("output_vcf", help="Path to filtered VCF file.")
    parser.add_argument("output_stats", help="Path to VCF filtering stats.")
    parser.add_argument("-m", "--min_qual", type=float, default=0, help="Minimum mapping quality to retain.")
    parser.add_argument("-d", "--min_depth", type=float, default=0, help="Minimum depth of coverage to retain.")
    parser.add_argument("-a", "--min_allele_count", type=int, default=0, help="Minimum allele count to consider heterozygous")
    parser.add_argument("-f", "--min_allele_freq", type=float, default=0, help="Minimum allele frequency to consider heterozygous")
    parser.add_argument("-p", "--remove_phase", default=False, action="store_true",
                        help="Remove phasing information from the VCF file.")
    args = parser.parse_args()

    logging.basicConfig(format='%(levelname)s:%(asctime)s %(message)s', level=logging.INFO)

    valid_positions = {}
    logging.info("Loading known variants")
    with open(args.known_variants, "r") as ipf:
        for line in ipf:
            items = line.strip().split("\t")
            valid_positions[items[1]] = items[0]

    logging.info("Filtering called variants")
    dropped = kept = homozygous = qual = total = depth = 0
    with pysam.VariantFile(args.input_vcf) as vcf_in, pysam.VariantFile(args.output_vcf, "w", header=vcf_in.header) as vcf_out:
        for r in vcf_in:
            total += 1
            sn = probe_io.get_sample_name(vcf_in)
            sd = r.samples[sn]
            gt = sd["GT"]
            ac = r.info["AC"]

            if (gt[0] == 1 and gt[1] == 1) or (min(ac) < args.min_allele_count) or (min(ac) / sum(ac) < args.min_allele_freq):
                homozygous += 1
                continue
            if r.qual < args.min_qual:
                qual += 1
                continue
            if sum(ac) < args.min_depth:
                depth += 1
                continue

            if args.remove_phase:
                sd['PS'] = None
                sd.phased = False

            ref = r.ref.upper()
            alt = r.alts[0].upper()

            key = "{0}-{1}-{2}-{3}".format(r.chrom, r.pos, ref, alt)
            if key in valid_positions:
                r.id = valid_positions[key]
                vcf_out.write(r)
                kept += 1
            else:
                dropped += 1

    pysam.tabix_index(args.output_vcf, preset='vcf', force=True)

    with open(args.output_stats, "w") as op_stats:
        op_stats.write("total_gnomad_var_count\t{0}\n".format(len(valid_positions.keys())))
        op_stats.write("initial_var_count\t{0}\n".format(total))
        op_stats.write("homozygous_var_count\t{0}\n".format(homozygous))
        op_stats.write("low_quality_var_count\t{0}\n".format(qual))
        op_stats.write("low_depth_var_count\t{0}\n".format(depth))
        op_stats.write("non_gnomad_var_count\t{0}\n".format(dropped))
        op_stats.write("final_var_count\t{0}\n".format(kept))


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(1)

