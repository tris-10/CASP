#!/usr/bin/env python

"""
Filter ONT reference alignment by removing reads with low mapping quality or that fall completely within the
RCCX CNV.
"""

import pysam
import argparse
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_alignment', help="Path to the reference sequence alignment in bam format.")
    parser.add_argument('output_alignment', help="Path the the filtered alignment output in bam format.")
    parser.add_argument('output_stats', help="Path the the filtered alignment statistics.")
    parser.add_argument("-m", "--min_mapq", type=int, default=0, help="Minimum mapping quality.")
    parser.add_argument("-s", "--rccx_start", type=int, default=31981000, help="Starting coordinate compliment")
    parser.add_argument("-e", "--rccx_end", type=int, default=32047000, help="Ending coordinate compliment")
    args = parser.parse_args()

    total = comp = qual = 0
    with pysam.AlignmentFile(args.input_alignment, 'rb') as bam_in, \
            pysam.AlignmentFile(args.output_alignment, 'wb', template=bam_in) as bam_out:
        for read in bam_in:
            if read.is_unmapped:
                continue

            total += 1
            if read.mapping_quality < args.min_mapq:
                qual += 1
                continue

            if read.reference_start > args.rccx_start and read.reference_end < args.rccx_end:
                comp += 1
                continue

            bam_out.write(read)

    with open(args.output_stats, "w") as op_stats:
        op_stats.write("total_align\t{0}\n".format(total))
        op_stats.write("rccx_align\t{0}\n".format(comp))
        op_stats.write("low_mapq\t{0}\n".format(qual))


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(1)
