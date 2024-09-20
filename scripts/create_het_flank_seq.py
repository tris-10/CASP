#!/usr/bin/env python


"""
Read through phased VCF and write out the flanking sequence for each het given a reference sequence
"""

from Bio import SeqIO
import argparse
import sys
import pysam
import logging


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("vcf_file", help="Path to phased VCF file.")
    parser.add_argument("ref_file", help="Path to reference sequence in fasta format.")
    parser.add_argument("probe_front", help="Path to upstream flanking sequence output file.")
    parser.add_argument("probe_back", help="Path to downstream flanking sequence output file.")
    parser.add_argument("-l", "--flank_length", help="Length of the flanking sequence.", default=60, type=int)
    args = parser.parse_args()

    logging.basicConfig(format='%(levelname)s:%(asctime)s %(message)s', level=logging.INFO)

    logging.info("Loading reference file")
    chrom_dict = SeqIO.to_dict(SeqIO.parse(args.ref_file, "fasta"))

    logging.info("Creating flanking sequences")
    with pysam.VariantFile(args.vcf_file, "r") as vcf_in, open(args.probe_front, 'w') as probe_front_out, \
            open(args.probe_back, 'w') as probe_back_out:
        for var in vcf_in:
            sd = var.samples["SAMPLE"]
            gt = sd["GT"]
            if sd.phased and ((gt[0] == 0 and gt[1] == 1) or (gt[0] == 1 and gt[1] == 0)):
                probe_front_out.write(">{0}\n{1}\n".format(var.id, str(chrom_dict[var.contig][var.pos-1 - args.flank_length: var.pos-1].seq)))
                probe_back_out.write(">{0}\n{1}\n".format(var.id, str(chrom_dict[var.contig][var.pos: var.pos + args.flank_length].seq)))


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logging.error("User interrupted, exiting")
        sys.exit(1)