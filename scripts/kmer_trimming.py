#!/usr/bin/env python

"""
Trim corrected ONT reads using k-mer counts.  For each k-mer in a given read, look up the k-mer depth in the full dataset.
Starting at the ends of each read, check if the k-mer count is --min_frac_cov the median read k-mer depth.  Reads are
cut to the starting position of the first k-mer above the --min_frac_cov threshold.  If more than --min_frac_low_kmer
positions are below --min_frac_cov, the entire read is discarded.  Reads shorter than --min_read_length are discarded.
"""

from Bio import SeqIO
from collections import defaultdict
import sys
import argparse
import numpy as np
import logging

trans_map = str.maketrans('ACGT', 'TGCA')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fasta_file', help='Path to corrected ONT reads in fasta format')
    parser.add_argument('kmer_file', help='Path to kmer count file generated by jellyfish')
    parser.add_argument('trim_file', help='Path to trimmed reads output file in fasta format')
    parser.add_argument("-c", "--min_frac_cov", type=float, default=0.1, help='Minimum allowed k-mer depth relative '
                                                                              'to the median read k-mer depth')
    parser.add_argument("-f", "--max_frac_low_kmer", type=float, default=0.01, help='Max fraction of the read that can '
                                                                        'be below --min_frac_cov before discarding')
    parser.add_argument("-l", "--min_read_length", type=int, default=5000, help='Minimum read length retained')
    args = parser.parse_args()

    logging.basicConfig(format='%(levelname)s:%(asctime)s %(message)s', level=logging.INFO)

    logging.info("Loading k-mer depth dictionary")
    kmer_dict = load_kmer_counts(args.kmer_file)
    if len(kmer_dict.keys()) == 0:
        logging.error("k-mer depth file is empty")
        sys.exit(1)
    kmer_flank = int(len(list(kmer_dict.keys())[0]) / 2)

    logging.info("Start read trimming")
    trim_reads(args.fasta_file, args.trim_file, kmer_dict, kmer_flank, args.min_frac_cov, args.max_frac_low_kmer,
               args.min_read_length)


def trim_reads(input_fasta, output_fasta, kmer_dict, kmer_flank, min_frac_cov, max_frac_low_kmer, min_read_length):
    """
    Calculate the median k-mer depth across each read.  Trim the ends of the reads if they contain k-mers with depths
    'min_frac_cov' of the median k-mer depth.  Discard reads with 'min_frac_low_kmer' positions below the threshold or
    shorter than 'min_read_length'

    :param input_fasta: Path to corrected ONT reads in fasta format
    :param output_fasta: Path to trimmed ONT reads in fasta format
    :param kmer_dict: k-mer depth dictionary
    :param kmer_flank: 1/2 k-mer length
    :param min_frac_cov: Minimum fraction of median k-mer depth allowed at read ends
    :param max_frac_low_kmer: Minimum fraction of read bases below threshold allowed
    :params min_read_length: Minimum retained read length

    """
    short = low_kmer = trimmed = 0
    with open(input_fasta, "r") as fasta_in, open(output_fasta, "w") as fasta_out:
        for record in SeqIO.parse(fasta_in, "fasta"):
            if len(record.seq) < min_read_length:
                continue

            kmer_counts = kmer_counts_across_read(str(record.seq), kmer_flank, kmer_dict)
            median_cov = np.median(kmer_counts)
            min_depth = int(median_cov * min_frac_cov + 1)

            bounds = [0, 0]
            for b, counts in zip([0, 1], [kmer_counts, kmer_counts[::-1]]):
                for idx, val in enumerate(counts):
                    if val >= min_depth:
                        bounds[b] = idx
                        break

            start = bounds[0]
            end = len(record.seq) - bounds[1]
            read_length = end - start

            if read_length >= min_read_length:
                f_low = len([x for x in kmer_counts[start:end] if x < min_depth]) / read_length
                if f_low <= max_frac_low_kmer:
                    record.seq = record.seq[start:end]
                    SeqIO.write(record, fasta_out, "fasta")
                    trimmed += 1
                else:
                    low_kmer += 1
            else:
                short += 1
    logging.info("Trimmed reads {0}, Dropped reads {1}, Short reads {2}".format(trimmed, low_kmer, short))


def load_kmer_counts(kmer_count_file):
    """
    Load in kmer counts generated by jellyfish. Each line contains the kmer sequence followed by the count of the kmer
    in the dateset. Low count kmers are not included in the output to reduce filesize.

    :param kmer_count_file: Path to a kmer count file generated by jellyfish
    :return: Dictionary of kmer counts, key is the kmer sequence, value is the count of the kmer in the dataset.
    """

    kmer_count_dict = defaultdict(int)
    with open(kmer_count_file, "r") as ipf:
        for line in ipf:
            items = line.strip().split(" ")
            kmer_count_dict[items[0]] = int(items[1])
    return kmer_count_dict


def generate_kmer(seq, pos, kmer_flank):
    """
    A kmer is generated by retrieving the sequence kmer_flank upstream and downstream of the position.  The
    reverse complement sequence is generated and the two are compared, retaining the lowest lexicographically.

    :param seq: Raw PacBio sequence
    :param pos: Potential breakpoint position
    :param kmer_flank: kmer sizes are odd, flank is kmer size / 2 rounded down the nearest integer
    :return:
    """
    kmer_for = seq[pos - kmer_flank: pos + kmer_flank + 1]
    kmer_rev = reverse_complement(kmer_for)
    return kmer_for if kmer_for < kmer_rev else kmer_rev


def kmer_counts_across_read(seq, kmer_flank, kmer_dict):
    """
    Generate a dictionary containing the count of every kmer across a read

    :param seq: Raw PacBio sequence
    :param kmer_flank: kmer sizes are odd, flank is kmer size / 2 rounded down the nearest integer
    :return: dictionary containing the count of every kmer in the read
    """
    kmer_self = []
    for i in range(kmer_flank, len(seq) - kmer_flank):
        kmer_self.append(kmer_dict[generate_kmer(seq, i, kmer_flank)])
    return kmer_self


def reverse_complement(seq):
    """
    Reverse complement the input sequence.  Assumes no ambiguous bases and uppercase letters.

    :param seq: DNA sequence
    :return: reverse complement sequence
    """

    return seq.translate(trans_map)[::-1]


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(1)