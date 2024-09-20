#!/usr/bin/env python

"""
Load a fasta file and split sequences into subsequences of size --target_length. The final subsequence for each read
will typically be shorter than the target length and is always written.
"""

from Bio import SeqIO
import sys
import argparse


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_reads", help="Path to input reads in fasta format")
    parser.add_argument("output_reads", help="Path to split output reads in fasta format")
    parser.add_argument("-l", "--target_length", type=int, default=30000, help="Target subsequence length")
    args = parser.parse_args()

    with open(args.input_reads, "r") as fasta_in, open(args.output_reads, "w") as fasta_out:
        for record in SeqIO.parse(fasta_in, "fasta"):
            start = 0
            while True:
                end = start + args.target_length if start + args.target_length < len(record) else len(record)
                sub_record = record[start:end]
                sub_record.id = "{0}_{1}_{2}".format(sub_record.id, str(start), str(end))
                sub_record.description = ""
                SeqIO.write(sub_record, fasta_out, "fasta")
                start = end
                if end == len(record):
                    break


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(1)
