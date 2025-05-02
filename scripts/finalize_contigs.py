#!/usr/bin/env python

"""
Finalize the patched assembly.  Find short contigs in haplotype A that are very strong matches to haplotype B and
remove. Trim assembies so that contig start at capture boundries (ie. MOG and DAXX). Stitch together contigs
within each overlap with overlap ends and high identity
"""

from pathlib import Path
from Bio import SeqIO
from utils import paf_io as pio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import argparse
import os
import subprocess
import logging

stl = []


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('hap1_contig_path', help="Path to hap1 patched assembly in fasta format.")
    parser.add_argument('hap2_contig_path', help="Path to hap2 patched assembly in fasta format.")
    parser.add_argument('boundary_sequences', help="Path to the boundary sequences fasta.")
    parser.add_argument('final_hap1_output', help="Path to hap1 final assembly output.")
    parser.add_argument('final_hap2_output', help="Path to hap2 final assembly output.")
    parser.add_argument('stats_output', help="Path to flitering/trimming stats output file.")
    parser.add_argument("-i", "--min_identity", type=float, default=0.999,
                        help="Minimum identity between two contigs to consider for removal or stitching")
    parser.add_argument("-p", "--minimap2_path", default="minimap2", help="Path to minimap2 binary")
    parser.add_argument('-t', '--threads', type=int, default=4, help="Minimap2 threads")
    parser.add_argument('-f', '--min_overlap_frac', type=float, default=0.999,
                        help="Minimum overlap required to remove contig or identify boundaries.")
    parser.add_argument("-m", "--max_contig_size", type=int, default=200000,
                        help="Maximum contig length that can be removed.")
    parser.add_argument("-e", "--max_edge_distance", type=int, default=1000,
                        help="Maximum distance from contig end to consider for overlap stitching")
    parser.add_argument('-l', "--min_overlap_length", type=int, default=5000,
                        help="Minimum overlap between two contigs required for stitiching")

    args = parser.parse_args()

    logging.basicConfig(format='%(levelname)s:%(asctime)s %(message)s', level=logging.INFO)

    logging.info("Removing bubbles")
    stl.append("Bubbles\nHaplotype\tTigName")
    contigs_to_remove = [set(), set()]
    remove_bubbles(args.hap1_contig_path, 0, contigs_to_remove)
    remove_bubbles(args.hap2_contig_path, 1, contigs_to_remove)

    logging.info("Align haplotype assemblies to each other")
    stl.append("\nOverlap Removal\nHap\tFiltTigName\tFitTigLength\tOther\tOvlTigName\tOvlTigLength\tIdentity\tOverlap")
    remove_hap_overlaps(args.hap1_contig_path, args.hap2_contig_path, 0, 1, contigs_to_remove, args.final_hap1_output,
                        args.max_contig_size, args.min_identity, args.min_overlap_frac, args.minimap2_path,
                        args.threads)
    remove_hap_overlaps(args.hap2_contig_path, args.hap1_contig_path, 1, 0, contigs_to_remove, args.final_hap2_output,
                        args.max_contig_size, args.min_identity, args.min_overlap_frac, args.minimap2_path,
                        args.threads)

    logging.info("Align haplotype assemblies to self")
    stl.append("\nOverlap Removal\nHap\tFiltTigName\tFitTigLength\tSelf\tOvlTigName\tOvlTigLength\tIdentity\tOverlap")
    remove_hap_overlaps(args.hap1_contig_path, args.hap1_contig_path, 0, 0, contigs_to_remove, args.final_hap1_output,
                        args.max_contig_size, args.min_identity, args.min_overlap_frac, args.minimap2_path,
                        args.threads)
    remove_hap_overlaps(args.hap2_contig_path, args.hap2_contig_path, 1, 1, contigs_to_remove, args.final_hap2_output,
                        args.max_contig_size, args.min_identity, args.min_overlap_frac, args.minimap2_path,
                        args.threads)

    logging.info("Remove filtered contigs")
    removed_filtered(args.hap1_contig_path, args.final_hap1_output, contigs_to_remove, 0)
    removed_filtered(args.hap2_contig_path, args.final_hap2_output, contigs_to_remove, 1)

    logging.info("Stitching together assemblies")
    stl.append("\nStitched Contigs\nHaplotype\tTig1\tTig2\tOverlap\tIdentity")
    stitch_contigs(args.final_hap1_output, 0, args.max_edge_distance, args.min_overlap_length, args.min_identity,
                   args.minimap2_path, args.threads)
    stitch_contigs(args.final_hap2_output, 1, args.max_edge_distance, args.min_overlap_length, args.min_identity,
                   args.minimap2_path, args.threads)

    logging.info("Trimming assemblies to the capture boundaries")
    stl.append("\nBoundary Trimming\nHaplotype\tTigName\tBoundary\tRemovedBases\tOrient")
    trim_contigs_to_boundaries(args.final_hap1_output, 0, args.boundary_sequences, args.minimap2_path, args.threads)
    trim_contigs_to_boundaries(args.final_hap2_output, 1, args.boundary_sequences, args.minimap2_path, args.threads)

    with open(args.stats_output, "w") as op_stats:
        op_stats.write("\n".join(stl))


def removed_filtered(original_contigs, filtered_contigs, contigs_to_remove, idx):
    """
    Remove filtered contigs

    :params original_contigs: original contig file
    :params filtered_contigs: filtered contig file
    :params contigs_to-remove: sets of contigs to remove
    :params idx: haplotype index
    """

    with open(original_contigs, "r") as ipf, open(filtered_contigs, "w") as opf:
        for record in SeqIO.parse(ipf, "fasta"):
            if record.id not in contigs_to_remove[idx]:
                SeqIO.write(record, opf, "fasta")


def remove_bubbles(contigs, idx, contigs_to_remove):
    """
    Flag contigs for removal that are marked as bubbles by canu

    :params contigs: contig file
    :params idx: haplotype index
    :params contigs_to_remove: sets of contigs to remove
    """
    hap = get_hap_name(idx)
    with open(contigs, "r") as contig_in:
        for record in SeqIO.parse(contig_in, "fasta"):
            if "suggestBubble=yes" in record.description:
                stl.append("{0}\t{1}".format(hap, record.id))
                contigs_to_remove[idx].add(record.id)


def remove_hap_overlaps(hap1_contigs, hap2_contigs, idx1, idx2, contigs_to_remove, output_path, max_size, min_identity,
                        min_overlap_frac, minimap_path, threads):
    """
    Find overlaps between haplotype assemblies, remove contigs in hapA that completely overlap a longer contig in
    haplotype B with a high identity.

    :param hap1_contigs: Path to hap1 assembly in fasta format
    :param hap2_contigs: Path to hap2 assembly in fasta format
    :param idx1: Query idx
    :param idx2: Target idx
    :param output_path: Path to filtered output file in fasta format
    :param max_size: Maximum contig length that can be filtered out.
    :param min_identity: Minimum identity between hapA/hapB contigs required for removing shorter hapA contig.
    :param min_overlap_frac: Minimum overlap between hapA/hapB contigs required for removing shorter hapA contig.
    :param minimap_path: Path to minimap2 binary.
    :param threads: Number of minimap2 threads.
    """

    hap1 = get_hap_name(idx1)
    hap2 = get_hap_name(idx2)
    out_path = Path(output_path)
    align_file = out_path.parent.joinpath(out_path.stem + "_temp_ovl.paf")

    try:
        subprocess.run([minimap_path, '-x', 'map-hifi', '-c', '-t', str(threads), '-o', align_file, hap2_contigs,
                        hap1_contigs], check=True)

        align_paf = pio.load_single_alignments(align_file)

        for a in align_paf:
            ovl = (a.query_end - a.query_start) / a.query_length
            size_ratio = a.target_length / a.query_length

            if ovl >= min_overlap_frac and size_ratio > 2 and a.target_name not in contigs_to_remove[idx2] and \
                    ((a.totalNM < 5 and a.query_name.startswith("h1")) or (a.query_length < max_size and a.identityNM >= min_identity)):
                stl.append("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.4f}\t{7:.4f}".format(hap1, a.query_name, a.query_length,
                                                                                   hap2, a.target_name, a.target_length,
                                                                                   a.identityNM, ovl))
                contigs_to_remove[idx1].add(a.query_name)

    except subprocess.CalledProcessError as cpe:
        logging.error("Error overlapping assemblies: {0}".format(cpe.cmd))
        sys.exit(1)
    finally:
        if os.path.exists(align_file):
            os.remove(align_file)


def trim_contigs_to_boundaries(contig_file, idx, boundary_file, minimap2_path, threads):
    """
    Locate MOG and DAXX on the assembly and trim sequences that extend beyond these genes.  Write out all
    passing reads to a new file

    :param contig_file: Path to the assembly in fasta format
    :param haplotype: Name of the processed haplotype
    :param boundary_file: Path to the sequences of MOG and DAXX
    :param minimap2_path: Path to the minimap2 binary
    :param threads: number of threads to use for alignment
    """

    contig_path = Path(contig_file)
    align_file = contig_path.parent.joinpath(contig_path.stem + "_temp_bound.paf")
    hap = get_hap_name(idx)

    try:
        subprocess.run([minimap2_path, '-x', 'asm10', '-c', '-t', str(threads), '-o',
                        align_file, contig_file, boundary_file], check=True)

        updated_starts = {}
        updated_ends = {}
        updated_orient = set()
        for a in pio.load_single_alignments(align_file):
            ovl = (a.query_end - a.query_start) / a.query_length
            if ovl >= 0.99:
                if a.query_name == "START":
                    if a.strand == "-":
                        updated_orient.add(a.target_name)
                        updated_ends[a.target_name] = a.target_end
                    else:
                        updated_starts[a.target_name] = a.target_start

                elif a.query_name == "END":
                    if a.strand == "-":
                        updated_orient.add(a.target_name)
                        updated_starts[a.target_name] = a.target_start
                    else:
                        updated_ends[a.target_name] = a.target_end

        trimmed_contigs = []
        with open(contig_file) as ipf:
            for record in SeqIO.parse(ipf, "fasta"):
                start = 0
                end = len(record.seq)
                if record.id in updated_starts:
                    start = updated_starts[record.id]
                    stl.append("{0}\t{1}\t{2}\t{3}\t{4}".format(hap, record.id, 'END' if record.id in updated_orient else 'START',
                                                                start, '-' if record.id in updated_orient else '+'))

                if record.id in updated_ends:
                    end = updated_ends[record.id]
                    stl.append("{0}\t{1}\t{2}\t{3}\t{4}".format(hap, record.id, 'START' if record.id in updated_orient else 'END',
                                                                len(record.seq) - end, '-' if record.id in updated_orient else '+'))

                record.seq = record.seq[start:end]
                if record.id in updated_orient:
                    record.seq = record.seq.reverse_complement()
                trimmed_contigs.append(record)
        SeqIO.write(trimmed_contigs, contig_path, 'fasta')
    except subprocess.CalledProcessError as cpe:
        logging.error("Error aligning boundary sequences to contigs: {0} see stderr for details".format(cpe.cmd))
        sys.exit(1)
    finally:
        if os.path.exists(align_file):
            os.remove(align_file)


def stitch_contigs(contig_path, idx, end_distance, overlap_length, overlap_identity, minimap_path, threads):
    """
    Align contigs within each haplotype, if two ends overlap by overlap_length and at overlap_identity,
    stitch the two contigs together

    :param contig_path: Path to trimmed contigs
    :param idx: haplotype index
    :param end_distance: Maximum distance from contig end to allow overlap
    :param overlap_length: Minimum overlap between contig ends
    :param overlap_identity: Minimum identity within overlap to consider for stitching
    :param minimap_path: Path to minimap2
    :param threads: Number of threads for minimap2
    """
    out_path = Path(contig_path)
    align_file = out_path.parent.joinpath(out_path.stem + "_stitch.paf")
    hap = get_hap_name(idx)

    try:
        continue_merge = True
        iter = 1

        while continue_merge:
            logging.info("Contig stitching for {0} iteration {1}".format(hap, iter))
            continue_merge = False
            contig_dict = SeqIO.to_dict(SeqIO.parse(contig_path, "fasta"))
            subprocess.run([minimap_path, '-x', 'asm10', '-D', '-c', '-t', str(threads), '-o',
                            align_file, contig_path, contig_path], check=True)

            align_paf = pio.load_single_alignments(align_file)

            merged_seq = None
            for a in align_paf:
                if a.target_name != a.query_name and a.identityNM >= overlap_identity and a.align_length >= overlap_length and a.target_name:
                    if abs(a.query_end - a.query_length) < end_distance and a.target_start < end_distance:
                        merged_seq = str(contig_dict[a.query_name].seq)[:a.query_end] + str(
                            contig_dict[a.target_name].seq)[a.target_end:]
                    elif abs(a.query_end - a.query_length) < end_distance and abs(
                            a.target_end - a.target_length) < end_distance:
                        merged_seq = str(contig_dict[a.query_name].seq)[:a.query_end] + str(
                            contig_dict[a.target_name].seq[:a.target_start].reverse_complement())
                    elif a.query_start < end_distance and a.target_start < end_distance:
                        merged_seq = str(contig_dict[a.target_name].seq[a.target_end:].reverse_complement()) + str(
                            contig_dict[a.query_name].seq)[a.query_start:]
                    elif a.query_start < end_distance and abs(a.target_end - a.target_length) < end_distance:
                        merged_seq = str(contig_dict[a.target_name].seq)[:a.target_start] + str(
                            contig_dict[a.query_name].seq)[a.query_start:]

                if merged_seq is not None:
                    continue_merge = True
                    iter += 1
                    logging.info("Found overlap between {0} and {1}".format(a.target_name, a.query_name))
                    stl.append("{0}\t{1}\t{2}\t{3}\t{4:.4f}".format(hap, a.query_name, a.target_name, a.align_length,a.identityNM))
                    record = SeqRecord(Seq(merged_seq), id='merged_{0}_{1}'.format(a.query_name, a.target_name), description="")
                    updated_contigs = [record]

                    for contig, record in contig_dict.items():
                        if contig != a.query_name and contig != a.target_name:
                            updated_contigs.append(record)
                    SeqIO.write(updated_contigs, contig_path, "fasta")
                    break

    except subprocess.CalledProcessError as cpe:
        logging.error("Error overlapping assemblies: {0}".format(cpe.cmd))
        sys.exit(1)
    finally:
        if os.path.exists(align_file):
            os.remove(align_file)


def get_hap_name(idx):
    return "hap1" if idx == 0 else "hap2"


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(1)
