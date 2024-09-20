#!/usr/bin/env python

"""
Find overlaps between contigs in a canu assembly.  Remove contigs with a high degree of overlap and identity with one or
more longer contigs.  Trim contigs ends to remove redundant sequences.
"""

from utils import paf_io as pio, contig_io as cio
from pathlib import Path
import sys
import argparse
import subprocess
import logging

info_list = []


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_contig_path', help="Path to canu assembly in fasta format.")
    parser.add_argument('output_contig_path', help="Path to filtered/trimmed contig fasta output file.")
    parser.add_argument('output_stats', help="Path to RCCX replacement statistics output file.")
    parser.add_argument('-c', '--min_contig_length', type=int, default=20000, help="Minimum contig length to retain.")
    parser.add_argument('-t', '--threads', type=int, default=4, help="Number of minimap2 threads.")
    parser.add_argument("-n", '--minimap2_path', default="minimap2", help="Path to minimap2 executable.")
    g1 = parser.add_argument_group("Contig Overlap Filtration")
    g1.add_argument('-i', '--min_identity', type=float, default=0.99,
                    help="Minimum overlap identity required for filtration.")
    g1.add_argument('-p', '--min_overlap', type=float, default=0.99,
                    help="Minimum fraction overlap with a single larger contigs required for filtration.")
    g1.add_argument('-l', '--max_align_length', type=int, default=200000,
                    help="Maximum alignment length allowed for filtration. Avoids removing extremely long contigs.")
    g2 = parser.add_argument_group("Contig Trimming")
    g2.add_argument('-m', '--max_diff', type=float, default=1e-3,
                 help="Max allowed difference between contigs to consider identical")

    args = parser.parse_args()

    logging.basicConfig(format='%(levelname)s:%(asctime)s %(message)s', level=logging.INFO)

    logging.info("Loading contig information")
    tig_dict = cio.create_tiginfo_dict(args.input_contig_path)

    logging.info("Finding contig overlaps")
    align_file = align_contigs_to_self(args.input_contig_path, args.output_contig_path, args.minimap2_path, args.threads)

    logging.info("Filter contigs by overlap")
    filter_by_overlap(align_file, tig_dict, args.min_identity, args.max_align_length, args.min_overlap)

    logging.info("Trim contigs")
    trim_contig_ends(align_file, tig_dict, args.max_diff)

    logging.info("Writing updated contig list")
    cio.write_trimmed_contigs(args.output_contig_path, tig_dict, args.min_contig_length)

    logging.info("Writing filter/trimming information")
    with open(args.output_stats, "w") as op_stats:
        op_stats.write("\n".join(info_list) + "\n")


def align_contigs_to_self(contig_file, output_path, minimap_path, threads):
    """
    Align contigs to themselves using minimap2

    :param contig_file: Path to the location/read count filtered contig file
    :param output_path: Path to output alignment file
    :param minimap_path: Path the minimap2 executable
    :param threads: Number of threads to use for alignment
    :return: Path to the alignment file in PAF format.
    """

    align_file = Path(output_path).parent.joinpath("temp_self.paf")

    try:
        subprocess.run([minimap_path, '-x', 'asm20', '-p', '0.20', '-D', '-c', '--eqx', '-t', str(threads), '-o',
                        align_file, contig_file, contig_file], check=True)
    except subprocess.CalledProcessError as error:
        logging.error("minimap2 self-alignment failed: %s:", error.cmd)
        sys.exit(1)

    return align_file


def trim_contig_ends(self_align, tig_dict, max_diff):
    """
    Check for long regions of homology between contigs and the edges of the shorter contig. Homology is determined
    by scoring +1 for matches and mm_penalty for mismatches. Find the location with the best positive score and
    trim the contig if it at least grace_bp away from the contig end.

    :param self_align: All-to-all contig alignment file in PAF format
    :param tig_dict: Dictionary of all contig objects reported by Canu
    :param grace_bp: Padding to leave upstream of the discordant region when trimming
    """

    info_list.append("\nContig Trimming\nTrimTig\tTrimLoc\tTrimLength\tTigLength\tOvlTig")
    self_paf = pio.load_query_alignments(self_align)
    for align_list in self_paf:
        for a in align_list:
            if a.query_length < a.target_length and not tig_dict[a.query_name].filtered and not tig_dict[a.target_name].filtered:
                if a.query_start < 100:
                    best_loc = process_cigar(a, True, max_diff)
                    if best_loc > 0:
                        info_list.append("{0}\ttrim_start\t{1}\t{2}\t{3}".format(a.query_name, best_loc,
                                                                                 a.query_length, a.target_name))
                        tig_dict[a.query_name].trim_start(best_loc)

                if (a.query_length - a.query_end) < 100:
                    best_loc = process_cigar(a, False, max_diff)
                    if best_loc > 0:
                        info_list.append("{0}\ttrim_end\t{1}\t{2}\t{3}".format(a.query_name, best_loc,
                                                                                 a.query_length, a.target_name))
                        tig_dict[a.query_name].trim_end(best_loc)


def filter_by_overlap(self_align, tig_dict, min_identity, max_align_length, min_overlap):
    """
    Check contigs for high overlap and identity with longer contigs.  These contigs are set to be filtered
    and not written to output.

    :param self_align: All-to-all contig alignment file in PAF format
    :param tig_dict: Dictionary of all contig objects reported by Canu
    :param min_identity: Minimum identity allowed to retain alignment
    :param max_align_length: Maximum overlap length allowed to retain alignment
    :param min_overlap: Minimum allowed fraction overlap to a single longer contig
    :param min_combined_overlap: Minimum allowed fraction overlap to all longer contigs
    """

    info_list.append("Single Contig Overlap\nFiltTig\tStatus\tOverlap\tIdentity\tOvlTig")

    # Remove contigs with high overlap to other contigs
    self_paf = pio.load_query_alignments(self_align)
    for align_list in self_paf:
        for a in align_list:
            if a.identityNM >= min_identity and a.query_length < a.target_length and a.query_length < max_align_length:
                if a.query_overlap >= min_overlap:
                    tig_dict[a.query_name].filter_contig()
                    info_list.append("{0}\toverlap_single\t{1:.4f}\t{2:.4f}\t{3}".format(a.query_name, a.query_overlap,
                                                                           a.identityNM, a.target_name))


def process_cigar(a, start, max_diff):
    """
    Scan the cigar string from the edge of the contig inwards.  Score matches + 1 and mismatches mm_penalty.  Find
    the location with the max non-zero score.

    :param a: alignment info
    :param mm_penalty: penalty for a mismatch
    :param start: True if overlap is at the start of the contig, False if the overlap is at the end
    :param grace_bp: Number of bases to skip prior to scoring trimming
    :return: Number of base pairs to remove from the end of the contig.
    """

    if a.strand == "+":
        if start:
            cigar = list(a.cigar.items())
            curr_loc = a.query_start
        else:
            cigar = list(a.cigar.items())[::-1]
            curr_loc = (a.query_length - a.query_end)
    else:
        if start:
            cigar = list(a.cigar.items())[::-1]
            curr_loc = a.query_start
        else:
            cigar = list(a.cigar.items())
            curr_loc = (a.query_length - a.query_end)

    loc = []
    score = []
    opl = []
    curr_score = -1

    for op in cigar:
        if op[1] == "=":
            curr_loc += op[0]
            curr_score += max_diff * op[0]
        elif op[1] == "X" or op[1] == "I":
            curr_loc += op[0]
            curr_score += -1 * op[0]
            score.append(curr_score)
            loc.append(curr_loc)
            opl.append("{0}:{1}".format(op[1], op[0]))
        elif op[1] == "D":
            curr_score += -1 * op[0]
            score.append(curr_score)
            loc.append(curr_loc)
            opl.append("{0}:{1}".format(op[1], op[0]))
        else:
            continue

    best_score = best_loc = 0
    for i in range(len(score)):
        if score[i] > best_score:
            best_score = score[i]
            best_loc = loc[i]
    return best_loc - 1000


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(1)
