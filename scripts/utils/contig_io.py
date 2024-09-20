"""
Container for contig data
"""

from Bio import SeqIO
import pysam


class TigInfo:
    """
    Container for contig information and used to track various filtering metrics.
    """

    def __init__(self, tig_name, tig_sequence):
        """
        Initialize TigInfo object.

        :param tig_name: Contig name, which is the ID padded with leading zeros out to 8 positions
        :param tig_sequence: Contig sequence as SeqRecord object
        """
        self.tig_name = tig_name
        self.tig_sequence = tig_sequence
        self.tig_length = len(tig_sequence)
        self.updated_start = 0
        self.updated_end = self.tig_length
        self.filtered = False

    def filter_contig(self):
        """
        Set filtered to True
        """
        self.filtered = True

    def trim_start(self, pred_start):
        """
        Adjust the start of the contig based on observed self-overlaps.

        :param pred_start: Number of base pairs from the start of the contig to trim.
        """
        if pred_start > self.updated_start:
            self.updated_start = pred_start

    def trim_end(self, pred_end):
        """
        Adjust the end of the contig based on observed self-overlaps

        :param pred_end: Number of base pairs from the end of the contig to trim.
        """
        updated_end = self.tig_length - pred_end
        if updated_end < self.updated_end:
            self.updated_end = updated_end

    def get_updated_sequence(self):
        """
        Adjust a SeqRecord object using updated start/end coordinates

        :param record: Original contig sequence
        :return: Trimmed contig sequence
        """
        return self.tig_sequence[self.updated_start:self.updated_end]


def create_tiginfo_dict(tig_fasta):
    """
    Load contig data into a dictionary of tigInfo objects

    :param tig_fasta: Path to contig fasta
    :return: dictionary of TigInfo objects, keyed on contig name
    """
    contig_dict = {}
    with open(tig_fasta, 'r') as ip_fasta:
        for contig in SeqIO.parse(ip_fasta, 'fasta'):
            contig_dict[contig.id] = TigInfo(contig.id, contig)
    return contig_dict


def write_trimmed_contigs(trimmed_assembly, tig_dict, min_length):
    """
    Write out trimmed sequences to file

    :param trimmed_assembly: Path to trimmed contig objects
    :param tig_dict: Dictionary of TigInfo objects
    :param min_length: Minimum allow contig length
    """
    with open(trimmed_assembly, 'w') as fasta_out:
        for tig in tig_dict.values():
            if tig.filtered or len(tig.get_updated_sequence()) <= min_length:
                continue
            SeqIO.write(tig.get_updated_sequence(), fasta_out, 'fasta')
    pysam.faidx(trimmed_assembly)
