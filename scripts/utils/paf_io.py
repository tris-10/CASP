from cigar import Cigar

"""
Classes and methods used to read, store and manipulate alignment information generated by minimap2 and
written in PAF format
"""


class PAFAlignment:
    """Simple container for alignment information"""

    def __init__(self, align):
        """
        Takes a string representing an alignment in PAF format. Each PAF field is stored by type in the object.

        :param align: String representing one alignment in PAF format.
        """

        align_parts = align.strip().split("\t")
        self.query_name = align_parts[0]
        self.query_length = int(align_parts[1])
        self.query_start = int(align_parts[2])
        self.query_end = int(align_parts[3])
        self.strand = align_parts[4]
        self.target_name = align_parts[5]
        self.target_length = int(align_parts[6])
        self.target_start = int(align_parts[7])
        self.target_end = int(align_parts[8])
        self.match_count = int(align_parts[9])
        self.align_length = int(align_parts[10])
        self.map_qual = int(align_parts[11])
        self.identityNM = None
        self.totalNM = None
        self.cigar = None

        for tag in align_parts[12:]:
            key, type, value = tag.split(":")
            if key == "NM":
                self.identityNM = (self.align_length - float(value)) / self.align_length
                self.totalNM = float(value)
            elif key == "cg":
                c = Cigar(value)
                self.cigar = c

        self.query_overlap = self.align_length / self.query_length


def load_query_alignments(paf_file):
    """
    Reads through an alignment file in PAF format. Returns all alignments for each query sequence. Assumes
    the PAF file isn't sorted differently than the default.

    :param paf_file: Path to alignment file in PAF format.
    :return: list of PAFAlignment objects for each query sequence
    """

    last_query_name = None
    query_align_list = []
    with open(paf_file, 'r') as paf_in:
        for line in paf_in:
            pa = PAFAlignment(line)
            if pa.query_name != last_query_name:
                if last_query_name is not None:
                    yield query_align_list
                query_align_list = []
                last_query_name = pa.query_name
            query_align_list.append(pa)
    yield query_align_list


def load_single_alignments(paf_file):
    """
     Reads through an alignment file in PAF format. Returns all alignments for each query sequence. Assumes
     the PAF file isn't sorted differently than the default.

     :param paf_file: Path to alignment file in PAF format.
     :return: list of PAFAlignment objects for each query sequence
     """
    
    with open(paf_file, 'r') as paf_in:
        for line in paf_in:
            pa = PAFAlignment(line)
            yield pa
