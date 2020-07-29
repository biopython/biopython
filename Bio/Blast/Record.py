# Copyright 1999-2000 by Jeffrey Chang.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Record classes to hold BLAST output.

Classes:
Blast              Holds all the information from a blast search.
PSIBlast           Holds all the information from a psi-blast search.

Header             Holds information from the header.
Description        Holds information about one hit description.
Alignment          Holds information about one alignment hit.
HSP                Holds information about one HSP.
MultipleAlignment  Holds information about a multiple alignment.
DatabaseReport     Holds information from the database report.
Parameters         Holds information from the parameters.

"""
# XXX finish printable BLAST output

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment


def fmt_(value, format_spec="%s", default_str="<unknown>"):
    """Ensure the given value formats to a string correctly."""
    if value is None:
        return default_str
    return format_spec % value


class Header:
    """Saves information from a blast header.

    Members:
    application         The name of the BLAST flavor that generated this data.
    version             Version of blast used.
    date                Date this data was generated.
    reference           Reference for blast.

    query               Name of query sequence.
    query_letters       Number of letters in the query sequence.  (int)

    database            Name of the database.
    database_sequences  Number of sequences in the database.  (int)
    database_letters    Number of letters in the database.  (int)

    """

    def __init__(self):
        """Initialize the class."""
        self.application = ""
        self.version = ""
        self.date = ""
        self.reference = ""

        self.query = ""
        self.query_letters = None

        self.database = ""
        self.database_sequences = None
        self.database_letters = None


class Description:
    """Stores information about one hit in the descriptions section.

    Members:
    title           Title of the hit.
    score           Number of bits.  (int)
    bits            Bit score. (float)
    e               E value.  (float)
    num_alignments  Number of alignments for the same subject.  (int)
    """

    def __init__(self):
        """Initialize the class."""
        self.title = ""
        self.score = None
        self.bits = None
        self.e = None
        self.num_alignments = None

    def __str__(self):
        """Return the description as a string."""
        return "%-66s %5s  %s" % (self.title, self.score, self.e)


class DescriptionExt(Description):
    """Extended description record for BLASTXML version 2.

    Members:
    items           List of DescriptionExtItem
    """

    def __init__(self):
        """Initialize the class."""
        super().__init__()

        self.items = []

    def append_item(self, item):
        """Add a description extended record."""
        if len(self.items) == 0:
            self.title = str(item)
        self.items.append(item)


class DescriptionExtItem:
    """Stores information about one record in hit description for BLASTXML version 2.

    Members:
    id              Database identifier
    title           Title of the hit.
    """

    def __init__(self):
        """Initialize the class."""
        self.id = None
        self.title = None
        self.accession = None
        self.taxid = None
        self.sciname = None

    def __str__(self):
        """Return the description identifier and title as a string."""
        return "%s %s" % (self.id, self.title)


class Alignment:
    """Stores information about one hit in the alignments section.

    Members:
    title      Name.
    hit_id     Hit identifier. (str)
    hit_def    Hit definition. (str)
    length     Length.  (int)
    hsps       A list of HSP objects.

    """

    def __init__(self):
        """Initialize the class."""
        self.title = ""
        self.hit_id = ""
        self.hit_def = ""
        self.length = None
        self.hsps = []

    def __str__(self):
        """Return the BLAST alignment as a formatted string."""
        lines = self.title.split("\n")
        lines.append("Length = %s\n" % self.length)
        return "\n           ".join(lines)


class HSP:
    """Stores information about one hsp in an alignment hit.

    Members:
        - score           BLAST score of hit.  (float)
        - bits            Number of bits for that score.  (float)
        - expect          Expect value.  (float)
        - num_alignments  Number of alignments for same subject.  (int)
        - identities      Number of identities (int) if using the XML parser.
          Tuple of number of identities/total aligned (int, int)
          if using the (obsolete) plain text parser.
        - positives       Number of positives (int) if using the XML parser.
          Tuple of number of positives/total aligned (int, int)
          if using the (obsolete) plain text parser.
        - gaps            Number of gaps (int) if using the XML parser.
          Tuple of number of gaps/total aligned (int, int) if
          using the (obsolete) plain text parser.
        - align_length    Length of the alignment. (int)
        - strand          Tuple of (query, target) strand.
        - frame           Tuple of 1 or 2 frame shifts, depending on the flavor.

        - query           The query sequence.
        - query_start     The start residue for the query sequence.  (1-based)
        - query_end       The end residue for the query sequence.  (1-based)
        - match           The match sequence.
        - sbjct           The sbjct sequence.
        - sbjct_start     The start residue for the sbjct sequence.  (1-based)
        - sbjct_end       The end residue for the sbjct sequence.  (1-based)

    Not all flavors of BLAST return values for every attribute::

                  score     expect     identities   positives    strand  frame
        BLASTP     X          X            X            X
        BLASTN     X          X            X            X          X
        BLASTX     X          X            X            X                  X
        TBLASTN    X          X            X            X                  X
        TBLASTX    X          X            X            X                 X/X

    Note: for BLASTX, the query sequence is shown as a protein sequence,
    but the numbering is based on the nucleotides.  Thus, the numbering
    is 3x larger than the number of amino acid residues.  A similar effect
    can be seen for the sbjct sequence in TBLASTN, and for both sequences
    in TBLASTX.

    Also, for negative frames, the sequence numbering starts from
    query_start and counts down.

    """

    def __init__(self):
        """Initialize the class."""
        self.score = None
        self.bits = None
        self.expect = None
        self.num_alignments = None
        self.identities = (None, None)
        self.positives = (None, None)
        self.gaps = (None, None)
        self.align_length = None
        self.strand = (None, None)
        self.frame = ()

        self.query = ""
        self.query_start = None
        self.query_end = None
        self.match = ""
        self.sbjct = ""
        self.sbjct_start = None
        self.sbjct_end = None

    def __str__(self):
        """Return the BLAST HSP as a formatted string."""
        lines = [
            "Score %s (%s bits), expectation %s, alignment length %s"
            % (
                fmt_(self.score, "%i"),
                fmt_(self.bits, "%i"),
                fmt_(self.expect, "%0.1e"),
                fmt_(self.align_length, "%i"),
            )
        ]
        if self.align_length is None:
            return "\n".join(lines)
        if self.align_length < 50:
            lines.append(
                "Query:%s %s %s"
                % (str(self.query_start).rjust(8), str(self.query), str(self.query_end))
            )
            lines.append("               %s" % (str(self.match)))
            lines.append(
                "Sbjct:%s %s %s"
                % (str(self.sbjct_start).rjust(8), str(self.sbjct), str(self.sbjct_end))
            )
        else:
            lines.append(
                "Query:%s %s...%s %s"
                % (
                    str(self.query_start).rjust(8),
                    str(self.query)[:45],
                    str(self.query)[-3:],
                    str(self.query_end),
                )
            )
            lines.append(
                "               %s...%s" % (str(self.match)[:45], str(self.match)[-3:])
            )
            lines.append(
                "Sbjct:%s %s...%s %s"
                % (
                    str(self.sbjct_start).rjust(8),
                    str(self.sbjct)[:45],
                    str(self.sbjct)[-3:],
                    str(self.sbjct_end),
                )
            )
        return "\n".join(lines)


class MultipleAlignment:
    """Holds information about a multiple alignment.

    Members:
    alignment  A list of tuples (name, start residue, sequence, end residue).

    The start residue is 1-based.  It may be blank, if that sequence is
    not aligned in the multiple alignment.

    """

    def __init__(self):
        """Initialize the class."""
        self.alignment = []

    def to_generic(self):
        """Retrieve generic alignment object for the given alignment.

        Instead of the tuples, this returns a MultipleSeqAlignment object
        from Bio.Align, through which you can manipulate and query
        the object.

        Thanks to James Casbon for the code.
        """
        # TODO - Switch to new Bio.Align.MultipleSeqAlignment class?
        seq_parts = []
        seq_names = []
        parse_number = 0
        n = 0
        for name, start, seq, end in self.alignment:
            if name == "QUERY":  # QUERY is the first in each alignment block
                parse_number += 1
                n = 0

            if parse_number == 1:  # create on first_parse, append on all others
                seq_parts.append(seq)
                seq_names.append(name)
            else:
                seq_parts[n] += seq
                n += 1

        generic = MultipleSeqAlignment([])
        for (name, seq) in zip(seq_names, seq_parts):
            generic.append(SeqRecord(Seq(seq), name))

        return generic


class Round:
    """Holds information from a PSI-BLAST round.

    Members:
    number       Round number.  (int)
    reused_seqs  Sequences in model, found again.  List of Description objects.
    new_seqs     Sequences not found, or below threshold.  List of Description.
    alignments          A list of Alignment objects.
    multiple_alignment  A MultipleAlignment object.
    """

    def __init__(self):
        """Initialize the class."""
        self.number = None
        self.reused_seqs = []
        self.new_seqs = []
        self.alignments = []
        self.multiple_alignment = None


class DatabaseReport:
    """Holds information about a database report.

    Members:
    database_name              List of database names.  (can have multiple dbs)
    num_letters_in_database    Number of letters in the database.  (int)
    num_sequences_in_database  List of number of sequences in the database.
    posted_date                List of the dates the databases were posted.
    ka_params                  A tuple of (lambda, k, h) values.  (floats)
    gapped                     # XXX this isn't set right!
    ka_params_gap              A tuple of (lambda, k, h) values.  (floats)

    """

    def __init__(self):
        """Initialize the class."""
        self.database_name = []
        self.posted_date = []
        self.num_letters_in_database = []
        self.num_sequences_in_database = []
        self.ka_params = (None, None, None)
        self.gapped = 0
        self.ka_params_gap = (None, None, None)


class Parameters:
    """Holds information about the parameters.

    Members:
    matrix              Name of the matrix.
    gap_penalties       Tuple of (open, extend) penalties.  (floats)
    sc_match            Match score for nucleotide-nucleotide comparison
    sc_mismatch         Mismatch penalty for nucleotide-nucleotide comparison
    num_hits            Number of hits to the database.  (int)
    num_sequences       Number of sequences.  (int)
    num_good_extends    Number of extensions.  (int)
    num_seqs_better_e   Number of sequences better than e-value.  (int)
    hsps_no_gap         Number of HSP's better, without gapping.  (int)
    hsps_prelim_gapped  Number of HSP's gapped in prelim test.  (int)
    hsps_prelim_gapped_attemped  Number of HSP's attempted in prelim.  (int)
    hsps_gapped         Total number of HSP's gapped.  (int)
    query_length        Length of the query.  (int)
    query_id            Identifier of the query sequence. (str)
    database_length     Number of letters in the database.  (int)
    effective_hsp_length         Effective HSP length.  (int)
    effective_query_length       Effective length of query.  (int)
    effective_database_length    Effective length of database.  (int)
    effective_search_space       Effective search space.  (int)
    effective_search_space_used  Effective search space used.  (int)
    frameshift          Frameshift window.  Tuple of (int, float)
    threshold           Threshold.  (int)
    window_size         Window size.  (int)
    dropoff_1st_pass    Tuple of (score, bits).  (int, float)
    gap_x_dropoff       Tuple of (score, bits).  (int, float)
    gap_x_dropoff_final Tuple of (score, bits).  (int, float)
    gap_trigger         Tuple of (score, bits).  (int, float)
    blast_cutoff        Tuple of (score, bits).  (int, float)
    """

    def __init__(self):
        """Initialize the class."""
        self.matrix = ""
        self.gap_penalties = (None, None)
        self.sc_match = None
        self.sc_mismatch = None
        self.num_hits = None
        self.num_sequences = None
        self.num_good_extends = None
        self.num_seqs_better_e = None
        self.hsps_no_gap = None
        self.hsps_prelim_gapped = None
        self.hsps_prelim_gapped_attemped = None
        self.hsps_gapped = None
        self.query_id = None
        self.query_length = None
        self.database_length = None
        self.effective_hsp_length = None
        self.effective_query_length = None
        self.effective_database_length = None
        self.effective_search_space = None
        self.effective_search_space_used = None
        self.frameshift = (None, None)
        self.threshold = None
        self.window_size = None
        self.dropoff_1st_pass = (None, None)
        self.gap_x_dropoff = (None, None)
        self.gap_x_dropoff_final = (None, None)
        self.gap_trigger = (None, None)
        self.blast_cutoff = (None, None)


# TODO - Add a friendly __str__ method to BLAST results
class Blast(Header, DatabaseReport, Parameters):
    """Saves the results from a blast search.

    Members:
    descriptions        A list of Description objects.
    alignments          A list of Alignment objects.
    multiple_alignment  A MultipleAlignment object.
    + members inherited from base classes

    """

    def __init__(self):
        """Initialize the class."""
        Header.__init__(self)
        DatabaseReport.__init__(self)
        Parameters.__init__(self)
        self.descriptions = []
        self.alignments = []
        self.multiple_alignment = None


class PSIBlast(Header, DatabaseReport, Parameters):
    """Saves the results from a blastpgp search.

    Members:
    rounds       A list of Round objects.
    converged    Whether the search converged.
    + members inherited from base classes

    """

    def __init__(self):
        """Initialize the class."""
        Header.__init__(self)
        DatabaseReport.__init__(self)
        Parameters.__init__(self)
        self.rounds = []
        self.converged = 0
