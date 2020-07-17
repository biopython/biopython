# Copyright 2000 by Bertrand Frottier.  All rights reserved.
# Revisions 2005-2006 copyright Michiel de Hoon
# Revisions 2006-2009 copyright Peter Cock
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Code to work with the BLAST XML output.

The BLAST XML DTD file is on the NCBI FTP site at:
ftp://ftp.ncbi.nlm.nih.gov/blast/documents/xml/NCBI_BlastOutput.dtd
"""

from Bio.Blast import Record
import xml.sax
from xml.sax.handler import ContentHandler


class _XMLparser(ContentHandler):
    """Generic SAX Parser (PRIVATE).

    Just a very basic SAX parser.

    Redefine the methods startElement, characters and endElement.
    """

    def __init__(self, debug=0):
        """Initialize the parser.

        Arguments:
         - debug - integer, amount of debug information to print

        """
        self._tag = []
        self._value = ""
        self._debug = debug
        self._debug_ignore_list = []
        self._method_name_level = 1
        self._method_map = None

    def startElement(self, name, attr):
        """Found XML start tag.

        No real need of attr, BLAST DTD doesn't use them

        Arguments:
         - name -- name of the tag
         - attr -- tag attributes

        """
        self._tag.append(name)

        if len(self._tag) == 1:
            # root node
            self._on_root_node(name)
            return

        # Try to call a method (defined in subclasses)
        method = "start_" + self._node_method_name(name)

        # Note could use try / except AttributeError
        # BUT I found often triggered by nested errors...
        if method in self._method_map:
            self._method_map[method]()
            if self._debug > 4:
                print("NCBIXML: Parsed:  " + method)
        elif self._debug > 3:
            # Doesn't exist (yet) and may want to warn about it
            if method not in self._debug_ignore_list:
                print("NCBIXML: Ignored: " + method)
                self._debug_ignore_list.append(method)

        # We don't care about white space in parent tags like Hsp,
        # but that white space doesn't belong to child tags like Hsp_midline
        if self._value.strip():
            raise ValueError(
                "What should we do with %s before the %r tag?" % (self._value, name)
            )
        self._value = ""

    def characters(self, ch):
        """Found some text.

        Arguments:
         - ch -- characters read

        """
        self._value += ch  # You don't ever get the whole string

    def endElement(self, name):
        """Found XML end tag.

        Arguments:
         - name -- tag name

        """
        # DON'T strip any white space, we may need it e.g. the hsp-midline

        # Try to call a method (defined in subclasses)
        method = "end_" + self._node_method_name(name)

        # Note could use try / except AttributeError
        # BUT I found often triggered by nested errors...
        if method in self._method_map:
            self._method_map[method]()
            if self._debug > 2:
                print("NCBIXML: Parsed:  %s %s" % (method, self._value))
        elif self._debug > 1:
            # Doesn't exist (yet) and may want to warn about it
            if method not in self._debug_ignore_list:
                print("NCBIXML: Ignored: %s %s" % (method, self._value))
                self._debug_ignore_list.append(method)

        # Reset character buffer
        self._value = ""

        self._tag.pop()

    def _node_method_name(self, name):
        if self._method_name_level == 1:
            return name
        return "/".join(self._tag[-self._method_name_level :])


class BlastParser(_XMLparser):
    """Parse XML BLAST data into a Record.Blast object.

    Parses XML output from BLAST (direct use discouraged).
    This (now) returns a list of Blast records.
    Historically it returned a single Blast record.
    You are expected to use this via the parse or read functions.

    All XML 'action' methods are private methods and may be:

    - ``_start_TAG`` called when the start tag is found
    - ``_end_TAG`` called when the end tag is found

    """

    def __init__(self, debug=0):
        """Initialize the parser.

        Arguments:
         - debug - integer, amount of debug information to print

        """
        # Calling superclass method
        _XMLparser.__init__(self, debug)

        self._parser = xml.sax.make_parser()
        self._parser.setContentHandler(self)

        # To avoid ValueError: unknown url type: NCBI_BlastOutput.dtd
        self._parser.setFeature(xml.sax.handler.feature_validation, 0)
        self._parser.setFeature(xml.sax.handler.feature_namespaces, 0)
        self._parser.setFeature(xml.sax.handler.feature_external_pes, 0)
        self._parser.setFeature(xml.sax.handler.feature_external_ges, 0)

        self._xml_version = 1

        self.reset()

    def reset(self):
        """Reset all the data allowing reuse of the BlastParser() object."""
        self._records = []
        self._header = Record.Header()
        self._parameters = Record.Parameters()
        self._parameters.filter = None  # Maybe I should update the class?

    def _on_root_node(self, name):
        if name == "BlastOutput":
            self._setup_blast_v1()
        elif name == "BlastXML2":
            self._setup_blast_v2()
        else:
            raise ValueError(
                "Invalid root node name: %s. Root node should be either"
                " BlastOutput or BlastXML2" % name
            )

    def _setup_blast_v1(self):
        self._method_map = {
            "start_Iteration": self._start_blast_record,
            "end_Iteration": self._end_blast_record,
            "end_BlastOutput_program": self._set_header_application,
            "end_BlastOutput_version": self._set_header_version,
            "end_BlastOutput_reference": self._set_header_reference,
            "end_BlastOutput_db": self._set_header_database,
            "end_BlastOutput_query-ID": self._set_header_query_id,
            "end_BlastOutput_query-def": self._set_header_query,
            "end_BlastOutput_query-len": self._set_header_query_letters,
            "end_Iteration_query-ID": self._set_record_query_id,
            "end_Iteration_query-def": self._set_record_query_def,
            "end_Iteration_query-len": self._set_record_query_letters,
            "end_BlastOutput_hits": self._set_record_hits,
            "end_Parameters_matrix": self._set_parameters_matrix,
            "end_Parameters_expect": self._set_parameters_expect,
            "end_Parameters_sc-match": self._set_parameters_sc_match,
            "end_Parameters_sc-mismatch": self._set_parameters_sc_mismatch,
            "end_Parameters_gap-open": self._set_parameters_gap_penalties,
            "end_Parameters_gap-extend": self._set_parameters_gap_extend,
            "end_Parameters_filter": self._set_parameters_filter,
            "start_Hit": self._start_hit,
            "end_Hit": self._end_hit,
            "end_Hit_id": self.set_hit_id,
            "end_Hit_def": self.set_hit_def,
            "end_Hit_accession": self.set_hit_accession,
            "end_Hit_len": self.set_hit_len,
            "start_Hsp": self._start_hsp,
            "end_Hsp_score": self._set_hsp_score,
            "end_Hsp_bit-score": self._set_hsp_bit_score,
            "end_Hsp_evalue": self._set_hsp_e_value,
            "end_Hsp_query-from": self._set_hsp_query_start,
            "end_Hsp_query-to": self._set_hsp_query_end,
            "end_Hsp_hit-from": self._set_hsp_hit_from,
            "end_Hsp_hit-to": self._set_hsp_hit_to,
            "end_Hsp_query-frame": self._set_hsp_query_frame,
            "end_Hsp_hit-frame": self._set_hsp_hit_frame,
            "end_Hsp_identity": self._set_hsp_identity,
            "end_Hsp_positive": self._set_hsp_positive,
            "end_Hsp_gaps": self._set_hsp_gaps,
            "end_Hsp_align-len": self._set_hsp_align_len,
            "end_Hsp_qseq": self._set_hsp_query_seq,
            "end_Hsp_hseq": self._set_hsp_subject_seq,
            "end_Hsp_midline": self._set_hsp_midline,
            "end_Statistics_db-num": self._set_statistics_db_num,
            "end_Statistics_db-len": self._set_statistics_db_len,
            "end_Statistics_hsp-len": self._set_statistics_hsp_len,
            "end_Statistics_eff-space": self._set_statistics_eff_space,
            "end_Statistics_kappa": self._set_statistics_kappa,
            "end_Statistics_lambda": self._set_statistics_lambda,
            "end_Statistics_entropy": self._set_statistics_entropy,
        }

    def _setup_blast_v2(self):
        self._method_name_level = 2
        self._xml_version = 2
        self._method_map = {
            "start_report/Report": self._start_blast_record,
            "end_report/Report": self._end_blast_record,
            "end_Report/program": self._set_header_application,
            "end_Report/version": self._set_header_version,
            "end_Report/reference": self._set_header_reference,
            "end_Target/db": self._set_header_database,
            "end_Search/query-id": self._set_record_query_id,
            "end_Search/query-title": self._set_record_query_def,
            "end_Search/query-len": self._set_record_query_letters,
            "end_BlastOutput_hits": self._set_record_hits,
            "end_Parameters/matrix": self._set_parameters_matrix,
            "end_Parameters/expect": self._set_parameters_expect,
            "end_Parameters/sc-match": self._set_parameters_sc_match,
            "end_Parameters/sc-mismatch": self._set_parameters_sc_mismatch,
            "end_Parameters/gap-open": self._set_parameters_gap_penalties,
            "end_Parameters/gap-extend": self._set_parameters_gap_extend,
            "end_Parameters/filter": self._set_parameters_filter,
            "start_hits/Hit": self._start_hit,
            "end_hits/Hit": self._end_hit,
            "start_description/HitDescr": self._start_hit_descr_item,
            "end_description/HitDescr": self._end_hit_descr_item,
            "end_HitDescr/id": self._end_description_id,
            "end_HitDescr/accession": self._end_description_accession,
            "end_HitDescr/title": self._end_description_title,
            "end_HitDescr/taxid": self._end_description_taxid,
            "end_HitDescr/sciname": self._end_description_sciname,
            "end_Hit/len": self.set_hit_len,
            "start_hsps/Hsp": self._start_hsp,
            "end_hsps/Hsp": self._end_hsp,
            "end_Hsp/score": self._set_hsp_score,
            "end_Hsp/bit-score": self._set_hsp_bit_score,
            "end_Hsp/evalue": self._set_hsp_e_value,
            "end_Hsp/query-from": self._set_hsp_query_start,
            "end_Hsp/query-to": self._set_hsp_query_end,
            "end_Hsp/hit-from": self._set_hsp_hit_from,
            "end_Hsp/hit-to": self._set_hsp_hit_to,
            "end_Hsp/query-frame": self._set_hsp_query_frame,
            "end_Hsp/hit-frame": self._set_hsp_hit_frame,
            "end_Hsp/query-strand": self._set_hsp_query_strand,
            "end_Hsp/hit-strand": self._set_hsp_hit_strand,
            "end_Hsp/identity": self._set_hsp_identity,
            "end_Hsp/positive": self._set_hsp_positive,
            "end_Hsp/gaps": self._set_hsp_gaps,
            "end_Hsp/align-len": self._set_hsp_align_len,
            "end_Hsp/qseq": self._set_hsp_query_seq,
            "end_Hsp/hseq": self._set_hsp_subject_seq,
            "end_Hsp/midline": self._set_hsp_midline,
            "end_Statistics/db-num": self._set_statistics_db_num,
            "end_Statistics/db-len": self._set_statistics_db_len,
            "end_Statistics/hsp-len": self._set_statistics_hsp_len,
            "end_Statistics/eff-space": self._set_statistics_eff_space,
            "end_Statistics/kappa": self._set_statistics_kappa,
            "end_Statistics/lambda": self._set_statistics_lambda,
            "end_Statistics/entropy": self._set_statistics_entropy,
        }

    def _start_blast_record(self):
        """Start interaction (PRIVATE)."""
        self._blast = Record.Blast()
        pass

    def _end_blast_record(self):
        """End interaction (PRIVATE)."""
        # We stored a lot of generic "top level" information
        # in self._header (an object of type Record.Header)
        self._blast.reference = self._header.reference
        self._blast.date = self._header.date
        self._blast.version = self._header.version
        self._blast.database = self._header.database
        self._blast.application = self._header.application

        # These are required for "old" pre 2.2.14 files
        # where only <BlastOutput_query-ID>, <BlastOutput_query-def>
        # and <BlastOutput_query-len> were used.  Now they
        # are supplemented/replaced by <Iteration_query-ID>,
        # <Iteration_query-def> and <Iteration_query-len>
        if not hasattr(self._blast, "query") or not self._blast.query:
            self._blast.query = self._header.query
        if not hasattr(self._blast, "query_id") or not self._blast.query_id:
            self._blast.query_id = self._header.query_id
        if not hasattr(self._blast, "query_letters") or not self._blast.query_letters:
            self._blast.query_letters = self._header.query_letters

        # Hack to record the query length as both the query_letters and
        # query_length properties (as in the plain text parser, see
        # Bug 2176 comment 12):
        self._blast.query_length = self._blast.query_letters
        # Perhaps in the long term we should deprecate one, but I would
        # prefer to drop query_letters - so we need a transition period
        # with both.

        # Hack to record the claimed database size as database_length
        # (as well as in num_letters_in_database, see Bug 2176 comment 13):
        self._blast.database_length = self._blast.num_letters_in_database
        # TODO? Deprecate database_letters next?

        # Hack to record the claimed database sequence count as database_sequences
        self._blast.database_sequences = self._blast.num_sequences_in_database

        # Apply the "top level" parameter information
        self._blast.matrix = self._parameters.matrix
        self._blast.num_seqs_better_e = self._parameters.num_seqs_better_e
        self._blast.gap_penalties = self._parameters.gap_penalties
        self._blast.filter = self._parameters.filter
        self._blast.expect = self._parameters.expect
        self._blast.sc_match = self._parameters.sc_match
        self._blast.sc_mismatch = self._parameters.sc_mismatch

        # Add to the list
        self._records.append(self._blast)
        # Clear the object (a new empty one is create in _start_Iteration)
        self._blast = None

        if self._debug:
            print("NCBIXML: Added Blast record to results")

    # Header
    def _set_header_application(self):
        """BLAST program, e.g., blastp, blastn, etc. (PRIVATE).

        Save this to put on each blast record object
        """
        self._header.application = self._value.upper()

    def _set_header_version(self):
        """Version number and date of the BLAST engine (PRIVATE).

        e.g. "BLASTX 2.2.12 [Aug-07-2005]" but there can also be
        variants like "BLASTP 2.2.18+" without the date.

        Save this to put on each blast record object
        """
        parts = self._value.split()
        # TODO - Check the first word starts with BLAST?

        # The version is the second word (field one)
        self._header.version = parts[1]

        # Check there is a third word (the date)
        if len(parts) >= 3:
            if parts[2][0] == "[" and parts[2][-1] == "]":
                self._header.date = parts[2][1:-1]
            else:
                # Assume this is still a date, but without the
                # square brackets
                self._header.date = parts[2]

    def _set_header_reference(self):
        """Record any article reference describing the algorithm (PRIVATE).

        Save this to put on each blast record object
        """
        self._header.reference = self._value

    def _set_header_database(self):
        """Record the database(s) searched (PRIVATE).

        Save this to put on each blast record object
        """
        self._header.database = self._value

    def _set_header_query_id(self):
        """Record the identifier of the query (PRIVATE).

        Important in old pre 2.2.14 BLAST, for recent versions
        <Iteration_query-ID> is enough
        """
        self._header.query_id = self._value

    def _set_header_query(self):
        """Record the definition line of the query (PRIVATE).

        Important in old pre 2.2.14 BLAST, for recent versions
        <Iteration_query-def> is enough
        """
        self._header.query = self._value

    def _set_header_query_letters(self):
        """Record the length of the query (PRIVATE).

        Important in old pre 2.2.14 BLAST, for recent versions
        <Iteration_query-len> is enough
        """
        self._header.query_letters = int(self._value)

    def _set_record_query_id(self):
        """Record the identifier of the query (PRIVATE)."""
        self._blast.query_id = self._value

    def _set_record_query_def(self):
        """Record the definition line of the query (PRIVATE)."""
        self._blast.query = self._value

    def _set_record_query_letters(self):
        """Record the length of the query (PRIVATE)."""
        self._blast.query_letters = int(self._value)

    # def _end_BlastOutput_query_seq(self):
    #     """The query sequence (PRIVATE)."""
    #     pass # XXX Missing in Record.Blast ?

    # def _end_BlastOutput_iter_num(self):
    #     """The psi-blast iteration number (PRIVATE)."""
    #     pass # XXX TODO PSI

    def _set_record_hits(self):
        """Hits to the database sequences, one for every sequence (PRIVATE)."""
        self._blast.num_hits = int(self._value)

    # def _end_BlastOutput_message(self):
    #     """error messages (PRIVATE)."""
    #     pass # XXX What to do ?

    # Parameters
    def _set_parameters_matrix(self):
        """Matrix used (-M on legacy BLAST) (PRIVATE)."""
        self._parameters.matrix = self._value

    def _set_parameters_expect(self):
        """Expect values cutoff (PRIVATE)."""
        # NOTE: In old text output there was a line:
        # Number of sequences better than 1.0e-004: 1
        # As far as I can see, parameters.num_seqs_better_e
        # would take the value of 1, and the expectation
        # value was not recorded.
        #
        # Anyway we should NOT record this against num_seqs_better_e
        self._parameters.expect = self._value

    # def _end_Parameters_include(self):
    #     """Inclusion threshold for a psi-blast iteration (-h) (PRIVATE)."""
    #     pass # XXX TODO PSI

    def _set_parameters_sc_match(self):
        """Match score for nucleotide-nucleotide comparison (-r) (PRIVATE)."""
        self._parameters.sc_match = int(self._value)

    def _set_parameters_sc_mismatch(self):
        """Mismatch penalty for nucleotide-nucleotide comparison (-r) (PRIVATE)."""
        self._parameters.sc_mismatch = int(self._value)

    def _set_parameters_gap_penalties(self):
        """Gap existence cost (-G) (PRIVATE)."""
        self._parameters.gap_penalties = int(self._value)

    def _set_parameters_gap_extend(self):
        """Gap extension cose (-E) (PRIVATE)."""
        self._parameters.gap_penalties = (
            self._parameters.gap_penalties,
            int(self._value),
        )

    def _set_parameters_filter(self):
        """Record filtering options (-F) (PRIVATE)."""
        self._parameters.filter = self._value

    # def _end_Parameters_pattern(self):
    #     """Pattern used for phi-blast search (PRIVATE).
    #     """
    #     pass # XXX TODO PSI

    # def _end_Parameters_entrez_query(self):
    #     """Entrez query used to limit search (PRIVATE).
    #     """
    #     pass # XXX TODO PSI

    # Hits
    def _start_hit(self):
        """Start filling records (PRIVATE)."""
        self._blast.alignments.append(Record.Alignment())
        self._descr = (
            Record.Description() if self._xml_version == 1 else Record.DescriptionExt()
        )
        self._blast.descriptions.append(self._descr)
        self._blast.multiple_alignment = []
        self._hit = self._blast.alignments[-1]

        self._descr.num_alignments = 0

    def _end_hit(self):
        """Clear variables (PRIVATE)."""
        # Cleanup
        self._blast.multiple_alignment = None
        self._hit = None
        self._descr = None

    def set_hit_id(self):
        """Record the identifier of the database sequence (PRIVATE)."""
        self._hit.hit_id = self._value
        self._hit.title = self._value + " "

    def set_hit_def(self):
        """Record the definition line of the database sequence (PRIVATE)."""
        self._hit.hit_def = self._value
        self._hit.title += self._value
        self._descr.title = self._hit.title

    def set_hit_accession(self):
        """Record the accession value of the database sequence (PRIVATE)."""
        self._hit.accession = self._value
        self._descr.accession = self._value

    def set_hit_len(self):
        """Record the length of the hit."""
        self._hit.length = int(self._value)

    # HSPs
    def _start_hsp(self):
        # Note that self._start_Hit() should have been called
        # to setup things like self._blast.multiple_alignment
        self._hsp = Record.HSP()
        self._hsp.positives = None
        self._hit.hsps.append(self._hsp)
        self._descr.num_alignments += 1
        self._blast.multiple_alignment.append(Record.MultipleAlignment())
        self._mult_al = self._blast.multiple_alignment[-1]

    def _end_hsp(self):
        if self._hsp.frame and len(self._hsp.frame) == 1:
            self._hsp.frame += (0,)

    # Hsp_num is useless
    def _set_hsp_score(self):
        """Record the raw score of HSP (PRIVATE)."""
        self._hsp.score = float(self._value)
        if self._descr.score is None:
            self._descr.score = float(self._value)

    def _set_hsp_bit_score(self):
        """Record the Bit score of HSP (PRIVATE)."""
        self._hsp.bits = float(self._value)
        if self._descr.bits is None:
            self._descr.bits = float(self._value)

    def _set_hsp_e_value(self):
        """Record the expect value of the HSP (PRIVATE)."""
        self._hsp.expect = float(self._value)
        if self._descr.e is None:
            self._descr.e = float(self._value)

    def _set_hsp_query_start(self):
        """Offset of query at the start of the alignment (one-offset) (PRIVATE)."""
        self._hsp.query_start = int(self._value)

    def _set_hsp_query_end(self):
        """Offset of query at the end of the alignment (one-offset) (PRIVATE)."""
        self._hsp.query_end = int(self._value)

    def _set_hsp_hit_from(self):
        """Offset of the database at the start of the alignment (one-offset) (PRIVATE)."""
        self._hsp.sbjct_start = int(self._value)

    def _set_hsp_hit_to(self):
        """Offset of the database at the end of the alignment (one-offset) (PRIVATE)."""
        self._hsp.sbjct_end = int(self._value)

    # def _end_Hsp_pattern_from(self):
    #     """Start of phi-blast pattern on the query (one-offset) (PRIVATE)."""
    #     pass # XXX TODO PSI

    # def _end_Hsp_pattern_to(self):
    #     """End of phi-blast pattern on the query (one-offset) (PRIVATE)."""
    #     pass # XXX TODO PSI

    def _set_hsp_query_frame(self):
        """Frame of the query if applicable (PRIVATE)."""
        v = int(self._value)
        self._hsp.frame = (v,)
        if self._header.application == "BLASTN":
            self._hsp.strand = ("Plus" if v > 0 else "Minus",)

    def _set_hsp_hit_frame(self):
        """Frame of the database sequence if applicable (PRIVATE)."""
        v = int(self._value)
        if len(self._hsp.frame) == 0:
            self._hsp.frame = (0, v)
        else:
            self._hsp.frame += (v,)
        if self._header.application == "BLASTN":
            self._hsp.strand += ("Plus" if v > 0 else "Minus",)

    def _set_hsp_query_strand(self):
        """Frame of the query if applicable (PRIVATE)."""
        self._hsp.strand = (self._value,)
        if self._header.application == "BLASTN":
            self._hsp.frame = (1 if self._value == "Plus" else -1,)

    def _set_hsp_hit_strand(self):
        """Frame of the database sequence if applicable (PRIVATE)."""
        self._hsp.strand += (self._value,)
        if self._header.application == "BLASTN":
            self._hsp.frame += (1 if self._value == "Plus" else -1,)

    def _set_hsp_identity(self):
        """Record the number of identities in the alignment (PRIVATE)."""
        v = int(self._value)
        self._hsp.identities = v
        if self._hsp.positives is None:
            self._hsp.positives = v

    def _set_hsp_positive(self):
        """Record the number of positive (conservative) substitutions in the alignment (PRIVATE)."""
        self._hsp.positives = int(self._value)

    def _set_hsp_gaps(self):
        """Record the number of gaps in the alignment (PRIVATE)."""
        self._hsp.gaps = int(self._value)

    def _set_hsp_align_len(self):
        """Record the length of the alignment (PRIVATE)."""
        self._hsp.align_length = int(self._value)

    # def _en_Hsp_density(self):
    #     """Score density (PRIVATE)."""
    #     pass # XXX ???

    def _set_hsp_query_seq(self):
        """Record the alignment string for the query (PRIVATE)."""
        self._hsp.query = self._value

    def _set_hsp_subject_seq(self):
        """Record the alignment string for the database (PRIVATE)."""
        self._hsp.sbjct = self._value

    def _set_hsp_midline(self):
        """Record the middle line as normally seen in BLAST report (PRIVATE)."""
        self._hsp.match = self._value  # do NOT strip spaces!
        assert len(self._hsp.match) == len(self._hsp.query)
        assert len(self._hsp.match) == len(self._hsp.sbjct)

    # Statistics
    def _set_statistics_db_num(self):
        """Record the number of sequences in the database (PRIVATE)."""
        self._blast.num_sequences_in_database = int(self._value)

    def _set_statistics_db_len(self):
        """Record the number of letters in the database (PRIVATE)."""
        self._blast.num_letters_in_database = int(self._value)

    def _set_statistics_hsp_len(self):
        """Record the effective HSP length (PRIVATE)."""
        self._blast.effective_hsp_length = int(self._value)

    def _set_statistics_eff_space(self):
        """Record the effective search space (PRIVATE)."""
        self._blast.effective_search_space = float(self._value)

    def _set_statistics_kappa(self):
        """Karlin-Altschul parameter K (PRIVATE)."""
        self._blast.ka_params = float(self._value)

    def _set_statistics_lambda(self):
        """Karlin-Altschul parameter Lambda (PRIVATE)."""
        self._blast.ka_params = (float(self._value), self._blast.ka_params)

    def _set_statistics_entropy(self):
        """Karlin-Altschul parameter H (PRIVATE)."""
        self._blast.ka_params = self._blast.ka_params + (float(self._value),)

    def _start_hit_descr_item(self):
        """XML v2. Start hit description item."""
        self._hit_descr_item = Record.DescriptionExtItem()

    def _end_hit_descr_item(self):
        """XML v2. Start hit description item."""
        self._descr.append_item(self._hit_descr_item)
        if not self._hit.title:
            self._hit.title = str(self._hit_descr_item)
        self._hit_descr_item = None

    def _end_description_id(self):
        """XML v2. The identifier of the database sequence(PRIVATE)."""
        self._hit_descr_item.id = self._value
        if not self._hit.hit_id:
            self._hit.hit_id = self._value

    def _end_description_accession(self):
        """XML v2. The accession value of the database sequence (PRIVATE)."""
        self._hit_descr_item.accession = self._value
        if not getattr(self._hit, "accession", None):
            self._hit.accession = self._value

    def _end_description_title(self):
        """XML v2. The hit description title (PRIVATE)."""
        self._hit_descr_item.title = self._value

    def _end_description_taxid(self):
        try:
            self._hit_descr_item.taxid = int(self._value)
        except ValueError:
            pass

    def _end_description_sciname(self):
        self._hit_descr_item.sciname = self._value


def read(handle, debug=0):
    """Return a single Blast record (assumes just one query).

    Uses the BlastParser internally.

    This function is for use when there is one and only one BLAST
    result in your XML file.

    Use the Bio.Blast.NCBIXML.parse() function if you expect more than
    one BLAST record (i.e. if you have more than one query sequence).
    """
    iterator = parse(handle, debug)
    try:
        record = next(iterator)
    except StopIteration:
        raise ValueError("No records found in handle") from None
    try:
        next(iterator)
        raise ValueError("More than one record found in handle")
    except StopIteration:
        pass
    return record


def parse(handle, debug=0):
    """Return an iterator a Blast record for each query.

    Incremental parser, this is an iterator that returns
    Blast records.  It uses the BlastParser internally.

    handle - file handle to and XML file to parse
    debug - integer, amount of debug information to print

    This is a generator function that returns multiple Blast records
    objects - one for each query sequence given to blast.  The file
    is read incrementally, returning complete records as they are read
    in.

    Should cope with new BLAST 2.2.14+ which gives a single XML file
    for multiple query records.

    Should also cope with XML output from older versions BLAST which
    gave multiple XML files concatenated together (giving a single file
    which strictly speaking wasn't valid XML).
    """
    from xml.parsers import expat

    BLOCK = 1024
    MARGIN = 10  # must be at least length of newline + XML start
    XML_START = "<?xml"
    NEW_LINE = "\n"
    NULL = ""

    pending = ""
    text = handle.read(BLOCK)
    if isinstance(text, bytes):
        # Not a text handle, raw bytes mode
        XML_START = b"<?xml"
        NEW_LINE = b"\n"
        NULL = b""
        pending = b""

    if not text:
        # NO DATA FOUND!
        raise ValueError("Your XML file was empty")

    while text:
        # We are now starting a new XML file
        if not text.startswith(XML_START):
            raise ValueError(
                "Your XML file did not start with %r... but instead %r"
                % (XML_START, text[:20])
            )

        expat_parser = expat.ParserCreate()
        blast_parser = BlastParser(debug)
        expat_parser.StartElementHandler = blast_parser.startElement
        expat_parser.EndElementHandler = blast_parser.endElement
        expat_parser.CharacterDataHandler = blast_parser.characters

        expat_parser.Parse(text, False)
        while blast_parser._records:
            record = blast_parser._records[0]
            blast_parser._records = blast_parser._records[1:]
            yield record

        while True:
            # Read in another block of the file...
            text, pending = pending + handle.read(BLOCK), ""
            if not text:
                # End of the file!
                expat_parser.Parse(NULL, True)  # End of XML record
                break

            # Now read a little bit more so we can check for the
            # start of another XML file...
            pending = handle.read(MARGIN)

            if (NEW_LINE + XML_START) not in (text + pending):
                # Good - still dealing with the same XML file
                expat_parser.Parse(text, False)
                while blast_parser._records:
                    yield blast_parser._records.pop(0)
            else:
                # This is output from pre 2.2.14 BLAST,
                # one XML file for each query!

                # Finish the old file:
                text, pending = (text + pending).split(NEW_LINE + XML_START, 1)
                pending = XML_START + pending

                expat_parser.Parse(text, True)  # End of XML record
                while blast_parser._records:
                    yield blast_parser._records.pop(0)

                # Now we are going to re-loop, reset the
                # parsers and start reading the next XML file
                text, pending = pending, NULL
                break

        # this was added because it seems that the Jython expat parser
        # was adding records later then the Python one
        while blast_parser._records:
            yield blast_parser._records.pop(0)

        # At this point we have finished the first XML record.
        # If the file is from an old version of blast, it may
        # contain more XML records (check if text=="").
        assert not pending, pending
        assert len(blast_parser._records) == 0, len(blast_parser._records)

    # We should have finished the file!
    assert not text, text
    assert not pending, pending
    assert len(blast_parser._records) == 0, len(blast_parser._records)
