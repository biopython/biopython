# Copyright 2001 Brad Chapman.
# Revisions copyright 2009-2010 by Peter Cock.
# Revisions copyright 2010 by Phillip Garland.
# All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Definitions for interacting with BLAST related applications (OBSOLETE).

Wrappers for the new NCBI BLAST+ tools (written in C++):

 - NcbiblastpCommandline - Protein-Protein BLAST
 - NcbiblastnCommandline - Nucleotide-Nucleotide BLAST
 - NcbiblastxCommandline - Translated Query-Protein Subject BLAST
 - NcbitblastnCommandline - Protein Query-Translated Subject BLAST
 - NcbitblastxCommandline - Translated Query-Protein Subject BLAST
 - NcbipsiblastCommandline - Position-Specific Initiated BLAST
 - NcbirpsblastCommandline - Reverse Position Specific BLAST
 - NcbirpstblastnCommandline - Translated Reverse Position Specific BLAST
 - NcbideltablastCommandline - Protein-Protein domain enhanced lookup time accelerated blast
 - NcbiblastformatterCommandline - Convert ASN.1 to other BLAST output formats
 - NcbimakeblastdbCommandline - Application to create BLAST databases

For further details, see:

Camacho et al. BLAST+: architecture and applications
BMC Bioinformatics 2009, 10:421
https://doi.org/10.1186/1471-2105-10-421

We have decided to remove this module in future, and instead recommend
building your command and invoking it via the subprocess module directly.
"""

from Bio.Application import _Option, AbstractCommandline, _Switch


class _NcbibaseblastCommandline(AbstractCommandline):
    """Base Commandline object for (new) NCBI BLAST+ wrappers (PRIVATE).

    This is provided for subclassing, it deals with shared options
    common to all the BLAST tools (blastn, rpsblast, rpsblast, etc
    AND blast_formatter).
    """

    def __init__(self, cmd=None, **kwargs):
        assert cmd is not None
        extra_parameters = [
            # Core:
            _Switch(
                ["-h", "h"], "Print USAGE and DESCRIPTION;  ignore other arguments."
            ),
            _Switch(
                ["-help", "help"],
                "Print USAGE, DESCRIPTION and ARGUMENTS description; "
                "ignore other arguments.",
            ),
            _Switch(
                ["-version", "version"],
                "Print version number;  ignore other arguments.",
            ),
            # Output configuration options
            _Option(
                ["-out", "out"],
                "Output file for alignment.",
                filename=True,
                equate=False,
            ),
            # Formatting options:
            _Option(
                ["-outfmt", "outfmt"],
                "Alignment view.  Typically an integer 0-14 but for some "
                "formats can be named columns like '6 qseqid sseqid'.  "
                "Use 5 for XML output (differs from classic BLAST which "
                "used 7 for XML).",
                filename=True,  # to ensure spaced inputs are quoted
                equate=False,
            ),
            # TODO - Document and test the column options
            _Switch(["-show_gis", "show_gis"], "Show NCBI GIs in deflines?"),
            _Option(
                ["-num_descriptions", "num_descriptions"],
                "Number of database sequences to show one-line descriptions for.\n\n"
                "Integer argument (at least zero). Default is 500. "
                "See also num_alignments.",
                equate=False,
            ),
            _Option(
                ["-num_alignments", "num_alignments"],
                "Number of database sequences to show num_alignments for.\n\n"
                "Integer argument (at least zero). Default is 200. "
                "See also num_alignments.",
                equate=False,
            ),
            _Option(
                ["-line_length", "line_length"],
                "Line length for formatting alignments "
                "(integer, at least 1, default 60).\n\n"
                "Not applicable for outfmt > 4. Added in BLAST+ 2.2.30.",
                equate=False,
            ),
            _Switch(
                ["-html", "html"], "Produce HTML output? See also the outfmt option."
            ),
            # Miscellaneous options
            _Switch(
                ["-parse_deflines", "parse_deflines"],
                "Should the query and subject defline(s) be parsed?",
            ),
        ]
        try:
            # Insert extra parameters - at the start just in case there
            # are any arguments which must come last:
            self.parameters = extra_parameters + self.parameters
        except AttributeError:
            # Should we raise an error?  The subclass should have set this up!
            self.parameters = extra_parameters
        AbstractCommandline.__init__(self, cmd, **kwargs)

    def _validate_incompatibilities(self, incompatibles):
        """Validate parameters for incompatibilities (PRIVATE).

        Used by the _validate method.
        """
        for a in incompatibles:
            if self._get_parameter(a):
                for b in incompatibles[a]:
                    if self._get_parameter(b):
                        raise ValueError("Options %s and %s are incompatible." % (a, b))


class _NcbiblastCommandline(_NcbibaseblastCommandline):
    """Base Commandline object for (new) NCBI BLAST+ wrappers (PRIVATE).

    This is provided for subclassing, it deals with shared options
    common to all the BLAST tools (blastn, rpsblast, rpsblast, etc).
    """

    def __init__(self, cmd=None, **kwargs):
        assert cmd is not None
        extra_parameters = [
            # Input query options:
            _Option(
                ["-query", "query"],
                "The sequence to search with.",
                filename=True,
                equate=False,
            ),  # Should this be required?
            _Option(
                ["-query_loc", "query_loc"],
                "Location on the query sequence (Format: start-stop).",
                equate=False,
            ),
            # General search options:
            _Option(["-db", "db"], "The database to BLAST against.", equate=False),
            _Option(["-evalue", "evalue"], "Expectation value cutoff.", equate=False),
            _Option(
                ["-word_size", "word_size"],
                "Word size for wordfinder algorithm.\n\nInteger. Minimum 2.",
                equate=False,
            ),
            # BLAST-2-Sequences options:
            # - see subclass
            # Formatting options:
            # - see baseclass
            # Query filtering options
            _Option(
                ["-soft_masking", "soft_masking"],
                "Apply filtering locations as soft masks (Boolean, Default = true).",
                equate=False,
            ),
            _Switch(
                ["-lcase_masking", "lcase_masking"],
                "Use lower case filtering in query and subject sequence(s)?",
            ),
            # Restrict search or results
            _Option(
                ["-gilist", "gilist"],
                "Restrict search of database to list of GI's.\n\n"
                "Incompatible with: negative_gilist, seqidlist, negative_seqidlist, "
                "remote, subject, subject_loc",
                filename=True,
                equate=False,
            ),
            _Option(
                ["-negative_gilist", "negative_gilist"],
                "Restrict search of database to everything except the listed GIs.\n\n"
                "Incompatible with: gilist, seqidlist, remote, subject, subject_loc",
                filename=True,
                equate=False,
            ),
            _Option(
                ["-seqidlist", "seqidlist"],
                "Restrict search of database to list of SeqID's.\n\n"
                "Incompatible with: gilist, negative_gilist, remote, subject, "
                "subject_loc",
                filename=True,
                equate=False,
            ),
            _Option(
                ["-negative_seqidlist", "negative_seqidlist"],
                "Restrict search of database to everything except listed SeqID's.\n\n"
                "Incompatible with: gilist, seqidlist, remote, subject, subject_loc",
                filename=True,
                equate=False,
            ),
            _Option(
                ["-entrez_query", "entrez_query"],
                "Restrict search with the given Entrez query (requires remote).",
                equate=False,
            ),
            _Option(
                ["-qcov_hsp_perc", "qcov_hsp_perc"],
                "Percent query coverage per hsp (float, 0 to 100).\n\n"
                "Added in BLAST+ 2.2.30.",
                equate=False,
            ),
            _Option(
                ["-max_target_seqs", "max_target_seqs"],
                "Maximum number of aligned sequences to keep (integer, at least one).",
                equate=False,
            ),
            # Statistical options
            _Option(
                ["-dbsize", "dbsize"],
                "Effective length of the database (integer).",
                equate=False,
            ),
            _Option(
                ["-searchsp", "searchsp"],
                "Effective length of the search space (integer).",
                equate=False,
            ),
            _Option(
                ["-max_hsps_per_subject", "max_hsps_per_subject"],
                "Override max number of HSPs per subject saved for ungapped searches "
                "(integer).",
                equate=False,
            ),
            _Option(
                ["-max_hsps", "max_hsps"],
                "Set max number of HSPs saved per subject sequence\n\n"
                "Ddefault 0 means no limit.",
                equate=False,
            ),
            _Switch(["-sum_statistics", "sum_statistics"], "Use sum statistics."),
            # Is -sum_stats a BLAST+ bug, why not use -sum_statistics switch?
            _Option(
                ["-sum_stats", "sum_stats"],
                "Use sum statistics (boolean).\n\nAdded in BLAST+ 2.2.30.",
                equate=False,
            ),
            # Extension options
            _Option(
                ["-xdrop_ungap", "xdrop_ungap"],
                "X-dropoff value (in bits) for ungapped extensions (float).",
                equate=False,
            ),
            _Option(
                ["-xdrop_gap", "xdrop_gap"],
                "X-dropoff value (in bits) for preliminary gapped extensions (float).",
                equate=False,
            ),
            _Option(
                ["-xdrop_gap_final", "xdrop_gap_final"],
                "X-dropoff value (in bits) for final gapped alignment (float).",
                equate=False,
            ),
            _Option(
                ["-window_size", "window_size"],
                "Multiple hits window size, use 0 to specify 1-hit algorithm "
                "(integer).",
                equate=False,
            ),
            # Search strategy options
            _Option(
                ["-import_search_strategy", "import_search_strategy"],
                "Search strategy to use.\n\n"
                "Incompatible with: export_search_strategy",
                filename=True,
                equate=False,
            ),
            _Option(
                ["-export_search_strategy", "export_search_strategy"],
                "File name to record the search strategy used.\n\n"
                "Incompatible with: import_search_strategy",
                filename=True,
                equate=False,
            ),
            # Miscellaneous options
            _Option(
                ["-num_threads", "num_threads"],
                "Number of threads to use in the BLAST search.\n\n"
                "Integer, at least one. Default is one. Incompatible with: remote",
                equate=False,
            ),
            _Switch(
                ["-remote", "remote"],
                "Execute search remotely?\n\n"
                "Incompatible with: gilist, negative_gilist, subject_loc, "
                "num_threads, ...",
            ),
        ]
        try:
            # Insert extra parameters - at the start just in case there
            # are any arguments which must come last:
            self.parameters = extra_parameters + self.parameters
        except AttributeError:
            # Should we raise an error?  The subclass should have set this up!
            self.parameters = extra_parameters
        _NcbibaseblastCommandline.__init__(self, cmd, **kwargs)

    def _validate(self):
        incompatibles = {
            "remote": ["gilist", "negative_gilist", "num_threads"],
            "import_search_strategy": ["export_search_strategy"],
            "gilist": ["negative_gilist"],
            "seqidlist": ["gilist", "negative_gilist", "remote"],
        }
        self._validate_incompatibilities(incompatibles)
        if self.entrez_query and not self.remote:
            raise ValueError("Option entrez_query requires remote option.")
        AbstractCommandline._validate(self)


class _Ncbiblast2SeqCommandline(_NcbiblastCommandline):
    """Base Commandline object for (new) NCBI BLAST+ wrappers (PRIVATE).

    This is provided for subclassing, it deals with shared options
    common to all the BLAST tools supporting two-sequence BLAST
    (blastn, psiblast, etc) but not rpsblast or rpstblastn.
    """

    def __init__(self, cmd=None, **kwargs):
        assert cmd is not None
        extra_parameters = [
            # General search options:
            _Option(
                ["-gapopen", "gapopen"], "Cost to open a gap (integer).", equate=False
            ),
            _Option(
                ["-gapextend", "gapextend"],
                "Cost to extend a gap (integer).",
                equate=False,
            ),
            # BLAST-2-Sequences options:
            _Option(
                ["-subject", "subject"],
                "Subject sequence(s) to search.\n\n"
                "Incompatible with: db, gilist, seqidlist, negative_gilist, "
                "negative_seqidlist, db_soft_mask, db_hard_mask\n\n"
                "See also subject_loc.",
                filename=True,
                equate=False,
            ),
            _Option(
                ["-subject_loc", "subject_loc"],
                "Location on the subject sequence (Format: start-stop).\n\n"
                "Incompatible with: db, gilist, seqidlist, negative_gilist, "
                "negative_seqidlist, db_soft_mask, db_hard_mask, remote.\n\n"
                "See also subject.",
                equate=False,
            ),
            # Restrict search or results:
            _Option(
                ["-culling_limit", "culling_limit"],
                "Hit culling limit (integer).\n\n"
                "If the query range of a hit is enveloped by that of at "
                "least this many higher-scoring hits, delete the hit.\n\n"
                "Incompatible with: best_hit_overhang, best_hit_score_edge.",
                equate=False,
            ),
            _Option(
                ["-best_hit_overhang", "best_hit_overhang"],
                "Best Hit algorithm overhang value (float, recommended value: 0.1)\n\n"
                "Float between 0.0 and 0.5 inclusive. "
                "Incompatible with: culling_limit.",
                equate=False,
            ),
            _Option(
                ["-best_hit_score_edge", "best_hit_score_edge"],
                "Best Hit algorithm score edge value (float).\n\n"
                "Float between 0.0 and 0.5 inclusive. Recommended value: 0.1\n\n"
                "Incompatible with: culling_limit.",
                equate=False,
            ),
        ]
        try:
            # Insert extra parameters - at the start just in case there
            # are any arguments which must come last:
            self.parameters = extra_parameters + self.parameters
        except AttributeError:
            # Should we raise an error?  The subclass should have set this up!
            self.parameters = extra_parameters
        _NcbiblastCommandline.__init__(self, cmd, **kwargs)

    def _validate(self):
        incompatibles = {
            "subject_loc": ["db", "gilist", "negative_gilist", "seqidlist", "remote"],
            "culling_limit": ["best_hit_overhang", "best_hit_score_edge"],
            "subject": ["db", "gilist", "negative_gilist", "seqidlist"],
        }
        self._validate_incompatibilities(incompatibles)
        _NcbiblastCommandline._validate(self)


class _NcbiblastMain2SeqCommandline(_Ncbiblast2SeqCommandline):
    """Base Commandline object for (new) NCBI BLAST+ wrappers (PRIVATE).

    This is provided for subclassing, it deals with shared options
    common to the main BLAST tools blastp, blastn, blastx, tblastx, tblastn
    but not psiblast, rpsblast or rpstblastn.
    """

    def __init__(self, cmd=None, **kwargs):
        assert cmd is not None
        extra_parameters = [
            # Restrict search or results:
            _Option(
                ["-db_soft_mask", "db_soft_mask"],
                "Filtering algorithm for soft masking (integer).\n\n"
                "Filtering algorithm ID to apply to BLAST database as soft masking. "
                "Incompatible with: db_hard_mask, subject, subject_loc",
                equate=False,
            ),
            _Option(
                ["-db_hard_mask", "db_hard_mask"],
                "Filtering algorithm for hard masking (integer).\n\n"
                "Filtering algorithm ID to apply to BLAST database as hard masking. "
                "Incompatible with: db_soft_mask, subject, subject_loc",
                equate=False,
            ),
        ]
        try:
            # Insert extra parameters - at the start just in case there
            # are any arguments which must come last:
            self.parameters = extra_parameters + self.parameters
        except AttributeError:
            # Should we raise an error?  The subclass should have set this up!
            self.parameters = extra_parameters
        _Ncbiblast2SeqCommandline.__init__(self, cmd, **kwargs)

    def _validate(self):
        incompatibles = {
            "db_soft_mask": ["db_hard_mask", "subject", "subject_loc"],
            "db_hard_mask": ["db_soft_mask", "subject", "subject_loc"],
        }
        self._validate_incompatibilities(incompatibles)
        _Ncbiblast2SeqCommandline._validate(self)


class NcbiblastpCommandline(_NcbiblastMain2SeqCommandline):
    """Create a commandline for the NCBI BLAST+ program blastp (for proteins).

    With the release of BLAST+ (BLAST rewritten in C++ instead of C), the NCBI
    replaced the old blastall tool with separate tools for each of the searches.
    This wrapper therefore replaces BlastallCommandline with option -p blastp.

    >>> from Bio.Blast.Applications import NcbiblastpCommandline
    >>> cline = NcbiblastpCommandline(query="rosemary.pro", db="nr",
    ...                               evalue=0.001, remote=True, ungapped=True)
    >>> cline
    NcbiblastpCommandline(cmd='blastp', query='rosemary.pro', db='nr', evalue=0.001, remote=True, ungapped=True)
    >>> print(cline)
    blastp -query rosemary.pro -db nr -evalue 0.001 -remote -ungapped

    You would typically run the command line with cline() or via the Python
    subprocess module, as described in the Biopython tutorial.
    """

    def __init__(self, cmd="blastp", **kwargs):
        """Initialize the class."""
        self.parameters = [
            # General search options:
            _Option(
                ["-task", "task"],
                "Task to execute (string, blastp (default), blastp-fast or blastp-short).",
                checker_function=lambda value: value
                in ["blastp", "blastp-fast", "blastp-short"],
                equate=False,
            ),
            _Option(["-matrix", "matrix"], "Scoring matrix name (default BLOSUM62)."),
            _Option(
                ["-threshold", "threshold"],
                "Minimum score for words to be added to the BLAST lookup table (float).",
                equate=False,
            ),
            _Option(
                ["-comp_based_stats", "comp_based_stats"],
                "Use composition-based statistics (string, default 2, i.e. True).\n\n"
                "0, F or f: no composition-based statistics\n\n"
                "2, T or t, D or d : Composition-based score adjustment as in "
                "Bioinformatics 21:902-911, 2005, conditioned on sequence "
                "properties\n\n"
                "Note that tblastn also supports values of 1 and 3.",
                checker_function=lambda value: value in "0Ft2TtDd",
                equate=False,
            ),
            # Query filtering options:
            _Option(
                ["-seg", "seg"],
                "Filter query sequence with SEG (string).\n\n"
                'Format: "yes", "window locut hicut", or "no" to disable\n'
                'Default is "12 2.2 2.5"',
                equate=False,
            ),
            # Extension options:
            _Switch(["-ungapped", "ungapped"], "Perform ungapped alignment only?"),
            # Miscellaneous options:
            _Switch(
                ["-use_sw_tback", "use_sw_tback"],
                "Compute locally optimal Smith-Waterman alignments?",
            ),
        ]
        _NcbiblastMain2SeqCommandline.__init__(self, cmd, **kwargs)


class NcbiblastnCommandline(_NcbiblastMain2SeqCommandline):
    """Wrapper for the NCBI BLAST+ program blastn (for nucleotides).

    With the release of BLAST+ (BLAST rewritten in C++ instead of C), the NCBI
    replaced the old blastall tool with separate tools for each of the searches.
    This wrapper therefore replaces BlastallCommandline with option -p blastn.

    For example, to run a search against the "nt" nucleotide database using the
    FASTA nucleotide file "m_code.fasta" as the query, with an expectation value
    cut off of 0.001, saving the output to a file in XML format:

    >>> from Bio.Blast.Applications import NcbiblastnCommandline
    >>> cline = NcbiblastnCommandline(query="m_cold.fasta", db="nt", strand="plus",
    ...                               evalue=0.001, out="m_cold.xml", outfmt=5)
    >>> cline
    NcbiblastnCommandline(cmd='blastn', out='m_cold.xml', outfmt=5, query='m_cold.fasta', db='nt', evalue=0.001, strand='plus')
    >>> print(cline)
    blastn -out m_cold.xml -outfmt 5 -query m_cold.fasta -db nt -evalue 0.001 -strand plus

    You would typically run the command line with cline() or via the Python
    subprocess module, as described in the Biopython tutorial.
    """

    def __init__(self, cmd="blastn", **kwargs):
        """Initialize the class."""
        self.parameters = [
            # Input query options:
            _Option(
                ["-strand", "strand"],
                "Query strand(s) to search against database/subject.\n\n"
                'Values allowed are "both" (default), "minus", "plus".',
                checker_function=lambda value: value in ["both", "minus", "plus"],
                equate=False,
            ),
            # General search options:
            _Option(
                ["-task", "task"],
                "Task to execute (string, default 'megablast')\n\n"
                "Allowed values 'blastn', 'blastn-short', 'dc-megablast', 'megablast' "
                "(the default), or 'vecscreen'.",
                checker_function=lambda value: value
                in ["blastn", "blastn-short", "dc-megablast", "megablast", "vecscreen"],
                equate=False,
            ),
            _Option(
                ["-penalty", "penalty"],
                "Penalty for a nucleotide mismatch (integer, at most zero).",
                equate=False,
            ),
            _Option(
                ["-reward", "reward"],
                "Reward for a nucleotide match (integer, at least zero).",
                equate=False,
            ),
            _Option(
                ["-use_index", "use_index"],
                "Use MegaBLAST database index (Boolean, Default = False)",
                equate=False,
            ),
            _Option(
                ["-index_name", "index_name"],
                "MegaBLAST database index name.",
                equate=False,
            ),
            # Query filtering options:
            _Option(
                ["-dust", "dust"],
                "Filter query sequence with DUST (string).\n\n"
                "Format: 'yes', 'level window linker', or 'no' to disable.\n\n"
                "Default = '20 64 1'.",
                equate=False,
            ),
            _Option(
                ["-filtering_db", "filtering_db"],
                "BLAST database containing filtering elements (i.e. repeats).",
                equate=False,
            ),
            _Option(
                ["-window_masker_taxid", "window_masker_taxid"],
                "Enable WindowMasker filtering using a Taxonomic ID (integer).",
                equate=False,
            ),
            _Option(
                ["-window_masker_db", "window_masker_db"],
                "Enable WindowMasker filtering using this repeats database (string).",
                equate=False,
            ),
            # Restrict search or results:
            _Option(
                ["-perc_identity", "perc_identity"],
                "Percent identity (real, 0 to 100 inclusive).",
                equate=False,
            ),
            # Discontiguous MegaBLAST options
            _Option(
                ["-template_type", "template_type"],
                "Discontiguous MegaBLAST template type (string).\n\n"
                "Allowed values: 'coding', 'coding_and_optimal' or 'optimal'.\n"
                "Requires: template_length.",
                checker_function=lambda value: value
                in ["coding", "coding_and_optimal", "optimal"],
                equate=False,
            ),
            _Option(
                ["-template_length", "template_length"],
                "Discontiguous MegaBLAST template length (integer).\n\n"
                "Allowed values: 16, 18, 21.\n\n"
                "Requires: template_type.",
                checker_function=lambda value: value in [16, 18, 21, "16", "18", "21"],
                equate=False,
            ),
            # Extension options:
            _Switch(
                ["-no_greedy", "no_greedy"],
                "Use non-greedy dynamic programming extension",
            ),
            _Option(
                ["-min_raw_gapped_score", "min_raw_gapped_score"],
                "Minimum raw gapped score to keep an alignment in the "
                "preliminary gapped and traceback stages (integer).",
                equate=False,
            ),
            _Switch(["-ungapped", "ungapped"], "Perform ungapped alignment only?"),
            _Option(
                ["-off_diagonal_range", "off_diagonal_range"],
                "Number of off-diagonals to search for the 2nd hit (integer).\n\n"
                "Expects a positive integer, or 0 (default) to turn off."
                "Added in BLAST 2.2.23+",
                equate=False,
            ),
        ]
        _NcbiblastMain2SeqCommandline.__init__(self, cmd, **kwargs)

    def _validate(self):
        if (self.template_type and not self.template_length) or (
            self.template_length and not self.template_type
        ):
            raise ValueError(
                "Options template_type and template_type require each other."
            )
        _NcbiblastMain2SeqCommandline._validate(self)


class NcbiblastxCommandline(_NcbiblastMain2SeqCommandline):
    """Wrapper for the NCBI BLAST+ program blastx (nucleotide query, protein database).

    With the release of BLAST+ (BLAST rewritten in C++ instead of C), the NCBI
    replaced the old blastall tool with separate tools for each of the searches.
    This wrapper therefore replaces BlastallCommandline with option -p blastx.

    >>> from Bio.Blast.Applications import NcbiblastxCommandline
    >>> cline = NcbiblastxCommandline(query="m_cold.fasta", db="nr", evalue=0.001)
    >>> cline
    NcbiblastxCommandline(cmd='blastx', query='m_cold.fasta', db='nr', evalue=0.001)
    >>> print(cline)
    blastx -query m_cold.fasta -db nr -evalue 0.001

    You would typically run the command line with cline() or via the Python
    subprocess module, as described in the Biopython tutorial.
    """

    def __init__(self, cmd="blastx", **kwargs):
        """Initialize the class."""
        self.parameters = [
            # Input query options:
            _Option(
                ["-task", "task"],
                "Task to execute (string, blastx (default) or blastx-fast).",
                checker_function=lambda value: value in ["blastx", "blastx-fast"],
                equate=False,
            ),
            _Option(
                ["-strand", "strand"],
                "Query strand(s) to search against database/subject.\n\n"
                'Values allowed are "both" (default), "minus", "plus".',
                checker_function=lambda value: value in ["both", "minus", "plus"],
                equate=False,
            ),
            # Input query options:
            _Option(
                ["-query_gencode", "query_gencode"],
                "Genetic code to use to translate query (integer, default 1).",
                equate=False,
            ),
            # General search options:
            _Option(
                ["-frame_shift_penalty", "frame_shift_penalty"],
                "Frame shift penalty (integer, at least 1, default ignored) (OBSOLETE).\n\n"
                "This was removed in BLAST 2.2.27+",
                equate=False,
            ),
            _Option(
                ["-max_intron_length", "max_intron_length"],
                "Maximum intron length (integer).\n\n"
                "Length of the largest intron allowed in a translated nucleotide "
                "sequence when linking multiple distinct alignments (a negative "
                "value disables linking). Default zero.",
                equate=False,
            ),
            _Option(
                ["-matrix", "matrix"],
                "Scoring matrix name (default BLOSUM62).",
                equate=False,
            ),
            _Option(
                ["-threshold", "threshold"],
                "Minimum score for words to be added to the BLAST lookup table (float).",
                equate=False,
            ),
            _Option(
                ["-comp_based_stats", "comp_based_stats"],
                "Use composition-based statistics for blastp, blastx, or tblastn.\n\n"
                "D or d: default (equivalent to 2 )\n\n"
                "0 or F or f: no composition-based statistics\n\n"
                "1: Composition-based statistics as in NAR 29:2994-3005, 2001\n\n"
                "2 or T or t : Composition-based score adjustment as in "
                "Bioinformatics 21:902-911, 2005, conditioned on sequence "
                "properties\n\n"
                "3: Composition-based score adjustment as in Bioinformatics "
                "21:902-911, 2005, unconditionally.\n\n"
                "For programs other than tblastn, must either be absent or be "
                "D, F or 0\n\n"
                "Default = 2.",
                equate=False,
            ),
            # Query filtering options:
            _Option(
                ["-seg", "seg"],
                "Filter query sequence with SEG (string).\n\n"
                'Format: "yes", "window locut hicut", or "no" to disable.'
                'Default is "12 2.2 2.5"',
                equate=False,
            ),
            # Extension options:
            _Switch(["-ungapped", "ungapped"], "Perform ungapped alignment only?"),
            _Switch(
                ["-use_sw_tback", "use_sw_tback"],
                "Compute locally optimal Smith-Waterman alignments?",
            ),
        ]
        _NcbiblastMain2SeqCommandline.__init__(self, cmd, **kwargs)


class NcbitblastnCommandline(_NcbiblastMain2SeqCommandline):
    """Wrapper for the NCBI BLAST+ program tblastn.

    With the release of BLAST+ (BLAST rewritten in C++ instead of C), the NCBI
    replaced the old blastall tool with separate tools for each of the searches.
    This wrapper therefore replaces BlastallCommandline with option -p tblastn.

    >>> from Bio.Blast.Applications import NcbitblastnCommandline
    >>> cline = NcbitblastnCommandline(help=True)
    >>> cline
    NcbitblastnCommandline(cmd='tblastn', help=True)
    >>> print(cline)
    tblastn -help

    You would typically run the command line with cline() or via the Python
    subprocess module, as described in the Biopython tutorial.
    """

    def __init__(self, cmd="tblastn", **kwargs):
        """Initialize the class."""
        self.parameters = [
            # General search options:
            _Option(
                ["-task", "task"],
                "Task to execute (string, tblastn (default) or tblastn-fast).",
                checker_function=lambda value: value in ["tblastn", "tblastn-fast"],
                equate=False,
            ),
            _Option(
                ["-db_gencode", "db_gencode"],
                "Genetic code to use to translate query (integer, default 1).",
                equate=False,
            ),
            _Option(
                ["-frame_shift_penalty", "frame_shift_penalty"],
                "Frame shift penalty (integer, at least 1, default ignored) (OBSOLETE).\n\n"
                "This was removed in BLAST 2.2.27+",
                equate=False,
            ),
            _Option(
                ["-max_intron_length", "max_intron_length"],
                "Maximum intron length (integer).\n\n"
                "Length of the largest intron allowed in a translated nucleotide "
                "sequence when linking multiple distinct alignments (a negative "
                "value disables linking). Default zero.",
                equate=False,
            ),
            _Option(
                ["-matrix", "matrix"],
                "Scoring matrix name (default BLOSUM62).",
                equate=False,
            ),
            _Option(
                ["-threshold", "threshold"],
                "Minimum score for words to be added to the BLAST lookup table (float).",
                equate=False,
            ),
            _Option(
                ["-comp_based_stats", "comp_based_stats"],
                "Use composition-based statistics (string, default 2, i.e. True).\n\n"
                "0, F or f: no composition-based statistics\n\n"
                "1: Composition-based statistics as in NAR 29:2994-3005, 2001\n\n"
                "2, T or t, D or d : Composition-based score adjustment as in "
                "Bioinformatics 21:902-911, 2005, conditioned on sequence properties\n\n"
                "3: Composition-based score adjustment as in Bioinformatics 21:902-911, "
                "2005, unconditionally\n\n"
                "Note that only tblastn supports values of 1 and 3.",
                checker_function=lambda value: value in "0Ft12TtDd3",
                equate=False,
            ),
            # Query filtering options:
            _Option(
                ["-seg", "seg"],
                "Filter query sequence with SEG (string).\n\n"
                'Format: "yes", "window locut hicut", or "no" to disable.\n\n'
                'Default is "12 2.2 2.5"',
                equate=False,
            ),
            # Extension options:
            _Switch(["-ungapped", "ungapped"], "Perform ungapped alignment only?"),
            # Miscellaneous options:
            _Switch(
                ["-use_sw_tback", "use_sw_tback"],
                "Compute locally optimal Smith-Waterman alignments?",
            ),
            # PSI-TBLASTN options:
            _Option(
                ["-in_pssm", "in_pssm"],
                "PSI-BLAST checkpoint file.\n\nIncompatible with: remote, query",
                filename=True,
                equate=False,
            ),
        ]
        _NcbiblastMain2SeqCommandline.__init__(self, cmd, **kwargs)


class NcbitblastxCommandline(_NcbiblastMain2SeqCommandline):
    """Wrapper for the NCBI BLAST+ program tblastx.

    With the release of BLAST+ (BLAST rewritten in C++ instead of C), the NCBI
    replaced the old blastall tool with separate tools for each of the searches.
    This wrapper therefore replaces BlastallCommandline with option -p tblastx.

    >>> from Bio.Blast.Applications import NcbitblastxCommandline
    >>> cline = NcbitblastxCommandline(help=True)
    >>> cline
    NcbitblastxCommandline(cmd='tblastx', help=True)
    >>> print(cline)
    tblastx -help

    You would typically run the command line with cline() or via the Python
    subprocess module, as described in the Biopython tutorial.
    """

    def __init__(self, cmd="tblastx", **kwargs):
        """Initialize the class."""
        self.parameters = [
            # Input query options:
            _Option(
                ["-strand", "strand"],
                "Query strand(s) to search against database/subject.\n\n"
                'Values allowed are "both" (default), "minus", "plus".',
                checker_function=lambda value: value in ["both", "minus", "plus"],
                equate=False,
            ),
            # Input query options:
            _Option(
                ["-query_gencode", "query_gencode"],
                "Genetic code to use to translate query (integer, default 1).",
                equate=False,
            ),
            # General search options:
            _Option(
                ["-db_gencode", "db_gencode"],
                "Genetic code to use to translate query (integer, default 1).",
                equate=False,
            ),
            _Option(
                ["-max_intron_length", "max_intron_length"],
                "Maximum intron length (integer).\n\n"
                "Length of the largest intron allowed in a translated nucleotide "
                "sequence when linking multiple distinct alignments (a negative "
                "value disables linking). Default zero.",
                equate=False,
            ),
            _Option(
                ["-matrix", "matrix"],
                "Scoring matrix name (default BLOSUM62).",
                equate=False,
            ),
            _Option(
                ["-threshold", "threshold"],
                "Minimum score for words to be added to the BLAST lookup table (float).",
                equate=False,
            ),
            # Query filtering options:
            _Option(
                ["-seg", "seg"],
                "Filter query sequence with SEG (string).\n\n"
                'Format: "yes", "window locut hicut", or "no" to disable.\n\n'
                'Default is "12 2.2 2.5"',
                equate=False,
            ),
        ]
        _NcbiblastMain2SeqCommandline.__init__(self, cmd, **kwargs)


class NcbipsiblastCommandline(_Ncbiblast2SeqCommandline):
    """Wrapper for the NCBI BLAST+ program psiblast.

    With the release of BLAST+ (BLAST rewritten in C++ instead of C), the NCBI
    replaced the old blastpgp tool with a similar tool psiblast. This wrapper
    therefore replaces BlastpgpCommandline, the wrapper for blastpgp.

    >>> from Bio.Blast.Applications import NcbipsiblastCommandline
    >>> cline = NcbipsiblastCommandline(help=True)
    >>> cline
    NcbipsiblastCommandline(cmd='psiblast', help=True)
    >>> print(cline)
    psiblast -help

    You would typically run the command line with cline() or via the Python
    subprocess module, as described in the Biopython tutorial.
    """

    def __init__(self, cmd="psiblast", **kwargs):
        """Initialize the class."""
        self.parameters = [
            # General search options:
            _Option(
                ["-matrix", "matrix"],
                "Scoring matrix name (default BLOSUM62).",
                equate=False,
            ),
            _Option(
                ["-threshold", "threshold"],
                "Minimum score for words to be added to the BLAST lookup table (float).",
                equate=False,
            ),
            _Option(
                ["-comp_based_stats", "comp_based_stats"],
                "Use composition-based statistics (string, default 2, i.e. True).\n\n"
                "0, F or f: no composition-based statistics\n\n"
                "2, T or t, D or d : Composition-based score adjustment as in "
                "Bioinformatics 21:902-911, 2005, conditioned on sequence properties\n\n"
                "Note that tblastn also supports values of 1 and 3.",
                checker_function=lambda value: value in "0Ft2TtDd",
                equate=False,
            ),
            # Query filtering options:
            _Option(
                ["-seg", "seg"],
                "Filter query sequence with SEG (string).\n\n"
                'Format: "yes", "window locut hicut", or "no" to disable. '
                'Default is "12 2.2 2.5"',
                equate=False,
            ),
            # Extension options:
            _Option(
                ["-gap_trigger", "gap_trigger"],
                "Number of bits to trigger gapping (float, default 22).",
                equate=False,
            ),
            # Miscellaneous options:
            _Switch(
                ["-use_sw_tback", "use_sw_tback"],
                "Compute locally optimal Smith-Waterman alignments?",
            ),
            # PSI-BLAST options:
            _Option(
                ["-num_iterations", "num_iterations"],
                "Number of iterations to perform (integer, at least one).\n\n"
                "Default is one. Incompatible with: remote",
                equate=False,
            ),
            _Option(
                ["-out_pssm", "out_pssm"],
                "File name to store checkpoint file.",
                filename=True,
                equate=False,
            ),
            _Option(
                ["-out_ascii_pssm", "out_ascii_pssm"],
                "File name to store ASCII version of PSSM.",
                filename=True,
                equate=False,
            ),
            _Switch(
                ["-save_pssm_after_last_round", "save_pssm_after_last_round"],
                "Save PSSM after the last database search.",
            ),
            _Switch(
                ["-save_each_pssm", "save_each_pssm"],
                "Save PSSM after each iteration\n\n"
                "File name is given in -save_pssm or -save_ascii_pssm options.",
            ),
            _Option(
                ["-in_msa", "in_msa"],
                "File name of multiple sequence alignment to restart PSI-BLAST.\n\n"
                "Incompatible with: in_pssm, query",
                filename=True,
                equate=False,
            ),
            _Option(
                ["-msa_master_idx", "msa_master_idx"],
                "Index of sequence to use as master in MSA.\n\n"
                "Index (1-based) of sequence to use as the master in the multiple "
                "sequence alignment. If not specified, the first sequence is used.",
                equate=False,
            ),
            _Option(
                ["-in_pssm", "in_pssm"],
                "PSI-BLAST checkpoint file.\n\n"
                "Incompatible with: in_msa, query, phi_pattern",
                filename=True,
                equate=False,
            ),
            # PSSM engine options:
            _Option(
                ["-pseudocount", "pseudocount"],
                "Pseudo-count value used when constructing PSSM.\n\n"
                "Integer. Default is zero.",
                equate=False,
            ),
            _Option(
                ["-inclusion_ethresh", "inclusion_ethresh"],
                "E-value inclusion threshold for pairwise alignments (float, default 0.002).",
                equate=False,
            ),
            _Switch(
                ["-ignore_msa_master", "ignore_msa_master"],
                "Ignore the master sequence when creating PSSM.\n\n"
                "Requires: in_msa\n"
                "Incompatible with: msa_master_idx, in_pssm, query, query_loc, "
                "phi_pattern",
            ),
            # PHI-BLAST options:
            _Option(
                ["-phi_pattern", "phi_pattern"],
                "File name containing pattern to search.\n\n"
                "Incompatible with: in_pssm",
                filename=True,
                equate=False,
            ),
        ]
        _Ncbiblast2SeqCommandline.__init__(self, cmd, **kwargs)

    def _validate(self):
        incompatibles = {
            "num_iterations": ["remote"],
            "in_msa": ["in_pssm", "query"],
            "in_pssm": ["in_msa", "query", "phi_pattern"],
            "ignore_msa_master": [
                "msa_master_idx",
                "in_pssm",
                "query",
                "query_loc",
                "phi_pattern",
            ],
        }
        self._validate_incompatibilities(incompatibles)
        _Ncbiblast2SeqCommandline._validate(self)


class NcbirpsblastCommandline(_NcbiblastCommandline):
    """Wrapper for the NCBI BLAST+ program rpsblast.

    With the release of BLAST+ (BLAST rewritten in C++ instead of C), the NCBI
    replaced the old rpsblast tool with a similar tool of the same name. This
    wrapper replaces RpsBlastCommandline, the wrapper for the old rpsblast.

    >>> from Bio.Blast.Applications import NcbirpsblastCommandline
    >>> cline = NcbirpsblastCommandline(help=True)
    >>> cline
    NcbirpsblastCommandline(cmd='rpsblast', help=True)
    >>> print(cline)
    rpsblast -help

    You would typically run the command line with cline() or via the Python
    subprocess module, as described in the Biopython tutorial.
    """

    def __init__(self, cmd="rpsblast", **kwargs):
        """Initialize the class."""
        # TODO - remove the -word_size argument as per BLAST+ 2.2.30
        # (BLAST team say it should never have been included, since
        # the word size is set when building the domain database.)
        # This likely means reviewing the class hierarchy again.
        self.parameters = [
            # Query filtering options:
            _Option(
                ["-seg", "seg"],
                "Filter query sequence with SEG (string).\n\n"
                'Format: "yes", "window locut hicut", or "no" to disable.'
                'Default is "12 2.2 2.5"',
                equate=False,
            ),
            # Restrict search or results:
            _Option(
                ["-culling_limit", "culling_limit"],
                "Hit culling limit (integer).\n\n"
                "If the query range of a hit is enveloped by that of at "
                "least this many higher-scoring hits, delete the hit. "
                "Incompatible with: best_hit_overhang, best_hit_score_edge.",
                equate=False,
            ),
            _Option(
                ["-best_hit_overhang", "best_hit_overhang"],
                "Best Hit algorithm overhang value (recommended value: 0.1).\n\n"
                "Float between 0.0 and 0.5 inclusive. "
                "Incompatible with: culling_limit.",
                equate=False,
            ),
            _Option(
                ["-best_hit_score_edge", "best_hit_score_edge"],
                "Best Hit algorithm score edge value (recommended value: 0.1).\n\n"
                "Float between 0.0 and 0.5 inclusive. "
                "Incompatible with: culling_limit.",
                equate=False,
            ),
            # General search options:
            _Option(
                ["-comp_based_stats", "comp_based_stats"],
                "Use composition-based statistics.\n\n"
                "D or d: default (equivalent to 0)\n\n"
                "0 or F or f: Simplified Composition-based statistics as in "
                "Bioinformatics 15:1000-1011, 1999\n\n"
                "1 or T or t: Composition-based statistics as in NAR 29:2994-3005, "
                "2001\n\n"
                "Default = 0.",
                checker_function=lambda value: value in "Dd0Ff1Tt",
                equate=False,
            ),
            # Misc options:
            _Switch(
                ["-use_sw_tback", "use_sw_tback"],
                "Compute locally optimal Smith-Waterman alignments?",
            ),
        ]
        _NcbiblastCommandline.__init__(self, cmd, **kwargs)

    def _validate(self):
        incompatibles = {"culling_limit": ["best_hit_overhang", "best_hit_score_edge"]}
        self._validate_incompatibilities(incompatibles)
        _NcbiblastCommandline._validate(self)


class NcbirpstblastnCommandline(_NcbiblastCommandline):
    """Wrapper for the NCBI BLAST+ program rpstblastn.

    With the release of BLAST+ (BLAST rewritten in C++ instead of C), the NCBI
    replaced the old rpsblast tool with a similar tool of the same name, and a
    separate tool rpstblastn for Translated Reverse Position Specific BLAST.

    >>> from Bio.Blast.Applications import NcbirpstblastnCommandline
    >>> cline = NcbirpstblastnCommandline(help=True)
    >>> cline
    NcbirpstblastnCommandline(cmd='rpstblastn', help=True)
    >>> print(cline)
    rpstblastn -help

    You would typically run the command line with cline() or via the Python
    subprocess module, as described in the Biopython tutorial.
    """

    def __init__(self, cmd="rpstblastn", **kwargs):
        """Initialize the class."""
        # TODO - remove the -word_size argument as per BLAST+ 2.2.30
        # (BLAST team say it should never have been included, since
        # the word size is set when building the domain database.)
        # This likely means reviewing the class hierarchy again.
        self.parameters = [
            # Input query options:
            _Option(
                ["-strand", "strand"],
                "Query strand(s) to search against database/subject.\n\n"
                'Values allowed are "both" (default), "minus", "plus".',
                checker_function=lambda value: value in ["both", "minus", "plus"],
                equate=False,
            ),
            # Input query options:
            _Option(
                ["-query_gencode", "query_gencode"],
                "Genetic code to use to translate query (integer, default 1).",
                equate=False,
            ),
            # Query filtering options:
            _Option(
                ["-seg", "seg"],
                "Filter query sequence with SEG (string).\n\n"
                'Format: "yes", "window locut hicut", or "no" to disable. '
                'Default is "12 2.2 2.5"',
                equate=False,
            ),
            # General search options:
            _Option(
                ["-comp_based_stats", "comp_based_stats"],
                "Use composition-based statistics.\n\n"
                "D or d: default (equivalent to 0)\n\n"
                "0 or F or f: Simplified Composition-based statistics as in "
                "Bioinformatics 15:1000-1011, 1999\n\n"
                "1 or T or t: Composition-based statistics as in NAR 29:2994-3005, "
                "2001\n\n"
                "Default = 0.",
                checker_function=lambda value: value in "Dd0Ff1Tt",
                equate=False,
            ),
            # Extension options:
            _Switch(["-ungapped", "ungapped"], "Perform ungapped alignment only?"),
            # Miscellaneous options:
            _Switch(
                ["-use_sw_tback", "use_sw_tback"],
                "Compute locally optimal Smith-Waterman alignments?",
            ),
        ]
        _NcbiblastCommandline.__init__(self, cmd, **kwargs)


class NcbiblastformatterCommandline(_NcbibaseblastCommandline):
    """Wrapper for the NCBI BLAST+ program blast_formatter.

    With the release of BLAST 2.2.24+ (i.e. the BLAST suite rewritten in C++
    instead of C), the NCBI added the ASN.1 output format option to all the
    search tools, and extended the blast_formatter to support this as input.

    The blast_formatter command allows you to convert the ASN.1 output into
    the other output formats (XML, tabular, plain text, HTML).

    >>> from Bio.Blast.Applications import NcbiblastformatterCommandline
    >>> cline = NcbiblastformatterCommandline(archive="example.asn", outfmt=5, out="example.xml")
    >>> cline
    NcbiblastformatterCommandline(cmd='blast_formatter', out='example.xml', outfmt=5, archive='example.asn')
    >>> print(cline)
    blast_formatter -out example.xml -outfmt 5 -archive example.asn

    You would typically run the command line with cline() or via the Python
    subprocess module, as described in the Biopython tutorial.

    Note that this wrapper is for the version of blast_formatter from BLAST
    2.2.24+ (or later) which is when the NCBI first announced the inclusion
    this tool. There was actually an early version in BLAST 2.2.23+ (and
    possibly in older releases) but this did not have the -archive option
    (instead -rid is a mandatory argument), and is not supported by this
    wrapper.
    """

    def __init__(self, cmd="blast_formatter", **kwargs):
        """Initialize the class."""
        self.parameters = [
            # Input options
            _Option(
                ["-rid", "rid"],
                "BLAST Request ID (RID), not compatible with archive arg.",
                equate=False,
            ),
            _Option(
                ["-archive", "archive"],
                "Archive file of results, not compatible with rid arg.",
                filename=True,
                equate=False,
            ),
            # Restrict search or results
            _Option(
                ["-max_target_seqs", "max_target_seqs"],
                "Maximum number of aligned sequences to keep.",
                checker_function=lambda value: value >= 1,
                equate=False,
            ),
        ]
        _NcbibaseblastCommandline.__init__(self, cmd, **kwargs)

    def _validate(self):
        incompatibles = {"rid": ["archive"]}
        self._validate_incompatibilities(incompatibles)
        _NcbibaseblastCommandline._validate(self)


class NcbideltablastCommandline(_Ncbiblast2SeqCommandline):
    """Create a commandline for the NCBI BLAST+ program deltablast (for proteins).

    This is a wrapper for the deltablast command line command included in
    the NCBI BLAST+ software (not present in the original BLAST).

    >>> from Bio.Blast.Applications import NcbideltablastCommandline
    >>> cline = NcbideltablastCommandline(query="rosemary.pro", db="nr",
    ...                               evalue=0.001, remote=True)
    >>> cline
    NcbideltablastCommandline(cmd='deltablast', query='rosemary.pro', db='nr', evalue=0.001, remote=True)
    >>> print(cline)
    deltablast -query rosemary.pro -db nr -evalue 0.001 -remote

    You would typically run the command line with cline() or via the Python
    subprocess module, as described in the Biopython tutorial.
    """

    def __init__(self, cmd="deltablast", **kwargs):
        """Initialize the class."""
        self.parameters = [
            # General search options:
            _Option(["-matrix", "matrix"], "Scoring matrix name (default BLOSUM62)."),
            _Option(
                ["-threshold", "threshold"],
                "Minimum score for words to be added to the BLAST lookup table (float).",
                equate=False,
            ),
            _Option(
                ["-comp_based_stats", "comp_based_stats"],
                "Use composition-based statistics (string, default 2, i.e. True).\n\n"
                "0, F or f: no composition-based statistics.\n\n"
                "2, T or t, D or d : Composition-based score adjustment as in "
                "Bioinformatics 21:902-911, 2005, conditioned on sequence properties\n\n"
                "Note that tblastn also supports values of 1 and 3.",
                checker_function=lambda value: value in "0Ft2TtDd",
                equate=False,
            ),
            # Query filtering options:
            _Option(
                ["-seg", "seg"],
                "Filter query sequence with SEG (string).\n\n"
                'Format: "yes", "window locut hicut", or "no" to disable. '
                'Default is "12 2.2 2.5"',
                equate=False,
            ),
            # Extension options:
            _Option(
                ["-gap_trigger", "gap_trigger"],
                "Number of bits to trigger gapping. Default = 22.",
                equate=False,
            ),
            # Miscellaneous options:
            _Switch(
                ["-use_sw_tback", "use_sw_tback"],
                "Compute locally optimal Smith-Waterman alignments?",
            ),
            # PSI-BLAST options
            _Option(
                ["-num_iterations", "num_iterations"],
                "Number of iterations to perform. (integer >=1, Default is 1).\n\n"
                "Incompatible with: remote",
                equate=False,
            ),
            _Option(
                ["-out_pssm", "out_pssm"],
                "File name to store checkpoint file.",
                filename=True,
                equate=False,
            ),
            _Option(
                ["-out_ascii_pssm", "out_ascii_pssm"],
                "File name to store ASCII version of PSSM.",
                filename=True,
                equate=False,
            ),
            _Switch(
                ["-save_pssm_after_last_round", "save_pssm_after_last_round"],
                "Save PSSM after the last database search.",
            ),
            _Switch(
                ["-save_each_pssm", "save_each_pssm"],
                "Save PSSM after each iteration.\n\n"
                "File name is given in -save_pssm or -save_ascii_pssm options.",
            ),
            # PSSM engine options
            _Option(
                ["-pseudocount", "pseudocount"],
                "Pseudo-count value used when constructing PSSM (integer, default 0).",
                equate=False,
            ),
            _Option(
                ["-domain_inclusion_ethresh", "domain_inclusion_ethresh"],
                "E-value inclusion threshold for alignments with conserved domains.\n\n"
                "(float, Default is 0.05)",
                equate=False,
            ),
            _Option(
                ["-inclusion_ethresh", "inclusion_ethresh"],
                "Pairwise alignment e-value inclusion threshold (float, default 0.002).",
                equate=False,
            ),
            # DELTA-BLAST options
            _Option(
                ["-rpsdb", "rpsdb"],
                "BLAST domain database name (dtring, Default = 'cdd_delta').",
                equate=False,
            ),
            _Switch(
                ["-show_domain_hits", "show_domain_hits"],
                "Show domain hits?\n\nIncompatible with: remote, subject",
            ),
        ]
        _Ncbiblast2SeqCommandline.__init__(self, cmd, **kwargs)


class NcbimakeblastdbCommandline(AbstractCommandline):
    """Wrapper for the NCBI BLAST+ program makeblastdb.

    This is a wrapper for the NCBI BLAST+ makeblastdb application
    to create BLAST databases. By default, this creates a blast database
    with the same name as the input file. The default output location
    is the same directory as the input.

    >>> from Bio.Blast.Applications import NcbimakeblastdbCommandline
    >>> cline = NcbimakeblastdbCommandline(dbtype="prot",
    ...                                    input_file="NC_005816.faa")
    >>> cline
    NcbimakeblastdbCommandline(cmd='makeblastdb', dbtype='prot', input_file='NC_005816.faa')
    >>> print(cline)
    makeblastdb -dbtype prot -in NC_005816.faa

    You would typically run the command line with cline() or via the Python
    subprocess module, as described in the Biopython tutorial.
    """

    def __init__(self, cmd="makeblastdb", **kwargs):
        """Initialize the class."""
        self.parameters = [
            # Basic input options
            _Switch(
                ["-h", "h"], "Print USAGE and DESCRIPTION; ignore other arguments."
            ),
            _Switch(
                ["-help", "help"],
                "Print USAGE, DESCRIPTION and ARGUMENTS description; "
                "ignore other arguments.",
            ),
            _Switch(
                ["-version", "version"],
                "Print version number;  ignore other arguments.",
            ),
            # Output configuration options
            _Option(
                ["-out", "out"],
                "Output file for alignment.",
                filename=True,
                equate=False,
            ),
            # makeblastdb specific options
            _Option(
                ["-blastdb_version", "blastdb_version"],
                "Version of BLAST database to be created. "
                "Tip: use BLAST database version 4 on 32 bit CPU. "
                "Default = 5",
                equate=False,
                checker_function=lambda x: x == 4 or x == 5,
            ),
            _Option(
                ["-dbtype", "dbtype"],
                "Molecule type of target db ('nucl' or 'prot').",
                equate=False,
                is_required=True,
                checker_function=lambda x: x == "nucl" or x == "prot",
            ),
            _Option(
                ["-in", "input_file"],
                "Input file/database name.",
                filename=True,
                equate=False,
            ),
            _Option(
                ["-input_type", "input_type"],
                "Type of the data specified in input_file.\n\n"
                "Default = 'fasta'. Added in BLAST 2.2.26.",
                filename=False,
                equate=False,
                checker_function=self._input_type_checker,
            ),
            _Option(
                ["-title", "title"],
                "Title for BLAST database.",
                filename=False,
                equate=False,
            ),
            _Switch(
                ["-parse_seqids", "parse_seqids"],
                "Option to parse seqid for FASTA input if set.\n\n"
                "For all other input types, seqids are parsed automatically",
            ),
            _Switch(
                ["-hash_index", "hash_index"], "Create index of sequence hash values."
            ),
            _Option(
                ["-mask_data", "mask_data"],
                "Comma-separated list of input files containing masking "
                "data as produced by NCBI masking applications "
                "(e.g. dustmasker, segmasker, windowmasker).",
                filename=True,
                equate=False,
            ),
            _Option(
                ["-mask_id", "mask_id"],
                "Comma-separated list of strings to uniquely identify the "
                "masking algorithm.",
                filename=False,
                equate=False,
            ),
            _Option(
                ["-mask_desc", "mask_desc"],
                "Comma-separated list of free form strings to describe "
                "the masking algorithm details.",
                filename=False,
                equate=False,
            ),
            _Switch(["-gi_mask", "gi_mask"], "Create GI indexed masking data."),
            _Option(
                ["-gi_mask_name", "gi_mask_name"],
                "Comma-separated list of masking data output files.",
                filename=False,
                equate=False,
            ),
            _Option(
                ["-max_file_sz", "max_file_sz"],
                "Maximum file size for BLAST database files. Default = '1GB'.",
                filename=False,
                equate=False,
            ),
            _Option(
                ["-logfile", "logfile"],
                "File to which the program log should be redirected.",
                filename=True,
                equate=False,
            ),
            _Option(
                ["-taxid", "taxid"],
                "Taxonomy ID to assign to all sequences.",
                filename=False,
                equate=False,
                checker_function=lambda x: type(x)(int(x)) == x,
            ),
            _Option(
                ["-taxid_map", "taxid_map"],
                "Text file mapping sequence IDs to taxonomy IDs.\n\n"
                "Format:<SequenceId> <TaxonomyId><newline>",
                filename=True,
                equate=False,
            ),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)

    def _input_type_checker(self, command):
        return command in ("asn1_bin", "asn1_txt", "blastdb", "fasta")

    def _validate(self):
        incompatibles = {
            "mask_id": ["gi_mask"],
            "gi_mask": ["mask_id"],
            "taxid": ["taxid_map"],
        }

        # Copied from _NcbibaseblastCommandline class above.
        # Code repeated here for python2 and 3 compatibility,
        # because this is not a _NcbibaseblastCommandline subclass.
        for a in incompatibles:
            if self._get_parameter(a):
                for b in incompatibles[a]:
                    if self._get_parameter(b):
                        raise ValueError("Options %s and %s are incompatible." % (a, b))

        if self.mask_id and not self.mask_data:
            raise ValueError("Option mask_id requires mask_data to be set.")
        if self.mask_desc and not self.mask_id:
            raise ValueError("Option mask_desc requires mask_id to be set.")
        if self.gi_mask and not self.parse_seqids:
            raise ValueError("Option gi_mask requires parse_seqids to be set.")
        if self.gi_mask_name and not (self.mask_data and self.gi_mask):
            raise ValueError(
                "Option gi_mask_name requires mask_data and gi_mask to be set."
            )
        if self.taxid_map and not self.parse_seqids:
            raise ValueError("Option taxid_map requires parse_seqids to be set.")
        AbstractCommandline._validate(self)


def _test():
    """Run the Bio.Blast.Applications module's doctests (PRIVATE)."""
    import doctest

    doctest.testmod(verbose=1)


if __name__ == "__main__":
    # Run the doctests
    _test()
