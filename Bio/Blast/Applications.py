# Copyright 2001 Brad Chapman.
# Revisions copyright 2009-2010 by Peter Cock.
# Revisions copyright 2010 by Phillip Garland.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Definitions for interacting with BLAST related applications.

Obsolete wrappers for the old/classic NCBI BLAST tools (written in C):

- FastacmdCommandline
- BlastallCommandline
- BlastpgpCommandline
- RpsBlastCommandline

Wrappers for the new NCBI BLAST+ tools (written in C++):

- NcbiblastpCommandline - Protein-Protein BLAST
- NcbiblastnCommandline - Nucleotide-Nucleotide BLAST
- NcbiblastxCommandline - Translated Query-Protein Subject BLAST
- NcbitblastnCommandline - Protein Query-Translated Subject BLAST
- NcbitblastxCommandline - Translated Query-Protein Subject BLAST
- NcbipsiblastCommandline - Position-Specific Initiated BLAST
- NcbirpsblastCommandline - Reverse Position Specific BLAST
- NcbirpstblastnCommandline - Translated Reverse Position Specific BLAST
- NcbiblastformatterCommandline - Convert ASN.1 to other BLAST output formats

For further details, see:

Camacho et al. BLAST+: architecture and applications
BMC Bioinformatics 2009, 10:421
doi:10.1186/1471-2105-10-421
"""
from Bio.Application import _Option, AbstractCommandline, _Switch

class FastacmdCommandline(AbstractCommandline):
    """Create a commandline for the fasta program from NCBI (OBSOLETE).

    """
    def __init__(self, cmd="fastacmd", **kwargs):
        self.parameters = [
           _Option(["-d", "database"],
                   "The database to retrieve from.",
                   is_required=True,
                   equate=False),
           _Option(["-s", "search_string"],
                   "The id to search for.",
                   is_required=True,
                   equate=False)
          ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


class _BlastCommandLine(AbstractCommandline):
    """Base Commandline object for (classic) NCBI BLAST wrappers (PRIVATE).

    This is provided for subclassing, it deals with shared options
    common to all the BLAST tools (blastall, rpsblast, blastpgp).
    """
    def __init__(self, cmd=None, **kwargs):
        assert cmd is not None
        extra_parameters = [
           _Switch(["--help", "help"],
                    "Print USAGE, DESCRIPTION and ARGUMENTS description;  ignore other arguments."),
           _Option(["-d", "database"],
                   "The database to BLAST against.",
                   is_required=True,
                   equate=False),
           _Option(["-i", "infile"],
                   "The sequence to search with.",
                   filename=True,
                   is_required=True,
                   equate=False),
           _Option(["-e", "expectation"], 
                   "Expectation value cutoff.",
                   equate=False),
           _Option(["-m", "align_view"], 
                   "Alignment view.  Integer 0-11.  Use 7 for XML output.",
                   equate=False),
           _Option(["-o", "align_outfile", "outfile"],
                   "Output file for alignment.",
                   filename=True,
                   equate=False),
           _Option(["-y", "xdrop_extension"], 
                   "Dropoff for blast extensions.",
                   equate=False),
           _Option(["-F", "filter"],
                   "Filter query sequence with SEG?  T/F",
                   equate=False),
           _Option(["-X", "xdrop"], 
                   "Dropoff value (bits) for gapped alignments.",
                   equate=False),
           _Option(["-I", "show_gi"], 
                   "Show GI's in deflines?  T/F",
                   equate=False),
           _Option(["-J", "believe_query"], 
                   "Believe the query defline?  T/F",
                   equate=False),
           _Option(["-Z", "xdrop_final"], 
                   "X dropoff for final gapped alignment.",
                   equate=False),
           _Option(["-z", "db_length"], 
                   "Effective database length.",
                   equate=False),
           _Option(["-O", "seqalign_file"],
                   "seqalign file to output.",
                   filename=True,
                   equate=False),
           _Option(["-v", "descriptions"], 
                   "Number of one-line descriptions.",
                   equate=False),
           _Option(["-b", "alignments"], 
                   "Number of alignments.",
                   equate=False),
           _Option(["-Y", "search_length"], 
                   "Effective length of search space (use zero for the "
                   "real size).",
                   equate=False),
           _Option(["-T", "html"], 
                   "Produce HTML output?  T/F",
                   equate=False),
           _Option(["-U", "case_filter"],
                   "Use lower case filtering of FASTA sequence? T/F",
                   equate=False),
           _Option(["-a", "nprocessors"],
                   "Number of processors to use.",
                   equate=False),
           _Option(["-g", "gapped"], 
                   "Whether to do a gapped alignment.  T/F",
                   equate=False),
        ]
        try:
            #Insert extra parameters - at the start just in case there
            #are any arguments which must come last:
            self.parameters = extra_parameters + self.parameters
        except AttributeError:
            #Should we raise an error?  The subclass should have set this up!
            self.parameters = extra_parameters
        AbstractCommandline.__init__(self, cmd, **kwargs)

    def _validate(self):
        if self.help:
            #Don't want to check the normally mandatory arguments like db
            return
        AbstractCommandline._validate(self)


class _BlastAllOrPgpCommandLine(_BlastCommandLine):
    """Base Commandline object for NCBI BLAST wrappers (PRIVATE).

    This is provided for subclassing, it deals with shared options
    common to all the blastall and blastpgp tools (but not rpsblast).
    """
    def __init__(self, cmd=None, **kwargs):
        assert cmd is not None
        extra_parameters = [
           _Option(["-G", "gap_open"], 
                   "Gap open penalty",
                   equate=False),
           _Option(["-E", "gap_extend"], 
                    "Gap extension penalty",
                   equate=False),
           _Option(["-A", "window_size"],
                    "Multiple hits window size",
                   equate=False),
           _Option(["-f", "hit_extend"], 
                   "Threshold for extending hits.",
                   equate=False),
           _Option(["-K", "keep_hits"],
                   " Number of best hits from a region to keep.",
                   equate=False),
           _Option(["-W", "wordsize"], 
                   "Word size",
                   equate=False),
           _Option(["-P", "passes"],
                   "Hits/passes.  Integer 0-2. 0 for multiple hit, "
                   "1 for single hit (does not apply to blastn)",
                   equate=False),
        ]
        try:
            #Insert extra parameters - at the start just in case there
            #are any arguments which must come last:
            self.parameters = extra_parameters + self.parameters
        except AttributeError:
            #Should we raise an error?  The subclass should have set this up!
            self.parameters = extra_parameters
        _BlastCommandLine.__init__(self, cmd, **kwargs)


class BlastallCommandline(_BlastAllOrPgpCommandLine):
    """Create a commandline for the blastall program from NCBI (OBSOLETE).

    With the release of BLAST+ (BLAST rewritten in C++ instead of C), the NCBI
    are replacing blastall with separate tools blastn, blastp, blastx, tblastn
    and tblastx.

    Like blastall, this wrapper is now obsolete, and will be deprecated and
    removed in a future release of Biopython.

    >>> from Bio.Blast.Applications import BlastallCommandline
    >>> cline = BlastallCommandline(program="blastx", infile="m_cold.fasta",
    ...                             database="nr", expectation=0.001)
    >>> cline
    BlastallCommandline(cmd='blastall', database='nr', infile='m_cold.fasta', expectation=0.001, program='blastx')
    >>> print cline
    blastall -d nr -i m_cold.fasta -e 0.001 -p blastx

    You would typically run the command line with cline() or via the Python
    subprocess module, as described in the Biopython tutorial.
    """
    #TODO - This could use more checking for valid parameters to the program.
    def __init__(self, cmd="blastall",**kwargs):
        import warnings
        warnings.warn("Like blastall, this wrapper is now obsolete, and will be deprecated and removed in a future release of Biopython.", PendingDeprecationWarning)
        self.parameters = [
            #Sorted in the same order as the output from blastall --help
            #which should make it easier to keep them up to date in future.
            #Note that some arguments are defined the the base clases (above).
           _Option(["-p", "program"],
                   "The blast program to use (e.g. blastp, blastn).",
                   is_required=True,
                   equate=False),
           _Option(["-q", "nuc_mismatch"], 
                   "Penalty for a nucleotide mismatch (blastn only).",
                   equate=False),
           _Option(["-r", "nuc_match"], 
                   "Reward for a nucleotide match (blastn only).",
                   equate=False),
           _Option(["-Q", "query_genetic_code"],
                   "Query Genetic code to use.",
                   equate=False),
           _Option(["-D", "db_genetic_code"],
                   "DB Genetic code (for tblast[nx] only).",
                   equate=False),
           _Option(["-M", "matrix"], 
                   "Matrix to use",
                   equate=False),
           _Option(["-S", "strands"], 
                   "Query strands to search against database (for blast[nx], "
                   "and tblastx). 3 is both, 1 is top, 2 is bottom.",
                   equate=False),
           _Option(["-l", "restrict_gi"],
                   "Restrict search of database to list of GI's.",
                   equate=False),
           _Option(["-R", "checkpoint"],
                   "PSI-TBLASTN checkpoint input file.",
                   filename=True,
                   equate=False),
           _Option(["-n", "megablast"],
                   "MegaBlast search T/F.",
                   equate=False),
           #The old name "region_length" is for consistency with our
           #old blastall function wrapper:
           _Option(["-L", "region_length", "range_restriction"],
                   """Location on query sequence (string format start,end).

                   In older versions of BLAST, -L set the length of region
                   used to judge hits (see -K parameter).""",
                   equate=False),
           _Option(["-w", "frame_shit_penalty"],
                   "Frame shift penalty (OOF algorithm for blastx).",
                   equate=False),
           _Option(["-t", "largest_intron"],
                   "Length of the largest intron allowed in a translated "
                   "nucleotide sequence when linking multiple distinct "
                   "alignments. (0 invokes default behavior; a negative value "
                   "disables linking.)",
                   equate=False),
           _Option(["-B", "num_concatenated_queries"],
                   "Number of concatenated queries, for blastn and tblastn.",
                   equate=False),
           _Option(["-V", "oldengine"],
                   "Force use of the legacy BLAST engine.",
                   equate=False),
           _Option(["-C", "composition_based"],
                   """Use composition-based statistics for tblastn:
                   D or d: default (equivalent to F)
                   0 or F or f: no composition-based statistics
                   1 or T or t: Composition-based statistics as in NAR 29:2994-3005, 2001
                   2: Composition-based score adjustment as in Bioinformatics
                       21:902-911, 2005, conditioned on sequence properties
                   3: Composition-based score adjustment as in Bioinformatics
                       21:902-911, 2005, unconditionally
                   For programs other than tblastn, must either be absent or be
                   D, F or 0.""",
                   equate=False),
           _Option(["-s", "smith_waterman"],
                   "Compute locally optimal Smith-Waterman alignments (This "
                   "option is only available for gapped tblastn.) T/F",
                   equate=False),
        ] 
        _BlastAllOrPgpCommandLine.__init__(self, cmd, **kwargs)


class BlastpgpCommandline(_BlastAllOrPgpCommandLine):
    """Create a commandline for the blastpgp program from NCBI (OBSOLETE).

    With the release of BLAST+ (BLAST rewritten in C++ instead of C), the NCBI
    are replacing blastpgp with a renamed tool psiblast. This module provides
    NcbipsiblastCommandline as a wrapper for the new tool psiblast.
    
    Like blastpgp (and blastall), this wrapper is now obsolete, and will be
    deprecated and removed in a future release of Biopython.

    >>> from Bio.Blast.Applications import BlastpgpCommandline
    >>> cline = BlastpgpCommandline(help=True)
    >>> cline
    BlastpgpCommandline(cmd='blastpgp', help=True)
    >>> print cline
    blastpgp --help

    You would typically run the command line with cline() or via the Python
    subprocess module, as described in the Biopython tutorial.
    """
    def __init__(self, cmd="blastpgp",**kwargs):
        import warnings
        warnings.warn("Like blastpgp (and blastall), this wrapper is now obsolete, and will be deprecated and removed in a future release of Biopython.", PendingDeprecationWarning)
        self.parameters = [
           _Option(["-C", "checkpoint_outfile"],
                   "Output file for PSI-BLAST checkpointing.",
                   filename=True,
                   equate=False),
           _Option(["-R", "restart_infile"],
                   "Input file for PSI-BLAST restart.",
                   filename=True,
                   equate=False),
           _Option(["-k", "hit_infile"],
                   "Hit file for PHI-BLAST.",
                   filename=True,
                   equate=False),
           _Option(["-Q", "matrix_outfile"],
                   "Output file for PSI-BLAST matrix in ASCII.",
                   filename=True,
                   equate=False),
           _Option(["-B", "align_infile"],
                   "Input alignment file for PSI-BLAST restart.",
                   filename=True,
                   equate=False),
           _Option(["-S", "required_start"], 
                   "Start of required region in query.",
                   equate=False),
           _Option(["-H", "required_end"],
                   "End of required region in query.",
                   equate=False),
           _Option(["-j", "npasses"],
                   "Number of passes",
                   equate=False),
           _Option(["-N", "nbits_gapping"], 
                   "Number of bits to trigger gapping.",
                   equate=False),
           _Option(["-c", "pseudocounts"],
                   "Pseudocounts constants for multiple passes.",
                   equate=False),
           _Option(["-h", "model_threshold"], 
                   "E-value threshold to include in multipass model.",
                   equate=False),
           #Does the old name "region_length" for -L make sense?
           _Option(["-L", "region_length"], 
                   "Cost to decline alignment (disabled when zero).",
                   equate=False),
           _Option(["-M", "matrix"], 
                   "Matrix (string, default BLOSUM62).",
                   equate=False),
           _Option(["-p", "program"],
                   "The blast program to use (e.g blastpgp, patseedp or seedp).",
                   is_required=True,
                   equate=False),
        ] 
        _BlastAllOrPgpCommandLine.__init__(self, cmd, **kwargs)


class RpsBlastCommandline(_BlastCommandLine):
    """Create a commandline for the classic rpsblast program from NCBI (OBSOLETE).

    With the release of BLAST+ (BLAST rewritten in C++ instead of C), the NCBI
    are replacing the old rpsblast with a new version of the same name plus a
    second tool rpstblastn, both taking different command line arguments. This
    module provides NcbirpsblastCommandline and NcbirpsblasntCommandline as
    wrappers for the new tools.
    
    Like the old rpsblast (and blastall), this wrapper is now obsolete, and will
    be deprecated and removed in a future release of Biopython.

    >>> from Bio.Blast.Applications import RpsBlastCommandline
    >>> cline = RpsBlastCommandline(help=True)
    >>> cline
    RpsBlastCommandline(cmd='rpsblast', help=True)
    >>> print cline
    rpsblast --help

    You would typically run the command line with cline() or via the Python
    subprocess module, as described in the Biopython tutorial.
    """
    def __init__(self, cmd="rpsblast",**kwargs):
        import warnings
        warnings.warn("Like the old rpsblast (and blastall), this wrapper is now obsolete, and will be deprecated and removed in a future release of Biopython.", PendingDeprecationWarning)
        self.parameters = [
           #Note -N is also in blastpgp, but not blastall
           _Option(["-N", "nbits_gapping"], 
                   "Number of bits to trigger gapping.",
                   equate=False),
           #Note blastall and blastpgp wrappers have -P with name "passes".
           #If this is the same thing, we should be consistent!
           _Option(["-P", "multihit"],
                   "0 for multiple hit, 1 for single hit",
                   equate=False),
           _Option(["-l", "logfile"],
                   "Logfile name.",
                   filename=True,
                   equate=False),
           _Option(["-p", "protein"], 
                   "Query sequence is protein. T/F",
                   equate=False),
           _Option(["-L", "range_restriction"], 
                   "Location on query sequence (string format start,end).",
                   equate=False),
        ] 
        _BlastCommandLine.__init__(self, cmd, **kwargs)

##############################################################################
# Legacy BLAST wrappers above, (new) BLAST+ wrappers below
##############################################################################

class _NcbibaseblastCommandline(AbstractCommandline):
    """Base Commandline object for (new) NCBI BLAST+ wrappers (PRIVATE).

    This is provided for subclassing, it deals with shared options
    common to all the BLAST tools (blastn, rpsblast, rpsblast, etc
    AND blast_formatter).
    """
    def __init__(self, cmd=None, **kwargs):
        assert cmd is not None
        extra_parameters = [
            #Core:
            _Switch(["-h", "h"],
                    "Print USAGE and DESCRIPTION;  ignore other arguments."),
            _Switch(["-help", "help"],
                    "Print USAGE, DESCRIPTION and ARGUMENTS description; "
                    "ignore other arguments."),
            _Switch(["-version", "version"],
                    "Print version number;  ignore other arguments."),
            # Output configuration options
            _Option(["-out", "out"],
                    "Output file for alignment.",
                    filename=True,
                    equate=False),
            #Formatting options:
            _Option(["-outfmt", "outfmt"], 
                    "Alignment view.  Integer 0-11.  Use 5 for XML output "
                    "(differs from classic BLAST which used 7 for XML).",
                    equate=False),
                    #TODO - Document and test the column options
            _Switch(["-show_gis","show_gis"],
                    "Show NCBI GIs in deflines?"),
            _Option(["-num_descriptions","num_descriptions"],
                    """Number of database sequences to show one-line descriptions for.

                    Integer argument (at least zero). Default is 500.
                    See also num_alignments.""",
                    equate=False),
            _Option(["-num_alignments","num_alignments"],
                    """Number of database sequences to show num_alignments for.

                    Integer argument (at least zero). Default is 200.
                    See also num_alignments.""",
                    equate=False),
            _Switch(["-html", "html"],
                    "Produce HTML output? See also the outfmt option."),
            #Miscellaneous options
            _Switch(["-parse_deflines", "parse_deflines"],
                    "Should the query and subject defline(s) be parsed?"),
            ]
        try:
            #Insert extra parameters - at the start just in case there
            #are any arguments which must come last:
            self.parameters = extra_parameters + self.parameters
        except AttributeError:
            #Should we raise an error?  The subclass should have set this up!
            self.parameters = extra_parameters
        AbstractCommandline.__init__(self, cmd, **kwargs)

    def _validate_incompatibilities(self, incompatibles):
        """Used by the BLAST+ _validate method (PRIVATE)."""
        for a in incompatibles:
            if self._get_parameter(a):
                for b in incompatibles[a]:
                    if self._get_parameter(b):
                        raise ValueError("Options %s and %s are incompatible." \
                                         % (a,b))

        
class _NcbiblastCommandline(_NcbibaseblastCommandline):
    """Base Commandline object for (new) NCBI BLAST+ wrappers (PRIVATE).

    This is provided for subclassing, it deals with shared options
    common to all the BLAST tools (blastn, rpsblast, rpsblast, etc).
    """
    def __init__(self, cmd=None, **kwargs):
        assert cmd is not None
        extra_parameters = [
            #Input query options:
            _Option(["-query", "query"],
                    "The sequence to search with.",
                    filename=True,
                    equate=False), #Should this be required?
            _Option(["-query_loc", "query_loc"],
                    "Location on the query sequence (Format: start-stop)",
                    equate=False),
            #General search options:
            _Option(["-db", "db"],
                    "The database to BLAST against.",
                    equate=False),
            _Option(["-evalue", "evalue"], 
                    "Expectation value cutoff.",
                    equate=False),
            _Option(["-word_size","word_size"],
                    """Word size for wordfinder algorithm.

                    Integer. Minimum 2.""",
                    equate=False),
            #BLAST-2-Sequences options:
            # - see subclass
            #Formatting options:
            # - see baseclass
            #Query filtering options
            # TODO -soft_masking <Boolean>, is this a switch or an option?
            #_Switch(["-soft_masking", "soft_masking"],
            #        "Apply filtering locations as soft masks?"),
            _Switch(["-lcase_masking", "lcase_masking"],
                    "Use lower case filtering in query and subject sequence(s)?"),
            #Restrict search or results
            _Option(["-gilist", "gilist"],
                    """Restrict search of database to list of GI's.
 
                    Incompatible with: negative_gilist, seqidlist, remote, subject, subject_loc""",
                    filename=True,
                    equate=False),
            _Option(["-negative_gilist", "negative_gilist"],
                    """Restrict search of database to everything except the listed GIs.
 
                    Incompatible with: gilist, seqidlist, remote, subject, subject_loc""",
                    filename=True,
                    equate=False),
            _Option(["-seqidlist", "seqidlist"],
                    """Restrict search of database to list of SeqID's.
 
                    Incompatible with: gilist, negative_gilist, remote, subject, subject_loc""",
                    filename=True,
                    equate=False),
            _Option(["-entrez_query", "entrez_query"],
                    "Restrict search with the given Entrez query (requires remote).",
                    equate=False),
            _Option(["-max_target_seqs", "max_target_seqs"],
                    """Maximum number of aligned sequences to keep.

                    Integer argument (at least one).""",
                    equate=False),
            #Statistical options
            _Option(["-dbsize", "dbsize"],
                    "Effective length of the database (integer)",
                    equate=False),
            _Option(["-searchsp", "searchsp"],
                    "Effective length of the search space (integer)",
                    equate=False),
            #Extension options
            _Option(["-xdrop_ungap", "xdrop_ungap"],
                    "X-dropoff value (in bits) for ungapped extensions. Float.",
                    equate=False),
            _Option(["-xdrop_gap", "xdrop_gap"],
                    "X-dropoff value (in bits) for preliminary gapped extensions. Float.",
                    equate=False),
            _Option(["-xdrop_gap_final", "xdrop_gap_final"],
                    "X-dropoff value (in bits) for final gapped alignment. Float.",
                    equate=False),
            _Option(["-window_size", "window_size"],
                    "Multiple hits window size, use 0 to specify 1-hit algorithm. Integer.",
                    equate=False),
            # Search strategy options
            _Option(["-import_search_strategy", "import_search_strategy"],
                    """Search strategy to use.

                    Incompatible with: export_search_strategy""",
                    filename=True,
                    equate=False),
            _Option(["-export_search_strategy", "export_search_strategy"],
                    """File name to record the search strategy used.

                    Incompatible with: import_search_strategy""",
                    filename=True,
                    equate=False),
            #Miscellaneous options
            _Option(["-num_threads", "num_threads"],
                    """Number of threads to use in the BLAST search.

                    Integer of at least one. Default is one.
                    Incompatible with: remote""",
                    equate=False),
            _Switch(["-remote", "remote"],
                    """Execute search remotely?

                    Incompatible with: gilist, negative_gilist, subject_loc, num_threads, ..."""),
            ]
        try:
            #Insert extra parameters - at the start just in case there
            #are any arguments which must come last:
            self.parameters = extra_parameters + self.parameters
        except AttributeError:
            #Should we raise an error?  The subclass should have set this up!
            self.parameters = extra_parameters
        _NcbibaseblastCommandline.__init__(self, cmd, **kwargs)

    def _validate(self):
        incompatibles = {"remote":["gilist", "negative_gilist", "num_threads"],
                         "import_search_strategy" : ["export_search_strategy"],
                         "gilist":["negative_gilist"],
                         "seqidlist":["gilist", "negative_gilist", "remote"]}
        self._validate_incompatibilities(incompatibles)
        if self.entrez_query and not self.remote :
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
            #General search options:
            _Option(["-gapopen", "gapopen"],
                    "Cost to open a gap (integer).",
                    equate=False),
            _Option(["-gapextend", "gapextend"],
                    "Cost to extend a gap (integer).",
                    equate=False),
            #BLAST-2-Sequences options:
            _Option(["-subject", "subject"],
                    """Subject sequence(s) to search.

                    Incompatible with: db, gilist, negative_gilist.
                    See also subject_loc.""",
                    filename=True,
                    equate=False),
            _Option(["-subject_loc", "subject_loc"],
                    """Location on the subject sequence (Format: start-stop)

                    Incompatible with: db, gilist, seqidlist, negative_gilist,
                    db_soft_mask, db_hard_mask, remote.
                    
                    See also subject.""",
                    equate=False),
            #Restrict search or results:
            _Option(["-culling_limit", "culling_limit"],
                    """Hit culling limit (integer).

                    If the query range of a hit is enveloped by that of at
                    least this many higher-scoring hits, delete the hit.

                    Incompatible with: best_hit_overhang, best_hit_score_edge.
                    """,
                    equate=False),
            _Option(["-best_hit_overhang", "best_hit_overhang"],
                    """Best Hit algorithm overhang value (recommended value: 0.1)

                    Float between 0.0 and 0.5 inclusive.

                    Incompatible with: culling_limit.""",
                    equate=False),
            _Option(["-best_hit_score_edge", "best_hit_score_edge"],
                    """Best Hit algorithm score edge value (recommended value: 0.1)

                    Float between 0.0 and 0.5 inclusive.

                    Incompatible with: culling_limit.""",
                    equate=False),
            ]
        try:
            #Insert extra parameters - at the start just in case there
            #are any arguments which must come last:
            self.parameters = extra_parameters + self.parameters
        except AttributeError:
            #Should we raise an error?  The subclass should have set this up!
            self.parameters = extra_parameters
        _NcbiblastCommandline.__init__(self, cmd, **kwargs)


    def _validate(self):
        incompatibles = {"subject_loc":["db", "gilist", "negative_gilist", "seqidlist", "remote"],
                         "culling_limit":["best_hit_overhang","best_hit_score_edge"],
                         "subject":["db", "gilist", "negative_gilist", "seqidlist"]}
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
            #Restrict search or results:
            _Option(["-db_soft_mask", "db_soft_mask"],
                    """Filtering algorithm for soft masking (integer).

                    Filtering algorithm ID to apply to the BLAST database as soft masking.

                    Incompatible with: db_hard_mask, subject, subject_loc""",
                    equate=False),
            _Option(["-db_hard_mask", "db_hard_mask"],
                    """Filtering algorithm for hard masking (integer).

                    Filtering algorithm ID to apply to the BLAST database as hard masking.

                    Incompatible with: db_soft_mask, subject, subject_loc""",
                    equate=False),
            ]
        try:
            #Insert extra parameters - at the start just in case there
            #are any arguments which must come last:
            self.parameters = extra_parameters + self.parameters
        except AttributeError:
            #Should we raise an error?  The subclass should have set this up!
            self.parameters = extra_parameters
        _Ncbiblast2SeqCommandline.__init__(self, cmd, **kwargs)

    def _validate(self):
        incompatibles = {"db_soft_mask":["db_hard_mask", "subject", "subject_loc"],
                         "db_hard_mask":["db_soft_mask", "subject", "subject_loc"]}
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
    >>> print cline
    blastp -query rosemary.pro -db nr -evalue 0.001 -remote -ungapped

    You would typically run the command line with cline() or via the Python
    subprocess module, as described in the Biopython tutorial.
    """
    def __init__(self, cmd="blastp", **kwargs):
        self.parameters = [
            #General search options:
            _Option(["-task", "task"],
                    "Task to execute (string, blastp (default) or blastp-short).",
                    checker_function=lambda value : value in ["blastp",
                                                              "blastp-short"],
                    equate=False),
            _Option(["-matrix", "matrix"],
                    "Scoring matrix name (default BLOSUM62)."),
            _Option(["-threshold", "threshold"],
                    "Minimum word score such that the word is added to the "
                    "BLAST lookup table (float)",
                    equate=False),
            _Option(["-comp_based_stats", "comp_based_stats"],
                    """Use composition-based statistics (string, default 2, i.e. True).

                    0, F or f: no composition-based statistics
                    2, T or t, D or d : Composition-based score adjustment as in
                    Bioinformatics 21:902-911, 2005, conditioned on sequence properties

                    Note that tblastn also supports values of 1 and 3.""",
                    checker_function=lambda value : value in "0Ft2TtDd",
                    equate=False),
            #Query filtering options:
            _Option(["-seg", "seg"],
                    """Filter query sequence with SEG (string).

                    Format: "yes", "window locut hicut", or "no" to disable.
                    Default is "12 2.2 2.5""",
                    equate=False),
            #Extension options:
            _Switch(["-ungapped", "ungapped"],
                    "Perform ungapped alignment only?"),
            #Miscellaneous options:
            _Switch(["-use_sw_tback", "use_sw_tback"],
                    "Compute locally optimal Smith-Waterman alignments?"),
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
    >>> print cline
    blastn -out m_cold.xml -outfmt 5 -query m_cold.fasta -db nt -evalue 0.001 -strand plus

    You would typically run the command line with cline() or via the Python
    subprocess module, as described in the Biopython tutorial.
    """
    def __init__(self, cmd="blastn", **kwargs):
        self.parameters = [
            #Input query options:
            _Option(["-strand", "strand"],
                    """Query strand(s) to search against database/subject.

                    Values allowed are "both" (default), "minus", "plus".""",
                    checker_function=lambda value : value in ["both",
                                                              "minus",
                                                              "plus"],
                    equate=False),
            #General search options:
            _Option(["-task", "task"],
                    """Task to execute (string, default 'megablast')

                    Allowed values 'blastn', 'blastn-short', 'dc-megablast', 'megablast'
                    (the default), or 'vecscreen'.""",
                    checker_function=lambda value : value in ['blastn',
                                                              'blastn-short',
                                                              'dc-megablast',
                                                              'megablast',
                                                              'vecscreen'],
                    equate=False),
            _Option(["-penalty", "penalty"],
                    "Penalty for a nucleotide mismatch (integer, at most zero).",
                    equate=False),
            _Option(["-reward", "reward"],
                    "Reward for a nucleotide match (integer, at least zero).",
                    equate=False),
            #TODO - Does this need an argument or is it a switch?
            #_Option(["-use_index", "use_index"],
            #        "Use MegaBLAST database index (boolean).",
            #        equate=False),
            _Option(["-index_name", "index_name"],
                    "MegaBLAST database index name.",
                    equate=False),
            #Query filtering options:
            _Option(["-dust", "dust"],
                    """Filter query sequence with DUST (string).

                    Format: 'yes', 'level window linker', or 'no' to disable.
                    Default = '20 64 1'.
                    """,
                    equate=False),
            _Option(["-filtering_db", "filtering_db"],
                    "BLAST database containing filtering elements (i.e. repeats).",
                    equate=False),
            _Option(["-window_masker_taxid", "window_masker_taxid"],
                    "Enable WindowMasker filtering using a Taxonomic ID (integer).",
                    equate=False),
            _Option(["-window_masker_db", "window_masker_db"],
                    "Enable WindowMasker filtering using this repeats database (string).",
                    equate=False),
            #Restrict search or results:
            _Option(["-perc_identity", "perc_identity"],
                    "Percent identity (real, 0 to 100 inclusive).",
                    equate=False),
            #Discontiguous MegaBLAST options
            _Option(["-template_type", "template_type"],
                    """Discontiguous MegaBLAST template type (string).

                    Allowed values: 'coding', 'coding_and_optimal' or 'optimal'
                    Requires: template_length.""",
                    checker_function=lambda value : value in ['coding', 'coding_and_optimal','optimal'],
                    equate=False),
            _Option(["-template_length", "template_length"],
                    """Discontiguous MegaBLAST template length (integer).

                    Allowed values: 16, 18, 21
                    
                    Requires: template_type.""",
                    checker_function=lambda value : value in [16,18,21,'16','18','21'],
                    equate=False),
            #Extension options:
            _Switch(["-no_greedy", "no_greedy"],
                    "Use non-greedy dynamic programming extension"),
            _Option(["-min_raw_gapped_score", "min_raw_gapped_score"],
                    "Minimum raw gapped score to keep an alignment in the "
                    "preliminary gapped and traceback stages (integer).",
                    equate=False),
            _Switch(["-ungapped", "ungapped"],
                    "Perform ungapped alignment only?"),
            _Option(["-off_diagonal_range", "off_diagonal_range"],
                    """Number of off-diagonals to search for the 2nd hit (integer).
                    
                    Expects a positive integer, or 0 (default) to turn off.
                    
                    Added in BLAST 2.2.23+
                    """,
                    equate=False),
            ]
        _NcbiblastMain2SeqCommandline.__init__(self, cmd, **kwargs)

    def _validate(self):
        if (self.template_type and not self.template_length) \
        or (self.template_length and not self.template_type) :
            raise ValueError("Options template_type and template_type require each other.")
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
    >>> print cline
    blastx -query m_cold.fasta -db nr -evalue 0.001

    You would typically run the command line with cline() or via the Python
    subprocess module, as described in the Biopython tutorial.
    """
    def __init__(self, cmd="blastx", **kwargs):
        self.parameters = [
            #Input query options:
            _Option(["-strand", "strand"],
                    """Query strand(s) to search against database/subject.

                    Values allowed are "both" (default), "minus", "plus".""",
                    checker_function=lambda value : value in ["both", "minus", "plus"],
                    equate=False),
            #Input query options:
            _Option(["-query_gencode", "query_gencode"],
                    """Genetic code to use to translate query

                    Integer. Default is one.""",
                    equate=False),
            #General search options:
            _Option(["-frame_shift_penalty", "frame_shift_penalty"],
                    "Frame shift penalty (integer, at least 1, default ignored).",
                    equate=False),
            _Option(["-max_intron_length", "max_intron_length"],
                    """Maximum intron length (integer).

                    Length of the largest intron allowed in a translated nucleotide
                    sequence when linking multiple distinct alignments (a negative
                    value disables linking). Default zero.""",
                    equate=False),
            _Option(["-matrix", "matrix"],
                    "Scoring matrix name (default BLOSUM62).",
                    equate=False),
            _Option(["-threshold", "threshold"],
                    "Minimum word score such that the word is added to the "
                    "BLAST lookup table (float)",
                    equate=False),
            #Query filtering options:
            _Option(["-seg", "seg"],
                    """Filter query sequence with SEG (string).

                    Format: "yes", "window locut hicut", or "no" to disable.
                    Default is "12 2.2 2.5""",
                    equate=False),
            #Extension options:
            _Switch(["-ungapped", "ungapped"],
                    "Perform ungapped alignment only?"),
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
    >>> print cline
    tblastn -help

    You would typically run the command line with cline() or via the Python
    subprocess module, as described in the Biopython tutorial.
    """
    def __init__(self, cmd="tblastn", **kwargs):
        self.parameters = [
            #General search options:
            _Option(["-db_gencode", "db_gencode"],
                    """Genetic code to use to translate query

                    Integer. Default is one.""",
                    equate=False),
            _Option(["-frame_shift_penalty", "frame_shift_penalty"],
                    "Frame shift penalty (integer, at least 1, default ignored).",
                    equate=False),
            _Option(["-max_intron_length", "max_intron_length"],
                    """Maximum intron length (integer).

                    Length of the largest intron allowed in a translated nucleotide
                    sequence when linking multiple distinct alignments (a negative
                    value disables linking). Default zero.""",
                    equate=False),
            _Option(["-matrix", "matrix"],
                    "Scoring matrix name (default BLOSUM62).",
                    equate=False),
            _Option(["-threshold", "threshold"],
                    "Minimum word score such that the word is added to the BLAST lookup table (float)",
                    equate=False),
            _Option(["-comp_based_stats", "comp_based_stats"],
                    """Use composition-based statistics (string, default 2, i.e. True).

                    0, F or f: no composition-based statistics
                    1: Composition-based statistics as in NAR 29:2994-3005, 2001
                    2, T or t, D or d : Composition-based score adjustment as in
                       Bioinformatics 21:902-911, 2005, conditioned on sequence properties
                    3: Composition-based score adjustment as in Bioinformatics 21:902-911,
                       2005, unconditionally

                    Note that only tblastn supports values of 1 and 3.""",
                    checker_function=lambda value : value in "0Ft12TtDd3",
                    equate=False),
            #Query filtering options:
            _Option(["-seg", "seg"],
                    """Filter query sequence with SEG (string).

                    Format: "yes", "window locut hicut", or "no" to disable.
                    Default is "12 2.2 2.5""",
                    equate=False),
            #Extension options:
            _Switch(["-ungapped", "ungapped"],
                    "Perform ungapped alignment only?"),
            #Miscellaneous options:
            _Switch(["-use_sw_tback", "use_sw_tback"],
                    "Compute locally optimal Smith-Waterman alignments?"),
            #PSI-TBLASTN options:
            _Option(["-in_pssm", "in_pssm"],
                    """PSI-BLAST checkpoint file

                    Incompatible with: remote, query""",
                    filename=True,
                    equate=False),
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
    >>> print cline
    tblastx -help

    You would typically run the command line with cline() or via the Python
    subprocess module, as described in the Biopython tutorial.
    """
    def __init__(self, cmd="tblastx", **kwargs):
        self.parameters = [
            #Input query options:
            _Option(["-strand", "strand"],
                    """Query strand(s) to search against database/subject.

                    Values allowed are "both" (default), "minus", "plus".""",
                    checker_function=lambda value : value in ["both", "minus", "plus"],
                    equate=False),
            #Input query options:
            _Option(["-query_gencode", "query_gencode"],
                    """Genetic code to use to translate query

                    Integer. Default is one.""",
                    equate=False),
            #General search options:
            _Option(["-db_gencode", "db_gencode"],
                    """Genetic code to use to translate query

                    Integer. Default is one.""",
                    equate=False),
            _Option(["-max_intron_length", "max_intron_length"],
                    """Maximum intron length (integer).

                    Length of the largest intron allowed in a translated nucleotide
                    sequence when linking multiple distinct alignments (a negative
                    value disables linking). Default zero.""",
                    equate=False),
            _Option(["-matrix", "matrix"],
                    "Scoring matrix name (default BLOSUM62).",
                    equate=False),
            _Option(["-threshold", "threshold"],
                    "Minimum word score such that the word is added to the "
                    "BLAST lookup table (float)",
                    equate=False),
            #Query filtering options:
            _Option(["-seg", "seg"],
                    """Filter query sequence with SEG (string).

                    Format: "yes", "window locut hicut", or "no" to disable.
                    Default is "12 2.2 2.5""",
                    equate=False),
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
    >>> print cline
    psiblast -help

    You would typically run the command line with cline() or via the Python
    subprocess module, as described in the Biopython tutorial.
    """
    def __init__(self, cmd="psiblast", **kwargs):
        self.parameters = [
            #General search options:
            _Option(["-matrix", "matrix"],
                    "Scoring matrix name (default BLOSUM62).",
                    equate=False),
            _Option(["-threshold", "threshold"],
                    "Minimum word score such that the word is added to the "
                    "BLAST lookup table (float)",
                    equate=False),
            _Option(["-comp_based_stats", "comp_based_stats"],
                    """Use composition-based statistics (string, default 2, i.e. True).

                    0, F or f: no composition-based statistics
                    2, T or t, D or d : Composition-based score adjustment
                    as in Bioinformatics 21:902-911, 2005, conditioned on
                    sequence properties

                    Note that tblastn also supports values of 1 and 3.""",
                    checker_function=lambda value : value in "0Ft2TtDd",
                    equate=False),
            #Query filtering options:
            _Option(["-seg", "seg"],
                    """Filter query sequence with SEG (string).

                    Format: "yes", "window locut hicut", or "no" to disable.
                    Default is "12 2.2 2.5""",
                    equate=False),
            #Extension options:
            _Option(["-gap_trigger", "gap_trigger"],
                    "Number of bits to trigger gapping (float, default 22)",
                    equate=False),
            #Miscellaneous options:
            _Switch(["-use_sw_tback", "use_sw_tback"],
                    "Compute locally optimal Smith-Waterman alignments?"),
            #PSI-BLAST options:
            _Option(["-num_iterations", "num_iterations"],
                    """Number of iterations to perform, integer

                    Integer of at least one. Default is one.
                    Incompatible with: remote""",
                    equate=False),
            _Option(["-out_pssm", "out_pssm"],
                    "File name to store checkpoint file",
                    filename=True,
                    equate=False),
            _Option(["-out_ascii_pssm", "out_ascii_pssm"],
                    "File name to store ASCII version of PSSM",
                    filename=True,
                    equate=False),
            _Option(["-in_msa", "in_msa"],
                    """File name of multiple sequence alignment to restart
                    PSI-BLAST

                    Incompatible with: in_pssm, query""",
                    filename=True,
                    equate=False),
            _Option(["-msa_master_idx", "msa_master_idx"],
                    """Index of sequence to use as master in MSA.

                    Index (1-based) of sequence to use as the master in the
                    multiple sequence alignment. If not specified, the first
                    sequence is used.""",
                    equate=False),
            _Option(["-in_pssm", "in_pssm"],
                    """PSI-BLAST checkpoint file

                    Incompatible with: in_msa, query, phi_pattern""",
                    filename=True,
                    equate=False),
            #PSSM engine options:
            _Option(["-pseudocount", "pseudocount"],
                    """Pseudo-count value used when constructing PSSM

                    Integer. Default is zero.""",
                    equate=False),
            _Option(["-inclusion_ethresh", "inclusion_ethresh"],
                    """E-value inclusion threshold for pairwise alignments

                    Float. Default is 0.002.""",
                    equate=False),
            #PHI-BLAST options:
            _Option(["-phi_pattern", "phi_pattern"],
                    """File name containing pattern to search

                    Incompatible with: in_pssm""",
                    filename=True,
                    equate=False),
            ]
        _Ncbiblast2SeqCommandline.__init__(self, cmd, **kwargs)

    def _validate(self):
        incompatibles = {"num_iterations":["remote"],
                         "in_msa":["in_pssm", "query"],
                         "in_pssm":["in_msa","query","phi_pattern"]}
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
    >>> print cline
    rpsblast -help

    You would typically run the command line with cline() or via the Python
    subprocess module, as described in the Biopython tutorial.
    """
    def __init__(self, cmd="rpsblast", **kwargs):
        self.parameters = [
            #Query filtering options:
            _Option(["-seg", "seg"],
                    """Filter query sequence with SEG (string).

                    Format: "yes", "window locut hicut", or "no" to disable.
                    Default is "12 2.2 2.5""",
                    equate=False),
            #Restrict search or results:
            _Option(["-culling_limit", "culling_limit"],
                    """Hit culling limit (integer).

                    If the query range of a hit is enveloped by that of at
                    least this many higher-scoring hits, delete the hit.

                    Incompatible with: best_hit_overhang, best_hit_score_edge.
                    """,
                    equate=False),
            _Option(["-best_hit_overhang", "best_hit_overhang"],
                    """Best Hit algorithm overhang value (recommended value: 0.1)

                    Float between 0.0 and 0.5 inclusive.

                    Incompatible with: culling_limit.""",
                    equate=False),
            _Option(["-best_hit_score_edge", "best_hit_score_edge"],
                    """Best Hit algorithm score edge value (recommended value: 0.1)

                    Float between 0.0 and 0.5 inclusive.

                    Incompatible with: culling_limit.""",
                    equate=False),
            ]
        _NcbiblastCommandline.__init__(self, cmd, **kwargs)

    def _validate(self):
        incompatibles = {"culling_limit":["best_hit_overhang","best_hit_score_edge"]}
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
    >>> print cline
    rpstblastn -help

    You would typically run the command line with cline() or via the Python
    subprocess module, as described in the Biopython tutorial.
    """
    def __init__(self, cmd="rpstblastn", **kwargs):
        self.parameters = [
            #Input query options:
            _Option(["-strand", "strand"],
                    """Query strand(s) to search against database/subject.

                    Values allowed are "both" (default), "minus", "plus".""",
                    checker_function=lambda value : value in ["both",
                                                              "minus",
                                                              "plus"],
                    equate=False),
            #Input query options:
            _Option(["-query_gencode", "query_gencode"],
                    """Genetic code to use to translate query

                    Integer. Default is one.""",
                    equate=False),
            #Query filtering options:
            _Option(["-seg", "seg"],
                    """Filter query sequence with SEG (string).

                    Format: "yes", "window locut hicut", or "no" to disable.
                    Default is "12 2.2 2.5""",
                    equate=False),
            #Extension options:
            _Switch(["-ungapped", "ungapped"],
                    "Perform ungapped alignment only?"),
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
    >>> print cline
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
        self.parameters = [
            # Input options
            _Option(["-rid", "rid"],
                    "BLAST Request ID (RID), not compatiable with archive arg",
                    equate=False),
            _Option(["-archive", "archive"],
                    "Archive file of results, not compatiable with rid arg.",
                    filename=True,
                    equate=False),
            # Restrict search or results
            _Option(["-max_target_seqs", "max_target_seqs"],
                    "Maximum number of aligned sequences to keep",
                    checker_function=lambda value: value >= 1,
                    equate=False),
            ]
        _NcbibaseblastCommandline.__init__(self, cmd, **kwargs)

    def _validate(self):
        incompatibles = {"rid":["archive"]}
        self._validate_incompatibilities(incompatibles)
        _NcbibaseblastCommandline._validate(self)
        

def _test():
    """Run the Bio.Blast.Applications module's doctests."""
    import doctest
    doctest.testmod(verbose=1)

if __name__ == "__main__":
    #Run the doctests
    _test()
