# Copyright 2001 Brad Chapman.
# Revisions copyright 2009 by Peter Cock.
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

"""
from Bio.Application import _Option, AbstractCommandline, _Switch

class FastacmdCommandline(AbstractCommandline):
    """Create a commandline for the fasta program from NCBI (OBSOLETE).

    """
    def __init__(self, cmd="fastacmd", **kwargs):
        self.parameters = \
          [
           _Option(["-d", "database"], ["input"], None, 1,
                   "The database to retrieve from."),
           _Option(["-s", "search_string"], ["input"], None, 1,
                   "The id to search for.")
          ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


class _BlastCommandLine(AbstractCommandline):
    """Base Commandline object for (classic) NCBI BLAST wrappers (PRIVATE).

    This is provided for subclassing, it deals with shared options
    common to all the BLAST tools (blastall, rpsblast, blastpgp).
    """
    def __init__(self, cmd=None, **kwargs):
        assert cmd is not None
        extra_parameters = [\
           _Switch(["--help", "help"], ["input"],
                    "Print USAGE, DESCRIPTION and ARGUMENTS description;  ignore other arguments."),
           _Option(["-d", "database"], ["input"], None, 1,
                   "The database to BLAST against.", False),
           _Option(["-i", "infile"], ["input", "file"], None, 1,
                   "The sequence to search with.", False),
           _Option(["-e", "expectation"], ["input"], None, 0, 
                   "Expectation value cutoff.", False),
           _Option(["-m", "align_view"], ["input"], None, 0, 
                   "Alignment view.  Integer 0-11.  Use 7 for XML output.",
                   False),
           _Option(["-o", "align_outfile", "outfile"], ["output", "file"], None, 0,
                   "Output file for alignment.", False),
           _Option(["-y", "xdrop_extension"], ["input"], None, 0, 
                   "Dropoff for blast extensions.", False),
           _Option(["-F", "filter"], ["input"], None, 0,
                   "Filter query sequence with SEG?  T/F", False),
           _Option(["-X", "xdrop"], ["input"], None, 0, 
                   "Dropoff value (bits) for gapped alignments."),
           _Option(["-I", "show_gi"], ["input"], None, 0, 
                   "Show GI's in deflines?  T/F", False),
           _Option(["-J", "believe_query"], ["input"], None, 0, 
                   "Believe the query defline?  T/F", False),
           _Option(["-Z", "xdrop_final"], ["input"], None, 0, 
                   "X dropoff for final gapped alignment.", False),
           _Option(["-z", "db_length"], ["input"], None, 0, 
                   "Effective database length.", False),
           _Option(["-O", "seqalign_file"], ["output", "file"], None, 0,
                   "seqalign file to output.", False),
           _Option(["-v", "descriptions"], ["input"], None, 0, 
                   "Number of one-line descriptions.", False),
           _Option(["-b", "alignments"], ["input"], None, 0, 
                   "Number of alignments.", False),
           _Option(["-Y", "search_length"], ["input"], None, 0, 
                   "Effective length of search space (use zero for the " + \
                   "real size).", False),
           _Option(["-T", "html"], ["input"], None, 0, 
                   "Produce HTML output?  T/F", False),
           _Option(["-U", "case_filter"], ["input"], None, 0,
                   "Use lower case filtering of FASTA sequence? T/F", False),

           _Option(["-a", "nprocessors"], ["input"], None, 0,
                   "Number of processors to use.", False),
           _Option(["-g", "gapped"], ["input"], None, 0, 
                   "Whether to do a gapped alignment.  T/F", False),
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
        extra_parameters = [\
           _Option(["-G", "gap_open"], ["input"], None, 0, 
                   "Gap open penalty", False),
           _Option(["-E", "gap_extend"], ["input"], None, 0, 
                    "Gap extension penalty", False),
           _Option(["-A", "window_size"], ["input"], None, 0,
                    "Multiple hits window size", False),
           _Option(["-f", "hit_extend"], ["input"], None, 0, 
                   "Threshold for extending hits.", False),
           _Option(["-K", "keep_hits"], ["input"], None, 0,
                   " Number of best hits from a region to keep.", False),
           _Option(["-W", "wordsize"], ["input"], None, 0, 
                   "Word size", False),
           _Option(["-P", "passes"], ["input"], None, 0,
                   "Hits/passes.  Integer 0-2. 0 for multiple hit, "
                   "1 for single hit (does not apply to blastn)", False),
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

    You would typically run the command line with the Python subprocess module,
    as described in the Biopython tutorial.
    """
    #TODO - This could use more checking for valid parameters to the program.
    def __init__(self, cmd="blastall",**kwargs):
        self.parameters = [ \
            #Sorted in the same order as the output from blastall --help
            #which should make it easier to keep them up to date in future.
            #Note that some arguments are defined the the base clases (above).
           _Option(["-p", "program"], ["input"], None, 1, 
                   "The blast program to use (e.g. blastp, blastn).", False),
           _Option(["-q", "nuc_mismatch"], ["input"], None, 0, 
                   "Penalty for a nucleotide mismatch (blastn only).", False),
           _Option(["-r", "nuc_match"], ["input"], None, 0, 
                   "Reward for a nucleotide match (blastn only).", False),
           _Option(["-Q", "query_genetic_code"], ["input"], None, 0,
                   "Query Genetic code to use.", False),
           _Option(["-D", "db_genetic_code"], ["input"], None, 0,
                   "DB Genetic code (for tblast[nx] only).", False),
           _Option(["-M", "matrix"], ["input"], None, 0, 
                   "Matrix to use", False),
           _Option(["-S", "strands"], ["input"], None, 0, 
                   "Query strands to search against database (for blast[nx], " + \
                   "and tblastx). 3 is both, 1 is top, 2 is bottom.", False),
           _Option(["-l", "restrict_gi"], ["input"], None, 0,
                   "Restrict search of database to list of GI's.", False),
           _Option(["-R"], ["input", "file"], None, 0,
                   "PSI-TBLASTN checkpoint input file.", False),
           _Option(["-n", "megablast"], ["input"], None, 0,
                   "MegaBlast search T/F.", False),
           #The old name "region_length" is for consistency with our
           #old blastall function wrapper:
           _Option(["-L", "region_length", "range_restriction"], ["input"],
                   None, 0, 
                   """Location on query sequence (string format start,end).

                   In older versions of BLAST, -L set the length of region
                   used to judge hits (see -K parameter).""", False),
           _Option(["-w"], ["input"], None, 0,
                   "Frame shift penalty (OOF algorithm for blastx).", False),
           _Option(["-t"], ["input"], None, 0,
                   "Length of the largest intron allowed in a translated " + \
                   "nucleotide sequence when linking multiple distinct " + \
                   "alignments. (0 invokes default behavior; a negative value " + \
                   "disables linking.)", False),
           _Option(["-B"], ["input"], None, 0,
                   "Number of concatenated queries, for blastn and tblastn.",
                   False),
           _Option(["-V", "oldengine"], ["input"], None, 0,
                   "Force use of the legacy BLAST engine.", False),
           _Option(["-C"], ["input"], None, 0,
                   """Use composition-based statistics for tblastn:
                   D or d: default (equivalent to F)
                   0 or F or f: no composition-based statistics
                   1 or T or t: Composition-based statistics as in NAR 29:2994-3005, 2001
                   2: Composition-based score adjustment as in Bioinformatics
                       21:902-911, 2005, conditioned on sequence properties
                   3: Composition-based score adjustment as in Bioinformatics
                       21:902-911, 2005, unconditionally
                   For programs other than tblastn, must either be absent or be
                   D, F or 0.""", False),
           _Option(["-s"], ["input"], None, 0,
                   "Compute locally optimal Smith-Waterman alignments (This " + \
                   "option is only available for gapped tblastn.) T/F", False),
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

    You would typically run the command line with the Python subprocess module,
    as described in the Biopython tutorial.
    """
    def __init__(self, cmd="blastpgp",**kwargs):
        self.parameters = [ \
           _Option(["-C", "checkpoint_outfile"], ["output", "file"], None, 0,
                   "Output file for PSI-BLAST checkpointing.", False),
           _Option(["-R", "restart_infile"], ["input", "file"], None, 0,
                   "Input file for PSI-BLAST restart.", False),
           _Option(["-k", "hit_infile"], ["input", "file"], None, 0,
                   "Hit file for PHI-BLAST.", False),
           _Option(["-Q", "matrix_outfile"], ["output", "file"], None, 0,
                   "Output file for PSI-BLAST matrix in ASCII.", False),
           _Option(["-B", "align_infile"], ["input", "file"], None, 0, 
                   "Input alignment file for PSI-BLAST restart.", False),
           _Option(["-S", "required_start"], ["input"], None, 0, 
                   "Start of required region in query.", False),
           _Option(["-H", "required_end"], ["input"], None, 0,
                   "End of required region in query.", False),
           _Option(["-j", "npasses"], ["input"], None, 0,
                    "Number of passes", False),
           _Option(["-N", "nbits_gapping"], ["input"], None, 0, 
                   "Number of bits to trigger gapping.", False),
           _Option(["-c", "pseudocounts"], ["input"], None, 0,
                   "Pseudocounts constants for multiple passes.", False),
           _Option(["-h", "model_threshold"], ["input"], None, 0, 
                   "E-value threshold to include in multipass model.", False),
           #Does the old name "region_length" for -L make sense?
           _Option(["-L", "region_length"], ["input"], None, 0, 
                   "Cost to decline alignment (disabled when zero).", False),
           _Option(["-M", "matrix"], ["input"], None, 0, 
                   "Matrix (string, default BLOSUM62).", False),
           _Option(["-p", "program"], ["input"], None, 1, 
                   "The blast program to use (e.g blastpgp, patseedp or seedp).", False),
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

    You would typically run the command line with the Python subprocess module,
    as described in the Biopython tutorial.
    """
    def __init__(self, cmd="rpsblast",**kwargs):
        self.parameters = [ \
           #Note -N is also in blastpgp, but not blastall
           _Option(["-N", "nbits_gapping"], ["input"], None, 0, 
                   "Number of bits to trigger gapping.", False),
           #Note blastall and blastpgp wrappers have -P with name "passes".
           #If this is the same thing, we should be consistent!
           _Option(["-P", "multihit"], ["input"], None, 0,
                   "0 for multiple hit, 1 for single hit", False),
           _Option(["-l", "logfile"], ["output", "file"], None, 0, 
                   "Logfile name.", False),
           _Option(["-p", "protein"], ["input"], None, 0, 
                   "Query sequence is protein. T/F", False),
           _Option(["-L", "range_restriction"], ["input"], None, 0, 
                   "Location on query sequence (string format start,end).",
                   False),
        ] 
        _BlastCommandLine.__init__(self, cmd, **kwargs)

   
class _NcbiblastCommandline(AbstractCommandline):
    """Base Commandline object for (classic) NCBI BLAST wrappers (PRIVATE).

    This is provided for subclassing, it deals with shared options
    common to all the BLAST tools (blastn, rpsblast, rpsblast, etc).
    """
    def __init__(self, cmd=None, **kwargs):
        assert cmd is not None
        extra_parameters = [ \
            #Core:
            _Switch(["-h", "h"], ["input"],
                    "Print USAGE and DESCRIPTION;  ignore other arguments."),
            _Switch(["-help", "help"], ["input"],
                    "Print USAGE, DESCRIPTION and ARGUMENTS description;  ignore other arguments."),
            _Switch(["-version", "version"], ["input"],
                    "Print version number;  ignore other arguments."),
            #Input query options:
            _Option(["-query", "query"], ["input", "file"], None, 0,
                    "The sequence to search with.", False), #Should this be required?
            _Option(["-query_loc", "query_loc"], ["input"], None, 0,
                    "Location on the query sequence (Format: start-stop)", False),
            #General search options:
            _Option(["-db", "db"], ["input"], None, 0,
                    "The database to BLAST against.", False), #Should this be required?
            _Option(["-out", "out"], ["output", "file"], None, 0,
                    "Output file for alignment.", False),
            _Option(["-evalue", "evalue"], ["input"], None, 0, 
                    "Expectation value cutoff.", False),
            _Option(["-word_size","word_size"], ["input"], None, 0,
                    """Word size for wordfinder algorithm.

                    Integer. Minimum 2.""", False),
            #BLAST-2-Sequences options:
            # - see subclass
            #Formatting options:
            _Option(["-outfmt", "outfmt"], ["input"], None, 0, 
                    "Alignment view.  Integer 0-10.  Use 5 for XML output (differs from classic BLAST which used 7 for XML).",
                    False), #Did not include old aliases as meaning has changed!
            _Switch(["-show_gis","show_gis"], ["input"],
                    "Show NCBI GIs in deflines?"),
            _Option(["-num_descriptions","num_descriptions"], ["input"], None, 0,
                    """Number of database sequences to show one-line descriptions for.

                    Integer argument (at least zero). Default is 500.
                    See also num_alignments.""", False),
            _Option(["-num_alignments","num_alignments"], ["input"], None, 0,
                    """Number of database sequences to show num_alignments for.

                    Integer argument (at least zero). Default is 200.
                    See also num_alignments.""", False),
            _Switch(["-html", "html"], ["input"],
                    "Produce HTML output? See also the outfmt option."),
            #Query filtering options
            # TODO -soft_masking <Boolean>, is this a switch or an option?
            #_Switch(["-soft_masking", "soft_masking"], ["input"],
            #        "Apply filtering locations as soft masks?"),
            _Switch(["-lcase_masking", "lcase_masking"], ["input"],
                    "Use lower case filtering in query and subject sequence(s)?"),
            #Restrict search or results
            _Option(["-gilist", "gilist"], ["input", "file"], None, 0,
                    """Restrict search of database to list of GI's.
 
                    Incompatible with: negative_gilist, remote, subject, subject_loc""",
                    False),
            _Option(["-negative_gilist", "negative_gilist"], ["input", "file"], None, 0,
                    """Restrict search of database to everything except the listed GIs.
 
                    Incompatible with: gilist, remote, subject, subject_loc""",
                    False),
            _Option(["-entrez_query", "entrez_query"], ["input"], None, 0,
                    "Restrict search with the given Entrez query (requires remote).", False),
            _Option(["-max_target_seqs", "max_target_seqs"], ["input"], None, 0,
                    """Maximum number of aligned sequences to keep.

                    Integer argument (at least one).""", False),
            #Statistical options
            _Option(["-dbsize", "dbsize"], ["input"], None, 0,
                    "Effective length of the database (integer)", False),
            _Option(["-searchsp", "searchsp"], ["input"], None, 0,
                    "Effective length of the search space (integer)", False),
            #Extension options
            _Option(["-xdrop_ungap", "xdrop_ungap"], ["input"], None, 0,
                    "X-dropoff value (in bits) for ungapped extensions. Float.",
                    False),
            _Option(["-xdrop_gap", "xdrop_gap"], ["input"], None, 0,
                    "X-dropoff value (in bits) for preliminary gapped extensions. Float.",
                    False),
            _Option(["-xdrop_gap_final", "xdrop_gap_final"], ["input"], None, 0,
                    "X-dropoff value (in bits) for final gapped alignment. Float.",
                    False),
            _Option(["-window_size", "window_size"], ["input"], None, 0,
                    "Multiple hits window size, use 0 to specify 1-hit algorithm. Integer.",
                    False),
            # Search strategy options
            _Option(["-import_search_strategy", "import_search_strategy"],
                    ["input", "file"], None, 0,
                    """Search strategy to use.

                    Incompatible with: export_search_strategy""", False),
            _Option(["-export_search_strategy", "export_search_strategy"],
                    ["output", "file"], None, 0,
                    """File name to record the search strategy used.

                    Incompatible with: import_search_strategy""", False),
            #Miscellaneous options
            _Switch(["-parse_deflines", "parse_deflines"], ["input"],
                    "Should the query and subject defline(s) be parsed?"),
            _Option(["-num_threads", "num_threads"], ["input"], None, 0,
                    """Number of threads to use in the BLAST search.

                    Integer of at least one. Default is one.
                    Incompatible with: remote""", False),
            _Switch(["-remote", "remote"], ["input"],
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
        AbstractCommandline.__init__(self, cmd, **kwargs)

    def _validate(self):
        incompatibles = {"remote":["gilist", "negative_gilist", "num_threads"],
                         "import_search_strategy" : ["export_search_strategy"],
                         "gilist":["negative_gilist"]}
        for a in incompatibles:
            if self._get_parameter(a):
                for b in incompatibles[a]:
                    if self._get_parameter(b):
                        raise ValueError("Options %s and %s are incompatible." \
                                         % (a,b))
        if self.entrez_query and not self.remote :
            raise ValueError("Option entrez_query requires remote option.")
        AbstractCommandline._validate(self)

class _Ncbiblast2SeqCommandline(_NcbiblastCommandline):
    """Base Commandline object for (classic) NCBI BLAST wrappers (PRIVATE).

    This is provided for subclassing, it deals with shared options
    common to all the BLAST tools supporting two-sequence BLAST
    (blastn, psiblast, etc) but not rpsblast or rpstblastn.
    """
    def __init__(self, cmd=None, **kwargs):
        assert cmd is not None
        extra_parameters = [ \
            #General search options:
            _Option(["-gapopen", "gapopen"], ["input"], None, 0,
                    "Cost to open a gap (integer).", False),
            _Option(["-gapextend", "gapextend"], ["input"], None, 0,
                    "Cost to extend a gap (integer).", False),
            #BLAST-2-Sequences options:
            _Option(["-subject", "subject"], ["input", "file"], None, 0,
                    """Subject sequence(s) to search.

                    Incompatible with: db, gilist, negative_gilist.
                    See also subject_loc.""", False),
            _Option(["-subject_loc", "subject_loc"], ["input"], None, 0,
                    """Location on the subject sequence (Format: start-stop)

                    Incompatible with: db, gilist, negative_gilist, remote.
                    See also subject.""", False),
            #Restrict search or results:
            _Option(["-culling_limit", "culling_limit"], ["input"], None, 0,
                    """Hit culling limit (integer).

                    If the query range of a hit is enveloped by that of at least this many
                    higher-scoring hits, delete the hit.

                    Incompatible with: best_hit_overhang, best_hit_score_edge.""", False),
            _Option(["-best_hit_overhang", "best_hit_overhang"], ["input"], None, 0,
                    """Best Hit algorithm overhang value (recommended value: 0.1)

                    Float between 0.0 and 0.5 inclusive.

                    Incompatible with: culling_limit.""", False),
            _Option(["-best_hit_score_edge", "best_hit_score_edge"], ["input"], None, 0,
                    """Best Hit algorithm score edge value (recommended value: 0.1)

                    Float between 0.0 and 0.5 inclusive.

                    Incompatible with: culling_limit.""", False),            ]
        try:
            #Insert extra parameters - at the start just in case there
            #are any arguments which must come last:
            self.parameters = extra_parameters + self.parameters
        except AttributeError:
            #Should we raise an error?  The subclass should have set this up!
            self.parameters = extra_parameters
        _NcbiblastCommandline.__init__(self, cmd, **kwargs)


    def _validate(self):
        incompatibles = {"subject_loc":["db, gilist, negative_gilist, remote"],
                         "culling_limit":["best_hit_overhang","best_hit_score_edge"],
                         "subject":["db", "gilist", "negative_gilist"]}
        for a in incompatibles:
            if self._get_parameter(a):
                for b in incompatibles[a]:
                    if self._get_parameter(b):
                        raise ValueError("Options %s and %s are incompatible." \
                                         % (a,b))
        _NcbiblastCommandline._validate(self)

class NcbiblastpCommandline(_Ncbiblast2SeqCommandline):
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

    You would typically run the command line with the Python subprocess module,
    as described in the Biopython tutorial.
    """
    def __init__(self, cmd="blastp", **kwargs):
        self.parameters = [ \
            #General search options:
            _Option(["-task", "task"], ["input"],
                    lambda value : value in ["blastp", "blastp-short"], 0,
                    "Task to execute (string, blastp (default) or blastp-short).", False),
            _Option(["-matrix", "matrix"], ["input"], None, 0,
                    "Scoring matrix name (default BLOSUM62).", False),
            _Option(["-threshold", "threshold"], ["input"], None, 0,
                    "Minimum word score such that the word is added to the BLAST lookup table (float)", False),
            _Option(["-comp_based_stats", "comp_based_stats"], ["input"],
                    lambda value : value in "0Ft2TtDd", 0,
                    """Use composition-based statistics (string, default 2, i.e. True).

                    0, F or f: no composition-based statistics
                    2, T or t, D or d : Composition-based score adjustment as in
                    Bioinformatics 21:902-911, 2005, conditioned on sequence properties

                    Note that tblastn also supports values of 1 and 3.""", False),
            #Query filtering options:
            _Option(["-seg", "seg"], ["input"], None, 0,
                    """Filter query sequence with SEG (string).

                    Format: "yes", "window locut hicut", or "no" to disable.
                    Default is "12 2.2 2.5""", False),
            #Restrict search or results:
            _Option(["-db_soft_mask", "db_soft_mask"], ["input"], None, 0,
                    """Filtering algorithm for soft masking (integer).

                    Filtering algorithm ID to apply to the BLAST database as soft masking.

                    Incompatible with: subject, subject_loc""", False),
            #Extension options:
            _Switch(["-ungapped", "ungapped"], ["input"],
                    "Perform ungapped alignment only?"),
            #Miscellaneous options:
            _Switch(["-use_sw_tback", "use_sw_tback"], ["input"],
                    "Compute locally optimal Smith-Waterman alignments?"),
            ]
        _Ncbiblast2SeqCommandline.__init__(self, cmd, **kwargs)

    def _validate(self):
        incompatibles = {"db_soft_mask":["subject", "subject_loc"]}
        for a in incompatibles:
            if self._get_parameter(a):
                for b in incompatibles[a]:
                    if self._get_parameter(b):
                        raise ValueError("Options %s and %s are incompatible." \
                                         % (a,b))
        _Ncbiblast2SeqCommandline._validate(self)


class NcbiblastnCommandline(_Ncbiblast2SeqCommandline):
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
    NcbiblastnCommandline(cmd='blastn', query='m_cold.fasta', db='nt', out='m_cold.xml', evalue=0.001, outfmt=5, strand='plus')
    >>> print cline
    blastn -query m_cold.fasta -db nt -out m_cold.xml -evalue 0.001 -outfmt 5 -strand plus

    You would typically run the command line with the Python subprocess module,
    as described in the Biopython tutorial.
    """
    def __init__(self, cmd="blastn", **kwargs):
        self.parameters = [ \
            #Input query options:
            _Option(["-strand", "strand"], ["input"],
                    lambda value : value in ["both", "minus", "plus"],0,
                    """Query strand(s) to search against database/subject.

                    Values allowed are "both" (default), "minus", "plus".""", False),
            #General search options:
            _Option(["-task", "task"], ["input"],
                    lambda value : value in ['blastn', 'blastn-short', 'dc-megablast',
                                             'megablast', 'vecscreen'], 0,
                    """Task to execute (string, default 'megablast')

                    Allowed values 'blastn', 'blastn-short', 'dc-megablast', 'megablast'
                    (the default), or 'vecscreen'.""", False),
            _Option(["-penalty", "penalty"], ["input"], None, 0,
                    "Penalty for a nucleotide mismatch (integer, at most zero).", False),
            _Option(["-reward", "reward"], ["input"], None, 0,
                    "Reward for a nucleotide match (integer, at least zero).", False),
            #TODO - Does this need an argument or is it a switch?
            #_Option(["-use_index", "use_index"], ["input"], None, 0,
            #        "Use MegaBLAST database index (boolean).", False),
            _Option(["-index_name", "index_name"], ["input"], None, 0,
                    "MegaBLAST database index name.", False),
            #Query filtering options:
            _Option(["-dust", "dust"], ["input"], None, 0,
                    """Filter query sequence with DUST (string).

                    Format: 'yes', 'level window linker', or 'no' to disable.
                    Default = '20 64 1'.
                    """, False),
            _Option(["-filtering_db", "filtering_db"], ["input"], None, 0,
                    "BLAST database containing filtering elements (i.e. repeats).", False),
            _Option(["-window_masker_taxid", "window_masker_taxid"], ["input"], None, 0,
                    "Enable WindowMasker filtering using a Taxonomic ID (integer).", False),
            _Option(["-window_masker_db", "window_masker_db"], ["input"], None, 0,
                    "Enable WindowMasker filtering using this repeats database (string).", False),
            #Restrict search or results:
            _Option(["-db_soft_mask", "db_soft_mask"], ["input"], None, 0,
                    """Filtering algorithm for soft masking (integer).

                    Filtering algorithm ID to apply to the BLAST database as soft masking.

                    Incompatible with: subject, subject_loc""", False),
            _Option(["-perc_identity", "perc_identity"], ["input"], None, 0,
                    "Percent identity (real, 0 to 100 inclusive).", False),
            #Discontiguous MegaBLAST options
            _Option(["-template_type", "template_type"], ["input"],
                    lambda value : value in ['coding', 'coding_and_optimal','optimal'], 0,
                    """Discontiguous MegaBLAST template type (string).

                    Allowed values: 'coding', 'coding_and_optimal' or 'optimal'
                    Requires: template_length.""", False),
            _Option(["-template_length", "template_length"], ["input"],
                    lambda value : value in [16,18,21,'16','18','21'], 0,
                    """Discontiguous MegaBLAST template length (integer).

                    Allowed values: 16, 18, 21
                    
                    Requires: template_type.""", False),
            #Extension options:
            _Switch(["-no_greedy", "no_greedy"], ["input"],
                    "Use non-greedy dynamic programming extension"),
            _Option(["-min_raw_gapped_score", "min_raw_gapped_score"], ["input"], None, 0,
                    "Minimum raw gapped score to keep an alignment in the preliminary gapped and traceback stages (integer).", False),
            _Switch(["-ungapped", "ungapped"], ["input"],
                    "Perform ungapped alignment only?"),
            ]
        _Ncbiblast2SeqCommandline.__init__(self, cmd, **kwargs)

    def _validate(self):
        incompatibles = {"db_soft_mask":["subject", "subject_loc"]}
        for a in incompatibles:
            if self._get_parameter(a):
                for b in incompatibles[a]:
                    if self._get_parameter(b):
                        raise ValueError("Options %s and %s are incompatible." \
                                         % (a,b))
        if (self.template_type and not self.template_length) \
        or (self.template_length and not self.template_type) :
            raise ValueError("Options template_type and template_type require each other.")
        _Ncbiblast2SeqCommandline._validate(self)


class NcbiblastxCommandline(_Ncbiblast2SeqCommandline):
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

    You would typically run the command line with the Python subprocess module,
    as described in the Biopython tutorial.
    """
    def __init__(self, cmd="blastx", **kwargs):
        self.parameters = [ \
            #Input query options:
            _Option(["-strand", "strand"], ["input"],
                    lambda value : value in ["both", "minus", "plus"],0,
                    """Query strand(s) to search against database/subject.

                    Values allowed are "both" (default), "minus", "plus".""", False),
            #Input query options:
            _Option(["-query_gencode", "query_gencode"], ["input"], None, 0,
                    """Genetic code to use to translate query

                    Integer. Default is one.""", False),
            #General search options:
            _Option(["-frame_shift_penalty", "frame_shift_penalty"], ["input"], None, 0,
                    "Frame shift penalty (integer, at least 1, default ignored).", False),
            _Option(["-max_intron_length", "max_intron_length"], ["input"], None, 0,
                    """Maximum intron length (integer).

                    Length of the largest intron allowed in a translated nucleotide
                    sequence when linking multiple distinct alignments (a negative
                    value disables linking). Default zero.""", False),
            _Option(["-matrix", "matrix"], ["input"], None, 0,
                    "Scoring matrix name (default BLOSUM62).", False),
            _Option(["-threshold", "threshold"], ["input"], None, 0,
                    "Minimum word score such that the word is added to the BLAST lookup table (float)", False),
            #Query filtering options:
            _Option(["-seg", "seg"], ["input"], None, 0,
                    """Filter query sequence with SEG (string).

                    Format: "yes", "window locut hicut", or "no" to disable.
                    Default is "12 2.2 2.5""", False),
            #Restrict search or results:
            _Option(["-db_soft_mask", "db_soft_mask"], ["input"], None, 0,
                    """Filtering algorithm for soft masking (integer).

                    Filtering algorithm ID to apply to the BLAST database as soft masking.

                    Incompatible with: subject, subject_loc""", False),
            #Extension options:
            _Switch(["-ungapped", "ungapped"], ["input"],
                    "Perform ungapped alignment only?"),
            ]
        _Ncbiblast2SeqCommandline.__init__(self, cmd, **kwargs)

    def _validate(self):
        incompatibles = {"db_soft_mask":["subject", "subject_loc"]}
        for a in incompatibles:
            if self._get_parameter(a):
                for b in incompatibles[a]:
                    if self._get_parameter(b):
                        raise ValueError("Options %s and %s are incompatible." \
                                         % (a,b))
        _Ncbiblast2SeqCommandline._validate(self)


class NcbitblastnCommandline(_Ncbiblast2SeqCommandline):
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

    You would typically run the command line with the Python subprocess module,
    as described in the Biopython tutorial.
    """
    def __init__(self, cmd="tblastn", **kwargs):
        self.parameters = [ \
            #General search options:
            _Option(["-db_gencode", "db_gencode"], ["input"], None, 0,
                    """Genetic code to use to translate query

                    Integer. Default is one.""", False),
            _Option(["-frame_shift_penalty", "frame_shift_penalty"], ["input"], None, 0,
                    "Frame shift penalty (integer, at least 1, default ignored).", False),
            _Option(["-max_intron_length", "max_intron_length"], ["input"], None, 0,
                    """Maximum intron length (integer).

                    Length of the largest intron allowed in a translated nucleotide
                    sequence when linking multiple distinct alignments (a negative
                    value disables linking). Default zero.""", False),
            _Option(["-matrix", "matrix"], ["input"], None, 0,
                    "Scoring matrix name (default BLOSUM62).", False),
            _Option(["-threshold", "threshold"], ["input"], None, 0,
                    "Minimum word score such that the word is added to the BLAST lookup table (float)", False),
            _Option(["-comp_based_stats", "comp_based_stats"], ["input"],
                    lambda value : value in "0Ft12TtDd3", 0,
                    """Use composition-based statistics (string, default 2, i.e. True).

                    0, F or f: no composition-based statistics
                    1: Composition-based statistics as in NAR 29:2994-3005, 2001
                    2, T or t, D or d : Composition-based score adjustment as in
                       Bioinformatics 21:902-911, 2005, conditioned on sequence properties
                    3: Composition-based score adjustment as in Bioinformatics 21:902-911,
                       2005, unconditionally

                    Note that only tblastn supports values of 1 and 3.""", False),
            #Query filtering options:
            _Option(["-seg", "seg"], ["input"], None, 0,
                    """Filter query sequence with SEG (string).

                    Format: "yes", "window locut hicut", or "no" to disable.
                    Default is "12 2.2 2.5""", False),
            #Extension options:
            _Switch(["-ungapped", "ungapped"], ["input"],
                    "Perform ungapped alignment only?"),
            #Miscellaneous options:
            _Switch(["-use_sw_tback", "use_sw_tback"], ["input"],
                    "Compute locally optimal Smith-Waterman alignments?"),
            #PSI-TBLASTN options:
            _Option(["-in_pssm", "in_pssm"], ["input", "file"], None, 0,
                    """PSI-BLAST checkpoint file

                    Incompatible with: remote, query""", False),
            ]
        _Ncbiblast2SeqCommandline.__init__(self, cmd, **kwargs)

    def _validate(self):
        incompatibles = {"in_pssm":["remote", "query"]}
        for a in incompatibles:
            if self._get_parameter(a):
                for b in incompatibles[a]:
                    if self._get_parameter(b):
                        raise ValueError("Options %s and %s are incompatible." \
                                         % (a,b))
        _Ncbiblast2SeqCommandline._validate(self)


class NcbitblastxCommandline(_Ncbiblast2SeqCommandline):
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

    You would typically run the command line with the Python subprocess module,
    as described in the Biopython tutorial.
    """
    def __init__(self, cmd="tblastx", **kwargs):
        self.parameters = [ \
            #Input query options:
            _Option(["-strand", "strand"], ["input"],
                    lambda value : value in ["both", "minus", "plus"],0,
                    """Query strand(s) to search against database/subject.

                    Values allowed are "both" (default), "minus", "plus".""", False),
            #Input query options:
            _Option(["-query_gencode", "query_gencode"], ["input"], None, 0,
                    """Genetic code to use to translate query

                    Integer. Default is one.""", False),
            #General search options:
            _Option(["-db_gencode", "db_gencode"], ["input"], None, 0,
                    """Genetic code to use to translate query

                    Integer. Default is one.""", False),
            _Option(["-max_intron_length", "max_intron_length"], ["input"], None, 0,
                    """Maximum intron length (integer).

                    Length of the largest intron allowed in a translated nucleotide
                    sequence when linking multiple distinct alignments (a negative
                    value disables linking). Default zero.""", False),
            _Option(["-matrix", "matrix"], ["input"], None, 0,
                    "Scoring matrix name (default BLOSUM62).", False),
            _Option(["-threshold", "threshold"], ["input"], None, 0,
                    "Minimum word score such that the word is added to the BLAST lookup table (float)", False),
            #Query filtering options:
            _Option(["-seg", "seg"], ["input"], None, 0,
                    """Filter query sequence with SEG (string).

                    Format: "yes", "window locut hicut", or "no" to disable.
                    Default is "12 2.2 2.5""", False),
           ]
        _Ncbiblast2SeqCommandline.__init__(self, cmd, **kwargs)


    def _validate(self):
        if self.remote and self.in_pssm:
            raise ValueError("The remote option cannot be used with in_pssm")
        if self.query and self.in_pssm:
            raise ValueError("The query option cannot be used with in_pssm")
        _Ncbiblast2SeqCommandline._validate(self)


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

    You would typically run the command line with the Python subprocess module,
    as described in the Biopython tutorial.
    """
    def __init__(self, cmd="psiblast", **kwargs):
        self.parameters = [ \
            #General search options:
            _Option(["-matrix", "matrix"], ["input"], None, 0,
                    "Scoring matrix name (default BLOSUM62).", False),
            _Option(["-threshold", "threshold"], ["input"], None, 0,
                    "Minimum word score such that the word is added to the BLAST lookup table (float)", False),
            _Option(["-comp_based_stats", "comp_based_stats"], ["input"],
                    lambda value : value in "0Ft2TtDd", 0,
                    """Use composition-based statistics (string, default 2, i.e. True).

                    0, F or f: no composition-based statistics
                    2, T or t, D or d : Composition-based score adjustment as in
                    Bioinformatics 21:902-911, 2005, conditioned on sequence properties

                    Note that tblastn also supports values of 1 and 3.""", False),
            #Query filtering options:
            _Option(["-seg", "seg"], ["input"], None, 0,
                    """Filter query sequence with SEG (string).

                    Format: "yes", "window locut hicut", or "no" to disable.
                    Default is "12 2.2 2.5""", False),
            #Extension options:
            _Option(["-gap_trigger", "gap_trigger"], ["input"], None, 0,
                    "Number of bits to trigger gapping (float, default 22)", False),
            #Miscellaneous options:
            _Switch(["-use_sw_tback", "use_sw_tback"], ["input"],
                    "Compute locally optimal Smith-Waterman alignments?"),
            #PSI-BLAST options:
            _Option(["-num_iterations", "num_iterations"], ["input"], None, 0,
                    """Number of iterations to perform, integer

                    Integer of at least one. Default is one.
                    Incompatible with: remote""", False),
            _Option(["-out_pssm", "out_pssm"], ["output", "file"], None, 0,
                    "File name to store checkpoint file", False),
            _Option(["-out_ascii_pssm", "out_ascii_pssm"], ["output", "file"], None, 0,
                    "File name to store ASCII version of PSSM", False),
            _Option(["-in_msa", "in_msa"], ["input", "file"], None, 0,
                    """File name of multiple sequence alignment to restart PSI-BLAST

                    Incompatible with: in_pssm, query""", False),
            _Option(["-in_pssm", "in_pssm"], ["input", "file"], None, 0,
                    """PSI-BLAST checkpoint file

                    Incompatible with: in_msa, query, phi_pattern""", False),
            #PSSM engine options:
            _Option(["-pseudocount", "pseudocount"], ["input"], None, 0,
                    """Pseudo-count value used when constructing PSSM

                    Integer. Default is zero.""", False),
            _Option(["-inclusion_ethresh", "inclusion_ethresh"], ["input"], None, 0,
                    """E-value inclusion threshold for pairwise alignments

                    Float. Default is 0.002.""", False),
            #PHI-BLAST options:
            _Option(["-phi_pattern", "phi_pattern"], ["input", "file"], None, 0,
                    """File name containing pattern to search

                    Incompatible with: in_pssm""", False),
            ]
        _Ncbiblast2SeqCommandline.__init__(self, cmd, **kwargs)

    def _validate(self):
        incompatibles = {"num_iterations":["remote"],
                         "in_msa":["in_pssm", "query"],
                         "in_pssm":["in_msa","query","phi_pattern"]}
        for a in incompatibles:
            if self._get_parameter(a):
                for b in incompatibles[a]:
                    if self._get_parameter(b):
                        raise ValueError("Options %s and %s are incompatible." \
                                         % (a,b))
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

    You would typically run the command line with the Python subprocess module,
    as described in the Biopython tutorial.
    """
    def __init__(self, cmd="rpsblast", **kwargs):
        self.parameters = [ \
            #Query filtering options:
            _Option(["-seg", "seg"], ["input"], None, 0,
                    """Filter query sequence with SEG (string).

                    Format: "yes", "window locut hicut", or "no" to disable.
                    Default is "12 2.2 2.5""", False),
            ]
        _NcbiblastCommandline.__init__(self, cmd, **kwargs)


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

    You would typically run the command line with the Python subprocess module,
    as described in the Biopython tutorial.
    """
    def __init__(self, cmd="rpstblastn", **kwargs):
        self.parameters = [ \
            #Input query options:
            _Option(["-strand", "strand"], ["input"],
                    lambda value : value in ["both", "minus", "plus"],0,
                    """Query strand(s) to search against database/subject.

                    Values allowed are "both" (default), "minus", "plus".""", False),
            #Input query options:
            _Option(["-query_gencode", "query_gencode"], ["input"], None, 0,
                    """Genetic code to use to translate query

                    Integer. Default is one.""", False),
            #Query filtering options:
            _Option(["-seg", "seg"], ["input"], None, 0,
                    """Filter query sequence with SEG (string).

                    Format: "yes", "window locut hicut", or "no" to disable.
                    Default is "12 2.2 2.5""", False),
            #Extension options:
            _Switch(["-ungapped", "ungapped"], ["input"],
                    "Perform ungapped alignment only?"),
            ]
        _NcbiblastCommandline.__init__(self, cmd, **kwargs)


def _test():
    """Run the Bio.Blast.Applications module's doctests."""
    import doctest
    doctest.testmod(verbose=1)

if __name__ == "__main__":
    #Run the doctests
    _test()
