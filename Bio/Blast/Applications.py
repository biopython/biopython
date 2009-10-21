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

- NcbiblastpCommandline
- NcbiblastnCommandline
- NcbiblastxCommandline
- NcbitblastnCommandline
- NcbitblastxCommandline
- Others pending

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


class _BlastCommandLine(AbstractCommandline) :
    """Base Commandline object for (classic) NCBI BLAST wrappers (PRIVATE).

    This is provided for subclassing, it deals with shared options
    common to all the BLAST tools (blastall, rpsblast, pgpblast).
    """
    def __init__(self, cmd=None, **kwargs):
        assert cmd is not None
        extra_parameters = [\
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
        try :
            #Insert extra parameters - at the start just in case there
            #are any arguments which must come last:
            self.parameters = extra_parameters + self.parameters
        except AttributeError:
            #Should we raise an error?  The subclass should have set this up!
            self.parameters = extra_parameters
        AbstractCommandline.__init__(self, cmd, **kwargs)


class _BlastAllOrPgpCommandLine(_BlastCommandLine) :
    """Base Commandline object for NCBI BLAST wrappers (PRIVATE).

    This is provided for subclassing, it deals with shared options
    common to all the blastall and pgpblast tools (but not rpsblast).
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
        try :
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
    are replacing blastpgp with a renamed tool psiblast.
    
    Like blastpgp (and blastall), this wrapper is now obsolete, and will be
    deprecated and removed in a future release of Biopython.
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
    are replacing the old rpsblast with a new version of the same name, taking
    different command line arguments.
    
    Like the old rpsblast (and blastall), this wrapper is now obsolete, and will be
    deprecated and removed in a future release of Biopython."""
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

#Further classes in the heirachy will probably be needed...
class _NcbiblastCommandline(AbstractCommandline) :
    """Base Commandline object for (classic) NCBI BLAST wrappers (PRIVATE).

    This is provided for subclassing, it deals with shared options
    common to all the BLAST tools (blastall, rpsblast, pgpblast).
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
            #BLAST-2-Sequences options:
            _Option(["-subject", "subject"], ["input", "file"], None, 0,
                    """Subject sequence(s) to search.

                    Incompatible with:  db, gilist, negative_gilist.
                    See also subject_loc.""", False),
            _Option(["-subject_loc", "subject_loc"], ["input"], None, 0,
                    """Location on the subject sequence (Format: start-stop)

                    Incompatible with:  db, gilist, negative_gilist, remote.
                    See also subject.""", False),
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
            _Switch(["-lcase_masking", "lcase_masking"], ["input"],
                    "Use lower case filtering in query and subject sequence(s)?"),
            #Restrict search or results
            _Option(["-gilist", "gilist"], ["input", "file"], None, 0,
                    """Restrict search of database to list of GI's.
 
                    Incompatible with:  negative_gilist, remote, subject, subject_loc""",
                    False),
            _Option(["-negative_gilist", "negative_gilist"], ["input", "file"], None, 0,
                    """Restrict search of database to everything except the listed GIs.
 
                    Incompatible with:  gilist, remote, subject, subject_loc""",
                    False),
            #Statistical options
            _Option(["-dbsize", "dbsize"], ["input"], None, 0,
                    "Effective length of the database (integer)", False),
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
            #Miscellaneous options
            _Option(["-num_threads", "num_threads"], ["input"], None, 0,
                    """Number of threads to use in the BLAST search.

                    Integer of at least one. Default is one.
                    Incompatible with:  remote""", False),
            _Switch(["-remote", "remote"], ["input"],
                    """Execute search remotely?

                    Incompatible with:  gilist, negative_gilist, subject_loc, num_threads, ..."""),
            ]
        try :
            #Insert extra parameters - at the start just in case there
            #are any arguments which must come last:
            self.parameters = extra_parameters + self.parameters
        except AttributeError:
            #Should we raise an error?  The subclass should have set this up!
            self.parameters = extra_parameters
        AbstractCommandline.__init__(self, cmd, **kwargs)

    def _validate(self) :
        incompatibles = {"remote":["gilist", "negative_gilist", "subject_loc", "num_threads"],
                         "subject":["db", "gilist", "negative_gilist"],
                         "gilist":["negative_gilist"]}
        for a in incompatibles :
            if self._get_parameter(a) :
                for b in incompatibles[a] :
                    if self._get_parameter(b) :
                        raise ValueError("Options %s and %s are incompatible." \
                                         % (a,b))
        AbstractCommandline._validate(self)

class NcbiblastpCommandline(_NcbiblastCommandline) :
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
            _Switch(["-ungapped", "ungapped"], ["input"],
                    "Perform ungapped alignment only?"),
            ]
        _NcbiblastCommandline.__init__(self, cmd, **kwargs)

class NcbiblastnCommandline(_NcbiblastCommandline) :
    """Wrapper for the NCBI BLAST+ program blastn (for nucleotides).

    With the release of BLAST+ (BLAST rewritten in C++ instead of C), the NCBI
    replaced the old blastall tool with separate tools for each of the searches.
    This wrapper therefore replaces BlastallCommandline with option -p blastn.

    For example, to run a search against the "nt" nucleotide database using the
    FASTA nucleotide file "m_code.fasta" as the query, with an expectation value
    cut off of 0.001, saving the output to a file in XML format:

    >>> from Bio.Blast.Applications import NcbiblastnCommandline
    >>> cline = NcbiblastnCommandline(query="m_cold.fasta", db="nt",
    ...                               evalue=0.001, out="m_cold.xml", outfmt=5)
    >>> cline
    NcbiblastnCommandline(cmd='blastn', query='m_cold.fasta', db='nt', out='m_cold.xml', evalue=0.001, outfmt=5)
    >>> print cline
    blastn -query m_cold.fasta -db nt -out m_cold.xml -evalue 0.001 -outfmt 5

    You would typically run the command line with the Python subprocess module,
    as described in the Biopython tutorial.
    """
    def __init__(self, cmd="blastn", **kwargs):
        self.parameters = [ \
            _Switch(["-ungapped", "ungapped"], ["input"],
                    "Perform ungapped alignment only?"),
            ]
        _NcbiblastCommandline.__init__(self, cmd, **kwargs)


class NcbiblastxCommandline(_NcbiblastCommandline) :
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
            _Switch(["-ungapped", "ungapped"], ["input"],
                    "Perform ungapped alignment only?"),
            ]
        _NcbiblastCommandline.__init__(self, cmd, **kwargs)


class NcbitblastnCommandline(_NcbiblastCommandline) :
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
            _Switch(["-ungapped", "ungapped"], ["input"],
                    "Perform ungapped alignment only?"),
            ]
        _NcbiblastCommandline.__init__(self, cmd, **kwargs)


class NcbitblastxCommandline(_NcbiblastCommandline) :
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
            #PSI-TBLASTN options:
            _Option(["-in_pssm", "in_pssm"], ["input", "file"], None, 0,
                    """PSI-TBLASTN checkpoint file

                    Incompatible with:  remote, query."""),
            ]
        _NcbiblastCommandline.__init__(self, cmd, **kwargs)


    def _validate(self) :
        if self.remote and self.in_pssm :
            raise ValueError("The remote option cannot be used with in_pssm")
        if self.query and self.in_pssm :
            raise ValueError("The query option cannot be used with in_pssm")
        _NcbiblastCommandline._validate(self)


def _test():
    """Run the Bio.Blast.Applications module's doctests."""
    import doctest
    doctest.testmod(verbose=1)

if __name__ == "__main__":
    #Run the doctests
    _test()
