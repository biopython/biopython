# Copyright 2001 Brad Chapman.
# Revisions copyright 2009 by Peter Cock.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Definitions for interacting with Blast related applications.
"""
from Bio.Application import _Option, AbstractCommandline

class FastacmdCommandline(AbstractCommandline):
    """Create a commandline for the fasta program from NCBI.

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
    """Base Commandline object for NCBI BLAST wrappers (PRIVATE).

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
    """Create a commandline for the blastall program from NCBI."""
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
    """Create a commandline for the blastpgp program from NCBI."""
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
    """Create a commandline for the rpsblast program from NCBI."""
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
