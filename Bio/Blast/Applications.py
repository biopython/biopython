"""Definitions for interacting with Blast related applications.
"""
from Bio import Application
from Bio.Application import _Option

class FastacmdCommandline(Application.AbstractCommandline):
    """Create a commandline for the fasta program from NCBI.

    """

    def __init__(self, fastacmd = "fastacmd"):
        Application.AbstractCommandline.__init__(self)
        
        self.program_name = fastacmd

        self.parameters = \
          [
           _Option(["-d", "database"], ["input"], None, 1,
                   "The database to retrieve from."),
           _Option(["-s", "search_string"], ["input"], None, 1,
                   "The id to search for.")
          ]
  
class BlastallCommandline(Application.AbstractCommandline):
    """Create a commandline for the blastall program from NCBI.

    XXX This could use more checking for valid paramters to the program.
    """
    def __init__(self, blastcmd = "blastall"):
        Application.AbstractCommandline.__init__(self)

        self.program_name = blastcmd

        self.parameters = \
          [# Scoring options
           _Option(["-M", "matrix"], ["input"], None, 0, 
                   "Matrix to use"),
           _Option(["-G", "gap_open"], ["input"], None, 0, 
                   "Gap open penalty"),
           _Option(["-E", "gap_extend"], ["input"], None, 0, 
                    "Gap extension penalty"),
           _Option(["-A", "window_size"], ["input"], None, 0,
                    "Multiple hits window size"),
           _Option(["-j", "npasses"], ["input"], None, 0,
                    "Number of passes"),
           _Option(["-p", "passes"], ["input"], None, 0,
                   "Hits/passes.  Integer 0-2."),

            # Algorithm options
           _Option(["-g", "gapped"], ["input"], None, 0, 
                   "Whether to do a gapped alignment.  T/F"),
           _Option(["-e", "expectation"], ["input"], None, 0, 
                   "Expectation value cutoff."),
           _Option(["-W", "wordsize"], ["input"], None, 0, 
                   "Word size"),
           _Option(["-K", "keep_hits"], ["input"], None, 0,
                   " Number of best hits from a region to keep."),
           _Option(["-X", "xdrop"], ["input"], None, 0, 
                   "Dropoff value (bits) for gapped alignments."),
           _Option(["-f", "hit_extend"], ["input"], None, 0, 
                   "Threshold for extending hits."),
           _Option(["-L", "region_length"], ["input"], None, 0, 
                   "Length of region used to judge hits."),
           _Option(["-Z", "db_length"], ["input"], None, 0, 
                   "Effective database length."),
           _Option(["-Y", "search_length"], ["input"], None, 0, 
                   "Effective length of search space."),
           _Option(["-N", "nbits_gapping"], ["input"], None, 0, 
                   "Number of bits to trigger gapping."),
           _Option(["-c", "pseudocounts"], ["input"], None, 0,
                   "Pseudocounts constants for multiple passes."),
           _Option(["-Z", "xdrop_final"], ["input"], None, 0, 
                   "X dropoff for final gapped alignment."),
           _Option(["-y", "xdrop_extension"], ["input"], None, 0, 
                   "Dropoff for blast extensions."),
           _Option(["-h", "model_threshold"], ["input"], None, 0, 
                   "E-value threshold to include in multipass model."),
           _Option(["-S", "required_start"], ["input"], None, 0, 
                   "Start of required region in query."),
           _Option(["-H", "required_end"], ["input"], None, 0,
                   "End of required region in query."),

            # Processing options
           _Option(["-p", "program"], ["input"], None, 1, 
                   "The blast program to use."),
           _Option(["-d", "database"], ["input"], None, 1,
                   "The database to BLAST against."),
           _Option(["-i", "infile"], ["input", "file"], None, 1,
                   "The sequence to search with."),
           _Option(["-F", "filter"], ["input"], None, 0,
                   "Filter query sequence with SEG?  T/F"),
           _Option(["-J", "believe_query"], ["input"], None, 0,
                   "Believe the query defline?  T/F"),
           _Option(["-a", "nprocessors"], ["input"], None, 0,
                   "Number of processors to use."),

           # Formatting options
           _Option(["-T", "html"], ["input"], None, 0, 
                   "Produce HTML output?  T/F"),
           _Option(["-v", "descriptions"], ["input"], None, 0, 
                   "Number of one-line descriptions."),
           _Option(["-b", "alignments"], ["input"], None, 0, 
                   "Number of alignments."),
           _Option(["-m", "align_view"], ["input"], None, 0, 
                   "Alignment view.  Integer 0-6."),
           _Option(["-I", "show_gi"], ["input"], None, 0, 
                   "Show GI's in deflines?  T/F"),
           _Option(["-O", "seqalign_file"], ["output", "file"], None, 0,
                   "seqalign file to output."),
           _Option(["-o", "align_outfile"], ["output", "file"], None, 1,
                   "Output file for alignment."),
           _Option(["-C", "checkpoint_outfile"], ["output", "file"], None, 0,
                   "Output file for PSI-BLAST checkpointing."),
           _Option(["-R", "restart_infile"], ["input", "file"], None, 0,
                   "Input file for PSI-BLAST restart."),
           _Option(["-k", "hit_infile"], ["input", "file"], None, 0,
                   "Hit file for PHI-BLAST."),
           _Option(["-Q", "matrix_outfile"], ["output", "file"], None, 0,
                   "Output file for PSI-BLAST matrix in ASCII."),
           _Option(["-B", "align_infile"], ["input", "file"], None, 0, 
                   "Input alignment file for PSI-BLAST restart.")
          ] 
