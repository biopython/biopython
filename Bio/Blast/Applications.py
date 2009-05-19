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

  
class BlastallCommandline(AbstractCommandline):
    """Create a commandline for the blastall program from NCBI."""
    #TODO - This could use more checking for valid parameters to the program.
    def __init__(self, cmd="blastall",**kwargs):
        self.parameters = [ \
            #Sorted in the same order as the output from blastall --help
            #which should make it easier to keep them up to date in future.
           _Option(["-p", "program"], ["input"], None, 1, 
                   "The blast program to use."),
           _Option(["-d", "database"], ["input"], None, 1,
                   "The database to BLAST against."),
           _Option(["-i", "infile"], ["input", "file"], None, 1,
                   "The sequence to search with."),
           _Option(["-e", "expectation"], ["input"], None, 0, 
                   "Expectation value cutoff."),
           _Option(["-m", "align_view"], ["input"], None, 0, 
                   "Alignment view.  Integer 0-11.  Use 7 for XML output."),
           _Option(["-o", "align_outfile", "outfile"], ["output", "file"], None, 1,
                   "Output file for alignment."),
           _Option(["-F", "filter"], ["input"], None, 0,
                   "Filter query sequence with SEG?  T/F"),
           _Option(["-G", "gap_open"], ["input"], None, 0, 
                   "Gap open penalty"),
           _Option(["-E", "gap_extend"], ["input"], None, 0, 
                    "Gap extension penalty"),
           _Option(["-X", "xdrop"], ["input"], None, 0, 
                   "Dropoff value (bits) for gapped alignments."),
           _Option(["-I", "show_gi"], ["input"], None, 0, 
                   "Show GI's in deflines?  T/F"),
           _Option(["-q", "nuc_mismatch"], ["input"], None, 0, 
                   "Penalty for a nucleotide mismatch (blastn only)."),
           _Option(["-r", "nuc_match"], ["input"], None, 0, 
                   "Reward for a nucleotide match (blastn only)."),
           _Option(["-v", "descriptions"], ["input"], None, 0, 
                   "Number of one-line descriptions."),
           _Option(["-b", "alignments"], ["input"], None, 0, 
                   "Number of alignments."),
           _Option(["-f", "hit_extend"], ["input"], None, 0, 
                   "Threshold for extending hits."),
           _Option(["-g", "gapped"], ["input"], None, 0, 
                   "Whether to do a gapped alignment.  T/F"),
           _Option(["-Q", "query_genetic_code"], ["input"], None, 0,
                   "Query Genetic code to use."),
           _Option(["-D", "db_genetic_code"], ["input"], None, 0,
                   "DB Genetic code (for tblast[nx] only)."),
           _Option(["-a", "nprocessors"], ["input"], None, 0,
                   "Number of processors to use."),
           _Option(["-O", "seqalign_file"], ["output", "file"], None, 0,
                   "seqalign file to output."),
           _Option(["-J", "believe_query"], ["input"], None, 0, 
                   "Believe the query defline?  T/F"),
           _Option(["-M", "matrix"], ["input"], None, 0, 
                   "Matrix to use"),
           _Option(["-W", "wordsize"], ["input"], None, 0, 
                   "Word size"),
           _Option(["-z", "db_length"], ["input"], None, 0, 
                   "Effective database length."),
           _Option(["-K", "keep_hits"], ["input"], None, 0,
                   " Number of best hits from a region to keep."),
           _Option(["-P", "passes"], ["input"], None, 0,
                   "Hits/passes.  Integer 0-2. 0 for multiple hit, "
                   "1 for single hit (does not apply to blastn)"),
           _Option(["-Y", "search_length"], ["input"], None, 0, 
                   "Effective length of search space (use zero for the real size)."),
           _Option(["-S", "strands"], ["input"], None, 0, 
                   "Query strands to search against database (for blast[nx], " + \
                   "and tblastx). 3 is both, 1 is top, 2 is bottom."),
           _Option(["-T", "html"], ["input"], None, 0, 
                   "Produce HTML output?  T/F"),
           _Option(["-l", "restrict_gi"], ["input"], None, 0,
                   "Restrict search of database to list of GI's."),
           _Option(["-U"], ["input"], None, 0,
                   "Use lower case filtering of FASTA sequence? T/F"),
           _Option(["-y", "xdrop_extension"], ["input"], None, 0, 
                   "Dropoff for blast extensions."),
           _Option(["-Z", "xdrop_final"], ["input"], None, 0, 
                   "X dropoff for final gapped alignment."),
           _Option(["-R"], ["input", "file"], None, 0,
                   "PSI-TBLASTN checkpoint input file."),
           _Option(["-n", "megablast"], ["input"], None, 0,
                   "MegaBlast search T/F."),
           _Option(["-L", "region_length"], ["input"], None, 0, 
                   "Length of region used to judge hits."),
           _Option(["-A", "window_size"], ["input"], None, 0,
                    "Multiple hits window size"),
           _Option(["-w"], ["input"], None, 0,
                   "Frame shift penalty (OOF algorithm for blastx)."),
           _Option(["-t"], ["input"], None, 0,
                   "Length of the largest intron allowed in a translated " + \
                   "nucleotide sequence when linking multiple distinct " + \
                   "alignments. (0 invokes default behavior; a negative value " + \
                   "disables linking.)"),
           _Option(["-B"], ["input"], None, 0,
                   "Number of concatenated queries, for blastn and tblastn."),
           _Option(["-V", "oldengine"], ["input"], None, 0,
                   "Force use of the legacy BLAST engine."),
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
                   D, F or 0."""),
           _Option(["-s"], ["input"], None, 1,
                   "Compute locally optimal Smith-Waterman alignments (This " + \
                   "option is only available for gapped tblastn.) T/F"),
           #These are blastpgp options - not blastall options!
           #_Option(["-C", "checkpoint_outfile"], ["output", "file"], None, 0,
           #        "Output file for PSI-BLAST checkpointing."),
           #_Option(["-R", "restart_infile"], ["input", "file"], None, 0,
           #        "Input file for PSI-BLAST restart."),
           #_Option(["-k", "hit_infile"], ["input", "file"], None, 0,
           #        "Hit file for PHI-BLAST."),
           #_Option(["-Q", "matrix_outfile"], ["output", "file"], None, 0,
           #        "Output file for PSI-BLAST matrix in ASCII."),
           #_Option(["-B", "align_infile"], ["input", "file"], None, 0, 
           #        "Input alignment file for PSI-BLAST restart.")
           #_Option(["-S", "required_start"], ["input"], None, 0, 
           #        "Start of required region in query."),
           #_Option(["-H", "required_end"], ["input"], None, 0,
           #        "End of required region in query."),
           #_Option(["-j", "npasses"], ["input"], None, 0,
           #         "Number of passes"),
           #_Option(["-N", "nbits_gapping"], ["input"], None, 0, 
           #        "Number of bits to trigger gapping."),
           #_Option(["-c", "pseudocounts"], ["input"], None, 0,
           #        "Pseudocounts constants for multiple passes."),
           #_Option(["-h", "model_threshold"], ["input"], None, 0, 
           #        "E-value threshold to include in multipass model."),
          ] 
        AbstractCommandline.__init__(self, cmd, **kwargs)

