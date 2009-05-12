# Copyright 2009 by Cymon J. Cox.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Command line wrapper for the multiple alignment program PROBCONS.

http://probcons.stanford.edu/

Citations:
Do, C.B., Mahabhashyam, M.S.P., Brudno, M., and Batzoglou, S. 2005. PROBCONS:
Probabilistic Consistency-based Multiple Sequence Alignment. Genome Research 15:
330-340. 

Last checked agains version: 1.12
"""
import types
from Bio.Application import _Option, _Switch, _Argument, AbstractCommandline

class ProbconsCommandline(AbstractCommandline):
    """Command line wrapper for the multiple alignment program PROBCONS."""
    def __init__(self, cmd="probcons", **kwargs):
        self.parameters = \
            [
            #Note that some options cannot be assigned via properties using the
            #original documented option:
            #e.g cmdline.pre-training = 3 will not work
            #In these cases the shortened option name should be used
            #cmdline.pre = 3
            
            #-clustalw
            #       use CLUSTALW output format instead of MFA
            _Switch(["-clustalw", "clustalw"], ["input"],
                    "Use CLUSTALW output format instead of MFA"),
            #-c, --consistency REPS
            #       use 0 <= REPS <= 5 (default: 2) passes of consistency transformation
            _Option(["-c", "c", "--consistency", "consistency" ], ["input"],
                    lambda x: x in range(0,6),
                    0,
                    "Use 0 <= REPS <= 5 (default: 2) passes of consistency transformation",
                    0),
            #-ir, --iterative-refinement REPS
            #       use 0 <= REPS <= 1000 (default: 100) passes of iterative-refinement
            _Option(["-ir", "--iterative-refinement", "iterative-refinement", "ir"], ["input"],
                    lambda x: x in range(0,1001),
                    0,
                    "Use 0 <= REPS <= 1000 (default: 100) passes of iterative-refinement",
                    0),
            #-pre, --pre-training REPS
            #       use 0 <= REPS <= 20 (default: 0) rounds of pretraining
            _Option(["-pre", "--pre-training", "pre-training", "pre"], ["input"],
                    lambda x: x in range(0,21),
                    0,
                    "Use 0 <= REPS <= 20 (default: 0) rounds of pretraining",
                    0),
            #-pairs
            #       generate all-pairs pairwise alignments
            _Switch(["-pairs", "pairs"], ["input"],
                    "Generate all-pairs pairwise alignments"),
            #-viterbi
            #       use Viterbi algorithm to generate all pairs (automatically enables -pairs)
            _Switch(["-viterbi", "viterbi"], ["input"],
                    "Use Viterbi algorithm to generate all pairs (automatically enables -pairs)"),
            #-v, --verbose
            #       report progress while aligning (default: off)
            _Switch(["-verbose", "verbose"], ["input"],
                    "Report progress while aligning (default: off)"),
            #-annot FILENAME
            #       write annotation for multiple alignment to FILENAME
            _Option(["-annot", "annot"], ["input"],
                    None,
                    0,
                    "Write annotation for multiple alignment to FILENAME",
                    0),
            #-t, --train FILENAME
            #       compute EM transition probabilities, store in FILENAME (default: no training)
            _Option(["-t", "t", "--train", "train"], ["input"],
                    None,
                    0,
                    "Compute EM transition probabilities, store in FILENAME (default: no training)",
                    0),
            #-e, --emissions
            #       also reestimate emission probabilities (default: off)
            _Switch(["-e", "e", "--emissions", "emissions"], ["input"],
                    "Also reestimate emission probabilities (default: off)"),
            #-p, --paramfile FILENAME
            #       read parameters from FILENAME (default: )
            _Option(["-p", "p", "--paramfile", "paramfile"], ["input"],
                    None,
                    0,
                    "Read parameters from FILENAME (default: )",
                    0),
            #-a, --alignment-order
            #       print sequences in alignment order rather than input order (default: off)
            _Switch(["-a", "--alignment-order", "alignment-order", "a"], ["input"],
                    "Print sequences in alignment order rather than input order (default: off)"),
            #Input file name
            _Argument(["input"], ["input", "file"], None, 1,
                      "Input file name. Must be multiple FASTA alignment "+ \
                      "(MFA) format"),
            ]
        AbstractCommandline.__init__(self, cmd, **kwargs)
