# -*- coding: utf-8 -*-
# Copyright 2017 by Rob Gilmore and Shaurita Hutchins. All rights reserved.
# Based on ClustalOmega wrapper copyright 2011 by Andreas Wilm.
#
# Wrapper for Guidance2 by Rob Gilmore (2017).  http://guidance.tau.ac.il/ver2/
# Used _ClustalOmega.py as template.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Command line wrapper for the multiple alignment program, GUIDANCE2.
It weights, filters or masks unreliably aligned positions in multiple sequence alignments.
"""

from __future__ import print_function
from pathlib import Path
from Bio.Application import _Option, AbstractCommandline


class Guidance2Commandline(AbstractCommandline):
    u""""Command line wrapper for GUIDANCE2.

    http://guidance.tau.ac.il/ver2/

    Example:
    --------

    >>> from Bio.Align.Applications import Guidance2Commandline
    >>> import os
    >>> seqFile = "HTR1A.ffn"
    >>> msaProgram = "CLUSTALW"
    >>> seqType = "codon"
    >>> outDir = Path(os.getcwd)
    >>> Guidance2_cline = Guidance2Commandline(seqFile, msaProgram, seqType, str(outDir), bootstraps=20, seqCutoff=0.63,\
                                    colCutoff=0.9, outOrder='as_input', dataset='HTR1A')
    >>>> print(Guidance2_cline)
    perl guidance.pl --seqFile HTR1A.ffn --msaProgram CLUSTALW --seqType codon --outDir C://Users/rgilmore/HTR1A
    --bootstraps 20 --seqCutoff 0.63 --colCutoff 0.9 --outOrder as_input --dataset HTR1A


    You would typically run the command line with Guidance2_cline() or via
        the Python subprocess module, as described in the Biopython tutorial.

    Citation:
    ---------
        Sela, I., Ashkenazy, H., Katoh, K. and Pupko, T. (2015)
        GUIDANCE2: accurate detection of unreliable alignment regions accounting for the uncertainty of multiple parameters.
        Nucleic Acids Research, 2015 Jul 1; 43 (Web Server issue): W7-W14.; doi: 10.1093/nar/gkq443

        Landan, G., and D. Graur. (2008).
        Local reliability measures from sets of co-optimal multiple sequence alignments.
        Pac Symp Biocomput 13:15-24
    """

    def __init__(self, cmd="perl guidance.pl", **kwargs):
            # TODO-ROB:  Add command lines for the alternative scripts in the guidance package
            self.parameters = \
                [
                    # Required Parameters
                    _Option(['--seqFile', 'seqFile'],
                            "Input sequence file in FASTA format",
                            filename=True, equate=False, is_required=True,
                            checker_function=lambda x: str(Path(x).suffix) in ['.fasta', 'fna', '.ffn', '.faa', '.fra']
                                                       and Path(x).is_file()),
                    _Option(['--msaProgram', 'msaProgram'],
                            "Which MSA program to use",
                            equate=False, is_required=True,
                            checker_function=lambda x: x in ['MAFFT', 'PRANK', 'CLUSTALW', 'MUSCLE']),
                    _Option(['--seqType', 'seqType'],
                            "Type of sequences for alignment (amino acids, nucleotides, or codons)",
                            equate=False, is_required=True,
                            checker_function=lambda x: x in ['aa', 'nuc', 'codon']),
                    _Option(['--outDir', 'outDir'],
                            "Output directory that will be created "
                            "automatically and hold all output files [please provid full (and not relative) path]",
                            filename=True, equate=False, is_required=True),

                    # Optional Parameters
                    _Option(['--program', 'program'],
                            "[GUIDANCE2|GUIDANCE|HoT] Default=GUIDANCE2",
                            equate=False,
                            checker_function=lambda x: x in ["GUIDANCE2", "GUIDANCE", "HoT"]),
                    _Option(['--bootstraps', 'bootstraps'],
                            "Number of bootstrap iterations (only for GUIDQANCE). Defaut=100",
                            equate=False,
                            checker_function=lambda x: isinstance(x, int)),
                    _Option(['--genCode', 'genCode'],
                            "Genetic code identifier (only for codon sequences). Default=1 \
                                1) Nuclear Standard\
                                15) Nuclear Blepharisma\
                                6) Nuclear Ciliate\
                                10) Nuclear Euplotid\
                                2) Mitochondria Vertebrate\
                                5) Mitochondria Invertebrate\
                                3) Mitochondria Yeast\
                                13) Mitochondria Ascidian\
                                9) Mitochondria Echinoderm\
                                14) Mitochondria Flatworm\
                                4) Mitochondria Protozoan",
                            equate=False,
                            checker_function=lambda x: isinstance(x, int)),
                    _Option(['--outOrder', 'outOrder'],
                            "[aligned|as_input] default=aligned",
                            equate=False,
                            checker_function=lambda x: x in ['aligned', 'as_input']),
                    _Option(['--msaFile', 'msaFile'],
                            "Input alignment file - not recommended",
                            filename=True, equate=False,
                            checker_function=lambda x: Path(x).is_file()),
                    # Confidence scores
                    _Option(['--seqCutoff', 'seqCutoff'],
                            "Confidence cutoff between 0 to 1. Default=0.6",
                            equate=False,
                            checker_function=lambda x: isinstance(x, (int, float))),
                    _Option(['--colCutoff', 'colCutoff'],
                            "Confidence cutoff between 0 to 1. Default=0.93",
                            equate=False,
                            checker_function=lambda x: isinstance(x, (int, float))),
                    # Alignment Programs
                    _Option(['--mafft', 'mafft'],
                            "path to mafft executable. Default=mafft",
                            filename=True, equate=False,
                            checker_function=lambda x: Path(x).is_file()),
                    _Option(['--prank', 'prank'],
                            "path to prank executable. Default=prank",
                            filename=True, equate=False,
                            checker_function=lambda x: Path(x).is_file()),
                    _Option(['--muscle', 'muscle'],
                            "path to muscle executable. default=muscle",
                            filename=True, equate=False,
                            checker_function=lambda x: Path(x).is_file()),
                    _Option(['--pagan', 'pagan'],
                            "path to pagan executable, default=pagan",
                            filename=True, equate=False,
                            checker_function=lambda x: Path(x).is_file()),
                    _Option(['--ruby', 'ruby'],
                            "path to ruby executable. default=ruby",
                            filename=True, equate=False,
                            checker_function=lambda x: Path(x).is_file()),
                    # Miscellaneous
                    _Option(['--dataset', 'dataset'],
                            "Unique name for the Dataset - will be used as prefix to outputs (default=MSA)",
                            equate=False),
                    _Option(['--MSA_Param', 'MSA_Param'],
                            "passing parameters for the alignment program e.g -F to prank. "
                            "To pass parameter containning '-' in it, add \ before each '-' e.g. \-F for PRANK",
                            equate=False),
                    _Option(['--proc_num', 'proc_num'],
                            "number of processors to use (default=1)",
                            equate=False,
                            checker_function=lambda x: isinstance(x, int))
                ]
            AbstractCommandline.__init__(self, cmd, **kwargs)


if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()
