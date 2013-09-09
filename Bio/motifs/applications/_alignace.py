# Copyright 2003-2009 by Bartek Wilczynski.  All rights reserved.
# Revisions copyright 2009 by Peter Cock.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""This module provides code to work with the standalone version of AlignACE,
for motif search in DNA sequences.

AlignACE homepage:

http://arep.med.harvard.edu/mrnadata/mrnasoft.html

AlignACE Citations:

Computational identification of cis-regulatory elements associated with
groups of functionally related genes in Saccharomyces cerevisiae,
Hughes, JD, Estep, PW, Tavazoie S, & GM Church, Journal of Molecular
Biology 2000 Mar 10;296(5):1205-14.

Finding DNA Regulatory Motifs within Unaligned Non-Coding Sequences
Clustered by Whole-Genome mRNA Quantitation,
Roth, FR, Hughes, JD, Estep, PE & GM Church, Nature Biotechnology
1998 Oct;16(10):939-45.

"""
from __future__ import print_function

from Bio.Application import AbstractCommandline, _Option, _Argument

import warnings
from Bio import BiopythonDeprecationWarning


class AlignAceCommandline(AbstractCommandline):
    """Create a commandline for the AlignAce program (DEPRECATED).

    Example:

    >>> from Bio.motifs.applications import AlignAceCommandline
    >>> in_file = "sequences.fasta"
    >>> alignace_cline = AlignAceCommandline(infile=in_file, gcback=0.55)
    >>> print(alignace_cline)
    AlignACE -i sequences.fasta -gcback 0.55

    You would typically run the command line with alignace_cline() or via
    the Python subprocess module, as described in the Biopython tutorial.
    """
    def __init__(self, cmd="AlignACE", **kwargs):
        warnings.warn("""The AlignACE application wrapper is deprecated and
                      is likely to be removed in a future release of Biopython,
                      since an up to date version of the AlignACE software
                      cannot be obtained anymore. If you have a copy of
                      AlignACE 4, please consider contacting the Biopython
                      developers.""", BiopythonDeprecationWarning)
        self.parameters = \
          [
            _Option(["-i", "infile"],
                    "Input Sequence file in FASTA format.",
                    checker_function=lambda x: isinstance(x, str),
                    equate=False,
                    filename=True),

            _Option(["-numcols", "numcols"],
                    "Number of columns to align",
                    equate=False,
                    checker_function=lambda x: isinstance(x, int)),

            _Option(["-expect", "expect"],
                    "number of sites expected in model",
                    equate=False,
                    checker_function=lambda x: isinstance(x, int)),

            _Option(["-gcback", "gcback"],
                    "background fractional GC content of input sequence",
                    equate=False,
                    checker_function=lambda x: isinstance(x, float)),

            _Option(["-minpass", "minpass"],
                    "minimum number of non-improved passes in phase 1",
                    equate=False,
                    checker_function=lambda x: isinstance(x, int)),

            _Option(["-seed", "seed"],
                    "set seed for random number generator (time)",
                    equate=False,
                    checker_function=lambda x: isinstance(x, int)),

            _Option(["-undersample", "undersample"],
                    "possible sites / (expect * numcols * seedings)",
                    equate=False,
                    checker_function=lambda x: isinstance(x, int)),

            _Option(["-oversample", "oversample"],
                    "1/undersample",
                    equate=False,
                    checker_function=lambda x: isinstance(x, int)),
          ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


class CompareAceCommandline(AbstractCommandline):
    """Create a commandline for the CompareAce program.

    Example:

    >>> from Bio.motifs.applications import CompareAceCommandline
    >>> m1_file = "sequences1.fasta"
    >>> m2_file = "sequences2.fasta"
    >>> compareace_cline = CompareAceCommandline(motif1=m1_file, motif2=m2_file)
    >>> print(compareace_cline)
    CompareACE sequences1.fasta sequences2.fasta

    You would typically run the command line with compareace_cline() or via
    the Python subprocess module, as described in the Biopython tutorial.
    """
    def __init__(self, cmd="CompareACE", **kwargs):
        warnings.warn("""The CompareACE application wrapper is deprecated and
                      is likely to be removed in a future release of Biopython,
                      since an up to date version of the AlignACE software
                      cannot be obtained anymore. If you have a copy of
                      AlignACE 4, please consider contacting the Biopython
                      developers.""", BiopythonDeprecationWarning)
        self.parameters = \
          [
            _Argument(["motif1"],
                        "name of file containing motif 1",
                        checker_function=lambda x: isinstance(x, str),
                        filename=True),
            _Argument(["motif2"],
                        "name of file containing motif 2",
                        checker_function=lambda x: isinstance(x, str),
                        filename=True),
          ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


def _test():
    """Run the module's doctests (PRIVATE)."""
    print("Running alignace doctests...")
    import doctest
    doctest.testmod()
    print("Done")


if __name__ == "__main__":
    _test()
