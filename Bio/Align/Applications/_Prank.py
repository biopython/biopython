# Copyright 2009 by Cymon J. Cox.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Command line wrapper for the multiple alignment program PRANK.
"""

from __future__ import print_function

__docformat__ = "restructuredtext en"  # Don't just use plain text in epydoc API pages!

from Bio.Application import _Option, _Switch, AbstractCommandline


class PrankCommandline(AbstractCommandline):
    """Command line wrapper for the multiple alignment program PRANK.

    http://www.ebi.ac.uk/goldman-srv/prank/prank/

    Example:
    --------

    To align a FASTA file (unaligned.fasta) with the output in aligned
    FASTA format with the output filename starting with "aligned" (you
    can't pick the filename explicitly), no tree output and no XML output,
    use:

    >>> from Bio.Align.Applications import PrankCommandline
    >>> prank_cline = PrankCommandline(d="unaligned.fasta",
    ...                                o="aligned", # prefix only!
    ...                                f=8, # FASTA output
    ...                                notree=True, noxml=True)
    >>> print(prank_cline)
    prank -d=unaligned.fasta -o=aligned -f=8 -noxml -notree

    You would typically run the command line with prank_cline() or via
    the Python subprocess module, as described in the Biopython tutorial.

    Citations:
    ----------

    Loytynoja, A. and Goldman, N. 2005. An algorithm for progressive
    multiple alignment of sequences with insertions. Proceedings of
    the National Academy of Sciences, 102: 10557--10562.

    Loytynoja, A. and Goldman, N. 2008. Phylogeny-aware gap placement
    prevents errors in sequence alignment and evolutionary analysis.
    Science, 320: 1632.

    Last checked against version: 081202
    """
    def __init__(self, cmd="prank", **kwargs):
        OUTPUT_FORMAT_VALUES = list(range(1, 18))
        self.parameters = [
            # ################# input/output parameters: ##################
            # -d=sequence_file
            _Option(["-d", "d"],
                    "Input filename",
                    filename=True,
                    is_required=True),
            # -t=tree_file [default: no tree, generate approximate NJ tree]
            _Option(["-t", "t"], "Input guide tree filename",
                    filename=True),
            # -tree="tree_string" [tree in newick format; in double quotes]
            _Option(["-tree", "tree"],
                    "Input guide tree as Newick string"),
            # -m=model_file [default: HKY2/WAG]
            _Option(["-m", "m"],
                    "User-defined alignment model filename. Default: "
                    "HKY2/WAG"),
            # -o=output_file [default: 'output']
            _Option(["-o", "o"],
                    "Output filenames prefix. Default: 'output'\n "
                    "Will write: output.?.fas (depending on requested "
                    "format), output.?.xml and output.?.dnd",
                    filename=True),
            # -f=output_format [default: 8]
            _Option(["-f", "f"],
                    "Output alignment format. Default: 8 FASTA\n"
                    "Option are:\n"
                    "1. IG/Stanford	8. Pearson/Fasta\n"
                    "2. GenBank/GB 	11. Phylip3.2\n"
                    "3. NBRF       	12. Phylip\n"
                    "4. EMBL       	14. PIR/CODATA\n"
                    "6. DNAStrider 	15. MSF\n"
                    "7. Fitch      	17. PAUP/NEXUS",
                    checker_function=lambda x: x in OUTPUT_FORMAT_VALUES),
            _Switch(["-noxml", "noxml"],
                    "Do not output XML files "
                    "(PRANK versions earlier than v.120626)"),
            _Switch(["-notree", "notree"],
                    "Do not output dnd tree files "
                    "(PRANK versions earlier than v.120626)"),
            _Switch(["-showxml", "showxml"],
                    "Output XML files (PRANK v.120626 and later)"),
            _Switch(["-showtree", "showtree"],
                    "Output dnd tree files (PRANK v.120626 and later)"),
            _Switch(["-shortnames", "shortnames"],
                    "Truncate names at first space"),
            _Switch(["-quiet", "quiet"],
                    "Reduce verbosity"),
            # ###################### model parameters: ######################
            # +F [force insertions to be always skipped]
            # -F [equivalent]
            _Switch(["-F", "+F", "F"],
                    "Force insertions to be always skipped: same as +F"),
            # -dots [show insertion gaps as dots]
            _Switch(["-dots", "dots"],
                    "Show insertion gaps as dots"),
            # -gaprate=# [gap opening rate; default: dna 0.025 / prot 0.0025]
            _Option(["-gaprate", "gaprate"],
                    "Gap opening rate. Default: dna 0.025 prot 0.0025",
                    checker_function=lambda x: isinstance(x, float)),
            # -gapext=# [gap extension probability; default: dna 0.5 / prot 0.5]
            _Option(["-gapext", "gapext"],
                    "Gap extension probability. Default: dna 0.5 "
                    "/ prot 0.5",
                    checker_function=lambda x: isinstance(x, float)),
            # -dnafreqs=#,#,#,# [ACGT; default: empirical]
            _Option(["-dnafreqs", "dnafreqs"],
                    "DNA frequencies - 'A,C,G,T'. eg '25,25,25,25' as a quote "
                    "surrounded string value. Default: empirical",
                    checker_function=lambda x: isinstance(x, bytes)),
            # -kappa=# [ts/tv rate ratio; default:2]
            _Option(["-kappa", "kappa"],
                    "Transition/transversion ratio. Default: 2",
                    checker_function=lambda x: isinstance(x, int)),
            # -rho=# [pur/pyr rate ratio; default:1]
            _Option(["-rho", "rho"],
                    "Purine/pyrimidine ratio. Default: 1",
                    checker_function=lambda x: isinstance(x, int)),
            # -codon [for DNA: use empirical codon model]
            # Assuming this is an input file as in -m
            _Option(["-codon", "codon"],
                    "Codon model filename. Default: empirical codon model"),
            # -termgap [penalise terminal gaps normally]
            _Switch(["-termgap", "termgap"],
                    "Penalise terminal gaps normally"),
            # ############### other parameters: ################################
            # -nopost [do not compute posterior support; default: compute]
            _Switch(["-nopost", "nopost"],
                    "Do not compute posterior support. Default: compute"),
            # -pwdist=# [expected pairwise distance for computing guidetree;
            # default: dna 0.25 / prot 0.5]
            _Option(["-pwdist", "pwdist"],
                    "Expected pairwise distance for computing guidetree. "
                    "Default: dna 0.25 / prot 0.5",
                    checker_function=lambda x: isinstance(x, float)),
            _Switch(["-once", "once"],
                    "Run only once. Default: twice if no guidetree given"),
            _Switch(["-twice", "twice"],
                    "Always run twice"),
            _Switch(["-skipins", "skipins"],
                    "Skip insertions in posterior support"),
            _Switch(["-uselogs", "uselogs"],
                    "Slower but should work for a greater number of sequences"),
            _Switch(["-writeanc", "writeanc"],
                    "Output ancestral sequences"),
            _Switch(["-printnodes", "printnodes"],
                    "Output each node; mostly for debugging"),
            # -matresize=# [matrix resizing multiplier]
            # Doesnt specify type but Float and Int work
            _Option(["-matresize", "matresize"],
                    "Matrix resizing multiplier",
                    checker_function=lambda x: isinstance(x, float) or
                                               isinstance(x, int)),
            # -matinitsize=# [matrix initial size multiplier]
            # Doesnt specify type but Float and Int work
            _Option(["-matinitsize", "matinitsize"],
                    "Matrix initial size multiplier",
                    checker_function=lambda x: isinstance(x, float) or
                                               isinstance(x, int)),
            _Switch(["-longseq", "longseq"],
                    "Save space in pairwise alignments"),
            _Switch(["-pwgenomic", "pwgenomic"],
                    "Do pairwise alignment, no guidetree"),
            # -pwgenomicdist=# [distance for pairwise alignment; default: 0.3]
            _Option(["-pwgenomicdist", "pwgenomicdist"],
                    "Distance for pairwise alignment. Default: 0.3",
                    checker_function=lambda x: isinstance(x, float)),
            # -scalebranches=# [scale branch lengths; default: dna 1 / prot 2]
            _Option(["-scalebranches", "scalebranches"],
                    "Scale branch lengths. Default: dna 1 / prot 2",
                    checker_function=lambda x: isinstance(x, int)),
            # -fixedbranches=# [use fixed branch lengths]
            # Assume looking for a float
            _Option(["-fixedbranches", "fixedbranches"],
                    "Use fixed branch lengths of input value",
                    checker_function=lambda x: isinstance(x, float)),
            # -maxbranches=# [set maximum branch length]
            # Assume looking for a float
            _Option(["-maxbranches", "maxbranches"],
                    "Use maximum branch lengths of input value",
                    checker_function=lambda x: isinstance(x, float)),
            # -realbranches [disable branch length truncation]
            _Switch(["-realbranches", "realbranches"],
                    "Disable branch length truncation"),
            _Switch(["-translate", "translate"],
                    "Translate to protein"),
            _Switch(["-mttranslate", "mttranslate"],
                    "Translate to protein using mt table"),
            # ##################### other: ####################
            _Switch(["-convert", "convert"],
                    "Convert input alignment to new format. Do "
                    "not perform alignment")
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


def _test():
    """Run the module's doctests (PRIVATE)."""
    print("Running modules doctests...")
    import doctest
    doctest.testmod()
    print("Done")

if __name__ == "__main__":
    _test()
