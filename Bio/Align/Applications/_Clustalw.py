# Copyright 2009 by Cymon J. Cox.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Command line wrapper for the multiple alignment program Clustal W.
"""

__docformat__ = "epytext en" #Don't just use plain text in epydoc API pages!

import os
from Bio.Application import _Option, _Switch, AbstractCommandline

class ClustalwCommandline(AbstractCommandline):
    """Command line wrapper for clustalw (version one or two).

    http://www.clustal.org/

    Example:

    >>> from Bio.Align.Applications import ClustalwCommandline
    >>> in_file = "unaligned.fasta"
    >>> clustalw_cline = ClustalwCommandline("clustalw2", infile=in_file)
    >>> print clustalw_cline
    clustalw2 -infile=unaligned.fasta

    You would typically run the command line with clustalw_cline() or via
    the Python subprocess module, as described in the Biopython tutorial.

    Citation:

    Larkin MA, Blackshields G, Brown NP, Chenna R, McGettigan PA,
    McWilliam H, Valentin F, Wallace IM, Wilm A, Lopez R, Thompson JD,
    Gibson TJ, Higgins DG. (2007). Clustal W and Clustal X version 2.0.
    Bioinformatics, 23, 2947-2948. 

    Last checked against versions: 1.83 and 2.0.10
    """
    #TODO - Should we default to cmd="clustalw2" now?
    def __init__(self, cmd="clustalw", **kwargs):
        self.parameters = \
            [
            _Option(["-infile", "-INFILE", "INFILE", "infile"],
                    "Input sequences.",
                    filename=True),
            _Option(["-profile1", "-PROFILE1", "PROFILE1", "profile1"],
                    "Profiles (old alignment).",
                    filename=True),
            _Option(["-profile2", "-PROFILE2", "PROFILE2", "profile2"],
                    "Profiles (old alignment).",
                    filename=True),
            ################## VERBS (do things) #############################
            _Switch(["-options", "-OPTIONS", "OPTIONS", "options"],
                    "List the command line parameters"),
            _Switch(["-help", "-HELP", "HELP", "help"],
                    "Outline the command line params."),
            _Switch(["-check", "-CHECK", "CHECK", "check"],
                    "Outline the command line params."),
            _Switch(["-fullhelp", "-FULLHELP", "FULLHELP", "fullhelp"],
                    "Output full help content."),
            _Switch(["-align", "-ALIGN", "ALIGN", "align"],
                    "Do full multiple alignment."),
            _Switch(["-tree", "-TREE", "TREE", "tree"],
                    "Calculate NJ tree."),
            _Option(["-bootstrap", "-BOOTSTRAP", "BOOTSTRAP", "bootstrap"],
                    "Bootstrap a NJ tree (n= number of bootstraps; def. = 1000).",
                    checker_function=lambda x: isinstance(x, int)),
            _Switch(["-convert", "-CONVERT", "CONVERT", "convert"],
                    "Output the input sequences in a different file format."),
            ##################### PARAMETERS (set things) #########################
            # ***General settings:****
            # Makes no sense in biopython
            #_Option(["-interactive", "-INTERACTIVE", "INTERACTIVE", "interactive"],
            #        [],
            #        lambda x: 0, #Does not take value
            #        False,
            #        "read command line, then enter normal interactive menus",
            #        False),
            _Switch(["-quicktree", "-QUICKTREE", "QUICKTREE", "quicktree"],
                    "Use FAST algorithm for the alignment guide tree"),
            _Option(["-type", "-TYPE", "TYPE", "type"],
                    "PROTEIN or DNA sequences",
                    checker_function=lambda x: x in ["PROTEIN", "DNA",
                                                     "protein", "dna"]),
            _Switch(["-negative", "-NEGATIVE", "NEGATIVE", "negative"],
                    "Protein alignment with negative values in matrix"),
            _Option(["-outfile", "-OUTFILE", "OUTFILE", "outfile"],
                    "Output sequence alignment file name",
                    filename=True),
            _Option(["-output", "-OUTPUT", "OUTPUT", "output"],
                    "Output format: GCG, GDE, PHYLIP, PIR or NEXUS",
                    checker_function=lambda x: x in ["GCG", "GDE", "PHYLIP",
                                                     "PIR", "NEXUS",
                                                     "gcg", "gde", "phylip",
                                                     "pir", "nexus"]),
            _Option(["-outorder", "-OUTORDER", "OUTORDER", "outorder"],
                    "Output taxon order: INPUT or ALIGNED",
                    checker_function=lambda x: x in ["INPUT", "input",
                                                     "ALIGNED", "aligned"]),
            _Option(["-case", "-CASE", "CASE", "case"],
                    "LOWER or UPPER (for GDE output only)",
                    checker_function=lambda x: x in ["UPPER", "upper",
                                                     "LOWER", "lower"]),
            _Option(["-seqnos", "-SEQNOS", "SEQNOS", "seqnos"],
                    "OFF or ON (for Clustal output only)",
                    checker_function=lambda x: x in ["ON", "on",
                                                     "OFF", "off"]),
            _Option(["-seqno_range", "-SEQNO_RANGE", "SEQNO_RANGE", "seqno_range"],
                    "OFF or ON (NEW- for all output formats)",
                    checker_function=lambda x: x in ["ON", "on",
                                                     "OFF", "off"]),
            _Option(["-range", "-RANGE", "RANGE", "range"],
                    "Sequence range to write starting m to m+n. "
                    "Input as string eg. '24,200'"),
            _Option(["-maxseqlen", "-MAXSEQLEN", "MAXSEQLEN", "maxseqlen"],
                    "Maximum allowed input sequence length",
                    checker_function=lambda x: isinstance(x, int)),
            _Switch(["-quiet", "-QUIET", "QUIET", "quiet"],
                    "Reduce console output to minimum"),
            _Switch(["-stats", "-STATS", "STATS", "stats"],
                    "Log some alignents statistics to file"),
            # ***Fast Pairwise Alignments:***
            _Option(["-ktuple", "-KTUPLE", "KTUPLE", "ktuple"],
                    "Word size",
                    checker_function=lambda x: isinstance(x, int) or \
                                               isinstance(x, float)),
            _Option(["-topdiags", "-TOPDIAGS", "TOPDIAGS", "topdiags"],
                    "Number of best diags.",
                    checker_function=lambda x: isinstance(x, int) or \
                                               isinstance(x, float)),
            _Option(["-window", "-WINDOW", "WINDOW", "window"],
                    "Window around best diags.",
                    checker_function=lambda x: isinstance(x, int) or \
                                               isinstance(x, float)),
            _Option(["-pairgap", "-PAIRGAP", "PAIRGAP", "pairgap"],
                    "Gap penalty",
                    checker_function=lambda x: isinstance(x, int) or \
                                               isinstance(x, float)),
            _Option(["-score", "-SCORE", "SCORE", "score"],
                    "Either: PERCENT or ABSOLUTE",
                    checker_function=lambda x: x in ["percent", "PERCENT",
                                                     "absolute","ABSOLUTE"]),
            # ***Slow Pairwise Alignments:***
            _Option(["-pwmatrix", "-PWMATRIX", "PWMATRIX", "pwmatrix"],
                    "Protein weight matrix=BLOSUM, PAM, GONNET, ID or filename",
                    checker_function=lambda x: x in ["BLOSUM", "PAM",
                                                     "GONNET", "ID",
                                                     "blosum", "pam",
                                                     "gonnet", "id"] or \
                                                os.path.exists(x),
                    filename=True),
            _Option(["-pwdnamatrix", "-PWDNAMATRIX", "PWDNAMATRIX", "pwdnamatrix"],
                    "DNA weight matrix=IUB, CLUSTALW or filename",
                    checker_function=lambda x: x in ["IUB", "CLUSTALW",
                                                     "iub", "clustalw"] or \
                                               os.path.exists(x),
                    filename=True),
            _Option(["-pwgapopen", "-PWGAPOPEN", "PWGAPOPEN", "pwgapopen"],
                    "Gap opening penalty",
                    checker_function=lambda x: isinstance(x, int) or \
                                               isinstance(x, float)),
            _Option(["-pwgapext", "-PWGAPEXT", "PWGAPEXT", "pwgapext"],
                    "Gap opening penalty",
                    checker_function=lambda x: isinstance(x, int) or \
                                               isinstance(x, float)),
            # ***Multiple Alignments:***
            _Option(["-newtree", "-NEWTREE", "NEWTREE", "newtree"],
                    "Output file name for newly created guide tree",
                    filename=True),
            _Option(["-usetree", "-USETREE", "USETREE", "usetree"],
                    "File name of guide tree",
                    checker_function=lambda x: os.path.exists,
                    filename=True),
            _Option(["-matrix", "-MATRIX", "MATRIX", "matrix"],
                    "Protein weight matrix=BLOSUM, PAM, GONNET, ID or filename",
                    checker_function=lambda x: x in ["BLOSUM", "PAM",
                                                     "GONNET", "ID",
                                                     "blosum", "pam",
                                                     "gonnet", "id"] or \
                                               os.path.exists(x),
                    filename=True),
            _Option(["-dnamatrix", "-DNAMATRIX", "DNAMATRIX", "dnamatrix"],
                    "DNA weight matrix=IUB, CLUSTALW or filename",
                    checker_function=lambda x: x in ["IUB", "CLUSTALW",
                                                     "iub", "clustalw"] or \
                                               os.path.exists(x),
                    filename=True),
            _Option(["-gapopen", "-GAPOPEN", "GAPOPEN", "gapopen"],
                    "Gap opening penalty",
                    checker_function=lambda x: isinstance(x, int) or \
                                               isinstance(x, float)),
            _Option(["-gapext", "-GAPEXT", "GAPEXT", "gapext"],
                    "Gap extension penalty",
                    checker_function=lambda x: isinstance(x, int) or \
                                               isinstance(x, float)),
            _Switch(["-endgaps", "-ENDGAPS", "ENDGAPS", "endgaps"],
                    "No end gap separation pen."),
            _Option(["-gapdist", "-GAPDIST", "GAPDIST", "gapdist"],
                    "Gap separation pen. range",
                    checker_function=lambda x: isinstance(x, int) or \
                                               isinstance(x, float)),
            _Switch(["-nopgap", "-NOPGAP", "NOPGAP", "nopgap"],
                    "Residue-specific gaps off"),
            _Switch(["-nohgap", "-NOHGAP", "NOHGAP", "nohgap"],
                    "Hydrophilic gaps off"),
            _Switch(["-hgapresidues", "-HGAPRESIDUES", "HGAPRESIDUES", "hgapresidues"],
                    "List hydrophilic res."),
            _Option(["-maxdiv", "-MAXDIV", "MAXDIV", "maxdiv"],
                    "% ident. for delay",
                    checker_function=lambda x: isinstance(x, int) or \
                                               isinstance(x, float)),
            _Option(["-transweight", "-TRANSWEIGHT", "TRANSWEIGHT", "transweight"],
                    "Transitions weighting",
                    checker_function=lambda x: isinstance(x, int) or \
                                               isinstance(x, float)),
            _Option(["-iteration", "-ITERATION", "ITERATION", "iteration"],
                    "NONE or TREE or ALIGNMENT",
                    checker_function=lambda x: x in ["NONE", "TREE",
                                                     "ALIGNMENT",
                                                     "none", "tree",
                                                     "alignment"]),
            _Option(["-numiter", "-NUMITER", "NUMITER", "numiter"],
                    "maximum number of iterations to perform",
                    checker_function=lambda x: isinstance(x, int)),
            _Switch(["-noweights", "-NOWEIGHTS", "NOWEIGHTS", "noweights"],
                    "Disable sequence weighting"),
            # ***Profile Alignments:***
            _Switch(["-profile", "-PROFILE", "PROFILE", "profile"],
                    "Merge two alignments by profile alignment"),
            _Option(["-newtree1", "-NEWTREE1", "NEWTREE1", "newtree1"],
                    "Output file name for new guide tree of profile1",
                    filename=True),
            _Option(["-newtree2", "-NEWTREE2", "NEWTREE2", "newtree2"],
                    "Output file for new guide tree of profile2",
                    filename=True),
            _Option(["-usetree1", "-USETREE1", "USETREE1", "usetree1"],
                    "File name of guide tree for profile1",
                    checker_function=lambda x: os.path.exists,
                    filename=True),
            _Option(["-usetree2", "-USETREE2", "USETREE2", "usetree2"],
                    "File name of guide tree for profile2",
                    checker_function=lambda x: os.path.exists,
                    filename=True),
            # ***Sequence to Profile Alignments:***
            _Switch(["-sequences", "-SEQUENCES", "SEQUENCES", "sequences"],
                    "Sequentially add profile2 sequences to profile1 alignment"),
            _Switch(["-nosecstr1", "-NOSECSTR1", "NOSECSTR1", "nosecstr1"],
                    "Do not use secondary structure-gap penalty mask for profile 1"),
            _Switch(["-nosecstr2", "-NOSECSTR2", "NOSECSTR2", "nosecstr2"],
                    "Do not use secondary structure-gap penalty mask for profile 2"),
            # ***Structure Alignments:***
            _Option(["-secstrout", "-SECSTROUT", "SECSTROUT", "secstrout"],
                    "STRUCTURE or MASK or BOTH or NONE output in alignment file",
                    checker_function=lambda x: x in ["STRUCTURE", "MASK",
                                                     "BOTH", "NONE",
                                                     "structure", "mask",
                                                     "both", "none"]),
            _Option(["-helixgap", "-HELIXGAP", "HELIXGAP", "helixgap"],
                    "Gap penalty for helix core residues",
                    checker_function=lambda x: isinstance(x, int) or \
                                               isinstance(x, float)),
            _Option(["-strandgap", "-STRANDGAP", "STRANDGAP", "strandgap"],
                    "gap penalty for strand core residues",
                    checker_function=lambda x: isinstance(x, int) or \
                                               isinstance(x, float)),
            _Option(["-loopgap", "-LOOPGAP", "LOOPGAP", "loopgap"],
                    "Gap penalty for loop regions",
                    checker_function=lambda x: isinstance(x, int) or \
                                               isinstance(x, float)),
            _Option(["-terminalgap", "-TERMINALGAP", "TERMINALGAP", "terminalgap"],
                    "Gap penalty for structure termini",
                    checker_function=lambda x: isinstance(x, int) or \
                                               isinstance(x, float)),
            _Option(["-helixendin", "-HELIXENDIN", "HELIXENDIN", "helixendin"],
                    "Number of residues inside helix to be treated as terminal",
                    checker_function=lambda x: isinstance(x, int)),
            _Option(["-helixendout", "-HELIXENDOUT", "HELIXENDOUT", "helixendout"],
                    "Number of residues outside helix to be treated as terminal",
                    checker_function=lambda x: isinstance(x, int)),
            _Option(["-strandendin", "-STRANDENDIN", "STRANDENDIN", "strandendin"],
                    "Number of residues inside strand to be treated as terminal",
                    checker_function=lambda x: isinstance(x, int)),
            _Option(["-strandendout", "-STRANDENDOUT", "STRANDENDOUT", "strandendout"],
                    "number of residues outside strand to be treated as terminal",
                    checker_function=lambda x: isinstance(x, int)),
            # ***Trees:***
            _Option(["-outputtree", "-OUTPUTTREE", "OUTPUTTREE", "outputtree"],
                    "nj OR phylip OR dist OR nexus",
                    checker_function=lambda x: x in ["NJ", "PHYLIP",
                                                     "DIST", "NEXUS",
                                                     "nj", "phylip",
                                                     "dist", "nexus"]),
            _Option(["-seed", "-SEED", "SEED", "seed"],
                    "Seed number for bootstraps.",
                    checker_function=lambda x: isinstance(x, int)),
            _Switch(["-kimura", "-KIMURA", "KIMURA", "kimura"],
                    "Use Kimura's correction."),
            _Switch(["-tossgaps", "-TOSSGAPS", "TOSSGAPS", "tossgaps"],
                    "Ignore positions with gaps."),
            _Option(["-bootlabels", "-BOOTLABELS", "BOOTLABELS", "bootlabels"],
                    "Node OR branch position of bootstrap values in tree display",
                    checker_function=lambda x: x in ["NODE", "BRANCH",
                                                     "node", "branch"]),
            _Option(["-clustering", "-CLUSTERING", "CLUSTERING", "clustering"],
                    "NJ or UPGMA",
                    checker_function=lambda x: x in ["NJ", "UPGMA", "nj", "upgma"])
            ]
        AbstractCommandline.__init__(self, cmd, **kwargs)

def _test():
    """Run the module's doctests (PRIVATE)."""
    print "Runing ClustalW doctests..."
    import doctest
    doctest.testmod()
    print "Done"

if __name__ == "__main__":
    _test()
