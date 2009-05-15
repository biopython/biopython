# Copyright 2009 by Cymon J. Cox.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Command line wrapper for the multiple alignment program PRANK.

http://www.ebi.ac.uk/goldman-srv/prank/prank/

Citations:

Loytynoja, A. and Goldman, N. 2005. An algorithm for progressive multiple
alignment of sequences with insertions. Proceedings of the National Academy
of Sciences, 102: 10557--10562.

Loytynoja, A. and Goldman, N. 2008. Phylogeny-aware gap placement prevents
errors in sequence alignment and evolutionary analysis. Science, 320: 1632.

Last checked agains version: 081202
"""
import types
from Bio.Application import _Option, _Switch, AbstractCommandline

class PrankCommandline(AbstractCommandline):
    """Command line wrapper for the multiple alignment program PRANK."""
    def __init__(self, cmd="prank", **kwargs):
        OUTPUT_FORMAT_VALUES = range(1,18)
        self.parameters = \
            [
            ################## input/output parameters: ##################
            #-d=sequence_file
            _Option(["-d", "d"], ["input", "file"],
                    None, 1, "Input filename"),
            #-t=tree_file [default: no tree, generate approximate NJ tree]
            _Option(["-t", "t"], ["input", "file"],
                    None, 0, "Input guide tree filename"),
            #-tree="tree_string" [tree in newick format; in double quotes]
            _Option(["-tree", "tree"], ["input"],
                    None, 0,
                    "Input guide tree as Newick string"),
            #-m=model_file [default: HKY2/WAG]
            _Option(["-m", "m"], ["input"],
                    None, 0,
                    "User-defined alignment model filename. Default: " + \
                    "HKY2/WAG"),
            #-o=output_file [default: 'output']
            _Option(["-o", "o"], ["output"],
                    None, 0,
                    "Output filenames prefix. Default: 'output'\n " + \
                    "Will write: output.?.fas (depending on requested " + \
                    "format), output.?.xml and output.?.dnd"),
            #-f=output_format [default: 8]
            _Option(["-f", "f"], ["input"],
                    lambda x: x in OUTPUT_FORMAT_VALUES, 0,
                    "Output alignment format. Default: 8 FASTA\n" + \
                    "Option are:\n" + \
                    "1. IG/Stanford	8. Pearson/Fasta\n" + \
                    "2. GenBank/GB 	11. Phylip3.2\n" + \
                    "3. NBRF       	12. Phylip\n" + \
                    "4. EMBL       	14. PIR/CODATA\n" + \
                    "6. DNAStrider 	15. MSF\n" + \
                    "7. Fitch      	17. PAUP/NEXUS"),
            _Switch(["-noxml", "noxml"], ["input"],
                    "Do not output XML files"),
            _Switch(["-notree", "notree"], ["input"],
                    "Do not output dnd tree files"),
            _Switch(["-shortnames", "shortnames"], ["input"],
                    "Truncate names at first space"),
            _Switch(["-quiet", "quiet"], ["input"],
                    "Reduce verbosity"),
            ####################### model parameters: ######################
            #+F [force insertions to be always skipped]
            #-F [equivalent]
            _Switch(["-F", "+F", "F"], ["input"],
                    "Force insertions to be always skipped: same as +F"),
            #-dots [show insertion gaps as dots]
            _Switch(["-dots", "dots"], ["input"],
                    "Show insertion gaps as dots"),
            #-gaprate=# [gap opening rate; default: dna 0.025 / prot 0.0025]
            _Option(["-gaprate", "gaprate"], ["input"],
                    lambda x: isinstance(x, types.FloatType), 
                    0,
                    "Gap opening rate. Default: dna 0.025 prot 0.0025"),
            #-gapext=# [gap extension probability; default: dna 0.5 / prot 0.5]
            _Option(["-gapext", "gapext"], ["input"],
                    lambda x: isinstance(x, types.FloatType), 
                    0,
                    "Gap extension probability. Default: dna 0.5 " + \
                    "/ prot 0.5"),
            #-dnafreqs=#,#,#,# [ACGT; default: empirical]
            _Option(["-dnafreqs", "dnafreqs"], ["input"],
                    lambda x: isinstance(x, types.StringType), 
                    0,
                    "DNA frequencies - 'A,C,G,T'. eg '25,25,25,25' as a quote " + \
                    "surrounded string value. Default: empirical"),
            #-kappa=# [ts/tv rate ratio; default:2]
            _Option(["-kappa", "kappa"], ["input"],
                    lambda x: isinstance(x, types.IntType), 
                    0,
                    "Transition/transversion ratio. Default: 2"),
            #-rho=# [pur/pyr rate ratio; default:1]
            _Option(["-rho", "rho"], ["input"],
                    lambda x: isinstance(x, types.IntType), 
                    0,
                    "Purine/pyrimidine ratio. Default: 1"),
            #-codon [for DNA: use empirical codon model]
            #Assuming this is an input file as in -m
            _Option(["-codon", "codon"], ["input"],
                    None, 
                    0,
                    "Codon model filename. Default: empirical codon model"),
            #-termgap [penalise terminal gaps normally]
            _Option(["-termgap", "termgap"], ["input"],
                    lambda x: 0, #Does not take a value 
                    0,
                    "Penalise terminal gaps normally",
                    0),
            ################ other parameters: ################################
            #-nopost [do not compute posterior support; default: compute]
            _Option(["-nopost", "nopost"], ["input"],
                    lambda x: 0, #Does not take a value 
                    0,
                    "Do not compute posterior support. Default: compute",
                    0),
            #-pwdist=# [expected pairwise distance for computing guidetree;
            #default: dna 0.25 / prot 0.5]
            _Option(["-pwdist", "pwdist"], ["input"],
                    lambda x: isinstance(x, types.FloatType),
                    0,
                    "Expected pairwise distance for computing guidetree. " + \
                    "Default: dna 0.25 / prot 0.5"),
            _Switch(["-once", "once"], ["input"],
                    "Run only once. Default: twice if no guidetree given"),
            _Switch(["-twice", "twice"], ["input"],
                    "Always run twice"),
            _Switch(["-skipins", "skipins"], ["input"],
                    "Skip insertions in posterior support"),
            _Switch(["-uselogs", "uselogs"], ["input"],
                    "Slower but should work for a greater number of sequences"),
            _Switch(["-writeanc", "writeanc"], ["input"],
                    "Output ancestral sequences"),
            _Switch(["-printnodes", "printnodes"], ["input"],
                    "Output each node; mostly for debugging"),
            #-matresize=# [matrix resizing multiplier]
            # Doesnt specify type but Float and Int work
            _Option(["-matresize", "matresize"], ["input"],
                    lambda x: isinstance(x, types.FloatType) or isinstance(x,
                              types.IntType),
                    0,
                    "Matrix resizing multiplier"),
            #-matinitsize=# [matrix initial size multiplier]
            # Doesnt specify type but Float and Int work
            _Option(["-matinitsize", "matinitsize"], ["input"],
                    lambda x: isinstance(x, types.FloatType) or isinstance(x,
                              types.IntType),
                    0,
                    "Matrix initial size multiplier"),
            _Switch(["-longseq", "longseq"], ["input"],
                    "Save space in pairwise alignments"),
            _Switch(["-pwgenomic", "pwgenomic"], ["input"],
                    "Do pairwise alignment, no guidetree"),
            #-pwgenomicdist=# [distance for pairwise alignment; default: 0.3]
            _Option(["-pwgenomicdist", "pwgenomicdist"], ["input"],
                    lambda x: isinstance(x, types.FloatType),
                    0,
                    "Distance for pairwise alignment. Default: 0.3"),
            #-scalebranches=# [scale branch lengths; default: dna 1 / prot 2]
            _Option(["-scalebranches", "scalebranches"], ["input"],
                    lambda x: isinstance(x, types.IntType),
                    0,
                    "Scale branch lengths. Default: dna 1 / prot 2"),
            #-fixedbranches=# [use fixed branch lengths]
            #Assume looking for a float
            _Option(["-fixedbranches", "fixedbranches"], ["input"],
                    lambda x: isinstance(x, types.FloatType),
                    0,
                    "Use fixed branch lengths of input value"),
            #-maxbranches=# [set maximum branch length]
            #Assume looking for a float
            _Option(["-maxbranches", "maxbranches"], ["input"],
                    lambda x: isinstance(x, types.FloatType),
                    0,
                    "Use maximum branch lengths of input value"),
            #-realbranches [disable branch length truncation]
            _Switch(["-realbranches", "realbranches"], ["input"],
                    "Disable branch length truncation"),
            _Switch(["-translate", "translate"], ["input"],
                    "Translate to protein"),
            _Switch(["-mttranslate", "mttranslate"], ["input"],
                    "Translate to protein using mt table"),
            ###################### other: ####################
            _Switch(["-convert", "convert"], ["input"],
                    "Convert input alignment to new format. Do " + \
                    "not perform alignment")
            ]
        AbstractCommandline.__init__(self, cmd, **kwargs)
