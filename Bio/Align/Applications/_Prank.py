# Copyright 2009 by Cymon J. Cox.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.Application command line for the multiple alignment program PRANK.

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
from Bio import Application
from Bio.Application import _Option
from Bio.Application import _Argument

class PrankCommandline(Application.AbstractCommandline):
    
    def __init__(self, cmd = "prank"):

        OUTPUT_FORMAT_VALUES = range(1,18)

        Application.AbstractCommandline.__init__(self)
        self.program_name = cmd
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

            #-noxml [do not output xml-files]
            _Option(["-noxml", "noxml"], ["input"],
                    lambda x: 0, # Does not take a value
                    0,
                    "Do not output XML files",
                    0), #No equate

            #-notree [do not output dnd-files]
            _Option(["-notree", "notree"], ["input"],
                    lambda x: 0, # Does not take a value
                    0,
                    "Do not output dnd tree files",
                    0), #No equate

            #-shortnames [truncate names at first space]
            _Option(["-shortnames", "shortnames"], ["input"],
                    lambda x: 0, # Does not take a value
                    0,
                    "Truncate names at first space",
                    0),

            #-quiet
            _Option(["-quiet", "quiet"], ["input"],
                    lambda x: 0, # Does not take a value
                    0,
                    "Reduce verbosity",
                    0),

            ####################### model parameters: ######################
            #+F [force insertions to be always skipped]
            #-F [equivalent]
            _Option(["-F", "+F", "F"], ["input"],
                    lambda x: 0, # Does not take a value
                    0,
                    "Force insertions to be always skipped: same as +F",
                    0),

            #-dots [show insertion gaps as dots]
            _Option(["-dots", "dots"], ["input"],
                    lambda x: 0, # Does not take a value
                    0,
                    "Show insertion gaps as dots",
                    0),

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

            #-once [run only once; default: twice if no guidetree given]
            _Option(["-once", "once"], ["input"],
                    lambda x: 0, #Does not take a value,
                    0,
                    "Run only once. Default: twice if no guidetree given",
                    0),

            #-twice [run always twice]
            _Option(["-twice", "twice"], ["input"],
                    lambda x: 0, #Does not take a value,
                    0,
                    "Always run twice",
                    0),

            #-skipins [skip insertions in posterior support]
            _Option(["-skipins", "skipins"], ["input"],
                    lambda x: 0, #Does not take a value,
                    0,
                    "Skip insertions in posterior support",
                    0),

            #-uselogs [slower but should work for a greater number of sequences]
            _Option(["-uselogs", "uselogs"], ["input"],
                    lambda x: 0, #Does not take a value,
                    0,
                    "Slower but should work for a greater number of sequences",
                    0),

            #-writeanc [output ancestral sequences]
            _Option(["-writeanc", "writeanc"], ["input"],
                    lambda x: 0, #Does not take a value,
                    0,
                    "Output ancestral sequences",
                    0),

            #-printnodes [output each node; mostly for debugging]
            _Option(["-printnodes", "printnodes"], ["input"],
                    lambda x: 0, #Does not take a value,
                    0,
                    "Output each node; mostly for debugging",
                    0),

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

            #-longseq [save space in pairwise alignments]
            _Option(["-longseq", "longseq"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Save space in pairwise alignments",
                    0),

            #-pwgenomic [do pairwise alignment, no guidetree]
            _Option(["-pwgenomic", "pwgenomic"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Do pairwise alignment, no guidetree",
                    0),

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
            _Option(["-realbranches", "realbranches"], ["input"],
                    lambda x: isinstance(x, types.FloatType),
                    0,
                    "Disable branch length truncation",
                    0),

            #-translate [translate to protein]
            _Option(["-translate", "translate"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Translate to protein",
                    0),

            #-mttranslate [translate to protein using mt table]
            _Option(["-mttranslate", "mttranslate"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Translate to protein using mt table",
                    0),

            ###################### other: ####################
            #-convert [no alignment, just convert to another format]
            _Option(["-convert", "convert"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Convert input alignment to new format. Do " + \
                    "not perform alignment",
                    0)
            ]

