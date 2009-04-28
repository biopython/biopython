# Copyright 2009 by Cymon J. Cox.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
Bio.Application command line for the multiple alignment programme MAFFT

Supports the commands: "mafft", "linsi", "ginsi", "einsi", "fftnsi",
                       "fftns", "nwns", "nwnsi", and "mafft-profile".

http://align.bmr.kyushu-u.ac.jp/mafft/software/

Citations:
    
Katoh, Toh (BMC Bioinformatics 9:212, 2008) Improved accuracy of
multiple ncRNA alignment by incorporating structural information into a
MAFFT-based framework (describes RNA structural alignment methods)

Katoh, Toh (Briefings in Bioinformatics 9:286-298, 2008) Recent developments in
the MAFFT multiple sequence alignment program (outlines version 6)

Katoh, Toh (Bioinformatics 23:372-374, 2007)  Errata PartTree: an algorithm to
build an approximate tree from a large number of unaligned sequences (describes
the PartTree algorithm)

Katoh, Kuma, Toh, Miyata (Nucleic Acids Res. 33:511-518, 2005) MAFFT version 5:
improvement in accuracy of multiple sequence alignment (describes [ancestral
versions of] the G-INS-i, L-INS-i and E-INS-i strategies) Katoh, Misawa, Kuma,
Miyata (Nucleic Acids Res. 30:3059-3066, 2002)

Last checked against version: 6.626b (2009/03/16)
"""

import os
import types
from Bio import Application
from Bio.Application import _Option
from Bio.Application import _Argument

class MafftCommandline(Application.AbstractCommandline):

    def __init__(self, cmd = "mafft"):

        BLOSUM_MATRICES = ["30","45","62","80"]

        Application.AbstractCommandline.__init__(self)
        self.program_name = cmd
        self.parameters = \
            [
            #**** Algorithm ****

            #Automatically selects an appropriate strategy from L-INS-i, FFT-NS-
            #i and FFT-NS-2, according to data size. Default: off (always FFT-NS-2)
            _Option(["--auto", "auto"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Automatically select strategy. Default off.",
                    0), #No equate

            #Distance is calculated based on the number of shared 6mers. Default: on
            _Option(["--6merpair", "6merpair"], ["input"],
                     lambda x: 0, #Does not take a value
                     0,
                     "Distance is calculated based on the number of shared " + \
                     "6mers. Default: on",
                     0), #No equate

            #All pairwise alignments are computed with the Needleman-Wunsch
            #algorithm. More accurate but slower than --6merpair. Suitable for a
            #set of globally alignable sequences. Applicable to up to ~200
            #sequences. A combination with --maxiterate 1000 is recommended (G-
            #INS-i). Default: off (6mer distance is used)
            _Option(["--globalpair", "globalpair"], ["input"],
                     lambda x: 0, #Does not take a value
                     0,
                     "All pairwise alignments are computed with the " + \
                     "Needleman-Wunsch algorithm. Default: off",
                     0),

            #All pairwise alignments are computed with the Smith-Waterman
            #algorithm. More accurate but slower than --6merpair. Suitable for a
            #set of locally alignable sequences. Applicable to up to ~200
            #sequences. A combination with --maxiterate 1000 is recommended (L-
            #INS-i). Default: off (6mer distance is used)
            _Option(["--localpair", "localpair"], ["input"],
                     lambda x: 0, #Does not take a value,
                     0,
                     "All pairwise alignments are computed with the " + \
                     "Smith-Waterman algorithm. Default: off",
                     0),

            #All pairwise alignments are computed with a local algorithm with
            #the generalized affine gap cost (Altschul 1998). More accurate but
            #slower than --6merpair. Suitable when large internal gaps are
            #expected. Applicable to up to ~200 sequences. A combination with --
            #maxiterate 1000 is recommended (E-INS-i). Default: off (6mer
            #distance is used)
            _Option(["--genafpair", "genafpair"], ["input"],
                     lambda x: 0, #Does not take a value
                     0,
                     "All pairwise alignments are computed with a local " + \
                     "algorithm with the generalized affine gap cost " + \
                     "(Altschul 1998). Default: off",
                     0),

            #All pairwise alignments are computed with FASTA (Pearson and Lipman
            #1988). FASTA is required. Default: off (6mer distance is used)
            _Option(["--fastapair", "fastapair"], ["input"],
                     lambda x: 0, #Does not take a value
                     0,
                     "All pairwise alignments are computed with FASTA " + \
                     "(Pearson and Lipman 1988). Default: off",
                     0),

            #Weighting factor for the consistency term calculated from pairwise
            #alignments. Valid when either of --blobalpair, --localpair, --
            #genafpair, --fastapair or --blastpair is selected. Default: 2.7
            _Option(["--weighti", "weighti"], ["input"],
                     lambda x: isinstance(x, types.FloatType), 0,
                     "Weighting factor for the consistency term calculated " + \
                     "from pairwise alignments. Default: 2.7",
                     0),

            #Guide tree is built number times in the progressive stage. Valid
            #with 6mer distance. Default: 2
            _Option(["--retree", "retree"], ["input"],
                     lambda x: isinstance(x, types.IntType), 0,
                     "Guide tree is built number times in the progressive " + \
                     "stage. Valid with 6mer distance. Default: 2",
                     0),

            #Number cycles of iterative refinement are performed. Default: 0
            _Option(["--maxiterate", "maxiterate"], ["input"],
                     lambda x: isinstance(x, types.IntType), 0,
                     "Number cycles of iterative refinement are performed. " + \
                     "Default: 0",
                     0),

            #Use FFT approximation in group-to-group alignment. Default: on
            _Option(["--fft", "fft"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Use FFT approximation in group-to-group alignment. " + \
                    "Default: on",
                    0),

            #Do not use FFT approximation in group-to-group alignment. Default:
            #off
            _Option(["--nofft", "nofft"], ["input"],
                     lambda x: 0, #Does not take a value
                     0,
                     "Do not use FFT approximation in group-to-group " + \
                     "alignment. Default: off",
                     0),
            
            #Alignment score is not checked in the iterative refinement stage.
            #Default: off (score is checked)
            _Option(["--noscore", "noscore"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Alignment score is not checked in the iterative " + \
                    "refinement stage. Default: off (score is checked)",
                    0),

            #Use the Myers-Miller (1988) algorithm. Default: automatically
            #turned on when the alignment length exceeds 10,000 (aa/nt).
            _Option(["--memsave", "memsave"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Use the Myers-Miller (1988) algorithm. Default: " + \
                    "automatically turned on when the alignment length " + \
                    "exceeds 10,000 (aa/nt).",
                    0),

            #Use a fast tree-building method (PartTree, Katoh and Toh 2007) with
            #the 6mer distance. Recommended for a large number (> ~10,000) of
            #sequences are input. Default: off
            _Option(["--parttree", "parttree"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Use a fast tree-building method with the 6mer " + \
                    "distance. Default: off",
                    0),

            #The PartTree algorithm is used with distances based on DP. Slightly
            #more accurate and slower than --parttree. Recommended for a large
            #number (> ~10,000) of sequences are input. Default: off
            _Option(["--dpparttree", "dpparttree"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "The PartTree algorithm is used with distances " + \
                    "based on DP. Default: off",
                    0),
              
            #The PartTree algorithm is used with distances based on FASTA.
            #Slightly more accurate and slower than --parttree. Recommended for
            #a large number (> ~10,000) of sequences are input. FASTA is
            #required. Default: off
            _Option(["--fastaparttree", "fastaparttree"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "The PartTree algorithm is used with distances based " + \
                    "on FASTA. Default: off",
                    0),

            #The number of partitions in the PartTree algorithm. Default: 50
            _Option(["--partsize", "partsize"], ["input"],
                    lambda x: isinstance(x, types.IntType), 0,
                    "The number of partitions in the PartTree algorithm. " + \
                    "Default: 50",
                    0),

            #Do not make alignment larger than number sequences. Valid only with
            #the --*parttree options. Default: the number of input sequences
            _Option(["--groupsize", "groupsize"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Do not make alignment larger than number sequences. " + \
                    "Default: the number of input sequences",
                    0),

            #**** Parameter ****
            #Gap opening penalty at group-to-group alignment. Default: 1.53
            _Option(["--op", "op"], ["input"],
                    lambda x: isinstance(x, types.FloatType), 0,
                    "Gap opening penalty at group-to-group alignment. " + \
                    "Default: 1.53",
                    0),

            #Offset value, which works like gap extension penalty, for group-to-
            #group alignment. Deafult: 0.123
            _Option(["--ep", "ep"], ["input"],
                    lambda x: isinstance(x, types.FloatType), 0,
                    "Offset value, which works like gap extension penalty, " + \
                    "for group-to- group alignment. Default: 0.123",
                    0),

            #Gap opening penalty at local pairwise alignment. Valid when the --
            #localpair or --genafpair option is selected. Default: -2.00
            _Option(["--lop", "lop"], ["input"],
                    lambda x: isinstance(x, types.FloatType), 0,
                    "Gap opening penalty at local pairwise alignment. " + \
                    "Default: 0.123",
                    0),

            #Offset value at local pairwise alignment. Valid when the --
            #localpair or --genafpair option is selected. Default: 0.1
            _Option(["--lep", "lep"], ["input"],
                    lambda x: isinstance(x, types.FloatType), 0,
                    "Offset value at local pairwise alignment. " + \
                    "Default: 0.1",
                    0),

            #Gap extension penalty at local pairwise alignment. Valid when the -
            #-localpair or --genafpair option is selected. Default: -0.1
            _Option(["--lexp", "lexp"], ["input"],
                    lambda x: isinstance(x, types.FloatType), 0,
                    "Gap extension penalty at local pairwise alignment. " + \
                    "Default: -0.1",
                    0),

            #Gap opening penalty to skip the alignment. Valid when the --
            #genafpair option is selected. Default: -6.00
            _Option(["--LOP", "LOP"], ["input"],
                    lambda x: isinstance(x, types.FloatType), 0,
                    "Gap opening penalty to skip the alignment. " + \
                    "Default: -6.00",
                    0),

            #Gap extension penalty to skip the alignment. Valid when the --
            #genafpair option is selected. Default: 0.00
            _Option(["--LEXP", "LEXP"], ["input"],
                    lambda x: isinstance(x, types.FloatType),
                    0,
                    "Gap extension penalty to skip the alignment. " + \
                    "Default: 0.00",
                    0),


            #BLOSUM number matrix (Henikoff and Henikoff 1992) is used.
            #number=30, 45, 62 or 80. Default: 62
            _Option(["--bl", "bl"], ["input"],
                    lambda x: x in BLOSUM_MATRICES, 0,
                    "BLOSUM number matrix is used. Default: 62",
                    0),
            
            #JTT PAM number (Jones et al. 1992) matrix is used. number>0.
            #Default: BLOSUM62
            _Option(["--jtt", "jtt"], ["input"], None, 0,
                    "JTT PAM number (Jones et al. 1992) matrix is used. " + \
                    "number>0. Default: BLOSUM62",
                    0),

            #Transmembrane PAM number (Jones et al. 1994) matrix is used.
            #number>0. Default: BLOSUM62
            _Option(["--tm", "tm"], ["input"],
                    os.path.exists, 0,
                    "Transmembrane PAM number (Jones et al. 1994) " + \
                    "matrix is used. number>0. Default: BLOSUM62",
                    0),

            #Use a user-defined AA scoring matrix. The format of matrixfile is
            #the same to that of BLAST. Ignored when nucleotide sequences are
            #input. Default: BLOSUM62
            _Option(["--aamatrix", "aamatrix"], ["input"],
                    os.path.exists, 0,
                    "Use a user-defined AA scoring matrix. " + \
                    "Default: BLOSUM62",
                    0),

            #Incorporate the AA/nuc composition information into the scoring
            #matrix. Default: off
            _Option(["--fmodel", "fmodel"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Incorporate the AA/nuc composition information " + \
                    "into the scoring matrix. Default: off",
                    0),
         
            #**** Output ****
            #Output format: clustal format. Default: off (fasta format)
            _Option(["--clustalout", "clustalout"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Output format: clustal format. Default: off (fasta" + \
                    "format)",
                    0),

            #Output order: same as input. Default: on
            _Option(["--inputorder", "inputorder"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Output order: same as input. Default: on",
                    0),

            #Output order: aligned. Default: off (inputorder)
            _Option(["--reorder", "reorder"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Output order: aligned. Default: off (inputorder)",
                    0),
            
            #Guide tree is output to the input.tree file. Default: off
            _Option(["--treeout", "treeout"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Guide tree is output to the input.tree file. Default: off",
                    0),

            #Do not report progress. Default: off
            _Option(["--quiet", "quiet"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Do not report progress. Default: off",
                    0),
         
            #**** Input ****

            #Assume the sequences are nucleotide. Deafult: auto
            _Option(["--nuc", "nuc"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Assume the sequences are nucleotide. Default: auto",
                    0),

            #Assume the sequences are amino acid. Deafult: auto
            _Option(["--amino", "amino"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Assume the sequences are amino acid. Default: auto",
                    0),

            ###################### SEEDS #####################################
            # MAFFT has multiple --seed commands where the unaligned input is
            # aligned to the seed alignment. There can be multiple seeds in the
            # form: "mafft --seed align1 --seed align2 [etc] input"
            # Effectively for n number of seed alignments. Here we're going to
            # assume 6 extra are enough
            _Option(["--seed", "seed"], ["input"], os.path.exists, 0,
                    "Seed alignments given in alignment_n (fasta format) " + \
                    "are aligned with sequences in input.",
                    0),

            #SEED 1
            _Option(["--seed", "seed1"], ["input"], os.path.exists, 0,
                    "Seed alignments given in alignment_n (fasta format) " + \
                    "are aligned with sequences in input.",
                    0),

            #SEED 2
            _Option(["--seed", "seed2"], ["input"], os.path.exists, 0,
                    "Seed alignments given in alignment_n (fasta format) " + \
                    "are aligned with sequences in input.",
                    0),

            #SEED 3
            _Option(["--seed", "seed3"], ["input"], os.path.exists, 0,
                    "Seed alignments given in alignment_n (fasta format) " + \
                    "are aligned with sequences in input.",
                    0),

            #SEED 4
            _Option(["--seed", "seed4"], ["input"], os.path.exists, 0,
                    "Seed alignments given in alignment_n (fasta format) " + \
                    "are aligned with sequences in input.",
                     False),

            #SEED 5
            _Option(["--seed", "seed5"], ["input"], os.path.exists, 0,
                    "Seed alignments given in alignment_n (fasta format) " + \
                    "are aligned with sequences in input.",
                    0),

            #SEED 6
            _Option(["--seed", "seed6"], ["input"], os.path.exists, 0,
                    "Seed alignments given in alignment_n (fasta format) " + \
                    "are aligned with sequences in input.",
                    0),
            ####################### END SEEDS  ################################

            #The input (must be FASTA format)
            _Argument(["input"], ["input"], os.path.exists, 1,
                      "Input file name"),

            ###################################################################
            #mafft-profile takes a second alignment input as an argument:
            #mafft-profile align1 align2
            _Argument(["input1"], ["input"], os.path.exists, 0,
                      "Second input file name for the mafft-profile command")
            ]



