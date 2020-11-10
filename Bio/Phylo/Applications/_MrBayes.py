# Based on code in _Phyml.py by Eric Talevich.
# All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Command-line wrapper for Bayesian Inference of Phylogeny program MrBayes"""

from Bio.Application import _Option, _Switch, _Argument, AbstractCommandline

class MrBayesCommandline(AbstractCommandline):
    """Command-line wrapper for the Bayesian Inference of Phylogeny program MrBayes.

        Homepage: http://nbisweden.github.io/MrBayes/index.html

        Huelsenbeck, J.P., and F. Ronquist. 2001. 
        MRBAYES: Bayesian inference of phylogeny. Bioinformatics 17:754-755.

        Ronquist, F., and J.P. Huelsenbeck. 2003. 
        MRBAYES 3: Bayesian phylogenetic inference under mixed models. 
        Bioinformatics 19:1572-1574.

        Ronquist, F., M. Teslenko, P. van der Mark, D.L. Ayres, A. Darling,
        S. Höhna, B. Larget, L. Liu, M.A. Suchard, and J.P. Huelsenbeck. 2012. 
        MRBAYES 3.2: Efficient Bayesian phylogenetic inference and model selection across a large model space. 
        Syst. Biol. 61:539-542.

    """

    def _init_(self, cmd="mb", **kwargs):
        """initialize the class"""
        self.parameters = [
            _Switch(
                ["About", "about"],
                """Provides general information about the program""",
                ),
            _Switch(    
                ["Acknoledgments", "acknoledgements"],
                """Shows the authors acknoledgments"""
                ),
            _Option(
                ["Calibrate", "calibrate"],
                """Dates a terminal or interior node in the tree"""

#              checker_function=lambda x: x in(
#                    "fixed",
#                    "uniform",
#                    "offsetexponential",
#                    "truncatednormal",
#                    "lognormal",
#                    "offsetlognormal",
#                    "gamma",
#                    "offsetgamma") 
            
                ),
            _Option(
                ["Charset", "charset"],
                """Defines a character set, this command is best used
                    not from the command line but rather as a line in the
                    MrBayes block of a file.""",

                ),
            _Switch(
                ["Citations", "citations"],
                """Shows a thorough list of citations you may consider using when
                    publishing the results of a MrBayes analysis""",
                ),
            _Switch(
                ["Plot", "plot"],
                """Plots parameters from MCMC analysis""",
                ),


            ]

        AbstractCommandline.__init__(self, cmd, **kwargs)

