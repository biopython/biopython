#Based on test_Fasttree_tool.py by Nate Sutton
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for MrBayes tool."""

from Bio import MissingExternalDependencyError

import sys
import os
import itertools
import unittest

from io import StringIO

from Bio import SeqIO
from Bio import Phylo
from Bio.Phylo.Applications import MrBayesCommandline
from Bio.Phylo.Applications import _MrBayes
from Bio.Application import ApplicationError


os.environ["LANG"] = "C"

mb_exe = None

if sys.platform == "win32":
    try:
        # This can vary depending on the Windows language.
        prog_files = os.environ["PROGRAMFILES"]
    except KeyError:
        prog_files = r"C:\Program Files (x86)"

        # Default MrBayes file path of "C:\Program Files (x86)\mb.exe

        likely_dirs = ["", "mrbayes"]
        likely_exes = ["mb.exe"]     # might be mb.exe need to verify

        for folder in likely_dirs:
            if os.path.isdir(os.path.join(program_files, folder)):
                for filename in likely_exes:
                    if os.path.isfile(os.path.join(prog_file, folder, filename)):
                        mb_exe = os.path.joinJ(prog_files, folder, filename)
                        break
                if mb_exe:
                    break

    
if not mb_exe:
    raise MissingExternalDependencyError(
        "Install MrBayes and correctly set the life path to the program"
        "if you want to use it from Biopython"
        )

class MrBayesTestCase(unittest.TestCase):
    """Test for application wrappers"""

    def test_mb(self):
        """Run MrBayes using the wrapper"""
        cmd = MrBayesCommandline(mb_exe)