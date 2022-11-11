# Copyright 2013 by Nate Sutton.  All rights reserved.
# Based on test_Clustalw_tool.py by Peter Cock.
# Example code used from Biopython's Phylo cookbook by Eric Talevich.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for Fasttree tool."""


from Bio import MissingExternalDependencyError

import sys
import os
import itertools
import unittest

from io import StringIO

from Bio import SeqIO
from Bio import Phylo
from Bio.Phylo.Applications import FastTreeCommandline
from Bio.Phylo.Applications import _Fasttree
from Bio.Application import ApplicationError

#################################################################

# Try to avoid problems when the OS is in another language
os.environ["LANG"] = "C"

fasttree_exe = None
if sys.platform == "win32":
    try:
        # This can vary depending on the Windows language.
        prog_files = os.environ["PROGRAMFILES"]
    except KeyError:
        prog_files = r"C:\Program Files (x86)"

    # A default fasttree file path of "C:\Program Files (x86)\Fasttree.exe"
    # was chosen here but users can alter the path according to where
    # fasttree is located on their systems

    likely_dirs = ["", "FastTree"]
    likely_exes = ["FastTree.exe"]
    for folder in likely_dirs:
        if os.path.isdir(os.path.join(prog_files, folder)):
            for filename in likely_exes:
                if os.path.isfile(os.path.join(prog_files, folder, filename)):
                    fasttree_exe = os.path.join(prog_files, folder, filename)
                    break
            if fasttree_exe:
                break
else:
    from subprocess import getoutput

    # Website uses 'FastTree', Nate's system had 'fasttree'
    likely_exes = ["FastTree", "fasttree"]
    for filename in likely_exes:
        # Checking the -help argument
        output = getoutput(f"{filename} -help")
        # Since "is not recognized" may be in another language, try and be sure this
        # is really the fasttree tool's output
        if (
            "is not recognized" not in output
            and "protein_alignment" in output
            and "nucleotide_alignment" in output
        ):
            fasttree_exe = filename
            break

if not fasttree_exe:
    raise MissingExternalDependencyError(
        "Install FastTree and correctly set the file path to the program "
        "if you want to use it from Biopython."
    )


class FastTreeTestCase(unittest.TestCase):
    def check(self, path, length):
        input_records = SeqIO.to_dict(SeqIO.parse(path, "fasta"))
        self.assertEqual(len(input_records), length)
        # Any filenames with spaces should get escaped with quotes
        #  automatically.
        # Using keyword arguments here.
        cline = _Fasttree.FastTreeCommandline(fasttree_exe, input=path, nt=True)
        self.assertEqual(str(eval(repr(cline))), str(cline))
        out, err = cline()
        self.assertTrue(err.strip().startswith("FastTree"))
        tree = Phylo.read(StringIO(out), "newick")

        names = {}
        for clade in tree.find_clades():
            if clade.name:
                self.assertNotIn(clade.name, names)
                names[clade.name] = clade

        self.assertGreater(len(names), 0)

        def terminal_neighbor_dists(self):
            """Return a list of distances between adjacent terminals."""

            def generate_pairs(self):
                pairs = itertools.tee(self)
                next(pairs[1])  # Advance second iterator one step
                return zip(pairs[0], pairs[1])

            return [
                self.distance(*i)
                for i in generate_pairs(self.find_clades(terminal=True))
            ]

        for dist in terminal_neighbor_dists(tree):
            self.assertGreater(dist, 0.0)

    def test_normal(self):
        self.check("Quality/example.fasta", 3)

    def test_filename_spaces(self):
        path = "Clustalw/temp horses.fasta"  # note spaces in filename
        records = SeqIO.parse("Phylip/hennigian.phy", "phylip")
        with open(path, "w") as handle:
            length = SeqIO.write(records, handle, "fasta")
        self.assertEqual(length, 10)
        self.check(path, length)

    def test_invalid(self):
        path = "Medline/pubmed_result1.txt"
        cline = FastTreeCommandline(fasttree_exe, input=path)
        with self.assertRaises(ApplicationError) as cm:
            stdout, stderr = cline()
        message = str(cm.exception)
        self.assertTrue(
            "invalid format" in message
            or "not produced" in message
            or "No sequences in file" in message
            or "Error parsing header line:" in message
            or "Non-zero return code " in message,
            msg=f"Unknown ApplicationError raised: {message}",
        )

    def test_single(self):
        path = "Fasta/f001"
        records = list(SeqIO.parse(path, "fasta"))
        self.assertEqual(len(records), 1)
        cline = FastTreeCommandline(fasttree_exe, input=path)
        stdout, stderr = cline()
        self.assertIn("Unique: 1/1", stderr)

    def test_empty(self):
        path = "does_not_exist.fasta"
        cline = FastTreeCommandline(fasttree_exe, input=path)
        with self.assertRaises(ApplicationError) as cm:
            stdout, stderr = cline()
        message = str(cm.exception)
        self.assertTrue(
            "Cannot open sequence file" in message
            or "Cannot open sequence file" in message
            or f"Cannot read {path}" in message
            or "Non-zero return code " in message,
            msg=f"Unknown ApplicationError raised: {message}",
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
