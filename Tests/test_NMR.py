# Copyright 2019 by Sergio Valqui. All rights reserved.
# Adapted from Doc/examples/nmr/simplepredict.py by Robert Bussell, Jr
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Unit tests for the Bio.NMR Module."""

import unittest
import tempfile
import os

from Bio.NMR import xpktools
from Bio.NMR import NOEtools


class NmrTests(unittest.TestCase):
    """Tests for NMR module."""

    def test_xpktools(self):
        """Self test for NMR.xpktools."""
        self.xpk_file = "NMR/noed.xpk"

        self.peaklist = xpktools.Peaklist(self.xpk_file)
        # Peaklist attributes checks
        self.assertEqual(self.peaklist.firstline, "label dataset sw sf ")
        self.assertEqual(self.peaklist.axislabels, "H1 15N2 N15 ")
        self.assertEqual(self.peaklist.dataset, "test.nv")
        self.assertEqual(self.peaklist.sw, "{1571.86 } {1460.01 } {1460.00 }")
        self.assertEqual(self.peaklist.sf, "{599.8230 } { 60.7860 } { 60.7860 }")
        self.assertEqual(
            self.peaklist.datalabels,
            " H1.L  H1.P  H1.W  H1.B  H1.E  "
            "H1.J  15N2.L  15N2.P  15N2.W  15N2.B  "
            "15N2.E  15N2.J  N15.L  N15.P  N15.W  "
            "N15.B  N15.E  N15.J  vol  int  stat ",
        )
        # Peaklist data check
        self.assertEqual(len(self.peaklist.data), 8)
        self.assertEqual(
            self.peaklist.data[0],
            "0  3.hn   8.853   0.021   0.010   ++   "
            "0.000   3.n   120.104   0.344   0.010   PP"
            "   0.000   3.n   120.117   0.344   0.010   "
            "PP   0.000  1.18200 1.18200 0",
        )
        self.assertEqual(
            self.peaklist.data[7],
            "8  10.hn   7.663   0.021   0.010   ++   "
            "0.000   10.n   118.341   0.324   0.010   "
            "+E   0.000   10.n   118.476   0.324   "
            "0.010   +E   0.000  0.49840 0.49840 0",
        )

        # Peaklist residue dict check
        self.assertEqual(len(self.peaklist.residue_dict("H1")["10"]), 1)
        self.assertEqual(
            self.peaklist.residue_dict("H1")["10"][0],
            "8  10.hn   7.663   0.021"
            "   0.010   ++   0.000   "
            "10.n   118.341   0.324   "
            "0.010   +E   0.000   10.n"
            "   118.476   0.324   0.010"
            "   +E   0.000  0.49840 "
            "0.49840 0",
        )

    def test_noetools(self):
        """Self test for NMR.NOEtools.

        Calculate and compare crosspeak peaklist files
        Adapted from Doc/examples/nmr/simplepredict.py by Robert Bussell, Jr.
        """
        self.xpk_i_file = os.path.join("NMR", "noed.xpk")
        # out_example.xpk is created by running Doc/examples/nmr/simplepredict.py
        # with noed.xpk as an input file.
        self.xpk_expected = os.path.join("NMR", "out_example.xpk")

        self.f_number, self.f_predicted = tempfile.mkstemp()
        os.close(self.f_number)

        # Calculate crosspeak from a peaklist input file
        # and save to a temporal file for comparison.
        try:
            self.peaklist = xpktools.Peaklist(self.xpk_i_file)
            self.res_dict = self.peaklist.residue_dict("H1")
            max_res = self.res_dict["maxres"]
            min_res = self.res_dict["minres"]

            self.peaklist.write_header(self.f_predicted)

            inc = 1  # The NOE increment (n where i->i+n and i->i-n are noes)
            count = 0  # A counter that number the output data lines in order
            res = min_res  # minimum residue number in the set
            out_list = []  # Holds the output data

            while res <= max_res:
                noe1 = NOEtools.predictNOE(self.peaklist, "15N2", "H1", res, res + inc)
                noe2 = NOEtools.predictNOE(self.peaklist, "15N2", "H1", res, res - inc)

                if noe1 != "":
                    noe1 = noe1 + "\012"
                    noe1 = xpktools.replace_entry(noe1, 1, count)
                    out_list.append(noe1)
                    count += 1

                    if noe2 != "":
                        noe2 = noe2 + "\012"
                        noe2 = xpktools.replace_entry(noe2, 1, count)
                        out_list.append(noe2)
                        count += 1
                res += 1

            # Open the output file and write the data
            with open(self.f_predicted, "a") as outfile:
                outfile.writelines(out_list)  # Write the output lines to the file

            # Compare the content of the predicted output file with expected file
            pre_content = open(self.f_predicted).read()
            exp_content = open(self.xpk_expected).read()
            self.assertEqual(pre_content, exp_content)

        finally:
            os.remove(self.f_predicted)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
