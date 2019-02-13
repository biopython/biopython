# Copyright (C) 2011 by Brandon Invergo (b.invergo@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

import unittest
import os
import os.path
from Bio.Phylo.PAML import yn00
from Bio.Phylo.PAML._paml import PamlError
import glob


class ModTest(unittest.TestCase):

    align_dir = os.path.join("PAML", "Alignments")
    tree_dir = os.path.join("PAML", "Trees")
    ctl_dir = os.path.join("PAML", "Control_files")
    results_dir = os.path.join("PAML", "Results")
    working_dir = os.path.join("PAML", "yn00_test")

    align_file = os.path.join(align_dir, "alignment.phylip")
    out_file = os.path.join(results_dir, "test.out")
    results_file = os.path.join(results_dir, "bad_results.out")
    bad_ctl_file1 = os.path.join(ctl_dir, "bad1.ctl")
    bad_ctl_file2 = os.path.join(ctl_dir, "bad2.ctl")
    ctl_file = os.path.join(ctl_dir, "yn00", "yn00.ctl")

    def tearDown(self):
        """Just in case yn00 creates some junk files, do a clean-up."""
        del_files = [self.out_file, "2YN.dN", "2YN.dS", "2YN.t", "rst",
                     "rst1", "rub"]
        for filename in del_files:
            if os.path.exists(filename):
                os.remove(filename)
        if os.path.exists(self.working_dir):
            for filename in os.listdir(self.working_dir):
                filepath = os.path.join(self.working_dir, filename)
                os.remove(filepath)
            os.rmdir(self.working_dir)

    def setUp(self):
        self.yn00 = yn00.Yn00()

    def testAlignmentFileIsValid(self):
        self.assertRaises((AttributeError, TypeError, OSError),
                          yn00.Yn00, alignment=list())
        self.yn00.alignment = list()
        self.yn00.out_file = self.out_file
        self.assertRaises((AttributeError, TypeError, OSError),
                          self.yn00.run)

    def testAlignmentExists(self):
        self.assertRaises((EnvironmentError, IOError), yn00.Yn00,
                          alignment="nonexistent")
        self.yn00.alignment = "nonexistent"
        self.yn00.out_file = self.out_file
        self.assertRaises(IOError, self.yn00.run)

    def testWorkingDirValid(self):
        self.yn00.alignment = self.align_file
        self.yn00.out_file = self.out_file
        self.yn00.working_dir = list()
        self.assertRaises((AttributeError, TypeError, OSError),
                          self.yn00.run)

    def testOptionExists(self):
        self.assertRaises((AttributeError, KeyError),
                          self.yn00.set_options, xxxx=1)
        self.assertRaises((AttributeError, KeyError),
                          self.yn00.get_option, "xxxx")

    def testAlignmentSpecified(self):
        self.yn00.out_file = self.out_file
        self.assertRaises((AttributeError, ValueError),
                          self.yn00.run)

    def testOutputFileSpecified(self):
        self.yn00.alignment = self.align_file
        self.assertRaises((AttributeError, ValueError),
                          self.yn00.run)

    # def testPamlErrorsCaught(self):
    #     self.yn00.alignment = self.align_file
    #     self.yn00.out_file = self.out_file
    #     self.assertRaises((EnvironmentError, PamlError),
    #                        self.yn00.run)

    def testCtlFileValidOnRun(self):
        self.yn00.alignment = self.align_file
        self.yn00.out_file = self.out_file
        self.assertRaises((AttributeError, TypeError, OSError),
                          self.yn00.run, ctl_file=list())

    def testCtlFileExistsOnRun(self):
        self.yn00.alignment = self.align_file
        self.yn00.out_file = self.out_file
        self.assertRaises(IOError,
                          self.yn00.run, ctl_file="nonexistent")

    def testCtlFileValidOnRead(self):
        self.assertRaises((AttributeError, TypeError, OSError),
                          self.yn00.read_ctl_file, list())
        self.assertRaises((AttributeError, KeyError),
                          self.yn00.read_ctl_file, self.bad_ctl_file1)
        self.assertRaises(AttributeError,
                          self.yn00.read_ctl_file, self.bad_ctl_file2)
        target_options = {"verbose": 1,
                          "icode": 0,
                          "weighting": 0,
                          "commonf3x4": 0,
                          "ndata": 1}
        self.yn00.read_ctl_file(self.ctl_file)
        self.assertEqual(self.yn00._options, target_options)

    def testCtlFileExistsOnRead(self):
        self.assertRaises(IOError,
                          self.yn00.read_ctl_file, ctl_file="nonexistent")

    def testResultsValid(self):
        self.assertRaises((AttributeError, TypeError, OSError),
                          yn00.read, list())

    def testResultsExist(self):
        self.assertRaises((EnvironmentError, IOError),
                          yn00.read, "nonexistent")

    def testResultsParsable(self):
        self.assertRaises(ValueError, yn00.read, self.results_file)

    def testParseAllVersions(self):
        pattern = os.path.join(self.results_dir, "yn00", 'yn00-*')
        for results_file in glob.glob(pattern):
            results = yn00.read(results_file)
            self.assertEqual(len(results), 5)
            self.assertEqual(len(results["Homo_sapie"]), 4)
            self.assertEqual(len(results["Homo_sapie"]["Pan_troglo"]), 5)

    def testParseLongNames(self):
        pattern = os.path.join(self.results_dir, "yn00", 'yn00_long-*')
        for results_file in glob.glob(pattern):
            results = yn00.read(results_file)
            # Expect seven taxa...
            self.assertEqual(len(results), 7)
            # ...each of which is compared to the other six.
            self.assertEqual(set([len(v) for v in results.values()]), set([6]))
            # ...each of which has five measures.
            self.assertEqual(set([len(v) for taxa in results.values() for v in taxa.values()]), set([5]))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
