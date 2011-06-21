# Copyright (C) 2011 by Brandon Invergo (b.invergo@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

import unittest
import os
import os.path
import sys
from Bio.Phylo.PAML import baseml
from Bio.Phylo.PAML._paml import PamlError

class ModTest(unittest.TestCase):
    
    align_file = os.path.join("PAML", "alignment.phylip")
    tree_file = os.path.join("PAML", "species.tree")
    bad_tree_file = os.path.join("PAML", "bad.tree")
    out_file = os.path.join("PAML", "test.out")
    working_dir = os.path.join("PAML", "baseml_test")
    results_file = os.path.join("PAML", "bad_results.out")
    bad_ctl_file1 = os.path.join("PAML", "bad1.ctl")
    bad_ctl_file2 = os.path.join("PAML", "bad2.ctl")
    ctl_file = os.path.join("PAML", "baseml.ctl")
       
    def __del__(self):
        """Just in case BASEML creates some junk files, do a clean-up."""
        del_files = [self.out_file, "2base.t", 
            "in.basemlg", "baseml.ctl", "lnf", "rates", "rst", "rst1", 
            "rub"]
        for filename in del_files:
            if os.path.exists(filename):
                os.remove(filename)
        if os.path.exists(self.working_dir):
            for filename in os.listdir(self.working_dir):
                os.remove(filename)
            os.rmdir(self.working_dir)
    
    def setUp(self):
        self.bml = baseml.Baseml()
        
    def testAlignmentFileIsValid(self):
        self.assertRaises((AttributeError, TypeError),
            baseml.Baseml, alignment = 1)
        self.bml.alignment = 1
        self.bml.tree = self.tree_file
        self.bml.out_file = self.out_file
        self.assertRaises((AttributeError, TypeError),
            self.bml.run)
        
    def testAlignmentExists(self):
        self.assertRaises((EnvironmentError, IOError), baseml.Baseml, 
            alignment = "nonexistent")
        self.bml.alignment = "nonexistent"
        self.bml.tree = self.tree_file
        self.bml.out_file = self.out_file
        self.assertRaises((EnvironmentError, IOError), 
            self.bml.run)
    
    def testTreeFileValid(self):
        self.assertRaises((AttributeError, TypeError),
            baseml.Baseml, tree = 1)
        self.bml.alignment = self.align_file
        self.bml.tree = 1
        self.bml.out_file = self.out_file
        self.assertRaises((AttributeError, TypeError),
            self.bml.run)
        
    def testTreeExists(self):
        self.assertRaises((EnvironmentError, IOError), baseml.Baseml, 
            tree = "nonexistent")
        self.bml.alignment = self.align_file
        self.bml.tree = "nonexistent"
        self.bml.out_file = self.out_file
        self.assertRaises((EnvironmentError, IOError),
            self.bml.run)
    
    def testWorkingDirValid(self):
        self.bml.tree = self.tree_file
        self.bml.alignment = self.align_file
        self.bml.out_file = self.out_file
        self.bml.working_dir = 1
        self.assertRaises((AttributeError, TypeError),
            self.bml.run)
    
    def testOutputFileValid(self):
        self.bml.tree = self.tree_file
        self.bml.alignment = self.align_file
        self.bml.out_file = 1
        self.assertRaises((AttributeError, TypeError),
            self.bml.run)
    
    def testOptionExists(self):
        self.assertRaises((AttributeError, KeyError),
            self.bml.set_option, "xxxx", 1)
        self.assertRaises((AttributeError, KeyError),
            self.bml.get_option, "xxxx")
    
    def testAlignmentSpecified(self):
        self.bml.tree = self.tree_file
        self.bml.out_file = self.out_file
        self.assertRaises((AttributeError, ValueError),
            self.bml.run)
        
    def testTreeSpecified(self):
        self.bml.alignment = self.align_file
        self.bml.out_file = self.out_file
        self.assertRaises((AttributeError, ValueError),
            self.bml.run)
        
    def testOutputFileSpecified(self):
        self.bml.alignment = self.align_file
        self.bml.tree = self.tree_file
        self.assertRaises((AttributeError, ValueError),
            self.bml.run)
        
    def testPamlErrorsCaught(self):
        self.bml.alignment = self.align_file
        self.bml.tree = self.bad_tree_file
        self.bml.out_file = self.out_file
        self.assertRaises((EnvironmentError, PamlError),
            self.bml.run)
        
    def testCtlFileValidOnRun(self):
        self.bml.alignment = self.align_file
        self.bml.tree = self.tree_file
        self.bml.out_file = self.out_file
        self.assertRaises((AttributeError, TypeError),
            self.bml.run, ctl_file = 1)
        
    def testCtlFileExistsOnRun(self):
        self.bml.alignment = self.align_file
        self.bml.tree = self.tree_file
        self.bml.out_file = self.out_file
        self.assertRaises((EnvironmentError, IOError),
            self.bml.run, ctl_file = "nonexistent")
            
    def testCtlFileValidOnRead(self):
        self.assertRaises((AttributeError, TypeError),
            self.bml.read_ctl_file, 1)
        self.assertRaises((AttributeError, KeyError), 
            self.bml.read_ctl_file, self.bad_ctl_file1)
        self.assertRaises(AttributeError, 
            self.bml.read_ctl_file, self.bad_ctl_file2)
        target_options = {"noisy": 0,
                        "verbose": 0,
                        "runmode": 0,
                        "model": 6,
                        "model_options": None,
                        "Mgene": 1,
                        "ndata": None,
                        "clock": 0,
                        "fix_kappa": 0,
                        "kappa": 5,
                        "fix_alpha": 0,
                        "alpha": 0.5,
                        "Malpha": 1,
                        "ncatG": 5,
                        "fix_rho": 1,
                        "rho": 0,
                        "nparK": 0,
                        "nhomo": 0,
                        "getSE": 0,
                        "RateAncestor": 0,
                        "Small_Diff": 7e-6,
                        "cleandata": 1,
                        "icode": None,
                        "fix_blength": None,
                        "method": 0}
        self.bml.read_ctl_file(self.ctl_file)
        self.assertEqual(self.bml._options, target_options)
        
    def testCtlFileExistsOnRead(self):
        self.assertRaises(IOError,
            self.bml.read_ctl_file, ctl_file = "nonexistent")
        
    def testResultsValid(self):
        self.assertRaises((AttributeError, TypeError),
            baseml.read, 1)
    
    def testResultsExist(self):
        self.assertRaises(IOError, baseml.read, "nonexistent")
        
    def testResultsParsable(self):
        self.assertRaises(ValueError, baseml.read, self.results_file)
        
    def testParseAllVersions(self):
        folder = os.path.join("PAML","Results", "baseml", "versions")
        for results_file in os.listdir(folder):
            file_path = os.path.join(folder, results_file)
            if os.path.isfile(file_path) and results_file[:6] == "baseml":
                results = baseml.read(file_path)
                self.assertEqual(len(results), 6)
                self.assertEqual(len(results["parameters"]), 7)
    
    def testParseSEs(self):
        SE_results_file = os.path.join("PAML", "Results", "baseml",
                                       "baseml_SE.out")
        SE_results = baseml.read(SE_results_file)
        SE_parameters = SE_results.get("parameters")
        self.assertNotEqual(SE_parameters.get("SEs"), None)
        

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)

