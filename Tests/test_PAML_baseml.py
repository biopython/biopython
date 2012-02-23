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
    align_dir = os.path.join("PAML", "Alignments")
    tree_dir = os.path.join("PAML", "Trees")
    ctl_dir = os.path.join("PAML", "Control_files")
    results_dir = os.path.join("PAML", "Results")
    working_dir = os.path.join("PAML", "baseml_test")
 
    align_file = os.path.join(align_dir, "alignment.phylip")
    tree_file = os.path.join(tree_dir, "species.tree")
    bad_tree_file = os.path.join(tree_dir, "bad.tree")
    out_file = os.path.join(results_dir, "test.out")
    results_file = os.path.join(results_dir, "bad_results.out")
    bad_ctl_file1 = os.path.join(ctl_dir, "bad1.ctl")
    bad_ctl_file2 = os.path.join(ctl_dir, "bad2.ctl")
    ctl_file = os.path.join(ctl_dir, "baseml.ctl")
       
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
                          self.bml.set_options, xxxx=1)
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
                        "Malpha": 0,
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
        self.assertEqual(sorted(self.bml._options.keys()), sorted(target_options.keys()))
        for key in target_options:
            self.assertEqual(self.bml._options[key], target_options[key], \
                             "%s: %r vs %r" \
                             % (key, self.bml._options[key], target_options[key]))
        
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
    
    def testParseModel(self):
        res_dir = os.path.join(self.results_dir, "baseml", "model")
        for results_file in os.listdir(res_dir):
            version = results_file.split('-')[1].split('.')[0]
            model = results_file[5]
            version_msg = "Improper parsing for model %s version %s" \
                        % (model, version.replace('_', '.'))
            results_path = os.path.join(res_dir, results_file)
            results = baseml.read(results_path)
            # There are 6 top-levels: parameters, tree, lnL, version,
            # tree length and lnL max
            self.assertEqual(len(results), 6, version_msg)
            self.assertTrue("parameters" in results, version_msg)
            params = results["parameters"]
            self.assertTrue("alpha" in params, version_msg)
            self.assertTrue("rates" in params, version_msg)
            self.assertTrue("parameter list" in params, version_msg)
            self.assertTrue("rate frequencies" in params, version_msg)
            if model in ["1", "3", "4", "5", "6"]:
                self.assertTrue("kappa" in params, version_msg)
            if model in ["7", "8"]:
                self.assertTrue("base frequencies" in params, version_msg)
                self.assertTrue("rate parameters" in params, version_msg)
                self.assertTrue("Q matrix" in params, version_msg)
                qmat = params["Q matrix"]
                self.assertEqual(len(qmat), 2, version_msg)
                self.assertTrue("matrix" in qmat)
                matrix = qmat["matrix"]
                self.assertEqual(len(matrix), 4, version_msg)
                self.assertEqual(len(matrix[0]), 4, version_msg)
            
    def testParseAlpha1Rho1(self):
        # Test the auto-discrete gamma model
        # Cannot test for baseml 4.3-4.5 due to bug in the program which
        # prevents this analysis from completing
        res_dir = os.path.join(self.results_dir, "baseml", "alpha1rho1")
        for results_file in os.listdir(res_dir):
            version = results_file.split('-')[1].split('.')[0]
            model = results_file[5]
            version_msg = "Improper parsing for model %s version %s" \
                        % (model, version.replace('_', '.'))
            results_path = os.path.join(res_dir, results_file)
            results = baseml.read(results_path)
            # There are 6 top-levels: parameters, tree, lnL, version,
            # tree length and lnL max
            self.assertEqual(len(results), 6, version_msg)
            self.assertTrue("parameters" in results, version_msg)
            params = results["parameters"]
            self.assertTrue("rho" in params, version_msg)
            self.assertTrue("transition probs." in params, version_msg)
            trans_p = params["transition probs."]
            self.assertEqual(len(trans_p), 5, version_msg)
            self.assertEqual(len(trans_p[0]), 5, version_msg)
 
    def testParseNhomo(self):
        res_dir = os.path.join(self.results_dir, "baseml", "nhomo")
        for results_file in os.listdir(res_dir):
            version = results_file.split('-')[1].split('.')[0]
            n = results_file[5]
            version_msg = "Improper parsing for nhomo %s version %s" \
                        % (n, version.replace('_', '.'))
            results_path = os.path.join(res_dir, results_file)
            results = baseml.read(results_path)
            # There are 6 top-levels: parameters, tree, lnL, version,
            # tree length and lnL max
            self.assertEqual(len(results), 6, version_msg)
            self.assertTrue("parameters" in results, version_msg)
            params = results["parameters"]
            if n == "1":
                self.assertTrue("base frequencies" in params, version_msg)
            else:
                self.assertTrue("nodes" in params)
                nodes = params["nodes"]
                self.assertEqual(len(nodes), 8, version_msg)
                self.assertEqual(len(nodes[1]), 2, version_msg)

    def testParseSEs(self):
        res_dir = os.path.join(self.results_dir, "baseml", "SE")
        for results_file in os.listdir(res_dir):
            version = results_file.split('-')[1].split('.')[0]
            version_msg = "Improper parsing for version %s" \
                        % version.replace('_', '.')
            results_path = os.path.join(res_dir, results_file)
            results = baseml.read(results_path)
            # There are 6 top-levels: parameters, tree, lnL, version,
            # tree length and lnL max
            self.assertEqual(len(results), 6, version_msg)
            self.assertTrue("parameters" in results, version_msg)
            params = results["parameters"]
            self.assertTrue("SEs" in params, version_msg)
        

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)

