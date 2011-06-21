# Copyright (C) 2011 by Brandon Invergo (b.invergo@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

import unittest
import os
import os.path
from Bio.Phylo.PAML import codeml
from Bio.Phylo.PAML._paml import PamlError

class ModTest(unittest.TestCase):
    
    align_file = os.path.join("PAML", "alignment.phylip")
    tree_file = os.path.join("PAML", "species.tree")
    bad_tree_file = os.path.join("PAML", "bad.tree")
    out_file = os.path.join("PAML", "test.out")
    working_dir = os.path.join("PAML", "codeml_test")
    results_file = os.path.join("PAML", "bad_results.out")
    bad_ctl_file1 = os.path.join("PAML", "bad1.ctl")
    bad_ctl_file2 = os.path.join("PAML", "bad2.ctl")
    bad_ctl_file3 = os.path.join("PAML", "bad3.ctl")
    ctl_file = os.path.join("PAML", "codeml.ctl")
       
    def __del__(self):
        """Just in case CODEML creates some junk files, do a clean-up."""
        del_files = [self.out_file, "2NG.dN", 
            "2NG.dS", "2NG.t", "codeml.ctl", "lnf", "rst", "rst1", 
            "rub"]
        for filename in del_files:
            if os.path.exists(filename):
                os.remove(filename)
        if os.path.exists(self.working_dir):
            for filename in os.listdir(self.working_dir):
                os.remove(filename)
            os.rmdir(self.working_dir)
    
    def setUp(self):
        self.cml = codeml.Codeml()
        
    def testAlignmentFileIsValid(self):
        self.assertRaises((AttributeError, TypeError),
            codeml.Codeml, alignment = 1)
        self.cml.alignment = 1
        self.cml.tree = self.tree_file
        self.cml.out_file = self.out_file
        self.assertRaises((AttributeError, TypeError),
            self.cml.run)
        
    def testAlignmentExists(self):
        self.assertRaises((EnvironmentError, IOError), codeml.Codeml, 
            alignment = "nonexistent")
        self.cml.alignment = "nonexistent"
        self.cml.tree = self.tree_file
        self.cml.out_file = self.out_file
        self.assertRaises(IOError, self.cml.run)
    
    def testTreeFileValid(self):
        self.assertRaises((AttributeError, TypeError),
            codeml.Codeml, tree = 1)
        self.cml.alignment = self.align_file
        self.cml.tree = 1
        self.cml.out_file = self.out_file
        self.assertRaises((AttributeError, TypeError),
            self.cml.run)
        
    def testTreeExists(self):
        self.assertRaises((EnvironmentError, IOError), codeml.Codeml, 
            tree = "nonexistent")
        self.cml.alignment = self.align_file
        self.cml.tree = "nonexistent"
        self.cml.out_file = self.out_file
        self.assertRaises(IOError, self.cml.run)
    
    def testWorkingDirValid(self):
        self.cml.tree = self.tree_file
        self.cml.alignment = self.align_file
        self.cml.out_file = self.out_file
        self.cml.working_dir = 1
        self.assertRaises((AttributeError, TypeError),
            self.cml.run)
    
    def testOutputFileValid(self):
        self.cml.tree = self.tree_file
        self.cml.alignment = self.align_file
        self.cml.out_file = 1
        self.assertRaises((AttributeError, TypeError),
            self.cml.run)
    
    def testOptionExists(self):
        self.assertRaises((AttributeError, KeyError),
            self.cml.set_option, "xxxx", 1)
        self.assertRaises((AttributeError, KeyError),
            self.cml.get_option, "xxxx")
    
    def testAlignmentSpecified(self):
        self.cml.tree = self.tree_file
        self.cml.out_file = self.out_file
        self.assertRaises((AttributeError, ValueError),
            self.cml.run)
        
    def testTreeSpecified(self):
        self.cml.alignment = self.align_file
        self.cml.out_file = self.out_file
        self.assertRaises((AttributeError, ValueError),
            self.cml.run)
        
    def testOutputFileSpecified(self):
        self.cml.alignment = self.align_file
        self.cml.tree = self.tree_file
        self.assertRaises((AttributeError, ValueError),
            self.cml.run)
        
    def testPamlErrorsCaught(self):
        self.cml.alignment = self.align_file
        self.cml.tree = self.bad_tree_file
        self.cml.out_file = self.out_file
        self.assertRaises((EnvironmentError, PamlError),
            self.cml.run)
        
    def testCtlFileValidOnRun(self):
        self.cml.alignment = self.align_file
        self.cml.tree = self.tree_file
        self.cml.out_file = self.out_file
        self.assertRaises((AttributeError, TypeError),
            self.cml.run, ctl_file = 1)
        
    def testCtlFileExistsOnRun(self):
        self.cml.alignment = self.align_file
        self.cml.tree = self.tree_file
        self.cml.out_file = self.out_file
        self.assertRaises((EnvironmentError, IOError),
            self.cml.run, ctl_file = "nonexistent")
            
    def testCtlFileValidOnRead(self):
        self.assertRaises((AttributeError, TypeError),
            self.cml.read_ctl_file, 1)
        self.assertRaises((AttributeError, KeyError), 
            self.cml.read_ctl_file, self.bad_ctl_file1)
        self.assertRaises(AttributeError, 
            self.cml.read_ctl_file, self.bad_ctl_file2)
        self.assertRaises(TypeError, 
            self.cml.read_ctl_file, self.bad_ctl_file3)
        target_options = {"noisy": 0,
                        "verbose": 0,
                        "runmode": 0,
                        "seqtype": 1, 
                        "CodonFreq": 2,
                        "ndata": None,
                        "clock": 0,
                        "aaDist": None,
                        "aaRatefile": None,
                        "model": 0,
                        "NSsites": [0, 1, 2, 7, 8],
                        "icode": 0,
                        "Mgene": 0,
                        "fix_kappa": 0,
                        "kappa": 4.54006,
                        "fix_omega": 0,
                        "omega": 1,
                        "fix_alpha": 1,
                        "alpha": 0,
                        "Malpha": 0,
                        "ncatG": None,
                        "getSE": 0,
                        "RateAncestor": 0,
                        "Small_Diff": None,
                        "cleandata": 1,
                        "fix_blength": 1,
                        "method": 0}
        self.cml.read_ctl_file(self.ctl_file)
        self.assertEqual(self.cml._options, target_options)
        
    def testCtlFileExistsOnRead(self):
        self.assertRaises((EnvironmentError, IOError),
            self.cml.read_ctl_file, ctl_file = "nonexistent")
        
    def testResultsValid(self):
        self.assertRaises((AttributeError, TypeError),
            codeml.read, 1)
    
    def testResultsExist(self):
        self.assertRaises((EnvironmentError, IOError),
            codeml.read, "nonexistent")
        
    def testResultsParsable(self):
        self.assertRaises(ValueError, codeml.read, self.results_file)
        
    def testParseAllVersions(self):
        for results_file in os.listdir(os.path.join("PAML",
                "Results","codeml","versions")):
            if os.path.isfile(results_file) and results_file[:6] == "codeml":
                results = codeml.read(os.path.join("PAML",
                    "Results", results_file))
                self.assertEqual(len(results["NSsites"]), 6)
                self.assertEqual(len(results["NSsites"][0]), 7)
                self.assertEqual(len(results["NSsites"][1]), 5)
                self.assertEqual(len(results["NSsites"][2]), 5)
                self.assertEqual(len(results["NSsites"][3]), 5)
                self.assertEqual(len(results["NSsites"][7]), 6)
                self.assertEqual(len(results["NSsites"][8]), 6)
    
    def testParseSEs(self):
        SE_results_file = os.path.join("PAML", "Results", "codeml",
            "codeml_SE.out")
        SE_results = codeml.read(SE_results_file)
        SE_models = SE_results.get("NSsites")
        for model in SE_models:
            SE_model = SE_models.get(model)
            SE_parameters = SE_model.get("parameters")
            self.assertNotEqual(SE_parameters.get("SEs"), None)
    
    def testParseAllNSsites(self):
        results_file = os.path.join("PAML", "Results", "codeml",
            "codeml_NSsites_all.out")
        results = codeml.read(results_file)
        models = results.get("NSsites")
        self.assertEqual(len(models), 6)
        for model in models:
            self.assertEqual(len(models.get(model)), 5)
        
    def testParseBranchSiteA(self):
        results_file = os.path.join("PAML", "Results", "codeml",
            "codeml_branchsiteA.out")
        results = codeml.read(results_file)
        self.assertEqual(len(results), 5)
        site_classes = results["NSsites"][2]["parameters"]["site classes"]
        self.assertEqual(len(site_classes), 4)        
        
    def testParseCladeModelC(self):
        results_file = os.path.join("PAML", "Results", "codeml",
            "codeml_clademodelC.out")
        results = codeml.read(results_file)
        self.assertEqual(len(results), 5)
        site_classes = results["NSsites"][2]["parameters"]["site classes"]
        self.assertEqual(len(site_classes), 3)        
    
    def testParseNgene2Mgene012(self):
        results_file = os.path.join("PAML", "Results", "codeml",
            "codeml_ngene2_mgene012.out")
        results = codeml.read(results_file)
        self.assertEqual(len(results), 4)
        site_classes = results["NSsites"][0]["parameters"]["rates"]
        self.assertEqual(len(site_classes), 2)        
    
    def testParseNgene2Mgene34(self):
        results_file = os.path.join("PAML", "Results", "codeml",
            "codeml_ngene2_mgene34.out")
        results = codeml.read(results_file)
        self.assertEqual(len(results), 4)
        site_classes = results["NSsites"][0]["parameters"]["genes"]
        self.assertEqual(len(site_classes), 2)   
    
    def testParseFreeBranch(self):
        results_file = os.path.join("PAML", "Results", "codeml",
            "codeml_freebranch.out")
        results = codeml.read(results_file)
        self.assertEqual(len(results), 4)
        branches = results["NSsites"][0]["parameters"]["branches"]
        self.assertEqual(len(branches), 7) 
    
    def testParsePairwise(self):
        results_file = os.path.join("PAML", "Results", "codeml",
            "codeml_pairwise.out")
        results = codeml.read(results_file)
        self.assertEqual(len(results), 5)
        pairwise = results["pairwise"]
        self.assertEqual(len(pairwise), 5) 
    
    def testParseAA(self):
        results_file = os.path.join("PAML", "Results", "codeml",
            "codeml_aa_model0.out")
        results = codeml.read(results_file)
        self.assertEqual(len(results), 5)
        distances = results["distances"]
        self.assertEqual(len(distances), 1) 

    def testParseAAPairwise(self):
        results_file = os.path.join("PAML", "Results", "codeml",
            "codeml_aa_pairwise.out")
        results = codeml.read(results_file)
        self.assertEqual(len(results), 4)
        distances = results["distances"]
        self.assertEqual(len(distances), 2) 
        
        
if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
