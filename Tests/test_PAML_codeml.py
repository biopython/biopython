# Copyright (C) 2011 by Brandon Invergo (b.invergo@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

import unittest
import os
import os.path
from Bio.Phylo.PAML import codeml
from Bio.Phylo.PAML._paml import PamlError


# Some constants to assist with testing:
# This is the number of parameters that should be parsed for each 
# NSsites site class model
SITECLASS_PARAMS = {0: 6, 1: 4, 2: 4, 3: 4, 7: 5, 8: 9}
# This is the default number of site classes per NSsites site
# class model
SITECLASSES = {0: None, 1: 2, 2: 3, 3: 3, 7: 10, 8: 11}


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
                          self.cml.set_options, xxxx=1)
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
                        "NSsites": [0],
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
        self.assertEqual(sorted(self.cml._options.keys()), sorted(target_options.keys()))
        for key in target_options:
            self.assertEqual(self.cml._options[key], target_options[key], \
                             "%s: %r vs %r" \
                             % (key, self.cml._options[key], target_options[key]))
        
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
        versions_path = os.path.join("PAML", "Results", "codeml", "versions")
        for results_file in versions_path:
            results_path = os.path.join(versions_path, results_file)
            if os.path.isfile(results_path) and results_file[:6] == "codeml":
                results = codeml.read(results_path)
                # There should be 4 top-level items: 'codon model', 'model', 
                # 'version', & 'NSsites'
                self.assertEqual(len(results), 4)
                self.assertIn("NSsites", results)
                # There should be 6 NSsites classes: 0, 1, 2, 3, 7 & 8
                self.assertEqual(len(results["NSsites"]), 6)
                # Each site class model should have 5 sub-items: 'lnL', 'tree', 
                # 'description', 'parameters', & 'tree length'. It should
                # have the correct number of parameters also.
                for model_num in [0, 1, 2, 3, 7, 8]:
                    model = results["NSsites"][model_num]
                    self.assertEqual(len(model), 5)
                    self.assertIn("parameters", model)
                    params = model["parameters"]
                    self.assertEqual(len(params), SITECLASS_PARAMS[model_num])
                    self.assertIn("branches", params)
                    branches = params["branches"]
                    # There are 7 branches in the test case (specific to these
                    # test cases)
                    self.assertEqual(len(branches), 7)
                    if "site classes" in params:
                        self.assertEqual(len(params["site classes"]),
                                         SITECLASSES[model_num])
    
    def testParseSEs(self):
        results_file = os.path.join("PAML", "Results", "codeml",
            "codeml_SE.out")
        results = codeml.read(results_file)
        self.assertEqual(len(results), 4)
        self.assertIn("NSsites", results)
        models = results["NSsites"]
        # Only site class model 0 was simulated
        self.assertEqual(len(models), 1)
        self.assertIn(0, models)
        model = models[0]
        self.assertEqual(len(model), 5)
        self.assertIn("parameters", model)
        params = model["parameters"]
        # There should be one new item in the parameters, "SEs" 
        self.assertEqual(len(params), SITECLASS_PARAMS[0] + 1)
        self.assertIn("SEs", params)
    
    def testParseAllNSsites(self):
        # This test is redundant (see testParseAllVersions)
        results_file = os.path.join("PAML", "Results", "codeml",
            "codeml_NSsites_all.out")
        results = codeml.read(results_file)
        # There should be 4 top-level items: 'codon model', 'model', 
        # 'version', & 'NSsites'
        self.assertEqual(len(results), 4)
        self.assertIn("NSsites", results)
        # There should be 6 NSsites classes: 0, 1, 2, 3, 7 & 8
        self.assertEqual(len(results["NSsites"]), 6)
        # Each site class model should have 5 sub-items: 'lnL', 'tree', 
        # 'description', 'parameters', & 'tree length'. It should
        # have the correct number of parameters also.
        for model_num in [0, 1, 2, 3, 7, 8]:
            model = results["NSsites"][model_num]
            self.assertEqual(len(model), 5)
            self.assertIn("parameters", model)
            params = model["parameters"]
            self.assertEqual(len(params), SITECLASS_PARAMS[model_num])
            self.assertIn("branches", params)
            branches = params["branches"]
            # There are 7 branches in the test case (specific to these
            # test cases)
            self.assertEqual(len(branches), 7)
            if "site classes" in params:
                self.assertEqual(len(params["site classes"]),
                                 SITECLASSES[model_num])
        
    def testParseBranchSiteA(self):
        results_file = os.path.join("PAML", "Results", "codeml",
            "codeml_branchsiteA.out")
        results = codeml.read(results_file)
        # There are 5 top-level items in this case:
        # 'codon model', 'model', 'version', 'NSsites' & 'site-class model'
        self.assertEqual(len(results), 5)
        self.assertIn("NSsites", results)
        models = results["NSsites"]
        # Only site class model 2 is simulated for Branch Site A
        self.assertEqual(len(models), 1)
        self.assertIn(2, models)
        model = models[2]
        self.assertEqual(len(model), 5)
        self.assertIn("parameters", model)
        params = model["parameters"]
        # Branch Site A results lack a "branches" parameter     
        self.assertEqual(len(params), SITECLASS_PARAMS[2] - 1)
        self.assertIn("site classes", params)
        site_classes = params["site classes"]
        # Branch Site A adds another site class
        self.assertEqual(len(site_classes), SITECLASSES[2] + 1)        
        for class_num in [0, 1, 2, 3]:
            self.assertIn(class_num, site_classes)
            site_class = site_classes[class_num]
            self.assertEqual(len(site_class), 2)
            self.assertIn("branch types", site_class)
            branches = site_class["branch types"]
            self.assertEqual(len(branches), 2)
 
    def testParseCladeModelC(self):
        results_file = os.path.join("PAML", "Results", "codeml",
            "codeml_clademodelC.out")
        results = codeml.read(results_file)
        # 5 top-level items again in this case
        self.assertEqual(len(results), 5)
        self.assertIn("NSsites", results)
        models = results["NSsites"]
        # Only site class model 2 is simulated for Clade Model C
        self.assertEqual(len(models), 1)
        self.assertIn(2, models)
        model = models[2]
        self.assertEqual(len(model), 5)
        self.assertIn("parameters", model)
        params = model["parameters"]
        # Clade Model C results lack a "branches" parameter     
        self.assertEqual(len(params), SITECLASS_PARAMS[2] - 1)
        self.assertIn("site classes", params)
        site_classes = params["site classes"]
        self.assertEqual(len(site_classes), SITECLASSES[2])        
        for class_num in [0, 1, 2]:
            self.assertIn(class_num, site_classes)
            site_class = site_classes[class_num]
            self.assertEqual(len(site_class), 2)
            self.assertIn("branch types", site_class)
            branches = site_class["branch types"]
            self.assertEqual(len(branches), 2)
    
    def testParseNgene2Mgene012(self):
        results_file = os.path.join("PAML", "Results", "codeml",
            "codeml_ngene2_mgene012.out")
        results = codeml.read(results_file)
        self.assertEqual(len(results), 4)
        self.assertIn("NSsites", results)
        models = results["NSsites"]
        self.assertEqual(len(models), 1)
        self.assertIn(0, models)
        model = models[0]
        self.assertEqual(len(model), 5)
        self.assertIn("parameters", model)
        params = model["parameters"]
        # This type of model has fewer parameters for model 0
        self.assertEqual(len(params), 4)
        self.assertIn("rates", params)
        rates = params["rates"]
        self.assertEqual(len(rates), 2)        
    
    def testParseNgene2Mgene34(self):
        results_file = os.path.join("PAML", "Results", "codeml",
            "codeml_ngene2_mgene34.out")
        results = codeml.read(results_file)
        self.assertEqual(len(results), 4)
        self.assertIn("NSsites", results)
        models = results["NSsites"]
        self.assertEqual(len(models), 1)
        self.assertIn(0, models)
        model = models[0]
        self.assertEqual(len(model), 5)
        self.assertIn("parameters", model)
        params = model["parameters"]
        # This type of model has fewer parameters for model 0
        self.assertEqual(len(params), 3)
        self.assertIn("rates", params)
        rates = params["rates"]
        self.assertEqual(len(rates), 2)        
        self.assertIn("genes", params)
        genes = params["genes"]
        self.assertEqual(len(genes), 2)   
    
    def testParseFreeBranch(self):
        results_file = os.path.join("PAML", "Results", "codeml",
            "codeml_freebranch.out")
        results = codeml.read(results_file)
        self.assertEqual(len(results), 4)
        self.assertIn("NSsites", results)
        models = results["NSsites"]
        self.assertEqual(len(models), 1)
        self.assertIn(0, models)
        model = models[0]
        # With the free ratio model, you get 3 extra trees: dN tree,
        # dS tree and omega tree
        self.assertEqual(len(model), 8)
        self.assertIn("parameters", model)
        params = model["parameters"]
        self.assertEqual(len(params), SITECLASS_PARAMS[0])
        self.assertIn("branches", params)
        # There should be 7 branches
        branches = params["branches"]
        self.assertEqual(len(branches), 7) 
        self.assertIn("omega", params)
        omega = params["omega"]
        self.assertEqual(len(omega), 7)
    
    def testParsePairwise(self):
        results_file = os.path.join("PAML", "Results", "codeml",
            "codeml_pairwise.out")
        results = codeml.read(results_file)
        # Pairwise models have an extra top-level item: pairwise
        self.assertEqual(len(results), 5)
        self.assertIn("pairwise", results)
        pairwise = results["pairwise"]
        self.assertEqual(len(pairwise), 5) 
    
    def testParseAA(self):
        results_file = os.path.join("PAML", "Results", "codeml",
            "codeml_aa_model0.out")
        results = codeml.read(results_file)
        # Amino Acid analysis has different top-levels:
        # 'NSsites', 'model', 'version', 'lnL max', 'distances'
        self.assertEqual(len(results), 5)
        self.assertIn("lnL max", results)
        self.assertIn("distances", results)
        distances = results["distances"]
        # non-pairwise AA analysis only gives raw distances
        self.assertEqual(len(distances), 1) 

    def testParseAAPairwise(self):
        results_file = os.path.join("PAML", "Results", "codeml",
            "codeml_aa_pairwise.out")
        results = codeml.read(results_file)
        # Pairwise AA analysis has one top-level fewer than non-pairwise
        self.assertEqual(len(results), 4)
        self.assertIn("lnL max", results)
        self.assertIn("distances", results)
        distances = results["distances"]
        # Pairwise AA analysis has ML & raw distances
        self.assertEqual(len(distances), 2) 
        
        
if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
