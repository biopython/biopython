# Copyright (C) 2011, 2015, 2016 by Brandon Invergo (b.invergo@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

import unittest
import os
import os.path
import itertools

from Bio.Phylo.PAML import codeml
from Bio.Phylo.PAML._paml import PamlError


# Some constants to assist with testing:
# This is the number of parameters that should be parsed for each
# NSsites site class model
SITECLASS_PARAMS = {0: 6, 1: 4, 2: 4, 3: 4, 7: 5, 8: 8, 22: 4}
# This is the default number of site classes per NSsites site
# class model
SITECLASSES = {0: None, 1: 2, 2: 3, 3: 3, 7: 10, 8: 11, 22: 3}


class ModTest(unittest.TestCase):

    align_dir = os.path.join("PAML", "Alignments")
    tree_dir = os.path.join("PAML", "Trees")
    ctl_dir = os.path.join("PAML", "Control_files")
    results_dir = os.path.join("PAML", "Results")
    working_dir = os.path.join("PAML", "codeml_test")

    align_file = os.path.join(align_dir, "alignment.phylip")
    tree_file = os.path.join(tree_dir, "species.tree")
    bad_tree_file = os.path.join(tree_dir, "bad.tree")
    out_file = os.path.join(results_dir, "test.out")
    results_file = os.path.join(results_dir, "bad_results.out")
    bad_ctl_file1 = os.path.join(ctl_dir, "bad1.ctl")
    bad_ctl_file2 = os.path.join(ctl_dir, "bad2.ctl")
    bad_ctl_file3 = os.path.join(ctl_dir, "bad3.ctl")
    ctl_file = os.path.join(ctl_dir, "codeml", "codeml.ctl")

    def tearDown(self):
        """Just in case CODEML creates some junk files, do a clean-up."""
        del_files = [self.out_file, "2NG.dN", "2NG.dS", "2NG.t", "codeml.ctl",
                     "lnf", "rst", "rst1", "rub"]
        for filename in del_files:
            if os.path.exists(filename):
                os.remove(filename)
        if os.path.exists(self.working_dir):
            for filename in os.listdir(self.working_dir):
                filepath = os.path.join(self.working_dir, filename)
                os.remove(filepath)
            os.rmdir(self.working_dir)

    def setUp(self):
        self.cml = codeml.Codeml()

    def testAlignmentFileIsValid(self):
        self.assertRaises((AttributeError, TypeError, OSError),
                          codeml.Codeml, alignment=list())
        self.cml.alignment = list()
        self.cml.tree = self.tree_file
        self.cml.out_file = self.out_file
        self.assertRaises((AttributeError, TypeError, OSError),
                          self.cml.run)

    def testAlignmentExists(self):
        self.assertRaises((EnvironmentError, IOError), codeml.Codeml,
                          alignment="nonexistent")
        self.cml.alignment = "nonexistent"
        self.cml.tree = self.tree_file
        self.cml.out_file = self.out_file
        self.assertRaises(IOError, self.cml.run)

    def testTreeFileValid(self):
        self.assertRaises((AttributeError, TypeError, OSError),
                          codeml.Codeml, tree=list())
        self.cml.alignment = self.align_file
        self.cml.tree = list()
        self.cml.out_file = self.out_file
        self.assertRaises((AttributeError, TypeError, OSError),
                          self.cml.run)

    def testTreeExists(self):
        self.assertRaises((EnvironmentError, IOError), codeml.Codeml,
                          tree="nonexistent")
        self.cml.alignment = self.align_file
        self.cml.tree = "nonexistent"
        self.cml.out_file = self.out_file
        self.assertRaises(IOError, self.cml.run)

    def testWorkingDirValid(self):
        self.cml.tree = self.tree_file
        self.cml.alignment = self.align_file
        self.cml.out_file = self.out_file
        self.cml.working_dir = list()
        self.assertRaises((AttributeError, TypeError, OSError),
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
        self.assertRaises((AttributeError, TypeError, OSError),
                          self.cml.run, ctl_file=list())

    def testCtlFileExistsOnRun(self):
        self.cml.alignment = self.align_file
        self.cml.tree = self.tree_file
        self.cml.out_file = self.out_file
        self.assertRaises((EnvironmentError, IOError),
                          self.cml.run, ctl_file="nonexistent")

    def testCtlFileValidOnRead(self):
        self.assertRaises((AttributeError, TypeError, OSError),
                          self.cml.read_ctl_file, list())
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
                          "method": 0,
                          "rho": None,
                          "fix_rho": None}
        self.cml.read_ctl_file(self.ctl_file)
        # Compare the dictionary keys:
        self.assertEqual(sorted(self.cml._options), sorted(target_options))
        for key in target_options:
            self.assertEqual(self.cml._options[key], target_options[key],
                             "%s: %r vs %r"
                             % (key, self.cml._options[key], target_options[key]))

    def testCtlFileExistsOnRead(self):
        self.assertRaises((EnvironmentError, IOError),
                          self.cml.read_ctl_file, ctl_file="nonexistent")

    def testResultsValid(self):
        self.assertRaises((AttributeError, TypeError, OSError),
                          codeml.read, list())

    def testResultsExist(self):
        self.assertRaises((EnvironmentError, IOError),
                          codeml.read, "nonexistent")

    def testResultsParsable(self):
        self.assertRaises(ValueError, codeml.read, self.results_file)

    def testParseSEs(self):
        res_dir = os.path.join(self.results_dir, "codeml", "SE")
        for results_file in os.listdir(res_dir):
            version = results_file.split('-')[1].split('.')[0]
            version_msg = ("Improper parsing for version %s"
                           % version.replace('_', '.'))
            results_path = os.path.join(res_dir, results_file)
            results = codeml.read(results_path)
            self.assertEqual(len(results), 4, version_msg)
            self.assertIn("NSsites", results, version_msg)
            models = results["NSsites"]
            # Only site class model 0 was simulated
            self.assertEqual(len(models), 1, version_msg)
            self.assertIn(0, models, version_msg)
            model = models[0]
            self.assertEqual(len(model), 5, version_msg)
            self.assertIn("parameters", model, version_msg)
            params = model["parameters"]
            # There should be one new item in the parameters, "SEs"
            self.assertEqual(len(params), SITECLASS_PARAMS[0] + 1, version_msg)
            self.assertIn("SEs", params, version_msg)

    def testParseAllNSsites(self):
        res_dir = os.path.join(self.results_dir, "codeml", "all_NSsites")
        for results_file in os.listdir(res_dir):
            version = results_file.split('-')[1].split('.')[0]
            version_msg = ("Improper parsing for version %s"
                           % version.replace('_', '.'))
            results_path = os.path.join(res_dir, results_file)
            results = codeml.read(results_path)
            # There should be 4 top-level items: 'codon model', 'model',
            # 'version', & 'NSsites'
            self.assertEqual(len(results), 4, version_msg)
            self.assertIn("NSsites", results, version_msg)
            # There should be 6 NSsites classes: 0, 1, 2, 3, 7 & 8
            self.assertEqual(len(results["NSsites"]), 6, version_msg)
            # Each site class model should have 5 sub-items: 'lnL', 'tree',
            # 'description', 'parameters', & 'tree length'. It should
            # have the correct number of parameters also.
            for model_num in [0, 1, 2, 3, 7, 8]:
                model = results["NSsites"][model_num]
                self.assertEqual(len(model), 5, version_msg)
                self.assertIn("parameters", model, version_msg)
                params = model["parameters"]
                self.assertEqual(len(params), SITECLASS_PARAMS[model_num],
                                 version_msg)
                self.assertIn("branches", params, version_msg)
                branches = params["branches"]
                # There are 7 branches in the test case (specific to these
                # test cases)
                self.assertEqual(len(branches), 7, version_msg)
                if "site classes" in params:
                    self.assertEqual(len(params["site classes"]),
                                     SITECLASSES[model_num], version_msg)

    def testParseNSsite3(self):
        res_dir = os.path.join(self.results_dir, "codeml", "NSsite3")
        for results_file in os.listdir(res_dir):
            version = results_file.split('-')[1].split('.')[0]
            version_msg = ("Improper parsing for version %s"
                           % version.replace('_', '.'))
            results_path = os.path.join(res_dir, results_file)
            results = codeml.read(results_path)
            # There should be 5 top-level items: 'codon model', 'model',
            # 'version', 'NSsites' & site-class model, the last of which
            # is only there when only one NSsites class is used
            self.assertEqual(len(results), 5, version_msg)
            self.assertIn('site-class model', results, version_msg)
            self.assertEqual(results['site-class model'], 'discrete',
                             version_msg)
            self.assertIn("NSsites", results, version_msg)
            # There should be 1 NSsites classe: 3
            self.assertEqual(len(results["NSsites"]), 1, version_msg)
            # Each site class model should have 5 sub-items: 'lnL', 'tree',
            # 'description', 'parameters', & 'tree length'. It should
            # have the correct number of parameters also.
            model = results["NSsites"][3]
            self.assertEqual(len(model), 5, version_msg)
            self.assertIn("parameters", model, version_msg)
            params = model["parameters"]
            self.assertEqual(len(params), SITECLASS_PARAMS[3],
                             version)
            self.assertIn("site classes", params, version_msg)
            site_classes = params["site classes"]
            self.assertEqual(len(site_classes), 4, version_msg)

    def testParseBranchSiteA(self):
        res_dir = os.path.join(self.results_dir, "codeml", "branchsiteA")
        for results_file in os.listdir(res_dir):
            version = results_file.split('-')[1].split('.')[0]
            version_msg = ("Improper parsing for version %s"
                           % version.replace('_', '.'))
            results_path = os.path.join(res_dir, results_file)
            results = codeml.read(results_path)
            # There are 5 top-level items in this case:
            # 'codon model', 'model', 'version', 'NSsites' & 'site-class model'
            self.assertEqual(len(results), 5, version_msg)
            self.assertIn("NSsites", results, version_msg)
            models = results["NSsites"]
            # Only site class model 2 is simulated for Branch Site A
            self.assertEqual(len(models), 1, version_msg)
            self.assertIn(2, models, version_msg)
            model = models[2]
            self.assertEqual(len(model), 5, version_msg)
            self.assertIn("parameters", model, version_msg)
            params = model["parameters"]
            # Branch Site A results lack a "branches" parameter
            self.assertEqual(len(params), SITECLASS_PARAMS[2] - 1, version_msg)
            self.assertIn("site classes", params, version_msg)
            site_classes = params["site classes"]
            # Branch Site A adds another site class
            self.assertEqual(len(site_classes), SITECLASSES[2] + 1,
                             version)
            for class_num in [0, 1, 2, 3]:
                self.assertIn(class_num, site_classes, version_msg)
                site_class = site_classes[class_num]
                self.assertEqual(len(site_class), 2, version_msg)
                self.assertIn("branch types", site_class, version_msg)
                branches = site_class["branch types"]
                self.assertEqual(len(branches), 2, version_msg)

    def testParseCladeModelC(self):
        cladeC_res_dir = os.path.join(self.results_dir, "codeml",
                                      "clademodelC")
        for results_file in os.listdir(cladeC_res_dir):
            version = results_file.split('-')[1].split('.')[0]
            version_msg = ("Improper parsing for version %s"
                           % version.replace('_', '.'))
            results_path = os.path.join(cladeC_res_dir, results_file)
            results = codeml.read(results_path)
            # 5 top-level items again in this case
            self.assertEqual(len(results), 5, version_msg)
            self.assertIn("NSsites", results, version_msg)
            models = results["NSsites"]
            # Only site class model 2 is simulated for Clade Model C
            self.assertEqual(len(models), 1, version_msg)
            self.assertIn(2, models, version_msg)
            model = models[2]
            self.assertEqual(len(model), 5, version_msg)
            self.assertIn("parameters", model, version_msg)
            params = model["parameters"]
            # Clade Model C results lack a "branches" parameter
            self.assertEqual(len(params), SITECLASS_PARAMS[2] - 1, version_msg)
            self.assertIn("site classes", params, version_msg)
            site_classes = params["site classes"]
            self.assertEqual(len(site_classes), SITECLASSES[2],
                             version)
            for class_num in [0, 1, 2]:
                self.assertIn(class_num, site_classes, version_msg)
                site_class = site_classes[class_num]
                self.assertEqual(len(site_class), 2, version_msg)
                self.assertIn("branch types", site_class, version_msg)
                branches = site_class["branch types"]
                self.assertEqual(len(branches), 2, version_msg)

    def testParseNgene2Mgene02(self):
        res_dir = os.path.join(self.results_dir, "codeml", "ngene2_mgene02")
        for results_file in os.listdir(res_dir):
            version = results_file.split('-')[1].split('.')[0]
            version_msg = ("Improper parsing for version %s"
                           % version.replace('_', '.'))
            results_path = os.path.join(res_dir, results_file)
            results = codeml.read(results_path)
            self.assertEqual(len(results), 4, version_msg)
            self.assertIn("NSsites", results, version_msg)
            models = results["NSsites"]
            self.assertEqual(len(models), 1, version_msg)
            self.assertIn(0, models, version_msg)
            model = models[0]
            self.assertEqual(len(model), 5, version_msg)
            self.assertIn("parameters", model, version_msg)
            params = model["parameters"]
            # This type of model has fewer parameters for model 0
            self.assertEqual(len(params), 4, version_msg)
            self.assertIn("rates", params, version_msg)
            rates = params["rates"]
            self.assertEqual(len(rates), 2, version_msg)

    def testParseNgene2Mgene1(self):
        res_dir = os.path.join(self.results_dir, "codeml", "ngene2_mgene1")
        for results_file in os.listdir(res_dir):
            version = results_file.split('-')[1].split('.')[0]
            version_msg = ("Improper parsing for version %s"
                           % version.replace('_', '.'))
            results_path = os.path.join(res_dir, results_file)
            results = codeml.read(results_path)
            self.assertEqual(len(results), 4, version_msg)
            self.assertIn("genes", results, version_msg)
            genes = results["genes"]
            self.assertEqual(len(genes), 2, version_msg)
            model = genes[0]
            self.assertEqual(len(model), 5, version_msg)
            self.assertIn("parameters", model, version_msg)
            params = model["parameters"]
            self.assertEqual(len(params), SITECLASS_PARAMS[0], version_msg)

    def testParseNgene2Mgene34(self):
        res_dir = os.path.join(self.results_dir, "codeml", "ngene2_mgene34")
        for results_file in os.listdir(res_dir):
            version = results_file.split('-')[1].split('.')[0]
            version_msg = ("Improper parsing for version %s"
                           % version.replace('_', '.'))
            results_path = os.path.join(res_dir, results_file)
            results = codeml.read(results_path)
            self.assertEqual(len(results), 4, version_msg)
            self.assertIn("NSsites", results, version_msg)
            models = results["NSsites"]
            self.assertEqual(len(models), 1, version_msg)
            self.assertIn(0, models, version_msg)
            model = models[0]
            self.assertEqual(len(model), 5, version_msg)
            self.assertIn("parameters", model, version_msg)
            params = model["parameters"]
            # This type of model has fewer parameters for model 0
            self.assertEqual(len(params), 3, version_msg)
            self.assertIn("rates", params, version_msg)
            rates = params["rates"]
            self.assertEqual(len(rates), 2, version_msg)
            self.assertIn("genes", params, version_msg)
            genes = params["genes"]
            self.assertEqual(len(genes), 2, version_msg)

    def testParseFreeRatio(self):
        res_dir = os.path.join(self.results_dir, "codeml", "freeratio")
        for results_file in os.listdir(res_dir):
            version = results_file.split('-')[1].split('.')[0]
            version_msg = ("Improper parsing for version %s"
                           % version.replace('_', '.'))
            results_path = os.path.join(res_dir, results_file)
            results = codeml.read(results_path)
            self.assertEqual(len(results), 4, version_msg)
            self.assertIn("NSsites", results, version_msg)
            models = results["NSsites"]
            self.assertEqual(len(models), 1, version_msg)
            self.assertIn(0, models, version_msg)
            model = models[0]
            # With the free ratio model, you get 3 extra trees: dN tree,
            # dS tree and omega tree
            self.assertEqual(len(model), 8, version_msg)
            self.assertIn("parameters", model, version_msg)
            params = model["parameters"]
            self.assertEqual(len(params), SITECLASS_PARAMS[0], version_msg)
            self.assertIn("branches", params, version_msg)
            # There should be 7 branches
            branches = params["branches"]
            self.assertEqual(len(branches), 7, version_msg)
            self.assertIn("omega", params, version_msg)
            omega = params["omega"]
            self.assertEqual(len(omega), 7, version_msg)

    def testParsePairwise(self):
        res_dir = os.path.join(self.results_dir, "codeml", "pairwise")
        for results_file in os.listdir(res_dir):
            version = results_file.split('-')[1].split('.')[0]
            version_msg = ("Improper parsing for version %s"
                           % version.replace('_', '.'))
            results_path = os.path.join(res_dir, results_file)
            results = codeml.read(results_path)
            # Pairwise models have an extra top-level item: pairwise
            self.assertEqual(len(results), 5, version_msg)
            self.assertIn("pairwise", results, version_msg)
            pairwise = results["pairwise"]
            self.assertGreaterEqual(len(pairwise), 2, version_msg +
                                    ": should have at least two sequences")
            for seq1, seq2 in itertools.combinations(pairwise.keys(), 2):
                self.assertEqual(len(pairwise[seq1][seq2]), 7, version_msg +
                                 ": wrong number of parameters parsed")
                self.assertEqual(len(pairwise[seq2][seq1]), 7, version_msg +
                                 ": wrong number of parameters parsed")

    def testParseSitesParamsForPairwise(self):
        """Verify that pairwise site estimates are indeed parsed. Fixes #483"""
        res_dir = os.path.join(self.results_dir, "codeml", "pairwise")
        for results_file in os.listdir(res_dir):
            version = results_file.split('-')[1].split('.')[0]
            version_msg = ("Improper parsing for version %s"
                           % version.replace('_', '.'))
            results_path = os.path.join(res_dir, results_file)
            results = codeml.read(results_path)
            self.assertIn("pairwise", results)
            seqs = list(results["pairwise"].keys())
            self.assertGreaterEqual(len(seqs), 2, version_msg +
                                    ": should have at least two sequences")
            for seq1, seq2 in itertools.combinations(seqs, 2):
                params = results["pairwise"][seq1][seq2]
                self.assertEqual(len(params), 7,
                                 version_msg + ": wrong number of parsed parameters" +
                                 " for %s-%s" % (seq1, seq2))
                for param in ("t", "S", "N", "omega", "dN", "dS", "lnL"):
                    self.assertTrue(param in params, version_msg +
                                    ": '%s' not in parsed parameters"
                                    % (param))
                    self.assertTrue(isinstance(params[param], float))
                    if param != "lnL":
                        self.assertTrue(params[param] >= 0)

    def testParseAA(self):
        res_dir = os.path.join(self.results_dir, "codeml", "aa_model0")
        for results_file in os.listdir(res_dir):
            version = results_file.split('-')[1].split('.')[0]
            version_msg = ("Improper parsing for version %s"
                           % version.replace('_', '.'))
            results_path = os.path.join(res_dir, results_file)
            results = codeml.read(results_path)
            # Amino Acid analysis has different top-levels:
            # 'NSsites', 'model', 'version', 'lnL max', 'distances'
            # Version 4.1 doesn't seem to produce distances in the results
            if version == "4_1":
                self.assertEqual(len(results), 4, version_msg)
                self.assertIn("lnL max", results, version_msg)
            else:
                self.assertEqual(len(results), 5, version_msg)
                self.assertIn("lnL max", results, version_msg)
                self.assertIn("distances", results, version_msg)
                distances = results["distances"]
                # non-pairwise AA analysis only gives raw distances
                self.assertEqual(len(distances), 1, version_msg)

    def testParseAAPairwise(self):
        res_dir = os.path.join(self.results_dir, "codeml", "aa_pairwise")
        for results_file in os.listdir(res_dir):
            version = results_file.split('-')[1].split('.')[0]
            version_msg = ("Improper parsing for version %s"
                           % version.replace('_', '.'))
            results_path = os.path.join(res_dir, results_file)
            results = codeml.read(results_path)
            # Pairwise AA analysis has one top-level fewer than non-pairwise
            self.assertEqual(len(results), 4, version_msg)
            self.assertIn("lnL max", results, version_msg)
            self.assertIn("distances", results, version_msg)
            distances = results["distances"]
            # Pairwise AA analysis has ML & raw distances
            self.assertEqual(len(distances), 2, version_msg)

    def testTreeParseVersatility(self):
        """Test finding trees in the results.

        In response to bug #453, where trees like (A, (B, C)); weren't being caught.
        """
        res_file = os.path.join(self.results_dir, "codeml",
                                "tree_regexp_versatility.out")
        results = codeml.read(res_file)
        self.assertIn("NSsites", results)
        nssites = results["NSsites"]
        self.assertIn(0, nssites)
        m0 = nssites[0]
        self.assertIn("tree", m0)
        self.assertTrue(m0["tree"] is not None)
        self.assertNotEqual(len(m0["tree"]), 0)

    def testParseM2arel(self):
        res_dir = os.path.join(self.results_dir, "codeml", "m2a_rel")
        for results_file in os.listdir(res_dir):
            version = results_file.split('-')[1].split('.')[0]
            version_msg = ("Improper parsing for version %s"
                           % version.replace('_', '.'))
            results_path = os.path.join(res_dir, results_file)
            results = codeml.read(results_path)
            self.assertIn("NSsites", results)
            self.assertIn(22, results["NSsites"])
            model = results["NSsites"][22]
            self.assertEqual(len(model), 5, version_msg)
            params = model["parameters"]
            self.assertEqual(len(params), SITECLASS_PARAMS[22],
                             version_msg)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
