# Copyright 2009 by Tiago Antao <tiagoantao@gmail.com>.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
import unittest

from Bio import MissingExternalDependencyError
from Bio.PopGen.GenePop.EasyController import EasyController

#Tests genepop related code for easy contorller. Note: this requires genepop
#test_PopGen_GenePop_nodepend tests code that does not require genepop

found = False
for path in os.environ['PATH'].split(os.pathsep):
    try:
        for filename in os.listdir(path):
            if filename.startswith('Genepop'):
                found = True
    except os.error:
        pass #Path doesn't exist - correct to pass
if not found:
    raise MissingExternalDependencyError(\
        "Install GenePop if you want to use Bio.PopGen.GenePop.")


cur_dir = os.path.abspath(".") #Tests directory

class AppTest(unittest.TestCase):
    """Tests genepop execution via biopython using EasyController.
    """

    def setUp(self):
        #Genepop likes to be on the directory where the file is.
        os.chdir("PopGen")
        self.ctrl = EasyController("big.gen")

    def tearDown(self):
        os.chdir(cur_dir)

    def test_basic_info(self):
        """Test basic info.
        """
        pops, loci = self.ctrl.get_basic_info()
        self.assertEqual(len(pops), 10)
        self.assertEqual(len(loci), 37)

    def test_get_heterozygosity_info(self):
        """Test heterozygosity info.
        """
        hz_info = self.ctrl.get_heterozygosity_info(0, "Locus2")
        self.assertEqual(hz_info[1], 24)
        self.assertEqual(hz_info[3], 7)

    def test_get_alleles(self):
        """Test get alleles.
        """
        #Returns keys of a dict, so order is Python implementation dependent
        self.assertEqual(set(self.ctrl.get_alleles(0,"Locus3")), set([3, 20]))

    def test_get_alleles_all_pops(self):
        """Test get alleles for all populations.
        """
        self.assertEqual(self.ctrl.get_alleles_all_pops("Locus4"), [1, 3])

    def test_get_fis(self):
        """Test get Fis.
        """
        alleles, overall = self.ctrl.get_fis(0,"Locus2")
        self.assertEqual(alleles[3][0], 55)
        self.assertEqual(overall[0], 62)

    def test_get_allele_frequency(self):
        """Test allele frequency.
        """
        tot_genes, alleles = self.ctrl.get_allele_frequency(0,"Locus2")
        self.assertEqual(tot_genes, 62)
        self.assertTrue(abs(alleles[20] - 0.113) < 0.05)

    def test_get_genotype_count(self):
        """Test genotype count.
        """
        self.assertEqual(len(self.ctrl.get_genotype_count(0,"Locus2")), 3)

    def test_estimate_nm(self):
        """Test Nm estimation.
        """
        nms = self.ctrl.estimate_nm()
        self.assertEqual(nms[0], 28.0)

#These tests are frequently failing, possibly due to a Genepop problem.
#    def test_get_avg_fst_pair_locus(self):
#        """Test get average Fst for pairwise pops on a locus.
#        """
#        self.assertEqual(len(self.ctrl.get_avg_fst_pair_locus("Locus4")), 45)
#
#    def test_get_avg_fst_pair(self):
#        """Test get pairwise Fst.
#        """
#        pop_fis =  self.ctrl.get_avg_fst_pair()
#        self.assertEqual(len(pop_fis), 45)

    def test_get_avg_fis(self):
        """Test average Fis.
        """
        self.ctrl.get_avg_fis()

    def test_get_multilocus_f_stats(self):
        """Test multilocus F stats.
        """
        mf = self.ctrl.get_multilocus_f_stats()
        self.assertEqual(len(mf), 3)
        self.assertTrue(mf[0]<0.1)

    def test_get_f_stats(self):
        """Test F stats.
        """
        fs = self.ctrl.get_f_stats("Locus2")
        self.assertEqual(len(fs), 5)
        self.assertTrue(fs[0]<0)

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
