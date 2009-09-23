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


class AppTest(unittest.TestCase):
    """Tests genepop execution via biopython using EasyController.
    """

    def setUp(self):
        #Genepop likes to be on the directory where the file is.
        os.chdir("PopGen")
        self.ctrl = EasyController("big.gen")

    def tearDown(self):
        os.chdir("..")

    def test_basic_info(self):
        """Test basic info.
        """
        pops, loci = self.ctrl.get_basic_info()
        assert len(pops) == 10
        assert len(loci) == 37

    def test_get_heterozygosity_info(self):
        """Test heterozygosity info.
        """
        hz_info = self.ctrl.get_heterozygosity_info(0, "Locus2")
        assert hz_info[1] == 24
        assert hz_info[3] == 7

    def test_get_alleles(self):
        """Test get alleles.
        """
        assert self.ctrl.get_alleles(0,"Locus3") == [3, 20]

    def test_get_alleles_all_pops(self):
        """Test get alleles for all populations.
        """
        assert self.ctrl.get_alleles_all_pops("Locus4") == [1, 3]

    def test_get_fis(self):
        """Test get Fis.
        """
        alleles, overall = self.ctrl.get_fis(0,"Locus2")
        assert alleles[3][0] == 55
        assert overall[0] == 62

    def test_get_allele_frequency(self):
        """Test allele frequency.
        """
        tot_genes, alleles = self.ctrl.get_allele_frequency(0,"Locus2")
        assert tot_genes == 62
        assert abs(alleles[20] - 0.113) < 0.05

    def test_get_genotype_count(self):
        """Test genotype count.
        """
        assert len(self.ctrl.get_genotype_count(0,"Locus2")) == 3

    def test_estimate_nm(self):
        """Test Nm estimation.
        """
        nms = self.ctrl.estimate_nm()
        assert nms[0] == 28.0

    def test_get_avg_fst_pair_locus(self):
        """Test get average Fst for pairwise pops on a locus.
        """
        assert len(self.ctrl.get_avg_fst_pair_locus("Locus4")) == 45

    def test_get_avg_fst_pair(self):
        """Test get pairwise Fst.
        """
        pop_fis =  self.ctrl.get_avg_fst_pair()
        assert len(pop_fis) == 45

    def test_get_avg_fis(self):
        """Test average Fis.
        """
        self.ctrl.get_avg_fis()

    def test_get_multilocus_f_stats(self):
        """Test multilocus F stats.
        """
        mf = self.ctrl.get_multilocus_f_stats()
        assert len(mf) == 3
        assert mf[0]<0.1

    def test_get_f_stats(self):
        """Test F stats.
        """
        fs = self.ctrl.get_f_stats("Locus2")
        assert len(fs)==5
        assert fs[0]<0

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
