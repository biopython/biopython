# Copyright 2009 by Tiago Antao <tiagoantao@gmail.com>.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
import unittest

from Bio import MissingExternalDependencyError
from Bio.PopGen.GenePop.Controller import GenePopController

#Tests genepop related code. Note: this case requires genepop
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
    """Tests genepop execution via biopython.
    """

    def test_allele_genotype_frequencies(self):
        """Test genepop execution on basic allele and genotype frequencies.
        """
        ctrl = GenePopController()
        pop_iter, locus_iter = ctrl.calc_allele_genotype_freqs("PopGen" + os.sep + "big.gen")
        #print pop, loci
        #for popc in pop_iter:
        #    pop_name, loci_content = popc
        #    print pop_name
        #    for locus in loci_content.keys():
        #        geno_list, hets, freq_fis = loci_content[locus]
        #        print locus
        #        print hets
        #        print freq_fis
        #        print geno_list
        #    print

    def test_calc_diversities_fis_with_identity(self):
        """Test calculations of diversities ...
        """
        ctrl = GenePopController()
        iter, avg_fis, avg_Qintra = ctrl.calc_diversities_fis_with_identity(
            "PopGen" + os.sep + "big.gen")
        liter = list(iter)
        assert len(liter) == 37
        assert liter[0][0] == "Locus1"
        assert len(avg_fis)==10
        assert len(avg_Qintra)==10

    def test_estimate_nm(self):
        """Test Nm estimation.
        """
        ctrl = GenePopController()
        mean_sample_size, mean_priv_alleles, mig10, mig25, mig50, mig_corrected =\
            ctrl.estimate_nm("PopGen" + os.sep + "big.gen")
        assert (mean_sample_size, mean_priv_alleles, mig10, mig25, mig50, mig_corrected) == \
               (28.0, 0.016129, 52.5578, 15.3006, 8.94583, 13.6612)


    def test_fst_all(self):
        """Test genepop execution on all fst.
        """
        ctrl = GenePopController()
        (allFis, allFst, allFit), itr = ctrl.calc_fst_all("PopGen" + os.sep + "c2line.gen")
        results = list(itr)
        assert (len(results) == 3)
        assert (results[0][0] == "136255903")
        assert (results[1][3] - 0.33 < 0.01)

    def test_haploidy(self):
        """Test haploidy.
        """
        ctrl = GenePopController()
        (allFis, allFst, allFit), itr = ctrl.calc_fst_all("PopGen" + os.sep + "haplo.gen")
        litr = list(itr)
        assert not type(allFst) == int
        assert len(litr) == 37
        assert litr[36][0] == "Locus37" 

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
