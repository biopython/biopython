# Copyright 2009 by Tiago Antao <tiagoantao@gmail.com>.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Control GenePop through an easier interface.

This interface is less efficient than the standard GenePopControler

"""

from .Controller import GenePopController
from Bio.PopGen import GenePop


class EasyController:
    """Define a class for an easier interface with the GenePop program."""

    def __init__(self, fname, genepop_dir=None):
        """Initialize the controller.

        genepop_dir is the directory where GenePop is.

        The binary should be called Genepop (capital G)
        """
        self._fname = fname
        self._controller = GenePopController(genepop_dir)
        self.__fst_pair_locus = {}  # More caches like this needed!
        self.__allele_frequency = {}  # More caches like this needed!

    def get_basic_info(self):
        """Obtain the population list and loci list from the file."""
        with open(self._fname) as f:
            rec = GenePop.read(f)
        return rec.pop_list, rec.loci_list

    # 1.3
    def test_hw_pop(self, pop_pos, test_type="probability"):
        """Perform Hardy-Weinberg test on the given position."""
        if test_type == "deficiency":
            hw_res = self._controller.test_pop_hz_deficiency(self._fname)
        elif test_type == "excess":
            hw_res = self._controller.test_pop_hz_excess(self._fname)
        else:
            loci_res, hw_res, fisher_full = self._controller.test_pop_hz_prob(
                self._fname, ".P"
            )
        for i in range(pop_pos - 1):
            next(hw_res)
        return next(hw_res)

    # 1.4
    def test_hw_global(
        self,
        test_type="deficiency",
        enum_test=True,
        dememorization=10000,
        batches=20,
        iterations=5000,
    ):
        """Perform Hardy-Weinberg global Heterozygote test."""
        if test_type == "deficiency":
            pop_res, loc_res, all = self._controller.test_global_hz_deficiency(
                self._fname, enum_test, dememorization, batches, iterations
            )
        else:
            pop_res, loc_res, all = self._controller.test_global_hz_excess(
                self._fname, enum_test, dememorization, batches, iterations
            )
        return list(pop_res), list(loc_res), all

    # 2.1
    def test_ld_all_pair(
        self, locus1, locus2, dememorization=10000, batches=20, iterations=5000
    ):
        """Test for linkage disequilibrium for each pair of loci in each population."""
        all_ld = self._controller.test_ld(
            self._fname, dememorization, batches, iterations
        )[1]
        for ld_case in all_ld:
            (l1, l2), result = ld_case
            if (l1 == locus1 and l2 == locus2) or (l1 == locus2 and l2 == locus1):
                return result

    def estimate_nm(self):
        """Estimate Nm. Just a simple bridge."""
        return self._controller.estimate_nm(self._fname)

    def get_heterozygosity_info(self, pop_pos, locus_name):
        """Return the heterozygosity info for a certain locus on a population.

        Returns (Expected homozygotes, observed homozygotes,
        Expected heterozygotes, observed heterozygotes)
        """
        geno_freqs = self._controller.calc_allele_genotype_freqs(self._fname)
        pop_iter, loc_iter = geno_freqs
        pops = list(pop_iter)
        return pops[pop_pos][1][locus_name][1]

    def get_genotype_count(self, pop_pos, locus_name):
        """Return the genotype counts for a certain population and locus."""
        geno_freqs = self._controller.calc_allele_genotype_freqs(self._fname)
        pop_iter, loc_iter = geno_freqs
        pop_iter = list(pop_iter)
        return pop_iter[pop_pos][1][locus_name][0]

    def get_fis(self, pop_pos, locus_name):
        """Return the Fis for a certain population and locus.

        Below CW means Cockerham and Weir and RH means Robertson and Hill.

        Returns a pair:

        - dictionary [allele] = (repetition count, frequency, Fis CW )
          with information for each allele
        - a triple with total number of alleles, Fis CW, Fis RH

        """
        geno_freqs = self._controller.calc_allele_genotype_freqs(self._fname)
        pop_iter, loc_iter = geno_freqs
        pops = list(pop_iter)
        return pops[pop_pos][1][locus_name][2:]

    def get_alleles(self, pop_pos, locus_name):
        """Return the alleles for a certain population and locus."""
        geno_freqs = self._controller.calc_allele_genotype_freqs(self._fname)
        pop_iter, loc_iter = geno_freqs
        pop_iter = list(pop_iter)
        return list(pop_iter[pop_pos][1][locus_name][2].keys())

    def get_alleles_all_pops(self, locus_name):
        """Return the alleles for a certain population and locus."""
        geno_freqs = self._controller.calc_allele_genotype_freqs(self._fname)
        pop_iter, loc_iter = geno_freqs
        for locus_info in loc_iter:
            if locus_info[0] == locus_name:
                return locus_info[1]

    def get_allele_frequency(self, pop_pos, locus_name):
        """Calculate the allele frequency for a certain locus on a population."""
        if len(self.__allele_frequency) == 0:
            geno_freqs = self._controller.calc_allele_genotype_freqs(self._fname)
            pop_iter, loc_iter = geno_freqs
            for locus_info in loc_iter:
                if locus_info[0] is None:
                    self.__allele_frequency[locus_info[0]] = None, None
                else:
                    self.__allele_frequency[locus_info[0]] = locus_info[1:]
        info = self.__allele_frequency[locus_name]
        pop_name, freqs, total = info[1][pop_pos]
        allele_freq = {}
        alleles = info[0]
        for i, allele in enumerate(alleles):
            allele_freq[allele] = freqs[i]
        return total, allele_freq

    def get_multilocus_f_stats(self):
        """Return the multilocus F stats.

        Explain averaging.
        Returns Fis(CW), Fst, Fit
        """
        return self._controller.calc_fst_all(self._fname)[0]

    def get_f_stats(self, locus_name):
        """Return F stats for a locus.

        Returns Fis(CW), Fst, Fit, Qintra, Qinter
        """
        loci_iter = self._controller.calc_fst_all(self._fname)[1]
        for name, fis, fst, fit, qintra, qinter in loci_iter:
            if name == locus_name:
                return fis, fst, fit, qintra, qinter

    def get_avg_fis(self):
        """Calculate identity-base average Fis."""
        return self._controller.calc_diversities_fis_with_identity(self._fname)[1]

    def get_avg_fst_pair(self):
        """Calculate Allele size-base average Fis for all population pairs."""
        return self._controller.calc_fst_pair(self._fname)[1]

    def get_avg_fst_pair_locus(self, locus):
        """Calculate Allele size-base average Fis for all population pairs of the given locus."""
        if len(self.__fst_pair_locus) == 0:
            iter = self._controller.calc_fst_pair(self._fname)[0]
            for locus_info in iter:
                self.__fst_pair_locus[locus_info[0]] = locus_info[1]
        return self.__fst_pair_locus[locus]

    def calc_ibd(self, is_diplo=True, stat="a", scale="Log", min_dist=0.00001):
        """Calculate isolation by distance statistics for Diploid or Haploid."""
        if is_diplo:
            return self._controller.calc_ibd_diplo(self._fname, stat, scale, min_dist)
        else:
            return self._controller.calc_ibd_haplo(self._fname, stat, scale, min_dist)
