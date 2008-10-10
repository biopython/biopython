# Copyright 2007 by Tiago Antao <tiagoantao@gmail.com>.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


import os
from Bio.PopGen import GenePop
import Bio.PopGen.FDist

# Quite a few utility functions could be done (like remove pop,
# add locus, etc...). The recommended strategy is convert back
# and forth from/to GenePop and use GenePop Utils

def convert_genepop_to_fdist(gp_rec):
    """Converts a GenePop record to a FDist one.

       Parameters:
       gp_rec - Genepop Record

       Returns:
       FDist record.
    """
    fd_rec = Bio.PopGen.FDist.Record()
    
    fd_rec.data_org = 0
    fd_rec.num_loci = len(gp_rec.loci_list)
    fd_rec.num_pops = len(gp_rec.populations)
    for lc_i in range(len(gp_rec.loci_list)):
        alleles = []
        pop_data = []
        for pop_i in range(len(gp_rec.populations)):
            for indiv in gp_rec.populations[pop_i]:
                for al in indiv[1][lc_i]:
                    if al is not None and al not in alleles:
                        alleles.append(al)
        #here we go again (necessary...)
        for pop_i in range(len(gp_rec.populations)):
            allele_counts = {}
            for indiv in gp_rec.populations[pop_i]:
                for al in indiv[1][lc_i]:
                    if al is not None:
                        count = allele_counts.get(al, 0)
                        allele_counts[al] = count + 1
            allele_array = [] #We need the same order as in alleles
            for allele in alleles:
                allele_array.append(allele_counts.get(allele, 0))
            pop_data.append(allele_array)
            #if lc_i==3:
            #    print alleles, allele_counts#, pop_data
        fd_rec.loci_data.append((len(alleles), pop_data))
    return fd_rec

def approximate_fst(desired_fst, simulated_fst, parameter_fst,
           max_run_fst = 1, min_run_fst = 0, limit = 0.005):
    """Calculates the next Fst attempt in order to approximate a
       desired Fst.
    
    """
    if abs(simulated_fst - desired_fst) < limit:
        return parameter_fst, max_run_fst, min_run_fst
    if simulated_fst > desired_fst:
        max_run_fst = parameter_fst
        next_parameter_fst = (min_run_fst + parameter_fst)/2
    else:
        min_run_fst = parameter_fst
        next_parameter_fst = (max_run_fst + parameter_fst)/2
    return next_parameter_fst, max_run_fst, min_run_fst

