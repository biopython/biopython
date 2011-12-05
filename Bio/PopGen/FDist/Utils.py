# Copyright 2007 by Tiago Antao <tiagoantao@gmail.com>.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


import os
from Bio.PopGen import GenePop
from Bio.PopGen.GenePop import FileParser
import Bio.PopGen.FDist

# Quite a few utility functions could be done (like remove pop,
# add locus, etc...). The recommended strategy is convert back
# and forth from/to GenePop and use GenePop Utils

def convert_genepop_to_fdist(gp_rec, report_pops = None):
    """Converts a GenePop record to a FDist one.

       Parameters:
       gp_rec - Genepop Record (either standard or big)

       Returns:
       FDist record.
    """
    if hasattr(gp_rec, "populations"):
        return _convert_genepop_to_fdist(gp_rec)
    else:
        return _convert_genepop_to_fdist_big(gp_rec, report_pops)

def _convert_genepop_to_fdist(gp_rec):
    """Converts a standard GenePop record to a FDist one.

       Parameters:
       gp_rec - Genepop Record (Standard)

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
        alleles.sort() #Dominance requires this

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

def _convert_genepop_to_fdist_big(gp_rec, report_pops = None):
    """Converts a big GenePop record to a FDist one.

       Parameters:
       gp_rec - Genepop Record (Big)

       Returns:
       FDist record.
    """
    fd_rec = Bio.PopGen.FDist.Record()

    fd_rec.data_org = 1
    fd_rec.num_loci = len(gp_rec.loci_list)
    num_loci = len(gp_rec.loci_list)
    loci = []
    for i in range(num_loci):
        loci.append(set())
    pops = []
    work_rec = FileParser.read(gp_rec.fname)
    lParser = work_rec.get_individual()
    def init_pop():
        my_pop = []
        for i in range(num_loci):
            my_pop.append({})
        return my_pop
    curr_pop = init_pop()
    num_pops = 1
    if report_pops:
        report_pops(num_pops)
    while lParser:
        if lParser != True:
            for loci_pos in range(num_loci):
                for al in lParser[1][loci_pos]:
                    if al is not None:
                        loci[loci_pos].add(al)
                        curr_pop[loci_pos][al]= curr_pop[loci_pos].get(al,0)+1
        else:
            pops.append(curr_pop)
            num_pops += 1
            if report_pops:
                report_pops(num_pops)
            curr_pop = init_pop()
        lParser = work_rec.get_individual()
    work_rec._handle.close() #TODO - Needs a proper fix
    pops.append(curr_pop)
    fd_rec.num_pops = num_pops
    for loci_pos in range(num_loci):
        alleles = list(loci[loci_pos])
        alleles.sort()
        loci_rec = [len(alleles), []]
        for pop in pops:
            pop_rec = []
            for allele in alleles:
                pop_rec.append(pop[loci_pos].get(allele, 0))
            loci_rec[1].append(pop_rec)
        fd_rec.loci_data.append(tuple(loci_rec))
    return fd_rec


def _convert_genepop_to_fdist_big_old(gp_rec, report_loci = None):
    """Converts a big GenePop record to a FDist one.

       Parameters:
       gp_rec - Genepop Record (Big)

       Returns:
       FDist record.
    """
    fd_rec = Bio.PopGen.FDist.Record()

    def countPops(rec):
        f2 = FileParser.read(rec.fname)
        popCnt = 1
        while f2.skip_population():
            popCnt += 1
        return popCnt

    
    fd_rec.data_org = 0
    fd_rec.num_loci = len(gp_rec.loci_list)
    work_rec0 = FileParser.read(gp_rec.fname)
    fd_rec.num_pops = countPops(work_rec0)

    num_loci = len(gp_rec.loci_list)
    for lc_i in range(num_loci):
        if report_loci:
            report_loci(lc_i, num_loci)
        work_rec = FileParser.read(gp_rec.fname)
        work_rec2 = FileParser.read(gp_rec.fname)

        alleles = []
        pop_data = []
        lParser = work_rec.get_individual()
        while lParser:
            if lParser != True:
                for al in lParser[1][lc_i]:
                    if al is not None and al not in alleles:
                        alleles.append(al)
            lParser = work_rec.get_individual()
        #here we go again (necessary...)
        alleles.sort()
        def process_pop(pop_data, alleles, allele_counts):
            allele_array = [] #We need the same order as in alleles
            for allele in alleles:
                allele_array.append(allele_counts.get(allele, 0))
            pop_data.append(allele_array)
        lParser = work_rec2.get_individual()
        allele_counts = {}
        for allele in alleles:
            allele_counts[allele] = 0
        allele_counts[None]=0
        while lParser:
            if lParser == True:
                process_pop(pop_data, alleles, allele_counts)
                allele_counts = {}
                for allele in alleles:
                    allele_counts[allele] = 0
                allele_counts[None]=0
            else:
                for al in lParser[1][lc_i]:
                    allele_counts[al] += 1
            lParser = work_rec2.get_individual()
        process_pop(pop_data, alleles, allele_counts)
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

