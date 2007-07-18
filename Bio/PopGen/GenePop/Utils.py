# Copyright 2007 by Tiago Antao.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


from Bio.PopGen import GenePop
from copy import deepcopy

"""
Utility functions to deal with GenePop files

Functions:
split_in_pops            - Splits a GP file in a list of 1 per pop (with all loci)
split_in_loci            - Splits a GP file in a list of 1 per locus (with all pops)
remove_population        - Removes a population (by position)
remove_locus_by_name     - Removes locus by name
remove_locus_by_position - Removes locus by position

"""
def split_in_pops(gp, pop_names):
    gp_pops = {}
    for i in range(len(gp.populations)):
        gp_pop = GenePop.Record()
        gp_pop.marker_len = gp.marker_len
        gp_pop.comment_line = gp.comment_line
        gp_pop.loci_list = deepcopy(gp.loci_list)
        gp_pop.populations = [deepcopy(gp.populations[i])]
        gp_pops[pop_names[i]] = gp_pop
    return gp_pops

def split_in_loci(gp):
    gp_loci = {}
    for i in range(len(gp.loci_list)):
        gp_pop = GenePop.Record()
        gp_pop.marker_len = gp.marker_len
        gp_pop.comment_line = gp.comment_line
        gp_pop.loci_list = [gp.loci_list[i]]
        gp_pop.populations = []
        for pop in gp.populations:
            my_pop = []
            for indiv in pop:
                my_pop.append((indiv[0], [indiv[1][i]]))
            gp_pop.populations.append(my_pop)
        gp_loci[gp_pop.loci_list[0]] = gp_pop
    return gp_loci


def remove_population(gp_rec, pos):
    """Removes a population (by position).
    """
    del gp_rec.populations[pos]
    
def remove_locus_by_position(gp_rec, pos):
    """Removes a locus by position.
    """
    del gp_rec.loci_list[pos]
    for pop in gp_rec.populations:
        for indiv in pop:
            name, loci = indiv
            del loci[pos]

def remove_locus_by_name(gp_rec, name):
    """Removes a locus by name.
    """
    for i in range(len(gp_rec.loci_list)):
        if gp_rec.loci_list[i] == name:
            remove_locus_by_position(gp_rec, i)
            return
    #If here than locus not existent... Maybe raise exception?
    #   Although it should be Ok... Just a boolean return, maybe?
