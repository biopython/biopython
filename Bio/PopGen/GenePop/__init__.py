# Copyright 2007 by Tiago Antao.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
This module provides code to work with GenePop.

See http://wbiomed.curtin.edu.au/genepop/ , the format is documented
here: http://wbiomed.curtin.edu.au/genepop/help_input.html .

Classes:
Record           Holds GenePop data.

Functions:
read             Parses a GenePop record (file) into a Record object.


Partially inspired on MedLine Code.

"""
from copy import deepcopy

__docformat__ = "restructuredtext en"

def get_indiv(line):
    def int_no_zero(val):
        v = int(val)
        if v == 0:
            return None
        return v
    indiv_name, marker_line = line.split(',')
    markers = marker_line.replace('\t', ' ').split(' ')
    markers = [marker for marker in markers if marker!='']
    if len(markers[0]) in [2, 4]:  # 2 digits per allele
        marker_len = 2
    else:
        marker_len = 3
    try:
        allele_list = [(int_no_zero(marker[0:marker_len]),
                       int_no_zero(marker[marker_len:]))
                   for marker in markers]
    except ValueError:  # Haploid
        allele_list = [(int_no_zero(marker[0:marker_len]),)
                   for marker in markers]
    return indiv_name, allele_list, marker_len


def read(handle):
    """Parses a handle containing a GenePop file.

       handle is a file-like object that contains a GenePop record.
    """
    record = Record()
    record.comment_line = str(next(handle)).rstrip()
    # We can now have one loci per line or all loci in a single line
    # separated by either space or comma+space...
    # We will remove all commas on loci... that should not be a problem
    sample_loci_line = str(next(handle)).rstrip().replace(',', '')
    all_loci = sample_loci_line.split(' ')
    record.loci_list.extend(all_loci)
    for line in handle:
        line = line.rstrip()
        if line.upper()=='POP':
            break
        record.loci_list.append(line)
    else:
        raise ValueError('No population data found, file probably not GenePop related')
    record.populations.append([])
    for line in handle:
        line = line.rstrip()
        if line.upper()=='POP':
            record.populations.append([])
        else:
            indiv_name, allele_list, record.marker_len = get_indiv(line)
            record.populations[-1].append((indiv_name, allele_list))
    loci = record.loci_list
    for pop in record.populations:
        record.pop_list.append(pop[-1][0])
        for indiv in pop:
            for mk_i in range(len(loci)):
                mk_orig = indiv[1][mk_i]
                mk_real = []
                for al in mk_orig:
                    if al == 0:
                        mk_real.append(None)
                    else:
                        mk_real.append(al)
                indiv[1][mk_i] = tuple(mk_real)
    return record


class Record(object):
    """Holds information from a GenePop record.

    Members:

        - marker_len         The marker length (2 or 3 digit code per allele).

        - comment_line       Comment line.

        - loci_list          List of loci names.

        - pop_list           List of population names.

        - populations        List of population data.

    In most genepop files, the population name is not trustable.
    It is strongly recommended that populations are referred by index.

    populations has one element per population. Each element is itself
    a list of individuals, each individual is a pair composed by individual
    name and a list of alleles (2 per marker or 1 for haploids):
    Example::

        [
            [
                ('Ind1', [(1,2),    (3,3), (200,201)],
                ('Ind2', [(2,None), (3,3), (None,None)],
            ],
            [
                ('Other1', [(1,1),  (4,3), (200,200)],
            ]
        ]

    """
    def __init__(self):
        self.marker_len = 0
        self.comment_line = ""
        self.loci_list = []
        self.pop_list = []
        self.populations = []

    def __str__(self):
        """Returns (reconstructs) a GenePop textual representation.
        """
        rep = [self.comment_line + '\n']
        rep.append('\n'.join(self.loci_list) + '\n')
        for pop in self.populations:
            rep.append('Pop\n')
            for indiv in pop:
                name, markers = indiv
                rep.append(name)
                rep.append(',')
                for marker in markers:
                    rep.append(' ')
                    for al in marker:
                        if al is None:
                            al = '0'
                        aStr = str(al)
                        while len(aStr)<self.marker_len:
                            aStr = "".join(['0', aStr])
                        rep.append(aStr)
                rep.append('\n')
        return "".join(rep)

    def split_in_pops(self, pop_names):
        """Splits a GP record in a dictionary with 1 pop per entry.

            Given a record with n pops and m loci returns a dictionary
            of records (key pop_name) where each item is a record
            with a single pop and m loci.

            Parameters:
            pop_names - Population names
        """
        gp_pops = {}
        for i in range(len(self.populations)):
            gp_pop = Record()
            gp_pop.marker_len = self.marker_len
            gp_pop.comment_line = self.comment_line
            gp_pop.loci_list = deepcopy(self.loci_list)
            gp_pop.populations = [deepcopy(self.populations[i])]
            gp_pops[pop_names[i]] = gp_pop
        return gp_pops

    def split_in_loci(self, gp):
        """Splits a GP record in a dictionary with 1 locus per entry.

            Given a record with n pops and m loci returns a dictionary
            of records (key locus name) where each item is a record
            with a single locus and n pops.
        """
        gp_loci = {}
        for i in range(len(self.loci_list)):
            gp_pop = Record()
            gp_pop.marker_len = self.marker_len
            gp_pop.comment_line = self.comment_line
            gp_pop.loci_list = [self.loci_list[i]]
            gp_pop.populations = []
            for pop in self.populations:
                my_pop = []
                for indiv in pop:
                    my_pop.append((indiv[0], [indiv[1][i]]))
                gp_pop.populations.append(my_pop)
            gp_loci[gp_pop.loci_list[0]] = gp_pop
        return gp_loci

    def remove_population(self, pos):
        """Removes a population (by position).
        """
        del self.populations[pos]

    def remove_locus_by_position(self, pos):
        """Removes a locus by position.
        """
        del self.loci_list[pos]
        for pop in self.populations:
            for indiv in pop:
                name, loci = indiv
                del loci[pos]

    def remove_locus_by_name(self, name):
        """Removes a locus by name.
        """
        for i in range(len(self.loci_list)):
            if self.loci_list[i] == name:
                self.remove_locus_by_position(i)
                return
        # If here than locus not existent... Maybe raise exception?
        #   Although it should be Ok... Just a boolean return, maybe?
