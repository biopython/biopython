# Copyright 2007 by Tiago Antao <tiagoantao@gmail.com>.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""This module provides code to work with FDist (DEPRECATED).

See http://www.rubic.rdg.ac.uk/~mab/software.html (old) and
http://www.maths.bris.ac.uk/~mamab/ (new) for downloading the
FDist tool by Mark Beaumont.

Classes:
Record           Holds FDist data.

Functions:
read             Parses a FDist record (file) into a Record object.


"""

import warnings

from Bio import BiopythonDeprecationWarning

warnings.warn("Bio.PopGen.FDist has been deprecated, and we intend to remove"
              " it in a future release of Biopython. If you would like to"
              " continue using it, please contact the Biopython developers"
              " via the mailing list.", BiopythonDeprecationWarning)


def read(handle):
    """Parses FDist data into a Record object.

       handle is a file-like object that contains a FDist record.
    """
    record = Record()
    record.data_org = int(str(next(handle)).rstrip())
    record.num_pops = int(str(next(handle)).rstrip())
    record.num_loci = int(str(next(handle)).rstrip())
    for i in range(record.num_loci):
        next(handle)
        num_alleles = int(str(next(handle)).rstrip())
        pops_data = []
        if record.data_org == 0:
            for j in range(record.num_pops):
                line_comp = str(next(handle)).rstrip().split(' ')
                pop_dist = [int(x) for x in line_comp]
                pops_data.append(pop_dist)
        else:
            raise NotImplementedError('1/alleles by rows not implemented')
        record.loci_data.append((num_alleles, pops_data))
    return record


class Record(object):
    """Holds information from a FDist record.

    Members:

        - data_org    Data organization (0 pops by rows, 1 alleles by rows).
          The Record will behave as if data was 0 (converting if needed)

        - num_pops       Number of populations

        - num_loci       Number of loci

        - loci_data      Loci data

    loci_data is a list, where each element represents a locus. Each element
    is a tuple, the first element is the number of alleles, the second
    element a list. Each element of the list is the count of each allele
    per population.
    """
    def __init__(self):
        self.data_org = 0
        self.num_pops = 0
        self.num_loci = 0
        self.loci_data = []

    def __str__(self):
        rep = ['0\n']  # We only export in 0 format, even if originally was 1
        rep.append(str(self.num_pops) + '\n')
        rep.append(str(self.num_loci) + '\n')
        rep.append('\n')
        for locus_data in self.loci_data:
            num_alleles, pops_data = locus_data
            rep.append(str(num_alleles) + '\n')
            for pop_data in pops_data:
                for allele_count in pop_data:
                    rep.append(str(allele_count) + ' ')
                rep.append('\n')
            rep.append('\n')
        return "".join(rep)
