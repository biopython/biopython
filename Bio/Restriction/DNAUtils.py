# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""DNA utilities for Bio.Restriction (DEPRECATED).

DNAUtils was written in C and therefore would not be available on Jython etc.
It offered three string based functions:
 - complement, duplicating the functionality of the Seq object
 - antiparallel, duplicating the functionality of the Seq object and
   the reverse_complement function in Bio.Seq
 - check_bases, a very odd validation routine unlikely to be of general use.
"""

import warnings
warnings.warn("Bio.Restriction.DNAUtils is deprecated, and will be "
              "removed in a future release of Biopython.")
del warnings

#expose these existing functions mimicking the old DNAUtils names:
from Bio.Seq import reverse_complement as antiparallel

#quick and dirty complement function, maybe we should add one to Bio.Seq?
def complement(seq_string) :
    return antiparallel(seq_string)[::-1]

#expose this re-implementation of the old C code function check_bases:
from Bio.Restriction.Restriction import _check_bases as check_bases
