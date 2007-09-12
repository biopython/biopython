# Copyright 2001 by Tarjei Mikkelsen. All rights reserved.
# Copyright 2007 by Michiel de Hoon. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
This module provides code to import KEGG Pathway maps for use with
the Biopython Pathway module.

The pathway maps are in the format:

RXXXXX:[X.X.X.X:] A + 2 B <=> C
RXXXXX:[X.X.X.X:] 3C <=> 2 D + E
...

where RXXXXX is a five-digit reaction id, and X.X.X.X is the optional
EC number of the enzyme that catalyze the reaction.


Classes:
Iterator             -- Iterates through a file of map file.

"""

from Bio.Pathway import Reaction

class Iterator:
    """Iterator to read a file of KEGG reactions one at a time.
    """
    def __init__(self, handle, parser = None):
        """Initialize the iterator.

        Arguments:
        o handle - A handle to a map file to iterate through.
        """
        import warnings
        warnings.warn("Bio.KEGG.Map.Iterator(handle, parser) is deprecated. Please use Bio.KEGG.Map.parse(handle) instead. It also returns an iterator.", DeprecationWarning)
        self.records = parse(handle)


    def next(self):
        """Return the next Pathway.Reaction object from the handle.

        Will return None if we ran out of records.
        """
        return self.records.next()
    
    def __iter__(self):
        return iter(self.next, None)


def parse(handle):
    for line in handle:
        data, catalysts, reaction = line.split(":")
        catalysts = [(catalysts,)]
        reactants = {}
        before, after = reaction.split("<=>")
        compounds = before.split(" + ")
        for compound in compounds:
            compound = compound.strip()
            try:
               number, compound = compound.split()
               number = -int(number)
            except ValueError:
               number = -1
            reactants[compound] = number
        compounds = after.split(" + ")
        for compound in compounds:
            compound = compound.strip()
            try:
               number, compound = compound.split()
               number = int(number)
            except ValueError:
               number = +1
            reactants[compound] = number
        yield Reaction(reactants, catalysts, True, data)
