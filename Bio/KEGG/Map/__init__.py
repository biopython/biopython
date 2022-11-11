# Copyright 2001 by Tarjei Mikkelsen. All rights reserved.
# Copyright 2007 by Michiel de Hoon. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Load KEGG Pathway maps for use with the Biopython Pathway module.

The pathway maps are in the format::

    RXXXXX:[X.X.X.X:] A + 2 B <=> C
    RXXXXX:[X.X.X.X:] 3C <=> 2 D + E
    ...

where RXXXXX is a five-digit reaction id, and X.X.X.X is the optional
EC number of the enzyme that catalyze the reaction.
"""

from Bio.Pathway import Reaction


def parse(handle):
    """Parse a KEGG pathway map."""
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
