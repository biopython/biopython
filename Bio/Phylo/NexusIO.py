# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""I/O function wrappers for ``Bio.Nexus`` trees."""

from itertools import chain

from Bio.Nexus import Nexus
from Bio.Phylo import Newick, NewickIO


# Structure of a Nexus tree-only file
NEX_TEMPLATE = """\
#NEXUS
Begin Taxa;
 Dimensions NTax=%(count)d;
 TaxLabels %(labels)s;
End;
Begin Trees;
 %(trees)s
End;
"""


def parse(handle):
    """Parse the trees in a Nexus file."""
    nex = Nexus.Nexus(handle)

    # Nexus.Trees has been modified to use Tree.Newick objects
    return iter(nex.trees)


def write(obj, handle, **kwargs):
    """Write a new Nexus file containing the given trees.

    Uses a simple Nexus template and the NewickIO writer to serialize just the
    trees and minimal supporting info needed for a valid Nexus file.
    """
    trees = list(obj)
    writer = NewickIO.Writer(trees)

    nexus_trees = list(writer.to_strings(plain=False, plain_newick=False, **kwargs))
    tax_labels = [str(x.name) for x in chain(*(t.get_terminals() for t in trees))]
    text = NEX_TEMPLATE % {
        "count": len(tax_labels),
        "labels": " ".join(tax_labels),
        "trees": "\n".join(nexus_trees),
    }
    handle.write(text)
    return len(nexus_trees)
