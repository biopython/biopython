# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""I/O function wrappers for `Bio.Nexus` trees."""
__docformat__ = "restructuredtext en"

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

# 'index' starts from 1; 'tree' is the Newick tree string
TREE_TEMPLATE = "Tree tree%(index)d=[&U]%(tree)s;"


def parse(handle):
    """Parse the trees in a Nexus file.

    Uses the old Nexus.Trees parser to extract the trees, converts them back to
    plain Newick trees, and feeds those strings through the new Newick parser.
    This way we don't have to modify the Nexus module yet. (Perhaps we'll
    eventually change Nexus to use the new NewickIO parser directly.)
    """
    nex = Nexus.Nexus(handle)
    # NB: Once Nexus.Trees is modified to use Tree.Newick objects, do this:
    # return iter(nex.trees)
    # Until then, convert the Nexus.Trees.Tree object hierarchy:
    def node2clade(nxtree, node):
        subclades = [node2clade(nxtree, nxtree.node(n)) for n in node.succ]
        return Newick.Clade(
                branch_length=node.data.branchlength,
                name=node.data.taxon,
                clades=subclades,
                confidence=node.data.support,
                comment=node.data.comment)

    for nxtree in nex.trees:
        newroot = node2clade(nxtree, nxtree.node(nxtree.root))
        yield Newick.Tree(root=newroot, rooted=nxtree.rooted, name=nxtree.name,
                          weight=nxtree.weight)

def write(obj, handle, **kwargs):
    """Write a new Nexus file containing the given trees.

    Uses a simple Nexus template and the NewickIO writer to serialize just the
    trees and minimal supporting info needed for a valid Nexus file.
    """
    trees = list(obj)
    writer = NewickIO.Writer(trees)
    nexus_trees = [TREE_TEMPLATE % {'index': idx+1, 'tree': nwk}
                   for idx, nwk in enumerate(
                        writer.to_strings(plain=False, plain_newick=True,
                                          **kwargs))]
    tax_labels = map(str, chain(*(t.get_terminals() for t in trees)))
    text = NEX_TEMPLATE % {
            'count':    len(tax_labels),
            'labels':   ' '.join(tax_labels),
            'trees':    '\n'.join(nexus_trees),
            }
    handle.write(text)
    return len(nexus_trees)
