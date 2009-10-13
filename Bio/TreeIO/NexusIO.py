# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""I/O function wrappers for Bio.Nexus trees.
"""
__docformat__ = "epytext en"

from itertools import chain

from Bio.Nexus import Nexus

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


def parse(file):
    nex = Nexus.Nexus(file)
    return iter(nex.trees)


# XXX Nexus/Newick-specific tree attributes: to_string, get_taxa
def write(obj, file, **kwargs):
    trees = list(obj)
    nexus_trees = [TREE_TEMPLATE % {
            'index': idx+1,
            'tree': tree.to_string(plain_newick=True, plain=False, **kwargs)
            } for idx, tree in enumerate(trees)]
    tax_labels = list(chain(*(t.get_taxa() for t in trees)))
    text = NEX_TEMPLATE % {
            'count':    len(tax_labels),
            'labels':   ' '.join(tax_labels),
            'trees':    '\n'.join(nexus_trees),
            }
    file.write(text)
    return len(nexus_trees)
