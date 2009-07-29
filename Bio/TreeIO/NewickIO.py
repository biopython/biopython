# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""I/O function wrappers for the Newick file format.

See: U{ http://evolution.genetics.washington.edu/phylip/newick_doc.html }
"""
__docformat__ = "epytext en"

from Bio.Nexus import Trees

def read(file):
    try:
        tree_gen = parse(file)
        tree = tree_gen.next()
    except StopIteration:
        raise RuntimeError("There are no trees in this file.")
    try:
        tree_gen.next()
    except StopIteration:
        return tree
    else:
        raise RuntimeError(
                "There are multiple trees in this file; use parse() instead.")


def parse(file):
    do_close = False
    if not hasattr(file, 'read'):
        file = open(file, 'r')
        do_close = True
    try:
        buf = None
        for line in file:
            buf = (buf and (buf + line) or line).strip()
            if buf.endswith(';'):
                yield Trees.Tree(tree=buf)
                buf = None
    finally:
        if do_close:
            file.close()


def write(trees, file, plain=False, **kwargs):
    do_close = False
    if not hasattr(file, 'write'):
        file = open(file, 'w+')
        do_close = True
    try:
        lines = (t.to_string(plain_newick=True, plain=plain, **kwargs)
                 for t in trees)
        file.write(';\n'.join(lines))
    finally:
        if do_close:
            file.close()
