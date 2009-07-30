# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""I/O function wrappers for the Newick file format.

See: U{ http://evolution.genetics.washington.edu/phylip/newick_doc.html }
"""
__docformat__ = "epytext en"

from Bio.Nexus import Trees


def parse(file):
    do_close = False
    if not hasattr(file, 'read'):
        file = open(file, 'r')
        do_close = True
    # Python 2.4 support: This should be in a try block, but yield is not OK
    buf = ''
    for line in file:
        buf += line.rstrip()
        if buf.endswith(';'):
            yield Trees.Tree(tree=buf)
            buf = ''
    if buf:
        # Last tree is missing a terminal ';' character
        yield Trees.Tree(tree=buf)
    # /py2.4
    if do_close:
        file.close()


def write(trees, file, plain=False, **kwargs):
    do_close = False
    if not hasattr(file, 'write'):
        file = open(file, 'w+')
        do_close = True
    try:
        for tree in trees:
            file.write(tree.to_string(plain_newick=True, plain=plain, **kwargs)
                        + ';\n')
    finally:
        if do_close:
            file.close()
