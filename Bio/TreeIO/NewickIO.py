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
    buf = ''
    for line in file:
        buf += line.rstrip()
        if buf.endswith(';'):
            yield Trees.Tree(tree=buf)
            buf = ''
    if buf:
        # Last tree is missing a terminal ';' character
        yield Trees.Tree(tree=buf)


def write(trees, file, plain=False, **kwargs):
    count = 0
    for tree in trees:
        file.write(tree.to_string(plain_newick=True, plain=plain, **kwargs)
                    + ';\n')
        count += 1
    return count
