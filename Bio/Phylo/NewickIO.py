# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# Based on Bio.Nexus, copyright 2005-2008 by Frank Kauff & Cymon J. Cox.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""I/O function wrappers for the Newick file format.

See: http://evolution.genetics.washington.edu/phylip/newick_doc.html
"""
__docformat__ = "restructuredtext en"

from cStringIO import StringIO

from Bio.Phylo import Newick

# Definitions retrieved from Bio.Nexus.Trees
NODECOMMENT_START = '[&'
NODECOMMENT_END = ']'


class NewickError(Exception):
    """Exception raised when Newick object construction cannot continue."""
    pass


# ---------------------------------------------------------
# Public API

def parse(handle, **kwargs):
    """Iterate over the trees in a Newick file handle.

    :returns: generator of Bio.Phylo.Newick.Tree objects.
    """
    return Parser(handle).parse(**kwargs)


def write(trees, handle, plain=False, **kwargs):
    """Write a trees in Newick format to the given file handle.

    :returns: number of trees written.
    """
    return Writer(trees).write(handle, plain=plain, **kwargs)


# ---------------------------------------------------------
# Input

class Parser(object):
    """Parse a Newick tree given a file handle.

    Based on the parser in `Bio.Nexus.Trees`.
    """

    def __init__(self, handle):
        self.handle = handle

    @classmethod
    def from_string(cls, treetext):
        handle = StringIO(treetext)
        return cls(handle)

    def parse(self, values_are_confidence=False, comments_are_confidence=True, rooted=False):
        """Parse the text stream this object was initialized with."""
        self.values_are_confidence = values_are_confidence
        self.comments_are_confidence = comments_are_confidence
        self.rooted = rooted
        buf = ''
        for line in self.handle:
            buf += line.rstrip()
            if buf.endswith(';'):
                yield self._parse_tree(buf)
                buf = ''
        if buf:
            # Last tree is missing a terminal ';' character -- that's OK
            yield self._parse_tree(buf)

    def _parse_tree(self, text):
        """Parses the text representation into an Tree object."""
        text = text.strip()

        def make_new_clade(parent=None):
            clade = Newick.Clade(name='')
            if parent: clade.parent = parent
            return clade
            
        root_clade = make_new_clade()
        
        def process_clade(clade):
            if not clade.name:
                clade.name = None
            if hasattr(clade, 'branch_length_string'):
                if self.values_are_confidence:
                    clade.confidence = float(clade.branch_length_string)
                else:
                    clade.branch_length = float(clade.branch_length_string)
                del clade.branch_length_string
            if self.comments_are_confidence and hasattr(clade, 'comment') and clade.comment:
                try:
                    clade.confidence = int(clade.comment) / 100.
                except ValueError, TypeError:
                    try:
                        clade.confidence = float(clade.comment)
                    except ValueError, TypeError:
                        pass
                if hasattr(clade, 'confidence'):
                    try:
                        assert 0 <= clade.confidence <= 1
                    except AssertionError:
                        del clade.confidence
            if hasattr(clade, 'parent'):
                parent = clade.parent
                parent.clades.append(clade)
                del clade.parent
                return parent
            if clade is root_clade:
                return root_clade
                
        current_clade = root_clade
        entering_quoted_string = False
        escaped = False
        entering_comment = False
        entering_branch_length = False
        
        for char in text:
            if entering_quoted_string:
                # add characters to clade name until closing quote is found
                if char == "'" and not escaped:
                    entering_quoted_string = False
                elif char == '\\':
                    escaped = True
                else:
                    current_clade.name += char
                    
            elif entering_comment:
                # add characters to comment until closing square bracket is found
                if char == ']':
                    entering_comment = False
                else:
                    current_clade.comment += char
                    
            else:
                if char == '(':
                    # start a new clade, which is a child of the current clade
                    current_clade = make_new_clade(current_clade)
                    entering_branch_length = False
                    
                elif char == ',':
                    # if the current clade is the root, then the external parentheses are missing
                    # and a new root should be created
                    if current_clade is root_clade:
                        root_clade = make_new_clade()
                        current_clade.parent = root_clade

                    # start a new child clade at the same level as the current clade
                    parent = process_clade(current_clade)
                    current_clade = make_new_clade(parent)
                    entering_branch_length = False
                    
                elif char == ')':
                    # done adding children for this parent clade
                    parent = process_clade(current_clade)
                    if not parent:
                        raise NewickError('Parenthesis mismatch.')
                    current_clade = parent
                    entering_branch_length = False
                    
                elif char == '[':
                    # start entering comment
                    current_clade.comment = ''
                    entering_comment = True
                    
                elif char == ':':
                    # start entering branch length
                    entering_branch_length = True
                    current_clade.branch_length_string = ''
                    
                elif char == "'":
                    # start entering quoted label
                    entering_quoted_string = True
                    
                elif char == ';': pass
                    
                elif entering_branch_length:
                    # add characters to branch length
                    current_clade.branch_length_string += char
                    
                else:
                    # add characters to node label
                    current_clade.name += char
            
            escape = False
            
        process_clade(current_clade)
        process_clade(root_clade)
        
        return Newick.Tree(root=root_clade, rooted=self.rooted)


# ---------------------------------------------------------
# Output

class Writer(object):
    """Based on the writer in Bio.Nexus.Trees (str, to_string)."""

    def __init__(self, trees):
        self.trees = trees

    def write(self, handle, **kwargs):
        """Write this instance's trees to a file handle."""
        count = 0
        for treestr in self.to_strings(**kwargs):
            handle.write(treestr + '\n')
            count += 1
        return count

    def to_strings(self, confidence_as_branch_length=False,
            branch_length_only=False, plain=False,
            plain_newick=True, ladderize=None, max_confidence=1.0,
            format_confidence='%1.2f', format_branch_length='%1.5f'):
        """Return an iterable of PAUP-compatible tree lines."""
        # If there's a conflict in the arguments, we override plain=True
        if confidence_as_branch_length or branch_length_only:
            plain = False
        make_info_string = self._info_factory(plain,
                confidence_as_branch_length, branch_length_only, max_confidence,
                format_confidence, format_branch_length)

        def newickize(clade):
            """Convert a node tree to a Newick tree string, recursively."""
            if clade.is_terminal():    # terminal
                return ((clade.name or '')
                        + make_info_string(clade, terminal=True))
            else:
                subtrees = (newickize(sub) for sub in clade)
                return '(%s)%s' % (','.join(subtrees),
                        (clade.name or '') + make_info_string(clade))

        # Convert each tree to a string
        for tree in self.trees:
            if ladderize in ('left', 'LEFT', 'right', 'RIGHT'):
                # Nexus compatibility shim, kind of
                tree.ladderize(reverse=(ladderize in ('right', 'RIGHT')))
            rawtree = newickize(tree.root) + ';'
            if plain_newick:
                yield rawtree
                continue
            # Nexus-style (?) notation before the raw Newick tree
            treeline = ['tree', (tree.name or 'a_tree'), '=']
            if tree.weight != 1:
                treeline.append('[&W%s]' % round(float(tree.weight), 3))
            if tree.rooted:
                treeline.append('[&R]')
            treeline.append(rawtree)
            yield ' '.join(treeline)

    def _info_factory(self, plain, confidence_as_branch_length,
            branch_length_only, max_confidence, format_confidence,
            format_branch_length):
        """Return a function that creates a nicely formatted node tag."""
        if plain:
            # Plain tree only. That's easy.
            def make_info_string(clade, terminal=False):
                return ''

        elif confidence_as_branch_length:
            # Support as branchlengths (eg. PAUP), ignore actual branchlengths
            def make_info_string(clade, terminal=False):
                if terminal:
                    # terminal branches have 100% support
                    return ':' + format_confidence % max_confidence
                else:
                    return ':' + format_confidence % clade.confidence

        elif branch_length_only:
            # write only branchlengths, ignore support
            def make_info_string(clade, terminal=False):
                return ':' + format_branch_length % clade.branch_length

        else:
            # write support and branchlengths (e.g. .con tree of mrbayes)
            def make_info_string(clade, terminal=False):
                if (terminal or
                        not hasattr(clade, 'confidence') or
                        clade.confidence is None):
                    return (':' + format_branch_length
                            ) % (clade.branch_length or 0.0)
                else:
                    return (format_confidence + ':' + format_branch_length
                            ) % (clade.confidence, clade.branch_length or 0.0)

        return make_info_string
