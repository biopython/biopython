# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Objects instantiated from elements in a parsed PhyloXML file.

"""

try:
    from xml.etree import cElementTree as ElementTree
except ImportError:
    try:
        from xml.etree import ElementTree as ElementTree
    except ImportError:
        # Python 2.4 -- check for 3rd-party implementations
        try:
            from lxml.etree import ElementTree
        except ImportError:
            try:
                import cElementTree as ElementTree
            except ImportError:
                try:
                    from elementtree import ElementTree
                except ImportError:
                    from Bio import MissingExternalDependencyError
                    raise MissingExternalDependencyError(
                            "No ElementTree module was found. " \
                            "Use Python 2.5+, lxml or elementtree if you " \
                            "want to use Bio.PhyloXML.")


def _dump_tags(source):
    """Extract tags from an XML document and print them to standard output.
    
    This function is meant for testing and debugging only.
    """
    events = ('start', 'end')
    for event, elem in ElementTree.iterparse(source, events=events):
        if event == 'start':
            print elem.tag
        else:
            elem.clear()


def parse(source):
    """Parse a phyloXML file or stream and build a tree of Biopython objects.
    """
    pass


# ---------------------------------------------------------------------
# Classes instantiated from phyloXML nodes

class PhyloElement(object):
    """Base class for all PhyloXML objects."""
    def __init__(self, attrib, text, children):
        self._attrib = attrib
        self._text = text
        self._children = children
        # Munge each of these into properties


class phyloxml(PhyloElement):
    """Root node of the PhyloXML document."""
    pass


# Tier 1:
class branch_length(PhyloElement):
    pass

class clade(PhyloElement):
    pass

class code(PhyloElement):
    pass

class confidence(PhyloElement):
    pass

class name(PhyloElement):
    pass

class phylogeny(PhyloElement):
    pass

class taxonomy(PhyloElement):
    pass


# Tier 2:
class domain(PhyloElement):
    pass

class domain_architecture(PhyloElement):
    pass

class duplications(PhyloElement):
    pass

class events(PhyloElement):
    pass

class scientific_name(PhyloElement):
    pass

class sequence(PhyloElement):
    pass

class speciations(PhyloElement):
    pass


