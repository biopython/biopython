# Copyright 2013 by Leighton Pritchard.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Classes and functions to parse a KGML pathway map.

The KGML pathway map is parsed into the object structure defined in
KGML_Pathway.py in this module.

Classes:
 - KGMLParser - Parses KGML file

Functions:
 - read - Returns a single Pathway object, using KGMLParser internally

"""

from __future__ import print_function

try:
    import xml.etree.cElementTree as ElementTree
except ImportError:
    import xml.etree.ElementTree as ElementTree

from Bio._py3k import StringIO

from Bio.KEGG.KGML.KGML_pathway import Component, Entry, Graphics
from Bio.KEGG.KGML.KGML_pathway import Pathway, Reaction, Relation


def read(handle, debug=0):
    """Parse a single KEGG Pathway from given file handle.

    Returns a single Pathway object.  There should be one and only
    one pathway in each file, but there may well be pathological
    examples out there.
    """
    iterator = parse(handle, debug)
    try:
        first = next(iterator)
    except StopIteration:
        first = None
    if first is None:
        raise ValueError("No pathways found in handle")
    try:
        second = next(iterator)
    except StopIteration:
        second = None
    if second is not None:
        raise ValueError("More than one pathway found in handle")
    return first


def parse(handle, debug=0):
    """Return an iterator over Pathway elements.

    Arguments:
     - handle - file handle to a KGML file for parsing
     - debug - integer for amount of debug information to print

    This is a generator for the return of multiple Pathway objects.

    """
    # Check handle
    if not hasattr(handle, "read"):
        if isinstance(handle, str):
            handle = StringIO(handle)
        else:
            exc_txt = "An XML-containing handle or an XML string must be provided"
            raise Exception(exc_txt)
    # Parse XML and return each Pathway
    for event, elem in ElementTree.iterparse(handle, events=("start", "end")):
        if event == "end" and elem.tag == "pathway":
            yield KGMLParser(elem).parse()
            elem.clear()


class KGMLParser(object):
    """Parses a KGML XML Pathway entry into a Pathway object.

    Example: Read and parse large metabolism file

    >>> from Bio.KEGG.KGML.KGML_parser import read
    >>> pathway = read(open('KEGG/ko01100.xml', 'r'))
    >>> print(len(pathway.entries))
    3628
    >>> print(len(pathway.reactions))
    1672
    >>> print(len(pathway.maps))
    149

    >>> pathway = read(open('KEGG/ko00010.xml', 'r'))
    >>> print(pathway) #doctest: +NORMALIZE_WHITESPACE
    Pathway: Glycolysis / Gluconeogenesis
    KEGG ID: path:ko00010
    Image file: http://www.kegg.jp/kegg/pathway/ko/ko00010.png
    Organism: ko
    Entries: 99
    Entry types:
        ortholog: 61
        compound: 31
        map: 7

    """

    def __init__(self, elem):
        """Initialize the class."""
        self.entry = elem

    def parse(self):
        """Parse the input elements."""
        # This comment stops black style adding a blank line here, which causes flake8 D202.
        def _parse_pathway(attrib):
            for k, v in attrib.items():
                self.pathway.__setattr__(k, v)

        def _parse_entry(element):
            new_entry = Entry()
            for k, v in element.attrib.items():
                new_entry.__setattr__(k, v)
            for subelement in element.getchildren():
                if subelement.tag == "graphics":
                    _parse_graphics(subelement, new_entry)
                elif subelement.tag == "component":
                    _parse_component(subelement, new_entry)
            self.pathway.add_entry(new_entry)

        def _parse_graphics(element, entry):
            new_graphics = Graphics(entry)
            for k, v in element.attrib.items():
                new_graphics.__setattr__(k, v)
            entry.add_graphics(new_graphics)

        def _parse_component(element, entry):
            new_component = Component(entry)
            for k, v in element.attrib.items():
                new_component.__setattr__(k, v)
            entry.add_component(new_component)

        def _parse_reaction(element):
            new_reaction = Reaction()
            for k, v in element.attrib.items():
                new_reaction.__setattr__(k, v)
            for subelement in element.getchildren():
                if subelement.tag == "substrate":
                    new_reaction.add_substrate(int(subelement.attrib["id"]))
                elif subelement.tag == "product":
                    new_reaction.add_product(int(subelement.attrib["id"]))
            self.pathway.add_reaction(new_reaction)

        def _parse_relation(element):
            new_relation = Relation()
            new_relation.entry1 = int(element.attrib["entry1"])
            new_relation.entry2 = int(element.attrib["entry2"])
            new_relation.type = element.attrib["type"]
            for subtype in element.getchildren():
                name, value = subtype.attrib["name"], subtype.attrib["value"]
                if name in ("compound", "hidden compound"):
                    new_relation.subtypes.append((name, int(value)))
                else:
                    new_relation.subtypes.append((name, value))
            self.pathway.add_relation(new_relation)

        # ==========
        # Initialize Pathway
        self.pathway = Pathway()
        # Get information about the pathway itself
        _parse_pathway(self.entry.attrib)
        for element in self.entry.getchildren():
            if element.tag == "entry":
                _parse_entry(element)
            elif element.tag == "reaction":
                _parse_reaction(element)
            elif element.tag == "relation":
                _parse_relation(element)
            # Parsing of some elements not implemented - no examples yet
            else:
                # This should warn us of any unimplemented tags
                import warnings
                from Bio import BiopythonParserWarning

                warnings.warn(
                    "Warning: tag %s not implemented in parser" % element.tag,
                    BiopythonParserWarning,
                )
        return self.pathway


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest(verbose=0)
