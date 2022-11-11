# Copyright 2015 by Gert Hulselmans.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Parse XMS motif files."""

from Bio import motifs


from xml.dom import minidom, Node
import re


class XMSScanner:
    """Class for scanning XMS XML file."""

    def __init__(self, doc):
        """Generate motif Record from xms document, an XML-like motif pfm file."""
        self.record = Record()
        for child in doc.getElementsByTagName("motif"):
            if child.nodeType == Node.ELEMENT_NODE:
                self.handle_motif(child)

    def handle_motif(self, node):
        """Read the motif's name and column from the node and add the motif record."""
        motif_name = self.get_text(node.getElementsByTagName("name"))
        nucleotide_counts = {"A": [], "C": [], "G": [], "T": []}

        for column in node.getElementsByTagName("column"):
            [
                nucleotide_counts[nucleotide].append(float(nucleotide_count))
                for nucleotide, nucleotide_count in zip(
                    ["A", "C", "G", "T"], self.get_acgt(column)
                )
            ]

        motif = motifs.Motif(alphabet="GATC", counts=nucleotide_counts)
        motif.name = motif_name

        self.record.append(motif)

    def get_property_value(self, node, key_name):
        """Extract the value of the motif's property named key_name from node."""
        for cur_property in node.getElementsByTagName("prop"):
            right_property = False
            cur_value = None
            for child in cur_property.childNodes:
                if child.nodeType != Node.ELEMENT_NODE:
                    continue
                if child.tagName == "key" and self.get_text([child]) == key_name:
                    right_property = True
                if child.tagName == "value":
                    cur_value = self.get_text([child])
            if right_property:
                return cur_value
        return None

    def get_acgt(self, node):
        """Get and return the motif's weights of A, C, G, T."""
        a, c, g, t = 0.0, 0.0, 0.0, 0.0
        for weight in node.getElementsByTagName("weight"):
            if weight.getAttribute("symbol") == "adenine":
                a = float(self.get_text([weight]))
            elif weight.getAttribute("symbol") == "cytosine":
                c = float(self.get_text([weight]))
            elif weight.getAttribute("symbol") == "guanine":
                g = float(self.get_text([weight]))
            elif weight.getAttribute("symbol") == "thymine":
                t = float(self.get_text([weight]))
        return a, c, g, t

    def get_text(self, nodelist):
        """Return a string representation of the motif's properties listed on nodelist ."""
        retlist = []
        for node in nodelist:
            if node.nodeType == Node.TEXT_NODE:
                retlist.append(node.wholeText)
            elif node.hasChildNodes:
                retlist.append(self.get_text(node.childNodes))

        return re.sub(r"\s+", " ", "".join(retlist))


class Record(list):
    """Class to store the information in a XMS matrix table.

    The record inherits from a list containing the individual motifs.
    """

    def __str__(self):
        """Return a string representation of the motifs in the Record object."""
        return "\n".join(str(motif) for motif in self)


def read(handle):
    """Read motifs in XMS matrix format from a file handle.

    XMS is an XML format for describing regulatory motifs and PSSMs.
    This format was defined by Thomas Down, and used in the NestedMICA and MotifExplorer programs.
    """
    xms_doc = minidom.parse(handle)
    record = XMSScanner(xms_doc).record

    return record
