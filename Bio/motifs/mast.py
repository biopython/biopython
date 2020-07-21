# Copyright 2008 by Bartek Wilczynski.
# Adapted from Bio.MEME.Parser by Jason A. Hackney.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Module for the support of Motif Alignment and Search Tool (MAST)."""

import xml.etree.ElementTree as ET

from Bio.motifs import meme


class Record(list):
    """The class for holding the results from a MAST run.

    A mast.Record holds data about matches between motifs and sequences.
    The motifs held by the Record are objects of the class meme.Motif.

    The mast.Record class inherits from list, so you can access individual
    motifs in the record by their index. Alternatively, you can find a motif
    by its name:

    >>> from Bio import motifs
    >>> with open("motifs/mast.crp0.de.oops.txt.xml") as f:
    ...     record = motifs.parse(f, 'MAST')
    >>> motif = record[0]
    >>> print(motif.name)
    1
    >>> motif = record['1']
    >>> print(motif.name)
    1
    """

    def __init__(self):
        """Initialize the class."""
        self.sequences = []
        self.version = ""
        self.database = ""
        self.diagrams = {}
        self.alphabet = None
        self.strand_handling = ""

    def __getitem__(self, key):
        """Return the motif of index key."""
        if isinstance(key, str):
            for motif in self:
                if motif.name == key:
                    return motif
        else:
            return list.__getitem__(self, key)


def read(handle):
    """Parse a MAST XML format handle as a Record object."""
    record = Record()
    try:
        xml_tree = ET.parse(handle)
    except ET.ParseError:
        raise ValueError(
            "Improper MAST XML input file. XML root tag should start with <mast version= ..."
        )
    __read_metadata(record, xml_tree)
    __read_sequences(record, xml_tree)
    return record


# Everything below is private


def __read_metadata(record, xml_tree):
    record.version = xml_tree.getroot().get("version")
    record.database = xml_tree.find("sequence_dbs").find("sequence_db").get("source")
    record.alphabet = xml_tree.find("alphabet").get("name")
    record.strand_handling = xml_tree.find("settings").get("strand_handling")
    # TODO - read other metadata
    for i, motif_tree in enumerate(xml_tree.find("motifs").findall("motif")):
        motif = meme.Motif(record.alphabet)
        # TODO - motif.name not in XML - always index?
        motif.name = str(i + 1)
        motif.id = motif_tree.get("id")
        motif.alt_id = motif_tree.get("alt")
        motif.length = int(motif_tree.get("length"))
        # TODO - add nsites, evalue
        record.append(motif)


def __read_sequences(record, xml_tree):
    """Read sequences from XML ElementTree object."""
    for sequence_tree in xml_tree.find("sequences").findall("sequence"):
        sequence_name = sequence_tree.get("name")
        record.sequences.append(sequence_name)
        diagram_str = __make_diagram(record, sequence_tree)
        record.diagrams[sequence_name] = diagram_str
        # TODO - add description, evalue, length, combined_pvalue


def __make_diagram(record, sequence_tree):
    """Make diagram string found in text file based on motif hit info."""
    sequence_length = int(sequence_tree.get("length"))
    hit_eles, hit_motifs, gaps = [], [], []
    for seg_tree in sequence_tree.findall("seg"):
        for hit_ele in seg_tree.findall("hit"):
            hit_pos = int(hit_ele.get("pos"))
            if not hit_eles:
                gap = hit_pos - 1
            else:
                gap = hit_pos - int(hit_eles[-1].get("pos")) - hit_motifs[-1].length
            gaps.append(gap)
            hit_motifs.append(record[int(hit_ele.get("idx"))])
            hit_eles.append(hit_ele)
    if not hit_eles:
        return str(sequence_length)
    if record.strand_handling == "combine":
        motif_strs = [
            f"[{'-' if hit_ele.get('rc') == 'y' else '+'}{hit_motif.name}]"
            for hit_ele, hit_motif in zip(hit_eles, hit_motifs)
        ]
    elif record.strand_handling == "unstranded":
        motif_strs = [
            f"[{hit_motif.name}]" for hit_ele, hit_motif in zip(hit_eles, hit_motifs)
        ]
    else:
        # TODO - more strand_handling possibilities?
        raise Exception(f"Strand handling option {record.strand_handling} not parsable")
    tail_length = (
        sequence_length - int(hit_eles[-1].get("pos")) - hit_motifs[-1].length + 1
    )
    motifs_with_gaps = [str(s) for pair in zip(gaps, motif_strs) for s in pair] + [
        str(tail_length)
    ]
    # remove 0-length gaps
    motifs_with_gaps = [s for s in motifs_with_gaps if s != "0"]
    return "-".join(motifs_with_gaps)
