# Copyright 2008 by Bartek Wilczynski
# Adapted from  Bio.MEME.Parser by Jason A. Hackney.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Module for the support of MEME motif format."""

from __future__ import print_function
import xml.etree.ElementTree as ET

from Bio import Seq
from Bio import motifs


def read(handle):
    """Parse the text output of the MEME program into a meme.Record object.

    Examples
    --------
    >>> from Bio.motifs import meme
    >>> with open("motifs/meme.out") as f:
    ...     record = meme.read(f)
    >>> for motif in record:
    ...     for instance in motif.instances:
    ...         print(instance.motif_name, instance.sequence_name, instance.strand, instance.pvalue)
    Motif 1 SEQ10; + 8.71e-07
    Motif 1 SEQ9; + 8.71e-07
    Motif 1 SEQ8; + 8.71e-07
    Motif 1 SEQ7; + 8.71e-07
    Motif 1 SEQ6; + 8.71e-07
    Motif 1 SEQ5; + 8.71e-07
    Motif 1 SEQ4; + 8.71e-07
    Motif 1 SEQ3; + 8.71e-07
    Motif 1 SEQ2; + 8.71e-07
    Motif 1 SEQ1; + 8.71e-07

    """
    record = Record()
    try:
        xml_tree = ET.parse(handle)
    except ET.ParseError:
        raise ValueError("Improper MAST XML input file. XML root tag should start with <mast version= ...")
    __read_metadata(record, xml_tree)
    __read_alphabet(record, xml_tree)
    __read_sequences(record, xml_tree)
    __read_motifs(record, xml_tree)
    return record


class Motif(motifs.Motif):
    """A subclass of Motif used in parsing MEME (and MAST) output.

    This subclass defines functions and data specific to MEME motifs.
    This includes the motif name, the evalue for a motif, and its number
    of occurrences.
    """

    def __init__(self, alphabet=None, instances=None):
        """Initialize the class."""
        motifs.Motif.__init__(self, alphabet, instances)
        self.evalue = 0.0
        self.num_occurrences = 0
        self.name = None
        self.id = None
        self.alt_id = None


class Instance(Seq.Seq):
    """A class describing the instances of a MEME motif, and the data thereof."""

    def __init__(self, *args, **kwds):
        """Initialize the class."""
        Seq.Seq.__init__(self, *args, **kwds)
        self.sequence_name = ""
        self.start = 0
        self.pvalue = 1.0
        self.strand = 0
        self.length = 0
        self.motif_name = ""


class Record(list):
    """A class for holding the results of a MEME run.

    A meme.Record is an object that holds the results from running
    MEME. It implements no methods of its own.

    The meme.Record class inherits from list, so you can access individual
    motifs in the record by their index. Alternatively, you can find a motif
    by its name:

    >>> from Bio import motifs
    >>> with open("motifs/meme.out") as f:
    ...     record = motifs.parse(f, 'MEME')
    >>> motif = record[0]
    >>> print(motif.name)
    Motif 1
    >>> motif = record['Motif 1']
    >>> print(motif.name)
    Motif 1
    """

    def __init__(self):
        """Initialize."""
        self.version = ""
        self.datafile = ""
        self.command = ""
        self.alphabet = ""
        self.sequences = []

    def __getitem__(self, key):
        """Return the motif of index key."""
        if isinstance(key, str):
            for motif in self:
                if motif.name == key:
                    return motif
        else:
            return list.__getitem__(self, key)


# Everything below is private


def __read_metadata(record, xml_tree):
    record.version = xml_tree.getroot().get("version")
    record.datafile = xml_tree.find("training_set").get("primary_sequences")
    record.command = xml_tree.find("model").find("command_line").text
    # TODO - background_frequencies, other metadata under model


def __read_alphabet(record, xml_tree):
    alphabet_tree = xml_tree.find("training_set").find("letter_frequencies").find("alphabet_array")
    for value in alphabet_tree.findall("value"):
        record.alphabet += value.get("letter_id")


def __read_sequences(record, xml_tree):
    for sequence_tree in xml_tree.find("training_set").findall("sequence"):
        sequence_name = sequence_tree.get('name')
        record.sequences.append(sequence_name)
        # TODO - sequence id, length, weight


def __read_motifs(record, xml_tree):
    for motif_tree in xml_tree.find("motifs").findall("motif"):
        instances = []
        for site_tree in motif_tree.find("contributing_sites").findall("contributing_site"):
            letters = [letter_ref.get("letter_id") for letter_ref in site_tree.find("site").findall("letter_ref")]
            sequence = ''.join(letters)
            instance = Instance(sequence, record.alphabet)
            instance.motif_name = motif_tree.get("name")
            instance.sequence_name = site_tree.get("sequence_id")  # TODO - rename to sequence_id, get sequence_name
            # TODO - left flank, right flank
            instance.start = int(site_tree.get("position")) + 1
            instance.pvalue = float(site_tree.get("pvalue"))
            instance.strand = __read_strand(site_tree.get("strand"))
            instance.length = len(sequence)
            instances.append(instance)
        instances = motifs.Instances(instances, record.alphabet)
        motif = Motif(record.alphabet, instances)
        motif.id = motif_tree.get("id")
        motif.name = motif_tree.get("name")
        motif.alt_id = motif_tree.get("alt")
        motif.length = int(motif_tree.get("width"))
        motif.num_occurrences = int(motif_tree.get("sites"))
        motif.evalue = float(motif_tree.get("e_value"))
        # TODO - ic, re, llr, pvalue, bayes_threshold, elapsed_time
        record.append(motif)


def __read_strand(strand):
    """Convert strand (+/-) from XML if present. Default: +"""
    if (strand == "minus"):
        return "-"
    if (strand == "plus" or strand == "none"):
        return "+"
