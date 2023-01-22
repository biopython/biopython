# Copyright 2017-2019 Damien Goutte-Gattat.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.SeqIO support for the SnapGene file format.

The SnapGene binary format is the native format used by the SnapGene program
from GSL Biotech LLC.
"""
from datetime import datetime
from re import sub
from struct import unpack
from xml.dom.minidom import parseString

from Bio.Seq import Seq
from Bio.SeqFeature import SimpleLocation
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from .Interfaces import SequenceIterator


def _iterate(handle):
    """Iterate over the packets of a SnapGene file.

    A SnapGene file is made of packets, each packet being a TLV-like
    structure comprising:

      - 1 single byte indicating the packet's type;
      - 1 big-endian long integer (4 bytes) indicating the length of the
        packet's data;
      - the actual data.
    """
    while True:
        packet_type = handle.read(1)
        if len(packet_type) < 1:  # No more packet
            return
        packet_type = unpack(">B", packet_type)[0]

        length = handle.read(4)
        if len(length) < 4:
            raise ValueError("Unexpected end of packet")
        length = unpack(">I", length)[0]

        data = handle.read(length)
        if len(data) < length:
            raise ValueError("Unexpected end of packet")

        yield (packet_type, length, data)


def _parse_dna_packet(length, data, record):
    """Parse a DNA sequence packet.

    A DNA sequence packet contains a single byte flag followed by the
    sequence itself.
    """
    if record.seq:
        raise ValueError("The file contains more than one DNA packet")

    flags, sequence = unpack(">B%ds" % (length - 1), data)
    record.seq = Seq(sequence.decode("ASCII"))
    record.annotations["molecule_type"] = "DNA"
    if flags & 0x01:
        record.annotations["topology"] = "circular"
    else:
        record.annotations["topology"] = "linear"


def _parse_notes_packet(length, data, record):
    """Parse a 'Notes' packet.

    This type of packet contains some metadata about the sequence. They
    are stored as a XML string with a 'Notes' root node.
    """
    xml = parseString(data.decode("UTF-8"))
    type = _get_child_value(xml, "Type")
    if type == "Synthetic":
        record.annotations["data_file_division"] = "SYN"
    else:
        record.annotations["data_file_division"] = "UNC"

    date = _get_child_value(xml, "LastModified")
    if date:
        record.annotations["date"] = datetime.strptime(date, "%Y.%m.%d")

    acc = _get_child_value(xml, "AccessionNumber")
    if acc:
        record.id = acc

    comment = _get_child_value(xml, "Comments")
    if comment:
        record.name = comment.split(" ", 1)[0]
        record.description = comment
        if not acc:
            record.id = record.name


def _parse_cookie_packet(length, data, record):
    """Parse a SnapGene cookie packet.

    Every SnapGene file starts with a packet of this type. It acts as
    a magic cookie identifying the file as a SnapGene file.
    """
    cookie, seq_type, exp_version, imp_version = unpack(">8sHHH", data)
    if cookie.decode("ASCII") != "SnapGene":
        raise ValueError("The file is not a valid SnapGene file")


def _parse_location(rangespec, strand, record, is_primer=False):
    start, end = (int(x) for x in rangespec.split("-"))
    # Account for SnapGene's 1-based coordinates
    start = start - 1
    if is_primer:
        # Primers' coordinates in SnapGene files are shifted by -1
        # for some reasons
        start += 1
        end += 1
    if start > end:
        # Range wrapping the end of the sequence
        l1 = SimpleLocation(start, len(record), strand=strand)
        l2 = SimpleLocation(0, end, strand=strand)
        location = l1 + l2
    else:
        location = SimpleLocation(start, end, strand=strand)
    return location


def _parse_features_packet(length, data, record):
    """Parse a sequence features packet.

    This packet stores sequence features (except primer binding sites,
    which are in a dedicated Primers packet). The data is a XML string
    starting with a 'Features' root node.
    """
    xml = parseString(data.decode("UTF-8"))
    for feature in xml.getElementsByTagName("Feature"):
        quals = {}

        type = _get_attribute_value(feature, "type", default="misc_feature")

        strand = +1
        directionality = int(
            _get_attribute_value(feature, "directionality", default="1")
        )
        if directionality == 2:
            strand = -1

        location = None
        subparts = []
        n_parts = 0
        for segment in feature.getElementsByTagName("Segment"):
            if _get_attribute_value(segment, "type", "standard") == "gap":
                continue
            rng = _get_attribute_value(segment, "range")
            n_parts += 1
            next_location = _parse_location(rng, strand, record)
            if not location:
                location = next_location
            elif strand == -1:
                # Reverse segments order for reverse-strand features
                location = next_location + location
            else:
                location = location + next_location

            name = _get_attribute_value(segment, "name")
            if name:
                subparts.append([n_parts, name])

        if len(subparts) > 0:
            # Add a "parts" qualifiers to represent "named subfeatures"
            if strand == -1:
                # Reverse segment indexes and order for reverse-strand features
                subparts = reversed([[n_parts - i + 1, name] for i, name in subparts])
            quals["parts"] = [";".join(f"{i}:{name}" for i, name in subparts)]

        if not location:
            raise ValueError("Missing feature location")

        for qualifier in feature.getElementsByTagName("Q"):
            qname = _get_attribute_value(
                qualifier, "name", error="Missing qualifier name"
            )
            qvalues = []
            for value in qualifier.getElementsByTagName("V"):
                if value.hasAttribute("text"):
                    qvalues.append(_decode(value.attributes["text"].value))
                elif value.hasAttribute("predef"):
                    qvalues.append(_decode(value.attributes["predef"].value))
                elif value.hasAttribute("int"):
                    qvalues.append(int(value.attributes["int"].value))
            quals[qname] = qvalues

        name = _get_attribute_value(feature, "name")
        if name:
            if "label" not in quals:
                # No explicit label attribute, use the SnapGene name
                quals["label"] = [name]
            elif name not in quals["label"]:
                # The SnapGene name is different from the label,
                # add a specific attribute to represent it
                quals["name"] = [name]

        feature = SeqFeature(location, type=type, qualifiers=quals)
        record.features.append(feature)


def _parse_primers_packet(length, data, record):
    """Parse a Primers packet.

    A Primers packet is similar to a Features packet but specifically
    stores primer binding features. The data is a XML string starting
    with a 'Primers' root node.
    """
    xml = parseString(data.decode("UTF-8"))
    for primer in xml.getElementsByTagName("Primer"):
        quals = {}

        name = _get_attribute_value(primer, "name")
        if name:
            quals["label"] = [name]

        locations = []
        for site in primer.getElementsByTagName("BindingSite"):
            rng = _get_attribute_value(
                site, "location", error="Missing binding site location"
            )
            strand = int(_get_attribute_value(site, "boundStrand", default="0"))
            if strand == 1:
                strand = -1
            else:
                strand = +1

            location = _parse_location(rng, strand, record, is_primer=True)
            simplified = int(_get_attribute_value(site, "simplified", default="0")) == 1
            if simplified and location in locations:
                # Duplicate "simplified" binding site, ignore
                continue

            locations.append(location)
            feature = SeqFeature(
                location,
                type="primer_bind",
                qualifiers=quals,
            )
            record.features.append(feature)


_packet_handlers = {
    0x00: _parse_dna_packet,
    0x05: _parse_primers_packet,
    0x06: _parse_notes_packet,
    0x0A: _parse_features_packet,
}

# Helper functions to process the XML data in
# some of the segments


def _decode(text):
    # Get rid of HTML tags in some values
    return sub("<[^>]+>", "", text)


def _get_attribute_value(node, name, default=None, error=None):
    if node.hasAttribute(name):
        return _decode(node.attributes[name].value)
    elif error:
        raise ValueError(error)
    else:
        return default


def _get_child_value(node, name, default=None, error=None):
    children = node.getElementsByTagName(name)
    if (
        children
        and children[0].childNodes
        and children[0].firstChild.nodeType == node.TEXT_NODE
    ):
        return _decode(children[0].firstChild.data)
    elif error:
        raise ValueError(error)
    else:
        return default


class SnapGeneIterator(SequenceIterator):
    """Parser for SnapGene files."""

    def __init__(self, source):
        """Parse a SnapGene file and return a SeqRecord object.

        Argument source is a file-like object or a path to a file.

        Note that a SnapGene file can only contain one sequence, so this
        iterator will always return a single record.
        """
        super().__init__(source, mode="b", fmt="SnapGene")

    def parse(self, handle):
        """Start parsing the file, and return a SeqRecord generator."""
        records = self.iterate(handle)
        return records

    def iterate(self, handle):
        """Iterate over the records in the SnapGene file."""
        record = SeqRecord(None)
        packets = _iterate(handle)
        try:
            packet_type, length, data = next(packets)
        except StopIteration:
            raise ValueError("Empty file.") from None

        if packet_type != 0x09:
            raise ValueError("The file does not start with a SnapGene cookie packet")
        _parse_cookie_packet(length, data, record)

        for (packet_type, length, data) in packets:
            handler = _packet_handlers.get(packet_type)
            if handler is not None:
                handler(length, data, record)

        if not record.seq:
            raise ValueError("No DNA packet in file")

        yield record
