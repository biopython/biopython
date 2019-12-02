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

from Bio import Alphabet
from Bio.File import as_handle
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord


class _PacketIterator:
    """Iterate over the packets of a SnapGene file.

    A SnapGene file is made of packets, each packet being a TLV-like
    structure comprising:

      - 1 single byte indicating the packet's type;
      - 1 big-endian long integer (4 bytes) indicating the length of the
        packet's data;
      - the actual data.
    """

    def __init__(self, handle):
        self.handle = handle

    def __iter__(self):
        return self

    def __next__(self):
        type = self.handle.read(1)
        if len(type) < 1:  # No more packet
            raise StopIteration
        type = unpack(">B", type)[0]

        length = self.handle.read(4)
        if len(length) < 4:
            raise ValueError("Unexpected end of packet")
        length = unpack(">I", length)[0]

        data = self.handle.read(length)
        if len(data) < length:
            raise ValueError("Unexpected end of packet")

        return (type, length, data)

    # Python2 compatibility
    def next(self):
        return self.__next__()


def _parse_dna_packet(length, data, record):
    """Parse a DNA sequence packet.

    A DNA sequence packet contains a single byte flag followed by the
    sequence itself.
    """
    if record.seq:
        raise ValueError("The file contains more than one DNA packet")

    flags, sequence = unpack(">B%ds" % (length - 1), data)
    record.seq = Seq(sequence.decode("ASCII"), alphabet=Alphabet.generic_dna)
    if flags & 0x01:
        record.annotations["topology"] = "circular"
    else:
        record.annotations["topology"] = "linear"


def _parse_notes_packet(length, data, record):
    """Parse a 'Notes' packet.

    This type of packet contains some metadata about the sequence. They
    are stored as a XML string with a 'Notes' root node.
    """
    xml = parseString(data.decode("ASCII"))
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


def _parse_features_packet(length, data, record):
    """Parse a sequence features packet.

    This packet stores sequence features (except primer binding sites,
    which are in a dedicated Primers packet). The data is a XML string
    starting with a 'Features' root node.
    """
    xml = parseString(data.decode("ASCII"))
    for feature in xml.getElementsByTagName("Feature"):
        quals = {}

        type = _get_attribute_value(feature, "type", default="misc_feature")
        label = _get_attribute_value(feature, "name")
        if label:
            quals["label"] = [label]

        strand = +1
        directionality = int(
            _get_attribute_value(feature, "directionality", default="1")
        )
        if directionality == 2:
            strand = -1

        location = None
        for segment in feature.getElementsByTagName("Segment"):
            rng = _get_attribute_value(segment, "range")
            start, end = [int(x) for x in rng.split("-")]
            # Account for SnapGene's 1-based coordinates
            start = start - 1
            if not location:
                location = FeatureLocation(start, end, strand=strand)
            else:
                location = location + FeatureLocation(start, end, strand=strand)
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

        feature = SeqFeature(location, type=type, qualifiers=quals)
        record.features.append(feature)


def _parse_primers_packet(length, data, record):
    """Parse a Primers packet.

    A Primers packet is similar to a Features packet but specifically
    stores primer binding features. The data is a XML string starting
    with a 'Primers' root node.
    """
    xml = parseString(data.decode("ASCII"))
    for primer in xml.getElementsByTagName("Primer"):
        quals = {}

        name = _get_attribute_value(primer, "name")
        if name:
            quals["label"] = [name]

        for site in primer.getElementsByTagName("BindingSite"):
            rng = _get_attribute_value(
                site, "location", error="Missing binding site location"
            )
            start, end = [int(x) for x in rng.split("-")]

            strand = int(_get_attribute_value(site, "boundStrand", default="0"))
            if strand == 1:
                strand = -1
            else:
                strand = +1

            feature = SeqFeature(
                FeatureLocation(start, end, strand=strand),
                type="primer_bind",
                qualifiers=quals,
            )
            record.features.append(feature)


_packet_handlers = {
    0x00: _parse_dna_packet,
    0x05: _parse_primers_packet,
    0x06: _parse_notes_packet,
    0x09: _parse_cookie_packet,
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


def SnapGeneIterator(handle):
    """Parse a SnapGene file and return a SeqRecord object.

    Note that a SnapGene file can only contain one sequence, so this
    iterator will always return a single record.
    """
    record = SeqRecord(None)
    n = 0

    # check if file is empty
    empty = True

    with as_handle(handle, "rb") as handle:

        for n, (type, length, data) in enumerate(_PacketIterator(handle)):
            empty = False
            if n == 0 and type != 0x09:
                raise ValueError(
                    "The file does not start with a SnapGene cookie packet"
                )

            if type in _packet_handlers:
                _packet_handlers[type](length, data, record)

        if empty:
            raise ValueError("Empty file.")

        if not record.seq:
            raise ValueError("No DNA packet in file")

        yield record
