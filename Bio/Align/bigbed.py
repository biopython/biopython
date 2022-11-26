# Copyright 2022 by Michiel de Hoon.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.Align support for alignment files in the bigBed format.

The bigBed format stores a series of pairwise alignments in a single indexed
binary file. Typically they are used for transcript to genome alignments. As
in the BED format, the alignment positions and alignment scores are stored,
but the aligned sequences are not.

See http://genome.ucsc.edu/goldenPath/help/bigBed.html for more information.

You are expected to use this module via the Bio.Align functions.
"""

# This parser was written based on the description of the bigBed file format in
# W. J. Kent, A. S. Zweig,* G. Barber, A. S. Hinrichs, and D. Karolchik:
# "BigWig and BigBed: enabling browsing of large distributed datasets:
# Bioinformatics 26(17): 2204–2207 (2010)
# in particular the tables in the supplemental materials listing the contents
# of a bigBed file byte-by-byte.


import numpy
import struct
import zlib
from collections import namedtuple


from Bio.Align import Alignment
from Bio.Align import interfaces
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


Field = namedtuple("Field", ("type", "name", "comment"))


class AutoSQLTable:
    """AutoSQL table describing the columns of an (possibly extended) BED format."""

    def __init__(self, text=None):
        """Create an AutoSQL table describing the columns of an (extended) BED format."""
        if text is None:
            self.name = None
            self.comment = None
            self.fields = []
        else:
            assert text.endswith(chr(0))  # NULL-terminated string
            word, text = text[:-1].split(None, 1)
            assert word == "table"
            name, text = text.split(None, 1)
            assert len(name.split()) == 1
            self.name = name
            assert text.startswith('"')
            i = text.find('"', 1)
            self.comment = text[1:i]
            text = text[i + 1 :].strip()
            assert text.startswith("(")
            assert text.endswith(")")
            text = text[1:-1].strip()
            fields = []
            while text:
                i = text.index('"')
                j = text.index('"', i + 1)
                field_comment = text[i + 1 : j]
                definition = text[:i].strip()
                assert definition.endswith(";")
                field_type, field_name = definition[:-1].rsplit(None, 1)
                if field_type.endswith("]"):
                    i = field_type.index("[")
                    data_type = field_type[:i]
                else:
                    data_type = field_type
                assert data_type in (
                    "int",
                    "uint",
                    "short",
                    "ushort",
                    "byte",
                    "ubyte",
                    "float",
                    "char",
                    "string",
                    "lstring",
                )
                field = Field(field_type, field_name, field_comment)
                fields.append(field)
                text = text[j + 1 :].strip()
            self.fields = fields

    def __str__(self):
        type_width = max(len(str(field.type)) for field in self.fields)
        name_width = max(len(field.name) for field in self.fields) + 1
        lines = []
        lines.append("table %s\n" % self.name)
        lines.append('"%s"\n' % self.comment)
        lines.append("(\n")
        for field in self.fields:
            name = field.name + ";"
            lines.append(
                '   %s %s    "%s"\n'
                % (field.type.ljust(type_width), name.ljust(name_width), field.comment)
            )
        lines.append(")\n")
        return "".join(lines)


class AlignmentIterator(interfaces.AlignmentIterator):
    """Alignment iterator for bigBed files.

    The pairwise alignments stored in the bigBed file are loaded and returned
    incrementally.  Additional alignment information is stored as attributes
    of each alignment.
    """

    fmt = "bigBed"
    mode = "b"

    def _read_header(self, stream):

        # Supplemental Table 5: Common header
        # magic                  4 bytes
        # version                2 bytes
        # zoomLevels             2 bytes
        # chromosomeTreeOffset   8 bytes
        # fullDataOffset         8 bytes; points to dataCount
        # fullIndexOffset        8 bytes
        # fieldCount             2 bytes
        # definedFieldCount      2 bytes
        # autoSqlOffset          8 bytes
        # totalSummaryOffset     8 bytes
        # uncompressBufSize      4 bytes
        # reserved               8 bytes
        signature = 0x8789F2EB
        magic = stream.read(4)
        for byteorder in ("little", "big"):
            if int.from_bytes(magic, byteorder=byteorder) == signature:
                break
        else:
            raise ValueError("not a bigBed file")
        self.byteorder = byteorder
        if byteorder == "little":
            byteorder_char = "<"
        elif byteorder == "big":
            byteorder_char = ">"
        else:
            raise ValueError("Unexpected byteorder '%s'" % byteorder)

        (
            version,
            zoomLevels,
            chromosomeTreeOffset,
            fullDataOffset,
            fullIndexOffset,
            fieldCount,
            definedFieldCount,
            autoSqlOffset,
            totalSummaryOffset,
            uncompressBufSize,
        ) = struct.unpack(byteorder_char + "hhqqqhhqqixxxxxxxx", stream.read(60))

        autoSqlSize = totalSummaryOffset - autoSqlOffset
        self.declaration = self._read_autosql(
            stream, autoSqlOffset, autoSqlSize, fieldCount, definedFieldCount
        )

        stream.seek(fullDataOffset)
        dataCount = int.from_bytes(stream.read(8), byteorder=byteorder)
        self._length = dataCount

        if uncompressBufSize > 0:
            self._compressed = True
        else:
            self._compressed = False

        self.targets = self._read_chromosomes(stream, chromosomeTreeOffset)

        self.tree = self._read_index(stream, fullIndexOffset)

        self._data = self._iterate_index(stream)

    def _read_autosql(self, stream, pos, size, fieldCount, definedFieldCount):
        if definedFieldCount < 3 or definedFieldCount > 12:
            raise ValueError(
                "expected between 3 and 12 columns, found %d" % definedFieldCount
            )
        self.bedN = definedFieldCount
        stream.seek(pos)
        data = stream.read(size)
        declaration = AutoSQLTable(data.decode())
        self._analyze_fields(declaration.fields, fieldCount, definedFieldCount)
        return declaration

    def _analyze_fields(self, fields, fieldCount, definedFieldCount):
        names = (
            "chrom",
            "chromStart",
            "chromEnd",
            "name",
            "score",
            "strand",
            "thickStart",
            "thickEnd",
            "reserved",
            "blockCount",
            "blockSizes",
            "chromStarts",
        )
        for i in range(self.bedN):
            name = fields[i].name
            if name != names[i]:
                raise ValueError(
                    "Expected field name '%s'; found '%s'" % (names[i], name)
                )
        if fieldCount > definedFieldCount:
            self._custom_fields = []
        for i in range(definedFieldCount, fieldCount):
            field_name = fields[i].name
            field_type = fields[i].type
            if "[" in field_type and "]" in field_type:
                make_array = True
                field_type, _ = field_type.split("[")
                field_type = field_type.strip()
            else:
                make_array = False
            if field_type in ("int", "uint", "short", "ushort"):
                converter = int
            elif field_type in ("byte", "ubyte"):
                converter = bytes
            elif field_type == "float":
                converter = float
            elif field_type in ("float", "char", "string", "lstring"):
                converter = str
            else:
                raise Exception("Unknown field type %s" % field_type)
            if make_array:
                item_converter = converter

                def converter(data, item_converter=item_converter):
                    values = data.rstrip(",").split(",")
                    return [item_converter(value) for value in values]

            self._custom_fields.append([field_name, converter])

    def _read_chromosomes(self, stream, pos):
        byteorder = self.byteorder
        if byteorder == "little":
            byteorder_char = "<"
        elif byteorder == "big":
            byteorder_char = ">"
        else:
            raise ValueError("Unexpected byteorder '%s'" % byteorder)
        # Supplemental Table 8: Chromosome B+ tree header
        # magic     4 bytes
        # blockSize 4 bytes
        # keySize   4 bytes
        # valSize   4 bytes
        # itemCount 8 bytes
        # reserved  8 bytes
        stream.seek(pos)
        signature = 0x78CA8C91
        magic = int.from_bytes(stream.read(4), byteorder=byteorder)
        assert magic == signature
        blockSize, keySize, valSize, itemCount = struct.unpack(
            byteorder_char + "iiiqxxxxxxxx", stream.read(28)
        )
        assert valSize == 8

        Node = namedtuple("Node", ["parent", "children"])
        targets = []
        node = None
        while True:
            # Supplemental Table 9: Chromosome B+ tree node
            # isLeaf    1 byte
            # reserved  1 byte
            # count     2 bytes
            isLeaf, count = struct.unpack(byteorder_char + "?xh", stream.read(4))
            if isLeaf:
                for i in range(count):
                    # Supplemental Table 10: Chromosome B+ tree leaf item format
                    # key        keySize bytes
                    # chromId    4 bytes
                    # chromSize  4 bytes
                    key = stream.read(keySize)
                    name = key.rstrip(b"\x00").decode()
                    chromId, chromSize = struct.unpack("<II", stream.read(valSize))
                    assert chromId == len(targets)
                    sequence = Seq(None, length=chromSize)
                    record = SeqRecord(sequence, id=name)
                    targets.append(record)
            else:
                children = []
                for i in range(count):
                    # Supplemental Table 11: Chromosome B+ tree non-leaf item
                    # key          keySize bytes
                    # childOffset  8 bytes
                    key = stream.read(keySize)
                    pos = int.from_bytes(stream.read(8), byteorder)
                    children.append(pos)
                parent = node
                node = Node(parent, children)
            while True:
                if node is None:
                    assert len(targets) == itemCount
                    return targets
                children = node.children
                try:
                    pos = children.pop(0)
                except IndexError:
                    node = node.parent
                else:
                    break
            stream.seek(pos)

    def _read_index(self, stream, pos):
        byteorder = self.byteorder
        if byteorder == "little":
            byteorder_char = "<"
        elif byteorder == "big":
            byteorder_char = ">"
        else:
            raise ValueError("Unexpected byteorder '%s'" % byteorder)
        Node = namedtuple(
            "Node",
            [
                "parent",
                "children",
                "startChromIx",
                "startBase",
                "endChromIx",
                "endBase",
            ],
        )
        Leaf = namedtuple(
            "Leaf",
            [
                "parent",
                "startChromIx",
                "startBase",
                "endChromIx",
                "endBase",
                "dataOffset",
                "dataSize",
            ],
        )
        # Supplemental Table 14: R tree index header
        # magic          4 bytes
        # blockSize      4 bytes
        # itemCount      8 bytes
        # startChromIx   4 bytes
        # startBase      4 bytes
        # endChromIx     4 bytes
        # endBase        4 bytes
        # endFileOffset  8 bytes
        # itemsPerSlot   4 bytes
        # reserved       4 bytes
        stream.seek(pos)
        signature = 0x2468ACE0
        magic = int.from_bytes(stream.read(4), byteorder=byteorder)
        assert magic == signature
        (
            blockSize,
            itemCount,
            startChromIx,
            startBase,
            endChromIx,
            endBase,
            endFileOffset,
            itemsPerSlot,
        ) = struct.unpack(byteorder_char + "iqiiiiqixxxx", stream.read(44))
        root = Node(None, [], startChromIx, startBase, endChromIx, endBase)
        node = root
        itemsCounted = 0
        dataOffsets = {}
        while True:
            # Supplemental Table 15: R tree node format
            # isLeaf    1 byte
            # reserved  1 byte
            # count     2 bytes
            isLeaf, count = struct.unpack(byteorder_char + "?xh", stream.read(4))
            if isLeaf:
                children = node.children
                for i in range(count):
                    # Supplemental Table 16: R tree leaf format
                    # startChromIx   4 bytes
                    # startBase      4 bytes
                    # endChromIx     4 bytes
                    # endBase        4 bytes
                    # dataOffset     8 bytes
                    # dataSize       8 bytes
                    (
                        startChromIx,
                        startBase,
                        endChromIx,
                        endBase,
                        dataOffset,
                        dataSize,
                    ) = struct.unpack(byteorder_char + "iiiiqq", stream.read(32))
                    child = Leaf(
                        node,
                        startChromIx,
                        startBase,
                        endChromIx,
                        endBase,
                        dataOffset,
                        dataSize,
                    )
                    children.append(child)
                itemsCounted += count
                while True:
                    parent = node.parent
                    if parent is None:
                        assert itemsCounted == itemCount
                        assert not dataOffsets
                        return node
                    for index, child in enumerate(parent.children):
                        if id(node) == id(child):
                            break
                    else:
                        raise RuntimeError("Failed to find child node")
                    try:
                        node = parent.children[index + 1]
                    except IndexError:
                        node = parent
                    else:
                        break
            else:
                children = node.children
                for i in range(count):
                    # Supplemental Table 17: R tree non-leaf format
                    # startChromIx   4 bytes
                    # startBase      4 bytes
                    # endChromIx     4 bytes
                    # endBase        4 bytes
                    # dataOffset     8 bytes
                    (
                        startChromIx,
                        startBase,
                        endChromIx,
                        endBase,
                        dataOffset,
                    ) = struct.unpack(byteorder_char + "iiiiq", stream.read(24))
                    child = Node(node, [], startChromIx, startBase, endChromIx, endBase)
                    dataOffsets[id(child)] = dataOffset
                    children.append(child)
                parent = node
                node = children[0]
            pos = dataOffsets.pop(id(node))
            stream.seek(pos)

    def _iterate_index(self, stream):
        byteorder = self.byteorder
        if byteorder == "little":
            byteorder_char = "<"
        elif byteorder == "big":
            byteorder_char = ">"
        else:
            raise ValueError("Unexpected byteorder '%s'" % byteorder)
        node = self.tree
        while True:
            try:
                children = node.children
            except AttributeError:
                stream.seek(node.dataOffset)
                data = stream.read(node.dataSize)
                if self._compressed > 0:
                    data = zlib.decompress(data)
                while data:
                    # Supplemental Table 12: Binary BED-data format
                    # chromId     4 bytes
                    # chromStart  4 bytes
                    # chromEnd    4 bytes
                    # rest        zero-terminated string in tab-separated format
                    chromId, chromStart, chromEnd = struct.unpack(
                        byteorder_char + "III", data[:12]
                    )
                    rest, data = data[12:].split(b"\00", 1)
                    yield (chromId, chromStart, chromEnd, rest)
                while True:
                    parent = node.parent
                    if parent is None:
                        return
                    for index, child in enumerate(parent.children):
                        if id(node) == id(child):
                            break
                    else:
                        raise RuntimeError("Failed to find child node")
                    try:
                        node = parent.children[index + 1]
                    except IndexError:
                        node = parent
                    else:
                        break
            else:
                node = children[0]

    def _search_index(self, stream, chromIx, start, end):
        byteorder = self.byteorder
        if byteorder == "little":
            byteorder_char = "<"
        elif byteorder == "big":
            byteorder_char = ">"
        else:
            raise ValueError("Unexpected byteorder '%s'" % byteorder)
        padded_start = start - 1
        padded_end = end + 1
        node = self.tree
        while True:
            try:
                children = node.children
            except AttributeError:
                stream.seek(node.dataOffset)
                data = stream.read(node.dataSize)
                if self._compressed > 0:
                    data = zlib.decompress(data)
                while data:
                    # Supplemental Table 12: Binary BED-data format
                    # chromId     4 bytes
                    # chromStart  4 bytes
                    # chromEnd    4 bytes
                    # rest        zero-terminated string in tab-separated format
                    child_chromIx, child_chromStart, child_chromEnd = struct.unpack(
                        byteorder_char + "III", data[:12]
                    )
                    rest, data = data[12:].split(b"\00", 1)
                    if child_chromIx != chromIx:
                        continue
                    if end <= child_chromStart or child_chromEnd <= start:
                        if child_chromStart != child_chromEnd:
                            continue
                        if child_chromStart != end and child_chromEnd != start:
                            continue
                    yield (child_chromIx, child_chromStart, child_chromEnd, rest)
            else:
                visit_child = False
                for child in children:
                    if (child.endChromIx, child.endBase) < (chromIx, padded_start):
                        continue
                    if (chromIx, padded_end) < (child.startChromIx, child.startBase):
                        continue
                    visit_child = True
                    break
                if visit_child:
                    node = child
                    continue
            while True:
                parent = node.parent
                if parent is None:
                    return
                for index, child in enumerate(parent.children):
                    if id(node) == id(child):
                        break
                else:
                    raise RuntimeError("Failed to find child node")
                try:
                    node = parent.children[index + 1]
                except IndexError:
                    node = parent
                else:
                    break

    def _read_next_alignment(self, stream):
        chunk = next(self._data)
        return self._create_alignment(chunk)

    def _create_alignment(self, chunk):
        chromId, chromStart, chromEnd, rest = chunk
        if rest:
            words = rest.decode().split("\t")
        else:
            words = []
        target_record = self.targets[chromId]
        if self.bedN > 3:
            name = words[0]
        else:
            name = None
        if self.bedN > 5:
            strand = words[2]
        else:
            strand = "+"
        if self.bedN > 9:
            blockCount = int(words[6])
            blockSizes = [
                int(blockSize) for blockSize in words[7].rstrip(",").split(",")
            ]
            blockStarts = [
                int(blockStart) for blockStart in words[8].rstrip(",").split(",")
            ]
            if len(blockSizes) != blockCount:
                raise ValueError(
                    "Inconsistent number of block sizes (%d found, expected %d)"
                    % (len(blockSizes), blockCount)
                )
            if len(blockStarts) != blockCount:
                raise ValueError(
                    "Inconsistent number of block start positions (%d found, expected %d)"
                    % (len(blockStarts), blockCount)
                )
            blockSizes = numpy.array(blockSizes)
            blockStarts = numpy.array(blockStarts)
            tPosition = 0
            qPosition = 0
            coordinates = [[tPosition, qPosition]]
            for blockSize, blockStart in zip(blockSizes, blockStarts):
                if blockStart != tPosition:
                    coordinates.append([blockStart, qPosition])
                    tPosition = blockStart
                tPosition += blockSize
                qPosition += blockSize
                coordinates.append([tPosition, qPosition])
            coordinates = numpy.array(coordinates).transpose()
            qSize = sum(blockSizes)
        else:
            blockSize = chromEnd - chromStart
            coordinates = numpy.array([[0, blockSize], [0, blockSize]])
            qSize = blockSize
        coordinates[0, :] += chromStart
        query_sequence = Seq(None, length=qSize)
        query_record = SeqRecord(query_sequence, id=name)
        records = [target_record, query_record]
        if strand == "-":
            coordinates[1, :] = qSize - coordinates[1, :]
        if chromStart != coordinates[0, 0]:
            raise ValueError(
                "Inconsistent chromStart found (%d, expected %d)"
                % (chromStart, coordinates[0, 0])
            )
        if chromEnd != coordinates[0, -1]:
            raise ValueError(
                "Inconsistent chromEnd found (%d, expected %d)"
                % (chromEnd, coordinates[0, -1])
            )
        alignment = Alignment(records, coordinates)
        if len(words) > self.bedN - 3:
            alignment.annotations = {}
            for word, custom_field in zip(words[self.bedN - 3 :], self._custom_fields):
                name, converter = custom_field
                alignment.annotations[name] = converter(word)
        if self.bedN <= 4:
            return alignment
        score = words[1]
        try:
            score = float(score)
        except ValueError:
            pass
        else:
            if score.is_integer():
                score = int(score)
        alignment.score = score
        if self.bedN <= 6:
            return alignment
        alignment.thickStart = int(words[3])
        if self.bedN <= 7:
            return alignment
        alignment.thickEnd = int(words[4])
        if self.bedN <= 8:
            return alignment
        alignment.itemRgb = words[5]
        return alignment

    def __len__(self):
        return self._length

    def search(self, chromosome=None, start=None, end=None):
        """Iterate over alignments overlapping the specified chromosome region..

        This method searches the index to find alignments to the specified
        chromosome that fully or partially overlap the chromosome region
        between start and end.

        Arguments:
         - chromosome - chromosome name. If None (default value), include all
           alignments.
         - start      - starting position on the chromosome. If None (default
           value), use 0 as the starting position.
         - end        - end position on the chromosome. If None (default value),
           use the length of the chromosome as the end position.

        """
        stream = self._stream
        if chromosome is None:
            if start is not None or end is not None:
                raise ValueError(
                    "start and end must both be None if chromosome is None"
                )
        else:
            for chromIx, target in enumerate(self.targets):
                if target.id == chromosome:
                    break
            else:
                raise ValueError("Failed to find %s in alignments" % chromosome)
            if start is None:
                if end is None:
                    start = 0
                    end = len(target)
                else:
                    raise ValueError("end must be None if start is None")
            elif end is None:
                end = start + 1
        data = self._search_index(stream, chromIx, start, end)
        for chunk in data:
            alignment = self._create_alignment(chunk)
            yield alignment
