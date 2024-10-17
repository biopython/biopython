# Copyright 2022 by Michiel de Hoon.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.Align support for the "bigmaf" multiple alignment format.

The bigMaf format stores multiple alignments in a format compatible with
the MAF (Multiple Alignment Format) format. BigMaf files are binary and are
indexed as a bigBed file.

See https://genome.ucsc.edu/goldenPath/help/bigMaf.html
"""

import re
import struct
import zlib

import numpy as np

from Bio.Align import _aligncore  # type: ignore
from Bio.Align import Alignment
from Bio.Align import Alignments
from Bio.Align import bigbed
from Bio.Align import maf
from Bio.Align.bigbed import AutoSQLTable
from Bio.Align.bigbed import Field
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

declaration = AutoSQLTable(
    "bedMaf",
    "Bed3 with MAF block",
    [
        Field(
            as_type="string",
            name="chrom",
            comment="Reference sequence chromosome or scaffold",
        ),
        Field(
            as_type="uint",
            name="chromStart",
            comment="Start position in chromosome",
        ),
        Field(
            as_type="uint",
            name="chromEnd",
            comment="End position in chromosome",
        ),
        Field(
            as_type="lstring",
            name="mafBlock",
            comment="MAF block",
        ),
    ],
)


class AlignmentWriter(bigbed.AlignmentWriter):
    """Alignment file writer for the bigMaf file format."""

    fmt = "bigMaf"

    def __init__(
        self,
        target,
        targets=None,
        compress=True,
        blockSize=256,
        itemsPerSlot=512,
    ):
        """Create an AlignmentWriter object.

        Arguments:
         - target       - output stream or file name.
         - targets      - A list of SeqRecord objects with the chromosomes in
                          the order as they appear in the alignments. The
                          sequence contents in each SeqRecord may be undefined,
                          but the sequence length must be defined, as in this
                          example:

                          SeqRecord(Seq(None, length=248956422), id="chr1")

                          If targets is None (the default value), the alignments
                          must have an attribute .targets providing the list of
                          SeqRecord objects.
         - compress     - If True (default), compress data using zlib.
                          If False, do not compress data.
                          Use compress=False for faster searching.
         - blockSize    - Number of items to bundle in r-tree.
                          See UCSC's bedToBigBed program for more information.
                          Default value is 256.
         - itemsPerSlot - Number of data points bundled at lowest level.
                          See UCSC's bedToBigBed program for more information.
                          Use itemsPerSlot=1 for faster searching.
                          Default value is 512.

        """
        super().__init__(
            target,
            bedN=3,
            declaration=declaration,
            targets=targets,
            compress=compress,
            blockSize=blockSize,
            itemsPerSlot=itemsPerSlot,
        )

    def write_file(self, stream, alignments):
        """Write the file."""
        fixed_alignments = Alignments()
        for alignment in alignments:
            if not isinstance(alignment, Alignment):
                raise TypeError("Expected an Alignment object")
            mafBlock = format(alignment, "maf")[:-1].replace("\n", ";")
            coordinates = alignment.coordinates
            if not coordinates.size:  # alignment consists of gaps only
                continue
            alignment = alignment[:2]
            reference, chromosome = alignment.target.id.split(".", 1)
            alignment.target.id = chromosome
            assert coordinates[0, 0] < coordinates[0, -1]
            alignment.annotations = {}
            alignment.annotations["mafBlock"] = mafBlock
            fixed_alignments.append(alignment)
        fixed_alignments.sort(
            key=lambda alignment: (alignment.target.id, alignment.coordinates[0, 0])
        )
        record = alignments.targets[0]
        reference, chromosome = record.id.split(".", 1)
        targets = list(alignments.targets)
        targets[0] = SeqRecord(record.seq, id=chromosome)
        fixed_alignments.targets = targets
        bigbed.AlignmentWriter(
            stream, bedN=3, declaration=declaration, compress=self.compress
        ).write(fixed_alignments)


class AlignmentIterator(bigbed.AlignmentIterator, maf.AlignmentIterator):
    """Alignment iterator for bigMaf files.

    The file may contain multiple alignments, which are loaded and returned
    incrementally.

    Alignment annotations are stored in the ``.annotations`` attribute of the
    ``Alignment`` object, except for the alignment score, which is stored as an
    attribute.  Sequence information of empty parts in the alignment block
    (sequences that connect the previous alignment block to the next alignment
    block, but do not align to the current alignment block) is stored in the
    alignment annotations under the ``"empty"`` key.  Annotations specific to
    each line in the alignment are stored in the ``.annotations`` attribute of
    the corresponding sequence record.
    """

    fmt = "bigMaf"
    mode = "b"

    def __init__(self, source):
        """Create an AlignmentIterator object.

        Arguments:
        - source - input file stream, or path to input file
        """
        self.reference = None
        super().__init__(source)

    def _read_reference(self, stream):
        # Supplemental Table 12: Binary BED-data format
        # chromId     4 bytes, unsigned
        # chromStart  4 bytes, unsigned
        # chromEnd    4 bytes, unsigned
        # rest        zero-terminated string in tab-separated format
        formatter = struct.Struct(self.byteorder + "III")
        size = formatter.size
        node = self.tree
        while True:
            try:
                children = node.children
            except AttributeError:
                break
            else:
                node = children[0]
        filepos = stream.tell()
        stream.seek(node.dataOffset)
        dataSize = 256
        data = b""
        compressed_data = b""
        while True:
            chunk = stream.read(dataSize)
            if self._compressed:
                compressed_data += chunk
                decompressor = zlib.decompressobj()
                data = decompressor.decompress(compressed_data)
            else:
                data += chunk
            try:
                i = data.index(b";s", size)
            except ValueError:
                continue
            words = data[i + 1 :].split()
            if len(words) > 2:
                break
        name = words[1]
        stream.seek(filepos)
        reference, chromosome = name.split(b".", 1)
        return reference.decode()

    def _read_header(self, stream):
        super()._read_header(stream)
        if self.reference is None:
            self.reference = self._read_reference(stream)
            self._index = 0
        self.targets[0].id = "%s.%s" % (self.reference, self.targets[0].id)

    def _create_alignment(
        self, chromId, chromStart, chromEnd, data, dataStart, dataEnd
    ):
        assert data[dataEnd - 1] == 0
        buffer = memoryview(data)
        buffer = buffer[dataStart:dataEnd]
        records = []
        strands = []
        annotations = {}
        score = None
        printed_alignment_parser = _aligncore.PrintedAlignmentParser(b";")
        j = -1
        sequences = []
        while True:
            i = j + 1
            prefix = buffer[i : i + 1]
            if prefix == b"#":
                m = re.match(b"^[^;]*", buffer[i:])
                j = i + m.span()[1]
            elif prefix == b"a":
                m = re.match(b"^[^;]*", buffer[i:])
                j = i + m.span()[1]
                line = buffer[i:j].tobytes()
                words = line[1:].split()
                for word in words:
                    key, value = word.split(b"=")
                    if key == b"score":
                        score = float(value)
                    elif key == b"pass":
                        value = int(value)
                        if value <= 0:
                            raise ValueError(
                                "pass value must be positive (found %d)" % value
                            )
                        annotations["pass"] = value
                    else:
                        raise ValueError(
                            "Unknown annotation variable '%s'" % key.decode()
                        )
            elif prefix == b"s":
                m = re.match(rb"^s\s*\S*\s*\d*\s*\d*\s*[+-]\s*\d*\s*", buffer[i:])
                j = i + m.span()[1]
                line = buffer[i:j].tobytes()
                words = line.split(None, 5)
                if len(words) != 6:
                    raise ValueError(
                        "Error parsing alignment - 's' line must have 7 fields"
                    )
                src = words[1].decode()
                start = int(words[2])
                size = int(words[3])
                strand = words[4]
                srcSize = int(words[5])
                i = j
                n, sequence = printed_alignment_parser.feed(data, dataStart + i)
                if len(sequence) != size:
                    raise ValueError(
                        "sequence size is incorrect (found %d, expected %d)"
                        % (len(sequence), size)
                    )
                seq = Seq({start: sequence}, length=srcSize)
                record = SeqRecord(seq, id=src, name="", description="")
                records.append(record)
                j = i + n
                i = j + 1
                strands.append(strand)
            elif prefix == b"i":
                m = re.match(b"^[^;]*", buffer[i:])
                j = i + m.span()[1]
                line = buffer[i:j].tobytes()
                words = line.split(None, 5)
                assert len(words) == 6
                assert words[1].decode() == src  # from the previous "s" line
                leftStatus = words[2].decode()
                leftCount = int(words[3])
                rightStatus = words[4].decode()
                rightCount = int(words[5])
                assert leftStatus in AlignmentIterator.status_characters
                assert rightStatus in AlignmentIterator.status_characters
                record.annotations["leftStatus"] = leftStatus
                record.annotations["leftCount"] = leftCount
                record.annotations["rightStatus"] = rightStatus
                record.annotations["rightCount"] = rightCount
            elif prefix == b"e":
                m = re.match(b"^[^;]*", buffer[i:])
                j = i + m.span()[1]
                line = buffer[i:j].tobytes()
                words = line.split(None, 6)
                assert len(words) == 7
                src = words[1].decode()
                start = int(words[2])
                size = int(words[3])
                strand = words[4]
                srcSize = int(words[5])
                status = words[6].decode()
                assert status in AlignmentIterator.empty_status_characters
                sequence = Seq(None, length=srcSize)
                record = SeqRecord(sequence, id=src, name="", description="")
                end = start + size
                if strand == b"+":
                    segment = (start, end)
                else:
                    segment = (srcSize - start, srcSize - end)
                empty = (record, segment, status)
                annotation = annotations.get("empty")
                if annotation is None:
                    annotation = []
                    annotations["empty"] = annotation
                annotation.append(empty)
            elif prefix == b"q":
                m = re.match(b"^[^;]*", buffer[i:])
                j = i + m.span()[1]
                line = buffer[i:j].tobytes()
                words = line.split(None, 2)
                assert len(words) == 3
                assert words[1].decode() == src  # from the previous "s" line
                value = words[2].replace(b"-", b"")
                record.annotations["quality"] = value.decode()
            elif prefix == b"\00":
                # reached the end of the alignment
                break
            else:
                raise ValueError(f"Error parsing alignment - unexpected line:\n{line}")
        shape = printed_alignment_parser.shape
        coordinates = np.empty(shape, np.int64)
        printed_alignment_parser.fill(coordinates)
        for record, strand, row in zip(records, strands, coordinates):
            if strand == b"+":
                row += record.seq.defined_ranges[0][0]
            else:  # strand == b"-"
                record.seq = record.seq.reverse_complement()
                row[:] = record.seq.defined_ranges[0][1] - row
        alignment = Alignment(records, coordinates)
        if annotations is not None:
            alignment.annotations = annotations
        if score is not None:
            alignment.score = score
        return alignment
