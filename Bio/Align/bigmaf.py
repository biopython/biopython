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

import struct
import zlib

import numpy as np

from Bio.Align import Alignment, Alignments
from Bio.Align import interfaces, bigbed, maf
from Bio.Align import _parser
from Bio.Align.bigbed import AutoSQLTable, Field
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
    ):
        """Create an AlignmentWriter object.

        Arguments:
         - target      - output stream or file name.
         - targets     - A list of SeqRecord objects with the chromosomes in the
                         order as they appear in the alignments. The sequence
                         contents in each SeqRecord may be undefined, but the
                         sequence length must be defined, as in this example:

                         SeqRecord(Seq(None, length=248956422), id="chr1")

                         If targets is None (the default value), the alignments
                         must have an attribute .targets providing the list of
                         SeqRecord objects.
         - compress    - If True (default), compress data using zlib.
                         If False, do not compress data.
        """
        super().__init__(
            target,
            bedN=3,
            declaration=declaration,
            targets=targets,
            compress=compress,
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

    def _create_alignment(self, chunk):
        chromId, chromStart, chromEnd, data = chunk
        records = []
        strands = []
        aligned_sequences = []
        annotations = {}
        j = -1
        # Avoid copying data around, as sequences in bigMaf files can be large
        starts = []
        sizes = []
        while True:
            i = j + 1
            try:
                j = data.index(b";", i)
            except ValueError:
                break
            prefix = data[i : i + 2]
            if prefix.startswith(b"#"):
                continue
            elif prefix == b"a ":
                words = data[i + 2 : j].split()
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
            elif prefix == b"s ":
                words = data[i + 2 : j].split(None, 5)
                if len(words) != 6:
                    raise ValueError(
                        "Error parsing alignment - 's' line must have 7 fields"
                    )
                src = words[0]
                start = int(words[1])
                size = int(words[2])
                strand = words[3]
                srcSize = int(words[4])
                text = words[5]
                aligned_sequences.append(text)
                seq = Seq(None, length=srcSize)
                sizes.append(size)
                if strand == b"+":
                    strands.append(+1)
                    starts.append(start)
                else:
                    strands.append(-1)
                    starts.append(srcSize - start)
                record = SeqRecord(seq, id=src.decode(), name="", description="")
                records.append(record)
            elif prefix == b"i ":
                words = data[i + 2 : j].split(None, 4)
                assert len(words) == 5
                assert words[0] == src  # from the previous "s" line
                leftStatus = words[1].decode()
                leftCount = int(words[2])
                rightStatus = words[3].decode()
                rightCount = int(words[4])
                assert leftStatus in AlignmentIterator.status_characters
                assert rightStatus in AlignmentIterator.status_characters
                record.annotations["leftStatus"] = leftStatus
                record.annotations["leftCount"] = leftCount
                record.annotations["rightStatus"] = rightStatus
                record.annotations["rightCount"] = rightCount
            elif prefix == b"e ":
                words = data[i + 2 : j].split(None, 5)
                assert len(words) == 6
                src = words[0]
                start = int(words[1])
                size = int(words[2])
                strand = words[3]
                srcSize = int(words[4])
                status = words[5].decode()
                assert status in AlignmentIterator.empty_status_characters
                sequence = Seq(None, length=srcSize)
                record = SeqRecord(sequence, id=src.decode(), name="", description="")
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
            elif prefix == b"q ":
                words = data[i + 2 : j].split(None, 1)
                assert len(words) == 2
                assert words[0] == src  # from the previous "s" line
                value = words[1].replace(b"-", b"")
                record.annotations["quality"] = value.decode()
            else:
                line = data[i:j]
                if line.strip():
                    raise ValueError(
                        f"Error parsing alignment - unexpected line:\n{line}"
                    )

        n = len(aligned_sequences)
        m = _parser.calculate_alignment_coordinates_columns(aligned_sequences)
        sizes = np.array(sizes)
        strands = np.array(strands)
        coordinates = np.empty((n, m), int)
        starts = np.array(starts)
        coordinates[:, 0] = starts
        sequences = _parser.fill_alignment_coordinates(
            aligned_sequences, coordinates, strands, sizes
        )

        for sequence, record, strand, start, size in zip(sequences, records, strands, starts, sizes):
            if len(sequence) != size:
                raise ValueError(
                    "sequence size is incorrect (found %d, expected %d)"
                    % (len(sequence), size)
                )
            srcSize = len(record.seq)
            if strand == +1:
                seq = Seq({start: sequence}, srcSize)
            else:
                seq = Seq({srcSize - start: sequence}, srcSize)
                seq = seq.reverse_complement()
            record.seq = seq
        alignment = Alignment(records, coordinates)
        if annotations is not None:
            alignment.annotations = annotations
        if score is not None:
            alignment.score = score
        return alignment
