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
from io import StringIO


from Bio.Align import Alignment, Alignments
from Bio.Align import interfaces, bigbed, maf
from Bio.Align.bigbed import AutoSQLTable, Field
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
            mafBlock = format(alignment, "maf")
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

    def _read_header(self, stream):
        super()._read_header(stream)
        if self.reference is None:
            alignment = next(self)
            self.reference, chromosome = alignment.sequences[0].id.split(".", 1)
            self.rewind()
        else:
            self.targets[0].id = "%s.%s" % (self.reference, self.targets[0].id)

    def _create_alignment(self, chunk):
        chromId, chromStart, chromEnd, rest = chunk
        data = rest.decode().replace(";", "\n")
        stream = StringIO()
        stream.write(data)
        stream.seek(0)
        line = next(stream)
        alignment = maf.AlignmentIterator._create_alignment(self, line, stream)
        return alignment
