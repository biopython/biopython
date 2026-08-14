# Copyright 2022 by Michiel de Hoon.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.Align support for alignment files in the bigPsl format.

A bigPsl file is a bigBed file with a BED12+13 format consisting of the 12
predefined BED fields and 13 custom fields defined in the autoSql file
bigPsl.as. This module uses the Bio.Align.bigbed module to parse the file,
but stores the data in a PSL-consistent manner as defined in bigPsl.as. As the
bigPsl format is a special case of the bigBed format, bigPsl files are binary
and are indexed as bigBed files.

See http://genome.ucsc.edu/goldenPath/help/bigPsl.html for more information.

You are expected to use this module via the Bio.Align functions.
"""

import numpy as np

from Bio.Align import Alignment
from Bio.Align import Alignments
from Bio.Align import bigbed
from Bio.Align.bigbed import AutoSQLTable
from Bio.Align.bigbed import Field
from Bio.Seq import reverse_complement
from Bio.Seq import Seq
from Bio.Seq import UndefinedSequenceError
from Bio.SeqFeature import Location
from Bio.SeqFeature import SeqFeature
from Bio.SeqIO.InsdcIO import _insdc_location_string
from Bio.SeqRecord import SeqRecord

declaration = AutoSQLTable(
    "bigPsl",
    "bigPsl pairwise alignment",
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
            as_type="string",
            name="name",
            comment="Name or ID of item, ideally both human readable and unique",
        ),
        Field(
            as_type="uint",
            name="score",
            comment="Score (0-1000)",
        ),
        Field(
            as_type="char[1]",
            name="strand",
            comment="+ or - indicates whether the query aligns to the + or - strand on the reference",
        ),
        Field(
            as_type="uint",
            name="thickStart",
            comment="Start of where display should be thick (start codon)",
        ),
        Field(
            as_type="uint",
            name="thickEnd",
            comment="End of where display should be thick (stop codon)",
        ),
        Field(
            as_type="uint",
            name="reserved",
            comment="RGB value (use R,G,B string in input file)",
        ),
        Field(
            as_type="int",
            name="blockCount",
            comment="Number of blocks",
        ),
        Field(
            as_type="int[blockCount]",
            name="blockSizes",
            comment="Comma separated list of block sizes",
        ),
        Field(
            as_type="int[blockCount]",
            name="chromStarts",
            comment="Start positions relative to chromStart",
        ),
        Field(
            as_type="uint",
            name="oChromStart",
            comment="Start position in other chromosome",
        ),
        Field(
            as_type="uint",
            name="oChromEnd",
            comment="End position in other chromosome",
        ),
        Field(
            as_type="char[1]",
            name="oStrand",
            comment="+ or -, - means that psl was reversed into BED-compatible coordinates",
        ),
        Field(
            as_type="uint",
            name="oChromSize",
            comment="Size of other chromosome.",
        ),
        Field(
            as_type="int[blockCount]",
            name="oChromStarts",
            comment="Start positions relative to oChromStart or from oChromStart+oChromSize depending on strand",
        ),
        Field(
            as_type="lstring",
            name="oSequence",
            comment="Sequence on other chrom (or edit list, or empty)",
        ),
        Field(
            as_type="string",
            name="oCDS",
            comment="CDS in NCBI format",
        ),
        Field(
            as_type="uint",
            name="chromSize",
            comment="Size of target chromosome",
        ),
        Field(
            as_type="uint",
            name="match",
            comment="Number of bases matched.",
        ),
        Field(
            as_type="uint",
            name="misMatch",
            comment="Number of bases that don't match",
        ),
        Field(
            as_type="uint",
            name="repMatch",
            comment="Number of bases that match but are part of repeats",
        ),
        Field(
            as_type="uint",
            name="nCount",
            comment="Number of 'N' bases",
        ),
        Field(
            as_type="uint",
            name="seqType",
            comment="0=empty, 1=nucleotide, 2=amino_acid",
        ),
    ],
)


class AlignmentWriter(bigbed.AlignmentWriter):
    """Alignment file writer for the bigPsl file format."""

    fmt = "bigPsl"

    def __init__(
        self,
        target,
        targets=None,
        compress=True,
        extraIndex=(),
        cds=False,
        fa=False,
        mask=None,
        wildcard="N",
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
         - extraIndex  - List of strings with the names of extra columns to be
                         indexed.
                         Default value is an empty list.
         - cds         - If True, look for a query feature of type CDS and write
                         it in NCBI style in the PSL file (default: False).
         - fa          - If True, include the query sequence in the PSL file
                         (default: False).
         - mask        - Specify if repeat regions in the target sequence are
                         masked and should be reported in the `repMatches` field
                         instead of in the `matches` field.
                         Acceptable values are
                         None   : no masking (default);
                         "lower": masking by lower-case characters;
                         "upper": masking by upper-case characters.
         - wildcard    - Report alignments to the wildcard character in the
                         target or query sequence in the `nCount` field instead
                         of in the `matches`, `misMatches`, or `repMatches`
                         fields.
                         Default value is 'N'.
        """
        super().__init__(
            target,
            bedN=12,
            declaration=declaration,
            targets=targets,
            compress=compress,
            extraIndex=extraIndex,
        )
        self.cds = cds
        self.fa = fa
        self.mask = mask
        self.wildcard = wildcard

    def _get_sequences(self, alignment):
        """Return target and query sequences from an alignment."""
        target, query = alignment.sequences

        try:
            query = query.seq
        except AttributeError:
            pass

        try:
            target = target.seq
        except AttributeError:
            pass

        return target, query

    def _prepare_coordinates(self, alignment, target, query):
        """Prepare coordinates and determine the alignment strand."""
        coordinates = alignment.coordinates
        t_size = len(target)
        q_size = len(query)
        dnax = None

        if coordinates[1, 0] > coordinates[1, -1]:
            strand = "-"
            query = reverse_complement(query)
            coordinates = coordinates.copy()
            coordinates[1, :] = q_size - coordinates[1, :]
        elif coordinates[0, 0] > coordinates[0, -1]:
            strand = "-"
            target = reverse_complement(target)
            coordinates = coordinates.copy()
            coordinates[0, :] = t_size - coordinates[0, :]
            dnax = True
        else:
            strand = "+"

        return coordinates, target, query, strand, dnax

    def _count_matches(self, target_sequence, query_sequence):
        """Count matches, mismatches, repeat matches and wildcards."""
        matches = 0
        mis_matches = 0
        rep_matches = 0
        n_count = 0

        if target_sequence is None or query_sequence is None:
            return 1, 0, 0, 0

        mask = self.mask
        wildcard = self.wildcard

        if mask == "lower":
            target_upper = target_sequence.upper()
            query_upper = query_sequence.upper()

            for u1, u2, c1 in zip(
                target_upper,
                query_upper,
                target_sequence,
            ):
                if wildcard in (u1, u2):
                    n_count += 1
                elif u1 == u2:
                    if u1 == c1:
                        matches += 1
                    else:
                        rep_matches += 1
                else:
                    mis_matches += 1

        elif mask == "upper":
            target_lower = target_sequence.lower()
            query_lower = query_sequence.lower()

            for u1, u2, c1 in zip(
                target_lower,
                query_lower,
                target_sequence,
            ):
                if wildcard in (u1, u2):
                    n_count += 1
                elif u1 == u2:
                    if u1 == c1:
                        matches += 1
                    else:
                        rep_matches += 1
                else:
                    mis_matches += 1

        else:
            for u1, u2 in zip(
                target_sequence.upper(),
                query_sequence.upper(),
            ):
                if wildcard in (u1, u2):
                    n_count += 1
                elif u1 == u2:
                    matches += 1
                else:
                    mis_matches += 1

        return matches, mis_matches, rep_matches, n_count

    def _get_block_sequences(
        self,
        target,
        query,
        t_start,
        t_end,
        q_start,
        q_end,
    ):
        """Return byte sequences for an alignment block."""
        target_sequence = target[t_start:t_end]
        query_sequence = query[q_start:q_end]

        try:
            target_sequence = bytes(target_sequence)
        except TypeError:
            target_sequence = bytes(target_sequence, "ASCII")
        except UndefinedSequenceError:
            target_sequence = None

        try:
            query_sequence = bytes(query_sequence)
        except TypeError:
            query_sequence = bytes(query_sequence, "ASCII")
        except UndefinedSequenceError:
            query_sequence = None

        return target_sequence, query_sequence

    def _process_blocks(
        self,
        coordinates,
        target,
        query,
        dnax,
    ):
        """Extract alignment blocks and calculate alignment statistics."""
        matches = 0
        mis_matches = 0
        rep_matches = 0
        n_count = 0

        block_sizes = []
        q_starts = []
        t_starts = []

        t_start, q_start = coordinates[:, 0]
        q_count = 0

        for t_end, q_end in coordinates[:, 1:].transpose():
            if t_start == t_end:
                q_start = q_end
                continue

            if q_start == q_end:
                t_start = t_end
                continue

            t_count = t_end - t_start
            q_count = q_end - q_start

            t_starts.append(t_start)
            q_starts.append(q_start)
            block_sizes.append(q_count)

            if t_count == q_count:
                assert dnax is not True
                dnax = False
            else:
                assert t_count == 3 * q_count
                assert dnax is not False
                dnax = True

            target_sequence, query_sequence = self._get_block_sequences(
                target,
                query,
                t_start,
                t_end,
                q_start,
                q_end,
            )

            if target_sequence is None or query_sequence is None:
                matches += q_count
            else:
                block_matches = self._count_matches(
                    target_sequence,
                    query_sequence,
                )
                matches += block_matches[0]
                mis_matches += block_matches[1]
                rep_matches += block_matches[2]
                n_count += block_matches[3]

            t_start = t_end
            q_start = q_end

        return (
            np.array(t_starts),
            np.array(q_starts),
            np.array(block_sizes),
            matches,
            mis_matches,
            rep_matches,
            n_count,
            dnax,
            q_count,
        )

    def _get_alignment_statistics(
        self,
        alignment,
        t_starts,
        q_starts,
        block_sizes,
        matches,
        mis_matches,
        rep_matches,
        n_count,
    ):
        """Use alignment-provided statistics when available."""
        try:
            matches = alignment.matches
        except AttributeError:
            pass

        try:
            mis_matches = alignment.misMatches
        except AttributeError:
            pass

        try:
            rep_matches = alignment.repMatches
        except AttributeError:
            pass

        try:
            n_count = alignment.nCount
        except AttributeError:
            pass

        return (
            t_starts,
            q_starts,
            block_sizes,
            matches,
            mis_matches,
            rep_matches,
            n_count,
        )

    def _prepare_query_coordinates(
        self,
        q_starts,
        block_sizes,
        q_size,
        q_start,
        q_end,
        strand,
        dnax,
        alignment,
    ):
        """Prepare query coordinates for the bigPsl representation."""
        o_strand = "+"

        if strand == "-":
            if dnax is True:
                o_strand = "-"
                q_starts = q_size - (q_starts + block_sizes)
                q_starts = q_starts[::-1]
                alignment.coordinates = alignment.coordinates[:, ::-1]
            else:
                q_start, q_end = q_size - q_end, q_size - q_start

        return q_starts, q_start, q_end, o_strand

    def _get_cds(self, alignment):
        """Return the CDS location in NCBI format."""
        if self.cds is False:
            return ""

        for feature in alignment.query.features:
            if feature.type == "CDS":
                return _insdc_location_string(
                    feature.location,
                    len(alignment.query),
                )

        return "n/a"

    def _get_sequence_type(self, alignment):
        """Return the bigPsl sequence type."""
        molecule_type = alignment.query.annotations.get("molecule_type")

        if molecule_type == "DNA":
            return "1"

        if molecule_type == "protein":
            return "2"

        return "0"

    def _set_alignment_annotations(
        self,
        alignment,
        q_start,
        q_end,
        o_strand,
        q_size,
        q_starts,
        o_sequence,
        o_cds,
        t_size,
        matches,
        mis_matches,
        rep_matches,
        n_count,
        seq_type,
    ):
        """Store bigPsl-specific information in the alignment."""
        alignment.annotations["oChromStart"] = str(q_start)
        alignment.annotations["oChromEnd"] = str(q_end)
        alignment.annotations["oStrand"] = o_strand
        alignment.annotations["oChromSize"] = str(q_size)
        alignment.annotations["oChromStarts"] = ",".join(
            map(str, q_starts)
        )
        alignment.annotations["oSequence"] = o_sequence
        alignment.annotations["oCDS"] = o_cds
        alignment.annotations["chromSize"] = str(t_size)
        alignment.annotations["match"] = str(matches)
        alignment.annotations["misMatch"] = str(mis_matches)
        alignment.annotations["repMatch"] = str(rep_matches)
        alignment.annotations["nCount"] = str(n_count)
        alignment.annotations["seqType"] = seq_type

    def _prepare_query_metadata(
        self,
        alignment,
        q_start,
        q_end,
        q_size,
        q_starts,
        strand,
        dnax,
    ):
        """Prepare query-related metadata for the bigPsl record."""
        q_starts, q_start, q_end, o_strand = (
            self._prepare_query_coordinates(
                q_starts,
                self._current_block_sizes,
                q_size,
                q_start,
                q_end,
                strand,
                dnax,
                alignment,
            )
        )

        if self.fa is True:
            o_sequence = str(alignment.query.seq)
        else:
            o_sequence = ""

        o_cds = self._get_cds(alignment)
        seq_type = self._get_sequence_type(alignment)

        return (
            q_starts,
            q_start,
            q_end,
            o_strand,
            o_sequence,
            o_cds,
            seq_type,
        )

    def write_file(self, stream, alignments):
        """Write the file."""
        fixed_alignments = Alignments()

        for alignment in alignments:
            if not isinstance(alignment, Alignment):
                raise TypeError("Expected an Alignment object")

            coordinates = alignment.coordinates

            if not coordinates.size:
                continue

            target, query = self._get_sequences(alignment)
            t_size = len(target)
            q_size = len(query)

            (
                coordinates,
                target,
                query,
                strand,
                dnax,
            ) = self._prepare_coordinates(
                alignment,
                target,
                query,
            )

            (
                t_starts,
                q_starts,
                block_sizes,
                matches,
                mis_matches,
                rep_matches,
                n_count,
                dnax,
                q_count,
            ) = self._process_blocks(
                coordinates,
                target,
                query,
                dnax,
            )

            (
                t_starts,
                q_starts,
                block_sizes,
                matches,
                mis_matches,
                rep_matches,
                n_count,
            ) = self._get_alignment_statistics(
                alignment,
                t_starts,
                q_starts,
                block_sizes,
                matches,
                mis_matches,
                rep_matches,
                n_count,
            )

            q_start = q_starts[0]
            q_end = q_starts[-1] + q_count

            self._current_block_sizes = block_sizes

            (
                q_starts,
                q_start,
                q_end,
                o_strand,
                o_sequence,
                o_cds,
                seq_type,
            ) = self._prepare_query_metadata(
                alignment,
                q_start,
                q_end,
                q_size,
                q_starts,
                strand,
                dnax,
            )

            self._set_alignment_annotations(
                alignment,
                q_start,
                q_end,
                o_strand,
                q_size,
                q_starts,
                o_sequence,
                o_cds,
                t_size,
                matches,
                mis_matches,
                rep_matches,
                n_count,
                seq_type,
            )

            fixed_alignments.append(alignment)

        fixed_alignments.sort(
            key=lambda alignment: (
                alignment.target.id,
                alignment.coordinates[0, 0],
            )
        )
        fixed_alignments.targets = alignments.targets

        bigbed.AlignmentWriter(
            stream,
            bedN=12,
            declaration=declaration,
            compress=self.compress,
        ).write(fixed_alignments)


class AlignmentIterator(bigbed.AlignmentIterator):
    """Alignment iterator for bigPsl files.

    The pairwise alignments stored in the bigPsl file are loaded and returned
    incrementally.  Additional alignment information is stored as attributes
    of each alignment.
    """

    fmt = "bigPsl"

    def _analyze_fields(self, fields, fieldCount, definedFieldCount):
        names = (
            "chrom",
            "chromStart",
            "chromEnd",
            "name",  # 0
            "score",  # 1
            "strand",  # 2
            "thickStart",  # 3
            "thickEnd",  # 4
            "reserved",  # 5
            "blockCount",  # 6
            "blockSizes",  # 7
            "chromStarts",  # 8
            "oChromStart",  # 9
            "oChromEnd",  # 10
            "oStrand",  # 11
            "oChromSize",  # 12
            "oChromStarts",  # 13
            "oSequence",  # 14
            "oCDS",  # 15
            "chromSize",  # 16
            "match",  # 17
            "misMatch",  # 18
            "repMatch",  # 19
            "nCount",  # 20
            "seqType",  # 21
        )
        for i, name in enumerate(names):
            if name != fields[i].name:
                raise ValueError(
                    f"Expected field name '{name}'; found '{fields[i].name}'"
                )

    def _create_alignment(self, chromId, tStart, tEnd, rest, dataStart, dataEnd):
        assert rest[dataEnd - 1] == 0
        words = rest[dataStart : dataEnd - 1].decode().split("\t")
        if len(words) != 22:
            raise ValueError(
                "Unexpected number of fields (%d, expected 22)" % len(words)
            )
        target_record = self.targets[chromId]
        tSize = int(words[16])
        if len(target_record) != tSize:
            raise ValueError(
                "Unexpected chromosome size %d (expected %d)"
                % (tSize, len(target_record))
            )
        strand = words[2]
        qName = words[0]
        qSize = int(words[12])
        blockCount = int(words[6])
        blockSizes = [int(blockSize) for blockSize in words[7].rstrip(",").split(",")]
        tStarts = [int(start) for start in words[8].rstrip(",").split(",")]
        qStarts = [int(start) for start in words[13].rstrip(",").split(",")]
        if len(blockSizes) != blockCount:
            raise ValueError(
                "Inconsistent number of blocks (%d found, expected %d)"
                % (len(blockSizes), blockCount)
            )
        if len(qStarts) != blockCount:
            raise ValueError(
                "Inconsistent number of query start positions (%d found, expected %d)"
                % (len(qStarts), blockCount)
            )
        if len(tStarts) != blockCount:
            raise ValueError(
                "Inconsistent number of target start positions (%d found, expected %d)"
                % (len(qStarts), blockCount)
            )
        qStarts = np.array(qStarts)
        tStarts = np.array(tStarts)
        tBlockSizes = np.array(blockSizes)
        query_sequence = words[14]
        if query_sequence == "":
            query_sequence = Seq(None, length=qSize)
        else:
            query_sequence = Seq(query_sequence)
            if len(query_sequence) != qSize:
                raise ValueError(
                    "Inconsistent query sequence length (%d, expected %d)"
                    % (len(query_sequence), qSize)
                )
        query_record = SeqRecord(query_sequence, id=qName)
        cds = words[15]
        if cds and cds != "n/a":
            location = Location.fromstring(cds)
            feature = SeqFeature(location, type="CDS")
            query_record.features.append(feature)
        seqType = words[21]
        if seqType == "0":
            qBlockSizes = tBlockSizes
        elif seqType == "1":
            query_record.annotations["molecule_type"] = "DNA"
            qBlockSizes = tBlockSizes
        elif seqType == "2":
            query_record.annotations["molecule_type"] = "protein"
            qBlockSizes = tBlockSizes // 3
        else:
            raise ValueError("Unexpected sequence type '%s'" % seqType)
        tStarts += tStart
        qStrand = words[11]
        if qStrand == "-" and strand == "-":
            tStart, tEnd = tEnd, tStart
            qStarts = qSize - qStarts - qBlockSizes
            tStarts = tSize - tStarts - tBlockSizes
            qStarts = qStarts[::-1]
            tStarts = tStarts[::-1]
            qBlockSizes = qBlockSizes[::-1]
            tBlockSizes = tBlockSizes[::-1]
        qPosition = qStarts[0]
        tPosition = tStarts[0]
        coordinates = [[tPosition, qPosition]]
        for tB, qB, tS, qS in zip(tBlockSizes, qBlockSizes, tStarts, qStarts):
            if tS != tPosition:
                coordinates.append([tS, qPosition])
                tPosition = tS
            if qS != qPosition:
                coordinates.append([tPosition, qS])
                qPosition = qS
            tPosition += tB
            qPosition += qB
            coordinates.append([tPosition, qPosition])
        coordinates = np.array(coordinates, np.intp).transpose()
        qStart = int(words[9])
        qEnd = int(words[10])
        if strand == "-":
            if qStrand == "-":
                coordinates[0, :] = tSize - coordinates[0, :]
            else:
                qStart, qEnd = qEnd, qStart
                coordinates[1, :] = qSize - coordinates[1, :]
        if tStart != coordinates[0, 0]:
            raise ValueError(
                "Inconsistent tStart found (%d, expected %d)"
                % (tStart, coordinates[0, 0])
            )
        if tEnd != coordinates[0, -1]:
            raise ValueError(
                "Inconsistent tEnd found (%d, expected %d)" % (tEnd, coordinates[0, -1])
            )
        if qStart != coordinates[1, 0]:
            raise ValueError(
                "Inconsistent qStart found (%d, expected %d)"
                % (qStart, coordinates[1, 0])
            )
        if qEnd != coordinates[1, -1]:
            raise ValueError(
                "Inconsistent qEnd found (%d, expected %d)" % (qEnd, coordinates[1, -1])
            )
        records = [target_record, query_record]
        alignment = Alignment(records, coordinates)
        alignment.annotations = {}
        score = words[1]
        try:
            score = float(score)
        except ValueError:
            pass
        else:
            if score.is_integer():
                score = int(score)
        alignment.score = score
        alignment.thickStart = int(words[3])
        alignment.thickEnd = int(words[4])
        alignment.itemRgb = words[5]
        alignment.matches = int(words[17])
        alignment.misMatches = int(words[18])
        alignment.repMatches = int(words[19])
        alignment.nCount = int(words[20])
        return alignment
