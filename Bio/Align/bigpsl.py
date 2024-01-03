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


from Bio.Align import Alignment, Alignments
from Bio.Align import bigbed, psl
from Bio.Align.bigbed import AutoSQLTable, Field
from Bio.Seq import Seq, reverse_complement, UndefinedSequenceError
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, Location
from Bio.SeqIO.InsdcIO import _insdc_location_string


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

    def write_file(self, stream, alignments):
        """Write the file."""
        fixed_alignments = Alignments()
        cds = self.cds
        fa = self.fa
        for alignment in alignments:
            if not isinstance(alignment, Alignment):
                raise TypeError("Expected an Alignment object")
            coordinates = alignment.coordinates
            if not coordinates.size:  # alignment consists of gaps only
                continue
            target, query = alignment.sequences
            try:
                query = query.seq
            except AttributeError:
                pass
            try:
                target = target.seq
            except AttributeError:
                pass
            tSize = len(target)
            qSize = len(query)
            # fmt: off
            dnax = None  # set to True for translated DNA aligned to protein,
            # and to False for DNA/RNA aligned to DNA/RNA  # noqa: E114, E116
            # fmt: on
            if coordinates[1, 0] > coordinates[1, -1]:
                # DNA/RNA mapped to reverse strand of DNA/RNA
                strand = "-"
                query = reverse_complement(query)
                coordinates = coordinates.copy()
                coordinates[1, :] = qSize - coordinates[1, :]
            elif coordinates[0, 0] > coordinates[0, -1]:
                # protein mapped to reverse strand of DNA
                strand = "-"
                target = reverse_complement(target)
                coordinates = coordinates.copy()
                coordinates[0, :] = tSize - coordinates[0, :]
                dnax = True
            else:
                # mapped to forward strand
                strand = "+"
            wildcard = self.wildcard
            mask = self.mask
            # variable names follow those in the PSL file format specification
            matches = 0
            misMatches = 0
            repMatches = 0
            nCount = 0
            blockSizes = []
            qStarts = []
            tStarts = []
            tStart, qStart = coordinates[:, 0]
            for tEnd, qEnd in coordinates[:, 1:].transpose():
                if tStart == tEnd:
                    qStart = qEnd
                elif qStart == qEnd:
                    tStart = tEnd
                else:
                    tCount = tEnd - tStart
                    qCount = qEnd - qStart
                    tStarts.append(tStart)
                    qStarts.append(qStart)
                    blockSizes.append(qCount)
                    if tCount == qCount:
                        assert dnax is not True
                        dnax = False
                    else:
                        # translated DNA aligned to protein, typically generated by
                        # blat -t=dnax -q=prot
                        assert tCount == 3 * qCount
                        assert dnax is not False
                        dnax = True
                    tSeq = target[tStart:tEnd]
                    qSeq = query[qStart:qEnd]
                    try:
                        tSeq = bytes(tSeq)
                    except TypeError:  # string
                        tSeq = bytes(tSeq, "ASCII")
                    except UndefinedSequenceError:  # sequence contents is unknown
                        tSeq = None
                    try:
                        qSeq = bytes(qSeq)
                    except TypeError:  # string
                        qSeq = bytes(qSeq, "ASCII")
                    except UndefinedSequenceError:  # sequence contents is unknown
                        qSeq = None
                    if tSeq is None or qSeq is None:
                        # contents of at least one sequence is unknown;
                        # count all aligned letters as matches:
                        matches += qCount
                    else:
                        if mask == "lower":
                            for u1, u2, c1 in zip(tSeq.upper(), qSeq.upper(), tSeq):
                                if u1 == wildcard or u2 == wildcard:
                                    nCount += 1
                                elif u1 == u2:
                                    if u1 == c1:
                                        matches += 1
                                    else:
                                        repMatches += 1
                                else:
                                    misMatches += 1
                        elif mask == "upper":
                            for u1, u2, c1 in zip(tSeq.lower(), qSeq.lower(), tSeq):
                                if u1 == wildcard or u2 == wildcard:
                                    nCount += 1
                                elif u1 == u2:
                                    if u1 == c1:
                                        matches += 1
                                    else:
                                        repMatches += 1
                                else:
                                    misMatches += 1
                        else:
                            for u1, u2 in zip(tSeq.upper(), qSeq.upper()):
                                if u1 == wildcard or u2 == wildcard:
                                    nCount += 1
                                elif u1 == u2:
                                    matches += 1
                                else:
                                    misMatches += 1
                    tStart = tEnd
                    qStart = qEnd
            tStarts = np.array(tStarts)
            qStarts = np.array(qStarts)
            blockSizes = np.array(blockSizes)
            try:
                matches = alignment.matches
            except AttributeError:
                pass
            try:
                misMatches = alignment.misMatches
            except AttributeError:
                pass
            try:
                repMatches = alignment.repMatches
            except AttributeError:
                pass
            try:
                nCount = alignment.nCount
            except AttributeError:
                pass
            qStart = qStarts[0]  # start of alignment in query
            qEnd = qStarts[-1] + qCount  # end of alignment in query
            oStrand = "+"
            if strand == "-":
                if dnax is True:
                    oStrand = "-"
                    qStarts = qSize - (qStarts + blockSizes)
                    qStarts = qStarts[::-1]
                    alignment.coordinates = alignment.coordinates[:, ::-1]
                else:
                    qStart, qEnd = qSize - qEnd, qSize - qStart
            if fa is True:
                oSequence = str(alignment.query.seq)
            else:
                oSequence = ""
            if cds is True:
                for feature in alignment.query.features:
                    if feature.type == "CDS":
                        oCDS = _insdc_location_string(
                            feature.location, len(alignment.query)
                        )
                        break
                else:
                    oCDS = "n/a"
            else:
                oCDS = ""
            seqType = 0
            molecule_type = alignment.query.annotations.get("molecule_type")
            if molecule_type == "DNA":
                seqType = "1"
            elif molecule_type == "protein":
                seqType = "2"
            else:
                seqType = "0"
            alignment.annotations["oChromStart"] = str(qStart)
            alignment.annotations["oChromEnd"] = str(qEnd)
            alignment.annotations["oStrand"] = oStrand
            alignment.annotations["oChromSize"] = str(qSize)
            alignment.annotations["oChromStarts"] = ",".join(map(str, qStarts))
            alignment.annotations["oSequence"] = oSequence
            alignment.annotations["oCDS"] = oCDS
            alignment.annotations["chromSize"] = str(tSize)
            alignment.annotations["match"] = str(matches)
            alignment.annotations["misMatch"] = str(misMatches)
            alignment.annotations["repMatch"] = str(repMatches)
            alignment.annotations["nCount"] = str(nCount)
            alignment.annotations["seqType"] = seqType
            fixed_alignments.append(alignment)
        fixed_alignments.sort(
            key=lambda alignment: (alignment.target.id, alignment.coordinates[0, 0])
        )
        fixed_alignments.targets = alignments.targets
        bigbed.AlignmentWriter(
            stream, bedN=12, declaration=declaration, compress=self.compress
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
                    "Expected field name '%s'; found '%s'" % (name, fields[i].name)
                )

    def _create_alignment(self, chunk):
        chromId, tStart, tEnd, rest = chunk
        words = rest.decode().split("\t")
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
        coordinates = np.array(coordinates).transpose()
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
