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

import numpy


from Bio.Align import Alignment
from Bio.Align import bigbed
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, Location


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
        qStarts = numpy.array(qStarts)
        tStarts = numpy.array(tStarts)
        tBlockSizes = numpy.array(blockSizes)
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
        coordinates = numpy.array(coordinates).transpose()
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
