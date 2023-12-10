# Copyright 2022 by Michiel de Hoon.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.Align support for the "sam" pairwise alignment format.

The Sequence Alignment/Map (SAM) format, created by Heng Li and Richard Durbin
at the Wellcome Trust Sanger Institute, stores a series of alignments to the
genome in a single file. Typically they are used for next-generation sequencing
data. SAM files store the alignment positions for mapped sequences, and may
also store the aligned sequences and other information associated with the
sequence.

See http://www.htslib.org/ for more information.

You are expected to use this module via the Bio.Align functions.

Coordinates in the SAM format are defined in terms of one-based start
positions; the parser converts these to zero-based coordinates to be consistent
with Python and other alignment formats.
"""
from itertools import chain
import copy

import numpy as np

from Bio.Align import Alignment
from Bio.Align import interfaces
from Bio.Seq import Seq, reverse_complement, UndefinedSequenceError
from Bio.SeqRecord import SeqRecord


class AlignmentWriter(interfaces.AlignmentWriter):
    """Alignment file writer for the Sequence Alignment/Map (SAM) file format."""

    fmt = "SAM"

    def __init__(self, target, md=False):
        """Create an AlignmentWriter object.

        Arguments:
         - md - If True, calculate the MD tag from the alignment and include it
                in the output.
                If False (default), do not include the MD tag in the output.

        """
        super().__init__(target)
        self.md = md

    def write_header(self, stream, alignments):
        """Write the SAM header."""
        try:
            metadata = alignments.metadata
        except AttributeError:
            metadata = {}
        try:
            targets = alignments.targets
        except AttributeError:
            targets = {}
        values = metadata.get("HD")
        if values is not None:
            # if HD is present, then VN is required and must come first
            fields = ["@HD", "VN:%s" % values["VN"]]
            for key, value in values.items():
                if key == "VN":
                    continue
                fields.append("%s:%s" % (key, value))
            line = "\t".join(fields) + "\n"
            stream.write(line)
        for record in targets:
            fields = ["@SQ"]
            fields.append("SN:%s" % record.id)
            length = len(record.seq)
            fields.append("LN:%d" % length)
            for key, value in record.annotations.items():
                if key == "alternate_locus":
                    fields.append("AH:%s" % value)
                elif key == "names":
                    fields.append("AN:%s" % ",".join(value))
                elif key == "assembly":
                    fields.append("AS:%s" % value)
                elif key == "MD5":
                    fields.append("M5:%s" % value)
                elif key == "species":
                    fields.append("SP:%s" % value)
                elif key == "topology":
                    assert value in ("linear", "circular")
                    fields.append("PP:%s" % value)
                elif key == "URI":
                    fields.append("UR:%s" % value)
                else:
                    fields.append("%s:%s" % (key[:2], value))
            try:
                description = record.description
            except AttributeError:
                pass
            else:
                if description != "<unknown description>":
                    fields.append("DS:%s" % description)
            line = "\t".join(fields) + "\n"
            stream.write(line)
        for tag, rows in metadata.items():
            if tag == "HD":  # already written
                continue
            for row in rows:
                fields = ["@" + tag]
                for key, value in row.items():
                    fields.append("%s:%s" % (key, value))
                line = "\t".join(fields) + "\n"
                stream.write(line)

    def format_alignment(self, alignment, md=None):
        """Return a string with a single alignment formatted as one SAM line."""
        if not isinstance(alignment, Alignment):
            raise TypeError("Expected an Alignment object")
        coordinates = alignment.coordinates.transpose()
        target, query = alignment.sequences
        hard_clip_left = None
        hard_clip_right = None
        try:
            qName = query.id
        except AttributeError:
            qName = "query"
            qual = "*"
        else:
            try:
                hard_clip_left = query.annotations["hard_clip_left"]
            except (AttributeError, KeyError):
                pass
            try:
                hard_clip_right = query.annotations["hard_clip_right"]
            except (AttributeError, KeyError):
                pass
            try:
                phred = query.letter_annotations["phred_quality"]
            except (AttributeError, KeyError):
                qual = "*"
            else:
                qual = "".join(chr(value + 33) for value in phred)
            query = query.seq
        qSize = len(query)
        try:
            rname = target.id
        except AttributeError:
            rname = "target"
        else:
            target = target.seq
        if coordinates[0, 0] > coordinates[-1, 0]:
            coordinates = coordinates[::-1, :]
        if coordinates[0, 1] < coordinates[-1, 1]:  # mapped to forward strand
            flag = 0
        else:  # mapped to reverse strand
            flag = 16
            query = reverse_complement(query)
            coordinates = np.array(coordinates)
            coordinates[:, 1] = qSize - coordinates[:, 1]
            hard_clip_left, hard_clip_right = hard_clip_right, hard_clip_left
        try:
            flag |= alignment.flag
        except AttributeError:
            pass
        try:
            query = bytes(query)
        except TypeError:  # string
            pass
        except UndefinedSequenceError:
            query = "*"
        else:
            query = str(query, "ASCII")
        tStart, qStart = coordinates[0, :]
        pos = tStart + 1  # 1-based coordinate
        cigar = ""
        if hard_clip_left is not None:
            cigar += "%dH" % hard_clip_left
        if qStart > 0:
            cigar += "%dS" % qStart
        try:
            operations = alignment.operations
        except AttributeError:
            operations = None
            for tEnd, qEnd in coordinates[1:, :]:
                tCount = tEnd - tStart
                qCount = qEnd - qStart
                if tCount == 0:
                    cigar += "%dI" % qCount  # insertion to the reference
                    qStart = qEnd
                elif qCount == 0:
                    cigar += "%dD" % tCount  # deletion from the reference
                    tStart = tEnd
                else:
                    if tCount != qCount:
                        raise ValueError("Unequal step sizes in alignment")
                    cigar += "%dM" % tCount
                    tStart = tEnd
                    qStart = qEnd
        else:
            for operation, (tEnd, qEnd) in zip(operations, coordinates[1:, :]):
                tCount = tEnd - tStart
                qCount = qEnd - qStart
                if tCount == 0:
                    assert operation == ord("I")
                    cigar += "%dI" % qCount  # insertion to the reference
                    qStart = qEnd
                elif qCount == 0:
                    if operation == ord("N"):
                        cigar += "%dN" % tCount  # skipped region from the reference
                    elif operation == ord("D"):
                        cigar += "%dD" % tCount  # deletion from the reference
                    else:
                        raise ValueError(f"Unexpected operation {operation}")
                    tStart = tEnd
                else:
                    if tCount != qCount:
                        raise ValueError("Unequal step sizes in alignment")
                    assert operation == ord("M")
                    cigar += "%dM" % tCount
                    tStart = tEnd
                    qStart = qEnd
        if qEnd < qSize:
            cigar += "%dS" % (qSize - qEnd)
        if hard_clip_right is not None:
            cigar += "%dH" % hard_clip_right
        try:
            mapq = alignment.mapq
        except AttributeError:
            mapq = 255  # not available
        try:
            rnext = alignment.rnext
        except AttributeError:
            rnext = "*"
        else:
            if rnext == rname:
                rnext = "="
        try:
            pnext = alignment.pnext
        except AttributeError:
            pnext = 0
        else:
            pnext += 1  # 1-based coordinates
        tLen = 0
        fields = [
            qName,
            str(flag),
            rname,
            str(pos),
            str(mapq),
            cigar,
            rnext,
            str(pnext),
            str(tLen),
            query,
            qual,
        ]
        if md is None:
            md = self.md
        if md is True:
            if query == "*":
                raise ValueError("requested MD tag with undefined sequence")
            # calculate the MD tag from the alignment coordinates and sequences
            tStart, qStart = coordinates[0, :]
            number = 0
            md = ""
            if operations is None:
                for tEnd, qEnd in coordinates[1:, :]:
                    tCount = tEnd - tStart
                    qCount = qEnd - qStart
                    if tCount == 0:
                        # insertion to the reference
                        qStart = qEnd
                    elif qCount == 0:
                        if True:
                            # deletion from the reference
                            if number:
                                md += str(number)
                                number = 0
                            md += "^" + target[tStart:tEnd]
                        tStart = tEnd
                    else:
                        # alignment match
                        if tCount != qCount:
                            raise ValueError("Unequal step sizes in alignment")
                        for tc, qc in zip(target[tStart:tEnd], query[qStart:qEnd]):
                            if tc == qc:
                                number += 1
                            else:
                                md += str(number) + tc
                                number = 0
                        tStart = tEnd
                        qStart = qEnd
                if number:
                    md += str(number)
            else:
                for operation, (tEnd, qEnd) in zip(operations, coordinates[1:, :]):
                    tCount = tEnd - tStart
                    qCount = qEnd - qStart
                    if tCount == 0:
                        # insertion to the reference
                        qStart = qEnd
                    elif qCount == 0:
                        if operation != ord("N"):
                            # deletion from the reference
                            if number:
                                md += str(number)
                                number = 0
                            md += "^" + target[tStart:tEnd]
                        tStart = tEnd
                    else:
                        # alignment match
                        if tCount != qCount:
                            raise ValueError("Unequal step sizes in alignment")
                        for tc, qc in zip(target[tStart:tEnd], query[qStart:qEnd]):
                            if tc == qc:
                                number += 1
                            else:
                                md += str(number) + tc
                                number = 0
                        tStart = tEnd
                        qStart = qEnd
                if number:
                    md += str(number)
            field = "MD:Z:%s" % md
            fields.append(field)
        try:
            score = alignment.score
        except AttributeError:
            pass
        else:
            field = "AS:i:%.0f" % score
            fields.append(field)
        try:
            annotations = alignment.annotations
        except AttributeError:
            pass
        else:
            for key, value in annotations.items():
                if isinstance(value, int):
                    datatype = "i"
                    value = str(value)
                elif isinstance(value, float):
                    datatype = "f"
                    value = str(value)
                elif isinstance(value, str):
                    if len(value) == 1:
                        datatype = "A"
                    else:
                        datatype = "Z"
                elif isinstance(value, bytes):
                    datatype = "H"
                    value = "".join(map(str, value))
                elif isinstance(value, np.array):
                    datatype = "B"
                    if np.issubdtype(value.dtype, np.integer):
                        pass
                    elif np.issubdtype(value.dtype, float):
                        pass
                    else:
                        raise ValueError(
                            f"Array of incompatible data type {value.dtype} in annotation '{key}'"
                        )
                    value = "".join(map(str, value))
                field = f"{key}:{datatype}:{value}"
                fields.append(field)
        line = "\t".join(fields) + "\n"
        return line


class AlignmentIterator(interfaces.AlignmentIterator):
    """Alignment iterator for Sequence Alignment/Map (SAM) files.

    Each line in the file contains one genomic alignment, which are loaded
    and returned incrementally.  The following columns are stored as attributes
    of the alignment:

      - flag: The FLAG combination of bitwise flags;
      - mapq: Mapping Quality (only stored if available)
      - rnext: Reference sequence name of the primary alignment of the next read
               in the alignment (only stored if available)
      - pnext: Zero-based position of the primary alignment of the next read in
               the template (only stored if available)
      - tlen: signed observed template length (only stored if available)

    Other information associated with the alignment by its tags are stored in
    the annotations attribute of each alignment.

    Any hard clipping (clipped sequences not present in the query sequence)
    are stored as 'hard_clip_left' and 'hard_clip_right' in the annotations
    dictionary attribute of the query sequence record.

    The sequence quality, if available, is stored as 'phred_quality' in the
    letter_annotations dictionary attribute of the query sequence record.
    """

    fmt = "SAM"

    def _read_header(self, stream):
        self.metadata = {}
        self.targets = []
        for line in stream:
            if not line.startswith("@"):
                self._line = line
                break
            fields = line[1:].strip().split("\t")
            tag = fields[0]
            values = {}
            if tag == "SQ":
                annotations = {}
                description = None
                for field in fields[1:]:
                    key, value = field.split(":", 1)
                    assert len(key) == 2
                    if key == "SN":
                        rname = value
                    elif key == "LN":
                        length = int(value)
                    elif key == "AH":
                        annotations["alternate_locus"] = value
                    elif key == "AN":
                        annotations["names"] = value.split(",")
                    elif key == "AS":
                        annotations["assembly"] = value
                    elif key == "DS":
                        description = value
                    elif key == "M5":
                        annotations["MD5"] = value
                    elif key == "SP":
                        annotations["species"] = value
                    elif key == "TP":
                        assert value in ("linear", "circular")
                        annotations["topology"] = value
                    elif key == "UR":
                        annotations["URI"] = value
                    else:
                        annotations[key] = value
                sequence = Seq(None, length=length)
                record = SeqRecord(
                    sequence, id=rname, description="", annotations=annotations
                )
                if description is not None:
                    record.description = description
                self.targets.append(record)
            else:
                for field in fields[1:]:
                    key, value = field.split(":", 1)
                    assert len(key) == 2
                    values[key] = value
                if tag == "HD":
                    self.metadata[tag] = values
                else:
                    if tag not in self.metadata:
                        self.metadata[tag] = []
                    self.metadata[tag].append(values)
        self._target_indices = {
            record.id: index for index, record in enumerate(self.targets)
        }

    def _read_next_alignment(self, stream):
        try:
            line = self._line
        except AttributeError:
            lines = stream
        else:
            lines = chain([line], stream)
            del self._line
        for line in lines:
            fields = line.split()
            if len(fields) < 11:
                raise ValueError(
                    "line has %d columns; expected at least 11" % len(fields)
                )
            qname = fields[0]
            flag = int(fields[1])
            rname = fields[2]
            target_pos = int(fields[3]) - 1
            mapq = int(fields[4])
            cigar = fields[5]
            rnext = fields[6]
            pnext = int(fields[7]) - 1
            tlen = int(fields[8])
            query = fields[9]
            qual = fields[10]
            md = None
            score = None
            annotations = {}
            for field in fields[11:]:
                tag, datatype, value = field.split(":", 2)
                if tag == "AS":
                    assert datatype == "i"
                    score = int(value)
                elif tag == "MD":
                    assert datatype == "Z"
                    md = value
                else:
                    if datatype == "i":
                        value = int(value)
                    elif datatype == "f":
                        value = float(value)
                    elif datatype in ("A", "Z"):  # string
                        pass
                    elif datatype == "H":
                        n = len(value)
                        value = bytes(int(value[i : i + 2]) for i in range(0, n, 2))
                    elif datatype == "B":
                        letter = value[0]
                        value = value[1:].split(",")
                        if letter in "cCsSiI":
                            dtype = int
                        elif letter == "f":
                            dtype = float
                        else:
                            raise ValueError(
                                f"Unknown number type '{letter}' in tag '{field}'"
                            )
                        value = np.array(value, dtype)
                    annotations[tag] = value
            if flag & 0x10:
                strand = "-"
            else:
                strand = "+"
            hard_clip_left = None
            hard_clip_right = None
            store_operations = False
            if flag & 0x4:  # unmapped
                target = None
                coordinates = None
            elif md is None:
                query_pos = 0
                coordinates = [[target_pos, query_pos]]
                number = ""
                operations = bytearray()
                for letter in cigar:
                    if letter == "M":
                        # M: alignment match
                        length = int(number)
                        target_pos += length
                        query_pos += length
                    elif letter in "=X":
                        # =: sequence match
                        # X: sequence mismatch
                        length = int(number)
                        target_pos += length
                        query_pos += length
                        store_operations = True
                    elif letter == "I":
                        # I: insertion to the reference
                        length = int(number)
                        query_pos += length
                    elif letter == "S":
                        # S: soft clipping
                        length = int(number)
                        if query_pos == 0:
                            coordinates[0][1] += length
                        query_pos += length
                        number = ""
                        continue
                    elif letter == "D":
                        # D: deletion from the reference
                        length = int(number)
                        target_pos += length
                    elif letter == "N":
                        # N: skipped region from the reference
                        length = int(number)
                        target_pos += length
                        store_operations = True
                    elif letter == "H":  # hard clipping
                        if query_pos == 0:
                            hard_clip_left = int(number)
                        else:
                            hard_clip_right = int(number)
                        number = ""
                        continue
                    elif letter == "P":  # padding
                        raise NotImplementedError(
                            "padding operator is not yet implemented"
                        )
                    else:
                        number += letter
                        continue
                    coordinates.append([target_pos, query_pos])
                    operations.append(ord(letter))
                    number = ""
                index = self._target_indices.get(rname)
                if index is None:
                    if self.targets:
                        raise ValueError(f"Found target {rname} missing from header")
                    target = SeqRecord(None, id=rname, description="")
                else:
                    target = self.targets[index]
            else:
                query_pos = 0
                coordinates = [[target_pos, query_pos]]
                seq = query
                target = ""
                starts = [target_pos]
                size = 0
                sizes = []
                number = ""
                operations = bytearray()
                for letter in cigar:
                    if letter in "M":
                        # M: alignment match
                        length = int(number)
                        target_pos += length
                        query_pos += length
                        target += seq[:length]
                        seq = seq[length:]
                        size += length
                    elif letter in "=X":
                        # =: sequence match
                        # X: sequence mismatch
                        length = int(number)
                        target_pos += length
                        query_pos += length
                        target += seq[:length]
                        seq = seq[length:]
                        size += length
                        store_operations = True
                    elif letter == "I":
                        # I: insertion to the reference
                        length = int(number)
                        query_pos += length
                        seq = seq[length:]
                    elif letter == "S":
                        # S: soft clipping
                        length = int(number)
                        if query_pos == 0:
                            coordinates[0][1] += length
                        query_pos += length
                        seq = seq[length:]
                        number = ""
                        continue
                    elif letter == "D":  # deletion from the reference
                        length = int(number)
                        target_pos += length
                        size += length
                        starts.append(target_pos)
                        sizes.append(size)
                        size = 0
                    elif letter == "N":  # skipped region from the reference
                        length = int(number)
                        target_pos += length
                        starts.append(target_pos)
                        sizes.append(size)
                        size = 0
                        store_operations = True
                    elif letter == "H":
                        # hard clipping (clipped sequences not present in sequence)
                        if query_pos == 0:
                            hard_clip_left = int(number)
                        else:
                            hard_clip_right = int(number)
                        number = ""
                        continue
                    elif letter == "P":  # padding
                        raise NotImplementedError(
                            "padding operator is not yet implemented"
                        )
                    else:
                        number += letter
                        continue
                    coordinates.append([target_pos, query_pos])
                    operations.append(ord(letter))
                    number = ""
                sizes.append(size)
                seq = target
                target = ""
                number = ""
                letters = iter(md)
                for letter in letters:
                    if letter in "ACGTNacgtn":
                        if number:
                            number = int(number)
                            target += seq[:number]
                            seq = seq[number:]
                            number = ""
                        target += letter
                        seq = seq[1:]
                    elif letter == "^":
                        if number:
                            number = int(number)
                            target += seq[:number]
                            seq = seq[number:]
                            number = ""
                        for letter in letters:
                            if letter not in "ACGTNacgtn":
                                break
                            target += letter
                        else:
                            break
                        number = letter
                    else:
                        number += letter
                if number:
                    number = int(number)
                    target += seq[:number]
                seq = target
                index = self._target_indices[rname]
                target = copy.deepcopy(self.targets[index])
                length = len(target.seq)
                data = {}
                index = 0
                for start, size in zip(starts, sizes):
                    data[start] = seq[index : index + size]
                    index += size
                target.seq = Seq(data, length=length)
            if coordinates is not None:
                coordinates = np.array(coordinates).transpose()
                if strand == "-":
                    coordinates[1, :] = query_pos - coordinates[1, :]
            if query == "*":
                length = query_pos
                sequence = Seq(None, length=length)
            else:
                sequence = Seq(query)
                if not (flag & 0x4):  # not unmapped
                    assert len(query) == query_pos
                    if strand == "-":
                        sequence = sequence.reverse_complement()
            query = SeqRecord(sequence, id=qname, description="")
            if strand == "-":
                hard_clip_left, hard_clip_right = hard_clip_right, hard_clip_left
            if hard_clip_left is not None:
                query.annotations["hard_clip_left"] = hard_clip_left
            if hard_clip_right is not None:
                query.annotations["hard_clip_right"] = hard_clip_right
            if qual != "*":
                phred = [ord(c) - 33 for c in qual]
                query.letter_annotations["phred_quality"] = phred
            records = [target, query]
            alignment = Alignment(records, coordinates)
            alignment.flag = flag
            if mapq != 255:
                alignment.mapq = mapq
            if rnext == "=":
                alignment.rnext = rname
            elif rnext != "*":
                alignment.rnext = rnext
            if pnext >= 0:
                alignment.pnext = pnext
            if tlen != 0:
                alignment.tlen = tlen
            if score is not None:
                alignment.score = score
            if annotations:
                alignment.annotations = annotations
            if hard_clip_left is not None:
                alignment.hard_clip_left = hard_clip_left
            if hard_clip_right is not None:
                alignment.hard_clip_right = hard_clip_right
            if store_operations:
                alignment.operations = operations
            return alignment
