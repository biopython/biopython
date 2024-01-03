# Copyright 2023 by Michiel de Hoon.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.Align support for the "chain" pairwise alignment format.

As described by UCSC, a chain file stores a series of pairwise alignments in a
single file. Typically they are used for genome to genome alignments. Chain
files store the lengths of the aligned segments, alignment gaps, and alignment
scores, but do not store the aligned sequences.

See https://genome.ucsc.edu/goldenPath/help/chain.html.

You are expected to use this module via the Bio.Align functions.

Coordinates in the chain file format are defined in terms of zero-based start
positions (like Python) and aligning region sizes.
"""
import numpy as np


from Bio.Align import Alignment
from Bio.Align import interfaces
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class AlignmentWriter(interfaces.AlignmentWriter):
    """Alignment file writer for the UCSC chain file format."""

    fmt = "chain"

    def format_alignment(self, alignment):
        """Return a string with one alignment formatted as a chain block."""
        if not isinstance(alignment, Alignment):
            raise TypeError("Expected an Alignment object")
        coordinates = alignment.coordinates
        if not coordinates.size:  # alignment consists of gaps only
            return ""
        target, query = alignment.sequences
        try:
            qName = query.id
        except AttributeError:
            qName = "query"
        try:
            query = query.seq
        except AttributeError:
            pass
        try:
            tName = target.id
        except AttributeError:
            tName = "target"
        try:
            target = target.seq
        except AttributeError:
            pass
        tSize = len(target)
        qSize = len(query)
        lines = []
        try:
            score = alignment.score
        except AttributeError:
            score = 0
        if coordinates[0, 0] > coordinates[0, -1]:
            tStrand = "-"
            coordinates = coordinates.copy()
            coordinates[0, :] = tSize - coordinates[0, :]
        else:
            tStrand = "+"
        if coordinates[1, 0] > coordinates[1, -1]:
            qStrand = "-"
            coordinates = coordinates.copy()
            coordinates[1, :] = qSize - coordinates[1, :]
        else:
            qStrand = "+"
        tStart = coordinates[0, 0]
        tEnd = coordinates[0, -1]
        qStart = coordinates[1, 0]
        qEnd = coordinates[1, -1]
        try:
            chainID = alignment.annotations["id"]
        except (AttributeError, KeyError):
            line = f"chain {score:g} {tName} {tSize} {tStrand} {tStart} {tEnd} {qName} {qSize} {qStrand} {qStart} {qEnd}"
        else:
            line = f"chain {score:g} {tName} {tSize} {tStrand} {tStart} {tEnd} {qName} {qSize} {qStrand} {qStart} {qEnd} {chainID}"
        lines.append(line)
        # variable names follow those in the UCSC chain file format description
        tStart, qStart = coordinates[:, 0]
        step = 0
        qGap = 0
        tGap = 0
        for tEnd, qEnd in coordinates[:, 1:].transpose():
            if tStart == tEnd:
                if tGap > 0:
                    qGap = qEnd - qStart
                    line = f"{step}\t{tGap}\t{qGap}"
                    lines.append(line)
                    step = 0
                    tGap = 0
                    qGap = 0
                elif qGap > 0:
                    line = f"{step}\t{tGap}\t{qGap}"
                    lines.append(line)
                    step = 0
                    qGap = qEnd - qStart
                else:
                    qGap = qEnd - qStart
                qStart = qEnd
            elif qStart == qEnd:
                if qGap > 0:
                    tGap = tEnd - tStart
                    line = f"{step}\t{tGap}\t{qGap}"
                    lines.append(line)
                    step = 0
                    tGap = 0
                    qGap = 0
                elif tGap > 0:
                    line = f"{step}\t{tGap}\t{qGap}"
                    lines.append(line)
                    step = 0
                    tGap = tEnd - tStart
                else:
                    tGap = tEnd - tStart
                tStart = tEnd
            else:
                if step == 0:
                    tStep = tEnd - tStart
                    qStep = qEnd - qStart
                    if tGap == 0 and qGap == 0:
                        step = tStep
                    else:
                        line = f"{tStep}\t{tGap}\t{qGap}"
                        lines.append(line)
                        tGap = 0
                        qGap = 0
                else:
                    line = f"{step}\t{tGap}\t{qGap}"
                    lines.append(line)
                    tStep = tEnd - tStart
                    qStep = qEnd - qStart
                    step = tStep
                    tGap = 0
                    qGap = 0
                if tStep != qStep:
                    raise ValueError(
                        f"Expected equal step size in target and query (found {tStep} and {qStep})"
                    )
                tStart = tEnd
                qStart = qEnd
        if tGap > 0 or qGap > 0:
            line = f"0\t{tGap}\t{qGap}"
            lines.append(line)
        line = f"{step}"
        lines.append(line)
        line = ""
        lines.append(line)
        block = "\n".join(lines) + "\n"
        return block


class AlignmentIterator(interfaces.AlignmentIterator):
    """Alignment iterator for UCSC chain files.

    Each chain block in the file contains one pairwise alignment, which are
    loaded and returned incrementally.  The alignment score is scored as an
    attribute of the alignment; the ID is stored under the key "id" in the
    dictionary referred to by the annotations attribute of the alignment.
    """

    fmt = "chain"

    def _read_header(self, stream):
        try:
            self._line = next(stream)
        except StopIteration:
            self._line = None

    def _read_next_alignment(self, stream):
        if self._line is None:
            return
        line = self._line
        self._line = None
        words = line.split()
        if len(words) == 12:
            chainID = None
        elif len(words) == 13:
            chainID = words[12]
        else:
            raise ValueError(
                "expected line with 12 or 13 columns at the beginning of a chain block"
            )
        if words[0] != "chain":
            raise ValueError(
                "expected first line of chain block to start with 'chain '"
            )
        score = float(words[1])
        tName = words[2]
        tSize = int(words[3])
        tStrand = words[4]
        tStart = int(words[5])
        tEnd = int(words[6])
        qName = words[7]
        qSize = int(words[8])
        qStrand = words[9]
        qStart = int(words[10])
        qEnd = int(words[11])
        qPosition = 0
        tPosition = 0
        qStarts = [qPosition]
        tStarts = [tPosition]
        for line in stream:
            words = line.split()
            if words[0] == "chain":
                self._line = line
                break
            step = int(words[0])
            qPosition += step
            tPosition += step
            qStarts.append(qPosition)
            tStarts.append(tPosition)
            if len(words) == 1:
                break
            tGap = int(words[1])
            qGap = int(words[2])
            if tGap > 0:
                tPosition += tGap
                qStarts.append(qPosition)
                tStarts.append(tPosition)
            if qGap > 0:
                qPosition += qGap
                qStarts.append(qPosition)
                tStarts.append(tPosition)
        target_sequence = Seq(None, length=tSize)
        query_sequence = Seq(None, length=qSize)
        target_record = SeqRecord(target_sequence, id=tName, description="")
        query_record = SeqRecord(query_sequence, id=qName, description="")
        records = [target_record, query_record]
        coordinates = np.array([tStarts, qStarts])
        coordinates[0, :] += tStart
        coordinates[1, :] += qStart
        if tStrand == "+":
            pass
        elif tStrand == "-":
            coordinates[0, :] = tSize - coordinates[0, :]
        else:
            raise ValueError(f"tStrand should be '+' or '-' (found '{tStrand}')")
        if qStrand == "+":
            pass
        elif qStrand == "-":
            coordinates[1, :] = qSize - coordinates[1, :]
        else:
            raise ValueError(f"qStrand should be '+' or '-' (found '{qStrand}')")
        alignment = Alignment(records, coordinates)
        alignment.score = score
        if chainID is not None:
            alignment.annotations = {"id": chainID}
        # There is supposed to be a blank line between chain blocks, but some
        # tools (e.g. pslToChain) do not include such a blank line in their
        # output.
        try:
            line = next(stream)
            if line.strip() == "":
                line = next(stream)
            self._line = line
        except StopIteration:
            self._line = None
        return alignment
