# Copyright 2022 by Michiel de Hoon.
# All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.Align support module (not for general use).

This module defines the Alignment class, the Alignments, LazyAlignments, and
ParsedAlignments class, and the AlignmentWriter class. Only the Alignment and
Alignments classes are made available to the user via Bio.Align; the other
classes are not intended to be directly used by a user.

The Lazyalignments and ParsedAlignments classes are abstract base classes
derived from the Alignments class. Concrete subclasses are implemented in the
alignment file parser modules (the AlignmentIterator class) and in the pairwise
alignment module (the PairwiseAlignments class). The inheritance relations are
shown in this diagram:

                                       .- ParsedAlignments <- AlignmentIterator
                                       |                      (in file parser
list <- Alignments <- LazyAlignments <-|                       modules)
                                       |
                                       .- PairwiseAlignments
                                          (in the pairwise
                                           alignment module)

AlignmentWriter is also an abstract base class, with concrete subclasses
implemented in the file parser modules)>

Unless you are writing a new parser or writer for Bio.Align, you should not
use this module directly.
"""

import sys
import warnings
import numbers
from itertools import zip_longest
from abc import ABC, abstractmethod

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Please install numpy if you want to use Bio.Align. "
        "See http://www.numpy.org/"
    ) from None

from Bio import StreamModeError
from Bio import BiopythonDeprecationWarning
from Bio.Align import substitution_matrices
from Bio.Seq import Seq, reverse_complement, UndefinedSequenceError


class Alignment:
    """An Alignment object represents a sequence alignment.

    An Alignment object has a `.sequences` attribute storing the sequences
    (Seq, MutableSeq, SeqRecord, or string objects) that were aligned, as well
    as a `.coordinates` attribute storing the sequence coordinates defining the
    alignment as a numpy array.

    Other commonly used attributes (which may or may not be present) are:
         - annotations        - A dictionary with annotations describing the
                                alignment;
         - column_annotations - A dictionary with annotations describing each
                                column in the alignment;
         - score              - The alignment score.
    """

    @classmethod
    def infer_coordinates(cls, lines, skipped_columns=None):
        """Infer the coordinates from a printed alignment.

        This method is primarily employed in Biopython's alignment parsers,
        though it may be useful for other purposes.

        For an alignment consisting of N sequences, printed as N lines with
        the same number of columns, where gaps are represented by dashes,
        this method will calculate the sequence coordinates that define the
        alignment. The coordinates are returned as a numpy array of integers,
        and can be used to create an Alignment object.

        The argument skipped columns should be None (the default) or an empty
        list. If skipped_columns is a list, then the indices of any columns in
        the alignment with a gap in all lines are appended to skipped_columns.

        This is an example for the alignment of three sequences TAGGCATACGTG,
        AACGTACGT, and ACGCATACTTG, with gaps in the second and third sequence:

        >>> from Bio.Align import Alignment
        >>> lines = ["TAGGCATACGTG",
        ...          "AACG--TACGT-",
        ...          "-ACGCATACTTG",
        ...         ]
        >>> sequences = [line.replace("-", "") for line in lines]
        >>> sequences
        ['TAGGCATACGTG', 'AACGTACGT', 'ACGCATACTTG']
        >>> coordinates = Alignment.infer_coordinates(lines)
        >>> coordinates
        array([[ 0,  1,  4,  6, 11, 12],
               [ 0,  1,  4,  4,  9,  9],
               [ 0,  0,  3,  5, 10, 11]])
        >>> alignment = Alignment(sequences, coordinates)
        """
        n = len(lines)
        m = len(lines[0])
        for line in lines:
            assert m == len(line)
        path = []
        if m > 0:
            indices = [0] * n
            current_state = [None] * n
            for i in range(m):
                next_state = [line[i] != "-" for line in lines]
                if not any(next_state):
                    # skip columns in which all rows have a gap
                    if skipped_columns is not None:
                        skipped_columns.append(i)
                elif next_state == current_state:
                    step += 1  # noqa: F821
                else:
                    indices = [
                        index + step if state else index
                        for index, state in zip(indices, current_state)
                    ]
                    path.append(indices)
                    step = 1
                    current_state = next_state
            indices = [
                index + step if state else index
                for index, state in zip(indices, current_state)
            ]
            path.append(indices)
        coordinates = numpy.array(path).transpose()
        return coordinates

    def __init__(self, sequences, coordinates=None):
        """Initialize a new Alignment object.

        Arguments:
         - sequences   - A list of the sequences (Seq, MutableSeq, SeqRecord,
                         or string objects)that were aligned.
         - coordinates - The sequence coordinates that define the alignment.
                         If None (the default value), assume that the sequences
                         align to each other without any gaps.
        """
        self.sequences = sequences
        if coordinates is None:
            try:
                lengths = {len(sequence) for sequence in sequences}
            except TypeError:
                # this may happen if sequences contain a SeqRecord where
                # the seq attribute is None, as neither the sequence nor
                # its length are known.
                pass
            else:
                if len(lengths) != 1:
                    raise ValueError(
                        "sequences must have the same length if coordinates is None"
                    )
                length = lengths.pop()
                coordinates = numpy.array([[0, length]] * len(sequences))
        self.coordinates = coordinates

    @property
    def target(self):
        """Return self.sequences[0] for a pairwise alignment."""
        n = len(self.sequences)
        if n != 2:
            raise ValueError(
                "self.target is defined for pairwise alignments only (found alignment of %d sequences)"
                % n
            )
        return self.sequences[0]

    @target.setter
    def target(self, value):
        """For a pairwise alignment, set self.sequences[0]."""
        n = len(self.sequences)
        if n != 2:
            raise ValueError(
                "self.target is defined for pairwise alignments only (found alignment of %d sequences)"
                % n
            )
        self.sequences[0] = value

    @property
    def query(self):
        """Return self.sequences[1] for a pairwise alignment."""
        n = len(self.sequences)
        if n != 2:
            raise ValueError(
                "self.query is defined for pairwise alignments only (found alignment of %d sequences)"
                % n
            )
        return self.sequences[1]

    @query.setter
    def query(self, value):
        """For a pairwise alignment, set self.sequences[1]."""
        n = len(self.sequences)
        if n != 2:
            raise ValueError(
                "self.query is defined for pairwise alignments only (found alignment of %d sequences)"
                % n
            )
        self.sequences[1] = value

    def __eq__(self, other):
        """Check if two Alignment objects specify the same alignment."""
        for left, right in zip_longest(self.sequences, other.sequences):
            try:
                left = left.seq
            except AttributeError:
                pass
            try:
                right = right.seq
            except AttributeError:
                pass
            if left != right:
                return False
        return numpy.array_equal(self.coordinates, other.coordinates)

    def __ne__(self, other):
        """Check if two Alignment objects have different alignments."""
        for left, right in zip_longest(self.sequences, other.sequences):
            try:
                left = left.seq
            except AttributeError:
                pass
            try:
                right = right.seq
            except AttributeError:
                pass
            if left != right:
                return True

        return not numpy.array_equal(self.coordinates, other.coordinates)

    def __lt__(self, other):
        """Check if self should come before other."""
        for left, right in zip_longest(self.sequences, other.sequences):
            try:
                left = left.seq
            except AttributeError:
                pass
            try:
                right = right.seq
            except AttributeError:
                pass
            if left < right:
                return True
            if left > right:
                return False
        for left, right in zip(
            self.coordinates.transpose(), other.coordinates.transpose()
        ):
            left, right = tuple(left), tuple(right)
            if left < right:
                return True
            if left > right:
                return False
        return False

    def __le__(self, other):
        """Check if self should come before or is equal to other."""
        for left, right in zip_longest(self.sequences, other.sequences):
            try:
                left = left.seq
            except AttributeError:
                pass
            try:
                right = right.seq
            except AttributeError:
                pass
            if left < right:
                return True
            if left > right:
                return False
        for left, right in zip(
            self.coordinates.transpose(), other.coordinates.transpose()
        ):
            left, right = tuple(left), tuple(right)
            if left < right:
                return True
            if left > right:
                return False
        return True

    def __gt__(self, other):
        """Check if self should come after other."""
        for left, right in zip_longest(self.sequences, other.sequences):
            try:
                left = left.seq
            except AttributeError:
                pass
            try:
                right = right.seq
            except AttributeError:
                pass
            if left < right:
                return False
            if left > right:
                return True
        for left, right in zip(
            self.coordinates.transpose(), other.coordinates.transpose()
        ):
            left, right = tuple(left), tuple(right)
            if left > right:
                return True
            if left < right:
                return False
        return False

    def __ge__(self, other):
        """Check if self should come after or is equal to other."""
        for left, right in zip_longest(self.sequences, other.sequences):
            try:
                left = left.seq
            except AttributeError:
                pass
            try:
                right = right.seq
            except AttributeError:
                pass
            if left < right:
                return False
            if left > right:
                return True
        for left, right in zip(
            self.coordinates.transpose(), other.coordinates.transpose()
        ):
            left, right = tuple(left), tuple(right)
            if left > right:
                return True
            if left < right:
                return False
        return True

    @property
    def path(self):
        """Return the path through the trace matrix."""
        warnings.warn(
            "The path attribute is deprecated; please use the coordinates "
            "attribute instead. The coordinates attribute is a numpy array "
            "containing the same values as the path attributes, after "
            "transposition.",
            BiopythonDeprecationWarning,
        )
        return tuple(tuple(row) for row in self.coordinates.transpose())

    @path.setter
    def path(self, value):
        warnings.warn(
            "The path attribute is deprecated; please use the coordinates "
            "attribute instead. The coordinates attribute is a numpy array "
            "containing the same values as the path attributes, after "
            "transposition.",
            BiopythonDeprecationWarning,
        )
        self.coordinates = numpy.array(value).transpose()

    def _get_row(self, index):
        """Return self[index], where index is an integer (PRIVATE).

        This method is called by __getitem__ for invocations of the form

        self[row]

        where row is an integer.
        """
        coordinates = self.coordinates.copy()
        for sequence, row in zip(self.sequences, coordinates):
            if row[0] > row[-1]:  # reverse strand
                row[:] = len(sequence) - row[:]
        steps = numpy.diff(coordinates, 1)
        gaps = steps.max(0)
        if not ((steps == gaps) | (steps == 0)).all():
            raise ValueError("Unequal step sizes in alignment")
        sequence = self.sequences[index]
        try:
            sequence = sequence.seq  # SeqRecord confusion
        except AttributeError:
            pass
        coordinates = coordinates[index]
        if self.coordinates[index, 0] > self.coordinates[index, -1]:
            # reverse strand
            sequence = reverse_complement(sequence, inplace=False)
        line = ""
        steps = steps[index]
        i = coordinates[0]
        for step, gap in zip(steps, gaps):
            if step:
                j = i + step
                line += str(sequence[i:j])
                i = j
            else:
                line += "-" * gap
        return line

    def _get_rows(self, key):
        """Return self[key], where key is a slice object (PRIVATE).

        This method is called by __getitem__ for invocations of the form

        self[rows]

        where rows is a slice object.
        """
        sequences = self.sequences[key]
        coordinates = self.coordinates[key].copy()
        alignment = Alignment(sequences, coordinates)
        if numpy.array_equal(self.coordinates, coordinates):
            try:
                alignment.score = self.score
            except AttributeError:
                pass
            try:
                alignment.column_annotations = self.column_annotations
            except AttributeError:
                pass
        return alignment

    def _get_row_col(self, j, col, steps, gaps, sequence):
        """Return the sequence contents at alignment column j (PRIVATE).

        This method is called by __getitem__ for invocations of the form

        self[row, col]

        where both row and col are integers.
        """
        indices = gaps.cumsum()
        index = indices.searchsorted(col, side="right")
        if steps[index]:
            offset = col - indices[index]
            j += sum(steps[: index + 1]) + offset
            return sequence[j]
        else:
            return "-"

    def _get_row_cols_slice(
        self, coordinate, start_index, stop_index, steps, gaps, sequence
    ):
        """Return the alignment contents of one row and consecutive columns (PRIVATE).

        This method is called by __getitem__ for invocations of the form

        self[row, cols]

        where row is an integer and cols is a slice object with step 1.
        """
        sequence_indices = coordinate + steps.cumsum()
        indices = gaps.cumsum()
        i = indices.searchsorted(start_index, side="right")
        j = i + indices[i:].searchsorted(stop_index, side="right")
        if i == j:
            length = stop_index - start_index
            if steps[i] == 0:
                line = "-" * length
            else:
                offset = start_index - indices[i]
                start = sequence_indices[i] + offset
                stop = start + length
                line = str(sequence[start:stop])
        else:
            length = indices[i] - start_index
            stop = sequence_indices[i]
            if steps[i] == 0:
                line = "-" * length
            else:
                start = stop - length
                line = str(sequence[start:stop])
            i += 1
            while i < j:
                step = gaps[i]
                if steps[i] == 0:
                    line += "-" * step
                else:
                    start = stop
                    stop = start + step
                    line += str(sequence[start:stop])
                i += 1
            length = stop_index - indices[j - 1]
            if length > 0:
                if steps[j] == 0:
                    line += "-" * length
                else:
                    start = stop
                    stop = start + length
                    line += str(sequence[start:stop])
        return line

    def _get_row_cols_iterable(self, i, cols, steps, gaps, sequence):
        """Return the alignment contents of one row and multiple columns (PRIVATE).

        This method is called by __getitem__ for invocations of the form

        self[row, cols]

        where row is an integer and cols is an iterable of integers.
        """
        line = ""
        for step, gap in zip(steps, gaps):
            if step:
                j = i + step
                line += str(sequence[i:j])
                i = j
            else:
                line += "-" * gap
        try:
            line = "".join(line[col] for col in cols)
        except IndexError:
            raise
        except Exception:
            raise TypeError(
                "second index must be an integer, slice, or iterable of integers"
            ) from None
        return line

    def _get_rows_col(self, coordinates, col, steps, gaps, sequences):
        """Return the alignment contents of multiple rows and one column (PRIVATE).

        This method is called by __getitem__ for invocations of the form

        self[rows, col]

        where rows is a slice object, and col is an integer.
        """
        indices = gaps.cumsum()
        j = indices.searchsorted(col, side="right")
        offset = indices[j] - col
        line = ""
        for sequence, coordinate, step in zip(sequences, coordinates, steps):
            if step[j] == 0:
                line += "-"
            else:
                index = coordinate[j] + step[j] - offset
                line += sequence[index]
        return line

    def _get_rows_cols_slice(
        self, coordinates, row, start_index, stop_index, steps, gaps, sequences
    ):
        """Return a subalignment of multiple rows and consecutive columns (PRIVATE).

        This method is called by __getitem__ for invocations of the form

        self[rows, cols]

        where rows is an arbitrary slice object, and cols is a slice object
        with step 1, allowing the alignment sequences to be reused in the
        subalignment.
        """
        indices = gaps.cumsum()
        i = indices.searchsorted(start_index, side="right")
        j = i + indices[i:].searchsorted(stop_index, side="left") + 1
        offset = steps[:, i] - indices[i] + start_index
        coordinates[:, i] += offset * (steps[:, i] > 0)
        offset = indices[j - 1] - stop_index
        coordinates[:, j] -= offset * (steps[:, j - 1] > 0)
        coordinates = coordinates[:, i : j + 1]
        reverse_complemented = self.coordinates[row, 0] > self.coordinates[row, -1]
        for i, sequence in enumerate(sequences):
            if reverse_complemented[i]:
                # mapped to reverse strand
                coordinates[i, :] = len(sequence) - coordinates[i, :]
        alignment = Alignment(self.sequences[row], coordinates)
        if numpy.array_equal(self.coordinates, coordinates):
            try:
                alignment.score = self.score
            except AttributeError:
                pass
        try:
            column_annotations = self.column_annotations
        except AttributeError:
            pass
        else:
            alignment.column_annotations = {}
            for key, value in column_annotations.items():
                value = value[start_index:stop_index]
                try:
                    value = value.copy()
                except AttributeError:
                    # immutable tuples like str, tuple
                    pass
                alignment.column_annotations[key] = value
        return alignment

    def _get_rows_cols_iterable(self, coordinates, col, steps, gaps, sequences):
        """Return a subalignment of multiple rows and columns (PRIVATE).

        This method is called by __getitem__ for invocations of the form

        self[rows, cols]

        where rows is a slice object and cols is an iterable of integers.
        This method will create new sequences for use by the subalignment
        object.
        """
        indices = tuple(col)
        lines = []
        for i, sequence in enumerate(sequences):
            line = ""
            k = coordinates[i, 0]
            for step, gap in zip(steps[i], gaps):
                if step:
                    j = k + step
                    line += str(sequence[k:j])
                    k = j
                else:
                    line += "-" * gap
            try:
                line = "".join(line[index] for index in indices)
            except IndexError:
                raise
            except Exception:
                raise TypeError(
                    "second index must be an integer, slice, or iterable of integers"
                ) from None
            lines.append(line)
        sequences = [line.replace("-", "") for line in lines]
        coordinates = self.infer_coordinates(lines)
        alignment = Alignment(sequences, coordinates)
        try:
            column_annotations = self.column_annotations
        except AttributeError:
            pass
        else:
            alignment.column_annotations = {}
            for key, value in column_annotations.items():
                values = (value[index] for index in indices)
                if isinstance(value, str):
                    value = "".join(values)
                else:
                    value = value.__class__(values)
                alignment.column_annotations[key] = value
        return alignment

    def __getitem__(self, key):
        """Return self[key].

        Indices of the form

        self[:, :]

        return a copy of the Alignment object;

        self[:, i:]
        self[:, :j]
        self[:, i:j]
        self[:, iterable] (where iterable returns integers)

        return a new Alignment object spanning the selected columns;

        self[k, i]
        self[k, i:]
        self[k, :j]
        self[k, i:j]
        self[k, iterable] (where iterable returns integers)
        self[k] (equivalent to self[k, :])

        return a string with the aligned sequence (including gaps) for the
        selected columns, where k = 0 represents the target and k = 1
        represents the query sequence; and

        self[:, i]

        returns a string with the selected column in the alignment.

        >>> from Bio.Align import PairwiseAligner
        >>> aligner = PairwiseAligner()
        >>> alignments = aligner.align("ACCGGTTT", "ACGGGTT")
        >>> alignment = alignments[0]
        >>> print(alignment)
        ACCGG-TTT
        ||-||-||-
        AC-GGGTT-
        <BLANKLINE>
        >>> alignment[0, :]
        'ACCGG-TTT'
        >>> alignment[1, :]
        'AC-GGGTT-'
        >>> alignment[0]
        'ACCGG-TTT'
        >>> alignment[1]
        'AC-GGGTT-'
        >>> alignment[0, 1:-2]
        'CCGG-T'
        >>> alignment[1, 1:-2]
        'C-GGGT'
        >>> alignment[0, (1, 5, 2)]
        'C-C'
        >>> alignment[1, ::2]
        'A-GT-'
        >>> alignment[1, range(0, 9, 2)]
        'A-GT-'
        >>> alignment[:, 0]
        'AA'
        >>> alignment[:, 5]
        '-G'
        >>> alignment[:, 1:]  # doctest:+ELLIPSIS
        <Alignment object (2 rows x 8 columns) at 0x...>
        >>> print(alignment[:, 1:])
        ACCGG-TTT
         |-||-||-
        AC-GGGTT-
        <BLANKLINE>
        >>> print(alignment[:, 2:])
        ACCGG-TTT
          -||-||-
        AC-GGGTT-
        <BLANKLINE>
        >>> print(alignment[:, 3:])
        ACCGG-TTT
           ||-||-
         ACGGGTT-
        <BLANKLINE>
        >>> print(alignment[:, 3:-1])
        ACCGG-TTT
           ||-||
         ACGGGTT
        <BLANKLINE>
        >>> print(alignment[:, ::2])
        ACGTT
        |-||-
        A-GT-
        <BLANKLINE>
        >>> print(alignment[:, range(1, 9, 2)])
        CG-T
        ||-|
        CGGT
        <BLANKLINE>
        >>> print(alignment[:, (2, 7, 3)])
        CTG
        -||
        -TG
        <BLANKLINE>
        """
        if isinstance(key, numbers.Integral):
            return self._get_row(key)
        if isinstance(key, slice):
            return self._get_rows(key)
        sequences = list(self.sequences)
        coordinates = self.coordinates.copy()
        for i, sequence in enumerate(sequences):
            try:
                sequence = sequence.seq  # SeqRecord confusion
            except AttributeError:
                pass
            if coordinates[i, 0] > coordinates[i, -1]:  # reverse strand
                coordinates[i, :] = len(sequence) - coordinates[i, :]
                sequence = reverse_complement(sequence, inplace=False)
            sequences[i] = sequence
        steps = numpy.diff(coordinates, 1)
        gaps = steps.max(0)
        if not ((steps == gaps) | (steps == 0)).all():
            raise ValueError("Unequal step sizes in alignment")
        m = sum(gaps)
        if isinstance(key, tuple):
            try:
                row, col = key
            except ValueError:
                raise ValueError("only tuples of length 2 can be alignment indices")
        else:
            raise TypeError("alignment indices must be integers, slices, or tuples")
        if isinstance(col, numbers.Integral):
            if col < 0:
                col += m
            if col < 0 or col >= m:
                raise IndexError(
                    "column index %d is out of bounds (%d columns)" % (col, m)
                )
        steps = steps[row]
        if isinstance(row, numbers.Integral):
            sequence = sequences[row]
            coordinate = coordinates[row, 0]
            if isinstance(col, numbers.Integral):
                return self._get_row_col(coordinate, col, steps, gaps, sequence)
            if isinstance(col, slice):
                start_index, stop_index, step = col.indices(m)
                if start_index < stop_index and step == 1:
                    return self._get_row_cols_slice(
                        coordinate, start_index, stop_index, steps, gaps, sequence
                    )
                # make an iterable if step != 1
                col = range(start_index, stop_index, step)
            return self._get_row_cols_iterable(coordinate, col, steps, gaps, sequence)
        if isinstance(row, slice):
            sequences = sequences[row]
            coordinates = coordinates[row]
            if isinstance(col, numbers.Integral):
                return self._get_rows_col(coordinates, col, steps, gaps, sequences)
            if isinstance(col, slice):
                start_index, stop_index, step = col.indices(m)
                if start_index < stop_index and step == 1:
                    return self._get_rows_cols_slice(
                        coordinates,
                        row,
                        start_index,
                        stop_index,
                        steps,
                        gaps,
                        sequences,
                    )
                # make an iterable if step != 1
                col = range(start_index, stop_index, step)
            # try if we can use col as an iterable
            return self._get_rows_cols_iterable(
                coordinates, col, steps, gaps, sequences
            )
        raise TypeError("first index must be an integer or slice")

    def _convert_sequence_string(self, sequence):
        """Convert given sequence to string using the appropriate method (PRIVATE)."""
        if isinstance(sequence, (bytes, bytearray)):
            return sequence.decode()
        if isinstance(sequence, str):
            return sequence
        if isinstance(sequence, Seq):
            return str(sequence)
        try:  # check if target is a SeqRecord
            sequence = sequence.seq
        except AttributeError:
            pass
        else:
            return str(sequence)
        try:
            view = memoryview(sequence)
        except TypeError:
            pass
        else:
            if view.format == "c":
                return str(sequence)
        return None

    def __format__(self, format_spec):
        """Return the alignment as a string in the specified file format.

        Wrapper for self.format() .
        """
        return self.format(format_spec)

    def format(self, fmt="", **kwargs):
        """Return the alignment as a string in the specified file format.

        Arguments:
         - fmt       - File format. Acceptable values are
                       ""   : create a human-readable representation of the
                              alignment (default);
                       "bed": create a line representing the alignment in
                              the Browser Extensible Data (BED) file format;
                       "psl": create a line representing the alignment in
                              the Pattern Space Layout (PSL) file format as
                              generated by BLAT;
                       "sam": create a line representing the alignment in
                              the Sequence Alignment/Map (SAM) format.
         - mask      - PSL format only. Specify if repeat regions in the target
                       sequence are masked and should be reported in the
                       `repMatches` field of the PSL file instead of in the
                       `matches` field. Acceptable values are
                       None   : no masking (default);
                       "lower": masking by lower-case characters;
                       "upper": masking by upper-case characters.
         - wildcard  - PSL format only. Report alignments to the wildcard
                       character in the target or query sequence in the
                       `nCount` field of the PSL file instead of in the
                       `matches`, `misMatches`, or `repMatches` fields.
                       Default value is 'N'.
        """
        if len(self.sequences) > 2:
            raise NotImplementedError(
                "format is currently implemented for pairwise alignments only"
            )
        if fmt == "":
            return self._format_pretty(**kwargs)
        elif fmt == "psl":
            return self._format_psl(**kwargs)
        elif fmt == "bed":
            return self._format_bed(**kwargs)
        elif fmt == "sam":
            return self._format_sam(**kwargs)
        else:
            raise ValueError("Unknown format %s" % fmt)

    def _format_pretty(self):
        """Return default string representation (PRIVATE).

        Helper for self.format() .
        """
        target, query = self.sequences
        seq1 = self._convert_sequence_string(target)
        if seq1 is None:
            return self._format_generalized()
        seq2 = self._convert_sequence_string(query)
        if seq2 is None:
            return self._format_generalized()
        n1 = len(seq1)
        n2 = len(seq2)
        aligned_seq1 = ""
        aligned_seq2 = ""
        pattern = ""
        coordinates = self.coordinates
        if coordinates[0, 0] > coordinates[0, -1]:  # mapped to reverse strand
            coordinates = coordinates.copy()
            coordinates[0, :] = n1 - coordinates[0, :]
            seq1 = reverse_complement(seq1, inplace=False)
        if coordinates[1, 0] > coordinates[1, -1]:  # mapped to reverse strand
            coordinates = coordinates.copy()
            coordinates[1, :] = n2 - coordinates[1, :]
            seq2 = reverse_complement(seq2, inplace=False)
        coordinates = coordinates.transpose()
        end1, end2 = coordinates[0, :]
        if end1 > 0 or end2 > 0:
            end = max(end1, end2)
            aligned_seq1 += " " * (end - end1) + seq1[:end1]
            aligned_seq2 += " " * (end - end2) + seq2[:end2]
            pattern += " " * end
        start1 = end1
        start2 = end2
        for end1, end2 in coordinates[1:]:
            if end1 == start1:
                gap = end2 - start2
                aligned_seq1 += "-" * gap
                aligned_seq2 += seq2[start2:end2]
                pattern += "-" * gap
            elif end2 == start2:
                gap = end1 - start1
                aligned_seq1 += seq1[start1:end1]
                aligned_seq2 += "-" * gap
                pattern += "-" * gap
            else:
                s1 = seq1[start1:end1]
                s2 = seq2[start2:end2]
                if len(s1) != len(s2):
                    raise ValueError("Unequal step sizes in alignment")
                aligned_seq1 += s1
                aligned_seq2 += s2
                for c1, c2 in zip(s1, s2):
                    if c1 == c2:
                        pattern += "|"
                    else:
                        pattern += "."
            start1 = end1
            start2 = end2
        aligned_seq1 += seq1[end1:]
        aligned_seq2 += seq2[end2:]
        return f"{aligned_seq1}\n{pattern}\n{aligned_seq2}\n"

    def _format_generalized(self):
        """Return generalized string representation (PRIVATE).

        Helper for self._format_pretty() .
        """
        seq1, seq2 = self.sequences
        aligned_seq1 = []
        aligned_seq2 = []
        pattern = []
        end1, end2 = self.coordinates[:, 0]
        if end1 > 0 or end2 > 0:
            if end1 <= end2:
                for c2 in seq2[: end2 - end1]:
                    s2 = str(c2)
                    s1 = " " * len(s2)
                    aligned_seq1.append(s1)
                    aligned_seq2.append(s2)
                    pattern.append(s1)
            else:  # end1 > end2
                for c1 in seq1[: end1 - end2]:
                    s1 = str(c1)
                    s2 = " " * len(s1)
                    aligned_seq1.append(s1)
                    aligned_seq2.append(s2)
                    pattern.append(s2)
        start1 = end1
        start2 = end2
        for end1, end2 in self.coordinates[:, 1:].transpose():
            if end1 == start1:
                for c2 in seq2[start2:end2]:
                    s2 = str(c2)
                    s1 = "-" * len(s2)
                    aligned_seq1.append(s1)
                    aligned_seq2.append(s2)
                    pattern.append(s1)
                start2 = end2
            elif end2 == start2:
                for c1 in seq1[start1:end1]:
                    s1 = str(c1)
                    s2 = "-" * len(s1)
                    aligned_seq1.append(s1)
                    aligned_seq2.append(s2)
                    pattern.append(s2)
                start1 = end1
            else:
                t1 = seq1[start1:end1]
                t2 = seq2[start2:end2]
                if len(t1) != len(t2):
                    raise ValueError("Unequal step sizes in alignment")
                for c1, c2 in zip(t1, t2):
                    s1 = str(c1)
                    s2 = str(c2)
                    m1 = len(s1)
                    m2 = len(s2)
                    if c1 == c2:
                        p = "|"
                    else:
                        p = "."
                    if m1 < m2:
                        space = (m2 - m1) * " "
                        s1 += space
                        pattern.append(p * m1 + space)
                    elif m1 > m2:
                        space = (m1 - m2) * " "
                        s2 += space
                        pattern.append(p * m2 + space)
                    else:
                        pattern.append(p * m1)
                    aligned_seq1.append(s1)
                    aligned_seq2.append(s2)
                start1 = end1
                start2 = end2
        aligned_seq1 = " ".join(aligned_seq1)
        aligned_seq2 = " ".join(aligned_seq2)
        pattern = " ".join(pattern)
        return f"{aligned_seq1}\n{pattern}\n{aligned_seq2}\n"

    def _format_bed(self, bedN=12):
        """Return BED file format string representation (PRIVATE).

        Helper for self.format().
        """
        from . import bed

        writer = bed.AlignmentWriter(None, bedN=bedN)
        return writer.format_alignment(self)

    def _format_psl(self, mask=False, wildcard="N"):
        """Return PSL file format string representation (PRIVATE).

        Helper for self.format() .
        """
        from . import psl

        writer = psl.AlignmentWriter(None, header=False, mask=mask, wildcard=wildcard)
        return writer.format_alignment(self)

    def _format_sam(self):
        """Return SAM file format string representation (PRIVATE).

        Helper for self.format() .
        """
        from . import sam

        writer = sam.AlignmentWriter(None)
        return writer.format_alignment(self)

    def __str__(self):
        """Return a string representation of the Alignment object.

        Wrapper for self.format().
        """
        if len(self.sequences) > 2:
            raise NotImplementedError(
                "__str__ is currently implemented for pairwise alignments only"
            )
        return self.format()

    def __repr__(self):
        """Return a representation of the alignment, including its shape.

        The representation cannot be used with eval() to recreate the object,
        which is usually possible with simple python objects.  For example:

        <Alignment object (2 rows x 14 columns) at 0x10403d850>

        The hex string is the memory address of the object and can be used to
        distinguish different Alignment objects.  See help(id) for more
        information.

        >>> import numpy
        >>> from Bio.Align import Alignment
        >>> alignment = Alignment(("ACCGT", "ACGT"),
        ...                       coordinates = numpy.array([[0, 2, 3, 5],
        ...                                                  [0, 2, 2, 4],
        ...                                                 ]))
        >>> print(alignment)
        ACCGT
        ||-||
        AC-GT
        <BLANKLINE>
        >>> alignment  # doctest:+ELLIPSIS
        <Alignment object (2 rows x 5 columns) at 0x...>
        """
        n, m = self.shape
        return "<%s object (%i rows x %i columns) at 0x%x>" % (
            self.__class__.__name__,
            n,
            m,
            id(self),
        )

    def __len__(self):
        """Return the number of sequences in the alignment."""
        return len(self.sequences)

    @property
    def shape(self):
        """Return the shape of the alignment as a tuple of two integer values.

        The first integer value is the number of sequences in the alignment as
        returned by len(alignment), which is always 2 for pairwise alignments.

        The second integer value is the number of columns in the alignment when
        it is printed, and is equal to the sum of the number of matches, number
        of mismatches, and the total length of gaps in the target and query.
        Sequence sections beyond the aligned segment are not included in the
        number of columns.

        For example,

        >>> from Bio import Align
        >>> aligner = Align.PairwiseAligner()
        >>> aligner.mode = "global"
        >>> alignments = aligner.align("GACCTG", "CGATCG")
        >>> alignment = alignments[0]
        >>> print(alignment)
        -GACCT-G
        -||--|-|
        CGA--TCG
        <BLANKLINE>
        >>> len(alignment)
        2
        >>> alignment.shape
        (2, 8)
        >>> aligner.mode = "local"
        >>> alignments = aligner.align("GACCTG", "CGATCG")
        >>> alignment = alignments[0]
        >>> print(alignment)
         GACCT-G
         ||--|-|
        CGA--TCG
        <BLANKLINE>
        >>> len(alignment)
        2
        >>> alignment.shape
        (2, 7)
        """
        coordinates = numpy.array(self.coordinates)
        n = len(coordinates)
        for i in range(n):
            if coordinates[i, 0] > coordinates[i, -1]:  # mapped to reverse strand
                k = len(self.sequences[i])
                coordinates[i, :] = k - coordinates[i, :]
        steps = numpy.diff(coordinates, 1)
        gaps = steps.max(0)
        if not ((steps == gaps) | (steps == 0)).all():
            raise ValueError("Unequal step sizes in alignment")
        m = sum(gaps)
        return (n, m)

    @property
    def aligned(self):
        """Return the indices of subsequences aligned to each other.

        This property returns the start and end indices of subsequences
        in the target and query sequence that were aligned to each other.
        If the alignment between target (t) and query (q) consists of N
        chunks, you get two tuples of length N:

            (((t_start1, t_end1), (t_start2, t_end2), ..., (t_startN, t_endN)),
             ((q_start1, q_end1), (q_start2, q_end2), ..., (q_startN, q_endN)))

        For example,

        >>> from Bio import Align
        >>> aligner = Align.PairwiseAligner()
        >>> alignments = aligner.align("GAACT", "GAT")
        >>> alignment = alignments[0]
        >>> print(alignment)
        GAACT
        ||--|
        GA--T
        <BLANKLINE>
        >>> alignment.aligned
        array([[[0, 2],
                [4, 5]],
        <BLANKLINE>
               [[0, 2],
                [2, 3]]])
        >>> alignment = alignments[1]
        >>> print(alignment)
        GAACT
        |-|-|
        G-A-T
        <BLANKLINE>
        >>> alignment.aligned
        array([[[0, 1],
                [2, 3],
                [4, 5]],
        <BLANKLINE>
               [[0, 1],
                [1, 2],
                [2, 3]]])

        Note that different alignments may have the same subsequences
        aligned to each other. In particular, this may occur if alignments
        differ from each other in terms of their gap placement only:

        >>> aligner.mismatch_score = -10
        >>> alignments = aligner.align("AAACAAA", "AAAGAAA")
        >>> len(alignments)
        2
        >>> print(alignments[0])
        AAAC-AAA
        |||--|||
        AAA-GAAA
        <BLANKLINE>
        >>> alignments[0].aligned
        array([[[0, 3],
                [4, 7]],
        <BLANKLINE>
               [[0, 3],
                [4, 7]]])
        >>> print(alignments[1])
        AAA-CAAA
        |||--|||
        AAAG-AAA
        <BLANKLINE>
        >>> alignments[1].aligned
        array([[[0, 3],
                [4, 7]],
        <BLANKLINE>
               [[0, 3],
                [4, 7]]])

        The property can be used to identify alignments that are identical
        to each other in terms of their aligned sequences.
        """
        if len(self.sequences) > 2:
            raise NotImplementedError(
                "aligned is currently implemented for pairwise alignments only"
            )
        coordinates = self.coordinates.copy()
        for i, sequence in enumerate(self.sequences):
            if coordinates[i, 0] > coordinates[i, -1]:  # mapped to reverse strand
                n = len(sequence)
                coordinates[i, :] = n - coordinates[i, :]
        coordinates = coordinates.transpose()
        steps = numpy.diff(coordinates, axis=0).min(1)
        indices = numpy.flatnonzero(steps)
        starts = coordinates[indices, :]
        ends = coordinates[indices + 1, :]
        segments = numpy.stack([starts, ends], axis=0).transpose()
        for i, sequence in enumerate(self.sequences):
            if self.coordinates[i, 0] > self.coordinates[i, -1]:
                # mapped to reverse strand
                n = len(sequence)
                segments[i, :] = n - segments[i, :]
        return segments

    def sort(self, key=None, reverse=False):
        """Sort the sequences of the alignment in place.

        By default, this sorts the sequences alphabetically using their id
        attribute if available, or by their sequence contents otherwise.
        For example,

        >>> from Bio.Align import PairwiseAligner
        >>> aligner = PairwiseAligner()
        >>> aligner.gap_score = -1
        >>> alignments = aligner.align("AATAA", "AAGAA")
        >>> len(alignments)
        1
        >>> alignment = alignments[0]
        >>> print(alignment)
        AATAA
        ||.||
        AAGAA
        <BLANKLINE>
        >>> alignment.sort()
        >>> print(alignment)
        AAGAA
        ||.||
        AATAA
        <BLANKLINE>

        Alternatively, a key function can be supplied that maps each sequence
        to a sort value.  For example, you could sort on the GC content of each
        sequence.

        >>> from Bio.SeqUtils import GC
        >>> alignment.sort(key=GC)
        >>> print(alignment)
        AATAA
        ||.||
        AAGAA
        <BLANKLINE>

        You can reverse the sort order by passing `reverse=True`:

        >>> alignment.sort(key=GC, reverse=True)
        >>> print(alignment)
        AAGAA
        ||.||
        AATAA
        <BLANKLINE>

        The sequences are now sorted by decreasing GC content value.
        """
        sequences = self.sequences
        if key is None:
            try:
                values = [sequence.id for sequence in sequences]
            except AttributeError:
                values = sequences
        else:
            values = [key(sequence) for sequence in sequences]
        indices = sorted(range(len(sequences)), key=values.__getitem__, reverse=reverse)
        self.sequences = [sequences[index] for index in indices]
        self.coordinates = self.coordinates.take(indices, 0)

    def map(self, alignment):
        r"""Map the alignment to self.target and return the resulting alignment.

        Here, self.query and alignment.target are the same sequence.

        A typical example is where self is the pairwise alignment between a
        chromosome and a transcript, the argument is the pairwise alignment
        between the transcript and a sequence (e.g., as obtained by RNA-seq),
        and we want to find the alignment of the sequence to the chromosome:

        >>> from Bio import Align
        >>> aligner = Align.PairwiseAligner()
        >>> aligner.mode = 'local'
        >>> aligner.open_gap_score = -1
        >>> aligner.extend_gap_score = 0
        >>> chromosome = "AAAAAAAACCCCCCCAAAAAAAAAAAGGGGGGAAAAAAAA"
        >>> transcript = "CCCCCCCGGGGGG"
        >>> alignments1 = aligner.align(chromosome, transcript)
        >>> len(alignments1)
        1
        >>> alignment1 = alignments1[0]
        >>> print(alignment1)
        AAAAAAAACCCCCCCAAAAAAAAAAAGGGGGGAAAAAAAA
                |||||||-----------||||||
                CCCCCCC-----------GGGGGG
        <BLANKLINE>
        >>> sequence = "CCCCGGGG"
        >>> alignments2 = aligner.align(transcript, sequence)
        >>> len(alignments2)
        1
        >>> alignment2 = alignments2[0]
        >>> print(alignment2)
        CCCCCCCGGGGGG
           ||||||||
           CCCCGGGG
        <BLANKLINE>
        >>> alignment = alignment1.map(alignment2)
        >>> print(alignment)
        AAAAAAAACCCCCCCAAAAAAAAAAAGGGGGGAAAAAAAA
                   ||||-----------||||
                   CCCC-----------GGGG
        <BLANKLINE>
        >>> format(alignment, "psl")
        '8\t0\t0\t0\t0\t0\t1\t11\t+\tquery\t8\t0\t8\ttarget\t40\t11\t30\t2\t4,4,\t0,4,\t11,26,\n'

        Mapping the alignment does not depend on the sequence contents. If we
        delete the sequence contents, the same alignment is found in PSL format
        (though we obviously lose the ability to print the sequence alignment):

        >>> alignment1.target = Seq(None, len(alignment1.target))
        >>> alignment1.query = Seq(None, len(alignment1.query))
        >>> alignment2.target = Seq(None, len(alignment2.target))
        >>> alignment2.query = Seq(None, len(alignment2.query))
        >>> alignment = alignment1.map(alignment2)
        >>> format(alignment, "psl")
        '8\t0\t0\t0\t0\t0\t1\t11\t+\tquery\t8\t0\t8\ttarget\t40\t11\t30\t2\t4,4,\t0,4,\t11,26,\n'
        """
        alignment1, alignment2 = self, alignment
        if len(alignment1.query) != len(alignment2.target):
            raise ValueError(
                "length of alignment1 query sequence (%d) != length of alignment2 target sequence (%d)"
                % (len(alignment1.query), len(alignment2.target))
            )
        target = alignment1.target
        query = alignment2.query
        coordinates1 = alignment1.coordinates
        coordinates2 = alignment2.coordinates
        n1 = len(alignment1.query)
        n2 = len(alignment2.query)
        if coordinates1[1, 0] < coordinates1[1, -1]:  # mapped to forward strand
            strand1 = "+"
        else:  # mapped to reverse strand
            strand1 = "-"
        if coordinates2[1, 0] < coordinates2[1, -1]:  # mapped to forward strand
            strand2 = "+"
        else:  # mapped to reverse strand
            strand2 = "-"
        if strand1 == "+":
            if strand2 == "-":  # mapped to reverse strand
                coordinates2 = coordinates2.copy()
                coordinates2[1, :] = n2 - coordinates2[1, :]
        else:  # mapped to reverse strand
            coordinates1 = coordinates1.copy()
            coordinates1[1, :] = n1 - coordinates1[1, :]
            coordinates2 = coordinates2.copy()
            coordinates2[0, :] = n1 - coordinates2[0, ::-1]
            if strand2 == "+":
                coordinates2[1, :] = n2 - coordinates2[1, ::-1]
            else:  # mapped to reverse strand
                coordinates2[1, :] = coordinates2[1, ::-1]
        steps1 = numpy.diff(coordinates1, 1)
        gaps1 = steps1.max(0)
        if not ((steps1 == gaps1) | (steps1 == 0)).all():
            raise ValueError("Unequal step sizes in first alignment")
        steps2 = numpy.diff(coordinates2, 1)
        gaps2 = steps2.max(0)
        if not ((steps2 == gaps2) | (steps2 == 0)).all():
            raise ValueError("Unequal step sizes in second alignment")
        path = []
        tEnd, qEnd = sys.maxsize, sys.maxsize
        coordinates1 = iter(coordinates1.transpose())
        tStart1, qStart1 = sys.maxsize, sys.maxsize
        for tEnd1, qEnd1 in coordinates1:
            if tStart1 < tEnd1 and qStart1 < qEnd1:
                break
            tStart1, qStart1 = tEnd1, qEnd1
        tStart2, qStart2 = sys.maxsize, sys.maxsize
        for tEnd2, qEnd2 in coordinates2.transpose():
            while qStart2 < qEnd2 and tStart2 < tEnd2:
                while True:
                    if tStart2 < qStart1:
                        if tEnd2 < qStart1:
                            size = tEnd2 - tStart2
                        else:
                            size = qStart1 - tStart2
                        break
                    elif tStart2 < qEnd1:
                        offset = tStart2 - qStart1
                        if tEnd2 > qEnd1:
                            size = qEnd1 - tStart2
                        else:
                            size = tEnd2 - tStart2
                        qStart = qStart2
                        tStart = tStart1 + offset
                        if tStart > tEnd and qStart > qEnd:
                            # adding a gap both in target and in query;
                            # add gap to target first:
                            path.append([tStart, qEnd])
                        qEnd = qStart2 + size
                        tEnd = tStart + size
                        path.append([tStart, qStart])
                        path.append([tEnd, qEnd])
                        break
                    tStart1, qStart1 = sys.maxsize, sys.maxsize
                    for tEnd1, qEnd1 in coordinates1:
                        if tStart1 < tEnd1 and qStart1 < qEnd1:
                            break
                        tStart1, qStart1 = tEnd1, qEnd1
                    else:
                        size = qEnd2 - qStart2
                        break
                qStart2 += size
                tStart2 += size
            tStart2, qStart2 = tEnd2, qEnd2
        coordinates = numpy.array(path).transpose()
        if strand1 != strand2:
            coordinates[1, :] = n2 - coordinates[1, :]
        sequences = [target, query]
        alignment = Alignment(sequences, coordinates)
        return alignment

    @property
    def substitutions(self):
        """Return an Array with the number of substitutions of letters in the alignment.

        As an example, consider a sequence alignment of two RNA sequences:

        >>> from Bio.Align import PairwiseAligner
        >>> target = "ATACTTACCTGGCAGGGGAGATACCATGATCACGAAGGTGGTTTTCCCAGGGCGAGGCTTATCCATTGCACTCCGGATGTGCTGACCCCTGCGATTTCCCCAAATGTGGGAAACTCGACTGCATAATTTGTGGTAGTGGGGGACTGCGTTCGCGCTTTCCCCTG"  # human spliceosomal small nuclear RNA U1
        >>> query = "ATACTTACCTGACAGGGGAGGCACCATGATCACACAGGTGGTCCTCCCAGGGCGAGGCTCTTCCATTGCACTGCGGGAGGGTTGACCCCTGCGATTTCCCCAAATGTGGGAAACTCGACTGTATAATTTGTGGTAGTGGGGGACTGCGTTCGCGCTATCCCCCG"  # sea lamprey spliceosomal small RNA U1
        >>> aligner = PairwiseAligner()
        >>> aligner.gap_score = -10
        >>> alignments = aligner.align(target, query)
        >>> len(alignments)
        1
        >>> alignment = alignments[0]
        >>> print(alignment)
        ATACTTACCTGGCAGGGGAGATACCATGATCACGAAGGTGGTTTTCCCAGGGCGAGGCTTATCCATTGCACTCCGGATGTGCTGACCCCTGCGATTTCCCCAAATGTGGGAAACTCGACTGCATAATTTGTGGTAGTGGGGGACTGCGTTCGCGCTTTCCCCTG
        |||||||||||.||||||||..|||||||||||..|||||||..|||||||||||||||..|||||||||||.|||..|.|.|||||||||||||||||||||||||||||||||||||||.||||||||||||||||||||||||||||||||||.|||||.|
        ATACTTACCTGACAGGGGAGGCACCATGATCACACAGGTGGTCCTCCCAGGGCGAGGCTCTTCCATTGCACTGCGGGAGGGTTGACCCCTGCGATTTCCCCAAATGTGGGAAACTCGACTGTATAATTTGTGGTAGTGGGGGACTGCGTTCGCGCTATCCCCCG
        <BLANKLINE>
        >>> m = alignment.substitutions
        >>> print(m)
             A    C    G    T
        A 28.0  1.0  2.0  1.0
        C  0.0 39.0  1.0  2.0
        G  2.0  0.0 45.0  0.0
        T  2.0  5.0  1.0 35.0
        <BLANKLINE>

        Note that the matrix is not symmetric: rows correspond to the target
        sequence, and columns to the query sequence.  For example, the number
        of T's in the target sequence that are aligned to a C in the query
        sequence is

        >>> m['T', 'C']
        5.0

        and the number of C's in the query sequence tat are aligned to a T in
        the query sequence is

        >>> m['C', 'T']
        2.0

        For some applications (for example, to define a scoring matrix from
        the substitution matrix), a symmetric matrix may be preferred, which
        can be calculated as follows:

        >>> m += m.transpose()
        >>> m /= 2.0
        >>> print(m)
             A    C    G    T
        A 28.0  0.5  2.0  1.5
        C  0.5 39.0  0.5  3.5
        G  2.0  0.5 45.0  0.5
        T  1.5  3.5  0.5 35.0
        <BLANKLINE>

        The matrix is now symmetric, with counts divided equally on both sides
        of the diagonal:

        >>> m['C', 'T']
        3.5
        >>> m['T', 'C']
        3.5

        The total number of substitutions between T's and C's in the alignment
        is 3.5 + 3.5 = 7.
        """
        sequences = self.sequences
        letters = set()
        for sequence in sequences:
            try:
                s = set(sequence)
            except UndefinedSequenceError:
                try:
                    sequence = sequence.seq  # SeqRecord confusion
                except AttributeError:
                    pass
                for start, end in sequence.defined_ranges:
                    s = set(sequence[start:end])
                    letters.update(s)
            else:
                letters.update(s)
        letters = "".join(sorted(letters))
        m = substitution_matrices.Array(letters, dims=2)
        n = len(sequences)
        for i1 in range(n):
            sequence1 = sequences[i1]
            coordinates1 = self.coordinates[i1, :]
            if coordinates1[0] > coordinates1[-1]:
                sequence1 = reverse_complement(sequence1)
                coordinates1 = len(sequence1) - coordinates1
            for i2 in range(i1 + 1, n):
                sequence2 = sequences[i2]
                coordinates2 = self.coordinates[i2, :]
                if coordinates2[0] > coordinates2[-1]:
                    sequence2 = reverse_complement(sequence2)
                    coordinates2 = len(sequence2) - coordinates2
                start1, start2 = sys.maxsize, sys.maxsize
                for end1, end2 in zip(coordinates1, coordinates2):
                    if start1 < end1 and start2 < end2:  # aligned
                        segment1 = sequence1[start1:end1]
                        segment2 = sequence2[start2:end2]
                        if len(segment1) != len(segment2):
                            raise ValueError("Unequal step sizes in alignment")
                        for c1, c2 in zip(segment1, segment2):
                            m[c1, c2] += 1.0
                    start1, start2 = end1, end2
        return m


class Alignments(list):
    """A list-like object storing sequence alignments.

    An `Alignments` object can be used as an iterator or as a list-like object.

    This is an example of an `Alignments` object used as an iterator:

    >>> from Bio.Align import emboss
    >>> alignments = emboss.AlignmentIterator("Emboss/needle.txt")
    >>> for alignment in alignments:
    ...     print("******************************")
    ...     print(alignment[0, :30])
    ...     print(alignment[1, :30])
    ...
    ******************************
    KILIVDD----QYGIRILLNEVFNKEGYQT
    -VLLADDHALVRRGFRLMLED--DPEIEIV
    ******************************
    KILIVDDQYGIRILLNEVFNKEGYQTFQAA
    -ILIVDDEANTLASLSRAFRLAGHEATVCD
    ******************************
    -KILIVDDQYGIRILLNEVFNKEGYQTFQA
    LHIVVVDDDPGTCVYIESVFAELGHTCKSF
    ******************************
    KILIVDDQYGIRILLNEVFNKEGYQTFQAA
    -VLLVEDEEALRAAAGDFLETRGYKIMTAR
    ******************************
    KILIVDDQYGIRILLNEVFNKEGYQTFQAA
    TVLLVEDEEGVRKLVRGILSRQGYHVLEAT

    In general, alignments obtained by parsing an alignment file can iterated
    over only once. However, an `Alignments` object is automatically converted
    to a list-like object when needed, for example if a specific alignment is
    accessed by index:

    >>> from Bio.Align import emboss
    >>> alignments = emboss.AlignmentIterator("Emboss/needle.txt")
    >>> alignment = alignments[2]
    >>> print(alignment[0, :30]); print(alignment[1, :30])
    -KILIVDDQYGIRILLNEVFNKEGYQTFQA
    LHIVVVDDDPGTCVYIESVFAELGHTCKSF

    Use `alignments[:]` to get all alignments as a list-like object:

    >>> from Bio.Align import emboss
    >>> alignments = emboss.AlignmentIterator("Emboss/needle.txt")
    >>> alignments = alignments[:]
    >>> len(alignments)
    5
    >>> alignment = alignments[2]
    >>> print(alignment[0, :30]); print(alignment[1, :30])
    -KILIVDDQYGIRILLNEVFNKEGYQTFQA
    LHIVVVDDDPGTCVYIESVFAELGHTCKSF
    >>> for alignment in alignments:
    ...     print("******************************")
    ...     print(alignment[0, :30])
    ...     print(alignment[1, :30])
    ...
    ******************************
    KILIVDD----QYGIRILLNEVFNKEGYQT
    -VLLADDHALVRRGFRLMLED--DPEIEIV
    ******************************
    KILIVDDQYGIRILLNEVFNKEGYQTFQAA
    -ILIVDDEANTLASLSRAFRLAGHEATVCD
    ******************************
    -KILIVDDQYGIRILLNEVFNKEGYQTFQA
    LHIVVVDDDPGTCVYIESVFAELGHTCKSF
    ******************************
    KILIVDDQYGIRILLNEVFNKEGYQTFQAA
    -VLLVEDEEALRAAAGDFLETRGYKIMTAR
    ******************************
    KILIVDDQYGIRILLNEVFNKEGYQTFQAA
    TVLLVEDEEGVRKLVRGILSRQGYHVLEAT

    Note that here, we are using `alignments` as an iterator after converting
    it to a list-like object. Importantly, using `alignments` as an iterator
    and then converting it to a list-like object will lose the alignments that
    were already extracted:

    >>> from Bio.Align import emboss
    >>> alignments = emboss.AlignmentIterator("Emboss/needle.txt")
    >>> alignment = next(alignments)
    >>> print(alignment[0, :30]); print(alignment[1, :30])
    KILIVDD----QYGIRILLNEVFNKEGYQT
    -VLLADDHALVRRGFRLMLED--DPEIEIV
    >>> alignment = next(alignments)
    >>> print(alignment[0, :30]); print(alignment[1, :30])
    KILIVDDQYGIRILLNEVFNKEGYQTFQAA
    -ILIVDDEANTLASLSRAFRLAGHEATVCD
    >>> alignments = alignments[:]
    >>> len(alignments)
    3
    >>> for alignment in alignments:
    ...     print("******************************")
    ...     print(alignment[0, :30])
    ...     print(alignment[1, :30])
    ...
    ******************************
    -KILIVDDQYGIRILLNEVFNKEGYQTFQA
    LHIVVVDDDPGTCVYIESVFAELGHTCKSF
    ******************************
    KILIVDDQYGIRILLNEVFNKEGYQTFQAA
    -VLLVEDEEALRAAAGDFLETRGYKIMTAR
    ******************************
    KILIVDDQYGIRILLNEVFNKEGYQTFQAA
    TVLLVEDEEGVRKLVRGILSRQGYHVLEAT

    """

    def __init__(self):
        """Initialize self."""
        super().__init__()
        self._index = 0

    def __repr__(self):
        """Return repr(self)."""
        return "<Alignments object at %s>" % hex(id(self))

    def __next__(self):
        """Return the next entry."""
        index = self._index
        length = super().__len__()
        if index == length:
            raise StopIteration from None
        self._index += 1
        return super().__getitem__(index)

    def __iter__(self):
        """Iterate over the alignments as Alignment objects."""
        self._index = 0
        return self


class LazyAlignments(Alignments, ABC):
    # The LazyAlignments class is an abstract base class for lazy loading of
    # sequence alignments. This class is a subclass of Alignments, which is a
    # subclass of list.
    #
    # A newly created LazyAlignments object will act as an iterator until a
    # method is used that requires list-like behavior. In that case, the _load
    # method will read in all alignments and store them in the grandparent
    # class list object. From that point on, the LazyAlignments class will act
    # as a list-like object.
    #
    # The PairwiseAlignments class in Bio.Align._pairwise is a concrete subclass
    # of the LazyAlignments class. The ParsedAlignments class is an abstract
    # subclass of LazyAlignments, with concrete subclasses in the alignment file
    # parser modules.

    def _load(self):
        for item in self:
            super().append(item)
        self._index = 0  # for use by the iterator
        self.__class__ = Alignments

    def __lt__(self, other):
        self._load()
        if isinstance(other, LazyAlignments):
            other._load()
        return self.__lt__(other)

    def __le__(self, other):
        self._load()
        if isinstance(other, LazyAlignments):
            other._load()
        return self.__le__(other)

    def __eq__(self, other):
        self._load()
        if isinstance(other, LazyAlignments):
            other._load()
        return self.__eq__(other)

    def __gt__(self, other):
        self._load()
        if isinstance(other, LazyAlignments):
            other._load()
        return self.__gt__(other)

    def __ge__(self, other):
        self._load()
        if isinstance(other, LazyAlignments):
            other._load()
        return self.__ge__(other)

    def __contains__(self, alignment):
        self._load()
        return self.__contains__(alignment)

    def __len__(self):
        self._load()
        return self.__len__()

    def __getitem__(self, i):
        self._load()
        if isinstance(i, slice):
            alignments = Alignments()
            items = self.__getitem__(i)
            alignments.extend(items)
            return alignments
        else:
            return self.__getitem__(i)

    def __setitem__(self, i, item):
        self._load()
        self.__setitem__(i, item)

    def __delitem__(self, i):
        self._load()
        self.__delitem__(i)

    def __add__(self, other):
        alignments = Alignments()
        alignments.extend(self)
        alignments.extend(other)
        return alignments

    def __radd__(self, other):
        alignments = Alignments()
        alignments.extend(other)
        alignments.extend(self)
        return alignments

    def __iadd__(self, other):
        self._load()
        self.extend(other)
        return self

    def __mul__(self, n):
        alignments = Alignments()
        items = list(self)
        for i in range(n):
            alignments.extend(items)
        return alignments

    def __rmul__(self, n):
        alignments = Alignments()
        items = list(self)
        for i in range(n):
            alignments.extend(items)
        return alignments

    def __imul__(self, n):
        self._load()
        return self.__imul__(n)

    def append(self, item):
        self._load()
        self.append(item)

    def insert(self, i, item):
        self._load()
        self.insert(i, item)

    def pop(self, i=-1):
        self._load()
        return self.pop(i)

    def remove(self, item):
        self._load()
        self.remove(item)

    def copy(self):
        self._load()
        alignments = Alignments()
        alignments.__dict__.update(self.__dict__)
        alignments.extend(self)
        return alignments

    def count(self, item):
        self._load()
        return self.count(item)

    def index(self, item, *args):
        self._load()
        return self.index(item, *args)

    def reverse(self):
        self._load()
        self.reverse()

    def sort(self, /, *args, **kwds):
        self._load()
        self.sort(*args, **kwds)

    def extend(self, other):
        self._load()
        self.extend(other)


class ParsedAlignments(LazyAlignments, ABC):
    # The ParsedAlignments class is an abstract base class for parsing sequence
    # alignment files. The alignment parser modules in Bio.Align implement
    # concrete subclasses of ParsedAlignments.
    #
    # To write a new parser, you would create a new private module in Bio.Align,
    # and define an AlignmentIterator class in it as a subclass of
    # ParsedAlignments. Typically, you would override the _read_header and
    # _read_next_alignment methods in this subclass, as well as the __init__
    # method to call the base class __init__ method with the appropriate
    # arguments.

    def __init__(self, source, mode="t", fmt=None):
        """Create an ParsedAlignments object.

        Arguments:
        - source - input file stream, or path to input file

        This method MAY be overridden by any subclass.

        Note when subclassing:
        - there should be a single non-optional argument, the source.
        - you can add additional optional arguments.
        """
        self.source = source
        if source is None:
            return
        try:
            self._stream = open(source, "r" + mode)
        except TypeError:  # not a path, assume we received a stream
            if mode == "t":
                if source.read(0) != "":
                    raise StreamModeError(
                        "%s files must be opened in text mode." % fmt
                    ) from None
            elif mode == "b":
                if source.read(0) != b"":
                    raise StreamModeError(
                        "%s files must be opened in binary mode." % fmt
                    ) from None
            else:
                raise ValueError("Unknown mode '%s'" % mode) from None
            self._stream = source
        try:
            self._read_header(self._stream)
        except Exception:
            self._close()
            raise

    def __next__(self):
        try:
            stream = self._stream
        except AttributeError:
            raise StopIteration from None
        try:
            alignment = self._read_next_alignment(stream)
            if alignment is None:
                raise StopIteration
        except Exception:
            self._close()
            raise
        return alignment

    def __iter__(self):
        return self

    def _read_header(self, stream):
        """Read the file header and store it in metadata."""

    @abstractmethod
    def _read_next_alignment(self, stream):
        """Read one Alignment from the stream, and return it."""

    def _close(self):
        try:
            stream = self._stream
        except AttributeError:
            return
        if stream is not self.source:
            stream.close()
        del self._stream

    def clear(self):
        self._close()


class AlignmentWriter:
    """Base class for alignment writers. This class should be subclassed.

    It is intended for sequential file formats with an (optional)
    header, one or more alignments, and an (optional) footer.

    The user may call the write_file() method to write a complete
    file containing the alignments.

    Alternatively, users may call the write_header(), followed
    by multiple calls to format_alignment() and/or write_alignments(),
    followed finally by write_footer().

    Note that write_header() cannot require any assumptions about
    the number of alignments.
    """

    def __init__(self, target, mode="w"):
        """Create the writer object."""
        if target is not None:
            # target is None if we only use the writer to format strings.
            if mode == "w":
                try:
                    target.write("")
                except TypeError:
                    # target was opened in binary mode
                    raise StreamModeError("File must be opened in text mode.") from None
                except AttributeError:
                    # target is a path
                    stream = open(target, mode)
                else:
                    stream = target
            elif mode == "wb":
                try:
                    target.write(b"")
                except TypeError:
                    # target was opened in text mode
                    raise StreamModeError(
                        "File must be opened in binary mode."
                    ) from None
                except AttributeError:
                    # target is a path
                    stream = open(target, mode)
                else:
                    stream = target
            else:
                raise RuntimeError("Unknown mode '%s'" % mode)
            self.stream = stream

        self._target = target

    def write_header(self, alignments):
        """Write the file header to the output file."""
        pass
        ##################################################
        # You MUST implement this method in the subclass #
        # if the file format defines a file header.      #
        ##################################################

    def write_footer(self):
        """Write the file footer to the output file."""
        pass
        ##################################################
        # You MUST implement this method in the subclass #
        # if the file format defines a file footer.      #
        ##################################################

    def format_alignment(self, alignment):
        """Format a single alignment as a string.

        alignment - an Alignment object
        """
        raise NotImplementedError("This method should be implemented")
        ###################################################
        # You MUST implement this method in the subclass. #
        ###################################################

    def write_alignments(self, alignments, maxcount=None):
        """Write alignments to the output file, and return the number of alignments.

        alignments - A list or iterator returning Alignment objects
        maxcount - The maximum number of alignments allowed by the
        file format, or None if there is no maximum.
        """
        count = 0
        if maxcount is None:
            for alignment in alignments:
                line = self.format_alignment(alignment)
                self.stream.write(line)
                count += 1
        else:
            for alignment in alignments:
                if count == maxcount:
                    if maxcount == 1:
                        raise ValueError("More than one alignment found")
                    else:
                        raise ValueError(
                            "Number of alignments is larger than %d" % maxcount
                        )
                line = self.format_alignment(alignment)
                self.stream.write(line)
                count += 1
        return count

    def write_file(self, alignments, mincount=0, maxcount=None):
        """Write a file with the alignments, and return the number of alignments.

        alignments - A list or iterator returning Alignment objects
        """
        try:
            self.write_header(alignments)
            count = self.write_alignments(alignments, maxcount)
            self.write_footer()
        finally:
            if self.stream is not self._target:
                self.stream.close()
        if count < mincount:
            if mincount == 1:  # Common case
                raise ValueError("Must have one alignment")
            elif mincount == maxcount:
                raise ValueError(
                    "Number of alignments is %d (expected %d)" % (count, mincount)
                )
            else:
                raise ValueError(
                    "Number of alignments is %d (expected at least %d)"
                    % (count, mincount)
                )
        return count


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
