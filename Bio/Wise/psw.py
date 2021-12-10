#!/usr/bin/env python
# Copyright 2004 by Michael Hoffman. All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Run and process output from the Wise2 package tool psw.

Bio.Wise contains modules for running and processing the output of
some of the models in the Wise2 package by Ewan Birney available from:
ftp://ftp.ebi.ac.uk/pub/software/unix/wise2/
http://www.ebi.ac.uk/Wise2/

Bio.Wise.psw is for protein Smith-Waterman alignments
Bio.Wise.dnal is for Smith-Waterman DNA alignments
"""


import os
import re
import sys

from Bio import Wise


_CMDLINE_PSW = ["psw", "-l", "-F"]
_OPTION_GAP_START = "-g"
_OPTION_GAP_EXTENSION = "-e"
_OPTION_SCORES = "-m"


class AlignmentColumnFullException(Exception):
    """Manage exception in the alignment output."""

    pass


class Alignment(list):
    """Define a container for all alignment Columns, output from running psw."""

    def append(self, column_unit):
        """Add an alignment Column to Alignment."""
        try:
            self[-1].append(column_unit)
        except AlignmentColumnFullException:
            list.append(self, AlignmentColumn(column_unit))
        except IndexError:
            list.append(self, AlignmentColumn(column_unit))


class AlignmentColumn(list):
    """Define a container for the units that made the Column."""

    def _set_kind(self, column_unit):
        if self.kind == "SEQUENCE":
            self.kind = column_unit.kind

    def __init__(self, column_unit):
        """Initialize the class."""
        assert column_unit.unit == 0
        self.kind = column_unit.kind
        list.__init__(self, [column_unit.column, None])

    def __repr__(self):
        """Represent the AlignmentColumn object as a string for debugging."""
        return f"{self.kind}({self[0]!r}, {self[1]!r})"

    def append(self, column_unit):
        """Add a unit to the Column."""
        if self[1] is not None:
            raise AlignmentColumnFullException

        assert column_unit.unit == 1

        self._set_kind(column_unit)
        self[1] = column_unit.column


class ColumnUnit:
    """Define a container for the details of each sequence alignment."""

    def __init__(self, unit, column, kind):
        """Initialize the class."""
        self.unit = unit
        self.column = column
        self.kind = kind

    def __repr__(self):
        """Represent the ColumnUnit object as a string for debugging."""
        return "ColumnUnit(unit=%r, column=%r, kind=%r)" % (
            self.unit,
            self.column,
            self.kind,
        )


_re_unit = re.compile(r"^Unit +([01])- \[ *(-?\d+)- *(-?\d+)\] \[(\w+)\]$")


def parse_line(line):
    """Parse a line from psw.

    >>> print(parse_line("Column 0:"))
    None
    >>> parse_line("Unit  0- [  -1-   0] [SEQUENCE]")
    ColumnUnit(unit=0, column=0, kind='SEQUENCE')
    >>> parse_line("Unit  1- [  85-  86] [SEQUENCE]")
    ColumnUnit(unit=1, column=86, kind='SEQUENCE')
    """
    match = _re_unit.match(line.rstrip())

    if not match:
        return

    return ColumnUnit(int(match.group(1)), int(match.group(3)), match.group(4))


def parse(iterable):
    """Parse a file.

    format

    Column 0:
    Unit  0- [  -1-   0] [SEQUENCE]
    Unit  1- [  85-  86] [SEQUENCE]

    means that seq1[0] == seq2[86] (0-based)
    """
    alignment = Alignment()
    for line in iterable:
        try:
            if os.environ["WISE_PY_DEBUG"]:
                print(line)
        except KeyError:
            pass

        column_unit = parse_line(line)
        if column_unit:
            alignment.append(column_unit)

    return alignment


def align(pair, scores=None, gap_start=None, gap_extension=None, *args, **keywds):
    """Align a pair of DNA files using Wise2 psw."""
    cmdline = _CMDLINE_PSW[:]
    if scores:
        cmdline.extend((_OPTION_SCORES, scores))
    if gap_start:
        cmdline.extend((_OPTION_GAP_START, str(gap_start)))
    if gap_extension:
        cmdline.extend((_OPTION_GAP_EXTENSION, str(gap_extension)))
    temp_file = Wise.align(cmdline, pair, *args, **keywds)
    return parse(temp_file)


def main():
    """Command line implementation."""
    print(align(sys.argv[1:3]))


def _test(*args, **keywds):
    import doctest

    doctest.testmod(sys.modules[__name__], *args, **keywds)


if __name__ == "__main__":
    if __debug__:
        _test()
        """Initialize the class."""
    main()
