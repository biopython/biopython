#!/usr/bin/env python
# Copyright 2004-2005 by Michael Hoffman. All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Run and process output from the Wise2 package tools.

Bio.Wise contains modules for running and processing the output of
some of the models in the Wise2 package by Ewan Birney available from:
ftp://ftp.ebi.ac.uk/pub/software/unix/wise2/
http://www.ebi.ac.uk/Wise2/

Bio.Wise.psw is for protein Smith-Waterman alignments
Bio.Wise.dnal is for Smith-Waterman DNA alignments
"""


import os
import sys
import tempfile

from Bio import SeqIO


def _build_align_cmdline(
    cmdline, pair, output_filename, kbyte=None, force_type=None, quiet=False
):
    """Build a command line string (PRIVATE).

    >>> os.environ["WISE_KBYTE"]="300000"
    >>> if os.isatty(sys.stderr.fileno()):
    ...    c = _build_align_cmdline(["dnal"], ("seq1.fna", "seq2.fna"),
    ...                             "/tmp/output", kbyte=100000)
    ...    assert c == 'dnal -kbyte 100000 seq1.fna seq2.fna > /tmp/output', c
    ...    c = _build_align_cmdline(["psw"], ("seq1.faa", "seq2.faa"),
    ...                             "/tmp/output_aa")
    ...    assert c == 'psw -kbyte 300000 seq1.faa seq2.faa > /tmp/output_aa', c
    ... else:
    ...    c = _build_align_cmdline(["dnal"], ("seq1.fna", "seq2.fna"),
    ...                             "/tmp/output", kbyte=100000)
    ...    assert c == 'dnal -kbyte 100000 -quiet seq1.fna seq2.fna > /tmp/output', c
    ...    c = _build_align_cmdline(["psw"], ("seq1.faa", "seq2.faa"),
    ...                             "/tmp/output_aa")
    ...    assert c == 'psw -kbyte 300000 -quiet seq1.faa seq2.faa > /tmp/output_aa', c

    """
    cmdline = cmdline[:]

    # XXX: force_type ignored

    if kbyte is None:
        try:
            cmdline.extend(("-kbyte", os.environ["WISE_KBYTE"]))
        except KeyError:
            pass
    else:
        cmdline.extend(("-kbyte", str(kbyte)))

    if not os.isatty(sys.stderr.fileno()):
        cmdline.append("-quiet")

    cmdline.extend(pair)
    cmdline.extend((">", output_filename))
    if quiet:
        cmdline.extend(("2>", "/dev/null"))
    return " ".join(cmdline)


def align(
    cmdline, pair, kbyte=None, force_type=None, dry_run=False, quiet=False, debug=False
):
    """Run an alignment. Returns a filehandle."""
    if not pair or len(pair) != 2:
        raise ValueError("Expected pair of filename, not %r" % pair)

    output_file = tempfile.NamedTemporaryFile(mode="r")
    input_files = (
        tempfile.NamedTemporaryFile(mode="w"),
        tempfile.NamedTemporaryFile(mode="w"),
    )

    if dry_run:
        print(
            _build_align_cmdline(
                cmdline, pair, output_file.name, kbyte, force_type, quiet
            )
        )
        return

    for filename, input_file in zip(pair, input_files):
        # Pipe the file through Biopython's Fasta parser/writer
        # to make sure it conforms to the Fasta standard (in particular,
        # Wise2 may choke on long lines in the Fasta file)
        records = SeqIO.parse(open(filename), "fasta")
        SeqIO.write(records, input_file, "fasta")
        input_file.flush()

    input_file_names = [input_file.name for input_file in input_files]

    cmdline_str = _build_align_cmdline(
        cmdline, input_file_names, output_file.name, kbyte, force_type, quiet
    )

    if debug:
        sys.stderr.write("%s\n" % cmdline_str)

    status = os.system(cmdline_str) >> 8

    # `status` here will be >1 for error codes >=256
    if status > 1:
        if kbyte != 0:  # possible memory problem; could be None
            sys.stderr.write("INFO trying again with the linear model\n")
            return align(cmdline, pair, 0, force_type, dry_run, quiet, debug)
        else:
            raise OSError("%s returned %s" % (" ".join(cmdline), status))

    return output_file


def all_pairs(singles):
    """Generate pairs list for all-against-all alignments.

    >>> all_pairs(range(4))
    [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    """
    pairs = []

    singles = list(singles)
    while singles:
        suitor = singles.pop(0)  # if sorted, stay sorted
        pairs.extend((suitor, single) for single in singles)

    return pairs


def main():
    """Provision for command line testing."""
    pass


def _test(*args, **keywds):
    import doctest

    doctest.testmod(sys.modules[__name__], *args, **keywds)


if __name__ == "__main__":
    if __debug__:
        _test()
    main()
