# Copyright 2008-2015 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
# This module is for reading and writing IntelliGenetics format files as
# SeqRecord objects.  This file format appears to be the same as the MASE
# multiple sequence alignment format.

"""Bio.SeqIO support for the "ig" (IntelliGenetics or MASE) file format.

You are expected to use this module via the Bio.SeqIO functions.
"""

from __future__ import print_function

from Bio.Alphabet import single_letter_alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def IgIterator(handle, alphabet=single_letter_alphabet):
    """Iterate over IntelliGenetics records (as SeqRecord objects).

    handle - input file
    alphabet - optional alphabet

    The optional free format file header lines (which start with two
    semi-colons) are ignored.

    The free format commentary lines at the start of each record (which
    start with a semi-colon) are recorded as a single string with embedded
    new line characters in the SeqRecord's annotations dictionary under the
    key 'comment'.

    Examples
    --------
    >>> with open("IntelliGenetics/TAT_mase_nuc.txt") as handle:
    ...     for record in IgIterator(handle):
    ...         print("%s length %i" % (record.id, len(record)))
    ...
    A_U455 length 303
    B_HXB2R length 306
    C_UG268A length 267
    D_ELI length 309
    F_BZ163A length 309
    O_ANT70 length 342
    O_MVP5180 length 348
    CPZGAB length 309
    CPZANT length 309
    A_ROD length 390
    B_EHOA length 420
    D_MM251 length 390
    STM_STM length 387
    VER_AGM3 length 354
    GRI_AGM677 length 264
    SAB_SAB1C length 219
    SYK_SYK length 330

    """
    # Skip any file header text before the first record (;; lines)
    while True:
        line = handle.readline()
        if not line:
            break  # Premature end of file, or just empty?
        if not line.startswith(";;"):
            break

    while line:
        # Now iterate over the records
        if line[0] != ";":
            raise ValueError(
                "Records should start with ';' and not:\n%r" % line)

        # Try and agree with SeqRecord convention from the GenBank parser,
        # (and followed in the SwissProt parser) which stores the comments
        # as a long string with newlines under annotations key 'comment'.

        # Note some examples use "; ..." and others ";..."
        comment_lines = []
        while line.startswith(";"):
            # TODO - Extract identifier from lines like "LOCUS\tB_SF2"?
            comment_lines.append(line[1:].strip())
            line = handle.readline()
        title = line.rstrip()

        seq_lines = []
        while True:
            line = handle.readline()
            if not line:
                break
            if line[0] == ";":
                break
            # Remove trailing whitespace, and any internal spaces
            seq_lines.append(line.rstrip().replace(" ", ""))
        seq_str = "".join(seq_lines)
        if seq_str.endswith("1"):
            # Remove the optional terminator (digit one)
            seq_str = seq_str[:-1]
        if "1" in seq_str:
            raise ValueError(
                "Potential terminator digit one found within sequence.")

        # Return the record and then continue...
        record = SeqRecord(Seq(seq_str, alphabet),
                           id=title, name=title)
        record.annotations['comment'] = "\n".join(comment_lines)
        yield record

    # We should be at the end of the file now
    assert not line


if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest(verbose=0)
