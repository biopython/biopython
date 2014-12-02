# Copyright 2008-2010 by Peter Cock.  All rights reserved.
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

__docformat__ = "restructuredtext en"


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
                "Records should start with ';' and not:\n%s" % repr(line))

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
    print("Running quick self test")

    import os
    path = "../../Tests/IntelliGenetics/"
    if os.path.isdir(path):
        for filename in os.listdir(path):
            if os.path.splitext(filename)[-1] == ".txt":
                print("")
                print(filename)
                print("-" * len(filename))
                with open(os.path.join(path, filename)) as handle:
                    for record in IgIterator(handle):
                        print("%s %i" % (record.id, len(record)))
        print("Done")
    else:
        print("Could not find input files")
