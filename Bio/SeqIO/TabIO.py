# Copyright 2008-2017,2020 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.SeqIO support for the "tab" (simple tab separated) file format.

You are expected to use this module via the Bio.SeqIO functions.

The "tab" format is an ad-hoc plain text file format where each sequence is
on one (long) line.  Each line contains the identifier/description, followed
by a tab, followed by the sequence.  For example, consider the following
short FASTA format file::

    >ID123456 possible binding site?
    CATCNAGATGACACTACGACTACGACTCAGACTAC
    >ID123457 random sequence
    ACACTACGACTACGACTCAGACTACAAN

Apart from the descriptions, this can be represented in the simple two column
tab separated format as follows::

    ID123456(tab)CATCNAGATGACACTACGACTACGACTCAGACTAC
    ID123457(tab)ACACTACGACTACGACTCAGACTACAAN

When reading this file, "ID123456" or "ID123457" will be taken as the record's
.id and .name property.  There is no other information to record.

Similarly, when writing to this format, Biopython will ONLY record the record's
.id and .seq (and not the description or any other information) as in the
example above.
"""
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .Interfaces import _clean
from .Interfaces import _get_seq_string
from .Interfaces import SequenceIterator
from .Interfaces import SequenceWriter


class TabIterator(SequenceIterator):
    """Parser for tab-delimited files."""

    def __init__(self, source):
        """Iterate over tab separated lines as SeqRecord objects.

        Each line of the file should contain one tab only, dividing the line
        into an identifier and the full sequence.

        Arguments:
         - source - file-like object opened in text mode, or a path to a file

        The first field is taken as the record's .id and .name (regardless of
        any spaces within the text) and the second field is the sequence.

        Any blank lines are ignored.

        Examples
        --------
        >>> with open("GenBank/NC_005816.tsv") as handle:
        ...     for record in TabIterator(handle):
        ...         print("%s length %i" % (record.id, len(record)))
        gi|45478712|ref|NP_995567.1| length 340
        gi|45478713|ref|NP_995568.1| length 260
        gi|45478714|ref|NP_995569.1| length 64
        gi|45478715|ref|NP_995570.1| length 123
        gi|45478716|ref|NP_995571.1| length 145
        gi|45478717|ref|NP_995572.1| length 357
        gi|45478718|ref|NP_995573.1| length 138
        gi|45478719|ref|NP_995574.1| length 312
        gi|45478720|ref|NP_995575.1| length 99
        gi|45478721|ref|NP_995576.1| length 90

        """
        super().__init__(source, mode="t", fmt="Tab-separated plain-text")

    def parse(self, handle):
        """Start parsing the file, and return a SeqRecord generator."""
        records = self.iterate(handle)
        return records

    def iterate(self, handle):
        """Parse the file and generate SeqRecord objects."""
        for line in handle:
            try:
                title, seq = line.split("\t")  # will fail if more than one tab!
            except ValueError:
                if line.strip() == "":
                    # It's a blank line, ignore it
                    continue
                raise ValueError(
                    "Each line should have one tab separating the"
                    + " title and sequence, this line has %i tabs: %r"
                    % (line.count("\t"), line)
                ) from None
            title = title.strip()
            seq = seq.strip()  # removes the trailing new line
            yield SeqRecord(Seq(seq), id=title, name=title, description="")


class TabWriter(SequenceWriter):
    """Class to write simple tab separated format files.

    Each line consists of "id(tab)sequence" only.

    Any description, name or other annotation is not recorded.

    This class is not intended to be used directly. Instead, please use
    the function ``as_tab``, or the top level ``Bio.SeqIO.write()`` function
    with ``format="tab"``.
    """

    def write_record(self, record):
        """Write a single tab line to the file."""
        assert self._header_written
        assert not self._footer_written
        self._record_written = True
        self.handle.write(as_tab(record))


def as_tab(record):
    """Return record as tab separated (id(tab)seq) string."""
    title = _clean(record.id)
    seq = _get_seq_string(record)  # Catches sequence being None
    assert "\t" not in title
    assert "\n" not in title
    assert "\r" not in title
    assert "\t" not in seq
    assert "\n" not in seq
    assert "\r" not in seq
    return f"{title}\t{seq}\n"


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest(verbose=0)
