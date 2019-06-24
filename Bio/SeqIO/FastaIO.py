# Copyright 2006-2017 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
# This module is for reading and writing FASTA format files as SeqRecord
# objects.  The code is partly inspired  by earlier Biopython modules,
# Bio.Fasta.* and the now deprecated Bio.SeqIO.FASTA

"""Bio.SeqIO support for the "fasta" (aka FastA or Pearson) file format.

You are expected to use this module via the Bio.SeqIO functions.
"""

from __future__ import print_function

from Bio.Alphabet import single_letter_alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.Interfaces import SequentialSequenceWriter
from Bio.SeqIO.Interfaces import _clean, _get_seq_string


def SimpleFastaParser(handle):
    """Iterate over Fasta records as string tuples.

    For each record a tuple of two strings is returned, the FASTA title
    line (without the leading '>' character), and the sequence (with any
    whitespace removed). The title line is not divided up into an
    identifier (the first word) and comment or description.

    >>> with open("Fasta/dups.fasta") as handle:
    ...     for values in SimpleFastaParser(handle):
    ...         print(values)
    ...
    ('alpha', 'ACGTA')
    ('beta', 'CGTC')
    ('gamma', 'CCGCC')
    ('alpha (again - this is a duplicate entry to test the indexing code)', 'ACGTA')
    ('delta', 'CGCGC')

    """
    # Skip any text before the first record (e.g. blank lines, comments)
    # This matches the previous implementation where .readline() was used
    for line in handle:
        if line[0] == '>':
            title = line[1:].rstrip()
            break
        elif isinstance(line[0], int):
            # Same exception as for FASTQ files
            raise ValueError("Is this handle in binary mode not text mode?")
    else:
        # no break encountered - probably an empty file
        return

    # Main logic
    # Note, remove trailing whitespace, and any internal spaces
    # (and any embedded \r which are possible in mangled files
    # when not opened in universal read lines mode)
    lines = []
    for line in handle:
        if line[0] == '>':
            yield title, ''.join(lines).replace(" ", "").replace("\r", "")
            lines = []
            title = line[1:].rstrip()
            continue
        lines.append(line.rstrip())

    yield title, ''.join(lines).replace(" ", "").replace("\r", "")


def FastaTwoLineParser(handle):
    """Iterate over no-wrapping Fasta records as string tuples.

    Functionally the same as SimpleFastaParser but with a strict
    interpretation of the FASTA format as exactly two lines per
    record, the greater-than-sign identifier with description,
    and the sequence with no line wrapping.

    Any line wrapping will raise an exception, as will excess blank
    lines (other than the special case of a zero-length sequence
    as the second line of a record).

    Examples
    --------
    This file uses two lines per FASTA record:

    >>> with open("Fasta/aster_no_wrap.pro") as handle:
    ...     for title, seq in FastaTwoLineParser(handle):
    ...         print("%s = %s..." % (title, seq[:3]))
    ...
    gi|3298468|dbj|BAA31520.1| SAMIPF = GGH...

    This equivalent file uses line wrapping:

    >>> with open("Fasta/aster.pro") as handle:
    ...     for title, seq in FastaTwoLineParser(handle):
    ...         print("%s = %s..." % (title, seq[:3]))
    ...
    Traceback (most recent call last):
       ...
    ValueError: Expected FASTA record starting with '>' character. Perhaps this file is using FASTA line wrapping? Got: 'MTFGLVYTVYATAIDPKKGSLGTIAPIAIGFIVGANI'

    """
    idx = -1  # for empty file
    for idx, line in enumerate(handle):
        if idx % 2 == 0:  # title line
            if line[0] != '>':
                raise ValueError("Expected FASTA record starting with '>' character. "
                                 "Perhaps this file is using FASTA line wrapping? "
                                 "Got: '{}'".format(line))
            title = line[1:].rstrip()
        else:  # sequence line
            if line[0] == '>':
                raise ValueError("Two '>' FASTA lines in a row. Missing sequence line "
                                 "if this is strict two-line-per-record FASTA format. "
                                 "Have '>{}' and '{}'".format(title, line))
            yield title, line.strip()

    if idx == -1:
        pass  # empty file
    elif idx % 2 == 0:  # on a title line
        raise ValueError("Missing sequence line at end of file "
                         "if this is strict two-line-per-record FASTA format. "
                         "Have title line '{}'".format(line))
    else:
        assert line[0] != '>', "line[0] == '>' ; this should be impossible!"


def FastaIterator(handle, alphabet=single_letter_alphabet, title2ids=None):
    """Iterate over Fasta records as SeqRecord objects.

    Arguments:
     - handle - input file
     - alphabet - optional alphabet
     - title2ids - A function that, when given the title of the FASTA
       file (without the beginning >), will return the id, name and
       description (in that order) for the record as a tuple of strings.
       If this is not given, then the entire title line will be used
       as the description, and the first word as the id and name.

    By default this will act like calling Bio.SeqIO.parse(handle, "fasta")
    with no custom handling of the title lines:

    >>> with open("Fasta/dups.fasta") as handle:
    ...     for record in FastaIterator(handle):
    ...         print(record.id)
    ...
    alpha
    beta
    gamma
    alpha
    delta

    However, you can supply a title2ids function to alter this:

    >>> def take_upper(title):
    ...     return title.split(None, 1)[0].upper(), "", title
    >>> with open("Fasta/dups.fasta") as handle:
    ...     for record in FastaIterator(handle, title2ids=take_upper):
    ...         print(record.id)
    ...
    ALPHA
    BETA
    GAMMA
    ALPHA
    DELTA

    """
    if title2ids:
        for title, sequence in SimpleFastaParser(handle):
            id, name, descr = title2ids(title)
            yield SeqRecord(Seq(sequence, alphabet),
                            id=id, name=name, description=descr)
    else:
        for title, sequence in SimpleFastaParser(handle):
            try:
                first_word = title.split(None, 1)[0]
            except IndexError:
                assert not title, repr(title)
                # Should we use SeqRecord default for no ID?
                first_word = ""
            yield SeqRecord(Seq(sequence, alphabet),
                            id=first_word, name=first_word, description=title)


def FastaTwoLineIterator(handle, alphabet=single_letter_alphabet):
    """Iterate over two-line Fasta records (as SeqRecord objects).

    Arguments:
     - handle - input file
     - alphabet - optional alphabet

    This uses a strict interpretation of the FASTA as requiring
    exactly two lines per record (no line wrapping).

    Only the default title to ID/name/description parsing offered
    by the relaxed FASTA parser is offered.
    """
    for title, sequence in FastaTwoLineParser(handle):
        try:
            first_word = title.split(None, 1)[0]
        except IndexError:
            assert not title, repr(title)
            # Should we use SeqRecord default for no ID?
            first_word = ""
        yield SeqRecord(Seq(sequence, alphabet),
                        id=first_word, name=first_word, description=title)


class FastaWriter(SequentialSequenceWriter):
    """Class to write Fasta format files (OBSOLETE).

    Please use the ``as_fasta`` function instead, or the top level
    ``Bio.SeqIO.write()`` function instead using ``format="fasta"``.
    """

    def __init__(self, handle, wrap=60, record2title=None):
        """Create a Fasta writer (OBSOLETE).

        Arguments:
         - handle - Handle to an output file, e.g. as returned
           by open(filename, "w")
         - wrap -   Optional line length used to wrap sequence lines.
           Defaults to wrapping the sequence at 60 characters
           Use zero (or None) for no wrapping, giving a single
           long line for the sequence.
         - record2title - Optional function to return the text to be
           used for the title line of each record.  By default
           a combination of the record.id and record.description
           is used.  If the record.description starts with the
           record.id, then just the record.description is used.

        You can either use::

            handle = open(filename, "w")
            writer = FastaWriter(handle)
            writer.write_file(myRecords)
            handle.close()

        Or, follow the sequential file writer system, for example::

            handle = open(filename, "w")
            writer = FastaWriter(handle)
            writer.write_header() # does nothing for Fasta files
            ...
            Multiple writer.write_record() and/or writer.write_records() calls
            ...
            writer.write_footer() # does nothing for Fasta files
            handle.close()

        """
        SequentialSequenceWriter.__init__(self, handle)
        self.wrap = None
        if wrap:
            if wrap < 1:
                raise ValueError
        self.wrap = wrap
        self.record2title = record2title

    def write_record(self, record):
        """Write a single Fasta record to the file."""
        assert self._header_written
        assert not self._footer_written
        self._record_written = True

        if self.record2title:
            title = self.clean(self.record2title(record))
        else:
            id = self.clean(record.id)
            description = self.clean(record.description)
            if description and description.split(None, 1)[0] == id:
                # The description includes the id at the start
                title = description
            elif description:
                title = "%s %s" % (id, description)
            else:
                title = id

        assert "\n" not in title
        assert "\r" not in title
        self.handle.write(">%s\n" % title)

        data = self._get_seq_string(record)  # Catches sequence being None

        assert "\n" not in data
        assert "\r" not in data

        if self.wrap:
            for i in range(0, len(data), self.wrap):
                self.handle.write(data[i:i + self.wrap] + "\n")
        else:
            self.handle.write(data + "\n")


class FastaTwoLineWriter(FastaWriter):
    """Class to write 2-line per record Fasta format files (OBSOLETE).

    This means we write the sequence information  without line
    wrapping, and will always write a blank line for an empty
    sequence.
    """

    def __init__(self, handle, record2title=None):
        """Create a 2-line per record Fasta writer (OBSOLETE).

        Arguments:
         - handle - Handle to an output file, e.g. as returned
           by open(filename, "w")
         - record2title - Optional function to return the text to be
           used for the title line of each record.  By default
           a combination of the record.id and record.description
           is used.  If the record.description starts with the
           record.id, then just the record.description is used.

        You can either use::

            handle = open(filename, "w")
            writer = FastaWriter(handle)
            writer.write_file(myRecords)
            handle.close()

        Or, follow the sequential file writer system, for example::

            handle = open(filename, "w")
            writer = FastaWriter(handle)
            writer.write_header() # does nothing for Fasta files
            ...
            Multiple writer.write_record() and/or writer.write_records() calls
            ...
            writer.write_footer() # does nothing for Fasta files
            handle.close()

        """
        super(FastaTwoLineWriter, self).__init__(handle, wrap=None,
                                                 record2title=record2title)


def as_fasta(record):
    """Turn a SeqRecord into a FASTA formated string.

    This is used internally by the SeqRecord's .format("fasta")
    method and by the SeqIO.write(..., ..., "fasta") function.
    """
    id = _clean(record.id)
    description = _clean(record.description)
    if description and description.split(None, 1)[0] == id:
        # The description includes the id at the start
        title = description
    elif description:
        title = "%s %s" % (id, description)
    else:
        title = id
    assert "\n" not in title
    assert "\r" not in title
    lines = [">%s\n" % title]

    data = _get_seq_string(record)  # Catches sequence being None
    assert "\n" not in data
    assert "\r" not in data
    for i in range(0, len(data), 60):
        lines.append(data[i:i + 60] + "\n")

    return "".join(lines)


def as_fasta_2line(record):
    """Turn a SeqRecord into a two-line FASTA formated string.

    This is used internally by the SeqRecord's .format("fasta-2line")
    method and by the SeqIO.write(..., ..., "fasta-2line") function.
    """
    id = _clean(record.id)
    description = _clean(record.description)
    if description and description.split(None, 1)[0] == id:
        # The description includes the id at the start
        title = description
    elif description:
        title = "%s %s" % (id, description)
    else:
        title = id
    assert "\n" not in title
    assert "\r" not in title

    data = _get_seq_string(record)  # Catches sequence being None
    assert "\n" not in data
    assert "\r" not in data

    return ">%s\n%s\n" % (title, data)


if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest(verbose=0)
