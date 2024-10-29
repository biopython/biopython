# Copyright 2006-2017,2020 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
#
# This module is for reading and writing FASTA format files as SeqRecord
# objects.  The code is partly inspired  by earlier Biopython modules,
# Bio.Fasta.* and the now removed module Bio.SeqIO.FASTA
"""Bio.SeqIO support for the "fasta" (aka FastA or Pearson) file format.

You are expected to use this module via the Bio.SeqIO functions.
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import BiopythonDeprecationWarning


from .Interfaces import _clean
from .Interfaces import _get_seq_string
from .Interfaces import _TextIOSource
from .Interfaces import SequenceIterator
from .Interfaces import SequenceWriter

import warnings


def SimpleFastaParser(handle):
    """Iterate over Fasta records as string tuples.

    Arguments:
     - handle - input stream opened in text mode

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
    for line in handle:
        if line[0] == ">":
            title = line[1:].rstrip()
            break
    else:
        # no break encountered - probably an empty file
        return

    # Main logic
    # Note, remove trailing whitespace, and any internal spaces
    # (and any embedded \r which are possible in mangled files
    # when not opened in universal read lines mode)
    lines = []
    for line in handle:
        if line[0] == ">":
            yield title, "".join(lines).replace(" ", "").replace("\r", "")
            lines = []
            title = line[1:].rstrip()
            continue
        lines.append(line.rstrip())

    yield title, "".join(lines).replace(" ", "").replace("\r", "")


def FastaTwoLineParser(handle):
    """Iterate over no-wrapping Fasta records as string tuples.

    Arguments:
     - handle - input stream opened in text mode

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
            if line[0] != ">":
                raise ValueError(
                    "Expected FASTA record starting with '>' character. "
                    "Perhaps this file is using FASTA line wrapping? "
                    f"Got: '{line}'"
                )
            title = line[1:].rstrip()
        else:  # sequence line
            if line[0] == ">":
                raise ValueError(
                    "Two '>' FASTA lines in a row. Missing sequence line "
                    "if this is strict two-line-per-record FASTA format. "
                    f"Have '>{title}' and '{line}'"
                )
            yield title, line.strip()

    if idx == -1:
        pass  # empty file
    elif idx % 2 == 0:  # on a title line
        raise ValueError(
            "Missing sequence line at end of file if this is strict "
            f"two-line-per-record FASTA format. Have title line '{line}'"
        )
    else:
        assert line[0] != ">", "line[0] == '>' ; this should be impossible!"


class FastaIterator(SequenceIterator):
    """Parser for plain Fasta files without comments."""

    modes = "t"

    def __init__(
        self,
        source: _TextIOSource,
        alphabet: None = None,
    ) -> None:
        """Iterate over Fasta records as SeqRecord objects.

        Arguments:
         - source - input stream opened in text mode, or a path to a file
         - alphabet - optional alphabet, not used. Leave as None.

        This parser expects a plain Fasta format without comments or header
        lines.

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

        If you want to modify the records before writing, for example to change
        the ID of each record, you can use a generator function as follows:

        >>> def modify_records(records):
        ...     for record in records:
        ...         record.id = record.id.upper()
        ...         yield record
        ...
        >>> with open('Fasta/dups.fasta') as handle:
        ...     for record in modify_records(FastaIterator(handle)):
        ...         print(record.id)
        ...
        ALPHA
        BETA
        GAMMA
        ALPHA
        DELTA

        """
        if alphabet is not None:
            raise ValueError("The alphabet argument is no longer supported")
        super().__init__(source, fmt="Fasta")
        try:
            line = next(self.stream)
        except StopIteration:
            line = None
        else:
            if not line.startswith(">"):
                warnings.warn(
                    "Previously, the FASTA parser silently ignored comments at the "
                    "beginning of the FASTA file (before the first sequence).\n"
                    "\n"
                    "Nowadays, the FASTA file format is usually understood not to "
                    "have any such comments, and most software packages do not allow "
                    "them. Therefore, the use of comments at the beginning of a FASTA "
                    "file is now deprecated in Biopython.\n"
                    "\n"
                    "In a future Biopython release, this deprecation warning will be "
                    "replaced by a ValueError. To avoid this, there are three "
                    "options:\n"
                    "\n"
                    "(1) Modify your FASTA file to remove such comments at the "
                    "beginning of the file.\n"
                    "\n"
                    "(2) Use SeqIO.parse with the 'fasta-pearson' format instead of "
                    "'fasta'. This format is consistent with the FASTA format defined "
                    "by William Pearson's FASTA aligner software. Thie format allows "
                    "for comments before the first sequence; lines starting with the "
                    "';' character anywhere in the file are also regarded as comment "
                    "lines and are ignored.\n"
                    "\n"
                    "(3) Use the 'fasta-blast' format. This format regards any lines "
                    "starting with '!', '#', or ';' as comment lines. The "
                    "'fasta-blast' format may be safer than the 'fasta-pearson' "
                    "format, as it explicitly indicates which lines are comments. ",
                    BiopythonDeprecationWarning,
                )
                for line in self.stream:
                    if line.startswith(">"):
                        break
                else:
                    line = None
        self._line = line

    def __next__(self):
        line = self._line
        if line is None:
            raise StopIteration
        title = line[1:].rstrip()
        # Main logic
        # Note, remove trailing whitespace, and any internal spaces
        # (and any embedded \r which are possible in mangled files
        # when not opened in universal read lines mode)
        lines = []
        for line in self.stream:
            if line[0] == ">":
                break
            lines.append(line)
        else:
            line = None
        self._line = line
        sequence = "".join(lines).encode().translate(None, b" \t\r\n")
        try:
            first_word = title.split(None, 1)[0]
        except IndexError:
            assert not title, repr(title)
            # Should we use SeqRecord default for no ID?
            first_word = ""
        return SeqRecord._from_validated(
            Seq(sequence), id=first_word, name=first_word, description=title
        )


class FastaTwoLineIterator(SequenceIterator):
    """Parser for Fasta files with exactly two lines per record."""

    modes = "t"

    def __init__(self, source):
        """Iterate over two-line Fasta records (as SeqRecord objects).

        Arguments:
         - source - input stream opened in text mode, or a path to a file

        This uses a strict interpretation of the FASTA as requiring
        exactly two lines per record (no line wrapping).

        Only the default title to ID/name/description parsing offered
        by the relaxed FASTA parser is offered.
        """
        super().__init__(source, fmt="FASTA")
        self._data = FastaTwoLineParser(self.stream)

    def __next__(self):
        try:
            title, sequence = next(self._data)
        except StopIteration:
            raise StopIteration from None
        try:
            first_word = title.split(None, 1)[0]
        except IndexError:
            assert not title, repr(title)
            # Should we use SeqRecord default for no ID?
            first_word = ""
        return SeqRecord(
            Seq(sequence), id=first_word, name=first_word, description=title
        )


class FastaBlastIterator(SequenceIterator):
    """Parser for Fasta files, allowing for comments as in BLAST."""

    modes = "t"

    def __init__(
        self,
        source: _TextIOSource,
        alphabet: None = None,
    ) -> None:
        """Iterate over Fasta records as SeqRecord objects.

        Arguments:
         - source - input stream opened in text mode, or a path to a file
         - alphabet - optional alphabet, not used. Leave as None.

        This parser expects the data to be in FASTA format. As in BLAST, lines
        starting with '#', '!', or ';' are interpreted as comments and ignored.

        This iterator acts like calling Bio.SeqIO.parse(handle, "fasta-blast")
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

        If you want to modify the records before writing, for example to change
        the ID of each record, you can use a generator function as follows:

        >>> def modify_records(records):
        ...     for record in records:
        ...         record.id = record.id.upper()
        ...         yield record
        ...
        >>> with open('Fasta/dups.fasta') as handle:
        ...     for record in modify_records(FastaIterator(handle)):
        ...         print(record.id)
        ...
        ALPHA
        BETA
        GAMMA
        ALPHA
        DELTA

        """
        if alphabet is not None:
            raise ValueError("The alphabet argument is no longer supported")
        super().__init__(source, fmt="FASTA")
        for line in self.stream:
            if line[0] not in "#!;":
                if not line.startswith(">"):
                    raise ValueError(
                        "Expected FASTA record starting with '>' character.\n"
                        "If this line is a comment, please use '#', '!', or ';' as "
                        "the first character, or use the 'fasta-pearson' "
                        "format for parsing.\n"
                        f"Got: '{line}'"
                    )
                self._line = line
                break
        else:
            self._line = None

    def __next__(self):
        line = self._line
        if line is None:
            raise StopIteration
        title = line[1:].rstrip()
        lines = []
        for line in self.stream:
            # Main logic
            # Note, remove trailing whitespace, and any internal spaces
            # (and any embedded \r which are possible in mangled files
            # when not opened in universal read lines mode)
            if line[0] in "#!;":
                pass
            elif line[0] == ">":
                self_line = line
                break
            else:
                lines.append(line.rstrip())
        else:
            self._line = None
        try:
            first_word = title.split(None, 1)[0]
        except IndexError:
            first_word = ""
        sequence = "".join(lines).replace(" ", "").replace("\r", "")
        return SeqRecord(
            Seq(sequence), id=first_word, name=first_word, description=title
        )


class FastaPearsonIterator(SequenceIterator):
    """Parser for Fasta files, allowing for comments as in the FASTA aligner."""

    modes = "t"

    def __init__(
        self,
        source: _TextIOSource,
        alphabet: None = None,
    ) -> None:
        """Iterate over Fasta records as SeqRecord objects.

        Arguments:
         - source - input stream opened in text mode, or a path to a file
         - alphabet - optional alphabet, not used. Leave as None.

        This parser expects a Fasta format allowing for a header (before the
        first sequence record) and comments (lines starting with ';') as in
        William Pearson's FASTA aligner software.

        This iterator acts as calling Bio.SeqIO.parse(handle, "fasta-pearson")
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

        If you want to modify the records before writing, for example to change
        the ID of each record, you can use a generator function as follows:

        >>> def modify_records(records):
        ...     for record in records:
        ...         record.id = record.id.upper()
        ...         yield record
        ...
        >>> with open('Fasta/dups.fasta') as handle:
        ...     for record in modify_records(FastaIterator(handle)):
        ...         print(record.id)
        ...
        ALPHA
        BETA
        GAMMA
        ALPHA
        DELTA

        """
        if alphabet is not None:
            raise ValueError("The alphabet argument is no longer supported")
        super().__init__(source, fmt="Fasta")
        for line in self.stream:
            if line.startswith(">"):
                self._line = line
                break
        else:
            self._line = None

    def __next__(self):
        line = self._line
        if line is None:
            raise StopIteration
        title = line[1:].rstrip()
        lines = []
        for line in self.stream:
            # Main logic
            # Note, remove trailing whitespace, and any internal spaces
            # (and any embedded \r which are possible in mangled files
            # when not opened in universal read lines mode)
            if line[0] == ";":
                pass
            elif line[0] == ">":
                self._line = line
                break
            else:
                lines.append(line.rstrip())
        else:
            self._line = None
        try:
            first_word = title.split(None, 1)[0]
        except IndexError:
            first_word = ""
        sequence = "".join(lines).replace(" ", "").replace("\r", "")
        return SeqRecord(
            Seq(sequence), id=first_word, name=first_word, description=title
        )


class FastaWriter(SequenceWriter):
    """Class to write Fasta format files (OBSOLETE).

    Please use the ``as_fasta`` function instead, or the top level
    ``Bio.SeqIO.write()`` function instead using ``format="fasta"``.
    """

    modes = "t"

    def __init__(self, target, wrap=60, record2title=None):
        """Create a Fasta writer (OBSOLETE).

        Arguments:
         - target - Output stream opened in text mode, or a path to a file.
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
            ...
            Multiple writer.write_record() and/or writer.write_records() calls
            ...
            handle.close()

        """
        super().__init__(target)
        if wrap:
            if wrap < 1:
                raise ValueError
        self.wrap = wrap
        self.record2title = record2title

    def write_record(self, record):
        """Write a single Fasta record to the file."""
        if self.record2title:
            title = self.clean(self.record2title(record))
        else:
            id = self.clean(record.id)
            description = self.clean(record.description)
            if description and description.split(None, 1)[0] == id:
                # The description includes the id at the start
                title = description
            elif description:
                title = f"{id} {description}"
            else:
                title = id

        assert "\n" not in title
        assert "\r" not in title
        self.handle.write(f">{title}\n")

        data = _get_seq_string(record)  # Catches sequence being None

        assert "\n" not in data
        assert "\r" not in data

        if self.wrap:
            for i in range(0, len(data), self.wrap):
                self.handle.write(data[i : i + self.wrap] + "\n")
        else:
            self.handle.write(data + "\n")


class FastaTwoLineWriter(FastaWriter):
    """Class to write 2-line per record Fasta format files (OBSOLETE).

    This means we write the sequence information  without line
    wrapping, and will always write a blank line for an empty
    sequence.

    Please use the ``as_fasta_2line`` function instead, or the top level
    ``Bio.SeqIO.write()`` function instead using ``format="fasta"``.
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
            ...
            Multiple writer.write_record() and/or writer.write_records() calls
            ...
            handle.close()

        """
        super().__init__(handle, wrap=None, record2title=record2title)


def as_fasta(record):
    """Turn a SeqRecord into a FASTA formatted string.

    This is used internally by the SeqRecord's .format("fasta")
    method and by the SeqIO.write(..., ..., "fasta") function.
    """
    id = _clean(record.id)
    description = _clean(record.description)
    if description and description.split(None, 1)[0] == id:
        # The description includes the id at the start
        title = description
    elif description:
        title = f"{id} {description}"
    else:
        title = id
    assert "\n" not in title
    assert "\r" not in title
    lines = [f">{title}\n"]

    data = _get_seq_string(record)  # Catches sequence being None
    assert "\n" not in data
    assert "\r" not in data
    for i in range(0, len(data), 60):
        lines.append(data[i : i + 60] + "\n")

    return "".join(lines)


def as_fasta_2line(record):
    """Turn a SeqRecord into a two-line FASTA formatted string.

    This is used internally by the SeqRecord's .format("fasta-2line")
    method and by the SeqIO.write(..., ..., "fasta-2line") function.
    """
    id = _clean(record.id)
    description = _clean(record.description)
    if description and description.split(None, 1)[0] == id:
        # The description includes the id at the start
        title = description
    elif description:
        title = f"{id} {description}"
    else:
        title = id
    assert "\n" not in title
    assert "\r" not in title

    data = _get_seq_string(record)  # Catches sequence being None
    assert "\n" not in data
    assert "\r" not in data

    return f">{title}\n{data}\n"


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest(verbose=0)
