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


# NCBI FASTA header identifier prefixes.
# Based on https://ncbi.github.io/cxx-toolkit/pages/ch_demo#ch_demo.id1_fetch.html_ref_fasta
# Cross-referenced with NCBI's dbxref list:
#   https://www.ncbi.nlm.nih.gov/genbank/collab/db_xref/
# and UniProt's dbxref list:
#   https://www.uniprot.org/docs/dbxref
ncbi_fasta_header_fields = {
    "bbs": (("id",), "GenInfo backbone seqid"),
    "bbm": (("id",), "GenInfo import ID"),
    "gb": (("accession", "locus"), "GenBank"),
    "emb": (("accession", "locus"), "EMBL"),
    "pir": (("accession", "name"), "PIR"),
    "sp": (("accession", "name"), "UniProtKB/Swiss-Prot"),
    "pat": (("country", "patent", "sequence"), "Patent"),
    "pgp": (("country", "application_number", "sequence"), "Pre-grant patent"),
    "ref": (("accession", "name"), "RefSeq"),
    "gnl": (("database", "id"), "General database reference"),
    "gi": (("id",), "GI"),
    "dbj": (("accession", "locus"), "DDBJ"),
    "prf": (("accession", "name"), "PRF"),
    "pdb": (("id", "chain"), "PDB"),
    "tpg": (("accession", "name"), "Third-party GenBank"),
    "tpe": (("accession", "name"), "Third-party EMBL"),
    "tpd": (("accession", "name"), "Third-party DDBJ"),
    "tr": (("accession", "name"), "UniProtKB/TrEMBL"),
    "gpp": (("accession", "name"), "Genome pipeline"),
    "nat": (("accession", "name"), "Named annotation track"),
}


def _parse_ncbi_fasta_header(title):
    """Parse NCBI-style pipe-delimited FASTA title lines into dbxrefs.

    Parses title lines formatted according to the NCBI FASTA standard:
    https://ncbi.github.io/cxx-toolkit/pages/ch_demo#ch_demo.id1_fetch.html_ref_fasta

    Returns a tuple of (id, name, dbxrefs) where:
     - id is the first whitespace-delimited word (the raw pipe string)
     - name is the first recognized accession, or same as id if none found
     - dbxrefs is a list of "DB:accession" strings

    Arguments:
     - title - title line as a stripped string without '>' as parsed by
       SimpleFastaParser

    >>> _parse_ncbi_fasta_header("gi|3298468|dbj|BAA31520.1| SAMIPF")
    ('gi|3298468|dbj|BAA31520.1|', '3298468', ['GI:3298468', 'DDBJ:BAA31520.1'])

    >>> _parse_ncbi_fasta_header("sp|P05698|LYSC_HUMAN Lysozyme C")
    ('sp|P05698|LYSC_HUMAN', 'P05698', ['UniProtKB/Swiss-Prot:P05698'])

    >>> _parse_ncbi_fasta_header("pdb|1A2B|C some PDB chain")
    ('pdb|1A2B|C', '1A2B', ['PDB:1A2B'])

    >>> _parse_ncbi_fasta_header("plain_id no pipes")
    ('plain_id', 'plain_id', [])

    >>> _parse_ncbi_fasta_header("")
    ('', '', [])

    """
    if not title:
        return "", "", []

    try:
        first_word = title.split(None, 1)[0]
    except IndexError:
        return "", "", []

    fields = first_word.split("|")

    # If there are no pipes, this is a plain FASTA header
    if len(fields) <= 1:
        return first_word, first_word, []

    dbxrefs = []
    first_accession = None
    i = 0
    while i < len(fields):
        field = fields[i]
        if field in ncbi_fasta_header_fields:
            field_names, db_name = ncbi_fasta_header_fields[field]
            # The accession/id is the next field after the prefix
            if i + 1 < len(fields) and fields[i + 1]:
                accession = fields[i + 1]
                dbxrefs.append(f"{db_name}:{accession}")
                if first_accession is None:
                    first_accession = accession
            # Skip over the fields consumed by this identifier
            i += 1 + len(field_names)
        else:
            i += 1

    name = first_accession if first_accession else first_word
    return first_word, name, dbxrefs


class FastaNcbiIterator(SequenceIterator):
    """Parser for FASTA files with NCBI-style pipe-delimited headers.

    This parser extracts database cross-references (dbxrefs) from NCBI
    FASTA title lines such as ``>gi|186972394|gb|EU490707.1| description``.

    The ``id`` and ``description`` fields are set identically to the plain
    ``fasta`` parser (for lossless FASTA roundtripping), but the ``name``
    field is set to the first recognized accession and the ``dbxrefs``
    list is populated with ``"DB:accession"`` strings.

    Use ``SeqIO.parse(handle, "fasta-ncbi")`` to access this parser.

    Examples
    --------
    Parse an NCBI FASTA file and access cross-references:

    >>> from Bio import SeqIO
    >>> for record in SeqIO.parse("Fasta/ncbi_headers.fasta", "fasta-ncbi"):
    ...     if record.dbxrefs:
    ...         print(f"{record.name}: {record.dbxrefs}")
    ...
    3298468: ['GI:3298468', 'DDBJ:BAA31520.1']
    186972394: ['GI:186972394', 'GenBank:EU490707.1']
    P05698: ['UniProtKB/Swiss-Prot:P05698']
    1A2B: ['PDB:1A2B']
    NM_001301717.2: ['RefSeq:NM_001301717.2']
    A0A0C5B5G6: ['UniProtKB/TrEMBL:A0A0C5B5G6']
    CAA12345.6: ['EMBL:CAA12345.6']
    US: ['Patent:US']
    taxon: ['General database reference:taxon']

    """

    modes = "t"

    def __init__(
        self,
        source: _TextIOSource,
        alphabet: None = None,
    ) -> None:
        """Iterate over FASTA records with NCBI dbxref parsing.

        Arguments:
         - source - input stream opened in text mode, or a path to a file
         - alphabet - optional alphabet, not used. Leave as None.

        """
        if alphabet is not None:
            raise ValueError("The alphabet argument is no longer supported")
        super().__init__(source, fmt="Fasta")
        line = self.stream.readline()
        if not line:
            line = None
        elif not line.startswith(">"):
            raise ValueError(
                "Expected FASTA record starting with '>' character. "
                f"Got: '{line.rstrip()}'"
            )
        self._line = line

    def __next__(self):
        line = self._line
        if line is None:
            raise StopIteration
        title = line[1:].rstrip()
        lines = []
        for line in self.stream:
            if line[0] == ">":
                break
            lines.append(line)
        else:
            line = None
        self._line = line
        sequence = "".join(lines).encode().translate(None, b" \t\r\n")

        first_word, name, dbxrefs = _parse_ncbi_fasta_header(title)

        return SeqRecord(
            Seq(sequence),
            id=first_word,
            name=name,
            description=title,
            dbxrefs=dbxrefs,
        )


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
        line = self.stream.readline()
        if not line:
            line = None
        else:
            if not line.startswith(">"):
                raise ValueError(
                    """\
This FASTA file contains comments at the beginning of the file, which are not
allowed by the 'fasta' parser.

To parse this file, you have three options:

(1) Modify your FASTA file to remove such comments at the beginning of the
file.

(2) Use SeqIO.parse with the 'fasta-pearson' format instead of 'fasta'. This
format is consistent with the FASTA format defined by William Pearson's FASTA
aligner software. This format allows for comments before the first sequence;
lines starting with the ';' character anywhere in the file are also regarded
as comment lines and are ignored.

(3) Use the 'fasta-blast' format. This format regards any lines "starting with
'!', '#', or ';' as comment lines. The 'fasta-blast' format may be safer than
the 'fasta-pearson' format, as it explicitly indicates which lines are comments.
"""
                )
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
    """FASTA file writer."""

    modes = "t"

    def __init__(self, target, wrap=60, record2title=None):
        """Create a Fasta writer.

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

    @classmethod
    def to_string(cls, record):
        """Turn a SeqRecord into a FASTA formatted string, and return it."""
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
    """Class to write 2-line per record Fasta format files.

    This means we write the sequence information  without line
    wrapping, and will always write a blank line for an empty
    sequence.
    """

    def __init__(self, handle, record2title=None):
        """Create a 2-line per record Fasta writer.

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

    @classmethod
    def to_string(cls, record):
        """Return a string in FASTA format with the sequence as one line."""
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


def as_fasta(record):
    """Turn a SeqRecord into a FASTA formatted string."""
    warnings.warn(
        """\
FastaIO.as_fasta is deprecated.

Instead of

FastaIO.as_fasta(record)

please use

format(record, "fasta")
""",
        DeprecationWarning,
    )
    return FastaWriter.to_string(record)


def as_fasta_2line(record):
    """Turn a SeqRecord into a two-line FASTA formatted string."""
    warnings.warn(
        """\
FastaIO.as_fasta_2line is deprecated.

Instead of

FastaIO.as_fasta_2line(record)

please use

format(record, "fasta-2line")
""",
        DeprecationWarning,
    )
    return FastaTwoLineWriter.to_string(record)


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest(verbose=0)
