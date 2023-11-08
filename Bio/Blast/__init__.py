# Copyright 1999 by Jeffrey Chang.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Code for dealing with BLAST programs and output."""


# fmt: off
# ruff: noqa
# flake8: noqa
# mypy: ignore-errors


from Bio import StreamModeError


class Records(list):
    def __init__(self, stream):
        from ._parser import DataHandler
        handler = DataHandler()
        self.header = handler.read_header(stream)
        self._handler = handler
        self._stream = stream

    def __next__(self):
        try:
            handler = self._handler
        except AttributeError:
            raise StopIteration from None
        stream = self._stream
        record = handler.read_next_record(stream)
        if record is None:
            del self._handler
            raise StopIteration
        return record


def parse(source):
    """Parse an XML file containing BLAST output and return a Bio.Blast.Records object.

    This is a generator function that returns BLAST records one by one.

    The source can be a file stream or the path to an XML file containing the
    BLAST output. If a file stream, source  must be in binary mode. This allows
    the parser to detect the encoding from the XML file,and to use it to convert
    any text in the XML to the correct Unicode string. The qblast function in
    Bio.Blast returns a file stream in binary mode. For files, please use mode
    "rb" when opening the file, as in

    >>> from Bio import Blast
    >>> stream = open("Blast/wnts.xml", "rb")  # opened in binary mode
    >>> records = Blast.parse(stream)
    >>> for record in records:
    ...     
    ...
    >>> stream.close()

    """
    try:
        stream = open(source, "rb")
    except TypeError:  # not a path, assume we received a stream
        if source.read(0) != b"":
            raise StreamModeError(
                f"BLAST output files must be opened in binary mode."
            ) from None
        stream = source

    try:
        return Records(stream)
    finally:
        if stream != source:
            stream.close()


def read(source):
    """Parse an XML file containing BLAST output for a single query and return it.

    Internally, this function uses Bio.Blast.parse to obtain an iterator over
    BLAST records.  The function then reads one record from the iterator,
    ensures that there are no further records, and returns the record it found
    as a Bio.Blast.Record object. An exception is raised if no records are
    found, or more than one record is found.

    The source can be a file stream or the path to an XML file containing the
    BLAST output. If a file stream, source  must be in binary mode. This allows
    the parser to detect the encoding from the XML file,and to use it to convert
    any text in the XML to the correct Unicode string. The qblast function in
    Bio.Blast returns a file stream in binary mode. For files, please use mode
    "rb" when opening the file, as in

    >>> from Bio import Blast
    >>> stream = open("Blast/wnts.xml", "rb")  # opened in binary mode
    >>> record = Blast.read(stream)
    >>> print(record)
    >>> stream.close()

    Use the Bio.Blast.parse function if you want to read a file containing
    BLAST output for more than one query.
    """
    records = parse(handle, fmt)
    try:
        record = next(records)
    except StopIteration:
        raise ValueError("No BLAST output found.") from None
    try:
        next(records)
        raise ValueError("BLAST output for more than one query found.")
    except StopIteration:
        pass
    return record
