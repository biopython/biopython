# Copyright 2023 copyright Michiel de Hoon
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Code to work with XML output from BLAST."""

# fmt: off
# flake8: noqa


from Bio import StreamModeError


class Record(list):
    def __init__(self):
        self.query = None


class Records:

    def __init__(self, stream):
        from Bio.Blast._parser import XMLHandler
        handler = XMLHandler(stream)
        handler.read_header(self)
        self.handler = handler

    def __iter__(self):
        return self

    def __next__(self):
        return next(self.handler)


def parse(source):
    """Parse an XML file containing BLAST output and return a Bio.Blast.Records object.

    This returns an iterator object; iterating over it returns Bio.Blast.Record
    objects one by one.

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
    ...     print(record.query.id, record.query.description)
    ... 
    Query_1 gi|195230749:301-1383 Homo sapiens wingless-type MMTV integration site family member 2 (WNT2), transcript variant 1, mRNA
    Query_2 gi|325053704:108-1166 Homo sapiens wingless-type MMTV integration site family, member 3A (WNT3A), mRNA
    Query_3 gi|156630997:105-1160 Homo sapiens wingless-type MMTV integration site family, member 4 (WNT4), mRNA
    Query_4 gi|371502086:108-1205 Homo sapiens wingless-type MMTV integration site family, member 5A (WNT5A), transcript variant 2, mRNA
    Query_5 gi|53729353:216-1313 Homo sapiens wingless-type MMTV integration site family, member 6 (WNT6), mRNA
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
    >>> stream = open("Blast/xml_2900_blastn_001.xml", "rb")  # opened in binary mode
    >>> record = Blast.read(stream)
    >>> record.query.id
    'G26684.1'
    >>> record.query.description
    'human STS STS_D11570, sequence tagged site'
    >>> len(record)
    10
    >>> stream.close()

    Use the Bio.Blast.parse function if you want to read a file containing
    BLAST output for more than one query.
    """
    records = parse(source)
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
