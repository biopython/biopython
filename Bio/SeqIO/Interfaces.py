# Copyright 2006-2013 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Bio.SeqIO support module (not for general use).

Unless you are writing a new parser or writer for Bio.SeqIO, you should not
use this module.  It provides base classes to try and simplify things.
"""

from __future__ import print_function

import sys  # for checking if Python 2

from Bio.Alphabet import generic_alphabet
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord

__docformat__ = "restructuredtext en"


class SequenceIterator(object):
    """Base class for building SeqRecord iterators.

    You should write a next() method to return SeqRecord
    objects.  You may wish to redefine the __init__
    method as well.
    """
    def __init__(self, handle, alphabet=generic_alphabet):
        """Create a SequenceIterator object.

            - handle - input file
            - alphabet - optional, e.g. Bio.Alphabet.generic_protein

        Note when subclassing:

            - there should be a single non-optional argument,
              the handle.
            - you do not have to require an alphabet.
            - you can add additional optional arguments."""
        self.handle = handle
        self.alphabet = alphabet
        #####################################################
        # You may want to subclass this, for example        #
        # to read through the file to find the first record,#
        # or if additional arguments are required.          #
        #####################################################

    def __next__(self):
        """Return the next record in the file.

        This method should be replaced by any derived class to do something useful."""
        raise NotImplementedError("This object should be subclassed")
        #####################################################
        # You SHOULD subclass this, to split the file up    #
        # into your individual records, and convert these   #
        # into useful objects, e.g. return SeqRecord object #
        #####################################################

    if sys.version_info[0] < 3:
        def next(self):
            """Python 2 style alias for Python 3 style __next__ method."""
            return self.__next__()

    def __iter__(self):
        """Iterate over the entries as a SeqRecord objects.

        Example usage for Fasta files::

            with open("example.fasta","r") as myFile:
                myFastaReader = FastaIterator(myFile)
                for record in myFastaReader:
                    print(record.id)
                    print(record.seq)
        """
        return iter(self.__next__, None)


class SequenceWriter(object):
    """This class should be subclassed.

    Interlaced file formats (e.g. Clustal) should subclass directly.

    Sequential file formats (e.g. Fasta, GenBank) should subclass
    the SequentialSequenceWriter class instead.
    """
    def __init__(self, handle):
        """Creates the writer object.

        Use the method write_file() to actually record your sequence records."""
        self.handle = handle

    def _get_seq_string(self, record):
        """Use this to catch errors like the sequence being None."""
        if not isinstance(record, SeqRecord):
            raise TypeError("Expected a SeqRecord object")
        if record.seq is None:
            raise TypeError("SeqRecord (id=%s) has None for its sequence."
                            % record.id)
        elif not isinstance(record.seq, (Seq, MutableSeq)):
            raise TypeError("SeqRecord (id=%s) has an invalid sequence."
                            % record.id)
        return str(record.seq)

    def clean(self, text):
        """Use this to avoid getting newlines in the output."""
        return text.replace("\n", " ").replace("\r", " ").replace("  ", " ")

    def write_file(self, records):
        """Use this to write an entire file containing the given records.

        records - A list or iterator returning SeqRecord objects

        Should return the number of records (as an integer).

        This method can only be called once."""
        # Note when implementing this, your writer class should NOT close the
        # file at the end, but the calling code should.
        raise NotImplementedError("This object should be subclassed")
        #####################################################
        # You SHOULD subclass this                          #
        #####################################################


class SequentialSequenceWriter(SequenceWriter):
    """This class should be subclassed.

    It is intended for sequential file formats with an (optional)
    header, repeated records, and an (optional) footer.

    In this case (as with interlaced file formats), the user may
    simply call the write_file() method and be done.

    However, they may also call the write_header(), followed
    by multiple calls to write_record() and/or write_records()
    followed finally by write_footer().

    Users must call write_header() and write_footer() even when
    the file format concerned doesn't have a header or footer.
    This is to try and make life as easy as possible when
    switching the output format.

    Note that write_header() cannot require any assumptions about
    the number of records.
    """
    def __init__(self, handle):
        self.handle = handle
        self._header_written = False
        self._record_written = False
        self._footer_written = False

    def write_header(self):
        assert not self._header_written, "You have aleady called write_header()"
        assert not self._record_written, "You have aleady called write_record() or write_records()"
        assert not self._footer_written, "You have aleady called write_footer()"
        self._header_written = True

    def write_footer(self):
        assert self._header_written, "You must call write_header() first"
        assert self._record_written, "You have not called write_record() or write_records() yet"
        assert not self._footer_written, "You have aleady called write_footer()"
        self._footer_written = True

    def write_record(self, record):
        """Write a single record to the output file.

        record - a SeqRecord object

        Once you have called write_header() you can call write_record()
        and/or write_records() as many times as needed.  Then call
        write_footer() and close()."""
        assert self._header_written, "You must call write_header() first"
        assert not self._footer_written, "You have already called write_footer()"
        self._record_written = True
        raise NotImplementedError("This object should be subclassed")
        #####################################################
        # You SHOULD subclass this                          #
        #####################################################

    def write_records(self, records):
        """Write multiple record to the output file.

        records - A list or iterator returning SeqRecord objects

        Once you have called write_header() you can call write_record()
        and/or write_records() as many times as needed.  Then call
        write_footer() and close().

        Returns the number of records written.
        """
        # Default implementation:
        assert self._header_written, "You must call write_header() first"
        assert not self._footer_written, "You have already called write_footer()"
        count = 0
        for record in records:
            self.write_record(record)
            count += 1
        # Mark as true, even if there where no records
        self._record_written = True
        return count

    def write_file(self, records):
        """Use this to write an entire file containing the given records.

        records - A list or iterator returning SeqRecord objects

        This method can only be called once.  Returns the number of records
        written.
        """
        self.write_header()
        count = self.write_records(records)
        self.write_footer()
        return count
