# Copyright 2006-2021 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.SeqIO support module (not for general use).

Unless you are writing a new parser or writer for Bio.SeqIO, you should not
use this module.  It provides base classes to try and simplify things.
"""
from abc import ABC
from abc import abstractmethod

from Bio import StreamModeError
from Bio.Seq import MutableSeq
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class SequenceIterator(ABC):
    """Base class for building SeqRecord iterators.

    You should write a parse method that returns a SeqRecord generator.  You
    may wish to redefine the __init__ method as well.
    """

    def __init__(self, source, alphabet=None, mode="t", fmt=None):
        """Create a SequenceIterator object.

        Arguments:
        - source - input file stream, or path to input file
        - alphabet - no longer used, should be None

        This method MAY be overridden by any subclass.

        Note when subclassing:
        - there should be a single non-optional argument, the source.
        - you do not have to require an alphabet.
        - you can add additional optional arguments.
        """
        if alphabet is not None:
            raise ValueError("The alphabet argument is no longer supported")
        try:
            self.stream = open(source, "r" + mode)
            self.should_close_stream = True
        except TypeError:  # not a path, assume we received a stream
            if mode == "t":
                if source.read(0) != "":
                    raise StreamModeError(
                        f"{fmt} files must be opened in text mode."
                    ) from None
            elif mode == "b":
                if source.read(0) != b"":
                    raise StreamModeError(
                        f"{fmt} files must be opened in binary mode."
                    ) from None
            else:
                raise ValueError(f"Unknown mode '{mode}'") from None
            self.stream = source
            self.should_close_stream = False
        try:
            self.records = self.parse(self.stream)
        except Exception:
            if self.should_close_stream:
                self.stream.close()
            raise

    def __next__(self):
        """Return the next entry."""
        try:
            return next(self.records)
        except Exception:
            if self.should_close_stream:
                self.stream.close()
            raise

    def __iter__(self):
        """Iterate over the entries as a SeqRecord objects.

        Example usage for Fasta files::

            with open("example.fasta","r") as myFile:
                myFastaReader = FastaIterator(myFile)
                for record in myFastaReader:
                    print(record.id)
                    print(record.seq)

        This method SHOULD NOT be overridden by any subclass. It should be
        left as is, which will call the subclass implementation of __next__
        to actually parse the file.
        """
        return self

    @abstractmethod
    def parse(self, handle):
        """Start parsing the file, and return a SeqRecord iterator."""


def _get_seq_string(record):
    """Use this to catch errors like the sequence being None (PRIVATE)."""
    if not isinstance(record, SeqRecord):
        raise TypeError("Expected a SeqRecord object")
    if record.seq is None:
        raise TypeError(f"SeqRecord (id={record.id}) has None for its sequence.")
    elif not isinstance(record.seq, (Seq, MutableSeq)):
        raise TypeError(f"SeqRecord (id={record.id}) has an invalid sequence.")
    return str(record.seq)


# Function variant of the SequenceWriter method.
def _clean(text):
    """Use this to avoid getting newlines in the output (PRIVATE)."""
    return text.replace("\n", " ").replace("\r", " ")


class SequenceWriter:
    """Base class for sequence writers. This class should be subclassed.

    It is intended for sequential file formats with an (optional)
    header, repeated records, and an (optional) footer, as well
    as for interlaced file formats such as Clustal.

    The user may call the write_file() method to write a complete
    file containing the sequences.

    Alternatively, users may call the write_header(), followed
    by multiple calls to write_record() and/or write_records(),
    followed finally by write_footer().

    Note that write_header() cannot require any assumptions about
    the number of records.
    """

    def __init__(self, target, mode="w"):
        """Create the writer object."""
        if mode == "w":
            try:
                target.write("")
            except TypeError:
                # target was opened in binary mode
                raise StreamModeError("File must be opened in text mode.") from None
            except AttributeError:
                # target is a path
                handle = open(target, mode)
            else:
                handle = target
        elif mode == "wb":
            try:
                target.write(b"")
            except TypeError:
                # target was opened in text mode
                raise StreamModeError("File must be opened in binary mode.") from None
            except AttributeError:
                # target is a path
                handle = open(target, mode)
            else:
                handle = target
        else:
            raise RuntimeError(f"Unknown mode '{mode}'")

        self._target = target
        self.handle = handle

    def clean(self, text):
        """Use this to avoid getting newlines in the output."""
        return text.replace("\n", " ").replace("\r", " ")

    def write_header(self):
        """Write the file header to the output file."""
        pass
        ##################################################
        # You MUST implement this method in the subclass #
        # if the file format defines a file header.      #
        ##################################################

    def write_footer(self):
        """Write the file footer to the output file."""
        pass
        ##################################################
        # You MUST implement this method in the subclass #
        # if the file format defines a file footer.      #
        ##################################################

    def write_record(self, record):
        """Write a single record to the output file.

        record - a SeqRecord object
        """
        raise NotImplementedError("This method should be implemented")
        ##################################################
        # You MUST implement this method in the subclass #
        # for sequential file formats.                   #
        ##################################################

    def write_records(self, records, maxcount=None):
        """Write records to the output file, and return the number of records.

        records - A list or iterator returning SeqRecord objects
        maxcount - The maximum number of records allowed by the
        file format, or None if there is no maximum.
        """
        count = 0
        if maxcount is None:
            for record in records:
                self.write_record(record)
                count += 1
        else:
            for record in records:
                if count == maxcount:
                    if maxcount == 1:
                        raise ValueError("More than one sequence found")
                    else:
                        raise ValueError(
                            "Number of sequences is larger than %d" % maxcount
                        )
                self.write_record(record)
                count += 1
        return count

    def write_file(self, records, mincount=0, maxcount=None):
        """Write a complete file with the records, and return the number of records.

        records - A list or iterator returning SeqRecord objects
        """
        ##################################################
        # You MUST implement this method in the subclass #
        # for interlaced file formats.                   #
        ##################################################
        try:
            self.write_header()
            count = self.write_records(records, maxcount)
            self.write_footer()
        finally:
            if self.handle is not self._target:
                self.handle.close()
        if count < mincount:
            if mincount == 1:  # Common case
                raise ValueError("Must have one sequence")
            elif mincount == maxcount:
                raise ValueError(
                    "Number of sequences is %d (expected %d)" % (count, mincount)
                )
            else:
                raise ValueError(
                    "Number of sequences is %d (expected at least %d)"
                    % (count, mincount)
                )
        return count
