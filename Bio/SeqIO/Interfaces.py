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
from os import PathLike
from typing import AnyStr
from typing import Generic
from typing import IO
from typing import Optional
from typing import Union

from Bio import StreamModeError
from Bio.Seq import MutableSeq
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# https://docs.python.org/3/glossary.html#term-path-like-object
_PathLikeTypes = (PathLike, str, bytes)
_IOSource = Union[IO[AnyStr], PathLike, str, bytes]
_TextIOSource = _IOSource[str]
_BytesIOSource = _IOSource[bytes]


class SequenceIterator(ABC, Generic[AnyStr]):
    """Base class for building SeqRecord iterators.

    You should write a __next__ method that returns the next SeqRecord.  You
    may wish to redefine the __init__ method as well.
    You must also create a class property `modes` specifying the allowable
    file stream modes.
    """

    @property
    @abstractmethod
    def modes(self):
        """File modes (binary or text) that the parser can handle.

        This property must be "t" (for text mode only), "b" (for binary mode
        only), "tb" (if both text and binary mode are accepted, but text mode
        is preferred), or "bt" (if both text and binary mode are accepted, but
        binary mode is preferred).
        """
        pass

    def __init__(
        self,
        source: _IOSource,
        alphabet: None = None,
        fmt: Optional[str] = None,
    ) -> None:
        """Create a SequenceIterator object.

        Arguments:
        - source - input file stream, or path to input file
        - alphabet - no longer used, should be None
        - fmt - string, mixed case format name for in error messages

        This method MAY be overridden by any subclass.

        Note when subclassing:
        - there should be a single non-optional argument, the source.
        - you do not have to require an alphabet.
        - you can add additional optional arguments.
        """
        if alphabet is not None:
            raise ValueError("The alphabet argument is no longer supported")
        modes = self.modes
        self.source = source
        if isinstance(source, _PathLikeTypes):
            mode = modes[0]
            self.stream = open(source, "r" + mode)
        else:
            value = source.read(0)
            if value == "":
                if modes == "b":
                    raise StreamModeError(
                        f"{fmt} files must be opened in binary mode."
                    ) from None
                mode = "t"
            elif value == b"":
                if modes == "t":
                    raise StreamModeError(
                        f"{fmt} files must be opened in text mode."
                    ) from None
                mode = "b"
            else:
                raise RuntimeError("Failed to read from input data") from None
            self.stream = source
        self.mode = mode

    @abstractmethod
    def __next__(self):
        """Return the next SeqRecord.

        This method must be implemented by the subclass.
        """

    def __iter__(self):
        """Iterate over the entries as a SeqRecord objects.

        Example usage for Fasta files::

            with open("example.fasta","r") as myFile:
                myFastaReader = FastaIterator(myFile)
                for record in myFastaReader:
                    print(record.id)
                    print(record.seq)

        This method SHOULD NOT be overridden by any subclass.
        """
        return self

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        try:
            stream = self.stream
        except AttributeError:
            return
        if self.stream is not self.source:
            self.stream.close()
        del self.stream
        return False


def _get_seq_string(record: SeqRecord) -> str:
    """Use this to catch errors like the sequence being None (PRIVATE)."""
    if not isinstance(record, SeqRecord):
        raise TypeError("Expected a SeqRecord object")
    if record.seq is None:
        raise TypeError(f"SeqRecord (id={record.id}) has None for its sequence.")
    elif not isinstance(record.seq, (Seq, MutableSeq)):
        raise TypeError(f"SeqRecord (id={record.id}) has an invalid sequence.")
    return str(record.seq)


# Function variant of the SequenceWriter method.
def _clean(text: str) -> str:
    """Use this to avoid getting newlines in the output (PRIVATE)."""
    return text.replace("\n", " ").replace("\r", " ")


class SequenceWriter(ABC, Generic[AnyStr]):
    """Base class for sequence writers. This class should be subclassed.

    The user may call the write_file() method to write a complete
    file containing the sequences.

    Most subclasses will only need to implement the write_record method.
    Subclasses must implement the write_records method to include a file
    header or footer, for file formats that only allow one record, or for
    file formats that cannot be written sequentially.
    """

    @property
    @abstractmethod
    def modes(self):
        """File modes (binary or text) that the writer can handle.

        This property must be "t" (for text mode only), "b" (for binary mode
        only), "tb" (if both text and binary mode are accepted, but text mode
        is preferred), or "bt" (if both text and binary mode are accepted, but
        binary mode is preferred).
        """
        pass

    def __init__(self, target: _IOSource) -> None:
        """Create the writer object."""
        if isinstance(target, _PathLikeTypes):
            mode = self.modes[0]
            stream = open(target, "w" + mode)
        else:
            stream = target
            modes = "tb"
            values = ("", b"")
            for mode, value in zip(modes, values):
                try:
                    stream.write(value)
                except TypeError:
                    continue
                else:
                    break
            else:
                raise RuntimeError("Failed to read from input data") from None
            if mode not in self.modes:
                if mode == "t":
                    # target was opened in text mode
                    raise StreamModeError("File must be opened in binary mode.")
                elif mode == "b":
                    # target was opened in binary mode
                    raise StreamModeError("File must be opened in text mode.")
        self.target = target
        self.handle = stream

    def clean(self, text: str) -> str:
        """Use this to avoid getting newlines in the output."""
        return text.replace("\n", " ").replace("\r", " ")

    def write_record(self, record):
        """Write a single record to the output file.

        record - a SeqRecord object
        """

    def write_records(self, records):
        """Write records to the output file, and return the number of records.

        records - A list or iterator returning SeqRecord objects
        """
        count = 0
        for record in records:
            self.write_record(record)
            count += 1
        return count

    def write_file(self, records):
        """Write a complete file with the records, and return the number of records.

        records - A list or iterator returning SeqRecord objects
        """
        try:
            count = self.write_records(records)
        finally:
            if self.handle is not self.target:
                self.handle.close()
        return count
