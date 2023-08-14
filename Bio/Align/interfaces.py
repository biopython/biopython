# Copyright 2006-2021 by Peter Cock.
# Copyright 2022 by Michiel de Hoon.
# All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.Align support module (not for general use).

Unless you are writing a new parser or writer for Bio.Align, you should not
use this module.  It provides base classes to try and simplify things.
"""
from abc import ABC
from abc import abstractmethod
from typing import Optional

from Bio import StreamModeError
from Bio.Align import AlignmentsAbstractBaseClass


class AlignmentIterator(AlignmentsAbstractBaseClass):
    """Base class for building Alignment iterators.

    You should write a parse method that returns an Alignment generator.  You
    may wish to redefine the __init__ method as well.

    Subclasses may define the following class attributes:
    - mode   - 't' or 'b' for text or binary files, respectively
    - fmt    - a human-readable name for the file format.
    """

    mode = "t"  # assume text files by default
    fmt: Optional[str] = None  # to be defined in the subclass

    def __init__(self, source):
        """Create an AlignmentIterator object.

        Arguments:
        - source - input file stream, or path to input file

        This method MAY be overridden by any subclass.

        Note when subclassing:
        - there should be a single non-optional argument, the source.
        - you can add additional optional arguments.
        """
        self.source = source
        try:
            self._stream = open(source, "r" + self.mode)
        except TypeError:  # not a path, assume we received a stream
            if self.mode == "t":
                if source.read(0) != "":
                    raise StreamModeError(
                        f"{self.fmt} files must be opened in text mode."
                    ) from None
            elif self.mode == "b":
                if source.read(0) != b"":
                    raise StreamModeError(
                        f"{self.fmt} files must be opened in binary mode."
                    ) from None
            else:
                raise ValueError(f"Unknown mode '{self.mode}'") from None
            self._stream = source
        self._read_header(self._stream)

    def __next__(self):
        """Return the next alignment."""
        try:
            stream = self._stream
        except AttributeError:
            raise StopIteration from None
        alignment = self._read_next_alignment(stream)
        if alignment is None:
            raise StopIteration
        return alignment

    def __len__(self):
        """Return the number of alignments.

        The number of alignments is cached. If not yet calculated, the iterator
        is rewound to the beginning, and the number of alignments is calculated
        by iterating over the alignments. The iterator is then returned to its
        original position in the file.
        """
        try:
            length = self._len
        except AttributeError:
            counter = 0
            for alignment in self:
                counter += 1
            self.rewind()
            length = 0
            for alignment in self:
                length += 1
            self.rewind()
            for i in range(length - counter):
                next(self)
            self._len = length
        return length

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        try:
            stream = self._stream
        except AttributeError:
            return
        if stream is not self.source:
            stream.close()
        del self._stream

    def _read_header(self, stream):
        """Read the file header and store it in metadata."""
        return

    @abstractmethod
    def _read_next_alignment(self, stream):
        """Read one Alignment from the stream, and return it."""

    def rewind(self):  # noqa: D102
        self._stream.seek(0)
        self._read_header(self._stream)


class AlignmentWriter(ABC):
    """Base class for alignment writers. This class should be subclassed.

    It is intended for alignment file formats with an (optional)
    header, one or more alignments, and an (optional) footer.

    The user may call the write_file() method to write a complete
    file containing the alignments.

    Alternatively, users may call the write_header(), followed
    by multiple calls to format_alignment() and/or write_alignments(),
    followed finally by write_footer().

    Subclasses may define the following class attributes:
    - mode   - 'w' or 'wb' for text or binary files, respectively
    - fmt    - a human-readable name for the file format.
    """

    mode = "w"  # assume text files by default
    fmt: Optional[str] = None  # to be defined in the subclass

    def __init__(self, target):
        """Create the writer object.

        Arguments:
        - target - output file stream, or path to output file

        This method MAY be overridden by any subclass.

        Note when subclassing:
        - there should be a single non-optional argument, the target.
        - you can add additional optional arguments.
        """
        if target is not None:
            # target is None if we only use the writer to format strings.
            if self.mode == "w":
                try:
                    target.write("")
                except TypeError:
                    # target was opened in binary mode
                    raise StreamModeError("File must be opened in text mode.") from None
                except AttributeError:
                    # target is a path
                    stream = open(target, self.mode)
                else:
                    stream = target
            elif self.mode == "wb":
                try:
                    target.write(b"")
                except TypeError:
                    # target was opened in text mode
                    raise StreamModeError(
                        "File must be opened in binary mode."
                    ) from None
                except AttributeError:
                    # target is a path
                    stream = open(target, self.mode)
                else:
                    stream = target
            else:
                raise RuntimeError("Unknown mode '%s'" % self.mode)
            self._stream = stream

        self._target = target

    def write_header(self, stream, alignments):
        """Write the file header to the output file."""
        return
        ##################################################
        # You MUST implement this method in the subclass #
        # if the file format defines a file header.      #
        ##################################################

    def write_footer(self, stream):
        """Write the file footer to the output file."""
        return
        ##################################################
        # You MUST implement this method in the subclass #
        # if the file format defines a file footer.      #
        ##################################################

    def format_alignment(self, alignment):
        """Format a single alignment as a string.

        alignment - an Alignment object
        """
        raise NotImplementedError("This method should be implemented")
        ###################################################
        # You MUST implement this method in the subclass. #
        ###################################################

    def write_single_alignment(self, stream, alignments):
        """Write a single alignment to the output file, and return 1.

        alignments - A list or iterator returning Alignment objects
        stream     - Output file stream.
        """
        count = 0
        for alignment in alignments:
            if count == 1:
                raise ValueError(
                    f"Alignment files in the {self.fmt} format can contain a single alignment only."
                )
            line = self.format_alignment(alignment)
            stream.write(line)
            count += 1
        return count

    def write_multiple_alignments(self, stream, alignments):
        """Write alignments to the output file, and return the number of alignments.

        alignments - A list or iterator returning Alignment objects
        stream     - Output file stream.
        """
        count = 0
        for alignment in alignments:
            line = self.format_alignment(alignment)
            stream.write(line)
            count += 1
        return count

    write_alignments = write_multiple_alignments

    def write_file(self, stream, alignments):
        """Write the alignments to the file strenm, and return the number of alignments.

        alignments - A list or iterator returning Alignment objects
        stream     - Output file stream.
        """
        self.write_header(stream, alignments)
        count = self.write_alignments(stream, alignments)
        self.write_footer(stream)
        return count

    def write(self, alignments):
        """Write a file with the alignments, and return the number of alignments.

        alignments - A list or iterator returning Alignment objects
        """
        stream = self._stream
        try:
            count = self.write_file(stream, alignments)
        finally:
            if stream is not self._target:
                stream.close()
        return count
