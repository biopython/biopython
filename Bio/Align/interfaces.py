# Copyright 2006-2021 by Peter Cock.  All rights reserved.
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

from Bio import StreamModeError


class AlignmentIterator(list, ABC):
    """Base class for building Alignment iterators.

    You should write a parse method that returns an Alignment generator.  You
    may wish to redefine the __init__ method as well.
    """

    def __init__(self, source, mode="t", fmt=None):
        """Create an AlignmentIterator object.

        Arguments:
        - source - input file stream, or path to input file

        This method MAY be overridden by any subclass.

        Note when subclassing:
        - there should be a single non-optional argument, the source.
        - you can add additional optional arguments.
        """
        self.source = source
        if source is None:
            return
        self._loaded = False
        try:
            self._stream = open(source, "r" + mode)
        except TypeError:  # not a path, assume we received a stream
            if mode == "t":
                if source.read(0) != "":
                    raise StreamModeError(
                        "%s files must be opened in text mode." % fmt
                    ) from None
            elif mode == "b":
                if source.read(0) != b"":
                    raise StreamModeError(
                        "%s files must be opened in binary mode." % fmt
                    ) from None
            else:
                raise ValueError("Unknown mode '%s'" % mode) from None
            self._stream = source
        try:
            self._read_header(self._stream)
        except Exception:
            self._close()
            raise

    def __next__(self):
        """Return the next entry."""
        if self._loaded:
            index = self._index
            length = super().__len__()
            if index == length:
                raise StopIteration from None
            self._index += 1
            return super().__getitem__(index)
        try:
            stream = self._stream
        except AttributeError:
            raise StopIteration from None
        try:
            alignment = self._read_next_alignment(stream)
            if alignment is None:
                raise StopIteration
        except Exception:
            self._close()
            raise
        return alignment

    def __iter__(self):
        """Iterate over the alignments as Alignment objects.

        This method SHOULD NOT be overridden by any subclass. It should be
        left as is, which will call the subclass implementation of __next__
        to actually parse the file.
        """
        if self._loaded:
            self._index = 0
        return self

    def _read_header(self, stream):
        """Read the file header and store it in metadata."""

    @abstractmethod
    def _read_next_alignment(self, stream):
        """Read one Alignment from the stream, and return it."""

    def _close(self):
        try:
            stream = self._stream
        except AttributeError:
            return
        if stream is not self.source:
            stream.close()
        del self._stream

    def _load(self):
        if not self._loaded:
            for item in self:
                super().append(item)
            self._loaded = True
            self._index = 0  # for use by the iterator

    def __repr__(self):
        return super(list, self).__repr__()

    def __lt__(self, other):
        self._load()
        if isinstance(other, AlignmentIterator):
            other._load()
        return super().__lt__(other)

    def __le__(self, other):
        self._load()
        if isinstance(other, AlignmentIterator):
            other._load()
        return super().__le__(other)

    def __eq__(self, other):
        self._load()
        if isinstance(other, AlignmentIterator):
            other._load()
        return super().__eq__(other)

    def __gt__(self, other):
        self._load()
        if isinstance(other, AlignmentIterator):
            other._load()
        return super().__gt__(other)

    def __ge__(self, other):
        self._load()
        if isinstance(other, AlignmentIterator):
            other._load()
        return super().__ge__(other)

    def __contains__(self, alignment):
        self._load()
        return super().__contains__(alignment)

    def __len__(self):
        self._load()
        return super().__len__()

    def __getitem__(self, i):
        self._load()
        if isinstance(i, slice):
            alignments = AlignmentIterator(None)
            items = super().__getitem__(i)
            list.extend(alignments, items)
            alignments._loaded = True
            return alignments
        else:
            return super().__getitem__(i)

    def __setitem__(self, i, item):
        self._load()
        super().__setitem__(i, item)

    def __delitem__(self, i):
        self._load()
        super().__delitem__(i)

    def __add__(self, other):
        alignments = AlignmentIterator(None)
        list.extend(alignments, self)
        list.extend(alignments, other)
        alignments._loaded = True
        return alignments

    def __radd__(self, other):
        alignments = AlignmentIterator(None)
        list.extend(alignments, other)
        list.extend(alignments, self)
        alignments._loaded = True
        return alignments

    def __iadd__(self, other):
        self.extend(other)
        return self

    def __mul__(self, n):
        alignments = AlignmentIterator(None)
        list.extend(alignments, list(self) * n)
        alignments._loaded = True
        return alignments

    def __rmul__(self, n):
        alignments = AlignmentIterator(None)
        list.extend(alignments, n * list(self))
        alignments._loaded = True
        return alignments

    def __imul__(self, n):
        self._load()
        return super().__imul__(n)

    def append(self, item):
        self._load()
        super().append(item)

    def insert(self, i, item):
        self._load()
        super().insert(i, item)

    def pop(self, i=-1):
        self._load()
        return super().pop(i)

    def remove(self, item):
        self._load()
        super().remove(item)

    def clear(self):
        self._close()
        super().clear()

    def copy(self):
        self._load()
        inst = self.__class__.__new__(self.__class__)
        inst.__dict__.update(self.__dict__)
        return inst

    def count(self, item):
        self._load()
        return super().count(item)

    def index(self, item, *args):
        self._load()
        return super().index(item, *args)

    def reverse(self):
        self._load()
        super().reverse()

    def sort(self, /, *args, **kwds):
        self._load()
        super().sort(*args, **kwds)

    def extend(self, other):
        self._load()
        for item in other:
            super().append(item)


class AlignmentWriter:
    """Base class for alignment writers. This class should be subclassed.

    It is intended for sequential file formats with an (optional)
    header, one or more alignments, and an (optional) footer.

    The user may call the write_file() method to write a complete
    file containing the alignments.

    Alternatively, users may call the write_header(), followed
    by multiple calls to format_alignment() and/or write_alignments(),
    followed finally by write_footer().

    Note that write_header() cannot require any assumptions about
    the number of alignments.
    """

    def __init__(self, target, mode="w"):
        """Create the writer object."""
        if target is not None:
            # target is None if we only use the writer to format strings.
            if mode == "w":
                try:
                    target.write("")
                except TypeError:
                    # target was opened in binary mode
                    raise StreamModeError("File must be opened in text mode.") from None
                except AttributeError:
                    # target is a path
                    stream = open(target, mode)
                else:
                    stream = target
            elif mode == "wb":
                try:
                    target.write(b"")
                except TypeError:
                    # target was opened in text mode
                    raise StreamModeError(
                        "File must be opened in binary mode."
                    ) from None
                except AttributeError:
                    # target is a path
                    stream = open(target, mode)
                else:
                    stream = target
            else:
                raise RuntimeError("Unknown mode '%s'" % mode)
            self.stream = stream

        self._target = target

    def write_header(self, alignments):
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

    def format_alignment(self, alignment):
        """Format a single alignment as a string.

        alignment - an Alignment object
        """
        raise NotImplementedError("This method should be implemented")
        ###################################################
        # You MUST implement this method in the subclass. #
        ###################################################

    def write_alignments(self, alignments, maxcount=None):
        """Write alignments to the output file, and return the number of alignments.

        alignments - A list or iterator returning Alignment objects
        maxcount - The maximum number of alignments allowed by the
        file format, or None if there is no maximum.
        """
        count = 0
        if maxcount is None:
            for alignment in alignments:
                line = self.format_alignment(alignment)
                self.stream.write(line)
                count += 1
        else:
            for alignment in alignments:
                if count == maxcount:
                    if maxcount == 1:
                        raise ValueError("More than one alignment found")
                    else:
                        raise ValueError(
                            "Number of alignments is larger than %d" % maxcount
                        )
                line = self.format_alignment(alignment)
                self.stream.write(line)
                count += 1
        return count

    def write_file(self, alignments, mincount=0, maxcount=None):
        """Write a file with the alignments, and return the number of alignments.

        alignments - A list or iterator returning Alignment objects
        """
        try:
            self.write_header(alignments)
            count = self.write_alignments(alignments, maxcount)
            self.write_footer()
        finally:
            if self.stream is not self._target:
                self.stream.close()
        if count < mincount:
            if mincount == 1:  # Common case
                raise ValueError("Must have one alignment")
            elif mincount == maxcount:
                raise ValueError(
                    "Number of alignments is %d (expected %d)" % (count, mincount)
                )
            else:
                raise ValueError(
                    "Number of alignmnets is %d (expected at least %d)"
                    % (count, mincount)
                )
        return count
