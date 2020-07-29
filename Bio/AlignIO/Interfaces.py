# Copyright 2008-2018 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""AlignIO support module (not for general use).

Unless you are writing a new parser or writer for Bio.AlignIO, you should not
use this module.  It provides base classes to try and simplify things.
"""


class AlignmentIterator:
    """Base class for building MultipleSeqAlignment iterators.

    You should write a next() method to return Alignment
    objects.  You may wish to redefine the __init__
    method as well.
    """

    def __init__(self, handle, seq_count=None):
        """Create an AlignmentIterator object.

        Arguments:
         - handle   - input file
         - count    - optional, expected number of records per alignment
           Recommend for fasta file format.

        Note when subclassing:
         - there should be a single non-optional argument, the handle,
           and optional count IN THAT ORDER.
         - you can add additional optional arguments.

        """
        self.handle = handle
        self.records_per_alignment = seq_count
        #####################################################
        # You may want to subclass this, for example        #
        # to read through the file to find the first record,#
        # or if additional arguments are required.          #
        #####################################################

    def __next__(self):
        """Return the next alignment in the file.

        This method should be replaced by any derived class to do something
        useful.
        """
        raise NotImplementedError("This object should be subclassed")
        #####################################################
        # You SHOULD subclass this, to split the file up    #
        # into your individual alignments and convert these #
        # into MultipleSeqAlignment objects.                #
        #####################################################

    def __iter__(self):
        """Iterate over the entries as MultipleSeqAlignment objects.

        Example usage for (concatenated) PHYLIP files::

            with open("many.phy","r") as myFile:
                for alignment in PhylipIterator(myFile):
                    print "New alignment:"
                    for record in alignment:
                        print record.id
                        print record.seq

        """
        return iter(self.__next__, None)


class AlignmentWriter:
    """Base class for building MultipleSeqAlignment writers.

    You should write a write_alignment() method.
    You may wish to redefine the __init__ method as well.
    """

    def __init__(self, handle):
        """Initialize the class."""
        self.handle = handle

    def write_file(self, alignments):
        """Use this to write an entire file containing the given alignments.

        Arguments:
         - alignments - A list or iterator returning MultipleSeqAlignment objects

        In general, this method can only be called once per file.

        This method should be replaced by any derived class to do something
        useful.  It should return the number of alignments..
        """
        raise NotImplementedError("This object should be subclassed")
        #####################################################
        # You SHOULD subclass this, to write the alignment  #
        # objects to the file handle                        #
        #####################################################

    def clean(self, text):
        """Use this to avoid getting newlines in the output."""
        return text.replace("\n", " ").replace("\r", " ")


class SequentialAlignmentWriter(AlignmentWriter):
    """Base class for building MultipleSeqAlignment writers.

    This assumes each alignment can be simply appended to the file.
    You should write a write_alignment() method.
    You may wish to redefine the __init__ method as well.
    """

    def __init__(self, handle):
        """Initialize the class."""
        self.handle = handle

    def write_file(self, alignments):
        """Use this to write an entire file containing the given alignments.

        Arguments:
         - alignments - A list or iterator returning MultipleSeqAlignment objects

        In general, this method can only be called once per file.
        """
        self.write_header()
        count = 0
        for alignment in alignments:
            self.write_alignment(alignment)
            count += 1
        self.write_footer()
        return count

    def write_header(self):
        """Use this to write any header.

        This method should be replaced by any derived class to do something
        useful.
        """
        pass

    def write_footer(self):
        """Use this to write any footer.

        This method should be replaced by any derived class to do something
        useful.
        """
        pass

    def write_alignment(self, alignment):
        """Use this to write a single alignment.

        This method should be replaced by any derived class to do something
        useful.
        """
        raise NotImplementedError("This object should be subclassed")
        #####################################################
        # You SHOULD subclass this, to write the alignment  #
        # objects to the file handle                        #
        #####################################################
