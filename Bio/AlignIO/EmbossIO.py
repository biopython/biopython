# Copyright 2008-2013 by Peter Cock.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Bio.AlignIO support for "emboss" alignment output from EMBOSS tools.

You are expected to use this module via the Bio.AlignIO functions (or the
Bio.SeqIO functions if you want to work directly with the gapped sequences).

This module contains a parser for the EMBOSS pairs/simple file format, for
example from the alignret, water and needle tools.
"""

from __future__ import print_function

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from .Interfaces import AlignmentIterator, SequentialAlignmentWriter

__docformat__ = "restructuredtext en"


class EmbossWriter(SequentialAlignmentWriter):
    """Emboss alignment writer (WORK IN PROGRESS).

    Writes a simplfied version of the EMBOSS pairs/simple file format.
    A lot of the information their tools record in their headers is not
    available and is omitted.
    """

    def write_header(self):
        handle = self.handle
        handle.write("########################################\n")
        handle.write("# Program: Biopython\n")
        try:
            handle.write("# Report_file: %s\n" % handle.name)
        except AttributeError:
            pass
        handle.write("########################################\n")

    def write_footer(self):
        handle = self.handle
        handle.write("#---------------------------------------\n")
        handle.write("#---------------------------------------\n")

    def write_alignment(self, alignment):
        """Use this to write (another) single alignment to an open file."""
        handle = self.handle
        handle.write("#=======================================\n")
        handle.write("#\n")
        handle.write("# Aligned_sequences: %i\n" % len(alignment))
        for i, record in enumerate(alignment):
            handle.write("# %i: %s\n" % (i + 1, record.id))
        handle.write("#\n")
        handle.write("# Length: %i\n" % alignment.get_alignment_length())
        handle.write("#\n")
        handle.write("#=======================================\n")
        handle.write("\n")
        # ...
        assert False


class EmbossIterator(AlignmentIterator):
    """Emboss alignment iterator.

    For reading the (pairwise) alignments from EMBOSS tools in what they
    call the "pairs" and "simple" formats.
    """

    def __next__(self):

        handle = self.handle

        try:
            # Header we saved from when we were parsing
            # the previous alignment.
            line = self._header
            del self._header
        except AttributeError:
            line = handle.readline()
        if not line:
            raise StopIteration

        while line.rstrip() != "#=======================================":
            line = handle.readline()
            if not line:
                raise StopIteration

        length_of_seqs = None
        number_of_seqs = None
        ids = []
        seqs = []

        while line[0] == "#":
            # Read in the rest of this alignment header,
            # try and discover the number of records expected
            # and their length
            parts = line[1:].split(":", 1)
            key = parts[0].lower().strip()
            if key == "aligned_sequences":
                number_of_seqs = int(parts[1].strip())
                assert len(ids) == 0
                # Should now expect the record identifiers...
                for i in range(number_of_seqs):
                    line = handle.readline()
                    parts = line[1:].strip().split(":", 1)
                    assert i + 1 == int(parts[0].strip())
                    ids.append(parts[1].strip())
                assert len(ids) == number_of_seqs
            if key == "length":
                length_of_seqs = int(parts[1].strip())

            # And read in another line...
            line = handle.readline()

        if number_of_seqs is None:
            raise ValueError("Number of sequences missing!")
        if length_of_seqs is None:
            raise ValueError("Length of sequences missing!")

        if self.records_per_alignment is not None \
        and self.records_per_alignment != number_of_seqs:
            raise ValueError("Found %i records in this alignment, told to expect %i"
                             % (number_of_seqs, self.records_per_alignment))

        seqs = ["" for id in ids]
        seq_starts = []
        index = 0

        # Parse the seqs
        while line:
            if len(line) > 21:
                id_start = line[:21].strip().split(None, 1)
                seq_end = line[21:].strip().split(None, 1)
                if len(id_start) == 2 and len(seq_end) == 2:
                    # identifier, seq start position, seq, seq end position
                    # (an aligned seq is broken up into multiple lines)
                    id, start = id_start
                    seq, end = seq_end
                    if start == end:
                        # Special case, either a single letter is present,
                        # or no letters at all.
                        if seq.replace("-", "") == "":
                            start = int(start)
                            end = int(end)
                        else:
                            start = int(start) - 1
                            end = int(end)
                    else:
                        assert seq.replace("-", "") != "", repr(line)
                        start = int(start) - 1  # python counting
                        end = int(end)

                    # The identifier is truncated...
                    assert 0 <= index and index < number_of_seqs, \
                           "Expected index %i in range [0,%i)" \
                           % (index, number_of_seqs)
                    assert id == ids[index] or id == ids[index][:len(id)]

                    if len(seq_starts) == index:
                        # Record the start
                        seq_starts.append(start)

                    # Check the start...
                    if start == end:
                        assert seq.replace("-", "") == "", line
                    else:
                        assert start - seq_starts[index] == len(seqs[index].replace("-", "")), \
                        "Found %i chars so far for sequence %i (%s, %s), line says start %i:\n%s" \
                            % (len(seqs[index].replace("-", "")), index, id, repr(seqs[index]),
                               start, line)

                    seqs[index] += seq

                    # Check the end ...
                    assert end == seq_starts[index] + len(seqs[index].replace("-", "")), \
                        "Found %i chars so far for sequence %i (%s, %s, start=%i), file says end %i:\n%s" \
                            % (len(seqs[index].replace("-", "")), index, id, repr(seqs[index]),
                               seq_starts[index], end, line)

                    index += 1
                    if index >= number_of_seqs:
                        index = 0
                else:
                    # just a start value, this is just alignment annotation (?)
                    # print "Skipping: " + line.rstrip()
                    pass
            elif line.strip() == "":
                # Just a spacer?
                pass
            else:
                print(line)
                assert False

            line = handle.readline()
            if line.rstrip() == "#---------------------------------------" \
            or line.rstrip() == "#=======================================":
                # End of alignment
                self._header = line
                break

        assert index == 0

        if self.records_per_alignment is not None \
        and self.records_per_alignment != len(ids):
            raise ValueError("Found %i records in this alignment, told to expect %i"
                             % (len(ids), self.records_per_alignment))

        records = []
        for id, seq in zip(ids, seqs):
            if len(seq) != length_of_seqs:
                # EMBOSS 2.9.0 is known to use spaces instead of minus signs
                # for leading gaps, and thus fails to parse.  This old version
                # is still used as of Dec 2008 behind the EBI SOAP webservice:
                # http://www.ebi.ac.uk/Tools/webservices/wsdl/WSEmboss.wsdl
                raise ValueError("Error parsing alignment - sequences of "
                                 "different length? You could be using an "
                                 "old version of EMBOSS.")
            records.append(SeqRecord(Seq(seq, self.alphabet),
                                     id=id, description=id))
        return MultipleSeqAlignment(records, self.alphabet)
