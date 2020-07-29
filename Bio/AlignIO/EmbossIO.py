# Copyright 2008-2016 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.AlignIO support for "emboss" alignment output from EMBOSS tools.

You are expected to use this module via the Bio.AlignIO functions (or the
Bio.SeqIO functions if you want to work directly with the gapped sequences).

This module contains a parser for the EMBOSS pairs/simple file format, for
example from the alignret, water and needle tools.
"""


from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.AlignIO.Interfaces import AlignmentIterator


class EmbossIterator(AlignmentIterator):
    """Emboss alignment iterator.

    For reading the (pairwise) alignments from EMBOSS tools in what they
    call the "pairs" and "simple" formats.
    """

    _header = None  # for caching lines between __next__ calls

    def __next__(self):
        """Parse the next alignment from the handle."""
        handle = self.handle

        if self._header is None:
            line = handle.readline()
        else:
            # Header we saved from when we were parsing
            # the previous alignment.
            line = self._header
            self._header = None

        if not line:
            raise StopIteration

        while line.rstrip() != "#=======================================":
            line = handle.readline()
            if not line:
                raise StopIteration

        length_of_seqs = None
        number_of_seqs = None
        ids = []
        header_dict = {}

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

            # Parse the rest of the header
            if key == "identity":
                header_dict["identity"] = int(parts[1].strip().split("/")[0])
            if key == "similarity":
                header_dict["similarity"] = int(parts[1].strip().split("/")[0])
            if key == "gaps":
                header_dict["gaps"] = int(parts[1].strip().split("/")[0])
            if key == "score":
                header_dict["score"] = float(parts[1].strip())

            # And read in another line...
            line = handle.readline()

        if number_of_seqs is None:
            raise ValueError("Number of sequences missing!")
        if length_of_seqs is None:
            raise ValueError("Length of sequences missing!")

        if (
            self.records_per_alignment is not None
            and self.records_per_alignment != number_of_seqs
        ):
            raise ValueError(
                "Found %i records in this alignment, told to expect %i"
                % (number_of_seqs, self.records_per_alignment)
            )

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
                    if start >= end:
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

                    if index < 0 or index >= number_of_seqs:
                        raise ValueError(
                            "Expected index %i in range [0,%i)"
                            % (index, number_of_seqs)
                        )
                    # The identifier is truncated...
                    assert id == ids[index] or id == ids[index][: len(id)]

                    if len(seq_starts) == index:
                        # Record the start
                        seq_starts.append(start)

                    # Check the start...
                    if start >= end:
                        assert seq.replace("-", "") == "", line
                    elif start - seq_starts[index] != len(seqs[index].replace("-", "")):
                        raise ValueError(
                            "Found %i chars so far for sequence %i (%s, %r), line says start %i:\n%s"
                            % (
                                len(seqs[index].replace("-", "")),
                                index,
                                id,
                                seqs[index],
                                start,
                                line,
                            )
                        )
                    seqs[index] += seq

                    # Check the end ...
                    if end != seq_starts[index] + len(seqs[index].replace("-", "")):
                        raise ValueError(
                            "Found %i chars so far for sequence %i (%s, %r, start=%i), file says end %i:\n%s"
                            % (
                                len(seqs[index].replace("-", "")),
                                index,
                                id,
                                seqs[index],
                                seq_starts[index],
                                end,
                                line,
                            )
                        )

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
                raise ValueError("Unrecognised EMBOSS pairwise line: %r\n" % line)

            line = handle.readline()
            if (
                line.rstrip() == "#---------------------------------------"
                or line.rstrip() == "#======================================="
            ):
                # End of alignment
                self._header = line
                break

        assert index == 0

        if (
            self.records_per_alignment is not None
            and self.records_per_alignment != len(ids)
        ):
            raise ValueError(
                "Found %i records in this alignment, told to expect %i"
                % (len(ids), self.records_per_alignment)
            )

        records = []
        for id, seq in zip(ids, seqs):
            if len(seq) != length_of_seqs:
                # EMBOSS 2.9.0 is known to use spaces instead of minus signs
                # for leading gaps, and thus fails to parse.  This old version
                # is still used as of Dec 2008 behind the EBI SOAP webservice:
                # http://www.ebi.ac.uk/Tools/webservices/wsdl/WSEmboss.wsdl
                raise ValueError(
                    "Error parsing alignment - sequences of "
                    "different length? You could be using an "
                    "old version of EMBOSS."
                )
            records.append(SeqRecord(Seq(seq), id=id, description=id))
        return MultipleSeqAlignment(records, annotations=header_dict)
