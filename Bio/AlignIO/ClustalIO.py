# Copyright 2006-2016 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.AlignIO support for "clustal" output from CLUSTAL W and other tools.

You are expected to use this module via the Bio.AlignIO functions (or the
Bio.SeqIO functions if you want to work directly with the gapped sequences).
"""
from Bio.Align import MultipleSeqAlignment
from Bio.AlignIO.Interfaces import AlignmentIterator
from Bio.AlignIO.Interfaces import SequentialAlignmentWriter
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class ClustalWriter(SequentialAlignmentWriter):
    """Clustalw alignment writer."""

    def write_alignment(self, alignment):
        """Use this to write (another) single alignment to an open file."""
        if len(alignment) == 0:
            raise ValueError("Must have at least one sequence")
        if alignment.get_alignment_length() == 0:
            # This doubles as a check for an alignment object
            raise ValueError("Non-empty sequences are required")

        # Old versions of the parser in Bio.Clustalw used a ._version property
        try:
            version = str(alignment._version)
        except AttributeError:
            version = ""
        if not version:
            version = "1.81"
        if version.startswith("2."):
            # e.g. 2.0.x
            output = f"CLUSTAL {version} multiple sequence alignment\n\n\n"
        else:
            # e.g. 1.81 or 1.83
            output = f"CLUSTAL X ({version}) multiple sequence alignment\n\n\n"

        cur_char = 0
        max_length = len(alignment[0])

        if max_length <= 0:
            raise ValueError("Non-empty sequences are required")

        if "clustal_consensus" in alignment.column_annotations:
            star_info = alignment.column_annotations["clustal_consensus"]
        else:
            try:
                # This was originally stored by Bio.Clustalw as ._star_info
                star_info = alignment._star_info
            except AttributeError:
                star_info = None

        # keep displaying sequences until we reach the end
        while cur_char != max_length:
            # calculate the number of sequences to show, which will
            # be less if we are at the end of the sequence
            if (cur_char + 50) > max_length:
                show_num = max_length - cur_char
            else:
                show_num = 50

            # go through all of the records and print out the sequences
            # when we output, we do a nice 80 column output, although this
            # may result in truncation of the ids.
            for record in alignment:
                # Make sure we don't get any spaces in the record
                # identifier when output in the file by replacing
                # them with underscores:
                line = record.id[0:30].replace(" ", "_").ljust(36)
                line += str(record.seq[cur_char : (cur_char + show_num)])
                output += line + "\n"

            # now we need to print out the star info, if we've got it
            if star_info:
                output += (
                    (" " * 36) + star_info[cur_char : (cur_char + show_num)] + "\n"
                )

            output += "\n"
            cur_char += show_num

        # Want a trailing blank new line in case the output is concatenated
        self.handle.write(output + "\n")


class ClustalIterator(AlignmentIterator):
    """Clustalw alignment iterator."""

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

        # Whitelisted headers we know about
        known_headers = [
            "CLUSTAL",
            "PROBCONS",
            "MUSCLE",
            "MSAPROBS",
            "Kalign",
            "Biopython",
        ]
        if line.strip().split()[0] not in known_headers:
            raise ValueError(
                "%s is not a known CLUSTAL header: %s"
                % (line.strip().split()[0], ", ".join(known_headers))
            )

        # find the clustal version in the header line
        version = None
        for word in line.split():
            if word[0] == "(" and word[-1] == ")":
                word = word[1:-1]
            if word[0] in "0123456789":
                version = word
                break

        # There should be two blank lines after the header line
        line = handle.readline()
        while line.strip() == "":
            line = handle.readline()

        # If the alignment contains entries with the same sequence
        # identifier (not a good idea - but seems possible), then this
        # dictionary based parser will merge their sequences.  Fix this?
        ids = []
        seqs = []
        consensus = ""
        seq_cols = None  # Used to extract the consensus

        # Use the first block to get the sequence identifiers
        while True:
            if line[0] != " " and line.strip() != "":
                # Sequences identifier...
                fields = line.rstrip().split()

                # We expect there to be two fields, there can be an optional
                # "sequence number" field containing the letter count.
                if len(fields) < 2 or len(fields) > 3:
                    raise ValueError(f"Could not parse line:\n{line}")

                ids.append(fields[0])
                seqs.append(fields[1])

                # Record the sequence position to get the consensus
                if seq_cols is None:
                    start = len(fields[0]) + line[len(fields[0]) :].find(fields[1])
                    end = start + len(fields[1])
                    seq_cols = slice(start, end)
                    del start, end
                assert fields[1] == line[seq_cols]

                if len(fields) == 3:
                    # This MAY be an old style file with a letter count...
                    try:
                        letters = int(fields[2])
                    except ValueError:
                        raise ValueError(
                            f"Could not parse line, bad sequence number:\n{line}"
                        ) from None
                    if len(fields[1].replace("-", "")) != letters:
                        raise ValueError(
                            f"Could not parse line, invalid sequence number:\n{line}"
                        )
            elif line[0] == " ":
                # Sequence consensus line...
                assert len(ids) == len(seqs)
                assert len(ids) > 0
                assert seq_cols is not None
                consensus = line[seq_cols]
                assert not line[: seq_cols.start].strip()
                assert not line[seq_cols.stop :].strip()
                # Check for blank line (or end of file)
                line = handle.readline()
                assert line.strip() == ""
                break
            else:
                # No consensus
                break
            line = handle.readline()
            if not line:
                break  # end of file

        assert line.strip() == ""
        assert seq_cols is not None

        # Confirm all same length
        for s in seqs:
            assert len(s) == len(seqs[0])
        if consensus:
            assert len(consensus) == len(seqs[0])

        # Loop over any remaining blocks...
        done = False
        while not done:
            # There should be a blank line between each block.
            # Also want to ignore any consensus line from the
            # previous block.
            while (not line) or line.strip() == "":
                line = handle.readline()
                if not line:
                    break  # end of file
            if not line:
                break  # end of file

            if line.split(None, 1)[0] in known_headers:
                # Found concatenated alignment.
                self._header = line
                break

            for i in range(len(ids)):
                if line[0] == " ":
                    raise ValueError(f"Unexpected line:\n{line!r}")
                fields = line.rstrip().split()

                # We expect there to be two fields, there can be an optional
                # "sequence number" field containing the letter count.
                if len(fields) < 2 or len(fields) > 3:
                    raise ValueError(f"Could not parse line:\n{line!r}")

                if fields[0] != ids[i]:
                    raise ValueError(
                        "Identifiers out of order? Got '%s' but expected '%s'"
                        % (fields[0], ids[i])
                    )

                if fields[1] != line[seq_cols]:
                    start = len(fields[0]) + line[len(fields[0]) :].find(fields[1])
                    if start != seq_cols.start:
                        raise ValueError("Old location %s -> %i:XX" % (seq_cols, start))
                    end = start + len(fields[1])
                    seq_cols = slice(start, end)
                    del start, end

                # Append the sequence
                seqs[i] += fields[1]
                assert len(seqs[i]) == len(seqs[0])

                if len(fields) == 3:
                    # This MAY be an old style file with a letter count...
                    try:
                        letters = int(fields[2])
                    except ValueError:
                        raise ValueError(
                            f"Could not parse line, bad sequence number:\n{line}"
                        ) from None
                    if len(seqs[i].replace("-", "")) != letters:
                        raise ValueError(
                            f"Could not parse line, invalid sequence number:\n{line}"
                        )

                # Read in the next line
                line = handle.readline()
            # There should now be a consensus line
            if consensus:
                assert line[0] == " "
                assert seq_cols is not None
                consensus += line[seq_cols]
                assert len(consensus) == len(seqs[0])
                assert not line[: seq_cols.start].strip()
                assert not line[seq_cols.stop :].strip()
                # Read in the next line
                line = handle.readline()

        assert len(ids) == len(seqs)
        if len(seqs) == 0 or len(seqs[0]) == 0:
            raise StopIteration

        if (
            self.records_per_alignment is not None
            and self.records_per_alignment != len(ids)
        ):
            raise ValueError(
                "Found %i records in this alignment, told to expect %i"
                % (len(ids), self.records_per_alignment)
            )

        records = (SeqRecord(Seq(s), id=i, description=i) for (i, s) in zip(ids, seqs))
        alignment = MultipleSeqAlignment(records)
        # TODO - Handle alignment annotation better, for now
        # mimic the old parser in Bio.Clustalw
        if version:
            alignment._version = version
        if consensus:
            alignment_length = len(seqs[0])
            if len(consensus) != alignment_length:
                raise ValueError(
                    "Alignment length is %i, consensus length is %i, '%s'"
                    % (alignment_length, len(consensus), consensus)
                )
            alignment.column_annotations["clustal_consensus"] = consensus
            # For backward compatibility prior to .column_annotations:
            alignment._star_info = consensus
        return alignment
