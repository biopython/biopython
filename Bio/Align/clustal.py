# Copyright 2006-2016 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.Align support for "clustal" output from CLUSTAL W and other tools.

You are expected to use this module via the Bio.Align functions (or the
Bio.SeqIO functions if you are interested in the sequences only).
"""
from Bio.Align import Alignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class ClustalWriter:
    """Clustalw alignment writer."""

    def write_alignment(self, alignment):
        """Use this to write (another) single alignment to an open file."""
        if len(alignment) == 0:
            raise ValueError("Must have at least one sequence")
        if alignment.get_alignment_length() == 0:
            # This doubles as a check for an alignment object
            raise ValueError("Non-empty sequences are required")

        try:
            version = str(alignment._version)
        except AttributeError:
            version = ""
        if not version:
            version = "1.81"
        if version.startswith("2."):
            # e.g. 2.0.x
            output = "CLUSTAL %s multiple sequence alignment\n\n\n" % version
        else:
            # e.g. 1.81 or 1.83
            output = "CLUSTAL X (%s) multiple sequence alignment\n\n\n" % version

        cur_char = 0
        max_length = len(alignment[0])

        if max_length <= 0:
            raise ValueError("Non-empty sequences are required")

        try:
            column_annotations = alignment.column_annotations
        except AttributeError:
            consensus = None
        else:
            consensus = column_annotations.get("clustal_consensus")

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
            if consensus is not None:
                output += (
                    (" " * 36) + consensus[cur_char : (cur_char + show_num)] + "\n"
                )

            output += "\n"
            cur_char += show_num

        # Want a trailing blank new line in case the output is concatenated
        self.stream.write(output + "\n")


class Iterator:
    """Clustalw alignment iterator."""

    _header = None  # for caching lines between __next__ calls

    def __init__(self, stream):
        """Create an Iterator object.

        Arguments:
         - stream   - input data or file name

        """
        self.stream = stream
        try:
            line = next(stream)
        except StopIteration:
            raise ValueError("Empty file.") from None

        # Whitelisted programs we know about
        words = line.split()
        known_programs = ["CLUSTAL", "PROBCONS", "MUSCLE", "MSAPROBS", "Kalign"]
        program = words[0]
        if program not in known_programs:
            raise ValueError(
                "%s is not known to generate CLUSTAL files: %s"
                % (program, ", ".join(known_programs))
            )
        self.program = program

        # find the clustal version in the header line
        for word in words:
            if word[0] == "(" and word[-1] == ")":
                word = word[1:-1]
            if word[0] in "0123456789":
                self.version = word
                break
        else:
            self.version = None

    def __iter__(self):
        return self

    def __next__(self):
        """Parse the next alignment from the stream."""
        stream = self.stream
        if stream is None:
            raise StopIteration

        # If the alignment contains entries with the same sequence
        # identifier (not a good idea - but seems possible), then this
        # dictionary based parser will merge their sequences.  Fix this?
        ids = []
        seqs = []
        aligned_seqs = []
        consensus = ""
        index = None  # Used to extract the consensus

        # Use the first block to get the sequence identifiers
        for line in stream:
            if line.startswith(" "):
                # Sequence consensus line...
                assert len(ids) > 0
                assert index is not None
                length = len(aligned_seq)
                consensus = line[index:index+length]
                break
            elif line.strip():
                # Sequences identifier...
                fields = line.split()

                # We expect there to be two fields, there can be an optional
                # "sequence number" field containing the letter count.
                if len(fields) < 2 or len(fields) > 3:
                    raise ValueError("Could not parse line:\n%s" % line)

                seqid, aligned_seq = fields[:2]
                ids.append(seqid)
                aligned_seqs.append(aligned_seq)
                seq = aligned_seq.replace("-", "")
                seqs.append(seq)

                # Record the sequence position to get the consensus
                if index is None:
                    index = line.find(aligned_seq, len(seqid))

                if len(fields) == 3:
                    # This MAY be an old style file with a letter count...
                    try:
                        letters = int(fields[2])
                    except ValueError:
                        raise ValueError(
                            "Could not parse line, bad sequence number:\n%s" % line
                        ) from None
                    if len(seq) != letters:
                        raise ValueError(
                            "Could not parse line, invalid sequence number:\n%s" % line
                        )
            else:
                # no consensus line
                if index:
                    break
        else:
            raise ValueError("Failed to find end of first block")

        assert index is not None

        # Confirm all same length
        length = len(aligned_seqs[0])
        for aligned_seq in aligned_seqs:
            assert len(aligned_seq) == length
        if consensus:
            assert len(consensus) == length

        n = len(seqs)
        i = 0
        # Loop over any remaining blocks...
        for line in stream:
            if line.startswith(" "):  # Sequence consensus line
                assert index is not None
                length = len(aligned_seq)
                consensus += line[index:index+length]
            elif not line.strip():  # Blank line
                continue
            else:
                seqid = ids[i]
                # Sequences identifier...
                fields = line.split()

                # We expect there to be two fields, there can be an optional
                # "sequence number" field containing the letter count.
                if len(fields) < 2 or len(fields) > 3:
                    raise ValueError("Could not parse line:\n%s" % line)

                assert seqid == fields[0]
                aligned_seq = fields[1]
                aligned_seqs[i] += aligned_seq
                seq = aligned_seq.replace("-", "")
                seqs[i] += seq

                if len(fields) == 3:
                    # This MAY be an old style file with a letter count...
                    try:
                        letters = int(fields[2])
                    except ValueError:
                        raise ValueError(
                            "Could not parse line, bad sequence number:\n%s" % line
                        ) from None
                    if len(seqs[i]) != letters:
                        raise ValueError(
                            "Could not parse line, invalid sequence number:\n%s" % line
                        )
                i += 1
                if i == n:
                    i = 0

        records = [SeqRecord(Seq(seq), id=seqid, description=seqid) for (seqid, seq) in zip(ids, seqs)]
        alignment = Alignment(records, None)
        coordinates = alignment._infer_coordinates(aligned_seqs)
        alignment.coordinates = coordinates
        # TODO - Handle alignment annotation better, for now
        # mimic the old parser in Bio.Clustalw
        if consensus:
            rows, columns = alignment.shape
            if len(consensus) != columns:
                for aligned_seq in aligned_seqs:
                    print(aligned_seq, len(aligned_seq))
                raise ValueError(
                    "Alignment has %i columns, consensus length is %i, '%s'"
                    % (columns, len(consensus), consensus)
                )
            alignment.column_annotations = {}
            alignment.column_annotations["clustal_consensus"] = consensus
        self.stream = None
        return alignment
