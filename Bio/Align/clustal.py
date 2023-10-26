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
import Bio
from Bio.Align import Alignment
from Bio.Align import interfaces
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class AlignmentWriter(interfaces.AlignmentWriter):
    """Clustalw alignment writer."""

    fmt = "Clustal"

    def write_header(self, stream, alignments):
        """Use this to write the file header."""
        try:
            metadata = alignments.metadata
            program = metadata["Program"]
        except (AttributeError, KeyError):
            program = "Biopython"
            version = Bio.__version__
        else:
            version = metadata.get("Version", "")
        line = f"{program} {version} multiple sequence alignment\n"
        stream.write(line)
        stream.write("\n")
        stream.write("\n")

    def format_alignment(self, alignment):
        """Return a string with a single alignment in the Clustal format."""
        nseqs, length = alignment.shape
        if nseqs == 0:
            raise ValueError("Must have at least one sequence")
        if length == 0:
            raise ValueError("Non-empty sequences are required")

        try:
            column_annotations = alignment.column_annotations
        except AttributeError:
            consensus = None
        else:
            consensus = column_annotations.get("clustal_consensus")

        gapped_sequences = list(alignment)
        names = []
        for i, sequence in enumerate(alignment.sequences):
            try:
                name = sequence.id
            except AttributeError:
                name = "sequence_%d" % i  # Clustal format doesn't allow an empty string
            else:
                # when we output, we do a nice 80 column output, although
                # this may result in truncation of the ids.  Also, make sure
                # we don't get any spaces in the record identifier when output
                # in the file by replacing them with underscores.
                name = name[:30].replace(" ", "_")
            name = name.ljust(36)
            names.append(name)

        lines = []
        start = 0
        while start != length:
            # calculate the number of letters to show, which will
            # be less if we are at the end of the alignment.
            stop = start + 50
            if stop > length:
                stop = length

            for name, gapped_sequence in zip(names, gapped_sequences):
                line = f"{name}{gapped_sequence[start:stop]}\n"
                lines.append(line)

            # now we need to print out the star info, if we've got it
            if consensus is not None:
                line = " " * 36 + consensus[start:stop] + "\n"
                lines.append(line)

            lines.append("\n")
            start = stop
        lines.append("\n")
        return "".join(lines)


class AlignmentIterator(interfaces.AlignmentIterator):
    """Clustalw alignment iterator."""

    fmt = "Clustal"

    def _read_header(self, stream):
        try:
            line = next(stream)
        except StopIteration:
            raise ValueError("Empty file.") from None

        self.metadata = {}
        # Whitelisted programs we know about
        words = line.split()
        known_programs = [
            "CLUSTAL",
            "PROBCONS",
            "MUSCLE",
            "MSAPROBS",
            "Kalign",
            "Biopython",
        ]
        program = words[0]
        if program not in known_programs:
            raise ValueError(
                "%s is not known to generate CLUSTAL files: %s"
                % (program, ", ".join(known_programs))
            )
        self.metadata["Program"] = program

        # find the clustal version in the header line
        for word in words:
            if word[0] == "(" and word[-1] == ")":
                word = word[1:-1]
            if word[0].isdigit():
                self.metadata["Version"] = word
                break

    def _read_next_alignment(self, stream):
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
                length = len(aligned_seq)  # noqa: F821
                consensus = line[index : index + length]
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
            return

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
                consensus += line[index : index + length]
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

        records = [
            SeqRecord(Seq(seq), id=seqid, description="")
            for (seqid, seq) in zip(ids, seqs)
        ]
        coordinates = Alignment.infer_coordinates(aligned_seqs)
        alignment = Alignment(records, coordinates)
        if consensus:
            columns = alignment.length
            if len(consensus) != columns:
                raise ValueError(
                    "Alignment has %i columns, consensus length is %i, '%s'"
                    % (columns, len(consensus), consensus)
                )
            alignment.column_annotations = {}
            alignment.column_annotations["clustal_consensus"] = consensus
        return alignment
