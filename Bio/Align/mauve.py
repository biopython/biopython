# Copyright 2015-2015 by Eric Rasche.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.AlignIO support for "xmfa" output from Mauve/ProgressiveMauve.

You are expected to use this module via the Bio.AlignIO functions (or the
Bio.SeqIO functions if you want to work directly with the gapped sequences).

For example, consider a progressiveMauve alignment file containing the following::

    #FormatVersion Mauve1
    #Sequence1File	a.fa
    #Sequence1Entry	1
    #Sequence1Format	FastA
    #Sequence2File	b.fa
    #Sequence2Entry	2
    #Sequence2Format	FastA
    #Sequence3File	c.fa
    #Sequence3Entry	3
    #Sequence3Format	FastA
    #BackboneFile	three.xmfa.bbcols
    > 1:0-0 + a.fa
    --------------------------------------------------------------------------------
    --------------------------------------------------------------------------------
    --------------------------------------------------------------------------------
    > 2:5417-5968 + b.fa
    TTTAAACATCCCTCGGCCCGTCGCCCTTTTATAATAGCAGTACGTGAGAGGAGCGCCCTAAGCTTTGGGAAATTCAAGC-
    --------------------------------------------------------------------------------
    CTGGAACGTACTTGCTGGTTTCGCTACTATTTCAAACAAGTTAGAGGCCGTTACCTCGGGCGAACGTATAAACCATTCTG
    > 3:9476-10076 - c.fa
    TTTAAACACCTTTTTGGATG--GCCCAGTTCGTTCAGTTGTG-GGGAGGAGATCGCCCCAAACGTATGGTGAGTCGGGCG
    TTTCCTATAGCTATAGGACCAATCCACTTACCATACGCCCGGCGTCGCCCAGTCCGGTTCGGTACCCTCCATGACCCACG
    ---------------------------------------------------------AAATGAGGGCCCAGGGTATGCTT
    =
    > 2:5969-6015 + b.fa
    -----------------------
    GGGCGAACGTATAAACCATTCTG
    > 3:9429-9476 - c.fa
    TTCGGTACCCTCCATGACCCACG
    AAATGAGGGCCCAGGGTATGCTT

This is a multiple sequence alignment with multiple aligned sections, so you
would probably load this using the Bio.AlignIO.parse() function:

    >>> from Bio import AlignIO
    >>> align = AlignIO.parse("Mauve/simple_short.xmfa", "mauve")
    >>> alignments = list(align)
    >>> for aln in alignments:
    ...     print(aln)
    ...
    Alignment with 3 rows and 240 columns
    --------------------------------------------...--- a.fa
    TTTAAACATCCCTCGGCCCGTCGCCCTTTTATAATAGCAGTACG...CTG b.fa/5416-5968
    TTTAAACACCTTTTTGGATG--GCCCAGTTCGTTCAGTTGTG-G...CTT c.fa/9475-10076
    Alignment with 2 rows and 46 columns
    -----------------------GGGCGAACGTATAAACCATTCTG b.fa/5968-6015
    TTCGGTACCCTCCATGACCCACGAAATGAGGGCCCAGGGTATGCTT c.fa/9428-9476

Additional information is extracted from the XMFA file and available through
the annotation attribute of each record::

    >>> for record in alignments[0]:
    ...     print(record.id, len(record))
    ...     print("  start: %d, end: %d, strand: %d" %(
    ...         record.annotations['start'], record.annotations['end'],
    ...         record.annotations['strand']))
    ...
    a.fa 240
      start: 0, end: 0, strand: 1
    b.fa/5416-5968 240
      start: 5416, end: 5968, strand: 1
    c.fa/9475-10076 240
      start: 9475, end: 10076, strand: -1

"""
import re

from Bio.Align import interfaces, Alignment
from Bio.Seq import Seq, reverse_complement
from Bio.SeqRecord import SeqRecord


ID_LINE_FMT = "> {seq_name}:{start}-{end} {strand} {filename} # {ugly_hack}"


class MauveWriter(interfaces.AlignmentWriter):
    """Mauve/XMFA alignment writer."""

    def __init__(self, *args, **kwargs):
        """Initialize the class."""
        super().__init__(*args, **kwargs)
        self._wrote_header = False
        self._wrote_first = False

    def write_alignment(self, alignment):
        """Use this to write (another) single alignment to an open file.

        Note that sequences and their annotation are recorded
        together (rather than having a block of annotation followed
        by a block of aligned sequences).
        """
        count = len(alignment)

        self._length_of_sequences = alignment.get_alignment_length()

        # NOTE - For now, the alignment object does not hold any per column
        # or per alignment annotation - only per sequence.

        if count == 0:
            raise ValueError("Must have at least one sequence")
        if self._length_of_sequences == 0:
            raise ValueError("Non-empty sequences are required")

        if not self._wrote_header:
            self._wrote_header = True
            self.handle.write("#FormatVersion Mauve1\n")
            # There are some more headers, but we ignore those for now.
            # Sequence1File	unknown.fa
            # Sequence1Entry	1
            # Sequence1Format	FastA
            for i in range(1, count + 1):
                self.handle.write("#Sequence%sEntry\t%s\n" % (i, i))

        for idx, record in enumerate(alignment):
            self._write_record(record, record_idx=idx)
        self.handle.write("=\n")

    def _write_record(self, record, record_idx=0):
        """Write a single SeqRecord to the file (PRIVATE)."""
        if self._length_of_sequences != len(record.seq):
            raise ValueError("Sequences must all be the same length")

        seq_name = record.name
        try:
            seq_name = str(int(record.name))
        except ValueError:
            seq_name = str(record_idx + 1)

        # We remove the "/{start}-{end}" before writing, as it cannot be part
        # of the produced XMFA file.
        if "start" in record.annotations and "end" in record.annotations:
            suffix0 = "/%s-%s" % (
                record.annotations["start"],
                record.annotations["end"],
            )
            suffix1 = "/%s-%s" % (
                record.annotations["start"] + 1,
                record.annotations["end"],
            )
            if seq_name[-len(suffix0) :] == suffix0:
                seq_name = seq_name[: -len(suffix0)]
            if seq_name[-len(suffix1) :] == suffix1:
                seq_name = seq_name[: -len(suffix1)]

        if (
            "start" in record.annotations
            and "end" in record.annotations
            and "strand" in record.annotations
        ):
            id_line = ID_LINE_FMT.format(
                seq_name=seq_name,
                start=record.annotations["start"] + 1,
                end=record.annotations["end"],
                strand=("+" if record.annotations["strand"] == 1 else "-"),
                filename=record.name + ".fa",
                ugly_hack=record.id,
            )
            lacking_annotations = False
        else:
            id_line = ID_LINE_FMT.format(
                seq_name=seq_name,
                start=0,
                end=0,
                strand="+",
                filename=record.name + ".fa",
                ugly_hack=record.id,
            )
            lacking_annotations = True

        # If the sequence is an empty one, skip writing it out
        if (":0-0 " in id_line or ":1-0 " in id_line) and not lacking_annotations:
            # Except in the first LCB
            if not self._wrote_first:
                self._wrote_first = True
                # The first LCB we write out is special, and must list ALL
                # sequences, for the Mauve GUI
                # http://darlinglab.org/mauve/user-guide/files.html#non-standard-xmfa-formatting-used-by-the-mauve-gui
                id_line = ID_LINE_FMT.format(
                    seq_name=seq_name,
                    start=0,
                    end=0,
                    strand="+",
                    filename=record.name + ".fa",
                    ugly_hack=record.id,
                )
                id_line = id_line.replace("\n", " ").replace("\r", " ")
                self.handle.write(id_line + "\n\n")
            # Alignments lacking a start/stop/strand were generated by
            # Biopython on load, and shouldn't exist according to XMFA
        else:
            # In other blocks, we only write sequences if they exist in a given
            # alignment.
            id_line = id_line.replace("\n", " ").replace("\r", " ")
            self.handle.write(id_line + "\n")
            for i in range(0, len(record.seq), 80):
                self.handle.write("%s\n" % record.seq[i : i + 80])


class AlignmentIterator(interfaces.AlignmentIterator):
    """Mauve xmfa alignment iterator."""

    _ids = []  # for caching IDs between __next__ calls

    def __init__(self, source):
        """Create an AlignmentIterator object.

        Arguments:
         - source   - input data or file name

        """
        super().__init__(source, mode="t", fmt="Mauve")
        stream = self.stream
        metadata = {}
        prefix = "Sequence"
        suffixes = ("File", "Entry", "Format")
        id_info = {}
        for suffix in suffixes:
            id_info[suffix] = []
        for line in stream:
            if not line.startswith("#"):
                self._line = line.strip()
                break
            key, value = line[1:].split()
            if key.startswith(prefix):
                for suffix in suffixes:
                    if key.endswith(suffix):
                        break
                else:
                    raise ValueError("Unexpected keyword '%s'" % key)
                seq_num = int(key[len(prefix):-len(suffix)])
                id_info[suffix].append(value)
                assert seq_num == len(id_info[suffix])  # Mauve uses 1-based counting
            else:
                metadata[key] = value.strip()
        else:
            if not metadata:
                raise ValueError("Empty file.") from None
        if len(set(id_info["File"])) == 1:
            self._identifiers = []
            # Use filename + entry number as ID
            for filename, entry in zip(id_info["File"], id_info["Entry"]):
                entry = int(entry) - 1  # Use 0-based counting
                identifier = "%s:%d" % (filename, entry)
                self._identifiers.append(identifier)
        else:
            assert len(set(id_info["File"])) == len(id_info["File"])
            # All filenames are unique
            self._identifiers = id_info["File"]
        self.metadata = metadata

    def _parse_description(self, line):
        assert line.startswith(">")
        locus, strand, comments = line[1:].split(None, 2)
        seq_num, start_end = locus.split(":")
        seq_num = int(seq_num) - 1  # python counting
        identifier = self._identifiers[seq_num]
        assert strand in "+-"
        start, end = start_end.split("-")
        start = int(start)
        end = int(end)
        if start == 0:
            assert end == 0  # unaligned sequence
        else:
            start -= 1  # python counting
        return (identifier, start, end, strand, comments)

    def parse(self, stream):
        """Parse the next alignment from the stream."""
        if stream is None:
            return

        descriptions = []
        seqs = []

        line = self._line
        del self._line
        description = self._parse_description(line)
        identifier, start, end, strand, comments = description
        if end > 0:
            descriptions.append(description)
            seqs.append("")

        for line in stream:
            line = line.strip()
            if line.startswith("="):
                # There may be more data, but we've reached the end of this
                # alignment
                coordinates = Alignment.infer_coordinates(seqs)
                records = []
                for index, (description, seq) in enumerate(zip(descriptions, seqs)):
                    identifier, start, end, strand, comments = description
                    length = end - start
                    seq = seq.replace("-", "")
                    assert len(seq) == end - start
                    if strand == "+":
                        pass
                    elif strand == "-":
                        seq = reverse_complement(seq, inplace=False)
                        coordinates[index, :] = len(seq) - coordinates[index, :]
                    else:
                        raise ValueError("Unexpected strand '%s'" % strand)
                    coordinates[index] += start
                    if start == 0:
                        seq = Seq(seq)
                    else:
                        seq = Seq({start: seq}, length=end)
                    record = SeqRecord(seq, id=identifier, description=comments)
                    records.append(record)

                yield Alignment(records, coordinates)

                descriptions = []
                seqs = []
            elif line.startswith(">"):
                description = self._parse_description(line)
                identifier, start, end, strand, comments = description
                if end > 0:
                    descriptions.append(description)
                    seqs.append("")
            elif end > 0:
                seqs[-1] += line
