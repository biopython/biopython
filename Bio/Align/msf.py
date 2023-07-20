# Copyright 2019, National Marrow Donor Program (NMPD).  All rights reserved.
# Written by Peter Cock, The James Hutton Institute, under contract to NMDP.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.Align support for GCG MSF format.

The file format was produced by the GCG PileUp and LocalPileUp tools, and later
tools such as T-COFFEE and MUSCLE support it as an optional output format.

You are expected to use this module via the Bio.Align functions.
"""
import warnings

from Bio.Align import Alignment
from Bio.Align import interfaces
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from Bio import BiopythonParserWarning


class AlignmentIterator(interfaces.AlignmentIterator):
    """GCG MSF alignment iterator."""

    fmt = "MSF"

    def _read_next_alignment(self, stream):
        try:
            line = next(stream)
        except StopIteration:
            if stream.tell() == 0:
                raise ValueError("Empty file.") from None
            return
        # Whitelisted headers we know about.
        known_headers = ["!!NA_MULTIPLE_ALIGNMENT", "!!AA_MULTIPLE_ALIGNMENT", "PileUp"]
        # Examples in "Molecular Biology Software Training Manual GCG version 10"
        # by BBSRC Bioscuences IT Services (BITS), Harpenden, UK, Copyright 1996-2001
        # would often start as follows:
        #
        # !!AA_MUTIPLE_ALIGNMENT 1.0
        # PileUp of: @/usr/users2/culhane/...
        #
        # etc with other seemingly free format text before getting to the
        # MSF/Type/Check line and the following Name: lines block and // line.
        #
        # MUSCLE just has a line "PileUp", while other sources just use the line
        # "!!AA_MULTIPLE_ALIGNMENT" (amino acid) or "!!NA_MULTIPLE_ALIGNMENT"
        # (nucleotide).
        if line.strip().split()[0] not in known_headers:
            raise ValueError(
                "%s is not a known GCG MSF header: %s"
                % (line.strip().split()[0], ", ".join(known_headers))
            )

        for line in stream:
            line = line.rstrip("\n")
            if "MSF: " in line and line.endswith(".."):
                break
        else:
            raise ValueError("Reached end of file without MSF/Type/Check header line")

        # Quoting from "Molecular Biology Software Training Manual GCG version 10"
        # by BBSRC Bioscuences IT Services (BITS), Harpenden, UK. Copyright 1996-2001.
        # Page 31:
        #
        # "Header information is before a .. (double dot) in a GCG format file.
        #  The file will also have a checksum specific for that file."
        #
        # This was followed by a single non-aligned sequence, but this convention
        # appears to also be used in the GCG MSF files. Quoting other examples in
        # this reference, page 31:
        #
        # localpileup_17.msf  MSF: 195  Type: P  January 6, 2000 15:41  Check: 4365 ..
        #
        # Except from page 148:
        #
        # localpileup_106.msf  MSF: 457  Type: P  November 28, 2000 16:09  Check: 2396 ..
        #
        # Quoting output from MUSCLE v3.8, have two leading spaces and a zero checksum:
        #
        #   MSF: 689  Type: N  Check: 0000  ..
        #
        # By observation, the MSF value is the column count, type is N (nucleotide)
        # or P (protein / amino acid).
        #
        # In a possible bug, EMBOSS v6.6.0.0 uses CompCheck: rather than Check: as shown,
        #
        # $ seqret -sequence Tests/Fasta/f002 -auto -stdout -osformat msf
        # !!NA_MULTIPLE_ALIGNMENT 1.0
        #
        #   stdout MSF: 633 Type: N 01/08/19 CompCheck: 8543 ..
        #
        #   Name: G26680     Len: 633  Check: 4334 Weight: 1.00
        #   Name: G26685     Len: 633  Check: 3818 Weight: 1.00
        #   Name: G29385     Len: 633  Check:  391 Weight: 1.00
        #
        # //
        #
        parts = line.split()
        offset = parts.index("MSF:")
        if parts[offset + 2] != "Type:" or parts[-3] not in ("Check:", "CompCheck:"):
            raise ValueError(
                "GCG MSF header line should be "
                "'<optional text> MSF: <int> Type: <letter> <optional date> Check: <int> ..', "
                " not: %r" % line
            )
        try:
            aln_length = int(parts[offset + 1])
        except ValueError:
            raise ValueError(
                "GCG MSF header line should have MSF: <int> for column count, not %r"
                % parts[offset + 1]
            ) from None
        seq_type = parts[offset + 3]
        if seq_type not in ["P", "N"]:
            raise ValueError(
                "GCG MSF header line should have 'Type: P' (protein) "
                "or 'Type: N' (nucleotide), not 'Type: %s'" % seq_type
            )

        # There should be a blank line after that header line, then the Name: lines
        #
        # The Name may be followed by 'oo', as shown here:
        #
        # PileUp
        #
        #
        #
        #    MSF:  628  Type: P    Check:   147   ..
        #
        #  Name: AK1H_ECOLI/1-378 oo  Len:  628  Check:  3643  Weight:  1.000
        #  Name: AKH_HAEIN/1-382 oo  Len:  628  Check:  6504  Weight:  1.000
        #
        # //
        names = []
        remaining = []
        checks = []
        weights = []
        for line in stream:
            line = line.strip()
            if line == "//":
                break
            if line.startswith("Name: "):
                words = line.split()
                try:
                    index_name = words.index("Name:")
                    index_len = words.index("Len:")
                    index_weight = words.index("Weight:")
                    index_check = words.index("Check:")
                except ValueError:
                    raise ValueError(f"Malformed GCG MSF name line: {line!r}") from None
                name = words[index_name + 1]
                length = int(words[index_len + 1])
                weight = float(words[index_weight + 1])
                check = words[index_check + 1]
                if name in names:
                    raise ValueError(f"Duplicated ID of {name!r}")
                names.append(name)
                remaining.append(length)
                checks.append(check)
                weights.append(weight)
        else:
            raise ValueError("End of file while looking for end of header // line.")

        try:
            line = next(stream)
        except StopIteration:
            raise ValueError("End of file after // line, expected sequences.") from None
        if line.strip():
            raise ValueError("After // line, expected blank line before sequences.")

        # Now load the sequences
        seqs = [""] * len(names)
        for line in stream:
            words = line.split()
            if not words:
                continue
            name = words[0]
            try:
                index = names.index(name)
            except ValueError:
                # This may be a coordinate line
                for word in words:
                    if not word.isdigit():
                        break
                else:
                    # all words are integers; assume this is a coordinate line
                    continue
                raise ValueError(f"Unexpected line '{line}' in input") from None
            seq = "".join(words[1:])
            length = remaining[index] - (len(seq) - seq.count("-"))
            if length < 0:
                raise ValueError("Received longer sequence than expected for %s" % name)
            seqs[index] += seq
            remaining[index] = length
            if all(length == 0 for length in remaining):
                break
        else:
            raise ValueError("End of file where expecting sequence data.")

        # skip any remaining empty lines
        for line in stream:
            assert line.strip() == ""

        length = max(len(seq) for seq in seqs)
        if length != aln_length:
            warnings.warn(
                "GCG MSF headers said alignment length %i, but found %i"
                % (aln_length, length),
                BiopythonParserWarning,
                stacklevel=2,
            )
            aln_length = length

        # Combine list of strings into single string, remap gaps
        for index, seq in enumerate(seqs):
            seq = "".join(seq).replace("~", "-").replace(".", "-")
            if len(seq) < aln_length:
                seq += "-" * (aln_length - len(seq))
            seqs[index] = seq

        coordinates = Alignment.infer_coordinates(seqs)
        seqs = (Seq(seq.replace("-", "")) for seq in seqs)
        records = [
            SeqRecord(
                seq,
                id=name,
                name=name,
                description=name,
                annotations={"weight": weight},
            )
            for (name, seq, weight) in zip(names, seqs, weights)
        ]

        alignment = Alignment(records, coordinates)
        # This will check alignment lengths are self-consistent:
        columns = alignment.length
        if columns != aln_length:
            raise ValueError(
                "GCG MSF headers said alignment length %i, but found %i"
                % (aln_length, columns)
            )
        return alignment
