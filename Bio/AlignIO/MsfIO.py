# Copyright 2019, National Marrow Donor Program (NMPD).  All rights reserved.
# Written by Peter Cock, The James Hutton Institute, under contract to NMDP.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.AlignIO support for GCG MSF format.

The file format was produced by the GCG PileUp and and LocalPileUp tools,
and later tools such as T-COFFEE and MUSCLE support it as an optional
output format.

The original GCG tool would write gaps at ends of each sequence which could
be missing data as tildes (``~``), whereas internal gaps were periods (``.``)
instead. This parser replaces both with minus signs (``-``) for consistency
with the rest of ``Bio.AlignIO``.

You are expected to use this module via the Bio.AlignIO functions (or the
Bio.SeqIO functions if you want to work directly with the gapped sequences).
"""
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .Interfaces import AlignmentIterator


class MsfIterator(AlignmentIterator):
    """GCG MSF alignment iterator."""

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

        while line and " MSF: " not in line:
            line = handle.readline()

        if not line:
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
        parts = line.strip("\n").split()
        offset = parts.index("MSF:")
        if (
            parts[offset + 2] != "Type:"
            or parts[-3] not in ("Check:", "CompCheck:")
            or parts[-1] != ".."
        ):
            raise ValueError(
                "GCG MSF header line should be "
                "'<optional text> MSF: <int> Type: <letter> <optional date> Check: <int> ..', "
                " not: %r" % line
            )
        try:
            aln_length = int(parts[offset + 1])
        except ValueError:
            aln_length = -1
        if aln_length < 0:
            raise ValueError(
                "GCG MSF header line should have MDF: <int> for column count, not %r"
                % parts[offset + 1]
            )
        seq_type = parts[offset + 3]
        if seq_type not in ["P", "N"]:
            raise ValueError(
                "GCG MSF header line should have 'Type: P' (protein) "
                "or 'Type: N' (nucleotide), not 'Type: %s'" % seq_type
            )

        # There should be a blank line after that header line, then the Name: lines
        #
        # In a possible bug, T-COFFEE v12.00 adds 'oo' after the names, as shown here,
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
        ids = []
        lengths = []
        checks = []
        weights = []
        line = handle.readline()
        while line and line.strip() != "//":
            line = handle.readline()
            if line.strip().startswith("Name: "):
                if " Len: " in line and " Check: " in line and " Weight: " in line:
                    rest = line[line.index("Name: ") + 6 :].strip()
                    name, rest = rest.split(" Len: ")
                    length, rest = rest.split(" Check: ")
                    check, weight = rest.split(" Weight: ")
                    name = name.strip()
                    if name.endswith(" oo"):
                        # T-COFFEE oddity, ignore this
                        name = name[:-3]
                    if name in ids:
                        raise ValueError("Duplicated ID of %r" % name)
                    if " " in name:
                        raise NotImplementedError("Space in ID %r" % name)
                    ids.append(name)
                    # Expect aln_length <= int(length.strip()), see below
                    lengths.append(int(length.strip()))
                    checks.append(int(check.strip()))
                    weights.append(float(weight.strip()))
                else:
                    raise ValueError("Malformed GCG MSF name line: %r" % line)
        if not line:
            raise ValueError("End of file while looking for end of header // line.")

        if aln_length != max(lengths):
            # In broken examples from IMGTHLA was possible to continue
            # https://github.com/ANHIG/IMGTHLA/issues/201
            max_length = max(lengths)
            max_count = sum(1 for _ in lengths if _ == max_length)
            raise ValueError(
                "GCG MSF header said alignment length %i, but %s of %i sequences said Len: %s"
                % (aln_length, max_count, len(ids), max_length)
            )

        line = handle.readline()
        if not line:
            raise ValueError("End of file after // line, expected sequences.")
        if line.strip():
            raise ValueError("After // line, expected blank line before sequences.")

        # Now load the sequences
        seqs = [[] for _ in ids]  # list of empty lists
        completed_length = 0
        while completed_length < aln_length:
            # Note might have a coordinate header line (seems to be optional)
            for idx, name in enumerate(ids):
                line = handle.readline()
                if idx == 0 and not line.strip():
                    # T-COFFEE uses two blank lines between blocks, rather than one
                    while line and not line.strip():
                        line = handle.readline()
                if not line:
                    raise ValueError("End of file where expecting sequence data.")
                # print("Looking for seq for %s in line: %r" % (name, line))
                words = line.strip().split()
                # Should we use column numbers, rather than assuming no spaces in names?
                if idx == 0 and words and words[0] != name:
                    # print("Actually have a coord line")
                    # Hopefully this is a coordinate header before the first seq
                    try:
                        i = int(words[0])
                    except ValueError:
                        i = -1
                    if i != completed_length + 1:
                        raise ValueError(
                            "Expected GCG MSF coordinate line starting %i, got: %r"
                            % (completed_length + 1, line)
                        )
                    if len(words) > 1:
                        # Final block usually not full 50 chars, so expect start only.
                        if len(words) != 2:
                            i = -1
                        else:
                            try:
                                i = int(words[1])
                            except ValueError:
                                i = -1
                        if i != (
                            completed_length + 50
                            if completed_length + 50 < aln_length
                            else aln_length
                        ):
                            raise ValueError(
                                "Expected GCG MSF coordinate line %i to %i, got: %r"
                                % (
                                    completed_length + 1,
                                    completed_length + 50
                                    if completed_length + 50 < aln_length
                                    else aln_length,
                                    line,
                                )
                            )
                    line = handle.readline()
                    words = line.strip().split()
                    # print("Still looking for seq for %s in line: %r" % (name, line))
                # Dealt with any coordinate header line, should now be sequence
                if not words:
                    # Should be sequence here, but perhaps its a short one?
                    if (
                        lengths[idx] < aln_length
                        and len("".join(seqs[idx])) == lengths[idx]
                    ):
                        # Is this actually allowed in the format? Personally I would
                        # expect a line with name and a block of trailing ~ here.
                        pass
                    else:
                        raise ValueError(
                            "Expected sequence for %s, got: %r" % (name, line)
                        )
                elif words[0] == name:
                    assert len(words) > 1, line
                    # print(i, name, repr(words))
                    seqs[idx].extend(words[1:])
                else:
                    raise ValueError("Expected sequence for %r, got: %r" % (name, line))
            # TODO - check the sequence lengths thus far are consistent
            # with blocks of 50?
            completed_length += 50
            line = handle.readline()
            if line.strip():
                raise ValueError("Expected blank line, got: %r" % line)

        # Skip over any whitespace at the end...
        while True:
            line = handle.readline()
            if not line:
                # End of file, no more alignments
                break
            elif not line.strip():
                # Blank line, ignore
                pass
            elif line.strip().split()[0] in known_headers:
                # Looks like the start of another alignment:
                self._header = line
                break
            else:
                raise ValueError("Unexpected line after GCG MSF alignment: %r" % line)

        # Combine list of strings into single string, remap gaps
        seqs = ["".join(s).replace("~", "-").replace(".", "-") for s in seqs]

        # Apply any trailing padding for short sequences
        padded = False
        for idx, (length, s) in enumerate(zip(lengths, seqs)):
            if len(s) < aln_length and len(s) == length:
                padded = True
                seqs[idx] = s + "-" * (aln_length - len(s))
        if padded:
            import warnings
            from Bio import BiopythonParserWarning

            warnings.warn(
                "One of more alignment sequences were truncated and have been gap padded",
                BiopythonParserWarning,
            )

        records = (
            SeqRecord(Seq(s), id=i, name=i, description=i, annotations={"weight": w},)
            for (i, s, w) in zip(ids, seqs, weights)
        )

        # This will check alignment lengths are self-consistent:
        align = MultipleSeqAlignment(records)
        # Check matches the header:
        if align.get_alignment_length() != aln_length:
            raise ValueError(
                "GCG MSF headers said alignment length %i, but have %i"
                % (aln_length, align.get_alignment_length())
            )
        return align
