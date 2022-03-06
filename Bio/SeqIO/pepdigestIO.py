# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
#
# This module is for reading .pepdigest format files as SeqRecord objects.
"""Bio.SeqIO support for the "pepdigest" file format.

You are expected to use this module via the Bio.SeqIO functions.

The "pepdigest" format is the default output format ("seqtable" format) for the
EMBOOS pepdigest tool, which performs in silico digestion of proteins.

See Also:
https://www.bioinformatics.nl/cgi-bin/emboss/pepdigest
https://www.bioinformatics.nl/cgi-bin/emboss/help/pepdigest
"""
import re

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .Interfaces import SequenceIterator


class PepdigestIterator(SequenceIterator):
    """Parser for EMBOSS pepdigest files."""

    main_header_pttrn = re.compile("#{2,}")
    digestion_title_pttrn = re.compile("#=+")
    max_iter = 100

    def __init__(self, source):
        """Break up a pepdigest seqtable file into SeqRecord objects.

        Argument source is a file-like object opened in text mode or a path to a file.
        """
        return super().__init__(source, mode="t")

    def parse(self, handle):
        """Parse the file and generate SeqRecord objects.

        Reads the file by reading each block of results, with each peptide
        produced by the in silico digestion being returned as a SeqRecord.
        """
        # Check for empty file
        line = handle.readline()
        if not line:
            raise ValueError("Empty file.")

        # Reads the main header
        header_sep_count = 0
        i = 0
        while header_sep_count < 2:
            if self.main_header_pttrn.search(line) is not None:
                header_sep_count += 1
            i += 1
            if i > self.max_iter:
                raise ValueError("Broken header in start of file.")
            line = handle.readline()

        while True:
            line = handle.readline()
            if not line:
                break

            if self.digestion_title_pttrn.search(line) is not None:
                try:
                    header = self.parse_digestion_header(handle)
                except IndexError:
                    raise ValueError("File has broken digestion results header")
                yield from self.parse_digestion_results(handle, header)

    def parse_digestion_header(self, handle):
        """Parse the header of a digestion.

        Return usefull information for parsing the results.
        """
        # Read sequence name and positions
        i = 0
        line = handle.readline()
        while not line.startswith("# Sequence"):
            line = handle.readline()
            i += 1
            if i > self.max_iter:
                raise ValueError(
                    "Can't find '# Sequence' section in digestion results header"
                )

        line = line.rsplit("# Sequence:")[1].strip()
        orig_seq_name, line = line.split(None, 1)
        line, seq_to = line.split("to:")
        seq_to = int(seq_to)
        seq_from = int(line.replace("from:", ""))

        # Read HitCount
        line = handle.readline()
        hitcount = int(line.replace("# HitCount:", ""))

        # Read enzyme utilized in digestion
        handle.readline()
        line = handle.readline()
        enzyme = line.split("with")[1].split("yields")[0].strip()

        # Read end of header
        i = 0
        while self.digestion_title_pttrn.search(line) is None:
            line = handle.readline()
            i += 1
            if i > self.max_iter:
                raise ValueError("Malformed digestion results header")

        return orig_seq_name, seq_to, seq_from, hitcount, enzyme

    def parse_digestion_results(self, handle, header):
        """Parse the digested peptides."""
        orig_seq_name, seq_to, seq_from, hitcount, enzyme = header
        # skip one line
        handle.readline()

        header = handle.readline().split()
        for i in range(hitcount):
            seq_line = handle.readline()
            start, end, mol_w, cterm, nterm, seq = seq_line.split()

            pep_name = f"{orig_seq_name}_{i+1}"
            yield SeqRecord(
                Seq(seq),
                id=pep_name,
                name=pep_name,
                annotations={
                    "Start": int(start),
                    "End": int(end),
                    "Mol_W": float(mol_w),
                    "Cterm": cterm,
                    "Nterm": nterm,
                    "source_sequence": orig_seq_name,
                    "source_seq_from": seq_from,
                    "source_seq_to": seq_to,
                    "digestion_enzyme": enzyme,
                },
                description="generated from pepdigest",
            )


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest(verbose=0)
