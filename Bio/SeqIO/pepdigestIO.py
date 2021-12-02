# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
#
# This module is for reading .pepdigest format files as SeqRecord objects.
"""Bio.SeqIO support for the "pepdigest" file format.

You are expected to use this module via the Bio.SeqIO functions.
"""
import re

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .Interfaces import SequenceIterator


class PepdigestIterator(SequenceIterator):
    """TODO: class docstring."""

    main_header_pttrn = re.compile("#{2,}")
    digestion_title_pttrn = re.compile("#=+")
    digestion_end_pttrn = re.compile("#-+")

    def __init__(self, source):
        """Break up a pepdigest seqtable file into SeqRecord objects.

        Argument source is a file-like object opened in text mode or a path to a file.
        """
        return super().__init__(source)

    def parse(self, handle):
        """Start parsing the file, and return a SeqRecord generator."""
        # Reads the main header
        header_sep_count = 0
        while header_sep_count < 2:
            line = handle.readline()
            if self.main_header_pttrn.search(line) is not None:
                header_sep_count += 1

        records = self.iterate(handle)

        return records

    def parse_digestion_header(self, handle):
        """TODO."""
        # Read sequence name and positions
        line = handle.readline()
        while not line.startswith("# Sequence"):
            line = handle.readline()

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
        while self.digestion_title_pttrn.search(line) is None:
            line = handle.readline()

        return orig_seq_name, seq_to, seq_from, hitcount, enzyme

    def parse_digestion_results(self, handle, header):
        """TODO."""
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

    def iterate(self, handle):
        """Parse the file and generate SeqRecord objects.

        Reads
        """
        while True:
            line = handle.readline()
            if not line:
                break

            if self.digestion_title_pttrn.search(line) is not None:
                header = self.parse_digestion_header(handle)
                yield from self.parse_digestion_results(handle, header)


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest(verbose=0)
