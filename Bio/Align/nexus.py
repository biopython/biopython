# Copyright 2006-2016 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.Align support for the "nexus" file format.

You are expected to use this module via the Bio.Align functions.

See also the Bio.Nexus module (which this code calls internally),
as this offers more than just accessing the alignment or its
sequences as SeqRecord objects.
"""
import Bio
from Bio.Align import Alignment
from Bio.Align import interfaces
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Nexus import Nexus


class AlignmentWriter(interfaces.AlignmentWriter):
    """Nexus alignment writer.

    Note that Nexus files are only expected to hold ONE alignment
    matrix.

    You are expected to call this class via Bio.Align.write().
    """

    def __init__(self, source):
        """Create an AlignmentWriter object.

        Arguments:
         - source   - input data or file name

        """
        super().__init__(source, mode="w")

    def write_file(self, alignments):
        count = super().write_file(alignments, mincount=1, maxcount=1)
        return count

    def write_header(self, alignments):
        """Use this to write the file header."""
        stream = self.stream
        stream.write("#NEXUS\n")

    def write_alignment(self, alignment, interleave=None):
        """Write an alignment to file.

        Creates an empty Nexus object, adds the sequences
        and then gets Nexus to prepare the output.
        Default interleave behavior: Interleave if columns > 1000
        --> Override with interleave=[True/False]
        """
        nseqs, length = alignment.shape
        if nseqs == 0:
            raise ValueError("Must have at least one sequence")
        if length == 0:
            raise ValueError("Non-empty sequences are required")

        stream = self.stream
        rows, columns = alignment.shape
        if rows == 0:
            raise ValueError("Must have at least one sequence")
        if columns == 0:
            raise ValueError("Non-empty sequences are required")
        datatype = self._classify_mol_type_for_nexus(alignment)
        minimal_record = (
            "begin data; dimensions ntax=0 nchar=0; format datatype=%s; end;" % datatype
        )
        n = Nexus.Nexus(minimal_record)
        for record, aligned_sequence in zip(alignment.sequences, alignment):
            # Sanity test sequences (should this be even stricter?)
            if datatype == "dna" and "U" in record.seq:
                raise ValueError(f"{record.id} contains U, but DNA alignment")
            elif datatype == "rna" and "T" in record.seq:
                raise ValueError(f"{record.id} contains T, but RNA alignment")
            n.add_sequence(record.id, aligned_sequence)

        # Note: MrBayes may choke on large alignments if not interleaved
        if interleave is None:
            interleave = columns > 1000
        n.write_nexus_data(stream, interleave=interleave)

    def _classify_mol_type_for_nexus(self, alignment):
        """Return 'protein', 'dna', or 'rna' based on records' molecule type (PRIVATE).

        All the records must have a molecule_type annotation, and they must
        agree.

        Raises an exception if this is not possible.
        """
        values = {
            sequence.annotations.get("molecule_type", None)
            for sequence in alignment.sequences
        }
        if all(_ and "DNA" in _ for _ in values):
            return "dna"  # could have been a mix of "DNA" and "gDNA"
        elif all(_ and "RNA" in _ for _ in values):
            return "rna"  # could have been a mix of "RNA" and "mRNA"
        elif all(_ and "protein" in _ for _ in values):
            return "protein"
        else:
            raise ValueError("Need the molecule type to be defined")


class AlignmentIterator(interfaces.AlignmentIterator):
    """Nexus alignment iterator."""

    def __init__(self, source):
        """Create an AlignmentIterator object.

        Arguments:
         - source   - input data or file name

        """
        super().__init__(source, mode="t", fmt="Nexus")
        stream = self.stream
        try:
            line = next(stream)
        except StopIteration:
            raise ValueError("Empty file.") from None

        if line.strip() != "#NEXUS":
            raise ValueError("File does not start with NEXUS header.")

    def parse(self, stream):
        """Parse the next alignment from the stream

        This uses the Bio.Nexus module to do the hard work.

        You are expected to call this function via Bio.Align
        (and not use it directly).

        NOTE - We only expect ONE alignment matrix per Nexus file,
        meaning this iterator will only yield one Alignment.
        """
        n = Nexus.Nexus(stream)
        if not n.matrix:
            # No alignment found
            return

        # Bio.Nexus deals with duplicated names by adding a '.copy' suffix.
        # The original names and the modified names are kept in these two lists:
        assert len(n.unaltered_taxlabels) == len(n.taxlabels)

        # TODO - Can we extract any annotation too?
        if n.datatype in ("dna", "nucleotide"):
            annotations = {"molecule_type": "DNA"}
        elif n.datatype == "rna":
            annotations = {"molecule_type": "RNA"}
        elif n.datatype == "protein":
            annotations = {"molecule_type": "protein"}
        else:
            annotations = None
        aligned_seqs = [str(n.matrix[new_name]) for new_name in n.taxlabels]
        records = [
            SeqRecord(
                n.matrix[new_name].replace("-", ""),
                id=old_name,
                annotations=annotations,
            )
            for old_name, new_name in zip(n.unaltered_taxlabels, n.taxlabels)
        ]
        coordinates = Alignment.infer_coordinates(aligned_seqs)
        alignment = Alignment(records, coordinates)
        yield alignment
