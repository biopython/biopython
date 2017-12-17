# Copyright 2008-2010, 2012-2014, 2016-2017 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.AlignIO support for the "nexus" file format.

You are expected to use this module via the Bio.AlignIO functions (or the
Bio.SeqIO functions if you want to work directly with the gapped sequences).

See also the Bio.Nexus module (which this code calls internally),
as this offers more than just accessing the alignment or its
sequences as SeqRecord objects.
"""

from __future__ import print_function

from Bio.SeqRecord import SeqRecord
from Bio.Nexus import Nexus
from Bio.Align import MultipleSeqAlignment
from .Interfaces import AlignmentWriter
from Bio import Alphabet


# You can get a couple of example files here:
# http://www.molecularevolution.org/resources/fileformats/


# This is a generator function!
def NexusIterator(handle, seq_count=None):
    """Return SeqRecord objects from a Nexus file.

    Thus uses the Bio.Nexus module to do the hard work.

    You are expected to call this function via Bio.SeqIO or Bio.AlignIO
    (and not use it directly).

    NOTE - We only expect ONE alignment matrix per Nexus file,
    meaning this iterator will only yield one MultipleSeqAlignment.
    """
    n = Nexus.Nexus(handle)
    if not n.matrix:
        # No alignment found
        return

    # Bio.Nexus deals with duplicated names by adding a '.copy' suffix.
    # The original names and the modified names are kept in these two lists:
    assert len(n.unaltered_taxlabels) == len(n.taxlabels)

    if seq_count and seq_count != len(n.unaltered_taxlabels):
        raise ValueError("Found %i sequences, but seq_count=%i"
                         % (len(n.unaltered_taxlabels), seq_count))

    # TODO - Can we extract any annotation too?
    records = (SeqRecord(n.matrix[new_name], id=new_name,
                         name=old_name, description="")
               for old_name, new_name
               in zip(n.unaltered_taxlabels, n.taxlabels))
    # All done
    yield MultipleSeqAlignment(records, n.alphabet)


class NexusWriter(AlignmentWriter):
    """Nexus alignment writer.

    Note that Nexus files are only expected to hold ONE alignment
    matrix.

    You are expected to call this class via the Bio.AlignIO.write() or
    Bio.SeqIO.write() functions.
    """

    def write_file(self, alignments):
        """Use this to write an entire file containing the given alignments.

        Arguments:
         - alignments - A list or iterator returning MultipleSeqAlignment objects.
           This should hold ONE and only one alignment.

        """
        align_iter = iter(alignments)  # Could have been a list
        try:
            first_alignment = next(align_iter)
        except StopIteration:
            first_alignment = None
        if first_alignment is None:
            # Nothing to write!
            return 0

        # Check there is only one alignment...
        try:
            second_alignment = next(align_iter)
        except StopIteration:
            second_alignment = None
        if second_alignment is not None:
            raise ValueError("We can only write one Alignment to a Nexus file.")

        # Good.  Actually write the single alignment,
        self.write_alignment(first_alignment)
        return 1  # we only support writing one alignment!

    def write_alignment(self, alignment, interleave=None):
        """Write an alignment to file.

        Creates an empty Nexus object, adds the sequences
        and then gets Nexus to prepare the output.
        Default interleave behaviour: Interleave if columns > 1000
        --> Override with interleave=[True/False]
        """
        if len(alignment) == 0:
            raise ValueError("Must have at least one sequence")
        columns = alignment.get_alignment_length()
        if columns == 0:
            raise ValueError("Non-empty sequences are required")
        minimal_record = "#NEXUS\nbegin data; dimensions ntax=0 nchar=0; " \
                         + "format datatype=%s; end;"  \
                         % self._classify_alphabet_for_nexus(alignment._alphabet)
        n = Nexus.Nexus(minimal_record)
        n.alphabet = alignment._alphabet
        for record in alignment:
            n.add_sequence(record.id, str(record.seq))

        # Note: MrBayes may choke on large alignments if not interleaved
        if interleave is None:
            interleave = (columns > 1000)
        n.write_nexus_data(self.handle, interleave=interleave)

    def _classify_alphabet_for_nexus(self, alphabet):
        """Return 'protein', 'dna', or 'rna' based on the alphabet (PRIVATE).

        Raises an exception if this is not possible.
        """
        # Get the base alphabet (underneath any Gapped or StopCodon encoding)
        a = Alphabet._get_base_alphabet(alphabet)

        if not isinstance(a, Alphabet.Alphabet):
            raise TypeError("Invalid alphabet")
        elif isinstance(a, Alphabet.ProteinAlphabet):
            return "protein"
        elif isinstance(a, Alphabet.DNAAlphabet):
            return "dna"
        elif isinstance(a, Alphabet.RNAAlphabet):
            return "rna"
        else:
            # Must be something like NucleotideAlphabet or
            # just the generic Alphabet (default for fasta files)
            raise ValueError("Need a DNA, RNA or Protein alphabet")


if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest(verbose=0)
