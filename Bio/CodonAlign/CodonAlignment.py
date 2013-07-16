# Copyright 2013 by Zheng Ruan.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Code for dealing with Codon Alignment.

CodonAlignment class is interited from MultipleSeqAlignment class. This is
the core class to deal with codon alignment in biopython.

"""
__docformat__ = "epytext en"  # Don't just use plain text in epydoc API pages!

from Bio.Align import MultipleSeqAlignment

from _utils import default_codon_table, default_codon_alphabet
from CodonSeq import CodonSeq

class CodonAlignment(MultipleSeqAlignment):
    """Codon Alignment class that inherits from MultipleSeqAlignment.

    >>> from Bio.Alphabet import generic_dna
    >>> from Bio.SeqRecord import SeqRecord
    >>> from Bio.Alphabet import IUPAC, Gapped
    >>> a = SeqRecord(CodonSeq("AAAACGTCG", alphabet=default_codon_alphabet), id="Alpha")
    >>> b = SeqRecord(CodonSeq("AAA---TCG", alphabet=default_codon_alphabet), id="Beta")
    >>> c = SeqRecord(CodonSeq("AAAAGGTGG", alphabet=default_codon_alphabet), id="Gamma")
    >>> print CodonAlignment([a, b, c])
    CodonAlphabet() CodonAlignment with 3 rows and 9 columns (3 codons)
    AAAACGTCG Alpha
    AAA---TCG Beta
    AAAAGGTGG Gamma

    """
    def __init__(self, records, name=None, alphabet=default_codon_alphabet):

        MultipleSeqAlignment.__init__(self, records, alphabet=alphabet)

        # check the type of the alignment to be nucleotide
        for rec in self:
            if not isinstance(rec.seq, CodonSeq):
                raise TypeError("CodonSeq object are expected in each SeqRecord in CodonAlignment")

        assert self.get_alignment_length() % 3 == 0, \
            "Alignment length is not a triple number"

    def __str__(self):
        """Return a multi-line string summary of the alignment.

        This output is indicated to be readable, but large alignment 
        is shown truncated. A maximum of 20 rows (sequences) and
        60 columns (20 codons) are shown, with the record identifiers.
        This should fit nicely on a single screen. e.g.

        """
        rows = len(self._records)
        lines = ["%s CodonAlignment with %i rows and %i columns (%i codons)"
                 % (str(self._alphabet), rows, \
                    self.get_alignment_length(), self.get_aln_length())]
        
        if rows <= 60:
            lines.extend([self._str_line(rec, length=60) for rec in self._records])
        else:
            lines.extend([self._str_line(rec, length=60) for rec in self._records[:18]])
            lines.append("...")
            lines.append(self._str_line(self._records[-1], length=60))
        return "\n".join(lines)


    def get_aln_length(self):
        return self.get_alignment_length() / 3

    def toMultipleSeqAlignment(self):
        """Return a MultipleSeqAlignment containing all the
        SeqRecord in the CodonAlignment using Seq to store 
        sequences
        """
        alignments = [SeqRecord(rec.seq.toSeq(), id=rec.id) for \
                rec in self._records]
        return MultipleSeqAlignment(alignments)


if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()

