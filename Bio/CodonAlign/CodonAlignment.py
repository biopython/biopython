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
from Bio.SeqRecord import SeqRecord

from CodonAlphabet import default_codon_table, default_codon_alphabet
from CodonSeq import CodonSeq, cal_dn_ds

class CodonAlignment(MultipleSeqAlignment):
    """Codon Alignment class that inherits from MultipleSeqAlignment.

    >>> from Bio.Alphabet import generic_dna
    >>> from Bio.SeqRecord import SeqRecord
    >>> from Bio.Alphabet import IUPAC, Gapped
    >>> a = SeqRecord(CodonSeq("AAAACGTCG", alphabet=default_codon_alphabet), id="Alpha")
    >>> b = SeqRecord(CodonSeq("AAA---TCG", alphabet=default_codon_alphabet), id="Beta")
    >>> c = SeqRecord(CodonSeq("AAAAGGTGG", alphabet=default_codon_alphabet), id="Gamma")
    >>> print CodonAlignment([a, b, c])
    CodonAlphabet(Standard) CodonAlignment with 3 rows and 9 columns (3 codons)
    AAAACGTCG Alpha
    AAA---TCG Beta
    AAAAGGTGG Gamma

    """
    def __init__(self, records, name=None, alphabet=default_codon_alphabet):

        MultipleSeqAlignment.__init__(self, records, alphabet=alphabet)

        # check the type of the alignment to be nucleotide
        for rec in self:
            if not isinstance(rec.seq, CodonSeq):
                raise TypeError("CodonSeq object are expected in each "
                                "SeqRecord in CodonAlignment")

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
            lines.extend([self._str_line(rec, length=60) \
                    for rec in self._records])
        else:
            lines.extend([self._str_line(rec, length=60) \
                    for rec in self._records[:18]])
            lines.append("...")
            lines.append(self._str_line(self._records[-1], length=60))
        return "\n".join(lines)

    def __getitem__(self, index):
        """Return a CodonAlignment object for single indexing
        """
        if isinstance(index, int):
            return self._records[index]
        elif isinstance(index, slice):
            return CodonAlignment(self._records[index], self._alphabet)
        elif len(index) != 2:
            raise TypeError("Invalid index type.")
        # Handle double indexing
        row_index, col_index = index
        if isinstance(row_index, int):
            return self._records[row_index][col_index]
        elif isinstance(col_index, int):
            return "".join(str(rec[col_index]) for rec in \
                                                    self._records[row_index])
        else:
            return MultipleSeqAlignment((rec[col_index] for rec in \
                                                    self._records[row_index]),
                                         self._alphabet)

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

    def get_dn_ds_matrix(self, method="NG86"):
        from numpy import zeros
        size = len(self._records)
        dN = zeros((size, size))
        dS = zeros((size, size))
        for i in range(size):
            for j in range(i):
                if i != j:
                    rec1 = self._records[i].seq
                    rec2 = self._records[j].seq
                    dn, ds = cal_dn_ds(rec1, rec2, method=method)
                    dN[i,j] = dn
                    dN[j,i] = dn
                    dS[i,j] = ds
                    dS[j,i] = ds
        return dN, dS


def toCodonAlignment(align, alphabet=default_codon_alphabet):
    """Function to convert a MultipleSeqAlignment to CodonAlignment.
    It is the user's responsibility to ensure all the requirement
    needed by CodonAlignment is met.

    """
    rec = [SeqRecord(CodonSeq(str(i.seq), alphabet=alphabet), id=i.id) \
             for i in align._records]
    return CodonAlignment(rec, alphabet=align._alphabet)


def mktest(codon_aln1, codon_aln2, codon_table=default_codon_table,
        alpha=0.05):
    """McDonald-Kreitman test for neutrality
    """
    if not isinstance(codon_aln1, CodonAlignment) or \
            not isinstance(codon_aln2, CodonAlignment):
        raise TypeError("mktest accept two CodonAlignment object ({0}, {1} "
                        "detected)".format(type(codon_aln1), type(codon_aln2))
                        )
    if codon_aln1.get_alignment_length() != codon_aln2.get_alignment_length():
        raise RuntimeError("Two CodonAlignment object for mktest should be of"
                           " equal length.")

if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()

