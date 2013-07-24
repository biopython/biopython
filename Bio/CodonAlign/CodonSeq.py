# Copyright 2013 by Zheng Ruan.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Code for dealing with Codon Seq.

CodonSeq class is interited from Seq class. This is the core class to
deal with sequences in CodonAlignment in biopython.

"""
__docformat__ = "epytext en"  # Don't just use plain text in epydoc API pages!


from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, Gapped, HasStopCodon, Alphabet, generic_dna, _ungap
from Bio.Data.CodonTable import generic_by_id

from CodonAlphabet import default_codon_alphabet, default_codon_table

class CodonSeq(Seq):
    """CodonSeq is designed to be within the SeqRecords of a 
    CodonAlignment class.

    CodonSeq is useful as it allows the user to specify
    reading frame when translate CodonSeq

    CodonSeq also accepts codon style slice by calling
    get_codon() method.

    Important: Ungapped CodonSeq can be any length if you
    specify the rf_table. Gapped CodonSeq should be a
    multiple of three.

    >>> codonseq = CodonSeq("AAATTTGGGCCAAATTT", rf_table=(0,3,6,8,11,14))
    >>> print codonseq.translate()
    KFGAKF

    test get_full_rf_table method

    >>> p = CodonSeq('AAATTTCCCGG-TGGGTTTAA', rf_table=(0, 3, 6, 9, 11, 14, 17))
    >>> full_rf_table = p.get_full_rf_table()
    >>> print full_rf_table
    [0, 3, 6, 9, 12, 15, 18]
    >>> print p.translate(rf_table=full_rf_table, ungap_seq=False)
    KFPPWV*
    >>> p = CodonSeq('AAATTTCCCGGGAA-TTTTAA', rf_table=(0, 3, 6, 9, 14, 17))
    >>> print p.get_full_rf_table()
    [0, 3, 6, 9, 15, 18]
    >>> p = CodonSeq('AAA------------TAA', rf_table=(0, 3)) 
    >>> print p.get_full_rf_table()
    [0, 3, 6, 9, 12, 15]

    """
    def __init__(self, data, alphabet=default_codon_alphabet, \
            gap_char="-", rf_table=None):
        # rf_table should be a tuple or list indicating the every
        # codon position along the sequence. For example:
        # sequence = 'AAATTTGGGCCAAATTT'
        # rf_table = (0, 3, 6, 8, 11, 14)
        # the translated protein sequences will be
        # AAA TTT GGG GCC AAA TTT
        #  K   F   G   A   K   F
        # Notice: rf_table applies to ungapped sequence. If there
        #   are gaps in the sequence, they will be discarded. This
        #   feature ensures the rf_table is independent of where the
        #   codon sequence appears in the alignment

        Seq.__init__(self, data.upper(), alphabet=alphabet)
        self.gap_char = gap_char

        # check the length of the alignment to be a triple
        if rf_table is None:
            seq_ungapped = self._data.replace(gap_char, "")
            if len(self) == 47:
                raise RuntimeError('trace back')
            assert len(self) % 3 == 0, "Sequence length is not a triple number"
            self.rf_table = filter(lambda x: x%3 == 0, range(len(seq_ungapped)))
            # check alphabet
            # Not use Alphabet._verify_alphabet function because it 
            # only works for single alphabet
            for i in self.rf_table:
                if self._data[i:i+3] not in alphabet.letters:
                    raise ValueError("Sequence contain undefined " \
                                  + "letters from alphabet (%s)! " \
                                  % self._data[i:i+3])
        else:
            if gap_char in self._data:
                assert  len(self) % 3 == 0, \
                        "Gapped sequence length is not a triple number"
            assert isinstance(rf_table, (tuple, list)), \
                    "rf_table should be a tuple or list object"
            assert all(isinstance(i, int) for i in rf_table), \
                    "elements in rf_table should be int that specify " \
                  + "the codon positions of the sequence"
            seq_ungapped = self._data.replace(gap_char, "")
            for i in rf_table:
                if seq_ungapped[i:i+3] not in alphabet.letters:
                    raise ValueError("Sequence contain undefined " \
                                  + "letters from alphabet (%s)! " \
                                  % seq_ungapped[i:i+3])
            self.rf_table = rf_table
    
    def __getitem__(self, index):
        #TODO: modify this to allow rf_table specification
        #      return a Seq object if index is not for codon (warning)
        seqlst = range(len(self._data))
        rf_table = [i for i in seqlst[index] if i in self.rf_table]
        rf_table = [i-rf_table[0] for i in rf_table]
        try:
            return CodonSeq(self._data[index], self.alphabet, rf_table=rf_table)
        except ValueError:
            # adjust alphabet
            return Seq(self._data[index], alphabet=generic_dna)

    def get_codon(self, index):
        """get the `index`-th codon in from the self.seq
        """
        if len(set([i % 3 for i in self.rf_table])) != 1:
            raise RuntimeError("frameshift detected. " \
                             + "CodonSeq object is not able to deal " \
                             + "with codon sequence with frameshift. " \
                             + "Plase use normal slice option ")
        if isinstance(index, int):
            if index != -1:
                return self._data[index*3:(index+1)*3]
            else:
                return self._data[index*3:]
        else:
        # This slice ensures that codon will always be the unit
        # in slicing (it won't change to other codon if you are 
        # using reverse slicing such as [::-1]).
        # The idea of the code below is to first map the slice
        # to amino acid sequence and then transform it into 
        # codon sequence.
            aa_index = range(len(self)/3)
            def cslice(p):
                aa_slice = aa_index[p]
                codon_slice = ''
                for i in aa_slice:
                    codon_slice += self._data[i*3:i*3+3]
                return codon_slice
            codon_slice = cslice(index)
            return CodonSeq(codon_slice, alphabet=self.alphabet)

    def get_codon_num(self):
        """Return the number of codons in the CodonSeq"""
        return len(self.rf_table)

    def translate(self, codon_table=default_codon_table, \
            stop_symbol="*", rf_table=None, ungap_seq=True):
        """Translate the CodonSeq based on the reading frame
        in rf_table. It is possible for the user to specify
        a rf_table at this point. If you want to include
        gaps in the translated sequence, this is the only
        way. ungap_seq should be set to true for this
        purpose.
        """
        amino_acids = []
        if ungap_seq is True:
            tr_seq = self._data.replace(self.gap_char, "")
        else:
            tr_seq = self._data
        if rf_table is None:
            rf_table = self.rf_table
        p = -1 #initiation
        for i in rf_table:
            if '---' == tr_seq[i:i+3]:
                amino_acids.append('-')
                continue
            elif '-' in tr_seq[i:i+3]:
                # considering two types of frameshift
                if p == -1 or p - i == 3:
                    p = i
                    codon = tr_seq[i:i+6].replace('-', '')[:3]
                elif p - i > 3:
                    codon = tr_seq[i:i+3]
                    p = i
            else:
                # normal condition without gaps
                codon = tr_seq[i:i+3]
                p = i
            if codon in codon_table.stop_codons:
                amino_acids.append(stop_symbol)
                continue
            try:
                amino_acids.append(codon_table.forward_table[codon])
            except KeyError:
                raise RuntimeError("Unknown codon detected (%s). Do you forget to speficy ungap_seq argument?" % codon)
        return "".join(amino_acids)

    def toSeq(self, alphabet=generic_dna):
        return Seq(self._data, generic_dna)

    def get_full_rf_table(self):
        """This function returns a full rf_table of the given
        CodonSeq records. A full rf_table is different from 
        normal rf_table in that it translate gaps in CodonSeq.
        It is helpful to construct alignment containing
        frameshift.
        """
        #assert self.gap_char in self._data, \
        #        "get_full_rf_table() method is only useful when your sequence contains gaps"
        full_rf_table = []
        accum = 0
        for i in filter(lambda x: x%3==0, range(len(self))):
            if self._data[i:i+3] == self.gap_char*3:
                full_rf_table.append(i)
            elif self._data[i:i+3] in self.alphabet.letters:
                full_rf_table.append(i)
                accum += 1
            else:
                # TODO: think about the last codon
                try:
                    nxt_shift = self.rf_table[accum+1]-self.rf_table[accum]-3
                except IndexError:
                    continue
                if nxt_shift < 0:
                    full_rf_table.append(i)
                    accum += 1
                elif nxt_shift == 0:
                    pre_shift = self.rf_table[accum]-self.rf_table[accum-1]-3
                    if pre_shift <= 0:
                        raise RuntimeError("Unexpected Codon %s", self._data[i:i+3])
                    else:
                        pass
                elif shift > 0:
                    pass
        return full_rf_table

    def ungap(self, gap=None):
        if hasattr(self.alphabet, "gap_char"):
            if not gap:
                gap = self.alphabet.gap_char
            elif gap != self.alphabet.gap_char:
                raise ValueError("Gap %s does not match %s from alphabet"
                        % (repr(gap), repr(self.alphabet.alphabet.gap_char)))
            alpha = _ungap(self.alphabet)
        elif not gap:
            raise ValueError("Gap character not given and not defined in alphabet")
        else:
            alpha = self.alphabet # modify!
        if len(gap) != 1 or not isinstance(gap, str):
            raise ValueError("Unexpected gap character, %s" % repr(gap))
        return CodonSeq(str(self._data).replace(gap, ""), alpha, rf_table=self.rf_table)


if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()

