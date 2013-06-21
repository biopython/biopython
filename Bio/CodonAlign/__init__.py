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

from Bio.Seq import Seq
from itertools import izip
from Bio.SeqRecord import SeqRecord
from Bio import Alphabet
from Bio.Alphabet import IUPAC, Gapped
from Bio.Align import MultipleSeqAlignment
from Bio.Data.CodonTable import generic_by_id


default_codon_table = generic_by_id[1]

def get_codon_alphabet(alphabet, gap="-", stop="*"):
    """function to get alignment alphabet for codon alignment. Only
    nucleotide alphabet is accepted. Raise an error when the type of 
    alphabet is incompatible.
    """
    from Bio.Alphabet import NucleotideAlphabet
    if isinstance(alphabet, NucleotideAlphabet):
        alpha = alphabet
        if gap:
            alpha = Alphabet.Gapped(alpha, gap_char=gap)
        if stop:
            alpha = Alphabet.HasStopCodon(alpha, stop_symbol=stop)
    else:
        raise TypeError("Only Nuclteotide Alphabet is accepted!")
    return alpha

default_alphabet = get_codon_alphabet(IUPAC.unambiguous_dna)


class CodonSeq(Seq):
    """CodonSeq is designed to be within the SeqRecords of a 
    CodonAlignment class. This most useful feature for Codon 
    Sequences bacause it slices three letters at once, i.e.
    codon slice.

    """
    def __init__(self, data, alphabet=default_alphabet):

        # Check the alphabet of the input sequences
        # TODO:
        # set up a codon alphabet to verify more rigorously
        Seq.__init__(self, data.upper(), alphabet=alphabet)
        # Not use Alphabet._verify_alphabet function because it 
        # only works for single alphabet
        for letter in self._data:
            if letter not in alphabet.letters:
                print letter
                raise ValueError("Sequence contain undefined letters from alphabet!")
    
        # check the length of the alignment to be a triple
        assert len(self) % 3 == 0, \
            "Sequence length is not a triple number"
    
    def __getitem__(self, index):
        """get the `index`-th codon in from the self.seq
        """
        if isinstance(index, int):
            if index != -1:
                return self._data[index*3:(index+1)*3]
            else:
                return self._data[index*3:]
        # is there a clever way to deal triple slice??
        # The following code is a little stupid.
        elif index.step:
            import warnings
            warnings.warn(
                "Slice method in CodonSeq object won't deal with step.\nReturn all codons from start(%s) to end(%s)" \
                        % (index.start, index.stop))
        def idx(p):
            return None if p is None else 3*p
        index = slice(idx(index.start), idx(index.stop), None)
        return CodonSeq(self._data[index], alphabet=self.alphabet)

    def get_codon_num(self):
        """Return the number of codons in the CodonSeq"""
        return len(self._data) / 3


class CodonAlignment(MultipleSeqAlignment):
    """Codon Alignment class that inherits from MultipleSeqAlignment.

    >>> from Bio.Alphabet import generic_dna
    >>> from Bio.SeqRecord import SeqRecord
    >>> from Bio.Alphabet import IUPAC, Gapped
    >>> a = SeqRecord(CodonSeq("AAAACGTCG", IUPAC.unambiguous_dna), id="Alpha")
    >>> b = SeqRecord(CodonSeq("AAA-CGTCG", Gapped(IUPAC.unambiguous_dna, gap_char="-")), id="Beta")
    >>> c = SeqRecord(CodonSeq("AAAAGGTGG", IUPAC.unambiguous_dna), id="Gamma")
    >>> print CodonAlignment([a, b, c])
    CodonAlignment Object
    HasStopCodon(Gapped(IUPACUnambiguousDNA(), '-'), '*') alignment with 3 rows and 9 columns (3 codons)
    AAAACGTCG Alpha
    AAA-CGTCG Beta
    AAAAGGTGG Gamma

    """
    def __init__(self, records, name=None, alphabet=default_alphabet):

        MultipleSeqAlignment.__init__(self, records, alphabet=alphabet)

        # check the type of the alignment to be nucleotide
        for rec in self:
            if not isinstance(rec.seq, CodonSeq):
                raise TypeError("CodonSeq object are expected in each SeqRecord in CodonAlignment")

        assert self.get_alignment_length() % 3 == 0, \
            "Alignment length is not a triple number"

    def _str_line(self, record):
        """Returns a truncated representation of SeqRecord storing 
        CodonSeq (PRIVATE).
        
        This is a PRIVATE function used by the __str__ method. The
        idea is the same with Alignment._str_line().
        """
        if len(record.seq) < 60:
            return "%s %s" % (record.seq, record.id)
        else:
            return "%s...%s %s" \
                    % (record.seq[:17], record.seq[-3:], record.id)

    def __str__(self):
        """Return a multi-line string summary of the alignment.

        This output is indicated to be readable, but large alignment 
        is shown truncated. A maximum of 20 rows (sequences) and
        60 columns (20 codons) are shown, with the record identifiers.
        This should fit nicely on a single screen. e.g.

        """
        rows = len(self._records)
        lines = ["CodonAlignment Object\n%s alignment with %i rows and %i columns (%i codons)"
                 % (str(self._alphabet), rows, \
                    self.get_alignment_length(), self.get_codon_num())]
        
        if rows <= 20:
            lines.extend([self._str_line(rec) for rec in self._records])
        else:
            lines.extend([self._str_line(rec) for rec in self._records[:18]])
            lines.append("...")
            lines.append(self._str_line(self._records[-1]))
        return "\n".join(lines)


    def get_codon_num(self):
        return self.get_alignment_length() / 3


def _get_aa_regex(codon_table, stop='*', unknown='X'):
    """Set up the regular expression of a given CodonTable for futher use.

    >>> from Bio.Data.CodonTable import generic_by_id
    >>> p = generic_by_id[1]
    >>> t = _get_aa_regex(p)
    >>> print t['A']
    (GC[ACUTG])
    >>> print t['L']
    ([CUT][UT][ACUTG])

    """
    from Bio.Data.CodonTable import CodonTable
    if not isinstance(codon_table, CodonTable):
        raise TypeError("Input table is not a instance of Bio.Data.CodonTable object")
    aa2codon = {}
    for codon, aa in codon_table.forward_table.iteritems():
        aa2codon.setdefault(aa, []).append(codon)

    def codons2re(codons):
        """Generate regular expression based on a given list of codons
        """
        reg = '('
        for i in izip(*codons):
            if len(set(i)) == 1:
                reg += ''.join(set(i))
            else:
                reg += '[' + ''.join(set(i)) + ']'
        return reg + ')'

    for aa, codons in aa2codon.items():
        aa2codon[aa] = codons2re(codons)
    aa2codon[stop] = codons2re(codon_table.stop_codons)
    aa2codon[unknown] = '...'
    return aa2codon


def _check_corr(pro, nucl, gap_char='-', codon_table=default_codon_table):
    """check if a give protein SeqRecord can be translated by another
    nucleotide SeqRecord.
    """
    import re
    from Bio.Alphabet import NucleotideAlphabet

    if (not isinstance(pro, SeqRecord)) or \
            (not isinstance(nucl, SeqRecord)):
        raise TypeError("_check_corr accept two SeqRecord object. Please check your input.")

    def get_alpha(alpha):
        if hasattr(alpha, 'alphabet'):
            return get_alpha(alpha.alphabet)
        else:
            return alpha

    if not isinstance(get_alpha(nucl.seq.alphabet), NucleotideAlphabet):
        raise TypeError("Alphabet for nucl should be an instance of NucleotideAlphabet, %s detected" \
                % str(nucl.seq.alphabet))

    aa2re  = _get_aa_regex(codon_table)
    pro_re = ""
    for aa in pro.seq:
        if aa != gap_char:
            pro_re += aa2re[aa][1:-1]
    #TODO:
    #  1) Allow mismatches between protein sequences and nucleotides
    #  2) Allow frameshift between protein sequences and nucleotides
    match = re.search(pro_re, str(nucl.seq.upper().ungap(gap_char)))
    if match:
        return match.span()


def _get_codon_aln(pro, nucl, span, gap_char="-", \
        codon_table=default_codon_table, mode=0):
    """Generate codon alignment based on regular re match (PRIVATE)

    mode represents the re match approach:
     - 0: direct match

    """
    nucl_seq = nucl.seq.ungap(gap_char)
    codon_seq = ""
    if mode == 0:
        if len(pro.seq.ungap(gap_char)) * 3 != (span[1] - span[0]):
            raise ValueError("Protein Record %s and Nucleotide Record %s do not match!" \
                    % (pro.id, nucl.id))
            #raise MissMatchError(pro.id, nucl.id)
        aa_num = 0
        for aa in pro.seq:
            if aa == "-":
                codon_seq += "---"
            else:
                codon_seq += nucl_seq._data[(span[0] + 3*aa_num):(span[0]+3*(aa_num+1))]
                aa_num += 1
    # If specialized codon alphabet are developed, the 
    # following code should be modified
    return SeqRecord(CodonSeq(codon_seq, \
            alphabet=Gapped(nucl_seq.alphabet, gap_char=gap_char)), \
            id=nucl.id)


def build(pro_align, nucl_seqs, gap_char='-', unknown='X', \
        codon_table=default_codon_table, alphabet=default_alphabet):
    """Build a codon alignment from a protein alignment and corresponding
    nucleotide sequences

    Arguments:
     - pro_align  - a MultipleSeqAlignment object that stores protein alignment
     - nucl_align - an object returned by SeqIO.parse or SeqIO.index or a 
                    colloction of SeqRecord.
     - alphabet   - alphabet for the returned codon alignment

    Return a CodonAlignment object
    """
    from Bio.Alphabet import ProteinAlphabet
    
    # check the type of object of pro_align
    if not isinstance(pro_align, MultipleSeqAlignment):
        raise TypeError("the first argument should be a MultipleSeqAlignment object")
    # check the alphabet of pro_align
    for pro in pro_align:
        if not isinstance(pro.seq.alphabet, ProteinAlphabet):
            raise TypeError("Alphabet Error!\nThe input alignment should be a *PROTEIN* alignment")
    # check whether the number of seqs in pro_align and nucl_seqs is the same
    pro_num = len(pro_align)
    if nucl_seqs.__class__.__name__ == "generator":
        # nucl_seqs will be a tuple if read by SeqIO.parse()
        nucl_seqs = tuple(nucl_seqs) 
    nucl_num = len(nucl_seqs)
    if pro_num > nucl_num:
        raise ValueError("More Number of SeqRecords in Protein Alignment (%d) than the Number of Nucleotide SeqRecords (%d) are found!" \
                % (pro_num, nucl_num))

    # Determine the protein sequences and nucl sequences correspondance. If
    # nucl_seqs is a list, tuple or read by SeqIO.parse(), we assume the order
    # of sequences in pro_align and nucl_seqs are the same. If nucl_seqs is a
    # dict or read by SeqIO.index(), we match seqs in pro_align and those in 
    # nucl_seq by their id.
    if   nucl_seqs.__class__.__name__ == "_IndexedSeqFileDict":
        corr_method = 1
    elif nucl_seqs.__class__.__name__ == "list":
        corr_method = 0
    elif nucl_seqs.__class__.__name__ == "tuple":
        corr_method = 0
    elif nucl_seqs.__class__.__name__ == "dict":
        corr_method = 1
    else:
        raise TypeError("Nucl Sequences Error, Unknown type to assign correspondance method")
    # set up pro-nucl correspondance based on corr_method
    # corr_method = 0, consecutive pairing
    if corr_method == 0:
        pro_nucl_pair = izip(pro_align, nucl_seqs)
    # corr_method = 1, keyword pairing
    elif corr_method == 1:
        nucl_id  = set(nucl_seqs.keys())
        pro_id = set([i.id for i in pro_align])
        # check if there is pro_id that does not have a nucleotide match
        if pro_id - nucl_id: 
            diff = pro_id - nucl_id
            raise ValueError("Protein Record %s cannot find a nucleotide sequence match, please check the id" % ', '.join(diff))
        else:
            pro_nucl_pair = []
            for pro_rec in pro_align:
                pro_nucl_pair.append((pro_rec, nucl_seqs[pro_rec.id]))

    codon_aln = []
    for pair in pro_nucl_pair:
        # Beaware that the following span corresponds to a ungapped 
        # nucleotide sequence.
        corr_span = _check_corr(pair[0], pair[1], gap_char=gap_char, \
                codon_table=codon_table)
        if not corr_span:
            raise ValueError("Protein Record %s and Nucleotide Record %s do not match!" \
                    % (pair[0].id, pair[1].id))
        else:
            codon_rec = _get_codon_aln(pair[0], pair[1], corr_span, \
                    codon_table=codon_table)
            codon_aln.append(codon_rec)
    return CodonAlignment(codon_aln, alphabet=alphabet)


if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()
