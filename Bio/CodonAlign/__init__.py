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

#from Bio.Seq import Seq
from itertools import izip
from Bio.SeqRecord import SeqRecord
from Bio import Alphabet
from Bio.Alphabet import IUPAC
from Bio.Align import MultipleSeqAlignment
from Bio.Data.CodonTable import generic_by_id

default_codon_table = generic_by_id[1]

def get_codon_alphabet(alphabet, gap="-", stop="*"):
    """function to get alignment alphabet for codon alignment. Only nucleotide
    alphabet is accepted. Raise an error with the type of alphabet is
    incompatible.
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


class CodonAlignment(MultipleSeqAlignment):
    """Codon Alignment class that inherits from MultipleSeqAlignment.

    >>> from Bio.Alphabet import generic_dna
    >>> from Bio.Seq import Seq
    >>> from Bio.SeqRecord import SeqRecord
    >>> a = SeqRecord(Seq("AAAACGTCG", generic_dna), id="Alpha")
    >>> b = SeqRecord(Seq("AAA-CGTCG", generic_dna), id="Beta")
    >>> c = SeqRecord(Seq("AAAAGGTGG", generic_dna), id="Gamma")
    >>> print CodonAlignment([a, b, c])
    HasStopCodon(Gapped(IUPACUnambiguousDNA(), '-'), '*') alignment with 3 rows and 9 columns
    AAAACGTCG Alpha
    AAA-CGTCG Beta
    AAAAGGTGG Gamma

    """
    def __init__(self, records, name=None, alphabet=default_alphabet):

        MultipleSeqAlignment.__init__(self, records, alphabet=alphabet)

        # check the type of the alignment to be nucleotide
        for rec in self:
            rec.seq.alphabet = alphabet
            rec.seq = rec.seq.upper()
            assert Alphabet._verify_alphabet(rec.seq), \
                "%s is incompatible with the %s\n" % (rec.id, str(alphabet))
    
        # check the length of the alignment to be a triple
        assert self.get_alignment_length() % 3 == 0, \
            "Alignment length is not a triple number"

class NumError(Exception):
    """NumError indicates the number of protein and nucleotide sequences
    is not the same.
    """
    def __init__(self, pro_num, nucl_num):
        self.pro_num  = pro_num
        self.nucl_num = nucl_num
    def __str__(self):
        return repr("Number of Seqs in Protein Alignment (%d) and Number of Nucleotide Seqs (%d) are not the same" \
                     % (self.pro_num, self.nucl_num))

class MissMatchError(Exception):
    """MissMatchError indicates the protein sequence and nucleotide
    sequence do not match.

    TODO:
      1) Modify this Exception to output a detailed information about 
         the discrepancy between prot and nucl
    """
    def __init__(self, pro_id, nucl_id):
        self.pro_id = pro_id
        self.nucl_id = nucl_id
    def __str__(self):
        return repr("Protein Record %s and Nucleotide Record %s do not match!" \
        % (self.pro_id, self.nucl_id))


def get_aa_regex(codon_table, stop='*', unknown='X'):
    """Set up the regular expression of a given CodonTable for futher use.

    >>> from __init__ import get_aa_regex
    >>> from Bio.Data.CodonTable import generic_by_id
    >>> p = generic_by_id[1]
    >>> t = get_aa_regex(p)
    >>> print t['A']
    (GC[ACUTG])
    >>> print t['L']
    ([CUT][UT][ACUTG])

    """
    from Bio.Data.CodonTable import CodonTable
    if not isinstance(codon_table, CodonTable):
        raise TypeError("""
        Input table is not a instance of Bio.Data.CodonTable object""")
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


def check_corr(pro, nucl, gap_char='-', codon_table=default_codon_table):
    """check if a give protein SeqRecord can be translated by another \
    nucleotide SeqRecord.
    """
    import re
    from Bio.Alphabet import NucleotideAlphabet

    if (not isinstance(pro, SeqRecord)) or \
            (not isinstance(nucl, SeqRecord)):
        raise TypeError("""
        check_corr accept two SeqRecord object. Please check your input.""")
    if not isinstance(nucl.seq.alphabet, NucleotideAlphabet):
        raise TypeError("Alphabet for nucl should be an instance of \
                NucleotideAlphabet, %s detected" % str(nucl.seq.alphabet))

    aa2re  = get_aa_regex(codon_table)
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

def get_codon_aln(pro, nucl, span, gap_char="-", \
        codon_table=default_codon_table, mode=0):
    """Generate codon alignment based on regular re match (PRIVATE)

    mode represents the re match approach:
     - 0: direct match

    """
    from Bio.Seq import Seq
    nucl_seq = nucl.seq.ungap(gap_char)
    codon_seq = ""
    if mode == 0:
        if len(pro.seq.ungap(gap_char)) * 3 != (span[1] - span[0]):
            raise MissMatchError(pro.id, nucl.id)
        aa_num = 0
        for aa in pro.seq:
            #p = len(codon_seq)
            if aa == "-":
                codon_seq += "---"
                #print aa
                #print len(codon_seq) - p
            else:
                codon_seq += nucl_seq[(span[0] + 3*aa_num):(span[0]+3*(aa_num+1))]
                aa_num += 1
                #print aa
                #print len(codon_seq) - p
        #p = 0
    return SeqRecord(codon_seq, id=nucl.id)


def build(pro_align, nucl_seqs, gap_char='-', unknown='X', \
        codon_table=default_codon_table, alphabet=default_alphabet):
    """Build a codon alignment from a protein alignment and corresponding
    nucleotide sequences

    Arguments:
     - pro_align  - a MultipleSeqAlignment object that stores protein alignment
     - nucl_align - an object returned by SeqIO.parse or SeqIO.index or a 
                    colloction of SeqRecord.

    Return a CodonAlignment object
    """
    from Bio.Alphabet import ProteinAlphabet
    #from Bio._utils import iterlen
    
    # check the type of object of pro_align
    if not isinstance(pro_align, MultipleSeqAlignment):
        raise TypeError("the first argument should be a MultipleSeqAlignment \
                object")
    # check the alphabet of pro_align
    for pro in pro_align:
        if not isinstance(pro.seq.alphabet, ProteinAlphabet):
            raise TypeError("""Alphabet Error!\nThe first argument should be a \
                    protein alignment""")
    # check whether the number of seqs in pro_align and nucl_seqs is the same
    pro_num = len(pro_align)
    if nucl_seqs.__class__.__name__ == "generator":
        nucl_seqs = tuple(nucl_seqs)
    nucl_num = len(nucl_seqs)
    if pro_num != nucl_num:
        raise NumError(pro_num, nucl_num)

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
        raise TypeError("Nucl Sequences Error, Unknown type to assign \
                correspondance method")
    
    # set up pro-nucl correspondance based on corr_method
    pro_nucl_pair = izip(pro_align, nucl_seqs)
    codon_aln = []
    for pair in pro_nucl_pair:
        # Beaware that the following span corresponds to a ungapped 
        # nucleotide sequence.
        corr_span = check_corr(pair[0], pair[1], gap_char=gap_char, \
                codon_table=codon_table)
        if not corr_span:
            raise MissMatchError(pair[0].id, pair[1].id)
        else:
            codon_rec = get_codon_aln(pair[0], pair[1], corr_span, \
                    codon_table=codon_table)
            codon_aln.append(codon_rec)
    return CodonAlignment(codon_aln, alphabet=alphabet)


if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()
