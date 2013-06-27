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
from Bio.Alphabet import IUPAC, Gapped, Alphabet, HasStopCodon, generic_dna
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
            alpha = Gapped(alpha, gap_char=gap)
        if stop:
            alpha = HasStopCodon(alpha, stop_symbol=stop)
    else:
        raise TypeError("Only Nuclteotide Alphabet is accepted!")
    return alpha

default_alphabet = get_codon_alphabet(IUPAC.unambiguous_dna)

class CodonAlphabet(Alphabet):
    """Generic Codon Alphabet with a size of three"""
    size = 3
    letters = None

generic_codon_alphabet = CodonAlphabet()


def get_codon_alphabet(codon_table, gap_char="-"):
    letters = codon_table.forward_table.keys()
    letters.extend(codon_table.stop_codons)
    letters.extend(codon_table.start_codons)
    letters.append(gap_char*3)
    generic_codon_alphabet.letters = letters
    return generic_codon_alphabet

default_codon_alphabet = get_codon_alphabet(default_codon_table)


class CodonSeq(Seq):
    """CodonSeq is designed to be within the SeqRecords of a 
    CodonAlignment class. This most useful feature for Codon 
    Sequences bacause it slices three letters at once, i.e.
    codon slice.

    """
    def __init__(self, data, alphabet=default_codon_alphabet):

        # Check the alphabet of the input sequences
        # TODO:
        # set up a codon alphabet to verify more rigorously
        Seq.__init__(self, data.upper(), alphabet=alphabet)

        # check the length of the alignment to be a triple
        assert len(self) % 3 == 0, "Sequence length is not a triple number"

        # check alphabet
        # Not use Alphabet._verify_alphabet function because it 
        # only works for single alphabet
        for i in range(self.get_codon_num()):
            if self[i] not in alphabet.letters:
                raise ValueError("Sequence contain undefined letters from alphabet (%s)!"\
                        % self[i])
    
    def __getitem__(self, index):
        """get the `index`-th codon in from the self.seq
        """
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
            aa_index = range(self.get_codon_num())
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
        return len(self._data) / 3

    def toSeq(self, alphabet=generic_dna):
        return Seq(self._data, generic_dna)


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
                    self.get_alignment_length(), self.get_codon_num())]
        
        if rows <= 20:
            lines.extend([self._str_line(rec, length=20) for rec in self._records])
        else:
            lines.extend([self._str_line(rec, length=20) for rec in self._records[:18]])
            lines.append("...")
            lines.append(self._str_line(self._records[-1], length=60))
        return "\n".join(lines)


    def get_codon_num(self):
        return self.get_alignment_length() / 3

    def toMultipleSeqAlignment(self):
        """Return a MultipleSeqAlignment containing all the
        SeqRecord in the CodonAlignment using Seq to store 
        sequences
        """
        alignments = [SeqRecord(rec.seq.toSeq(), id=rec.id) for \
                rec in self._records]
        return MultipleSeqAlignment(alignments)


def _codons2re(codons):
    """Generate regular expression based on a given list of codons
    """
    reg = ''
    for i in izip(*codons):
        if len(set(i)) == 1:
            reg += ''.join(set(i))
        else:
            reg += '[' + ''.join(set(i)) + ']'
    return reg

def _get_aa_regex(codon_table, stop='*', unknown='X'):
    """Set up the regular expression of a given CodonTable for futher use.

    >>> from Bio.Data.CodonTable import generic_by_id
    >>> p = generic_by_id[1]
    >>> t = _get_aa_regex(p)
    >>> print t['A']
    GC[ACUTG]
    >>> print t['L']
    [CUT][UT][ACUTG]

    """
    from Bio.Data.CodonTable import CodonTable
    if not isinstance(codon_table, CodonTable):
        raise TypeError("Input table is not a instance of Bio.Data.CodonTable object")
    aa2codon = {}
    for codon, aa in codon_table.forward_table.iteritems():
        aa2codon.setdefault(aa, []).append(codon)

    for aa, codons in aa2codon.items():
        aa2codon[aa] = _codons2re(codons)
    aa2codon[stop] = _codons2re(codon_table.stop_codons)
    aa2codon[unknown] = '...'
    return aa2codon


def _check_corr(pro, nucl, gap_char='-', \
        codon_table=default_codon_table, complete_protein=False):
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
            pro_re += aa2re[aa]
    #TODO:
    #  1) Allow frameshift between protein sequences and nucleotides
    nucl_seq = str(nucl.seq.upper().ungap(gap_char))
    match = re.search(pro_re, nucl_seq)
    if match:
        # mode = 0, direct match
        return (match.span(), 0)
    else:
        # Might caused by mismatches, using Anchors to have a try
        anchor_len = 10 # adjust this value to test performance
        anchors = [pro.seq[i:(i+anchor_len)] for i in \
                range(0, len(pro.seq), anchor_len)]
        # if the last anchor is less than the specified anchor
        # size, we combine the penultimate and the last anchor
        # together as the last one.
        if len(anchors[-1]) < anchor_len:
            anchors[-1] = anchors[-2] + anchors[-1]

        pro_re = ""
        for i, anchor in enumerate(anchors):
            this_anchor_len = len(anchor)
            qcodon  = ""
            fncodon = ""
            ## dirty code to deal with the last anchor ##
            # as the last anchor is combined in the steps
            # above, we need to get the true last anchor to
            # pro_re
            if this_anchor_len == 10:
                for aa in str(anchor).replace(gap_char, ""):
                    if complete_protein is True and i == 0:
                        qcodon += _codons2re(codon_table.start_codons)
                        print qcodon
                        fncodon += aa2re['X']
                        continue
                    qcodon += aa2re[aa]
                    fncodon += aa2re['X']
            elif this_anchor_len > 10:
                last_qcodon = ""
                pos = 0
                for aa in anchor:
                    if aa != gap_char:
                        qcodon += aa2re[aa]
                        fncodon += aa2re['X']
                        if pos >= 10:
                            last_qcodon += aa2re[aa]
                    pos += 1
            match = re.search(qcodon, nucl_seq)
            if match:
                if this_anchor_len == 10:
                    pro_re += qcodon
                else:
                    pro_re += last_qcodon
            else:
                if this_anchor_len == 10:
                    pro_re += fncodon
                else:
                    pro_re += fncodon[10:]
        match = re.search(pro_re, nucl_seq)
        if match:
            print '1'
            # mode = 1, mismatch
            return (match.span(), 1)
        else:
            raise RuntimeError("Protein SeqRecord (%s) and Nucleotide SeqRecord (%s) do not match!" \
                    % (pro.id, nucl.id))
            

def _get_codon_aln(pro, nucl, span_mode, alphabet, gap_char="-", \
        codon_table=default_codon_table, complete_protein=False):
    """Generate codon alignment based on regular re match (PRIVATE)

    span_mode is a tuple returned by _check_corr. The first element
    is the span of a re search, and the second element is the mode
    for the match.

    mode
     - 0: direct match
     - 1: mismatch (no indels)

    """
    import re, warnings

    nucl_seq = nucl.seq.ungap(gap_char)
    codon_seq = ""
    span = span_mode[0]
    mode = span_mode[1]
    aa2re = _get_aa_regex(codon_table)
    if mode == 0 or mode == 1:
        if len(pro.seq.ungap(gap_char)) * 3 != (span[1] - span[0]):
            raise ValueError("Protein Record %s and Nucleotide Record %s do not match!" \
                    % (pro.id, nucl.id))
        aa_num = 0
        for aa in pro.seq:
            if aa == "-":
                codon_seq += "---"
            elif complete_protein is True and aa_num == 0:
                this_codon = nucl_seq._data[span[0]:(span[0]+3)]
                if not re.search(_codons2re[codon_table.start_codons], this_codon_upper()):
                    warnings.warn("start codon of %s (%s %d) does not correspond to %s (%s)" \
                            % (pro.id, aa, aa_num, nucl.id, this_codon))
                codon_seq += this_codon
                aa_num += 1
            else:
                this_codon = nucl_seq._data[(span[0] + 3*aa_num):(span[0]+3*(aa_num+1))]
                if not re.search(aa2re[aa], this_codon.upper()):
                    warnings.warn("%s (%s %d) does not correspond to %s (%s)" % (pro.id, aa, aa_num, nucl.id, this_codon))
                codon_seq += this_codon
                aa_num += 1
    return SeqRecord(CodonSeq(codon_seq, alphabet=alphabet), id=nucl.id)


def build(pro_align, nucl_seqs, corr_dict=None, gap_char='-', unknown='X', \
        codon_table=default_codon_table, alphabet=None, complete_protein=False):
    """Build a codon alignment from a protein alignment and corresponding
    nucleotide sequences

    Arguments:
     - pro_align  - a MultipleSeqAlignment object that stores protein alignment
     - nucl_align - an object returned by SeqIO.parse or SeqIO.index or a 
                    colloction of SeqRecord.
     - alphabet   - alphabet for the returned codon alignment
     - corr_dict  - a dict that maps protein id to nucleotide id

    Return a CodonAlignment object

    """
    # TODO
    # add an option to allow the user to specify the returned object?

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
    if corr_dict is None:
        if nucl_seqs.__class__.__name__ == "generator":
            # nucl_seqs will be a tuple if read by SeqIO.parse()
            nucl_seqs = tuple(nucl_seqs) 
        nucl_num = len(nucl_seqs)
        if pro_num > nucl_num:
            raise ValueError("More Number of SeqRecords in Protein Alignment (%d) than the Number of Nucleotide SeqRecords (%d) are found!" \
                    % (pro_num, nucl_num))

        if alphabet is None:
            alphabet = get_codon_alphabet(codon_table, gap_char=gap_char)
        # Determine the protein sequences and nucl sequences correspondance. If
        # nucl_seqs is a list, tuple or read by SeqIO.parse(), we assume the order
        # of sequences in pro_align and nucl_seqs are the same. If nucl_seqs is a
        # dict or read by SeqIO.index(), we match seqs in pro_align and those in 
        # nucl_seq by their id.
        if   nucl_seqs.__class__.__name__ in ("_IndexedSeqFileDict", "dict"):
            corr_method = 1
        elif nucl_seqs.__class__.__name__ in ("list", "tuple"):
            corr_method = 0
        else:
            raise TypeError("Nucl Sequences Error, Unknown type to assign correspondance method")
    else:
        if not isinstance(corr_dict, dict):
            raise TypeError("corr_dict should be a dict that corresponds protein id to nucleotide id!")
        if len(corr_dict) >= pro_num:
            if nucl_seqs.__class__.__name__ == "generator": # read by SeqIO.parse()
                from Bio import SeqIO
                nucl_seqs = SeqIO.to_dict(nucl_seqs)
            elif nucl_seqs.__class__.__name__ in ("list", "tuple"):
                nucl_seqs = {i.id: i for i in nucl_seqs}
            elif nucl_seqs.__class__.__name__ in ("_IndexedSeqFileDict", "dict"):
                pass
            else:
                raise TypeError("Nucl Sequences Error, Unknown type of Nucleotide Records!")
            corr_method = 2
        else:
            raise RuntimeError("Number of items in corr_dict (%n) is less than number of protein records (%n)" \
                    % (len(corr_dict), pro_num))

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
            raise ValueError("Protein Record %s cannot find a nucleotide sequence match, please check the id" \
                    % ', '.join(diff))
        else:
            pro_nucl_pair = []
            for pro_rec in pro_align:
                pro_nucl_pair.append((pro_rec, nucl_seqs[pro_rec.id]))
    elif corr_method == 2:
        pro_nucl_pair = []
        for pro_rec in pro_align:
            try:
                nucl_id = corr_dict[pro_rec.id]
            except KeyError:
                print "Protein record (%s) is not in corr_dict!" % pro_rec.id
                exit(1)
            pro_nucl_pair.append((pro_rec, nucl_seqs[nucl_id]))

    codon_aln = []
    for pair in pro_nucl_pair:
        # Beaware that the following span corresponds to a ungapped 
        # nucleotide sequence.
        corr_span = _check_corr(pair[0], pair[1], gap_char=gap_char, \
                codon_table=codon_table, complete_protein=complete_protein)
        if not corr_span:
            raise ValueError("Protein Record %s and Nucleotide Record %s do not match!" \
                    % (pair[0].id, pair[1].id))
        else:
            codon_rec = _get_codon_aln(pair[0], pair[1], corr_span, \
                    alphabet=alphabet, complete_protein=False)
            codon_aln.append(codon_rec)
    return CodonAlignment(codon_aln, alphabet=alphabet)

def toCodonAlignment(align, alphabet=default_codon_alphabet):
    """Function to convert a MultipleSeqAlignment to CodonAlignment.
    It is the user's responsibility to ensure all the requirement
    needed by CodonAlignment is met.

    """
    m = align._records[0]
    print m.seq.upper()
    rec = [SeqRecord(CodonSeq(str(i.seq), alphabet=alphabet), id=i.id) \
            for i in align._records]
    return CodonAlignment(rec, alphabet=align._alphabet)



if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()

