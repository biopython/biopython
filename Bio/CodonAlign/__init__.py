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

from itertools import izip

from Bio.Seq import Seq
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
            assert len(self) % 3 == 0, "Sequence length is not a triple number"
            self.rf_table = filter(lambda x: x%3 == 0, range(len(seq_ungapped)))
            # check alphabet
            # Not use Alphabet._verify_alphabet function because it 
            # only works for single alphabet
            for i in self.rf_table:
                if self[i:i+3] not in alphabet.letters:
                    raise ValueError("Sequence contain undefined " \
                                  + "letters from alphabet (%s)! " \
                                  % self[i:i+3])
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
        return self._data[index]

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
        return len(self.rf_table) / 3

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
        assert self.gap_char in self._data, \
                "get_full_rf_table() method is only useful when your sequence contains gaps"
        full_rf_table = []
        accum = 0
        for i in filter(lambda x: x%3==0, range(len(self))):
            if self[i:i+3] == self.gap_char*3:
                full_rf_table.append(i)
            elif self[i:i+3] in self.alphabet.letters:
                full_rf_table.append(i)
                accum += 1
            else:
                # think about the last codon
                # code is dirty right now
                #print self[i-100:i]
                #print i
                nxt_shift = self.rf_table[accum+1]-self.rf_table[accum]-3
                if nxt_shift < 0:
                    full_rf_table.append(i)
                elif nxt_shift == 0:
                    #print self.rf_table[accum-1]
                    #print self.rf_table[accum]
                    pre_shift = self.rf_table[accum]-self.rf_table[accum-1]-3
                    if pre_shift <= 0:
                        raise RuntimeError("Unexpected Codon %s", self[i:i+3])
                    else:
                        pass
                elif shift > 0:
                    pass
        #print len(full_rf_table)
        return full_rf_table


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
        # Might caused by mismatches or frameshift, using anchors to
        # have a try
        anchor_len = 10 # adjust this value to test performance
        pro_seq = str(pro.seq).replace(gap_char, "")
        anchors = [pro_seq[i:(i+anchor_len)] for i in \
                range(0, len(pro_seq), anchor_len)]
        # if the last anchor is less than the specified anchor
        # size, we combine the penultimate and the last anchor
        # together as the last one.
        # TODO: modify this to deal with short sequence with only
        # one anchor.
        if len(anchors[-1]) < anchor_len:
            anchors[-1] = anchors[-2] + anchors[-1]

        pro_re = []
        anchor_distance = 0
        anchor_pos = []
        for i, anchor in enumerate(anchors):
            this_anchor_len = len(anchor)
            qcodon  = ""
            fncodon = ""
            ## dirty code to deal with the last anchor ##
            # as the last anchor is combined in the steps
            # above, we need to get the true last anchor to
            # pro_re
            if this_anchor_len == anchor_len:
                for aa in anchor:#str(anchor).replace(gap_char, ""):
                    if complete_protein is True and i == 0:
                        qcodon += _codons2re(codon_table.start_codons)
                        fncodon += aa2re['X']
                        continue
                    qcodon += aa2re[aa]
                    fncodon += aa2re['X']
                match = re.search(qcodon, nucl_seq)
            elif this_anchor_len > anchor_len:
                last_qcodon = ""
                last_fcodon = ""
                for j in range(anchor_len, len(anchor)):
                    last_qcodon += aa2re[anchor[j]]
                    last_fcodon += aa2re['X']
                match = re.search(last_qcodon, nucl_seq)
            # build full_pro_re from anchors
            if match:
                anchor_pos.append((match.start(), match.end(), i))
                if this_anchor_len == anchor_len:
                    pro_re.append(qcodon)
                else:
                    pro_re.append(last_qcodon)
            else:
                if this_anchor_len == anchor_len:
                    pro_re.append(fncodon)
                else:
                    pro_re.append(last_fcodon)
                #anchor_pos.append(-1)

        full_pro_re = "".join(pro_re)
        match = re.search(full_pro_re, nucl_seq)
        if match:
            # mode = 1, mismatch
            return (match.span(), 1)
        else:
            # check frames of anchors
            shift_id = [chr(i) for i in range(97,107)]
            shift_id_pos = 0
            for i in range(len(anchor_pos)-1):
                # TODO: think about the first and last anchor (mismatch)
                shift_val = (anchor_pos[i+1][0] - anchor_pos[i][0]) % anchor_len
                #if shift_val != 0:
                if shift_val > 0:
                    # obtain shifted anchor and corresponding nucl
                    sh_anc = "".join(anchors[anchor_pos[i][2]+1:anchor_pos[i+1][2]])
                    sh_nuc = nucl_seq[anchor_pos[i][1]:anchor_pos[i+1][0]]
                    # re substring matching doesn't allow number as id
                    # use ascii instead
                    qcodon = "^(?P<a>.*)"
                    id_dict = {'a': 0}
                    for j, aa in enumerate(sh_anc):
                        qcodon += aa2re[aa] + "(?P<" + chr(j+98) + ">.*)"
                        id_dict[chr(j+98)] = j+1
                    qcodon += "$"
                    match = re.search(qcodon, sh_nuc)
                    if match:
                        anc_shift_pos = [num for id, num in id_dict.iteritems() \
                                if match.group(id) != ""]
                        #print anc_shift_pos
                        if 0 in anc_shift_pos:
                            qcodon = "(?P<" + shift_id[shift_id_pos] + ">.*)"
                            shift_id_pos += 1
                        else:
                            qcodon = ""
                        for j, aa in enumerate(sh_anc):
                            qcodon +=  aa2re[aa]
                            if j+1 in anc_shift_pos:
                                #print shift_id
                                qcodon += "(?P<" + shift_id[shift_id_pos] + ">.*)"
                                shift_id_pos += 1
                        pro_re[anchor_pos[i][2]+1:anchor_pos[i+1][2]] = [qcodon]
                    else:
                        # failed to find a match (frameshift)
                        import warnings
                        warnings.warn("Frameshift detection failed")
                elif shift_val == -1:
                    print nucl_seq[anchor_pos[i]+3*anchor_len:anchor_pos[i+1]]
                elif shift_val == -2:
                    print nucl_seq[anchor_pos[i]+3*anchor_len:anchor_pos[i+1]]
            full_pro_re = "".join(pro_re)
            match = re.search(full_pro_re, nucl_seq)
            if match:
                return (match.span(), 2, match)
            else:
                raise RuntimeError("Protein SeqRecord (%s) and Nucleotide SeqRecord (%s) do not match!" \
                        % (pro.id, nucl.id))
            

def _get_codon_rec(pro, nucl, span_mode, alphabet, gap_char="-", \
        codon_table=default_codon_table, complete_protein=False):
    """Generate codon alignment based on regular re match (PRIVATE)

    span_mode is a tuple returned by _check_corr. The first element
    is the span of a re search, and the second element is the mode
    for the match.

    mode
     - 0: direct match
     - 1: mismatch (no indels)
     - 2: frameshift

    """
    import re, warnings

    nucl_seq = nucl.seq.ungap(gap_char)
    codon_seq = ""
    span = span_mode[0]
    mode = span_mode[1]
    aa2re = _get_aa_regex(codon_table)
    if mode in (0, 1):
        if len(pro.seq.ungap(gap_char))*3 != (span[1]-span[0]):
            raise ValueError("Protein Record %s and Nucleotide Record %s do not match!" \
                    % (pro.id, nucl.id))
        aa_num = 0
        for aa in pro.seq:
            if aa == "-":
                codon_seq += "---"
            elif complete_protein is True and aa_num == 0:
                this_codon = nucl_seq._data[span[0]:span[0]+3]
                if not re.search(_codons2re[codon_table.start_codons], this_codon_upper()):
                    warnings.warn("start codon of %s (%s %d) does not correspond to %s (%s)" \
                            % (pro.id, aa, aa_num, nucl.id, this_codon))
                codon_seq += this_codon
                aa_num += 1
            else:
                this_codon = nucl_seq._data[(span[0] + 3*aa_num):(span[0]+3*(aa_num+1))]
                if not re.search(aa2re[aa], this_codon.upper()):
                    warnings.warn("%s (%s %d) does not correspond to %s (%s)" \
                            % (pro.id, aa, aa_num, nucl.id, this_codon))
                codon_seq += this_codon
                aa_num += 1
        return SeqRecord(CodonSeq(codon_seq, alphabet=alphabet), id=nucl.id)
    elif mode == 2:
        from collections import deque
        shift_pos = deque([])
        match = span_mode[2]
        for i in match.groupdict():
            shift_pos.append(match.span(i))
        # this rf_table is relative to nucl_seq
        rf_table = []
        i = match.start()
        while True:
            rf_table.append(i)
            i += 3
            if len(shift_pos) != 0 and i == shift_pos[0][0]:
                i = shift_pos[0][1]
                shift_pos.popleft()
            if i >= match.end():
                break
        aa_num = 0
        for aa in pro.seq:
            if aa == "-":
                codon_seq += "---"
            elif complete_protein is True and aa_num == 0:
                this_codon = nucl_seq._data[rf_table[0]:rf_table[0]+3]
                if not re.search(_codons2re[codon_table.start_codons], this_codon_upper()):
                    warnings.warn("start codon of %s (%s %d) does not correspond to %s (%s)" \
                            % (pro.id, aa, aa_num, nucl.id, this_codon))
                    codon_seq += this_codon
                    aa_num += 1
            else:
                # two types of frameshift
                if aa_num < len(pro.seq.ungap('-'))-1 and \
                        rf_table[aa_num+1]-rf_table[aa_num]-3 < 0:
                    start = rf_table[aa_num]
                    end   = rf_table[aa_num+1]-3
                    ngap  = rf_table[aa_num+1]-rf_table[aa_num]-3
                    this_codon = nucl_seq._data[start:end] + '-' * ngap
                elif rf_table[aa_num]-rf_table[aa_num-1]-3 > 0:
                    start = rf_table[aa_num-1]+3
                    end   = rf_table[aa_num]
                    ngap  = 3-(rf_table[aa_num]-rf_table[aa_num-1]-3)
                    this_codon = nucl_seq._data[start:end] + '-'*ngap + \
                            nucl_seq._data[rf_table[aa_num]:rf_table[aa_num]+3]
                else:
                    this_codon = nucl_seq._data[rf_table[aa_num]:rf_table[aa_num]+3]
                codon_seq += this_codon
                aa_num += 1
        #print len(codon_seq)
        return SeqRecord(CodonSeq(codon_seq, alphabet=alphabet, \
                rf_table=rf_table), id=nucl.id)


def _align_shift_recs(recs):
    """This function is useful to build alignment according to the
    frameshift detected by _check_corr.

    Argument:
    - recs     - a list of SeqRecords containing a CodonSeq dictated
                 by a rf_table (with frameshift in some of them).
    """
    # first check the full_rf_table of the recs are the same
    frt_lst = [rec.seq.get_full_rf_table() for rec in recs]
    if len(set([len(i) for i in frt_lst])) != 1:
        raise RuntimeError("full_rf_table of given records are not the same!!")
    else:
        frt_len = len(frt_lst[0])
    curate = [0 for i in range(len(frt_lst))]
    for i in range(frt_len-1):
        codons_range = [k[i+1]-k[i] for k in frt_lst]
        if codons_range.count(3) == len(codons_range):
            pass
        elif set(codons_range) == set([3, 6]):
            for n, j in enumerate(codons_range):
                if j == 3:
                    gaps = '-'*3
                    insert_pos = frt_lst[n][i]+3
                    seq = recs[n].seq
                    seq = CodonSeq(recs[n].seq[:insert_pos] + gaps + \
                            recs[n].seq[insert_pos:], rf_table=seq.rf_table)
                    recs[n].seq = seq
    return recs


def build(pro_align, nucl_seqs, corr_dict=None, gap_char='-', unknown='X', \
        codon_table=default_codon_table, alphabet=None, \
        complete_protein=False):
    """Build a codon alignment from a protein alignment and corresponding
    nucleotide sequences

    Arguments:
     - pro_align  - a MultipleSeqAlignment object that stores protein alignment
     - nucl_align - an object returned by SeqIO.parse or SeqIO.index or a 
                    colloction of SeqRecord.
     - alphabet   - alphabet for the returned codon alignment
     - corr_dict  - a dict that maps protein id to nucleotide id
     - complete_protein - whether the sequence begins with a start codon
                TODO: check stop codon
     - frameshift - whether to appply frameshift detection

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
    shift = None
    for pair in pro_nucl_pair:
        # Beaware that the following span corresponds to a ungapped 
        # nucleotide sequence.
        corr_span = _check_corr(pair[0], pair[1], gap_char=gap_char, \
                codon_table=codon_table, complete_protein=complete_protein)
        if not corr_span:
            raise ValueError("Protein Record %s and Nucleotide Record %s do not match!" \
                    % (pair[0].id, pair[1].id))
        else:
            codon_rec = _get_codon_rec(pair[0], pair[1], corr_span, \
                    alphabet=alphabet, complete_protein=False)
            codon_aln.append(codon_rec)
            if corr_span[1] == 2:
                shift = True
    if shift is True:
        return CodonAlignment(_align_shift_recs(codon_aln), alphabet=alphabet)
    else:
        return CodonAlignment(codon_aln, alphabet=alphabet)


def toCodonAlignment(align, alphabet=default_codon_alphabet):
    """Function to convert a MultipleSeqAlignment to CodonAlignment.
    It is the user's responsibility to ensure all the requirement
    needed by CodonAlignment is met.

    """
    rec = [SeqRecord(CodonSeq(str(i.seq), alphabet=alphabet), id=i.id) \
            for i in align._records]
    return CodonAlignment(rec, alphabet=align._alphabet)


if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()

