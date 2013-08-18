# Copyright 2013 by Zheng Ruan.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Code for dealing with Codon Seq.

CodonSeq class is interited from Seq class. This is the core class to
deal with sequences in CodonAlignment in biopython.

"""
from __future__ import division

__docformat__ = "epytext en"  # Don't just use plain text in epydoc API pages!

from math import log

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
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
                             + "Plase use normal slice option.")
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
                elif nxt_shift > 0:
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


def _get_codon_list(codonseq):
    """get a list of codons according to full_rf_table for counting (PRIVATE).
    """
    if not isinstance(codonseq, CodonSeq):
        raise TypeError("_get_codon_list accept a CodonSeq object (%s detected)" \
                % type(codonseq))
    full_rf_table = codonseq.get_full_rf_table()
    codon_lst = []
    for i in range(len(full_rf_table)):
        start = full_rf_table[i]
        try:
            end = full_rf_table[i+1]
        except IndexError:
            end = start+3
        this_codon = str(codonseq[start:end])
        if len(this_codon) == 3:
            codon_lst.append(this_codon)
        else:
            codon_lst.append(str(this_codon.ungap()))
    return codon_lst


def cal_dn_ds(codon_seq1, codon_seq2, method="NG86", \
        codon_table=default_codon_table, w=1):
    """Function to calculate the dN and dS of the given two CodonSeq
    or SeqRecord that contain CodonSeq objects.

    Available methods:
        - NG86  - PMID: 3444411
        - LWL85 - PMID: 3916709
        - ML    - PMID: 7968486
    
    Arguments:
        - w  - transition/transvertion ratio
    """
    if all([isinstance(codon_seq1, CodonSeq), isinstance(codon_seq2, CodonSeq)]):
        pass
    elif all([isinstance(codon_seq1, SeqRecord), isinstance(codon_seq2, SeqRecord)]):
        assert isinstance(codon_seq1.seq, CodonSeq), "cal_dn_ds accepts SeqRecords that contain CodonSeq as its seq!"
        assert isinstance(codon_seq2.seq, CodonSeq), "cal_dn_ds accepts SeqRecords that contain CodonSeq as its seq!"
        codon_seq1 = codon_seq1.seq
        codon_seq2 = codon_seq2.seq
    else:
        raise TypeError("cal_dn_ds accepts two CodonSeq objects or SeqRecord that contains CodonSeq as its seq!")
    if len(codon_seq1.get_full_rf_table()) != len(codon_seq2.get_full_rf_table()):
        raise RuntimeError("full_rf_table length of seq1 (%d) and seq2 (%s) are not the same" \
                % (len(codon_seq1.get_full_rf_table()), len(codon_seq2.get_full_rf_table())))
    seq1_codon_lst = _get_codon_list(codon_seq1)
    seq2_codon_lst = _get_codon_list(codon_seq2)
    if method == "NG86":
        S_sites1, N_sites1 = _count_site(seq1_codon_lst, codon_table=codon_table, w=w)
        S_sites2, N_sites2 = _count_site(seq2_codon_lst, codon_table=codon_table, w=w)
        S_sites = (S_sites1 + S_sites2) / 2.0
        N_sites = (N_sites1 + N_sites2) / 2.0
        NS = [0, 0]
        for i, j in zip(seq1_codon_lst, seq2_codon_lst):
            NS = [m+n for m,n in zip(NS, _count_diff(i, j, codon_table=codon_table))]
        ps = NS[0] / S_sites
        pn = NS[1] / N_sites
        dN = -3.0/4*log(1-4.0/3*pn)
        dS = -3.0/4*log(1-4.0/3*ps)
        return dN, dS
    elif method == "LWL85":
        # Nomenclature is according to PMID (3916709)
        codon_fold_dict = _get_codon_fold(codon_table)
        # count number of sites in different degenerate classes
        fold0 = [0, 0]
        fold2 = [0, 0]
        fold4 = [0, 0]
        for codon in seq1_codon_lst:
            fold_num = codon_fold_dict[codon]
            for f in fold_num:
                if f == '0':
                    fold0[0] += 1
                elif f == '2':
                    fold2[0] += 1
                elif f == '4':
                    fold4[0] += 1
        for codon in seq2_codon_lst:
            fold_num = codon_fold_dict[codon]
            for f in fold_num:
                if f == '0':
                    fold0[1] += 1
                elif f == '2':
                    fold2[1] += 1
                elif f == '4':
                    fold4[1] += 1
        L = [sum(fold0)/2.0, sum(fold2)/2.0, sum(fold4)/2.0]
        # count number of differences in different degenerate classes
        PQ = [0] * 6 # with P0, P2, P4, Q0, Q2, Q4 in each position
        for codon1, codon2 in zip(seq1_codon_lst, seq2_codon_lst):
            if (codon1 == "---" or codon2 == "---") or codon1 == codon2:
                continue
            else:
                PQ = [i+j for i, j in zip(PQ, _diff_codon(codon1, codon2, fold_dict=codon_fold_dict))]
        PQ = [i/j for i, j in zip(PQ, L*2)]
        P = PQ[:3]
        Q = PQ[3:]
        A = [(1./2)*log(1./(1-2*i-j)) - (1./4)*log(1./(1-2*j)) for i, j in zip(P, Q)]
        B = [(1./2)*log(1./(1-2*i)) for i in Q]
        dS = 3*(L[2]*A[1]+L[2]*(A[2]+B[2]))/(L[1]+3*L[2])
        dN = 3*(L[2]*B[1]+L[0]*(A[0]+B[0]))/(2*L[1]+3*L[0])
        return dN, dS
    elif method == "ML":
        from collections import Counter
        from scipy.optimize import minimize
        codon_cnt = Counter()
        # three codon position
        fcodon = [{'A': 0, 'G': 0, 'C': 0, 'T': 0},
                  {'A': 0, 'G': 0, 'C': 0, 'T': 0},
                  {'A': 0, 'G': 0, 'C': 0, 'T': 0}]
        for i in seq1_codon_lst + seq2_codon_lst:
            if i != '---':
                fcodon[0][i[0]] += 1
                fcodon[1][i[1]] += 1
                fcodon[2][i[2]] += 1
        for i in range(3):
            tot = sum(fcodon[i].values())
            fcodon[i] = {j: k/tot for j, k in fcodon[i].items()}
        pi = {}
        for i in set(seq1_codon_lst+seq2_codon_lst):
            if i != '---':
                pi[i] = fcodon[0][i[0]]*fcodon[1][i[1]]*fcodon[2][i[2]]
        for i, j in zip(seq1_codon_lst, seq2_codon_lst):
            #if i != j and ('---' not in (i, j)):
            if '---' not in (i, j):
                codon_cnt[(i,j)] += 1
        codon_lst = [i for i in \
                codon_table.forward_table.keys() + codon_table.stop_codons if 'U' not in i]
        # apply optimization
        def func(params, pi=pi, codon_cnt=codon_cnt, codon_lst=codon_lst, \
                codon_table=codon_table):
            """params = [t, k, w]"""
            return -_likelihood_func(params[0], params[1], params[2], pi, \
                    codon_cnt, codon_lst=codon_lst, codon_table=codon_table)
        # count sites
        opt_res = minimize(func, [1, 0.1, 2], method='L-BFGS-B', \
                bounds=((1e-10, 20), (1e-10, 20), (1e-10, 10)), tol=1e-5)
        t, k, w = opt_res.x
        Q = _get_Q(pi, k, w, codon_lst, codon_table)
        Sd = Nd = 0
        for i, c1 in enumerate(codon_lst):
            for j, c2 in enumerate(codon_lst):
                if i != j:
                    try:
                        if codon_table.forward_table[c1] == codon_table.forward_table[c2]:
                            # synonymous count
                            Sd += pi[c1] * Q[i, j]
                        else:
                            # nonsynonymous count
                            Nd += pi[c1] * Q[i, j]
                    except:
                        # This is probably due to stop codons
                        pass
        Sd *= t
        Nd *= t
        # count differences (with w fixed to 1)
        opt_res = minimize(func, [1, 0.1, 2], method='L-BFGS-B', \
                bounds=((1e-10, 20), (1e-10, 20), (1, 1)), tol=1e-5)
        t, k, w = opt_res.x
        Q = _get_Q(pi, k, w, codon_lst, codon_table)
        rhoS = rhoN = 0
        for i, c1 in enumerate(codon_lst):
            for j, c2 in enumerate(codon_lst):
                if i != j:
                    try:
                        if codon_table.forward_table[c1] == codon_table.forward_table[c2]:
                            # synonymous count
                            rhoS += pi[c1] * Q[i, j]
                        else:
                            # nonsynonymous count
                            rhoN += pi[c1] * Q[i, j]
                    except:
                        # This is probably due to stop codons
                        pass
        rhoS *= 3
        rhoN *= 3
        dN = Nd/rhoN
        dS = Sd/rhoS
        return dN, dS



#################################################################
#  private functions for NG86 method
#################################################################
def _count_site(codon_lst, w=1, codon_table=default_codon_table):
    """count synonymous and non-synonymous sites of a list of codons (PRIVATE).
    Argument:
        - codon_lst - A three letter codon list from a CodonSeq object. This
                      can be returned from _get_codon_list method.
        - w         - transition/transversion rate ratio
    """
    S_site = 0 # synonymous sites
    N_site = 0 # non-synonymous sites
    purine     = ('A', 'G')
    pyrimidine = ('T', 'C')
    base_tuple = ('A', 'T', 'C', 'G')
    for codon in codon_lst:
        neighbor_codon = {'transition': [], 'transversion': []}
        # classify neighbor codons
        codon = codon.replace('U', 'T')
        if codon == '---': continue
        for n, i in enumerate(codon):
            for j in base_tuple:
                if  i == j:
                    pass
                elif i in purine and j in purine:
                    codon_chars = [c for c in codon]
                    codon_chars[n] = j
                    this_codon = ''.join(codon_chars)
                    neighbor_codon['transition'].append(this_codon)
                elif i in pyrimidine and j in pyrimidine:
                    codon_chars = [c for c in codon]
                    codon_chars[n] = j
                    this_codon = ''.join(codon_chars)
                    neighbor_codon['transition'].append(this_codon)
                else: 
                    codon_chars = [c for c in codon]
                    codon_chars[n] = j
                    this_codon = ''.join(codon_chars)
                    neighbor_codon['transversion'].append(this_codon)
        # count synonymous and non-synonymous sites
        aa = codon_table.forward_table[codon]
        for neighbor in neighbor_codon['transition']:
            if neighbor in codon_table.stop_codons:
                N_site += 1
            elif codon_table.forward_table[neighbor] == aa:
                S_site += 1
            else:
                N_site += 1
        for neighbor in neighbor_codon['transversion']:
            if neighbor in codon_table.stop_codons:
                N_site += w
            elif codon_table.forward_table[neighbor] == aa:
                S_site += w
            else:
                N_site += w
    return (S_site/3.0, N_site/3.0)


def _count_diff(codon1, codon2, codon_table=default_codon_table):
    """Count differences between two codons (three-letter string).
    The function will take multiple pathways from codon1 to codon2
    into account (PRIVATE).
    """
    if not all([isinstance(codon1, str), isinstance(codon2, str)]):
        raise TypeError("_count_diff accept string object to represent codon (%s, %s detected)" \
                % (type(codon1), type(codon2)))
    if len(codon1) != 3 or len(codon2) != 3:
        raise RuntimeError("codon should be three letter string (%d, %d detected)" \
                (len(codon1), len(codon2)))
    SN = [0, 0]
    Sd = 0 # synonymous differences
    Nd = 0 # non-synonymous differences
    if codon1 == '---' or codon2 == '---':
        return SN
    base_tuple = ('A', 'C', 'G', 'T')
    if not all([i in base_tuple for i in codon1]):
        raise RuntimeError("Unrecognized character detected in codon1 %s (Codon are consists of A, T, C or G)" % codon1)
    if not all([i in base_tuple for i in codon2]):
        raise RuntimeError("Unrecognized character detected in codon2 %s (Codon are consists of A, T, C or G)" % codon2)
    if codon1 == codon2:
        return SN
    else:
        diff_pos = []
        for i, k in enumerate(zip(codon1, codon2)):
            if k[0] != k[1]:
                diff_pos.append(i)
        def compare_codon(codon1, codon2, codon_table=default_codon_table, weight=1):
            """Method to compare two codon accounting for different pathways"""
            sd = nd = 0
            if len(set(map(codon_table.forward_table.get, [codon1, codon2]))) == 1:
                sd += weight
            else:
                nd += weight
            return (sd, nd)
        if len(diff_pos) == 1:
            SN = [i+j for i,j in zip(SN, \
                    compare_codon(codon1, codon2, codon_table=codon_table))]
        elif len(diff_pos) == 2:
            codon2_aa = codon_table.forward_table[codon2]
            for i in diff_pos:
                codon1_chars = [c for c in codon1]
                codon1_chars[i] = codon2[i]
                temp_codon = ''.join(codon1_chars)
                SN = [i+j for i,j in zip(SN, \
                        compare_codon(codon1, temp_codon, codon_table=codon_table, weight=0.5))]
                SN = [i+j for i,j in zip(SN, \
                        compare_codon(temp_codon, codon2, codon_table=codon_table, weight=0.5))]
        elif len(diff_pos) == 3:
            # we are now in the most complex situation
            # I don't want to think about this.
            # the substitution is considered non-synonymous (modify!!!)
            SN[1] += 1
    return SN
        

#################################################################
#  private functions for LWL85 method
#################################################################
def _get_codon_fold(codon_table):
    """function to classify different position in a codon into
    different fold (PRIVATE).
    """
    def find_fold_class(codon, forward_table):
        base = set(['A', 'T', 'C', 'G'])
        fold = ''
        codon_base_lst = [i for i in codon]
        for i, b in enumerate(codon_base_lst):
            other_base = base - set(b)
            aa = []
            for j in other_base:
                codon_base_lst[i] = j
                try:
                    aa.append(forward_table[''.join(codon_base_lst)])
                except KeyError:
                    aa.append('stop')
            if aa.count(forward_table[codon]) == 0:
                fold += '0'
            elif aa.count(forward_table[codon]) in (1,2):
                fold += '2'
            elif aa.count(forward_table[codon]) == 3:
                fold += '4'
            else:
                raise RuntimeError("Unknown Error, cannot assign the position to a fold")
            codon_base_lst[i] = b
        return fold
    fold_table = {}
    for codon in codon_table.forward_table:
        if 'U' not in codon:
            fold_table[codon] = find_fold_class(codon, \
                    codon_table.forward_table)
    fold_table["---"] = '000'
    return fold_table


def _diff_codon(codon1, codon2, fold_dict):
    """function to get the differences of two codon and return
    number of different types substitutions

    return (P0, P2, P4, Q0, Q2, Q4)
    Nomenclature is according to PMID (3916709)
    """
    P0 = P2 = P4 = Q0 = Q2 = Q4 = 0
    fold_num = fold_dict[codon1]
    purine = ('A', 'G')
    pyrimidine = ('T', 'C')
    for n, (i, j) in enumerate(zip(codon1, codon2)):
            if i!= j and (i in purine and j in purine):
                if fold_num[n] == '0':
                    P0 += 1
                elif fold_num[n] == '2':
                    P2 += 1
                elif fold_num[n] == '4':
                    P4 += 1
                else:
                    raise RuntimeError("Unexpected fold_num %d" % fold_num[n])
            if i!= j and (i in pyrimidine and j in pyrimidine):
                if fold_num[n] == '0':
                    P0 += 1
                elif fold_num[n] == '2':
                    P2 += 1
                elif fold_num[n] == '4':
                    P4 += 1
                else:
                    raise RuntimeError("Unexpected fold_num %d" % fold_num[n])
            if i != j and ((i in purine and j in pyrimidine) \
                    or (i in pyrimidine and j in purine)):
                if fold_num[n] == '0':
                    Q0 += 1
                elif fold_num[n] == '2':
                    Q2 += 1
                elif fold_num[n] == '4':
                    Q4 += 1
                else:
                    raise RuntimeError("Unexpected fold_num %d" % fold_num[n])
    return (P0, P2, P4, Q0, Q2, Q4)


#################################################################
#  private functions for Maximum Likelihood method
#################################################################

def _q(i, j, pi, k, w, codon_table=default_codon_table):
    """Q matrix for codon substitution.

    Arguments:
        - i, j  : three letter codon string
        - pi    : expected codon frequency
        - k     : transition/transversion ratio
        - w     : nonsynonymous/synonymous rate ratio
        - codon_table: Bio.Data.CodonTable object
    """
    if i == j:
        # diagonal elements is the sum of all other elements
        return 0
    if i in codon_table.stop_codons or j in codon_table.stop_codons:
        return 0
    if (i not in pi) or (j not in pi):
        return 0
    purine = ('A', 'G')
    pyrimidine = ('T', 'C')
    diff = []
    for n, (c1, c2) in enumerate(zip(i, j)):
        if c1 != c2:
            diff.append((n, c1, c2))
    if len(diff) >= 2:
        return 0
    if codon_table.forward_table[i] == codon_table.forward_table[j]:
        # synonymous substitution
        if diff[0][1] in purine and diff[0][2] in purine:
            # transition
            return k*pi[j]
        elif diff[0][1] in pyrimidine and diff[0][2] in pyrimidine:
            # transition
            return k*pi[j]
        else:
            # transversion
            return pi[j]
    else:
        # nonsynonymous substitution
        if diff[0][1] in purine and diff[0][2] in purine:
            # transition
            return w*k*pi[j]
        elif diff[0][1] in pyrimidine and diff[0][2] in pyrimidine:
            # transition
            return w*k*pi[j]
        else:
            # transversion
            return w*pi[j]

def _get_Q(pi, k, w, codon_lst, codon_table):
    """Q matrix for codon substitution"""
    import numpy as np
    codon_num = len(codon_lst)
    Q = np.zeros((codon_num, codon_num))
    for i in range(codon_num):
        for j in range(codon_num):
            if i != j:
                Q[i, j] = _q(codon_lst[i], codon_lst[j], pi, k, w, \
                        codon_table=codon_table)
    nucl_substitutions = 0
    for i in range(codon_num):
        Q[i,i] = -sum(Q[i,:])
        try:
            nucl_substitutions += pi[codon_lst[i]] * (-Q[i, i])
        except KeyError:
            pass
    Q = Q / nucl_substitutions
    return Q


def _likelihood_func(t, k, w, pi, codon_cnt, codon_lst, codon_table):
    """likelihood function for ML method
    """
    from scipy.linalg import expm
    Q = _get_Q(pi, k, w, codon_lst, codon_table)
    P = expm(Q*t)
    l = 0 # likelihood value
    for i, c1 in enumerate(codon_lst):
        for j, c2 in enumerate(codon_lst):
            #if i == j: print P[i, j]
            if (c1, c2) in codon_cnt:
                if P[i, j] * pi[c1] <= 0:
                    #print c1, c2
                    l += codon_cnt[(c1, c2)] * 0
                else:
                    l += codon_cnt[(c1, c2)] * log(pi[c1] * P[i, j])
    return l


if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()

