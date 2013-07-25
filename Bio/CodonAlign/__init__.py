# Copyright 2013 by Zheng Ruan.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Code for dealing with Codon Alignment.

"""
__docformat__ = "epytext en"  # Don't just use plain text in epydoc API pages!

from itertools import izip
from Bio.SeqRecord import SeqRecord

from CodonSeq import CodonSeq
from CodonAlignment import CodonAlignment, toCodonAlignment
from CodonAlphabet import default_codon_table, default_codon_alphabet, CodonAlphabet
from CodonAlphabet import get_codon_alphabet as _get_codon_alphabet


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
    
    >>> from Bio.Alphabet import IUPAC
    >>> from Bio.Seq import Seq
    >>> from Bio.Align import MultipleSeqAlignment
    >>> seq1 = SeqRecord(Seq('TCAGGGACTGCGAGAACCAAGCTACTGCTGCTGCTGGCTGCGCTCTGCGCCGCAGGTGGGGCGCTGGAG', 
    ...    alphabet=IUPAC.IUPACUnambiguousDNA()), id='pro1')
    >>> seq2 = SeqRecord(Seq('TCAGGGACTTCGAGAACCAAGCGCTCCTGCTGCTGGCTGCGCTCGGCGCCGCAGGTGGAGCACTGGAG', 
    ...    alphabet=IUPAC.IUPACUnambiguousDNA()), id='pro2')
    >>> pro1 = SeqRecord(Seq('SGTARTKLLLLLAALCAAGGALE', alphabet=IUPAC.protein),id='pro1')
    >>> pro2 = SeqRecord(Seq('SGTSRTKRLLLLAALGAAGGALE', alphabet=IUPAC.protein),id='pro2')
    >>> aln = MultipleSeqAlignment([pro1, pro2])
    >>> codon_aln = build(aln, [seq1, seq2])
    >>> print codon_aln
    CodonAlphabet() CodonAlignment with 2 rows and 69 columns (23 codons)
    TCAGGGACTGCGAGAACCAAGCTACTGCTGCTGCTGGCTGCGCTCTGCGCCGCAGGT...GAG pro1
    TCAGGGACTTCGAGAACCAAGCG-CTCCTGCTGCTGGCTGCGCTCGGCGCCGCAGGT...GAG pro2

    """
    # TODO
    # add an option to allow the user to specify the returned object?

    from Bio.Alphabet import ProteinAlphabet
    from Bio.Align import MultipleSeqAlignment
    
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
            alphabet = _get_codon_alphabet(codon_table, gap_char=gap_char)
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

    if not all([isinstance(pro, SeqRecord), isinstance(nucl, SeqRecord)]):
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
        print pro_re[0]
        full_pro_re = "".join(pro_re)
        match = re.search(full_pro_re, nucl_seq)
        if match:
            # mode = 1, mismatch
            return (match.span(), 1)
        else:
            # check frames of anchors
            # ten frameshift events are allowed in a sequence
            first_anchor = True
            shift_id_pos = 0
            for i in range(len(anchor_pos)-1):
                # TODO: think about last anchor (mismatch)
                if first_anchor is True and anchor_pos[0][2] != 0:
                    shift_val_lst = (1,2,anchor_len-2,anchor_len-1)
                    sh_anc = anchors[0]
                    for shift_val in shift_val_lst:
                        if shift_val in (1,2):
                            sh_nuc_len = anchor_len*3+shift_val
                        elif shift_val in (anchor_len-2,anchor_len-1):
                            sh_nuc_len = anchor_len*3-(anchor_len-shift_val)
                        if anchor_pos[0][0] >= sh_nuc_len:
                            sh_nuc = nucl_seq[anchor_pos[0][0]-sh_nuc_len:anchor_pos[0][0]]
                        else:
                            #this is unlikely to produce the correct output
                            sh_nuc = nucl_seq[:anchor_pos[0][0]]
                        qcodon = None
                        qcodon, shift_id_pos = _get_shift_anchor_re(sh_anc, sh_nuc, \
                                shift_val, aa2re, anchor_len, shift_id_pos)
                        if qcodon is not None and qcodon != -1:
                            # pro_re[0] should be '.'*anchor_len, therefore I replace it.
                            pro_re[0] = qcodon
                            break
                    if qcodon == -1:
                        import warnings
                        warnings.warn("frameshift detection failed for %s" % nucl.id)
                first_anchor = False
                shift_val = (anchor_pos[i+1][0]-anchor_pos[i][0]) % anchor_len
                sh_anc = "".join(anchors[anchor_pos[i][2]:anchor_pos[i+1][2]])
                sh_nuc = nucl_seq[anchor_pos[i][0]:anchor_pos[i+1][0]]
                qcodon = None
                if shift_val != 0:
                    qcodon, shift_id_pos = _get_shift_anchor_re(sh_anc, sh_nuc, \
                            shift_val, aa2re, anchor_len, shift_id_pos)
                if qcodon is not None and qcodon != -1:
                    pro_re[anchor_pos[i][2]:anchor_pos[i+1][2]] = [qcodon]
                    qcodon = None
                elif qcodon == -1:
                    import warnings
                    warnings.warn("frameshift detection failed for %s" % nucl.id)
            # try global match
            full_pro_re = "".join(pro_re)
            match = re.search(full_pro_re, nucl_seq)
            if match:
                return (match.span(), 2, match)
            else:
                raise RuntimeError("Protein SeqRecord (%s) and Nucleotide SeqRecord (%s) do not match!" \
                        % (pro.id, nucl.id))


def _get_shift_anchor_re(sh_anc, sh_nuc, shift_val, aa2re, \
        anchor_len, shift_id_pos):
    """This function tries all the best to come up with an re that
    matches a potentially shifted anchor.

    Arguments:
        - sh_anc    - shifted anchor sequence
        - sh_nuc    - potentially corresponding nucleotide sequence
                      of sh_anc
        - shift_val - 1 or 2 indicates forward frame shift, whereas
                      anchor_len-1 or anchor_len-2 indicates backward
                      shift
        - aa2re     - aa to codon re dict
        - anchor_len - length of the anchor
        - shift_id_pos - specify current shift name we are at
    """
    import re
    shift_id = [chr(i) for i in range(97,107)]
    if shift_val in (1, 2):
#        # obtain shifted anchor and corresponding nucl
#        # re substring matching doesn't allow number as id
#        # use ascii instead
#        qcodon = "^(?P<a>.*)"
#        # this stores the id of subre match, at least 'a' should be there
#        id_dict = {'a': 0}
#        # append subre between each codon
#        for j, aa in enumerate(sh_anc):
#            qcodon += aa2re[aa] + "(?P<" + chr(j+98) + ">.*)"
#            id_dict[chr(j+98)] = j+1
#        qcodon += "$"
#        match = re.search(qcodon, sh_nuc)
#        if match:
#            # where shift events happend (in the anchor)
#            anc_shift_pos = [num for id, num in id_dict.iteritems() \
#                    if match.group(id) != ""]
#            if 0 in anc_shift_pos:
#                qcodon = "(?P<" + shift_id[shift_id_pos] + ">.*)"
#                shift_id_pos += 1
#            else:
#                qcodon = ""
#            for j, aa in enumerate(sh_anc):
#                qcodon +=  aa2re[aa]
#                if j+1 in anc_shift_pos:
#                    qcodon += "(?P<" + shift_id[shift_id_pos] + ">.*)"
#                    shift_id_pos += 1
#            return qcodon, shift_id_pos
#        else:
#            # failed to find a match (frameshift)
#            return -1, shift_id_pos
        for j in range(len(sh_anc)):
            qcodon = "^"
            for k, aa in enumerate(sh_anc):
                if k == j:
                    qcodon += "(?P<" + shift_id[shift_id_pos] + ">.*)" \
                            + aa2re[aa]
                else:
                    qcodon += aa2re[aa]
            qcodon += "$"
            match = re.search(qcodon, sh_nuc)
            if match:
                qcodon = qcodon.replace('^','').replace('$','')
                shift_id_pos += 1
                return qcodon, shift_id_pos
        if not match:
            # failed to find a match (frameshift)
            return -1, shift_id_pos
    elif shift_val in (anchor_len-1, anchor_len-2):
        shift_val = anchor_len-shift_val
        # obtain shifted anchor and corresponding nucl
        # first check if the shifted pos is just at the end of the 
        # previous anchor.
        for j in range(1,len(sh_anc)):
            qcodon = "^"
            for k, aa in enumerate(sh_anc):
                if k == j-1:
                    # will be considered in the next step
                    pass
                elif k == j:
                    qcodon += _merge_aa2re(sh_anc[j-1], sh_anc[j], \
                            shift_val, aa2re, shift_id[shift_id_pos].upper())
                else:
                    qcodon += aa2re[aa]
            qcodon += '$'
            match = re.search(qcodon, sh_nuc)
            if match:
                qcodon = qcodon.replace('^','').replace('$','')
                shift_id_pos += 1
                return qcodon, shift_id_pos
        if not match:
            # failed to find a match (frameshift)
            return -1, shift_id_pos


def _merge_aa2re(aa1, aa2, shift_val, aa2re, reid):
    """Function to merge two amino acids based on detected frame shift
    value.
    """
    def get_aa_from_codonre(re_aa):
        aas = []
        m = 0
        for i in re_aa:
            if i == '[':
                m = -1
                aas.append('')
            elif i == ']':
                m = 0
                continue
            elif m == -1:
                aas[-1] = aas[-1] + i
            elif m == 0:
                aas.append(i)
        return aas
    scodon = map(get_aa_from_codonre, (aa2re[aa1], aa2re[aa2]))
    if shift_val == 1:
        intersect = ''.join(set(scodon[0][2]) & set(scodon[1][0]))
        scodonre = '(?P<' + reid + '>'
        scodonre += '[' + scodon[0][0] + ']' + \
                    '[' + scodon[0][1] + ']' + \
                    '[' + intersect    + ']' + \
                    '[' + scodon[1][1] + ']' + \
                    '[' + scodon[1][2] + ']'
    elif shift_val == 2:
        intersect1 = ''.join(set(scodon[0][1]) & set(scodon[1][0]))
        intersect2 = ''.join(set(scodon[0][2]) & set(scodon[1][1]))
        scodonre = '(?P<' + reid + '>'
        scodonre += '[' + scodon[0][0] + ']' + \
                    '[' + intersect1   + ']' + \
                    '[' + intersect2   + ']' + \
                    '[' + scodon[1][2] + ']'
    scodonre += ')'
    return scodonre


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
                if not re.search(_codons2re[codon_table.start_codons], this_codon.upper()):
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
        shift_start = []
        match = span_mode[2]
        m_groupdict = match.groupdict().keys()
        #if all([i.isupper() for i in m_groupdict]):
        # backward frameshift
        for i in m_groupdict:
            shift_pos.append(match.span(i))
            shift_start.append(match.start(i))
        rf_table = []
        i = match.start()
        while True:
            rf_table.append(i)
            i += 3
            #if len(shift_pos) != 0 and i in shift_start:
            if i in shift_start and m_groupdict[shift_start.index(i)].isupper():
                shift_index = shift_start.index(i)
                shift_val = 6-(shift_pos[shift_index][1]-shift_pos[shift_index][0])
                rf_table.append(i)
                rf_table.append(i+3-shift_val)
                i = shift_pos[shift_index][1]
            elif i in shift_start and m_groupdict[shift_start.index(i)].islower():
                i = shift_pos[shift_start.index(i)][1]
            if i >= match.end():
                break
        aa_num = 0
        for aa in pro.seq:
            if aa == "-":
                codon_seq += "---"
            elif complete_protein is True and aa_num == 0:
                this_codon = nucl_seq._data[rf_table[0]:rf_table[0]+3]
                if not re.search(_codons2re[codon_table.start_codons], this_codon.upper()):
                    warnings.warn("start codon of %s (%s %d) does not correspond to %s (%s)" \
                            % (pro.id, aa, aa_num, nucl.id, this_codon))
                    codon_seq += this_codon
                    aa_num += 1
            else:
                if aa_num < len(pro.seq.ungap('-'))-1 and \
                        rf_table[aa_num+1]-rf_table[aa_num]-3 < 0:
                    start = rf_table[aa_num]
                    end   = start + (3-shift_val)
                    ngap  = shift_val
                    this_codon = nucl_seq._data[start:end] + '-'*ngap
                elif rf_table[aa_num]-rf_table[aa_num-1]-3 > 0:
                    start = rf_table[aa_num-1]+3
                    end   = rf_table[aa_num]
                    ngap  = 3-(rf_table[aa_num]-rf_table[aa_num-1]-3)
                    this_codon = nucl_seq._data[start:end] + '-'*ngap + \
                            nucl_seq._data[rf_table[aa_num]:rf_table[aa_num]+3]
                else:
                    start = rf_table[aa_num]
                    end   = start + 3
                    this_codon = nucl_seq._data[start:end]
                codon_seq += this_codon
                aa_num += 1
        return SeqRecord(CodonSeq(codon_seq, alphabet=alphabet, \
                    rf_table=rf_table), id=nucl.id)

#        else:
#            #for i in match.groupdict():
#            for i in m_groupdict:
#                if i.islower(): shift_pos.append(match.span(i))
#            # this rf_table is relative to nucl_seq
#            rf_table = []
#            i = match.start()
#            while True:
#                rf_table.append(i)
#                i += 3
#                if len(shift_pos) != 0 and i == shift_pos[0][0]:
#                    i = shift_pos[0][1]
#                    shift_pos.popleft()
#                # TODO:
#                # if i > shift_pos[0][0], then an exception should raise up.
#                if i >= match.end():
#                    break
#            aa_num = 0
#            for aa in pro.seq:
#                if aa == "-":
#                    codon_seq += "---"
#                elif complete_protein is True and aa_num == 0:
#                    this_codon = nucl_seq._data[rf_table[0]:rf_table[0]+3]
#                    if not re.search(_codons2re[codon_table.start_codons], this_codon.upper()):
#                        warnings.warn("start codon of %s (%s %d) does not correspond to %s (%s)" \
#                                % (pro.id, aa, aa_num, nucl.id, this_codon))
#                        codon_seq += this_codon
#                        aa_num += 1
#                else:
#                    # two types of frameshift
#                    if aa_num < len(pro.seq.ungap('-'))-1 and \
#                            rf_table[aa_num+1]-rf_table[aa_num]-3 < 0:
#                        start = rf_table[aa_num]
#                        end   = rf_table[aa_num+1]-3
#                        ngap  = rf_table[aa_num+1]-rf_table[aa_num]-3
#                        this_codon = nucl_seq._data[start:end] + '-' * ngap
#                    elif rf_table[aa_num]-rf_table[aa_num-1]-3 > 0:
#                        start = rf_table[aa_num-1]+3
#                        end   = rf_table[aa_num]
#                        ngap  = 3-(rf_table[aa_num]-rf_table[aa_num-1]-3)
#                        this_codon = nucl_seq._data[start:end] + '-'*ngap + \
#                                nucl_seq._data[rf_table[aa_num]:rf_table[aa_num]+3]
#                    else:
#                        this_codon = nucl_seq._data[rf_table[aa_num]:rf_table[aa_num]+3]
#                    codon_seq += this_codon
#                    aa_num += 1
#            return SeqRecord(CodonSeq(codon_seq, alphabet=alphabet, \
#                    rf_table=rf_table), id=nucl.id)
#

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
                    seq = CodonSeq(recs[n].seq._data[:insert_pos] + gaps + \
                            recs[n].seq._data[insert_pos:], rf_table=seq.rf_table)
                    recs[n].seq = seq
    return recs


if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()

