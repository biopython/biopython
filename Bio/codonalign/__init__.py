# Copyright 2013 by Zheng Ruan (zruan1991@gmail.com).
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Code for dealing with Codon Alignments."""

import copy
from collections.abc import Mapping, Iterable

from Bio import BiopythonWarning
from Bio import BiopythonExperimentalWarning

from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable

from Bio.codonalign.codonseq import CodonSeq
from Bio.codonalign.codonalignment import CodonAlignment, mktest

import warnings

warnings.warn(
    "Bio.codonalign is an experimental module which may undergo "
    "significant changes prior to its future official release.",
    BiopythonExperimentalWarning,
)


def build(
    pro_align,
    nucl_seqs,
    corr_dict=None,
    gap_char="-",
    unknown="X",
    codon_table=None,
    complete_protein=False,
    anchor_len=10,
    max_score=10,
):
    """Build a codon alignment from protein alignment and corresponding nucleotides.

    Arguments:
     - pro_align  - a protein MultipleSeqAlignment object
     - nucl_seqs - an object returned by SeqIO.parse or SeqIO.index
       or a collection of SeqRecord.
     - corr_dict  - a dict that maps protein id to nucleotide id
     - complete_protein - whether the sequence begins with a start
       codon

    Return a CodonAlignment object.

    The example below answers this Biostars question: https://www.biostars.org/p/89741/

    >>> from Bio.Seq import Seq
    >>> from Bio.SeqRecord import SeqRecord
    >>> from Bio.Align import MultipleSeqAlignment
    >>> from Bio.codonalign import build
    >>> seq1 = SeqRecord(Seq('ATGTCTCGT'), id='pro1')
    >>> seq2 = SeqRecord(Seq('ATGCGT'), id='pro2')
    >>> pro1 = SeqRecord(Seq('MSR'), id='pro1')
    >>> pro2 = SeqRecord(Seq('M-R'), id='pro2')
    >>> aln = MultipleSeqAlignment([pro1, pro2])
    >>> codon_aln = build(aln, [seq1, seq2])
    >>> print(codon_aln)
    CodonAlignment with 2 rows and 9 columns (3 codons)
    ATGTCTCGT pro1
    ATG---CGT pro2

    """
    # TODO
    # add an option to allow the user to specify the returned object?

    from Bio.Align import MultipleSeqAlignment

    # check the type of object of pro_align
    if not isinstance(pro_align, MultipleSeqAlignment):
        raise TypeError("the first argument should be a MultipleSeqAlignment object")
    # check whether the number of seqs in pro_align and nucl_seqs is
    # the same
    pro_num = len(pro_align)
    if corr_dict is None:
        try:
            nucl_num = len(nucl_seqs)
        except TypeError:
            # nucl_seqs will be an iterator if returned by SeqIO.parse()
            nucl_seqs = tuple(nucl_seqs)
            nucl_num = len(nucl_seqs)
        if pro_num > nucl_num:
            raise ValueError(
                f"Higher Number of SeqRecords in Protein Alignment ({pro_num}) "
                f"than the Number of Nucleotide SeqRecords ({nucl_num}) are found!"
            )

        # Determine the protein sequences and nucl sequences
        # correspondence. If nucl_seqs is a list, tuple or read by
        # SeqIO.parse(), we assume the order of sequences in pro_align
        # and nucl_seqs are the same. If nucl_seqs is a dict or read by
        # SeqIO.index(), we match seqs in pro_align and those in
        # nucl_seq by their id.
        if isinstance(nucl_seqs, Mapping):
            corr_method = 1
        elif isinstance(nucl_seqs, Iterable):
            corr_method = 0
        else:
            raise TypeError(
                "Nucl Sequences Error, Unknown type to assign correspondence method"
            )
    else:
        if not isinstance(corr_dict, dict):
            raise TypeError(
                "corr_dict should be a dict that corresponds "
                "protein id to nucleotide id!"
            )
        if len(corr_dict) >= pro_num:
            if isinstance(nucl_seqs, Mapping):
                pass
            else:
                d = {}
                for record in nucl_seqs:
                    key = record.id
                    if key in d:
                        raise ValueError("Duplicate key '%s'" % key)
                    d[key] = record
                nucl_seqs = d
            corr_method = 2
        else:
            raise RuntimeError(
                f"Number of items in corr_dict ({len(corr_dict)}) "
                f"is less than number of protein records ({pro_num})"
            )

    # set up pro-nucl correspondence based on corr_method
    # corr_method = 0, consecutive pairing
    if corr_method == 0:
        pro_nucl_pair = zip(pro_align, nucl_seqs)
    # corr_method = 1, keyword pairing
    elif corr_method == 1:
        nucl_id = set(nucl_seqs.keys())
        pro_id = {i.id for i in pro_align}
        # check if there is pro_id that does not have a nucleotide match
        if pro_id - nucl_id:
            diff = pro_id - nucl_id
            raise ValueError(
                f"Protein Record {', '.join(diff)} cannot find a "
                "nucleotide sequence match, please check the id"
            )
        else:
            pro_nucl_pair = []
            for pro_rec in pro_align:
                pro_nucl_pair.append((pro_rec, nucl_seqs[pro_rec.id]))
    # corr_method = 2, dict pairing
    elif corr_method == 2:
        pro_nucl_pair = []
        for pro_rec in pro_align:
            try:
                nucl_id = corr_dict[pro_rec.id]
            except KeyError:
                print("Protein record (%s) is not in corr_dict!" % pro_rec.id)
                exit(1)
            pro_nucl_pair.append((pro_rec, nucl_seqs[nucl_id]))

    if codon_table is None:
        codon_table = CodonTable.generic_by_id[1]

    codon_aln = []
    shift = False
    for pair in pro_nucl_pair:
        # Beware that the following span corresponds to an ungapped
        # nucleotide sequence.
        corr_span = _check_corr(
            pair[0],
            pair[1],
            gap_char=gap_char,
            codon_table=codon_table,
            complete_protein=complete_protein,
            anchor_len=anchor_len,
        )
        if not corr_span:
            raise ValueError(
                f"Protein Record {pair[0].id} and "
                f"Nucleotide Record {pair[1].id} do not match!"
            )
        else:
            codon_rec = _get_codon_rec(
                pair[0],
                pair[1],
                corr_span,
                gap_char=gap_char,
                complete_protein=complete_protein,
                codon_table=codon_table,
                max_score=max_score,
            )
            codon_aln.append(codon_rec)
            if corr_span[1] == 2:
                shift = True
    if shift:
        return CodonAlignment(_align_shift_recs(codon_aln))
    else:
        return CodonAlignment(codon_aln)


def _codons2re(codons):
    """Generate regular expression based on a given list of codons (PRIVATE)."""
    reg = ""
    for i in zip(*codons):
        if len(set(i)) == 1:
            reg += "".join(set(i))
        else:
            reg += "[" + "".join(set(i)) + "]"
    return reg


def _get_aa_regex(codon_table, stop="*", unknown="X"):
    """Set up the regular expression of a given CodonTable (PRIVATE).

    >>> from Bio.Data.CodonTable import generic_by_id
    >>> p = generic_by_id[1]
    >>> t = _get_aa_regex(p)
    >>> print(t['A'][0])
    G
    >>> print(t['A'][1])
    C
    >>> print(sorted(list(t['A'][2:])))
    ['A', 'C', 'G', 'T', 'U', '[', ']']
    >>> print(sorted(list(t['L'][:5])))
    ['C', 'T', 'U', '[', ']']
    >>> print(sorted(list(t['L'][5:9])))
    ['T', 'U', '[', ']']
    >>> print(sorted(list(t['L'][9:])))
    ['A', 'C', 'G', 'T', 'U', '[', ']']

    """
    from Bio.Data.CodonTable import CodonTable

    if not isinstance(codon_table, CodonTable):
        raise TypeError("Input table is not a instance of Bio.Data.CodonTable object")
    aa2codon = {}
    for codon, aa in codon_table.forward_table.items():
        aa2codon.setdefault(aa, []).append(codon)
    for aa, codons in aa2codon.items():
        aa2codon[aa] = _codons2re(codons)
    aa2codon[stop] = _codons2re(codon_table.stop_codons)
    aa2codon[unknown] = "..."
    return aa2codon


def _check_corr(
    pro, nucl, gap_char, codon_table, complete_protein=False, anchor_len=10,
):
    """Check if the nucleotide can be translated into the protein (PRIVATE).

    Expects two SeqRecord objects.
    """
    import re

    if not isinstance(pro, SeqRecord) or not isinstance(nucl, SeqRecord):
        raise TypeError(
            "_check_corr accepts two SeqRecord object. Please check your input."
        )

    aa2re = _get_aa_regex(codon_table)
    pro_re = ""
    for aa in pro.seq:
        if aa != gap_char:
            pro_re += aa2re[aa]

    nucl_seq = str(nucl.seq.upper().replace(gap_char, ""))
    match = re.search(pro_re, nucl_seq)
    if match:
        # mode = 0, direct match
        return (match.span(), 0)
    else:
        # Might caused by mismatches or frameshift, using anchors to
        # have a try
        # anchor_len = 10 # adjust this value to test performance
        pro_seq = str(pro.seq).replace(gap_char, "")
        anchors = [
            pro_seq[i : (i + anchor_len)] for i in range(0, len(pro_seq), anchor_len)
        ]
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
            qcodon = ""
            fncodon = ""
            # dirty code to deal with the last anchor
            # as the last anchor is combined in the steps
            # above, we need to get the true last anchor to
            # pro_re
            if this_anchor_len == anchor_len:
                for aa in anchor:
                    if complete_protein and i == 0:
                        qcodon += _codons2re(codon_table.start_codons)
                        fncodon += aa2re["X"]
                        continue
                    qcodon += aa2re[aa]
                    fncodon += aa2re["X"]
                match = re.search(qcodon, nucl_seq)
            elif this_anchor_len > anchor_len:
                last_qcodon = ""
                last_fcodon = ""
                for j in range(anchor_len, len(anchor)):
                    last_qcodon += aa2re[anchor[j]]
                    last_fcodon += aa2re["X"]
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
            # check the first anchor
            if first_anchor and anchor_pos[0][2] != 0:
                shift_val_lst = [1, 2, 3 * anchor_len - 2, 3 * anchor_len - 1, 0]
                sh_anc = anchors[0]
                for shift_val in shift_val_lst:
                    if shift_val == 0:
                        qcodon = None
                        break
                    if shift_val in (1, 2):
                        sh_nuc_len = anchor_len * 3 + shift_val
                    elif shift_val in (3 * anchor_len - 2, 3 * anchor_len - 1):
                        sh_nuc_len = anchor_len * 3 - (3 * anchor_len - shift_val)
                    if anchor_pos[0][0] >= sh_nuc_len:
                        sh_nuc = nucl_seq[
                            anchor_pos[0][0] - sh_nuc_len : anchor_pos[0][0]
                        ]
                    else:
                        # this is unlikely to produce the correct output
                        sh_nuc = nucl_seq[: anchor_pos[0][0]]
                    qcodon, shift_id_pos = _get_shift_anchor_re(
                        sh_anc, sh_nuc, shift_val, aa2re, anchor_len, shift_id_pos
                    )
                    if qcodon is not None and qcodon != -1:
                        # pro_re[0] should be '.'*anchor_len, therefore I
                        # replace it.
                        pro_re[0] = qcodon
                        break
                if qcodon == -1:
                    warnings.warn(
                        f"first frameshift detection failed for {nucl.id}",
                        BiopythonWarning,
                    )
            # check anchors in the middle
            for i in range(len(anchor_pos) - 1):
                shift_val = (anchor_pos[i + 1][0] - anchor_pos[i][0]) % (3 * anchor_len)
                sh_anc = "".join(anchors[anchor_pos[i][2] : anchor_pos[i + 1][2]])
                sh_nuc = nucl_seq[anchor_pos[i][0] : anchor_pos[i + 1][0]]
                qcodon = None
                if shift_val != 0:
                    qcodon, shift_id_pos = _get_shift_anchor_re(
                        sh_anc, sh_nuc, shift_val, aa2re, anchor_len, shift_id_pos
                    )
                if qcodon is not None and qcodon != -1:
                    pro_re[anchor_pos[i][2] : anchor_pos[i + 1][2]] = [qcodon]
                    qcodon = None
                elif qcodon == -1:
                    warnings.warn(
                        f"middle frameshift detection failed for {nucl.id}",
                        BiopythonWarning,
                    )
            # check the last anchor
            if anchor_pos[-1][2] + 1 == len(anchors) - 1:
                sh_anc = anchors[-1]
                this_anchor_len = len(sh_anc)
                shift_val_lst = [
                    1,
                    2,
                    3 * this_anchor_len - 2,
                    3 * this_anchor_len - 1,
                    0,
                ]
                for shift_val in shift_val_lst:
                    if shift_val == 0:
                        qcodon = None
                        break
                    if shift_val in (1, 2):
                        sh_nuc_len = this_anchor_len * 3 + shift_val
                    elif shift_val in (
                        3 * this_anchor_len - 2,
                        3 * this_anchor_len - 1,
                    ):
                        sh_nuc_len = this_anchor_len * 3 - (
                            3 * this_anchor_len - shift_val
                        )
                    if len(nucl_seq) - anchor_pos[-1][0] >= sh_nuc_len:
                        sh_nuc = nucl_seq[
                            anchor_pos[-1][0] : anchor_pos[-1][0] + sh_nuc_len
                        ]
                    else:
                        # this is unlikely to produce the correct output
                        sh_nuc = nucl_seq[anchor_pos[-1][0] :]
                    qcodon, shift_id_pos = _get_shift_anchor_re(
                        sh_anc, sh_nuc, shift_val, aa2re, this_anchor_len, shift_id_pos
                    )
                    if qcodon is not None and qcodon != -1:
                        pro_re.pop()
                        pro_re[-1] = qcodon
                        break
                if qcodon == -1:
                    warnings.warn(
                        f"last frameshift detection failed for {nucl.id}",
                        BiopythonWarning,
                    )
            # try global match
            full_pro_re = "".join(pro_re)
            match = re.search(full_pro_re, nucl_seq)
            if match:
                return (match.span(), 2, match)
            else:
                raise RuntimeError(
                    f"Protein SeqRecord ({pro.id}) and "
                    f"Nucleotide SeqRecord ({nucl.id}) do not match!"
                )


def _get_shift_anchor_re(sh_anc, sh_nuc, shift_val, aa2re, anchor_len, shift_id_pos):
    """Find a regular expression matching a potentially shifted anchor (PRIVATE).

    Arguments:
     - sh_anc    - shifted anchor sequence
     - sh_nuc    - potentially corresponding nucleotide sequence
       of sh_anc
     - shift_val - 1 or 2 indicates forward frame shift, whereas
       3*anchor_len-1 or 3*anchor_len-2 indicates
       backward shift
     - aa2re     - aa to codon re dict
     - anchor_len - length of the anchor
     - shift_id_pos - specify current shift name we are at

    """
    import re

    shift_id = [chr(i) for i in range(97, 107)]
    if 0 < shift_val < 3 * anchor_len - 2:
        # if shift_val in (1, 2):
        for j in range(len(sh_anc)):
            qcodon = "^"
            for k, aa in enumerate(sh_anc):
                if k == j:
                    qcodon += aa2re[aa] + "(?P<" + shift_id[shift_id_pos] + ">..*)"
                else:
                    qcodon += aa2re[aa]
            qcodon += "$"
            match = re.search(qcodon, sh_nuc)
            if match:
                qcodon = qcodon.replace("^", "").replace("$", "")
                shift_id_pos += 1
                return qcodon, shift_id_pos
        if not match:
            # failed to find a match (frameshift)
            return -1, shift_id_pos
    elif shift_val in (3 * anchor_len - 1, 3 * anchor_len - 2):
        shift_val = 3 * anchor_len - shift_val
        # obtain shifted anchor and corresponding nucl
        # first check if the shifted pos is just at the end of the
        # previous anchor.
        for j in range(1, len(sh_anc)):
            qcodon = "^"
            for k, aa in enumerate(sh_anc):
                if k == j - 1:
                    # will be considered in the next step
                    pass
                elif k == j:
                    qcodon += _merge_aa2re(
                        sh_anc[j - 1],
                        sh_anc[j],
                        shift_val,
                        aa2re,
                        shift_id[shift_id_pos].upper(),
                    )
                else:
                    qcodon += aa2re[aa]
            qcodon += "$"
            match = re.search(qcodon, sh_nuc)
            if match:
                qcodon = qcodon.replace("^", "").replace("$", "")
                shift_id_pos += 1
                return qcodon, shift_id_pos
        if not match:
            # failed to find a match (frameshift)
            return -1, shift_id_pos


def _merge_aa2re(aa1, aa2, shift_val, aa2re, reid):
    """Merge two amino acids based on detected frame shift value (PRIVATE)."""

    def get_aa_from_codonre(re_aa):
        aas = []
        m = 0
        for i in re_aa:
            if i == "[":
                m = -1
                aas.append("")
            elif i == "]":
                m = 0
                continue
            elif m == -1:
                aas[-1] = aas[-1] + i
            elif m == 0:
                aas.append(i)
        return aas

    scodon = list(map(get_aa_from_codonre, (aa2re[aa1], aa2re[aa2])))
    if shift_val == 1:
        intersect = "".join(set(scodon[0][2]) & set(scodon[1][0]))
        scodonre = "(?P<" + reid + ">"
        scodonre += (
            "["
            + scodon[0][0]
            + "]"
            + "["
            + scodon[0][1]
            + "]"
            + "["
            + intersect
            + "]"
            + "["
            + scodon[1][1]
            + "]"
            + "["
            + scodon[1][2]
            + "]"
        )
    elif shift_val == 2:
        intersect1 = "".join(set(scodon[0][1]) & set(scodon[1][0]))
        intersect2 = "".join(set(scodon[0][2]) & set(scodon[1][1]))
        scodonre = "(?P<" + reid + ">"
        scodonre += (
            "["
            + scodon[0][0]
            + "]"
            + "["
            + intersect1
            + "]"
            + "["
            + intersect2
            + "]"
            + "["
            + scodon[1][2]
            + "]"
        )
    scodonre += ")"
    return scodonre


def _get_codon_rec(
    pro, nucl, span_mode, gap_char, codon_table, complete_protein=False, max_score=10,
):
    """Generate codon alignment based on regular re match (PRIVATE).

    span_mode is a tuple returned by _check_corr. The first element
    is the span of a re search, and the second element is the mode
    for the match.

    mode
     - 0: direct match
     - 1: mismatch (no indels)
     - 2: frameshift

    """
    import re
    from Bio.Seq import Seq

    nucl_seq = nucl.seq.replace(gap_char, "")
    span = span_mode[0]
    mode = span_mode[1]
    aa2re = _get_aa_regex(codon_table)
    if mode in (0, 1):
        if len(pro.seq.replace(gap_char, "")) * 3 != (span[1] - span[0]):
            raise ValueError(
                f"Protein Record {pro.id} and "
                f"Nucleotide Record {nucl.id} do not match!"
            )
        aa_num = 0
        codon_seq = CodonSeq()
        for aa in pro.seq:
            if aa == "-":
                codon_seq += "---"
            elif complete_protein and aa_num == 0:
                this_codon = nucl_seq[span[0] : span[0] + 3]
                if not re.search(
                    _codons2re(codon_table.start_codons), str(this_codon.upper())
                ):
                    max_score -= 1
                    warnings.warn(
                        f"start codon of {pro.id} ({aa} {aa_num}) does not "
                        f"correspond to {nucl.id} ({this_codon})",
                        BiopythonWarning,
                    )
                if max_score == 0:
                    raise RuntimeError(
                        f"max_score reached for {nucl.id}! Please raise up "
                        "the tolerance to get an alignment in anyway"
                    )
                codon_seq += this_codon
                aa_num += 1
            else:
                this_codon = nucl_seq[span[0] + 3 * aa_num : span[0] + 3 * (aa_num + 1)]
                if this_codon.upper().translate(table=codon_table) != aa:
                    max_score -= 1
                    warnings.warn(
                        "%s(%s %d) does not correspond to %s(%s)"
                        % (pro.id, aa, aa_num, nucl.id, this_codon),
                        BiopythonWarning,
                    )
                if max_score == 0:
                    raise RuntimeError(
                        f"max_score reached for {nucl.id}! Please raise up "
                        "the tolerance to get an alignment in anyway"
                    )
                codon_seq += this_codon
                aa_num += 1
        return SeqRecord(codon_seq, id=nucl.id)
    elif mode == 2:
        from collections import deque

        shift_pos = deque([])
        shift_start = []
        match = span_mode[2]
        m_groupdict = list(match.groupdict().keys())
        # backward frameshift
        for i in m_groupdict:
            shift_pos.append(match.span(i))
            shift_start.append(match.start(i))
        rf_table = []
        i = match.start()
        while True:
            rf_table.append(i)
            i += 3
            if i in shift_start and m_groupdict[shift_start.index(i)].isupper():
                shift_index = shift_start.index(i)
                shift_val = 6 - (shift_pos[shift_index][1] - shift_pos[shift_index][0])
                rf_table.append(i)
                rf_table.append(i + 3 - shift_val)
                i = shift_pos[shift_index][1]
            elif i in shift_start and m_groupdict[shift_start.index(i)].islower():
                i = shift_pos[shift_start.index(i)][1]
            if i >= match.end():
                break
        codon_seq = CodonSeq()
        aa_num = 0
        for aa in pro.seq:
            if aa == "-":
                codon_seq += "---"
            elif complete_protein and aa_num == 0:
                this_codon = nucl_seq[rf_table[0] : rf_table[0] + 3]
                if not re.search(
                    _codons2re(codon_table.start_codons), str(this_codon.upper())
                ):
                    max_score -= 1
                    warnings.warn(
                        f"start codon of {pro.id}({aa} {aa_num}) does not "
                        f"correspond to {nucl.id}({this_codon})",
                        BiopythonWarning,
                    )
                    codon_seq += this_codon
                    aa_num += 1
            else:
                if (
                    aa_num < len(pro.seq.replace("-", "")) - 1
                    and rf_table[aa_num + 1] - rf_table[aa_num] - 3 < 0
                ):
                    max_score -= 1
                    start = rf_table[aa_num]
                    end = start + (3 - shift_val)
                    ngap = shift_val
                    this_codon = nucl_seq[start:end] + "-" * ngap
                elif rf_table[aa_num] - rf_table[aa_num - 1] - 3 > 0:
                    max_score -= 1
                    start = rf_table[aa_num - 1] + 3
                    end = rf_table[aa_num]
                    ngap = 3 - (rf_table[aa_num] - rf_table[aa_num - 1] - 3)
                    this_codon = (
                        nucl_seq[start:end]
                        + "-" * ngap
                        + nucl_seq[rf_table[aa_num] : rf_table[aa_num] + 3]
                    )
                else:
                    start = rf_table[aa_num]
                    end = start + 3
                    this_codon = nucl_seq[start:end]
                    if this_codon.upper().translate(table=codon_table) != aa:
                        max_score -= 1
                        warnings.warn(
                            f"Codon of {pro.id}({aa} {aa_num}) does not "
                            f"correspond to {nucl.id}({this_codon})",
                            BiopythonWarning,
                        )
                if max_score == 0:
                    raise RuntimeError(
                        f"max_score reached for {nucl.id}! Please raise up "
                        "the tolerance to get an alignment in anyway"
                    )
                codon_seq += this_codon
                aa_num += 1
        codon_seq.rf_table = rf_table
        return SeqRecord(codon_seq, id=nucl.id)


def _align_shift_recs(recs):
    """Build alignment according to the frameshift detected by _check_corr (PRIVATE).

    Argument:
     - recs - a list of SeqRecords containing a CodonSeq dictated
       by a rf_table (with frameshift in some of them).

    """

    def find_next_int(k, lst):
        idx = lst.index(k)
        p = 0
        while True:
            if isinstance(lst[idx + p], int):
                return lst[idx + p], p
            p += 1

    full_rf_table_lst = [rec.seq.get_full_rf_table() for rec in recs]
    rf_num = [0] * len(recs)
    for k, rec in enumerate(recs):
        for i in rec.seq.get_full_rf_table():
            if isinstance(i, int):
                rf_num[k] += 1
            # isinstance(i, float) should be True
            elif rec.seq[int(i) : int(i) + 3] == "---":
                rf_num[k] += 1
    if len(set(rf_num)) != 1:
        raise RuntimeError("Number of alignable codons unequal in given records")
    i = 0
    rec_num = len(recs)
    while True:
        add_lst = []
        try:
            col_rf_lst = [k[i] for k in full_rf_table_lst]
        except IndexError:
            # we probably reached the last codon
            break
        for j, k in enumerate(col_rf_lst):
            add_lst.append((j, int(k)))
            if isinstance(k, float) and recs[j].seq[int(k) : int(k) + 3] != "---":
                m, p = find_next_int(k, full_rf_table_lst[j])
                if (m - k) % 3 != 0:
                    gap_num = 3 - (m - k) % 3
                else:
                    gap_num = 0
                if gap_num != 0:
                    gaps = "-" * int(gap_num)
                    seq = CodonSeq(rf_table=recs[j].seq.rf_table)
                    seq += recs[j].seq[: int(k)] + gaps + recs[j].seq[int(k) :]
                    full_rf_table = full_rf_table_lst[j]
                    bp = full_rf_table.index(k)
                    full_rf_table = full_rf_table[:bp] + [
                        v + int(gap_num) for v in full_rf_table[bp + 1 :]
                    ]
                    full_rf_table_lst[j] = full_rf_table
                    recs[j].seq = seq
                add_lst.pop()
                gap_num += m - k
                i += p - 1
        if len(add_lst) != rec_num:
            for j, k in add_lst:
                seq = CodonSeq(rf_table=recs[j].seq.rf_table)
                gaps = "-" * int(gap_num)
                seq += recs[j].seq[: int(k)] + gaps + recs[j].seq[int(k) :]
                full_rf_table = full_rf_table_lst[j]
                bp = full_rf_table.index(k)
                inter_rf = []
                for t in range(0, len(gaps), 3):
                    inter_rf.append(k + t + 3.0)
                full_rf_table = (
                    full_rf_table[:bp]
                    + inter_rf
                    + [v + int(gap_num) for v in full_rf_table[bp:]]
                )
                full_rf_table_lst[j] = full_rf_table
                recs[j].seq = seq
        i += 1
    return recs


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
