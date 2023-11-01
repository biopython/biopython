# Copyright 2013 by Zheng Ruan (zruan1991@gmail.com).
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Code for performing calculations on codon alignments."""

import sys
from math import sqrt, erfc, log, floor
from heapq import heapify, heappop, heappush
from itertools import permutations
from collections import defaultdict, Counter

import numpy as np

from Bio.Align import Alignment
from Bio.Data import CodonTable


def calculate_dn_ds(alignment, method="NG86", codon_table=None, k=1, cfreq=None):
    """Calculate dN and dS of the given two sequences.

    Available methods:
        - NG86  - `Nei and Gojobori (1986)`_ (PMID 3444411).
        - LWL85 - `Li et al. (1985)`_ (PMID 3916709).
        - ML    - `Goldman and Yang (1994)`_ (PMID 7968486).
        - YN00  - `Yang and Nielsen (2000)`_ (PMID 10666704).

    .. _`Nei and Gojobori (1986)`: http://www.ncbi.nlm.nih.gov/pubmed/3444411
    .. _`Li et al. (1985)`: http://www.ncbi.nlm.nih.gov/pubmed/3916709
    .. _`Goldman and Yang (1994)`: http://mbe.oxfordjournals.org/content/11/5/725
    .. _`Yang and Nielsen (2000)`: https://doi.org/10.1093/oxfordjournals.molbev.a026236

    Arguments:
     - k  - transition/transversion rate ratio
     - cfreq - Current codon frequency vector can only be specified
       when you are using ML method. Possible ways of
       getting cfreq are: F1x4, F3x4 and F61.

    """
    if cfreq is None:
        cfreq = "F3x4"
    elif cfreq is not None and method != "ML":
        raise ValueError("cfreq can only be specified when you are using ML method")
    elif cfreq not in ("F1x4", "F3x4", "F61"):
        raise ValueError("cfreq must be 'F1x4', 'F3x4', or 'F61'")
    if codon_table is None:
        codon_table = CodonTable.generic_by_id[1]
    codons1 = []
    codons2 = []
    sequence1, sequence2 = alignment.sequences
    try:
        sequence1 = sequence1.seq  # stupid SeqRecord
    except AttributeError:
        pass
    sequence1 = str(sequence1)
    try:
        sequence2 = sequence2.seq  # stupid SeqRecord
    except AttributeError:
        pass
    sequence2 = str(sequence2)
    aligned1, aligned2 = alignment.aligned
    for block1, block2 in zip(aligned1, aligned2):
        start1, end1 = block1
        start2, end2 = block2
        codons1.extend(sequence1[i : i + 3] for i in range(start1, end1, 3))
        codons2.extend(sequence2[i : i + 3] for i in range(start2, end2, 3))
    bases = {"A", "T", "C", "G"}
    for codon1 in codons1:
        if not all(nucleotide in bases for nucleotide in codon1):
            raise ValueError(
                f"Unrecognized character in {codon1} in the target sequence"
                " (Codons consist of A, T, C or G)"
            )
    for codon2 in codons2:
        if not all(nucleotide in bases for nucleotide in codon2):
            raise ValueError(
                f"Unrecognized character in {codon2} in the query sequence"
                " (Codons consist of A, T, C or G)"
            )
    if method == "ML":
        return _ml(codons1, codons2, cfreq, codon_table)
    elif method == "NG86":
        return _ng86(codons1, codons2, k, codon_table)
    elif method == "LWL85":
        return _lwl85(codons1, codons2, codon_table)
    elif method == "YN00":
        return _yn00(codons1, codons2, codon_table)
    else:
        raise ValueError(f"Unknown method '{method}'")


#################################################################
#              private functions for NG86 method
#################################################################


def _ng86(codons1, codons2, k, codon_table):
    """NG86 method main function (PRIVATE)."""
    S_sites1, N_sites1 = _count_site_NG86(codons1, codon_table=codon_table, k=k)
    S_sites2, N_sites2 = _count_site_NG86(codons2, codon_table=codon_table, k=k)
    S_sites = (S_sites1 + S_sites2) / 2.0
    N_sites = (N_sites1 + N_sites2) / 2.0
    SN = [0, 0]
    for codon1, codon2 in zip(codons1, codons2):
        SN = [
            m + n
            for m, n in zip(
                SN, _count_diff_NG86(codon1, codon2, codon_table=codon_table)
            )
        ]

    ps = SN[0] / S_sites
    pn = SN[1] / N_sites
    if ps < 3 / 4:
        dS = abs(-3.0 / 4 * log(1 - 4.0 / 3 * ps))
    else:
        dS = -1
    if pn < 3 / 4:
        dN = abs(-3.0 / 4 * log(1 - 4.0 / 3 * pn))
    else:
        dN = -1
    return dN, dS


def _count_site_NG86(codons, codon_table, k=1):
    """Count synonymous and non-synonymous sites of a list of codons (PRIVATE).

    Arguments:
     - codons - A list of three letter codons.
     - k - transition/transversion rate ratio.

    """
    S_site = 0  # synonymous sites
    N_site = 0  # non-synonymous sites
    purine = ("A", "G")
    pyrimidine = ("T", "C")
    bases = ("A", "T", "C", "G")
    for codon in codons:
        neighbor_codon = {"transition": [], "transversion": []}
        # classify neighbor codons
        codon = codon.replace("U", "T")
        for i, nucleotide in enumerate(codon):
            for base in bases:
                if nucleotide == base:
                    pass
                elif nucleotide in purine and base in purine:
                    codon_chars = list(codon)
                    codon_chars[i] = base
                    this_codon = "".join(codon_chars)
                    neighbor_codon["transition"].append(this_codon)
                elif nucleotide in pyrimidine and base in pyrimidine:
                    codon_chars = list(codon)
                    codon_chars[i] = base
                    this_codon = "".join(codon_chars)
                    neighbor_codon["transition"].append(this_codon)
                else:
                    codon_chars = list(codon)
                    codon_chars[i] = base
                    this_codon = "".join(codon_chars)
                    neighbor_codon["transversion"].append(this_codon)
        # count synonymous and non-synonymous sites
        aa = codon_table.forward_table[codon]
        this_codon_N_site = this_codon_S_site = 0
        for neighbor in neighbor_codon["transition"]:
            if neighbor in codon_table.stop_codons:
                this_codon_N_site += 1
            elif codon_table.forward_table[neighbor] == aa:
                this_codon_S_site += 1
            else:
                this_codon_N_site += 1
        for neighbor in neighbor_codon["transversion"]:
            if neighbor in codon_table.stop_codons:
                this_codon_N_site += k
            elif codon_table.forward_table[neighbor] == aa:
                this_codon_S_site += k
            else:
                this_codon_N_site += k
        norm_const = (this_codon_N_site + this_codon_S_site) / 3
        S_site += this_codon_S_site / norm_const
        N_site += this_codon_N_site / norm_const
    return (S_site, N_site)


def _count_diff_NG86(codon1, codon2, codon_table):
    """Count differences between two codons, three-letter string (PRIVATE).

    The function will take multiple pathways from codon1 to codon2
    into account.
    """
    SN = [0, 0]  # synonymous and nonsynonymous counts
    if codon1 == codon2:
        return SN
    else:
        diff_pos = [
            i
            for i, (nucleotide1, nucleotide2) in enumerate(zip(codon1, codon2))
            if nucleotide1 != nucleotide2
        ]

        def compare_codon(codon1, codon2, codon_table, weight=1):
            """Compare two codon accounting for different pathways."""
            sd = nd = 0
            if len(set(map(codon_table.forward_table.get, [codon1, codon2]))) == 1:
                sd += weight
            else:
                nd += weight
            return (sd, nd)

        if len(diff_pos) == 1:
            SN = [
                i + j
                for i, j in zip(
                    SN, compare_codon(codon1, codon2, codon_table=codon_table)
                )
            ]
        elif len(diff_pos) == 2:
            for i in diff_pos:
                temp_codon = codon1[:i] + codon2[i] + codon1[i + 1 :]
                SN = [
                    i + j
                    for i, j in zip(
                        SN,
                        compare_codon(
                            codon1, temp_codon, codon_table=codon_table, weight=0.5
                        ),
                    )
                ]
                SN = [
                    i + j
                    for i, j in zip(
                        SN,
                        compare_codon(
                            temp_codon, codon2, codon_table=codon_table, weight=0.5
                        ),
                    )
                ]
        elif len(diff_pos) == 3:
            paths = list(permutations([0, 1, 2], 3))
            tmp_codon = []
            for index1, index2, index3 in paths:
                tmp1 = codon1[:index1] + codon2[index1] + codon1[index1 + 1 :]
                tmp2 = tmp1[:index2] + codon2[index2] + tmp1[index2 + 1 :]
                tmp_codon.append((tmp1, tmp2))
                SN = [
                    i + j
                    for i, j in zip(
                        SN, compare_codon(codon1, tmp1, codon_table, weight=0.5 / 3)
                    )
                ]
                SN = [
                    i + j
                    for i, j in zip(
                        SN, compare_codon(tmp1, tmp2, codon_table, weight=0.5 / 3)
                    )
                ]
                SN = [
                    i + j
                    for i, j in zip(
                        SN, compare_codon(tmp2, codon2, codon_table, weight=0.5 / 3)
                    )
                ]
    return SN


#################################################################
#               private functions for LWL85 method
#################################################################


def _lwl85(codons1, codons2, codon_table):
    """LWL85 method main function (PRIVATE).

    Nomenclature is according to Li et al. (1985), PMID 3916709.
    """
    codon_fold_dict = _get_codon_fold(codon_table)
    # count number of sites in different degenerate classes
    fold0 = [0, 0]
    fold2 = [0, 0]
    fold4 = [0, 0]
    for codon in codons1 + codons2:
        fold_num = codon_fold_dict[codon]
        for f in fold_num:
            if f == "0":
                fold0[0] += 1
            elif f == "2":
                fold2[0] += 1
            elif f == "4":
                fold4[0] += 1
    L = [sum(fold0) / 2.0, sum(fold2) / 2.0, sum(fold4) / 2.0]
    # count number of differences in different degenerate classes
    PQ = [0] * 6  # with P0, P2, P4, Q0, Q2, Q4 in each position
    for codon1, codon2 in zip(codons1, codons2):
        if codon1 == codon2:
            continue
        PQ = [
            i + j
            for i, j in zip(PQ, _diff_codon(codon1, codon2, fold_dict=codon_fold_dict))
        ]
    PQ = [i / j for i, j in zip(PQ, L * 2)]
    P = PQ[:3]
    Q = PQ[3:]
    A = [
        (1.0 / 2) * log(1.0 / (1 - 2 * i - j)) - (1.0 / 4) * log(1.0 / (1 - 2 * j))
        for i, j in zip(P, Q)
    ]
    B = [(1.0 / 2) * log(1.0 / (1 - 2 * i)) for i in Q]
    dS = 3 * (L[2] * A[1] + L[2] * (A[2] + B[2])) / (L[1] + 3 * L[2])
    dN = 3 * (L[2] * B[1] + L[0] * (A[0] + B[0])) / (2 * L[1] + 3 * L[0])
    return dN, dS


def _get_codon_fold(codon_table):
    """Classify different position in a codon into different folds (PRIVATE)."""
    fold_table = {}
    forward_table = codon_table.forward_table
    bases = {"A", "T", "C", "G"}
    for codon in forward_table:
        if "U" in codon:
            continue
        fold = ""
        codon_base_lst = list(codon)
        for i, base in enumerate(codon_base_lst):
            other_bases = bases - set(base)
            aa = []
            for other_base in other_bases:
                codon_base_lst[i] = other_base
                try:
                    aa.append(forward_table["".join(codon_base_lst)])
                except KeyError:
                    aa.append("stop")
            if aa.count(forward_table[codon]) == 0:
                fold += "0"
            elif aa.count(forward_table[codon]) in (1, 2):
                fold += "2"
            elif aa.count(forward_table[codon]) == 3:
                fold += "4"
            else:
                raise RuntimeError(
                    "Unknown Error, cannot assign the position to a fold"
                )
            codon_base_lst[i] = base
        fold_table[codon] = fold
    return fold_table


def _diff_codon(codon1, codon2, fold_dict):
    """Count number of different substitution types between two codons (PRIVATE).

    returns tuple (P0, P2, P4, Q0, Q2, Q4)

    Nomenclature is according to Li et al. (1958), PMID 3916709.
    """
    P0 = P2 = P4 = Q0 = Q2 = Q4 = 0
    fold_num = fold_dict[codon1]
    purine = ("A", "G")
    pyrimidine = ("T", "C")
    for n, (nucleotide1, nucleotide2) in enumerate(zip(codon1, codon2)):
        if nucleotide1 == nucleotide2:
            pass
        elif nucleotide1 in purine and nucleotide2 in purine:
            if fold_num[n] == "0":
                P0 += 1
            elif fold_num[n] == "2":
                P2 += 1
            elif fold_num[n] == "4":
                P4 += 1
            else:
                raise RuntimeError("Unexpected fold_num %d" % fold_num[n])
        elif nucleotide1 in pyrimidine and nucleotide2 in pyrimidine:
            if fold_num[n] == "0":
                P0 += 1
            elif fold_num[n] == "2":
                P2 += 1
            elif fold_num[n] == "4":
                P4 += 1
            else:
                raise RuntimeError("Unexpected fold_num %d" % fold_num[n])
        else:
            # nucleotide1 in purine and nucleotide2 in pyrimidine, or
            # nucleotide1 in pyrimidine and nucleotide2 in purine
            if fold_num[n] == "0":
                Q0 += 1
            elif fold_num[n] == "2":
                Q2 += 1
            elif fold_num[n] == "4":
                Q4 += 1
            else:
                raise RuntimeError("Unexpected fold_num %d" % fold_num[n])
    return (P0, P2, P4, Q0, Q2, Q4)


#################################################################
#               private functions for YN00 method
#################################################################


def _yn00(codons1, codons2, codon_table):
    """YN00 method main function (PRIVATE).

    Nomenclature is according to Yang and Nielsen (2000), PMID 10666704.
    """
    from scipy.linalg import expm

    fcodon = [
        {"A": 0, "G": 0, "C": 0, "T": 0},
        {"A": 0, "G": 0, "C": 0, "T": 0},
        {"A": 0, "G": 0, "C": 0, "T": 0},
    ]
    codon_fold_dict = _get_codon_fold(codon_table)
    fold0_cnt = defaultdict(int)
    fold4_cnt = defaultdict(int)
    for codon in codons1 + codons2:
        # count sites at different codon position
        fcodon[0][codon[0]] += 1
        fcodon[1][codon[1]] += 1
        fcodon[2][codon[2]] += 1
        # count sites in different degenerate fold class
        fold_num = codon_fold_dict[codon]
        for i, f in enumerate(fold_num):
            if f == "0":
                fold0_cnt[codon[i]] += 1
            elif f == "4":
                fold4_cnt[codon[i]] += 1
    f0_total = sum(fold0_cnt.values())
    f4_total = sum(fold4_cnt.values())
    for i, j in zip(fold0_cnt, fold4_cnt):
        fold0_cnt[i] = fold0_cnt[i] / f0_total
        fold4_cnt[i] = fold4_cnt[i] / f4_total
    # TODO:
    # the initial kappa is different from what yn00 gives,
    # try to find the problem.
    TV = _get_TV(codons1, codons2, codon_table=codon_table)
    k04 = (_get_kappa_t(fold0_cnt, TV), _get_kappa_t(fold4_cnt, TV))
    kappa = (f0_total * k04[0] + f4_total * k04[1]) / (f0_total + f4_total)
    # kappa = 2.4285
    # count synonymous sites and non-synonymous sites
    for i in range(3):
        tot = sum(fcodon[i].values())
        fcodon[i] = {j: k / tot for j, k in fcodon[i].items()}
    pi = defaultdict(int)
    for codon in list(codon_table.forward_table.keys()) + codon_table.stop_codons:
        if "U" not in codon:
            pi[codon] = 0
    for codon in codons1 + codons2:
        pi[codon] += 1
    S_sites1, N_sites1, bfreqSN1 = _count_site_YN00(
        codons1, codons2, pi, k=kappa, codon_table=codon_table
    )
    S_sites2, N_sites2, bfreqSN2 = _count_site_YN00(
        codons2, codons1, pi, k=kappa, codon_table=codon_table
    )
    N_sites = (N_sites1 + N_sites2) / 2
    S_sites = (S_sites1 + S_sites2) / 2
    bfreqSN = [{"A": 0, "T": 0, "C": 0, "G": 0}, {"A": 0, "T": 0, "C": 0, "G": 0}]
    for i in range(2):
        for base in ("A", "T", "C", "G"):
            bfreqSN[i][base] = (bfreqSN1[i][base] + bfreqSN2[i][base]) / 2
    # use NG86 method to get initial t and w
    SN = [0, 0]
    for codon1, codon2 in zip(codons1, codons2):
        SN = [
            m + n
            for m, n in zip(
                SN, _count_diff_NG86(codon1, codon2, codon_table=codon_table)
            )
        ]
    ps = SN[0] / S_sites
    pn = SN[1] / N_sites
    p = sum(SN) / (S_sites + N_sites)
    w = log(1 - 4.0 / 3 * pn) / log(1 - 4.0 / 3 * ps)
    t = -3 / 4 * log(1 - 4 / 3 * p)
    tolerance = 1e-5
    dSdN_pre = [0, 0]
    for temp in range(20):
        # count synonymous and nonsynonymous differences under kappa, w, t
        codons = [
            codon
            for codon in list(codon_table.forward_table.keys())
            + codon_table.stop_codons
            if "U" not in codon
        ]
        Q = _get_Q(pi, kappa, w, codons, codon_table)
        P = expm(Q * t)
        TV = [0, 0, 0, 0]  # synonymous/nonsynonymous transition/transversion
        codon_npath = Counter(zip(codons1, codons2))
        for (nucleotide1, nucleotide2), count in codon_npath.items():
            tv = _count_diff_YN00(nucleotide1, nucleotide2, P, codons, codon_table)
            TV = [m + n * count for m, n in zip(TV, tv)]
        TV = (TV[0] / S_sites, TV[1] / S_sites), (TV[2] / N_sites, TV[3] / N_sites)
        # according to the DistanceF84() function of yn00.c in paml,
        # the t (e.q. 10) appears in PMID: 10666704 is dS and dN
        dSdN = []
        for f, tv in zip(bfreqSN, TV):
            dSdN.append(_get_kappa_t(f, tv, t=True))
        t = dSdN[0] * 3 * S_sites / (S_sites + N_sites) + dSdN[1] * 3 * N_sites / (
            S_sites + N_sites
        )
        w = dSdN[1] / dSdN[0]
        if all(abs(i - j) < tolerance for i, j in zip(dSdN, dSdN_pre)):
            return dSdN[1], dSdN[0]  # dN, dS
        dSdN_pre = dSdN


def _get_TV(codons1, codons2, codon_table):
    """Get TV (PRIVATE).

    Arguments:
     - T - proportions of transitional differences
     - V - proportions of transversional differences

    """
    purine = ("A", "G")
    pyrimidine = ("C", "T")
    TV = [0, 0]
    sites = 0
    for codon1, codon2 in zip(codons1, codons2):
        for nucleotide1, nucleotide2 in zip(codon1, codon2):
            if nucleotide1 == nucleotide2:
                pass
            elif nucleotide1 in purine and nucleotide2 in purine:
                TV[0] += 1
            elif nucleotide1 in pyrimidine and nucleotide2 in pyrimidine:
                TV[0] += 1
            else:
                TV[1] += 1
            sites += 1
    return (TV[0] / sites, TV[1] / sites)
    # return (TV[0], TV[1])


def _get_kappa_t(pi, TV, t=False):
    """Calculate kappa (PRIVATE).

    The following formula and variable names are according to PMID: 10666704
    """
    pi["Y"] = pi["T"] + pi["C"]
    pi["R"] = pi["A"] + pi["G"]
    A = (
        2 * (pi["T"] * pi["C"] + pi["A"] * pi["G"])
        + 2
        * (
            pi["T"] * pi["C"] * pi["R"] / pi["Y"]
            + pi["A"] * pi["G"] * pi["Y"] / pi["R"]
        )
        * (1 - TV[1] / (2 * pi["Y"] * pi["R"]))
        - TV[0]
    ) / (2 * (pi["T"] * pi["C"] / pi["Y"] + pi["A"] * pi["G"] / pi["R"]))
    B = 1 - TV[1] / (2 * pi["Y"] * pi["R"])
    a = -0.5 * log(A)  # this seems to be an error in YANG's original paper
    b = -0.5 * log(B)
    kappaF84 = a / b - 1
    if t is False:
        kappaHKY85 = 1 + (
            pi["T"] * pi["C"] / pi["Y"] + pi["A"] * pi["G"] / pi["R"]
        ) * kappaF84 / (pi["T"] * pi["C"] + pi["A"] * pi["G"])
        return kappaHKY85
    else:
        t = (
            4 * pi["T"] * pi["C"] * (1 + kappaF84 / pi["Y"])
            + 4 * pi["A"] * pi["G"] * (1 + kappaF84 / pi["R"])
            + 4 * pi["Y"] * pi["R"]
        ) * b
        return t


def _count_site_YN00(codons1, codons2, pi, k, codon_table):
    """Site counting method from Ina / Yang and Nielsen (PRIVATE).

    Method from `Ina (1995)`_ as modified by `Yang and Nielsen (2000)`_.
    This will return the total number of synonymous and nonsynonymous sites
    and base frequencies in each category. The function is equivalent to
    the ``CountSites()`` function in ``yn00.c`` of PAML.

    .. _`Ina (1995)`: https://doi.org/10.1007/BF00167113
    .. _`Yang and Nielsen (2000)`: https://doi.org/10.1093/oxfordjournals.molbev.a026236

    """
    length = len(codons1)
    assert length == len(codons2)
    purine = ("A", "G")
    pyrimidine = ("T", "C")
    bases = ("A", "T", "C", "G")
    codon_dict = codon_table.forward_table
    stop = codon_table.stop_codons
    codon_npath = Counter(zip(codons1, codons2))
    S_sites = N_sites = 0
    freqSN = [
        {"A": 0, "T": 0, "C": 0, "G": 0},  # synonymous
        {"A": 0, "T": 0, "C": 0, "G": 0},
    ]  # nonsynonymous
    for codon_pair, npath in codon_npath.items():
        codon = codon_pair[0]
        S = N = 0
        for pos in range(3):
            for base in bases:
                if codon[pos] == base:
                    continue
                neighbor_codon = codon[:pos] + base + codon[pos + 1 :]
                if neighbor_codon in stop:
                    continue
                weight = pi[neighbor_codon]
                if codon[pos] in pyrimidine and base in pyrimidine:
                    weight *= k
                elif codon[pos] in purine and base in purine:
                    weight *= k
                if codon_dict[codon] == codon_dict[neighbor_codon]:
                    S += weight
                    freqSN[0][base] += weight * npath
                else:
                    N += weight
                    freqSN[1][base] += weight * npath
        S_sites += S * npath
        N_sites += N * npath
    norm_const = 3 * length / (S_sites + N_sites)
    S_sites *= norm_const
    N_sites *= norm_const
    for i in freqSN:
        norm_const = sum(i.values())
        for b in i:
            i[b] /= norm_const
    return S_sites, N_sites, freqSN


def _count_diff_YN00(codon1, codon2, P, codons, codon_table):
    """Count differences between two codons (three-letter string; PRIVATE).

    The function will weighted multiple pathways from codon1 to codon2
    according to P matrix of codon substitution. The proportion
    of transition and transversion (TV) will also be calculated in
    the function.
    """
    TV = [
        0,
        0,
        0,
        0,
    ]  # transition and transversion counts (synonymous and nonsynonymous)
    if codon1 == codon2:
        return TV
    else:
        diff_pos = [
            i
            for i, (nucleotide1, nucleotide2) in enumerate(zip(codon1, codon2))
            if nucleotide1 != nucleotide2
        ]

        def count_TV(codon1, codon2, diff, codon_table, weight=1):
            purine = ("A", "G")
            pyrimidine = ("T", "C")
            dic = codon_table.forward_table
            stop = codon_table.stop_codons
            if codon1 in stop or codon2 in stop:
                # stop codon is always considered as nonsynonymous
                if codon1[diff] in purine and codon2[diff] in purine:
                    return [0, 0, weight, 0]
                elif codon1[diff] in pyrimidine and codon2[diff] in pyrimidine:
                    return [0, 0, weight, 0]
                else:
                    return [0, 0, 0, weight]
            elif dic[codon1] == dic[codon2]:
                if codon1[diff] in purine and codon2[diff] in purine:
                    return [weight, 0, 0, 0]
                elif codon1[diff] in pyrimidine and codon2[diff] in pyrimidine:
                    return [weight, 0, 0, 0]
                else:
                    return [0, weight, 0, 0]
            else:
                if codon1[diff] in purine and codon2[diff] in purine:
                    return [0, 0, weight, 0]
                elif codon1[diff] in pyrimidine and codon2[diff] in pyrimidine:
                    return [0, 0, weight, 0]
                else:
                    return [0, 0, 0, weight]

        if len(diff_pos) == 1:
            TV = [
                p + q
                for p, q in zip(TV, count_TV(codon1, codon2, diff_pos[0], codon_table))
            ]
        elif len(diff_pos) == 2:
            tmp_codons = [codon1[:i] + codon2[i] + codon1[i + 1 :] for i in diff_pos]
            path_prob = []
            for codon in tmp_codons:
                codon_idx = list(map(codons.index, [codon1, codon, codon2]))
                prob = (P[codon_idx[0], codon_idx[1]], P[codon_idx[1], codon_idx[2]])
                path_prob.append(prob[0] * prob[1])
            path_prob = [2 * i / sum(path_prob) for i in path_prob]
            for n, i in enumerate(diff_pos):
                codon = codon1[:i] + codon2[i] + codon1[i + 1 :]
                TV = [
                    p + q
                    for p, q in zip(
                        TV,
                        count_TV(
                            codon1, codon, i, codon_table, weight=path_prob[n] / 2
                        ),
                    )
                ]
                TV = [
                    p + q
                    for p, q in zip(
                        TV,
                        count_TV(
                            codon1, codon, i, codon_table, weight=path_prob[n] / 2
                        ),
                    )
                ]
        elif len(diff_pos) == 3:
            paths = list(permutations([0, 1, 2], 3))
            path_prob = []
            tmp_codons = []
            for index1, index2, index3 in paths:
                tmp1 = codon1[:index1] + codon2[index1] + codon1[index1 + 1 :]
                tmp2 = tmp1[:index2] + codon2[index2] + tmp1[index2 + 1 :]
                tmp_codons.append((tmp1, tmp2))
                codon_idx = list(map(codons.index, [codon1, tmp1, tmp2, codon2]))
                prob = (
                    P[codon_idx[0], codon_idx[1]],
                    P[codon_idx[1], codon_idx[2]],
                    P[codon_idx[2], codon_idx[3]],
                )
                path_prob.append(prob[0] * prob[1] * prob[2])
            path_prob = [3 * i / sum(path_prob) for i in path_prob]
            for codon, j, k in zip(tmp_codons, path_prob, paths):
                TV = [
                    p + q
                    for p, q in zip(
                        TV, count_TV(codon1, codon[0], k[0], codon_table, weight=j / 3)
                    )
                ]
                TV = [
                    p + q
                    for p, q in zip(
                        TV,
                        count_TV(codon[0], codon[1], k[1], codon_table, weight=j / 3),
                    )
                ]
                TV = [
                    p + q
                    for p, q in zip(
                        TV, count_TV(codon[1], codon2, k[1], codon_table, weight=j / 3)
                    )
                ]
    return TV


#################################################################
#        private functions for Maximum Likelihood method
#################################################################


def _ml(codons1, codons2, cmethod, codon_table):
    """ML method main function (PRIVATE)."""
    from scipy.optimize import minimize

    pi = _get_pi(codons1, codons2, cmethod, codon_table=codon_table)
    codon_cnt = Counter(zip(codons1, codons2))
    codons = [
        codon
        for codon in list(codon_table.forward_table.keys()) + codon_table.stop_codons
        if "U" not in codon
    ]

    # apply optimization
    def func(
        params, pi=pi, codon_cnt=codon_cnt, codons=codons, codon_table=codon_table
    ):
        """Temporary function, params = [t, k, w]."""
        return -_likelihood_func(
            params[0],
            params[1],
            params[2],
            pi,
            codon_cnt,
            codons=codons,
            codon_table=codon_table,
        )

    # count sites
    opt_res = minimize(
        func,
        [1, 0.1, 2],
        method="L-BFGS-B",
        bounds=((1e-10, 20), (1e-10, 20), (1e-10, 10)),
        tol=1e-5,
    )
    t, k, w = opt_res.x
    Q = _get_Q(pi, k, w, codons, codon_table)
    Sd = Nd = 0
    for i, codon1 in enumerate(codons):
        for j, codon2 in enumerate(codons):
            if i != j:
                try:
                    if (
                        codon_table.forward_table[codon1]
                        == codon_table.forward_table[codon2]
                    ):
                        # synonymous count
                        Sd += pi[codon1] * Q[i, j]
                    else:
                        # nonsynonymous count
                        Nd += pi[codon1] * Q[i, j]
                except KeyError:
                    # This is probably due to stop codons
                    pass
    Sd *= t
    Nd *= t

    # count differences (with w fixed to 1)
    def func_w1(
        params, pi=pi, codon_cnt=codon_cnt, codons=codons, codon_table=codon_table
    ):
        """Temporary function, params = [t, k]. w is fixed to 1."""
        return -_likelihood_func(
            params[0],
            params[1],
            1.0,
            pi,
            codon_cnt,
            codons=codons,
            codon_table=codon_table,
        )

    opt_res = minimize(
        func_w1,
        [1, 0.1],
        method="L-BFGS-B",
        bounds=((1e-10, 20), (1e-10, 20)),
        tol=1e-5,
    )
    t, k = opt_res.x
    w = 1.0
    Q = _get_Q(pi, k, w, codons, codon_table)
    rhoS = rhoN = 0
    for i, codon1 in enumerate(codons):
        for j, codon2 in enumerate(codons):
            if i != j:
                try:
                    if (
                        codon_table.forward_table[codon1]
                        == codon_table.forward_table[codon2]
                    ):
                        # synonymous count
                        rhoS += pi[codon1] * Q[i, j]
                    else:
                        # nonsynonymous count
                        rhoN += pi[codon1] * Q[i, j]
                except KeyError:
                    # This is probably due to stop codons
                    pass
    rhoS *= 3
    rhoN *= 3
    dN = Nd / rhoN
    dS = Sd / rhoS
    return dN, dS


def _get_pi(codons1, codons2, cmethod, codon_table):
    """Obtain codon frequency dict (pi) from two codon list (PRIVATE).

    This function is designed for ML method. Available counting methods
    (cfreq) are F1x4, F3x4 and F64.
    """
    # TODO:
    # Stop codon should not be allowed according to Yang.
    # Try to modify this!
    pi = {}
    if cmethod == "F1x4":
        fcodon = Counter(
            nucleotide for codon in codons1 + codons2 for nucleotide in codon
        )
        tot = sum(fcodon.values())
        fcodon = {j: k / tot for j, k in fcodon.items()}
        for codon in codon_table.forward_table.keys() + codon_table.stop_codons:
            if "U" not in codon:
                pi[codon] = fcodon[codon[0]] * fcodon[codon[1]] * fcodon[codon[2]]
    elif cmethod == "F3x4":
        # three codon position
        fcodon = [
            {"A": 0, "G": 0, "C": 0, "T": 0},
            {"A": 0, "G": 0, "C": 0, "T": 0},
            {"A": 0, "G": 0, "C": 0, "T": 0},
        ]
        for codon in codons1 + codons2:
            fcodon[0][codon[0]] += 1
            fcodon[1][codon[1]] += 1
            fcodon[2][codon[2]] += 1
        for i in range(3):
            tot = sum(fcodon[i].values())
            fcodon[i] = {j: k / tot for j, k in fcodon[i].items()}
        for codon in list(codon_table.forward_table.keys()) + codon_table.stop_codons:
            if "U" not in codon:
                pi[codon] = (
                    fcodon[0][codon[0]] * fcodon[1][codon[1]] * fcodon[2][codon[2]]
                )
    elif cmethod == "F61":
        for codon in codon_table.forward_table.keys() + codon_table.stop_codons:
            if "U" not in codon:
                pi[codon] = 0.1
        for codon in codons1 + codons2:
            pi[codon] += 1
        tot = sum(pi.values())
        pi = {j: k / tot for j, k in pi.items()}
    return pi


def _q(codon1, codon2, pi, k, w, codon_table):
    """Q matrix for codon substitution (PRIVATE).

    Arguments:
     - codon1, codon2  : three letter codon string
     - pi              : expected codon frequency
     - k               : transition/transversion ratio
     - w               : nonsynonymous/synonymous rate ratio
     - codon_table     : Bio.Data.CodonTable object

    """
    if codon1 == codon2:
        # diagonal elements is the sum of all other elements
        return 0
    if codon1 in codon_table.stop_codons or codon2 in codon_table.stop_codons:
        return 0
    if (codon1 not in pi) or (codon2 not in pi):
        return 0
    purine = ("A", "G")
    pyrimidine = ("T", "C")
    diff = [
        (i, nucleotide1, nucleotide2)
        for i, (nucleotide1, nucleotide2) in enumerate(zip(codon1, codon2))
        if nucleotide1 != nucleotide2
    ]
    if len(diff) >= 2:
        return 0
    if codon_table.forward_table[codon1] == codon_table.forward_table[codon2]:
        # synonymous substitution
        if diff[0][1] in purine and diff[0][2] in purine:
            # transition
            return k * pi[codon2]
        elif diff[0][1] in pyrimidine and diff[0][2] in pyrimidine:
            # transition
            return k * pi[codon2]
        else:
            # transversion
            return pi[codon2]
    else:
        # nonsynonymous substitution
        if diff[0][1] in purine and diff[0][2] in purine:
            # transition
            return w * k * pi[codon2]
        elif diff[0][1] in pyrimidine and diff[0][2] in pyrimidine:
            # transition
            return w * k * pi[codon2]
        else:
            # transversion
            return w * pi[codon2]


def _get_Q(pi, k, w, codons, codon_table):
    """Q matrix for codon substitution (PRIVATE)."""
    codon_num = len(codons)
    Q = np.zeros((codon_num, codon_num))
    for i1, codon1 in enumerate(codons):
        for i2, codon2 in enumerate(codons):
            if i1 != i2:
                Q[i1, i2] = _q(codon1, codon2, pi, k, w, codon_table=codon_table)
    nucl_substitutions = 0
    for i, codon in enumerate(codons):
        Q[i, i] = -sum(Q[i, :])
        try:
            nucl_substitutions += pi[codon] * (-Q[i, i])
        except KeyError:
            pass
    Q /= nucl_substitutions
    return Q


def _likelihood_func(t, k, w, pi, codon_cnt, codons, codon_table):
    """Likelihood function for ML method (PRIVATE)."""
    from scipy.linalg import expm

    Q = _get_Q(pi, k, w, codons, codon_table)
    P = expm(Q * t)
    likelihood = 0
    for i, codon1 in enumerate(codons):
        for j, codon2 in enumerate(codons):
            if (codon1, codon2) in codon_cnt:
                if P[i, j] * pi[codon1] <= 0:
                    likelihood += codon_cnt[(codon1, codon2)] * 0
                else:
                    likelihood += codon_cnt[(codon1, codon2)] * log(
                        pi[codon1] * P[i, j]
                    )
    return likelihood


def calculate_dn_ds_matrix(alignment, method="NG86", codon_table=None):
    """Calculate dN and dS pairwise for the multiple alignment, and return as matrices.

    Argument:
     - method       - Available methods include NG86, LWL85, YN00 and ML.
     - codon_table  - Codon table to use for forward translation.

    """
    from Bio.Phylo.TreeConstruction import DistanceMatrix

    if codon_table is None:
        codon_table = CodonTable.generic_by_id[1]
    sequences = alignment.sequences
    coordinates = alignment.coordinates
    names = [record.id for record in sequences]
    size = len(names)
    dn_matrix = []
    ds_matrix = []
    for i in range(size):
        dn_matrix.append([])
        ds_matrix.append([])
        for j in range(i):
            pairwise_sequences = [sequences[i], sequences[j]]
            pairwise_coordinates = coordinates[(i, j), :]
            pairwise_alignment = Alignment(pairwise_sequences, pairwise_coordinates)
            dn, ds = calculate_dn_ds(
                pairwise_alignment, method=method, codon_table=codon_table
            )
            dn_matrix[i].append(dn)
            ds_matrix[i].append(ds)
        dn_matrix[i].append(0.0)
        ds_matrix[i].append(0.0)
    dn_dm = DistanceMatrix(names, matrix=dn_matrix)
    ds_dm = DistanceMatrix(names, matrix=ds_matrix)
    return dn_dm, ds_dm


def mktest(alignment, species=None, codon_table=None):
    """McDonald-Kreitman test for neutrality.

    Implement the McDonald-Kreitman test for neutrality (PMID: 1904993)
    This method counts changes rather than sites
    (http://mkt.uab.es/mkt/help_mkt.asp).

    Arguments:
     - alignment    - Alignment of gene nucleotide sequences to compare.
     - species      - List of the species ID for each sequence in the alignment.
       Typically, the species ID is the species name as a string, or an integer.
     - codon_table  - Codon table to use for forward translation.

    Return the p-value of test result.
    """
    if codon_table is None:
        codon_table = CodonTable.generic_by_id[1]
    G, nonsyn_G = _get_codon2codon_matrix(codon_table=codon_table)
    unique_species = set(species)
    sequences = []
    for sequence in alignment.sequences:
        try:
            sequence = sequence.seq
        except AttributeError:
            pass
        sequence = str(sequence)
        sequences.append(sequence)
    syn_fix, nonsyn_fix, syn_poly, nonsyn_poly = 0, 0, 0, 0
    starts = sys.maxsize
    for ends in alignment.coordinates.transpose():
        step = min(ends - starts)
        for j in range(0, step, 3):
            codons = {key: [] for key in unique_species}
            for key, sequence, start in zip(species, sequences, starts):
                codon = sequence[start + j : start + j + 3]
                codons[key].append(codon)
            fixed = True
            all_codons = set()
            for value in codons.values():
                value = set(value)
                if len(value) > 1:
                    fixed = False
                all_codons.update(value)
            if len(all_codons) == 1:
                continue
            nonsyn = _count_replacement(all_codons, nonsyn_G)
            syn = _count_replacement(all_codons, G) - nonsyn
            if fixed is True:
                # fixed
                nonsyn_fix += nonsyn
                syn_fix += syn
            else:
                # not fixed
                nonsyn_poly += nonsyn
                syn_poly += syn
        starts = ends
    return _G_test([syn_fix, nonsyn_fix, syn_poly, nonsyn_poly])


def _get_codon2codon_matrix(codon_table):
    """Get codon codon substitution matrix (PRIVATE).

    Elements in the matrix are number of synonymous and nonsynonymous
    substitutions required for the substitution.
    """
    bases = ("A", "T", "C", "G")
    codons = [
        codon
        for codon in list(codon_table.forward_table.keys()) + codon_table.stop_codons
        if "U" not in codon
    ]
    # set up codon_dict considering stop codons
    codon_dict = codon_table.forward_table.copy()
    for stop in codon_table.stop_codons:
        codon_dict[stop] = "stop"
    # count site
    num = len(codons)
    G = {}  # graph for substitution
    nonsyn_G = {}  # graph for nonsynonymous substitution
    graph = {}
    graph_nonsyn = {}
    for i, codon in enumerate(codons):
        graph[codon] = {}
        graph_nonsyn[codon] = {}
        for p in range(3):
            for base in bases:
                tmp_codon = codon[0:p] + base + codon[p + 1 :]
                if codon_dict[codon] != codon_dict[tmp_codon]:
                    graph_nonsyn[codon][tmp_codon] = 1
                    graph[codon][tmp_codon] = 1
                else:
                    if codon != tmp_codon:
                        graph_nonsyn[codon][tmp_codon] = 0.1
                        graph[codon][tmp_codon] = 1
    for codon1 in codons:
        nonsyn_G[codon1] = {}
        G[codon1] = {}
        for codon2 in codons:
            if codon1 == codon2:
                nonsyn_G[codon1][codon2] = 0
                G[codon1][codon2] = 0
            else:
                nonsyn_G[codon1][codon2] = _dijkstra(graph_nonsyn, codon1, codon2)
                G[codon1][codon2] = _dijkstra(graph, codon1, codon2)
    return G, nonsyn_G


def _dijkstra(graph, start, end):
    """Dijkstra's algorithm Python implementation (PRIVATE).

    Algorithm adapted from
    http://thomas.pelletier.im/2010/02/dijkstras-algorithm-python-implementation/.
    However, an obvious bug in::

        if D[child_node] >(<) D[node] + child_value:

    is fixed.
    This function will return the distance between start and end.

    Arguments:
     - graph: Dictionary of dictionary (keys are vertices).
     - start: Start vertex.
     - end: End vertex.

    Output:
       List of vertices from the beginning to the end.

    """
    D = {}  # Final distances dict
    P = {}  # Predecessor dict
    # Fill the dicts with default values
    for node in graph.keys():
        D[node] = 100  # Vertices are unreachable
        P[node] = ""  # Vertices have no predecessors
    D[start] = 0  # The start vertex needs no move
    unseen_nodes = list(graph.keys())  # All nodes are unseen
    while len(unseen_nodes) > 0:
        # Select the node with the lowest value in D (final distance)
        shortest = None
        node = ""
        for temp_node in unseen_nodes:
            if shortest is None:
                shortest = D[temp_node]
                node = temp_node
            elif D[temp_node] < shortest:
                shortest = D[temp_node]
                node = temp_node
        # Remove the selected node from unseen_nodes
        unseen_nodes.remove(node)
        # For each child (ie: connected vertex) of the current node
        for child_node, child_value in graph[node].items():
            if D[child_node] > D[node] + child_value:
                D[child_node] = D[node] + child_value
                # To go to child_node, you have to go through node
                P[child_node] = node
        if node == end:
            break
    # Set a clean path
    path = []
    # We begin from the end
    node = end
    distance = 0
    # While we are not arrived at the beginning
    while not (node == start):
        if path.count(node) == 0:
            path.insert(0, node)  # Insert the predecessor of the current node
            node = P[node]  # The current node becomes its predecessor
        else:
            break
    path.insert(0, start)  # Finally, insert the start vertex
    for i in range(len(path) - 1):
        distance += graph[path[i]][path[i + 1]]
    return distance


def _count_replacement(codons, G):
    """Count replacement needed for a given codon_set (PRIVATE)."""
    if len(codons) == 1:
        return 0, 0
    elif len(codons) == 2:
        codons = list(codons)
        return floor(G[codons[0]][codons[1]])
    else:
        subgraph = {
            codon1: {codon2: G[codon1][codon2] for codon2 in codons if codon1 != codon2}
            for codon1 in codons
        }
        return _prim(subgraph)


def _prim(G):
    """Prim's algorithm to find minimum spanning tree (PRIVATE).

    Code is adapted from
    http://programmingpraxis.com/2010/04/09/minimum-spanning-tree-prims-algorithm/
    """
    nodes = []
    edges = []
    for i in G.keys():
        nodes.append(i)
        for j in G[i]:
            if (i, j, G[i][j]) not in edges and (j, i, G[i][j]) not in edges:
                edges.append((i, j, G[i][j]))
    conn = defaultdict(list)
    for n1, n2, c in edges:
        conn[n1].append((c, n1, n2))
        conn[n2].append((c, n2, n1))
    mst = []  # minimum spanning tree
    used = set(nodes[0])
    usable_edges = conn[nodes[0]][:]
    heapify(usable_edges)
    while usable_edges:
        cost, n1, n2 = heappop(usable_edges)
        if n2 not in used:
            used.add(n2)
            mst.append((n1, n2, cost))
            for e in conn[n2]:
                if e[2] not in used:
                    heappush(usable_edges, e)
    length = 0
    for p in mst:
        length += floor(p[2])
    return length


def _G_test(site_counts):
    """G test for 2x2 contingency table (PRIVATE).

    Arguments:
     - site_counts - [syn_fix, nonsyn_fix, syn_poly, nonsyn_poly]

    >>> print("%0.6f" % _G_test([17, 7, 42, 2]))
    0.004924
    """
    # TODO:
    #   Apply continuity correction for Chi-square test.
    G = 0
    tot = sum(site_counts)
    tot_syn = site_counts[0] + site_counts[2]
    tot_non = site_counts[1] + site_counts[3]
    tot_fix = sum(site_counts[:2])
    tot_poly = sum(site_counts[2:])
    exp = [
        tot_fix * tot_syn / tot,
        tot_fix * tot_non / tot,
        tot_poly * tot_syn / tot,
        tot_poly * tot_non / tot,
    ]
    for obs, ex in zip(site_counts, exp):
        G += obs * log(obs / ex)
    # with only 1 degree of freedom for a 2x2 table,
    # the cumulative chi-square distribution reduces to a simple form:
    return erfc(sqrt(G))


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
