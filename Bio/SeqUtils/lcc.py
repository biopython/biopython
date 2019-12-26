# Copyright 2003, 2007 by Sebastian Bassi. sbassi@genesdigitales.com
# All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Local Composition Complexity."""

import math


def lcc_mult(seq, wsize):
    """Calculate Local Composition Complexity (LCC) values over sliding window.

    Returns a list of floats, the LCC values for a sliding window over
    the sequence.

    seq - an unambiguous DNA sequence (a string or Seq object)
    wsize - window size, integer

    The result is the same as applying lcc_simp multiple times, but this
    version is optimized for speed. The optimization works by using the
    value of previous window as a base to compute the next one.
    """
    l2 = math.log(2)
    tamseq = len(seq)
    try:
        # Assume its a string
        upper = seq.upper()
    except AttributeError:
        # Should be a Seq object then
        upper = str(seq).upper()
    compone = [0]
    lccsal = [0]
    for i in range(wsize):
        compone.append(
            ((i + 1) / float(wsize)) * ((math.log((i + 1) / float(wsize))) / l2)
        )
    window = seq[0:wsize]
    cant_a = window.count("A")
    cant_c = window.count("C")
    cant_t = window.count("T")
    cant_g = window.count("G")
    term_a = compone[cant_a]
    term_c = compone[cant_c]
    term_t = compone[cant_t]
    term_g = compone[cant_g]
    lccsal.append(-(term_a + term_c + term_t + term_g))
    tail = seq[0]
    for x in range(tamseq - wsize):
        window = upper[x + 1 : wsize + x + 1]
        if tail == window[-1]:
            lccsal.append(lccsal[-1])
        elif tail == "A":
            cant_a -= 1
            if window.endswith("C"):
                cant_c += 1
                term_a = compone[cant_a]
                term_c = compone[cant_c]
                lccsal.append(-(term_a + term_c + term_t + term_g))
            elif window.endswith("T"):
                cant_t += 1
                term_a = compone[cant_a]
                term_t = compone[cant_t]
                lccsal.append(-(term_a + term_c + term_t + term_g))
            elif window.endswith("G"):
                cant_g += 1
                term_a = compone[cant_a]
                term_g = compone[cant_g]
                lccsal.append(-(term_a + term_c + term_t + term_g))
        elif tail == "C":
            cant_c -= 1
            if window.endswith("A"):
                cant_a += 1
                term_a = compone[cant_a]
                term_c = compone[cant_c]
                lccsal.append(-(term_a + term_c + term_t + term_g))
            elif window.endswith("T"):
                cant_t += 1
                term_c = compone[cant_c]
                term_t = compone[cant_t]
                lccsal.append(-(term_a + term_c + term_t + term_g))
            elif window.endswith("G"):
                cant_g += 1
                term_c = compone[cant_c]
                term_g = compone[cant_g]
                lccsal.append(-(term_a + term_c + term_t + term_g))
        elif tail == "T":
            cant_t -= 1
            if window.endswith("A"):
                cant_a += 1
                term_a = compone[cant_a]
                term_t = compone[cant_t]
                lccsal.append(-(term_a + term_c + term_t + term_g))
            elif window.endswith("C"):
                cant_c += 1
                term_c = compone[cant_c]
                term_t = compone[cant_t]
                lccsal.append(-(term_a + term_c + term_t + term_g))
            elif window.endswith("G"):
                cant_g += 1
                term_t = compone[cant_t]
                term_g = compone[cant_g]
                lccsal.append(-(term_a + term_c + term_t + term_g))
        elif tail == "G":
            cant_g -= 1
            if window.endswith("A"):
                cant_a += 1
                term_a = compone[cant_a]
                term_g = compone[cant_g]
                lccsal.append(-(term_a + term_c + term_t + term_g))
            elif window.endswith("C"):
                cant_c += 1
                term_c = compone[cant_c]
                term_g = compone[cant_g]
                lccsal.append(-(term_a + term_c + term_t + term_g))
            elif window.endswith("T"):
                cant_t += 1
                term_t = compone[cant_t]
                term_g = compone[cant_g]
                lccsal.append(-(term_a + term_c + term_t + term_g))
        tail = window[0]
    return lccsal


def lcc_simp(seq):
    """Calculate Local Composition Complexity (LCC) for a sequence.

    seq - an unambiguous DNA sequence (a string or Seq object)

    Returns the Local Composition Complexity (LCC) value for the entire
    sequence (as a float).

    Reference:
    Andrzej K Konopka (2005) Sequence Complexity and Composition
    https://doi.org/10.1038/npg.els.0005260
    """
    wsize = len(seq)
    try:
        # Assume its a string
        upper = seq.upper()
    except AttributeError:
        # Should be a Seq object then
        upper = str(seq).upper()
    l2 = math.log(2)
    if "A" not in seq:
        term_a = 0
        # Check to avoid calculating the log of 0.
    else:
        term_a = ((upper.count("A")) / float(wsize)) * (
            (math.log((upper.count("A")) / float(wsize))) / l2
        )
    if "C" not in seq:
        term_c = 0
    else:
        term_c = ((upper.count("C")) / float(wsize)) * (
            (math.log((upper.count("C")) / float(wsize))) / l2
        )
    if "T" not in seq:
        term_t = 0
    else:
        term_t = ((upper.count("T")) / float(wsize)) * (
            (math.log((upper.count("T")) / float(wsize))) / l2
        )
    if "G" not in seq:
        term_g = 0
    else:
        term_g = ((upper.count("G")) / float(wsize)) * (
            (math.log((upper.count("G")) / float(wsize))) / l2
        )
    return -(term_a + term_c + term_t + term_g)
