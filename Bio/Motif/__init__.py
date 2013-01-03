# Copyright 2003-2009 by Bartek Wilczynski.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
Module containing different tools for sequence motif analysis.

It contains the core Motif class containing various I/O methods
as well as methods for motif comparisons and motif searching in sequences.
It also includes functionality for parsing AlignACE and MEME programs.
"""
from Bio.Motif._Motif import Motif
from Bio.Motif.AlignAce import read as _AlignAce_read
from Bio.Motif.MEME import read as _MEME_read
from Bio.Motif import Jaspar
from Bio.Motif.Thresholds import ScoreDistribution
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq


def _from_pfm(handle):
    return Motif()._from_jaspar_pfm(handle)


def _from_sites(handle):
    return Motif()._from_jaspar_sites(handle)


def create(instances, alphabet=None):
    for instance in instances:
        try:
            a = instance.alphabet
        except AttributeError:
            # The instance is a plain string
            continue
        if alphabet is None:
            alphabet = a
        elif alphabet != a:
            raise ValueError("Alphabets are inconsistent")
    if alphabet is None or alphabet.letters is None:
        # If we didn't get a meaningful alphabet from the instances,
        # assume it is DNA.
        alphabet = IUPAC.unambiguous_dna
    seqs = []
    for instance in instances:
        seq = Seq(str(instance), alphabet=alphabet)
        seqs.append(seq)
    return Motif(instances=seqs, alphabet=alphabet)


def parse(handle,format):
    """Parses an output file of motif finding programs.

    Currently supported formats:
     - AlignAce
     - MEME

    You can also use single-motif formats, although the Bio.Motif.read()
    function is simpler to use in this situation.
     - jaspar-pfm
     - jaspar-sites

    For example:

    >>> from Bio import Motif
    >>> for motif in Motif.parse(open("Motif/alignace.out"),"AlignAce"):
    ...     print motif.consensus
    TCTACGATTGAG
    CTGCAGCTAGCTACGAGTGAG
    GTGCTCTAAGCATAGTAGGCG
    GCCACTAGCAGAGCAGGGGGC
    CGACTCAGAGGTT
    CCACGCTAAGAGAGGTGCCGGAG
    GCGCGTCGCTGAGCA
    GTCCATCGCAAAGCGTGGGGC
    GGGATCAGAGGGCCG
    TGGAGGCGGGG
    GACCAGAGCTTCGCATGGGGG
    GGCGTGCGTG
    GCTGGTTGCTGTTCATTAGG
    GCCGGCGGCAGCTAAAAGGG
    GAGGCCGGGGAT
    CGACTCGTGCTTAGAAGG
    """
    if format in ('pfm', 'sites'):
        yield Jaspar.read(handle, format)
    elif format=="AlignAce":
        record = _AlignAce_read(handle)
        for m in record.motifs:
            yield m
    elif format=="MEME":
        record = _MEME_read(handle)
        for m in record.motifs:
            yield m
    elif format=="jaspar-pfm":
        yield _from_pfm(handle)
    elif format=="jaspar-sites":
        yield _from_sites(handle)
    else:
        raise ValueError("Unknown format %s" % format)


def read(handle,format):
    """Reads a motif from a handle using a specified file-format.

    This supports the same formats as Bio.Motif.parse(), but
    only for files containing exactly one record.  For example,
    reading a pfm file:

    >>> from Bio import Motif
    >>> motif = Motif.read(open("Motif/SRF.pfm"), "pfm")
    >>> motif.consensus
    Seq('GCCCATATATGG', IUPACUnambiguousDNA())

    Or a single-motif MEME file,

    >>> from Bio import Motif
    >>> motif =  Motif.read(open("Motif/meme.out"),"MEME")
    >>> motif.consensus
    Seq('CTCAATCGTA', IUPACUnambiguousDNA())

    If the handle contains no records, or more than one record,
    an exception is raised:

    >>> from Bio import Motif
    >>> motif = Motif.read(open("Motif/alignace.out"),"AlignAce")
    Traceback (most recent call last):
        ...
    ValueError: More than one motif found in handle

    If however you want the first record from a file containing
    multiple records this function would raise an exception (as
    shown in the example above).  Instead use:

    >>> from Bio import Motif
    >>> motif = Motif.parse(open("Motif/alignace.out"),"AlignAce").next()
    >>> motif.consensus
    Seq('TCTACGATTGAG', IUPACUnambiguousDNA())

    Use the Bio.Motif.parse(handle, format) function if you want
    to read multiple records from the handle.
    """
    iterator = parse(handle, format)
    try:
        first = iterator.next()
    except StopIteration:
        first = None
    if first is None:
        raise ValueError("No motifs found in handle")
    try:
        second = iterator.next()
    except StopIteration:
        second = None
    if second is not None:
        raise ValueError("More than one motif found in handle")
    return first


if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()
