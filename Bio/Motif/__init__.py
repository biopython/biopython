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


def create(instances, alphabet=None):
    from Bio.Alphabet import IUPAC
    from Bio.Seq import Seq
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
     - AlignAce:      AlignAce output file format
     - MEME:          MEME output file motif
     - TRANSFAC:      TRANSFAC database file format
     - pfm:           JASPAR-style position-frequency matrix
     - sites:         JASPAR-style sites file
     - jaspar-pfm:    JASPAR-style position-frequency matrix [DEPRECATED]
     - jaspar-sites:  JASPAR-style sites file [DEPRECATED]
    As files in the pfm and sites formats contain only a single motif,
    it is easier to use Bio.Motif.read() instead of Bio.Motif.parse()
    for those.

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
    if format=="AlignAce":
        from Bio.Motif import AlignAce
        record = AlignAce.read(handle)
        return record
    elif format=="MEME":
        from Bio.Motif import MEME
        record = MEME.read(handle)
        return record
    elif format=="TRANSFAC":
        from Bio.Motif import TRANSFAC
        record = TRANSFAC.read(handle)
        return record
    elif format in ('pfm', 'sites'):
        from Bio.Motif import Jaspar
        motif = Jaspar.read(handle, format)
    elif format=="jaspar-pfm":
        motif = Motif()._from_jaspar_pfm(handle)
    elif format=="jaspar-sites":
        motif = Motif()._from_jaspar_sites(handle)
    else:
        raise ValueError("Unknown format %s" % format)
    # Treat the single-motif formats
    motifs = [motif]
    return motifs


def read(handle,format):
    """Reads a motif from a handle using a specified file-format.

    This supports the same formats as Bio.Motif.parse(), but
    only for files containing exactly one motif.  For example,
    reading a JASPAR-style pfm file:

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

    If however you want the first motif from a file containing
    multiple motifs this function would raise an exception (as
    shown in the example above).  Instead use:

    >>> from Bio import Motif
    >>> motifs = Motif.parse(open("Motif/alignace.out"),"AlignAce")
    >>> motif = motifs[0]
    >>> motif.consensus
    Seq('TCTACGATTGAG', IUPACUnambiguousDNA())

    Use the Bio.Motif.parse(handle, format) function if you want
    to read multiple records from the handle.
    """
    motifs = parse(handle, format)
    if len(motifs)==0:
        raise ValueError("No motifs found in handle")
    if len(motifs) > 1:
        raise ValueError("More than one motif found in handle")
    motif = motifs[0]
    return motif


if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()
