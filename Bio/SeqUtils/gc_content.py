"""GC content calculation."""

from Bio.SeqUtils import gc_fraction


def gc_content(seq, ambiguous="remove"):
    """Calculate GC content percentage.

    Parameters
    ----------
    seq : str, Seq, or SeqRecord
        DNA or RNA sequence.
    ambiguous : {"remove", "ignore", "weighted"}, optional
        How to handle ambiguous nucleotides.
        'remove' (default) removes ambiguous nucleotides from the sequence.
        'ignore' counts only G and C; ambiguous nucleotides are not counted.
        'weighted' weights ambiguous nucleotides based on their GC content.

    Returns
    -------
    float
        GC content percentage (0-100).
    """
    return gc_fraction(seq, ambiguous) * 100.0
