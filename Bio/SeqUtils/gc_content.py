"""GC content calculation for sequences."""

def gc_content(seq, ambiguous="remove"):
    """Calculate G+C percentage in seq (float between 0 and 100).

    This is a convenience wrapper around gc_fraction that returns
    percentage instead of fraction.

    Parameters
    ----------
    seq : str or Seq object
        DNA or RNA sequence.
    ambiguous : {"remove", "ignore", "weighted"}, optional
        How to handle ambiguous nucleotides. Default is "remove".
        See gc_fraction for details.

    Returns
    -------
    float
        G+C percentage (0 to 100).

    Examples
    --------
    >>> from Bio.SeqUtils import gc_content
    >>> gc_content("ACTG")
    50.0
    >>> gc_content("ACTG", ambiguous="weighted")
    50.0
    """
    from . import gc_fraction
    return gc_fraction(seq, ambiguous) * 100.0