# Copyright 2003-2009 by Bartek Wilczynski.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
Module containing different tools for sequence motif analysis.

it contains the core Motif class containing various I/O methods
as well as methods for motif comparisons and motif searching in sequences.
It also inlcudes functionality for parsing AlignACE and MEME programs
"""
from Bio.Motif._Motif import Motif
from Bio.Motif.Parsers.AlignAce import read as _AlignAce_read
from Bio.Motif.Parsers.MEME import read as _MEME_read
from Bio.Motif.Thresholds import ScoreDistribution

_parsers={"AlignAce" : _AlignAce_read,
          "MEME" : _MEME_read,
          }

def _from_pfm(handle):
    return Motif()._from_jaspar_pfm(handle)

def _from_sites(handle):
    return Motif()._from_jaspar_sites(handle)

_readers={"jaspar-pfm": _from_pfm,
          "jaspar-sites": _from_sites
          }


          
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
    ...     print motif.consensus()
    TCTACGATTGAG
    CTGCACCTAGCTACGAGTGAG
    GTGCCCTAAGCATACTAGGCG
    GCCACTAGCAGAGCAGGGGGC
    CGACTCAGAGGTT
    CCACGCTAAGAGAAGTGCCGGAG
    GCACGTCCCTGAGCA
    GTCCATCGCAAAGCGTGGGGC
    GAGATCAGAGGGCCG
    TGGACGCGGGG
    GACCAGAGCCTCGCATGGGGG
    AGCGCGCGTG
    GCCGGTTGCTGTTCATTAGG
    ACCGACGGCAGCTAAAAGGG
    GACGCCGGGGAT
    CGACTCGCGCTTACAAGG
    """
    try:
        parser=_parsers[format]
        
    except KeyError:
        try: #not a true parser, try reader formats
            reader=_readers[format]
        except:
            raise ValueError("Wrong parser format")
        else: #we have a proper reader 
            yield reader(handle)
    else: # we have a proper reader
        for m in parser(handle).motifs:
            yield m

def read(handle,format):
    """Reads a motif from a handle using a specified file-format.

    This supports the same formats as Bio.Motif.parse(), but
    only for files containing exactly one record.  For example,
    reading a pfm file:

    >>> from Bio import Motif
    >>> motif = Motif.read(open("Motif/SRF.pfm"),"jaspar-pfm")
    >>> motif.consensus()
    Seq('GCCCATATATGG', IUPACUnambiguousDNA())

    Or a single-motif MEME file,

    >>> from Bio import Motif
    >>> motif =  Motif.read(open("Motif/meme.out"),"MEME")
    >>> motif.consensus()
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
    >>> motif.consensus()
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


def _test():
    """Run the Bio.Motif module's doctests.

    This will try and locate the unit tests directory, and run the doctests
    from there in order that the relative paths used in the examples work.
    """
    import doctest
    import os
    if os.path.isdir(os.path.join("..","..","Tests")):
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join("..","..","Tests"))
        doctest.testmod()
        os.chdir(cur_dir)
        del cur_dir
        print "Done"

if __name__ == "__main__":
    #Run the doctests
    _test()
