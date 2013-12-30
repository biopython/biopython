#!/usr/bin/env python
# Created: Wed May 29 08:07:18 2002
# thomas@cbs.dtu.dk, Cecilia.Alsmark@ebc.uu.se
# Copyright 2001 by Thomas Sicheritz-Ponten and Cecilia Alsmark.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Miscellaneous functions for dealing with sequences."""

from __future__ import print_function

import re
from math import pi, sin, cos

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data import IUPACData


######################################
# DNA
######################
# {{{


def GC(seq):
    """Calculates G+C content, returns the percentage (float between 0 and 100).

    Copes mixed case sequences, and with the ambiguous nucleotide S (G or C)
    when counting the G and C content.  The percentage is calculated against
    the full length, e.g.:

    >>> from Bio.SeqUtils import GC
    >>> GC("ACTGN")
    40.0

    Note that this will return zero for an empty sequence.
    """
    try:
        gc = sum(seq.count(x) for x in ['G', 'C', 'g', 'c', 'S', 's'])
        return gc*100.0/len(seq)
    except ZeroDivisionError:
        return 0.0


def GC123(seq):
    """Calculates total G+C content plus first, second and third positions.

    Returns a tuple of four floats (percentages between 0 and 100) for the
    entire sequence, and the three codon positions.  e.g.

    >>> from Bio.SeqUtils import GC123
    >>> GC123("ACTGTN")
    (40.0, 50.0, 50.0, 0.0)

    Copes with mixed case sequences, but does NOT deal with ambiguous
    nucleotides.
    """
    d= {}
    for nt in ['A', 'T', 'G', 'C']:
        d[nt] = [0, 0, 0]

    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        if len(codon) < 3:
            codon += '  '
        for pos in range(0, 3):
            for nt in ['A', 'T', 'G', 'C']:
                if codon[pos] == nt or codon[pos] == nt.lower():
                    d[nt][pos] += 1
    gc = {}
    gcall = 0
    nall = 0
    for i in range(0, 3):
        try:
            n = d['G'][i] + d['C'][i] +d['T'][i] + d['A'][i]
            gc[i] = (d['G'][i] + d['C'][i])*100.0/n
        except:
            gc[i] = 0

        gcall = gcall + d['G'][i] + d['C'][i]
        nall = nall + n

    gcall = 100.0*gcall/nall
    return gcall, gc[0], gc[1], gc[2]


def GC_skew(seq, window=100):
    """Calculates GC skew (G-C)/(G+C) for multiple windows along the sequence.

    Returns a list of ratios (floats), controlled by the length of the sequence
    and the size of the window.

    Does NOT look at any ambiguous nucleotides.
    """
    # 8/19/03: Iddo: added lowercase
    values = []
    for i in range(0, len(seq), window):
        s = seq[i: i + window]
        g = s.count('G') + s.count('g')
        c = s.count('C') + s.count('c')
        skew = (g-c)/float(g+c)
        values.append(skew)
    return values


def xGC_skew(seq, window=1000, zoom=100,
                         r=300, px=100, py=100):
    """Calculates and plots normal and accumulated GC skew (GRAPHICS !!!)."""
    try:
        import Tkinter as tkinter # Python 2
    except ImportError:
        import tkinter # Python 3

    yscroll = tkinter.Scrollbar(orient=tkinter.VERTICAL)
    xscroll = tkinter.Scrollbar(orient=tkinter.HORIZONTAL)
    canvas = tkinter.Canvas(yscrollcommand=yscroll.set,
                            xscrollcommand=xscroll.set, background='white')
    win = canvas.winfo_toplevel()
    win.geometry('700x700')

    yscroll.config(command=canvas.yview)
    xscroll.config(command=canvas.xview)
    yscroll.pack(side=tkinter.RIGHT, fill=tkinter.Y)
    xscroll.pack(side=tkinter.BOTTOM, fill=tkinter.X)
    canvas.pack(fill=tkinter.BOTH, side=tkinter.LEFT, expand=1)
    canvas.update()

    X0, Y0 = r + px, r + py
    x1, x2, y1, y2 = X0 - r, X0 + r, Y0 - r, Y0 + r

    ty = Y0
    canvas.create_text(X0, ty, text='%s...%s (%d nt)' % (seq[:7], seq[-7:], len(seq)))
    ty += 20
    canvas.create_text(X0, ty, text='GC %3.2f%%' % (GC(seq)))
    ty += 20
    canvas.create_text(X0, ty, text='GC Skew', fill='blue')
    ty += 20
    canvas.create_text(X0, ty, text='Accumulated GC Skew', fill='magenta')
    ty += 20
    canvas.create_oval(x1, y1, x2, y2)

    acc = 0
    start = 0
    for gc in GC_skew(seq, window):
        r1 = r
        acc += gc
        # GC skew
        alpha = pi - (2*pi*start)/len(seq)
        r2 = r1 - gc*zoom
        x1 = X0 + r1 * sin(alpha)
        y1 = Y0 + r1 * cos(alpha)
        x2 = X0 + r2 * sin(alpha)
        y2 = Y0 + r2 * cos(alpha)
        canvas.create_line(x1, y1, x2, y2, fill='blue')
        # accumulated GC skew
        r1 = r - 50
        r2 = r1 - acc
        x1 = X0 + r1 * sin(alpha)
        y1 = Y0 + r1 * cos(alpha)
        x2 = X0 + r2 * sin(alpha)
        y2 = Y0 + r2 * cos(alpha)
        canvas.create_line(x1, y1, x2, y2, fill='magenta')

        canvas.update()
        start += window

    canvas.configure(scrollregion=canvas.bbox(tkinter.ALL))


def molecular_weight(seq):
    """Calculate the molecular weight of a DNA sequence."""
    if isinstance(seq, str):
        seq = Seq(seq, IUPAC.unambiguous_dna)
    weight_table = IUPACData.unambiguous_dna_weights
    return sum(weight_table[x] for x in seq)


def nt_search(seq, subseq):
    """Search for a DNA subseq in sequence.

    use ambiguous values (like N = A or T or C or G, R = A or G etc.)
    searches only on forward strand
    """
    pattern = ''
    for nt in subseq:
        value = IUPACData.ambiguous_dna_values[nt]
        if len(value) == 1:
            pattern += value
        else:
            pattern += '[%s]' % value

    pos = -1
    result = [pattern]
    l = len(seq)
    while True:
        pos += 1
        s = seq[pos:]
        m = re.search(pattern, s)
        if not m:
            break
        pos += int(m.start(0))
        result.append(pos)
    return result

# }}}

######################################
# Protein
######################
# {{{

def seq3(seq, custom_map={'*': 'Ter'}, undef_code='Xaa'):
    """Turn a one letter code protein sequence into one with three letter codes.

    The single input argument 'seq' should be a protein sequence using single
    letter codes, either as a python string or as a Seq or MutableSeq object.

    This function returns the amino acid sequence as a string using the three
    letter amino acid codes. Output follows the IUPAC standard (including
    ambiguous characters B for "Asx", J for "Xle" and X for "Xaa", and also U
    for "Sel" and O for "Pyl") plus "Ter" for a terminator given as an asterisk.
    Any unknown character (including possible gap characters), is changed into
    'Xaa'.

    e.g.

    >>> from Bio.SeqUtils import seq3
    >>> seq3("MAIVMGRWKGAR*")
    'MetAlaIleValMetGlyArgTrpLysGlyAlaArgTer'

    You can set a custom translation of the codon termination code using the
    "custom_map" argument, e.g.

    >>> seq3("MAIVMGRWKGAR*", custom_map={"*": "***"})
    'MetAlaIleValMetGlyArgTrpLysGlyAlaArg***'

    You can also set a custom translation for non-amino acid characters, such
    as '-', using the "undef_code" argument, e.g.

    >>> seq3("MAIVMGRWKGA--R*", undef_code='---')
    'MetAlaIleValMetGlyArgTrpLysGlyAla------ArgTer'

    If not given, "undef_code" defaults to "Xaa", e.g.

    >>> seq3("MAIVMGRWKGA--R*")
    'MetAlaIleValMetGlyArgTrpLysGlyAlaXaaXaaArgTer'

    This function was inspired by BioPerl's seq3.
    """
    # not doing .update() on IUPACData dict with custom_map dict
    # to preserve its initial state (may be imported in other modules)
    threecode = dict(list(IUPACData.protein_letters_1to3_extended.items()) +
                     list(custom_map.items()))
    #We use a default of 'Xaa' for undefined letters
    #Note this will map '-' to 'Xaa' which may be undesirable!
    return ''.join(threecode.get(aa, undef_code) for aa in seq)


def seq1(seq, custom_map={'Ter': '*'}, undef_code='X'):
    """Turns a three-letter code protein sequence into one with single letter codes.

    The single input argument 'seq' should be a protein sequence using three-
    letter codes, either as a python string or as a Seq or MutableSeq object.

    This function returns the amino acid sequence as a string using the one
    letter amino acid codes. Output follows the IUPAC standard (including
    ambiguous characters "B" for "Asx", "J" for "Xle", "X" for "Xaa", "U" for
    "Sel", and "O" for "Pyl") plus "*" for a terminator given the "Ter" code.
    Any unknown character (including possible gap characters), is changed into
    '-'.

    e.g.

    >>> from Bio.SeqUtils import seq3
    >>> seq1("MetAlaIleValMetGlyArgTrpLysGlyAlaArgTer")
    'MAIVMGRWKGAR*'

    The input is case insensitive, e.g.

    >>> from Bio.SeqUtils import seq3
    >>> seq1("METalaIlEValMetGLYArgtRplysGlyAlaARGTer")
    'MAIVMGRWKGAR*'

    You can set a custom translation of the codon termination code using the
    "custom_map" argument, e.g.

    >>> seq1("MetAlaIleValMetGlyArgTrpLysGlyAlaArg***", custom_map={"***": "*"})
    'MAIVMGRWKGAR*'

    You can also set a custom translation for non-amino acid characters, such
    as '-', using the "undef_code" argument, e.g.

    >>> seq1("MetAlaIleValMetGlyArgTrpLysGlyAla------ArgTer", undef_code='?')
    'MAIVMGRWKGA??R*'

    If not given, "undef_code" defaults to "X", e.g.

    >>> seq1("MetAlaIleValMetGlyArgTrpLysGlyAla------ArgTer")
    'MAIVMGRWKGAXXR*'

    """
    # reverse map of threecode
    # upper() on all keys to enable caps-insensitive input seq handling
    onecode = dict((k.upper(), v) for k, v in
                   IUPACData.protein_letters_3to1_extended.items())
    # add the given termination codon code and custom maps
    onecode.update((k.upper(), v) for (k, v) in custom_map.items())
    seqlist = [seq[3*i:3*(i+1)] for i in range(len(seq) // 3)]
    return ''.join(onecode.get(aa.upper(), undef_code) for aa in seqlist)


# }}}

######################################
# Mixed ???
######################
# {{{


def six_frame_translations(seq, genetic_code=1):
    """Formatted string showing the 6 frame translations and GC content.

    nice looking 6 frame translation with GC content - code from xbbtools
    similar to DNA Striders six-frame translation

    >>> from Bio.SeqUtils import six_frame_translations
    >>> print(six_frame_translations("AUGGCCAUUGUAAUGGGCCGCUGA"))
    GC_Frame: a:5 t:0 g:8 c:5 
    Sequence: auggccauug ... gggccgcuga, 24 nt, 54.17 %GC
    <BLANKLINE>
    <BLANKLINE>
    1/1
      G  H  C  N  G  P  L
     W  P  L  *  W  A  A
    M  A  I  V  M  G  R  *
    auggccauuguaaugggccgcuga   54 %
    uaccgguaacauuacccggcgacu
    A  M  T  I  P  R  Q 
     H  G  N  Y  H  A  A  S
      P  W  Q  L  P  G  S
    <BLANKLINE>
    <BLANKLINE>

    """
    from Bio.Seq import reverse_complement, translate
    anti = reverse_complement(seq)
    comp = anti[::-1]
    length = len(seq)
    frames = {}
    for i in range(0, 3):
        fragment_length = 3 * ((length-i) // 3)
        frames[i+1] = translate(seq[i:i+fragment_length], genetic_code)
        frames[-(i+1)] = translate(anti[i:i+fragment_length], genetic_code)[::-1]

    # create header
    if length > 20:
        short = '%s ... %s' % (seq[:10], seq[-10:])
    else:
        short = seq
    header = 'GC_Frame: '
    for nt in ['a', 't', 'g', 'c']:
        header += '%s:%d ' % (nt, seq.count(nt.upper()))

    header += '\nSequence: %s, %d nt, %0.2f %%GC\n\n\n' % (short.lower(), length, GC(seq))
    res = header

    for i in range(0, length, 60):
        subseq = seq[i:i+60]
        csubseq = comp[i:i+60]
        p = i//3
        res += '%d/%d\n' % (i+1, i/3+1)
        res += '  ' + '  '.join(frames[3][p:p+20]) + '\n'
        res += ' ' + '  '.join(frames[2][p:p+20]) + '\n'
        res += '  '.join(frames[1][p:p+20]) + '\n'
        # seq
        res += subseq.lower() + '%5d %%\n' % int(GC(subseq))
        res += csubseq.lower() + '\n'
        # - frames
        res += '  '.join(frames[-2][p:p+20]) +' \n'
        res += ' ' + '  '.join(frames[-1][p:p+20]) + '\n'
        res += '  ' + '  '.join(frames[-3][p:p+20]) + '\n\n'
    return res

# }}}

######################################
# FASTA file utilities
######################
# {{{


def quick_FASTA_reader(file):
    """Simple FASTA reader, returning a list of string tuples (OBSOLETE).

    The single argument 'file' should be the filename of a FASTA format file.
    This function will open and read in the entire file, constructing a list
    of all the records, each held as a tuple of strings (the sequence name or
    title, and its sequence).

    >>> seqs = quick_FASTA_reader("Fasta/dups.fasta")
    >>> for title, sequence in seqs:
    ...     print("%s %s" % (title, sequence))
    alpha ACGTA
    beta CGTC
    gamma CCGCC
    alpha (again - this is a duplicate entry to test the indexing code) ACGTA
    delta CGCGC

    This function was is fast, but because it returns the data as a single in
    memory list, is unsuitable for large files where an iterator approach is
    preferable.

    You are generally encouraged to use Bio.SeqIO.parse(handle, "fasta") which
    allows you to iterate over the records one by one (avoiding having all the
    records in memory at once).  Using Bio.SeqIO also makes it easy to switch
    between different input file formats.  However, please note that rather
    than simple strings, Bio.SeqIO uses SeqRecord objects for each record.

    If you want to use simple strings, use the function SimpleFastaParser
    added to Bio.SeqIO.FastaIO in Biopython 1.61 instead.
    """
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    with open(file) as handle:
        entries = list(SimpleFastaParser(handle))
    return entries


# }}}


def _test():
    """Run the module's doctests (PRIVATE)."""
    import os
    import doctest
    if os.path.isdir(os.path.join("..", "Tests")):
        print("Running doctests...")
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join("..", "Tests"))
        doctest.testmod()
        os.chdir(cur_dir)
        del cur_dir
        print("Done")
    elif os.path.isdir(os.path.join("Tests")):
        print("Running doctests...")
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join("Tests"))
        doctest.testmod()
        os.chdir(cur_dir)
        del cur_dir
        print("Done")


if __name__ == "__main__":
    _test()
