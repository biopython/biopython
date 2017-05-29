#!/usr/bin/env python
# Created: Wed May 29 08:07:18 2002
# thomas@cbs.dtu.dk, Cecilia.Alsmark@ebc.uu.se
# Copyright 2001 by Thomas Sicheritz-Ponten and Cecilia Alsmark.
# Revisions copyright 2014 by Markus Piotrowski.
# Revisions copyright 2014-2016 by Peter Cock.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Miscellaneous functions for dealing with sequences."""

from __future__ import print_function

import re
from math import pi, sin, cos
from random import random, randint, seed

from Bio.Seq import Seq, MutableSeq
from Bio import Alphabet
from Bio.Alphabet import IUPAC
from Bio.Data import IUPACData, CodonTable

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
    gc = sum(seq.count(x) for x in ['G', 'C', 'g', 'c', 'S', 's'])
    try:
        return gc * 100.0 / len(seq)
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
    d = {}
    for nt in ['A', 'T', 'G', 'C']:
        d[nt] = [0, 0, 0]

    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
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
            n = d['G'][i] + d['C'][i] + d['T'][i] + d['A'][i]
            gc[i] = (d['G'][i] + d['C'][i]) * 100.0 / n
        except Exception:  # TODO - ValueError?
            gc[i] = 0

        gcall = gcall + d['G'][i] + d['C'][i]
        nall = nall + n

    gcall = 100.0 * gcall / nall
    return gcall, gc[0], gc[1], gc[2]


def GC_skew(seq, window=100):
    """Calculates GC skew (G-C)/(G+C) for multiple windows along the sequence.

    Returns a list of ratios (floats), controlled by the length of the sequence
    and the size of the window.

    Returns 0 for windows without any G/C by handling zero division errors.

    Does NOT look at any ambiguous nucleotides.
    """
    # 8/19/03: Iddo: added lowercase
    values = []
    for i in range(0, len(seq), window):
        s = seq[i: i + window]
        g = s.count('G') + s.count('g')
        c = s.count('C') + s.count('c')
        try:
            skew = (g - c) / float(g + c)
        except ZeroDivisionError:
            skew = 0.0
        values.append(skew)
    return values


def xGC_skew(seq, window=1000, zoom=100,
                         r=300, px=100, py=100):
    """Calculates and plots normal and accumulated GC skew (GRAPHICS !!!)."""
    try:
        import Tkinter as tkinter  # Python 2
    except ImportError:
        import tkinter  # Python 3

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
        alpha = pi - (2 * pi * start) / len(seq)
        r2 = r1 - gc * zoom
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


def seq3(seq, custom_map=None, undef_code='Xaa'):
    """Turn a one letter code protein sequence into one with three letter codes.

    The single required input argument 'seq' should be a protein sequence using
    single letter codes, either as a python string or as a Seq or MutableSeq
    object.

    This function returns the amino acid sequence as a string using the three
    letter amino acid codes. Output follows the IUPAC standard (including
    ambiguous characters B for "Asx", J for "Xle" and X for "Xaa", and also U
    for "Sel" and O for "Pyl") plus "Ter" for a terminator given as an asterisk.
    Any unknown character (including possible gap characters), is changed into
    'Xaa' by default.

    e.g.

    >>> from Bio.SeqUtils import seq3
    >>> seq3("MAIVMGRWKGAR*")
    'MetAlaIleValMetGlyArgTrpLysGlyAlaArgTer'

    You can set a custom translation of the codon termination code using the
    dictionary "custom_map" argument (which defaults to {'*': 'Ter'}), e.g.

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
    if custom_map is None:
        custom_map = {'*': 'Ter'}
    # not doing .update() on IUPACData dict with custom_map dict
    # to preserve its initial state (may be imported in other modules)
    threecode = dict(list(IUPACData.protein_letters_1to3_extended.items()) +
                     list(custom_map.items()))
    # We use a default of 'Xaa' for undefined letters
    # Note this will map '-' to 'Xaa' which may be undesirable!
    return ''.join(threecode.get(aa, undef_code) for aa in seq)


def seq1(seq, custom_map=None, undef_code='X'):
    """Turns a three-letter code protein sequence into one with single letter codes.

    The single required input argument 'seq' should be a protein sequence
    using three-letter codes, either as a python string or as a Seq or
    MutableSeq object.

    This function returns the amino acid sequence as a string using the one
    letter amino acid codes. Output follows the IUPAC standard (including
    ambiguous characters "B" for "Asx", "J" for "Xle", "X" for "Xaa", "U" for
    "Sel", and "O" for "Pyl") plus "*" for a terminator given the "Ter" code.
    Any unknown character (including possible gap characters), is changed
    into '-' by default.

    e.g.

    >>> from Bio.SeqUtils import seq3
    >>> seq1("MetAlaIleValMetGlyArgTrpLysGlyAlaArgTer")
    'MAIVMGRWKGAR*'

    The input is case insensitive, e.g.

    >>> from Bio.SeqUtils import seq3
    >>> seq1("METalaIlEValMetGLYArgtRplysGlyAlaARGTer")
    'MAIVMGRWKGAR*'

    You can set a custom translation of the codon termination code using the
    dictionary "custom_map" argument (defaulting to {'Ter': '*'}), e.g.

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
    if custom_map is None:
        custom_map = {'Ter': '*'}
    # reverse map of threecode
    # upper() on all keys to enable caps-insensitive input seq handling
    onecode = dict((k.upper(), v) for k, v in
                   IUPACData.protein_letters_3to1_extended.items())
    # add the given termination codon code and custom maps
    onecode.update((k.upper(), v) for (k, v) in custom_map.items())
    seqlist = [seq[3 * i:3 * (i + 1)] for i in range(len(seq) // 3)]
    return ''.join(onecode.get(aa.upper(), undef_code) for aa in seqlist)


# }}}

######################################
# Mixed ???
######################
# {{{


def molecular_weight(seq, seq_type=None, double_stranded=False, circular=False,
                     monoisotopic=False):
    """Calculates the molecular weight of a DNA, RNA or protein sequence.

    Only unambiguous letters are allowed. Nucleotide sequences are assumed to
    have a 5' phosphate.

        - seq: String or Biopython sequence object.
        - seq_type: The default (None) is to take the alphabet from the seq argument,
          or assume DNA if the seq argument is a string. Override this with
          a string 'DNA', 'RNA', or 'protein'.
        - double_stranded: Calculate the mass for the double stranded molecule?
        - circular: Is the molecule circular (has no ends)?
        - monoisotopic: Use the monoisotopic mass tables?

    Note that for backwards compatibility, if the seq argument is a string,
    or Seq object with a generic alphabet, and no seq_type is specified
    (i.e. left as None), then DNA is assumed.

    >>> print("%0.2f" % molecular_weight("AGC"))
    949.61
    >>> print("%0.2f" % molecular_weight(Seq("AGC")))
    949.61

    However, it is better to be explicit - for example with strings:

    >>> print("%0.2f" % molecular_weight("AGC", "DNA"))
    949.61
    >>> print("%0.2f" % molecular_weight("AGC", "RNA"))
    997.61
    >>> print("%0.2f" % molecular_weight("AGC", "protein"))
    249.29

    Or, with the sequence alphabet:

    >>> from Bio.Seq import Seq
    >>> from Bio.Alphabet import generic_dna, generic_rna, generic_protein
    >>> print("%0.2f" % molecular_weight(Seq("AGC", generic_dna)))
    949.61
    >>> print("%0.2f" % molecular_weight(Seq("AGC", generic_rna)))
    997.61
    >>> print("%0.2f" % molecular_weight(Seq("AGC", generic_protein)))
    249.29

    Also note that contradictory sequence alphabets and seq_type will also
    give an exception:

    >>> from Bio.Seq import Seq
    >>> from Bio.Alphabet import generic_dna
    >>> print("%0.2f" % molecular_weight(Seq("AGC", generic_dna), "RNA"))
    Traceback (most recent call last):
      ...
    ValueError: seq_type='RNA' contradicts DNA from seq alphabet

    """
    # Rewritten by Markus Piotrowski, 2014

    # Find the alphabet type
    tmp_type = ''
    if isinstance(seq, (Seq, MutableSeq)):
        base_alphabet = Alphabet._get_base_alphabet(seq.alphabet)
        if isinstance(base_alphabet, Alphabet.DNAAlphabet):
            tmp_type = 'DNA'
        elif isinstance(base_alphabet, Alphabet.RNAAlphabet):
            tmp_type = 'RNA'
        elif isinstance(base_alphabet, Alphabet.ProteinAlphabet):
            tmp_type = 'protein'
        elif isinstance(base_alphabet, Alphabet.ThreeLetterProtein):
            tmp_type = 'protein'
            # Convert to one-letter sequence. Have to use a string for seq1
            seq = Seq(seq1(str(seq)), alphabet=Alphabet.ProteinAlphabet())
        elif not isinstance(base_alphabet, Alphabet.Alphabet):
            raise TypeError("%s is not a valid alphabet for mass calculations"
                             % base_alphabet)
        else:
            tmp_type = "DNA"  # backward compatibity
        if seq_type and tmp_type and tmp_type != seq_type:
            raise ValueError("seq_type=%r contradicts %s from seq alphabet"
                             % (seq_type, tmp_type))
        seq_type = tmp_type
    elif isinstance(seq, str):
        if seq_type is None:
            seq_type = "DNA"  # backward compatibity
    else:
        raise TypeError("Expected a string or Seq object, not seq=%r" % seq)

    seq = ''.join(str(seq).split()).upper()  # Do the minimum formatting

    if seq_type == 'DNA':
        if monoisotopic:
            weight_table = IUPACData.monoisotopic_unambiguous_dna_weights
        else:
            weight_table = IUPACData.unambiguous_dna_weights
    elif seq_type == 'RNA':
        if monoisotopic:
            weight_table = IUPACData.monoisotopic_unambiguous_rna_weights
        else:
            weight_table = IUPACData.unambiguous_rna_weights
    elif seq_type == 'protein':
        if monoisotopic:
            weight_table = IUPACData.monoisotopic_protein_weights
        else:
            weight_table = IUPACData.protein_weights
    else:
        raise ValueError("Allowed seq_types are DNA, RNA or protein, not %r"
                         % seq_type)

    if monoisotopic:
        water = 18.010565
    else:
        water = 18.0153

    try:
        weight = sum(weight_table[x] for x in seq) - (len(seq) - 1) * water
        if circular:
            weight -= water
    except KeyError as e:
        raise ValueError('%s is not a valid unambiguous letter for %s'
                         % (e, seq_type))
    except:
        raise

    if seq_type in ('DNA', 'RNA') and double_stranded:
        seq = str(Seq(seq).complement())
        weight += sum(weight_table[x] for x in seq) - (len(seq) - 1) * water
        if circular:
            weight -= water
    elif seq_type == 'protein' and double_stranded:
        raise ValueError('double-stranded proteins await their discovery')

    return weight


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

    """  # noqa for pep8 W291 trailing whitespace
    from Bio.Seq import reverse_complement, translate
    anti = reverse_complement(seq)
    comp = anti[::-1]
    length = len(seq)
    frames = {}
    for i in range(0, 3):
        fragment_length = 3 * ((length - i) // 3)
        frames[i + 1] = translate(seq[i:i + fragment_length], genetic_code)
        frames[-(i + 1)] = translate(anti[i:i + fragment_length], genetic_code)[::-1]

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
        subseq = seq[i:i + 60]
        csubseq = comp[i:i + 60]
        p = i // 3
        res += '%d/%d\n' % (i + 1, i / 3 + 1)
        res += '  ' + '  '.join(frames[3][p:p + 20]) + '\n'
        res += ' ' + '  '.join(frames[2][p:p + 20]) + '\n'
        res += '  '.join(frames[1][p:p + 20]) + '\n'
        # seq
        res += subseq.lower() + '%5d %%\n' % int(GC(subseq))
        res += csubseq.lower() + '\n'
        # - frames
        res += '  '.join(frames[-2][p:p + 20]) + ' \n'
        res += ' ' + '  '.join(frames[-1][p:p + 20]) + '\n'
        res += '  ' + '  '.join(frames[-3][p:p + 20]) + '\n\n'
    return res


shuffle_seed = randint(1, 9000)


def testseq(size=30, alphabet=IUPAC.unambiguous_dna, table=1, gc_target=None, persistent=True, from_start=True, to_stop=True, stop_symbol="*", truncate=True, messenger=False, shuffle=False, rand_seed=9001):
    """Generate and return a Seq object.

    This function will generate and return a custom Seq object
    using any IUPAC alphabet. These sequences are a faux representation
    of biological data and can be used for testing/demonstration purposes.

    Arguments:
        - size - The number of letters in the generated sequence.
        This preferably accepts an integer value and will attempt
        to convert any input to an integer.
        - alphabet - Any IUPAC alphabet can be used to generate
        the sequence. Defaults to unambiguous DNA.
        - table - Select which codon table you would like to use.  
        This preferably accepts an integer value and will attempt
        to convert any input to an integer. 
        - gc_target - The function will attempt to generate a sequence
        with a GC-content equal to the 'gc_target' argument. The argument will
        accept any integer between 0 and 100. Alternatively, if 'gc_target'
        is set to 'None', the 'gc_target' argument will be ignored when
        generating the sequence.
        - persistent - A boolean that, if set to True, will remove
        any stop codons that are generated by chance within the sequence.
        - from_start - A boolean that, if set to True, will ensure
        that the first codon in the generated sequence is a start codon.
        - to_stop - A boolean that, if set to True, will ensure
        that the last codon in the generated sequence is a stop codon.
        - stop_symbol - Single character string representing a translated
        stop codon.  This defaults to the asterisk, "*".
        - truncate - A boolean that, if set to True, will ensure that
        the size of any generated (non-protein) sequence is a multiple of three (3).
        - messenger - A boolean that, if set to True, will ensure that
        any RNA sequence generated will additionally have a 5'-UTR,
        3'-UTR, and a poly-A tail.
        - shuffle - A boolean that, if set to True, will ensure that every
        sequence generated is different from the last.
        - rand_seed - The seed used to generate the randomized sequence.
        This preferably accepts an integer value and will attempt
        to convert any input to an integer.

    Hey there! We can use the 'testseq' function to quickly generate sequences:

    >>> from Bio.SeqUtils import testseq
    >>> my_seq = testseq()
    >>> my_seq
    Seq('ATGGTGTTTGCACCAACAATACTATCGTAA', IUPACUnambiguousDNA())

    The default size of the generated sequence is 30 letters,
    but you can change that at any time, like so:

    >>> my_seq = testseq(1500)
    >>> len(my_seq)
    1500

    You can generate your sequence using any IUPAC alphabet, like so:

    >>> from Bio.Alphabet import IUPAC
    >>> my_seq = testseq(alphabet=IUPAC.extended_protein)
    >>> my_seq
    Seq('MLKDWKWYTCZXKYBVBUUNLWTDVRJUI*', HasStopCodon(ExtendedIUPACProtein(), '*'))

    You may have noticed that the above sequence starts with Methionine(M) and
    ends in an asterisk(*). That's because of the two arguments 'from_start'
    and 'to_stop' respectively. Curiously, there are no asterisks (or terminators)
    within the sequence either; this is due to the 'persistent' argument.
    The 'from_start', 'to_stop', and 'persistent' arguments are all set to True by default.
    You can read more about what they do in the "Arguments" section above. It's
    useful to note that all three of those arguments involve the use of codon tables!
    When generating your sequence, you can set which codon table you'd like to use:

    >>> my_seq = testseq(table=5)
    >>> my_seq
    Seq('ATGGTGTTTGCACCAACAATACTATCGTAA', IUPACUnambiguousDNA())

    Now we can translate our sequence with ease!

    >>> my_seq.translate(table=6)
    Seq('MVFAPTILSQ', IUPACProtein())

    Oops! We're missing a stop codon. We've generated a sequence using Table 5,
    but translated it using Table 6. Those tables don't share a common stop codon!
    Let's fix that...

    >>> my_seq.translate(table=5)
    Seq('MVFAPTMLS*', HasStopCodon(IUPACProtein(), '*'))

    That's better!

    The 'testseq' function can also attempt to generate sequences with a
    custom GC-content. You can alter your desired GC-content by declaring
    it in the 'gc_target' argument, like so:

    >>> my_seq = testseq(gc_target=60)
    >>> from Bio.SeqUtils import GC
    >>> GC(my_seq)
    46.666666666666664

    What happened? Note that in the above example, the sequence is at the
    default size of 30 letters. Since the sequence is generated letter by letter,
    larger sequences will have a tendency to be closer to the desired GC-content
    than smaller sequences. Let's try that again with a much larger sequence:

    >>> my_seq = testseq(size=10000, gc_target=60)
    >>> GC(my_seq)
    59.665966596659665

    Much better! It's also worth noting that the 'gc_target' argument is ignored
    when generating protein sequences.

    Lets revisit the 'size' argument for a moment. You will notice that when
    generating the sequence above, I declared a 'size' of 10000. Lets confirm
    whether that was generated as requested:

    >>> len(my_seq)
    9999

    That may seem like a bug, but it isn't! All non-protein sequences are
    truncated to a multiple of three (3). This is to allow for smooth translation
    from nucleotide to amino-acid alphabets. This behavior is controlled by the
    'truncate' argument, which is set to True by default. Let's see what happens
    when we set it to False:

    >>> my_seq = testseq(10000, truncate=False)
    >>> len(my_seq)
    10000

    The result above seems cleaner, but would result in a 'BiopythonWarning' if
    you translated that sequence. So please be careful when altering the
    'truncate' argument.

    Let us briefly discuss the 'messenger' argument. You can use it to add
    messenger RNA components to a generated RNA sequence. Though, it has
    two caveats. First, the messenger argument is ignored unless an RNA alphabet
    is declared. More importantly though, all mRNA compenents are added to
    the sequence addtionally. Let's look at an example:

    >>> my_seq = testseq(300, alphabet=IUPAC.unambiguous_rna, messenger=True)
    >>> len(my_seq)
    726

    Notice that the sequence requested was 300 letters, however the final length of
    the sequence is 726 letters. Those extra letters are the mRNA components. The
    generated sequence is buried in there, somewhere.

    Lastly, lets discuss the sequence generator itself. The sequence is created
    using a pseudo-random number generator which relies on a seed to process
    and spit out random numbers.

    Lets look at an example:

    >>> a = testseq()
    >>> b = testseq()
    >>> a == b
    True

    Since sequence "a" and sequence "b" were both generated using the same seed,
    they ended up being the exact same sequence. We can change that behavior
    to shuffle the seed every time the function is called by altering
    the 'shuffle' argument, like so:

    >>> a = testseq(shuffle=True)
    >>> b = testseq(shuffle=True)
    >>> a == b
    False

    If you'd like even more control over which sequence is generated, you
    can alter the 'rand_seed' arguement to manually change the seed:

    >>> a = testseq(rand_seed=12345)
    >>> b = testseq(rand_seed=67890)
    >>> a == b
    False
    
    Please note that if the 'shuffle' arguement is set to True,
    the 'rand_seed' arguement will be ignored.

    This concludes our discussion. Thanks again for using Biopython!
    Contribution by Adil Iqbal (2017).
    """
    # Set seed, gather data, clean-up logic, validate arguments.
    if shuffle:
        global shuffle_seed
        shuffle_seed += 1
        if shuffle_seed > 9000:
            shuffle_seed = 1
        seed(shuffle_seed)
    else:
        rand_seed = int(rand_seed)
        seed(rand_seed)
    typeof = _SeqType(alphabet)
    if not typeof.rna:
        messenger = False
    if typeof.rna and messenger:
        from_start = True
        to_stop = True
        persistent = True
    if persistent or from_start or to_stop:
        stop_symbol = str(stop_symbol)[0]
        codon_set = _CodonSet(alphabet, table, stop_symbol)
    if gc_target is not None and not typeof.protein:
        gc_target = int(gc_target)
        if gc_target > 100:
            gc_target = 100
        if gc_target < 0:
            gc_target = 0
        probability_table = _construct_probability_table(alphabet, gc_target)
    size = int(size)
    if not typeof.protein and truncate:
        size -= size % 3
    # Begin generating sequence.
    seq = ""
    for i in range(size):
        if gc_target is not None and not typeof.protein:
            seq += _pick_one(probability_table)
        else:
            roll = randint(0, len(alphabet.letters) - 1)
            seq += alphabet.letters[roll]
        if len(seq) >= 3 and len(seq) % 3 == 0 and not typeof.protein and persistent:
            # Replace stop codons with non-stop codons.
            this_codon = seq[-3:]
            if this_codon in codon_set.stop:
                roll = randint(0, len(codon_set.nonstop) - 1)
                new_codon = codon_set.nonstop[roll]
                seq = seq[:-3] + new_codon
    # Additional processing of generated sequence.
    x = 3
    if typeof.protein:
        x = 1
    if from_start:
        aug = None
        if typeof.dna:
            aug = "ATG"
        elif typeof.rna:
            aug = "AUG"
        elif typeof.protein:
            aug = "M"
        if aug is not None and aug in codon_set.start:
            start = aug
        else:
            roll = randint(0, len(codon_set.start) - 1)
            start = codon_set.start[roll]
        seq = start + seq[x:]
    if to_stop:
        roll = randint(0, len(codon_set.stop) - 1)
        stop = codon_set.stop[roll]
        seq = seq[:-x] + stop
    if messenger:
        seq = _add_messenger_parts(seq, size, alphabet, codon_set)
    if typeof.protein and stop_symbol in seq:
        alphabet = Alphabet.HasStopCodon(alphabet, stop_symbol)
    return Seq(seq, alphabet)


def _construct_probability_table(alphabet, gc_target):
    """Assign a probability of being chosen to each nucleotide based on desired GC-content. (PRIVATE)"""
    gc_nt_total = 2
    if alphabet == IUPAC.ambiguous_dna or alphabet == IUPAC.ambiguous_rna:
        gc_nt_total = 3
    probability_table = []
    total = len(alphabet.letters) - gc_nt_total
    for letter in alphabet.letters:
        if letter in ["G", "C", "S"]:
            value = gc_target / gc_nt_total
        else:
            value = (100 - gc_target) / total
        probability_table.append(_Letter(letter, value))
    sum_ = 0
    for letter in range(len(probability_table)):
        sum_ += probability_table[letter].probability_value
    for letter in range(len(probability_table)):
        probability_table[letter].probability_value /= sum_
    return probability_table


def _pick_one(probability_table):
    """Choose a nucleotide based on its probability of being chosen. (PRIVATE)"""
    roll = random()
    index = 0
    while roll > 0:
        roll -= probability_table[index].probability_value
        index += 1
    index -= 1
    return probability_table[index].letter


def _add_messenger_parts(seq, size, alphabet, codon_set):
    """Generate and add the 5' UTR, 3' UTR, and PolyA-Tail to an RNA sequence. (PRIVATE)"""
    utr_size = int(size / 3)
    utr5 = ""
    for i in range(utr_size):
        roll = randint(0, len(alphabet.letters) - 1)
        utr5 += alphabet.letters[roll]
        if len(utr5) >= 3 and len(utr5) % 3 == 0:
            # Replace start codons with non-start codons.
            this_codon = utr5[-3:]
            for startCodon in codon_set.start:
                if this_codon == startCodon:
                    roll = randint(0, len(codon_set.nonstart) - 1)
                    new_codon = codon_set.nonstart[roll]
                    utr5 = utr5[:-3] + new_codon
                    break
    utr3 = ""
    for i in range(utr_size):
        roll = randint(0, len(alphabet.letters) - 1)
        utr3 += alphabet.letters[roll]
    roll = randint(100, 251)
    poly_a_tail = "A" * roll
    seq = utr5 + seq + utr3 + poly_a_tail
    while len(seq) % 3 != 0:
        seq = seq[:-1]
    return seq


class _SeqType(object):
    """Evaluate alphabet to determine Seq type and return boolean values. (PRIVATE)"""

    def __init__(self, alphabet):
        self.dna = isinstance(alphabet, Alphabet.DNAAlphabet)
        self.rna = isinstance(alphabet, Alphabet.RNAAlphabet)
        self.protein = isinstance(alphabet, Alphabet.ProteinAlphabet)
        if self.dna is True and self.rna is False and self.protein is False:
            pass
        elif self.rna is True and self.dna is False and self.protein is False:
            pass
        elif self.protein is True and self.dna is False and self.rna is False:
            pass
        else:
            raise TypeError("Not a valid IUPAC alphabet.")


class _CodonSet(object):
    """Populate lists of codons from appropriate codon table. Return lists. (PRIVATE)"""

    def __init__(self, alphabet=IUPAC.unambiguous_dna, table=1, stop_symbol="*"):
        self.alphabet = alphabet
        self.table = table
        self.stop_symbol = stop_symbol
        self.typeof = _SeqType(self.alphabet)
        self.table = self._get_codon_table()
        self.start = self.table.start_codons
        self.stop = self.table.stop_codons
        self.nonstop = self._get_non_codons(self.stop)
        self.nonstart = self._get_non_codons(self.start)
        if self.typeof.protein:
            self._translate_codon_sets()
        del self.alphabet, self.table, self.stop_symbol, self.typeof

    def _get_codon_table(self):
        """Retrieve codon data from Bio.Data.CodonTable. (PRIVATE)"""
        if self.alphabet == IUPAC.unambiguous_dna or self.typeof.protein:
            codon_table = CodonTable.unambiguous_dna_by_id[self.table]
        elif self.alphabet == IUPAC.ambiguous_dna:
            codon_table = CodonTable.ambiguous_dna_by_id[self.table]
        elif self.alphabet == IUPAC.unambiguous_rna:
            codon_table = CodonTable.unambiguous_rna_by_id[self.table]
        elif self.alphabet == IUPAC.ambiguous_rna:
            codon_table = CodonTable.ambiguous_rna_by_id[self.table]
        else:
            codon_table = CodonTable.unambiguous_dna_by_id[self.table]
        return codon_table

    def _get_non_codons(self, exceptions):
        """Return a list of all codons in a table that are not exception codons. (PRIVATE)"""
        if self.typeof.rna:
            letters = IUPAC.unambiguous_rna.letters
        else:
            letters = IUPAC.unambiguous_dna.letters
        non_codons = []
        for nt1 in letters:
            for nt2 in letters:
                for nt3 in letters:
                    codon = nt1 + nt2 + nt3
                    for exception_codon in exceptions:
                        if codon == exception_codon:
                            break
                    else:
                        non_codons.append(codon)
        return non_codons

    def _translate_codon_sets(self):
        """Replace all codons in codon sets with corresponding amino acids. (PRIVATE)"""
        def translate_set(codon_set):
            amino_acids = []
            for codon in codon_set:
                if codon_set == self.stop:
                    translation = self.stop_symbol
                else:
                    translation = Seq(codon, IUPAC.unambiguous_dna).translate(table=self.table.id)._data
                this_amino_acid = translation
                amino_acids.append(this_amino_acid)
            return amino_acids

        self.start = translate_set(self.start)
        self.stop = translate_set(self.stop)
        self.nonstop = translate_set(self.nonstop)
        self.nonstart = translate_set(self.nonstart)


class _Letter(object):
    """Pair a letter with its probability of being chosen. (PRIVATE)"""

    def __init__(self, letter=None, value=None):
        self.letter = letter
        self.probability_value = value

    def __repr__(self):
        return "<" + str(self.letter) + " : " + str(self.probability_value) + ">"

# }}}


if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()
