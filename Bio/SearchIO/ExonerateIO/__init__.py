# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.SearchIO support for Exonerate output formats.

This module adds support for handling Exonerate outputs. Exonerate is a generic
tool for pairwise sequence comparison that allows you to align sequences using
several different models.

Bio.SearchIO.ExonerateIO was tested on the following Exonerate versions and
models:

    - version: 2.2
    - models:
      - affine:local                - cdna2genome
      - coding2coding               - est2genome
      - genome2genome               - ner
      - protein2dna                 - protein2genome
      - ungapped                    - ungapped:translated

Although model testing were not exhaustive, ExonerateIO should be able to cope
with all Exonerate models. Please file a bug report if you stumble upon an
unparsable file.

More information on Exonerate is available on its home page at
www.ebi.ac.uk/~guy/exonerate/


Supported Formats
=================

  - Plain text alignment - 'exonerate-text'   - parsing, indexing
  - Vulgar line          - 'exonerate-vulgar' - parsing, indexing
  - Cigar line           - 'exonerate-cigar'  - parsing, indexing

On Exonerate, these output formats are not exclusive to one another. For
example, you may have both plain text and vulgar output in the same file.
ExonerateIO can only handle one of these at a time, however. If you have a file
containing both plain text and vulgar lines, for example, you have to pick
either 'exonerate-text' or 'exonerate-vulgar' to parse it.

Due to the cigar format specification, many features of the alignments such as
introns or frameshifts may be collapsed into a single feature (in this case,
they are labelled 'D' for 'deletion'). The parser does not attempt to guess
whether the D label it encounters is a real deletion or a collapsed feature.
As such, parsing or indexing using 'exonerate-cigar' may yield different results
compared to 'exonerate-text' or 'exonerate-vulgar'.


exonerate-text
==============

The plain text output / C4 alignment is the output triggered by the
'--showalignemnt' flag. Compared to the two other output formats, this format
contains the most information, having the complete query and hit sequences of
the alignment.

Here are some examples of the C4 output alignment that ExonerateIO can handle
(coordinates not written in scale)::

    1. simple ungapped alignments

           1 : ATGGGCAATATCCTTCGGAAAGGTCAGCAAAT :      56
               ||||||||||||||||||||||||||||||||
     1319275 : ATGGGCAATATCCTTCGGAAAGGTCAGCAAAT : 1319220

    2. alignments with frameshifts:

         129 : -TGCCGTTACCAT----GACGAAAGTATTAAT : 160
               -CysArgTyrHis----AspGluSerIleAsn
               #||||||||||||####|||||||||||||||
               #CysArgTyrHis####AspGluSerIleAsn
     1234593 : GTGCCGTTACCATCGGTGACGAAAGTATTAAT : 1234630

    3. alignments with introns and split codons:

        382 :    {A}                             {CC}AAA                 :    358
              AAA{T}  >>>> Target Intron 3 >>>>  {hr}LysATGAGCGATGAAAATA
              || { }++         55423 bp        ++{  } !  |||  ||||||||||
              AAC{L}gt.........................ag{eu}AspTTGAATGATGAAAATA
      42322 :    {C}                             {TG}GAT                 :  97769

    4. alignments with NER blocks

        111 : CAGAAAA--<   31  >--CTGCCCAGAAT--<   10  >--AACGAGCGTTCCG- :    184
              | |||||--< NER 1 >--| ||||| | |--< NER 2 >--|||  | ||||||-
     297911 : CTGAAAA--<   29  >--CCGCCCAAAGT--<   13  >--AACTGGAGTTCCG- : 297993

ExonerateIO utilizes the HSPFragment model quite extensively to deal with non-
ungapped alignments. For any single HSPFragment, if ExonerateIO sees an intron,
a NER block, or a frameshift, it will break the fragment into two HSPFragment
objects and adjust each of their start and end coordinate appropriately.

You may notice that Exonerate always uses the three letter amino acid codes to
display protein sequences. If the protein itself is part of the query sequence,
such as in the protein2dna model, ExonerateIO will transform the protein
sequence into using one letter codes. This is because the SeqRecord objects that
store the sequences are designed for single-letter sequences only. If Exonerate
also outputs the underlying nucleotide sequence, it will be saved into an
``aln_annotation`` entry as a list of triplets.

If the protein sequence is not part of the actual alignment, such as in the
est2genome or genome2genome models, ExonerateIO will keep the three letter codes
and store them as ``aln_annotation`` entries. In these cases, the hit and
query sequences may be used directly as SeqRecord objects as they are one-letter
nucleotide codes. The three-letter protein sequences are then stored as entries
in the ``aln_annotation`` dictionary.


For 'exonerate-text', ExonerateIO provides the following object attributes:

+-----------------+-------------------------+----------------------------------+
| Object          | Attribute               | Value                            |
+=================+=========================+==================================+
| QueryResult     | description             | query sequence description       |
|                 +-------------------------+----------------------------------+
|                 | id                      | query sequence ID                |
|                 +-------------------------+----------------------------------+
|                 | model                   | alignment model                  |
|                 +-------------------------+----------------------------------+
|                 | program                 | 'exonerate'                      |
+-----------------+-------------------------+----------------------------------+
| Hit             | description             | hit sequence description         |
|                 +-------------------------+----------------------------------+
|                 | id                      | hit sequence ID                  |
+-----------------+-------------------------+----------------------------------+
| HSP             | hit_split_codons        | list of split codon coordinates  |
|                 |                         | in the hit sequence              |
|                 +-------------------------+----------------------------------+
|                 | score                   | alignment score                  |
|                 +-------------------------+----------------------------------+
|                 | query_split_codons      | list of split codon coordinates  |
|                 |                         | in the query sequence            |
+-----------------+-------------------------+----------------------------------+
| HSPFragment     | aln_annotation          | alignment similarity string, hit |
|                 |                         | sequence annotation, and/or      |
|                 |                         | query sequence annotation        |
|                 +-------------------------+----------------------------------+
|                 | hit                     | hit sequence                     |
|                 +-------------------------+----------------------------------+
|                 | hit_end                 | hit sequence end coordinate      |
|                 +-------------------------+----------------------------------+
|                 | hit_frame               | hit sequence reading frame       |
|                 +-------------------------+----------------------------------+
|                 | hit_start               | hit sequence start coordinate    |
|                 +-------------------------+----------------------------------+
|                 | hit_strand              | hit sequence strand              |
|                 +-------------------------+----------------------------------+
|                 | query                   | query sequence                   |
|                 +-------------------------+----------------------------------+
|                 | query_end               | query sequence end coordinate    |
|                 +-------------------------+----------------------------------+
|                 | query_frame             | query sequence reading frame     |
|                 +-------------------------+----------------------------------+
|                 | query_start             | query sequence start coordinate  |
|                 +-------------------------+----------------------------------+
|                 | query_strand            | query sequence strand            |
+-----------------+-------------------------+----------------------------------+

Note that you can also use the default HSP or HSPFragment properties. For
example, to check the intron coordinates of your result you can use the
``query_inter_ranges`` or ``hit_inter_ranges`` properties:

    >>> from Bio import SearchIO
    >>> fname = 'Exonerate/exn_22_m_genome2genome.exn'
    >>> all_qresult = list(SearchIO.parse(fname, 'exonerate-text'))
    >>> hsp = all_qresult[-1][-1][-1]   # last qresult, last hit, last hsp
    >>> hsp
    HSP(...)
    >>> hsp.query_inter_ranges
    [(388, 449), (284, 319), (198, 198), (114, 161)]
    >>> hsp.hit_inter_ranges
    [(487387, 641682), (386207, 487327), (208677, 386123), (71917, 208639)]

Here you can see that for both query and hit introns, the coordinates
in each tuple is always (start, end) where start <= end. But when you compare
each tuple to the next, the coordinates decrease. This is an indication that
both the query and hit sequences lie on the minus strand. Exonerate outputs
minus strand results in a decreasing manner; the start coordinate is always
bigger than the end coordinate. ExonerateIO preserves the fragment ordering as a
whole, but uses its own standard to store an individual fragment's start and end
coordinates.

You may also notice that the third tuple in ``query_inter_ranges`` is (198, 198),
two exact same numbers. This means that the query sequence does not have any
gaps at that position. The gap is only present in the hit sequence, where we see
that the third tuple contains (208677, 386123), a gap of about 177k bases.

Another example is to use the ``hit_frame_all`` and ``query_frame_all`` to see if
there are any frameshifts in your alignment:

    >>> from Bio import SearchIO
    >>> fname = 'Exonerate/exn_22_m_coding2coding_fshifts.exn'
    >>> qresult = next(SearchIO.parse(fname, 'exonerate-text'))
    >>> hsp = qresult[0][0]      # first hit, first hsp
    >>> hsp
    HSP(...)
    >>> hsp.query_frame_all
    [1, 2, 2, 2]
    >>> hsp.hit_frame_all
    [1, 1, 3, 1]

Here you can see that the alignment as a whole has three frameshifts. The first
one occurs in the query sequence, after the first fragment (1 -> 2 shift), the
second one occurs in the hit sequence, after the second fragment (1 -> 3 shift),
and the last one also occurs in the hit sequence, before the last fragment (3 ->
1 shift).

There are other default HSP properties that you can use to ease your workflow.
Please refer to the HSP object documentation for more details.


exonerate-vulgar
================

The vulgar format provides a compact way of representing alignments created by
Exonerate. In general, it contains the same information as the plain text output
except for the 'model' information and the actual sequences themselves. You can
expect that the coordinates obtained from using 'exonerate-text' and
'exonerate-vulgar' to be the same. Both formats also creates HSPFragment using
the same triggers: introns, NER blocks, and/or frameshifts.


exonerate-cigar
===============

The cigar format provides an even more compact representation of Exonerate
alignments. However, this comes with a cost of losing information. In the cigar
format, for example, introns are treated as simple deletions. This makes it
impossible for the parser to distinguish between simple deletions or intron
regions. As such, 'exonerate-cigar' may produce different sets of coordinates
and fragments compared to 'exonerate-vulgar' or 'exonerate-text'.

"""

# Known issues & gotchas:
# - The cigar parser does not use the extended cigar string; only supports MID
# - Cigar and vulgar parsing results will most likely be different, due to the
#   different type of data stored by both formats

from .exonerate_text import ExonerateTextParser, ExonerateTextIndexer
from .exonerate_vulgar import ExonerateVulgarParser, ExonerateVulgarIndexer
from .exonerate_cigar import ExonerateCigarParser, ExonerateCigarIndexer


# if not used as a module, run the doctest
if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
