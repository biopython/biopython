.. _`chapter:pairwise2`:

Pairwise alignments using pairwise2
===================================

**Please note that Bio.pairwise2 was deprecated in Release 1.80.** As an
alternative, please consider using ``Bio.Align.PairwiseAligner``
(described in Chapter :ref:`chapter:pairwise`).

``Bio.pairwise2`` contains essentially the same algorithms as ``water``
(local) and ``needle`` (global) from the
`EMBOSS <http://emboss.sourceforge.net/>`__ suite (see above) and should
return the same results. The ``pairwise2`` module has undergone some
optimization regarding speed and memory consumption recently (Biopython
versions >1.67) so that for short sequences (global alignments: ~2000
residues, local alignments ~600 residues) it’s faster (or equally fast)
to use ``pairwise2`` than calling EMBOSS’ ``water`` or ``needle`` via
the command line tools.

Suppose you want to do a global pairwise alignment between the same two
hemoglobin sequences from above (``HBA_HUMAN``, ``HBB_HUMAN``) stored in
``alpha.faa`` and ``beta.faa``:

.. doctest examples

.. code:: pycon

   >>> from Bio import pairwise2
   >>> from Bio import SeqIO
   >>> seq1 = SeqIO.read("alpha.faa", "fasta")
   >>> seq2 = SeqIO.read("beta.faa", "fasta")
   >>> alignments = pairwise2.align.globalxx(seq1.seq, seq2.seq)

As you see, we call the alignment function with ``align.globalxx``. The
tricky part are the last two letters of the function name (here:
``xx``), which are used for decoding the scores and penalties for
matches (and mismatches) and gaps. The first letter decodes the match
score, e.g. ``x`` means that a match counts 1 while mismatches have no
costs. With ``m`` general values for either matches or mismatches can be
defined (for full details see :py:mod:`Bio.pairwise2`). The
second letter decodes the cost for gaps; ``x`` means no gap costs at
all, with ``s`` different penalties for opening and extending a gap can
be assigned. So, ``globalxx`` means that only matches between both
sequences are counted.

Our variable ``alignments`` now contains a list of alignments (at least
one) which have the same optimal score for the given conditions. In our
example this are 80 different alignments with the score 72
(``Bio.pairwise2`` will return up to 1000 alignments). Have a look at
one of these alignments:

.. cont-doctest

.. code:: pycon

   >>> len(alignments)
   80
   >>> print(alignments[0])
   Alignment(seqA='MV-LSPADKTNV---K-A--A-WGKVGAHAG...YR-', seqB='MVHL-----T--PEEKSAVTALWGKV----...Y-H', score=72.0, start=0, end=217)

Each alignment is a named tuple consisting of the two aligned sequences,
the score, the start and the end positions of the alignment (in global
alignments the start is always 0 and the end the length of the
alignment). ``Bio.pairwise2`` has a function ``format_alignment`` for a
nicer printout:

.. cont-doctest

.. code:: pycon

   >>> print(pairwise2.format_alignment(*alignments[0]))
   MV-LSPADKTNV---K-A--A-WGKVGAHAG---EY-GA-EALE-RMFLSF----PTTK-TY--F...YR-
   || |     |     | |  | ||||        |  |  |||  |  |      |    |   |...|  
   MVHL-----T--PEEKSAVTALWGKV-----NVDE-VG-GEAL-GR--L--LVVYP---WT-QRF...Y-H
     Score=72
   <BLANKLINE>

Since Biopython 1.77 the required parameters can be supplied with
keywords. The last example can now also be written as:

.. cont-doctest

.. code:: pycon

   >>> alignments = pairwise2.align.globalxx(sequenceA=seq1.seq, sequenceB=seq2.seq)

Better alignments are usually obtained by penalizing gaps: higher costs
for opening a gap and lower costs for extending an existing gap. For
amino acid sequences match scores are usually encoded in matrices like
``PAM`` or ``BLOSUM``. Thus, a more meaningful alignment for our example
can be obtained by using the BLOSUM62 matrix, together with a gap open
penalty of 10 and a gap extension penalty of 0.5 (using ``globalds``):

.. code:: pycon

   >>> from Bio import pairwise2
   >>> from Bio import SeqIO
   >>> from Bio.Align import substitution_matrices
   >>> blosum62 = substitution_matrices.load("BLOSUM62")
   >>> seq1 = SeqIO.read("alpha.faa", "fasta")
   >>> seq2 = SeqIO.read("beta.faa", "fasta")
   >>> alignments = pairwise2.align.globalds(seq1.seq, seq2.seq, blosum62, -10, -0.5)
   >>> len(alignments)
   2
   >>> print(pairwise2.format_alignment(*alignments[0]))
   MV-LSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTY...KYR
   || |.|..|..|.|.|||| ......|............|.......||.
   MVHLTPEEKSAVTALWGKV-NVDEVGGEALGRLLVVYPWTQRFF...KYH
     Score=292.5

This alignment has the same score that we obtained earlier with EMBOSS
needle using the same sequences and the same parameters.

Local alignments are called similarly with the function
``align.localXX``, where again XX stands for a two letter code for the
match and gap functions:

.. doctest

.. code:: pycon

   >>> from Bio import pairwise2
   >>> from Bio.Align import substitution_matrices
   >>> blosum62 = substitution_matrices.load("BLOSUM62")
   >>> alignments = pairwise2.align.localds("LSPADKTNVKAA", "PEEKSAV", blosum62, -10, -1)
   >>> print(pairwise2.format_alignment(*alignments[0]))
   3 PADKTNV
     |..|..|
   1 PEEKSAV
     Score=16
   <BLANKLINE>

In recent Biopython versions, ``format_alignment`` will only print the
aligned part of a local alignment (together with the start positions in
1-based notation, as shown in the above example). If you are also
interested in the non- aligned parts of the sequences, use the
keyword-parameter ``full_sequences=True``:

.. doctest

.. code:: pycon

   >>> from Bio import pairwise2
   >>> from Bio.Align import substitution_matrices
   >>> blosum62 = substitution_matrices.load("BLOSUM62")
   >>> alignments = pairwise2.align.localds("LSPADKTNVKAA", "PEEKSAV", blosum62, -10, -1)
   >>> print(pairwise2.format_alignment(*alignments[0], full_sequences=True))
   LSPADKTNVKAA
     |..|..|   
   --PEEKSAV---
     Score=16
   <BLANKLINE>

Note that local alignments must, as defined by Smith & Waterman, have a
positive score (>0). Thus, ``pairwise2`` may return no alignments if no
score >0 has been obtained. Also, ``pairwise2`` will not report
alignments which are the result of the addition of zero-scoring
extensions on either site. In the next example, the pairs
serine/aspartic acid (S/D) and lysine/asparagine (K/N) both have a match
score of 0. As you see, the aligned part has not been extended:

.. doctest

.. code:: pycon

   >>> from Bio import pairwise2
   >>> from Bio.Align import substitution_matrices
   >>> blosum62 = substitution_matrices.load("BLOSUM62")
   >>> alignments = pairwise2.align.localds("LSSPADKTNVKKAA", "DDPEEKSAVNN", blosum62, -10, -1)
   >>> print(pairwise2.format_alignment(*alignments[0]))
   4 PADKTNV
     |..|..|
   3 PEEKSAV
     Score=16
   <BLANKLINE>

Instead of supplying a complete match/mismatch matrix, the match code
``m`` allows for easy defining general match/mismatch values. The next
example uses match/mismatch scores of 5/-4 and gap penalties
(open/extend) of 2/0.5 using ``localms``:

.. cont-doctest

.. code:: pycon

   >>> alignments = pairwise2.align.localms("AGAACT", "GAC", 5, -4, -2, -0.5)
   >>> print(pairwise2.format_alignment(*alignments[0]))
   2 GAAC
     | ||
   1 G-AC
     Score=13
   <BLANKLINE>

One useful keyword argument of the ``Bio.pairwise2.align`` functions is
``score_only``. When set to ``True`` it will only return the score of
the best alignment(s), but in a significantly shorter time. It will also
allow the alignment of longer sequences before a memory error is raised.
Another useful keyword argument is ``one_alignment_only=True`` which
will also result in some speed gain.

Unfortunately, ``Bio.pairwise2`` does not work with Biopython’s multiple
sequence alignment objects (yet). However, the module has some
interesting advanced features: you can define your own match and gap
functions (interested in testing affine logarithmic gap costs?), gap
penalties and end gaps penalties can be different for both sequences,
sequences can be supplied as lists (useful if you have residues that are
encoded by more than one character), etc. These features are hard (if at
all) to realize with other alignment tools. For more details see the
module's API documentation :py:mod:`Bio.pairwise2`.
