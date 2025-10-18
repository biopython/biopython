.. _`chapter:align`:

Sequence alignments
===================

Sequence alignments are a collection of two or more sequences that have
been aligned to each other – usually with the insertion of gaps, and the
addition of leading or trailing gaps – such that all the sequence
strings are the same length.

Alignments may extend over the full length of each sequence, or may be
limited to a subsection of each sequence. In Biopython, all sequence
alignments are represented by an ``Alignment`` object, described in
section :ref:`sec:alignmentobject`. ``Alignment`` objects can be
obtained by parsing the output of alignment software such as Clustal or
BLAT (described in section :ref:`sec:alignmentparsers`. or by using
Biopython’s pairwise sequence aligner, which can align two sequences to
each other (described in
Chapter :ref:`chapter:pairwise`).

See Chapter :ref:`chapter:msa` for a description of the
older ``MultipleSeqAlignment`` class and the parsers in ``Bio.AlignIO``
that parse the output of sequence alignment software, generating
``MultipleSeqAlignment`` objects.

.. _`sec:alignmentobject`:

Alignment objects
-----------------

The ``Alignment`` class is defined in ``Bio.Align``. Usually you would
get an ``Alignment`` object by parsing the output of alignment programs
(section :ref:`sec:alignmentparsers`) or by running Biopython’s
pairwise aligner (Chapter :ref:`chapter:pairwise`).
For the benefit of this section, however, we will create an
``Alignment`` object from scratch.

.. _`subsec:align_sequences_coordinates`:

Creating an Alignment object from sequences and coordinates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Suppose you have three sequences:

.. doctest ../Tests lib:numpy

.. code:: pycon

   >>> seqA = "CCGGTTTTT"
   >>> seqB = "AGTTTAA"
   >>> seqC = "AGGTTT"
   >>> sequences = [seqA, seqB, seqC]

To create an ``Alignment`` object, we also need the coordinates that
define how the sequences are aligned to each other. We use a NumPy array
for that:

.. cont-doctest

.. code:: pycon

   >>> import numpy as np
   >>> coordinates = np.array([[1, 3, 4, 7, 9], [0, 2, 2, 5, 5], [0, 2, 3, 6, 6]])

These coordinates define the alignment for the following sequence
segments:

-  ``SeqA[1:3]``, ``SeqB[0:2]``, and ``SeqC[0:2]`` are aligned to each
   other;

-  ``SeqA[3:4]`` and ``SeqC[2:3]`` are aligned to each other, with a gap
   of one nucleotide in ``seqB``;

-  ``SeqA[4:7]``, ``SeqB[2:5]``, and ``SeqC[3:6]`` are aligned to each
   other;

-  ``SeqA[7:9]`` is not aligned to ``seqB`` or ``seqC``.

Note that the alignment does not include the first nucleotide of
``seqA`` and last two nucleotides of ``seqB``.

Now we can create the ``Alignment`` object:

.. cont-doctest

.. code:: pycon

   >>> from Bio.Align import Alignment
   >>> alignment = Alignment(sequences, coordinates)
   >>> alignment  # doctest: +ELLIPSIS
   <Alignment object (3 rows x 8 columns) at ...>

The alignment object has an attribute ``sequences`` pointing to the
sequences included in this alignment:

.. cont-doctest

.. code:: pycon

   >>> alignment.sequences
   ['CCGGTTTTT', 'AGTTTAA', 'AGGTTT']

and an attribute ``coordinates`` with the alignment coordinates:

.. cont-doctest

.. code:: pycon

   >>> alignment.coordinates
   array([[1, 3, 4, 7, 9],
          [0, 2, 2, 5, 5],
          [0, 2, 3, 6, 6]])

Print the ``Alignment`` object to show the alignment explicitly:

.. cont-doctest

.. code:: pycon

   >>> print(alignment)
                     1 CGGTTTTT 9
                     0 AG-TTT-- 5
                     0 AGGTTT-- 6
   <BLANKLINE>

with the starting and end coordinate for each sequence are shown to the
left and right, respectively, of the alignment.

.. _`subsec:align_parse_printed_alignment`:

Creating an Alignment object from aligned sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you start out with the aligned sequences, with dashes representing
gaps, then you can calculate the coordinates using the
``parse_printed_alignment`` class method. This method is primarily employed in
Biopython’s alignment parsers (see
Section :ref:`sec:alignmentparsers`), but it may be useful for other
purposes. For example, you can construct the ``Alignment`` object from
aligned sequences as follows:

.. cont-doctest

.. code:: pycon

   >>> lines = ["CGGTTTTT", "AG-TTT--", "AGGTTT--"]
   >>> for line in lines:
   ...     print(line)
   ...
   CGGTTTTT
   AG-TTT--
   AGGTTT--
   >>> lines = [line.encode() for line in lines]  # convert to bytes
   >>> lines
   [b'CGGTTTTT', b'AG-TTT--', b'AGGTTT--']
   >>> sequences, coordinates = Alignment.parse_printed_alignment(lines)
   >>> sequences
   [b'CGGTTTTT', b'AGTTT', b'AGGTTT']
   >>> sequences = [sequence.decode() for sequence in sequences]
   >>> sequences
   ['CGGTTTTT', 'AGTTT', 'AGGTTT']
   >>> print(coordinates)
   [[0 2 3 6 8]
    [0 2 2 5 5]
    [0 2 3 6 6]]

The initial ``G`` nucleotide of ``seqA`` and the final ``CC``
nucleotides of ``seqB`` were not included in the alignment and is
therefore missing here. But this is easy to fix:

.. cont-doctest

.. code:: pycon

   >>> from Bio.Seq import Seq
   >>> sequences[0] = "C" + sequences[0]
   >>> sequences[1] = sequences[1] + "AA"
   >>> sequences
   ['CCGGTTTTT', 'AGTTTAA', 'AGGTTT']
   >>> coordinates[0, :] += 1
   >>> print(coordinates)
   [[1 3 4 7 9]
    [0 2 2 5 5]
    [0 2 3 6 6]]

Now we can create the ``Alignment`` object:

.. cont-doctest

.. code:: pycon

   >>> alignment = Alignment(sequences, coordinates)
   >>> print(alignment)
                     1 CGGTTTTT 9
                     0 AG-TTT-- 5
                     0 AGGTTT-- 6
   <BLANKLINE>

which identical to the ``Alignment`` object created above in
section :ref:`subsec:align_sequences_coordinates`.

By default, the ``coordinates`` argument to the ``Alignment``
initializer is ``None``, which assumes that there are no gaps in the
alignment. All sequences in an ungapped alignment must have the same
length. If the ``coordinates`` argument is ``None``, then the
initializer will fill in the ``coordinates`` attribute of the
``Alignment`` object for you:

.. cont-doctest

.. code:: pycon

   >>> ungapped_alignment = Alignment(["ACGTACGT", "AAGTACGT", "ACGTACCT"])
   >>> ungapped_alignment  # doctest: +ELLIPSIS
   <Alignment object (3 rows x 8 columns) at ...>
   >>> print(ungapped_alignment.coordinates)
   [[0 8]
    [0 8]
    [0 8]]
   >>> print(ungapped_alignment)
                     0 ACGTACGT 8
                     0 AAGTACGT 8
                     0 ACGTACCT 8
   <BLANKLINE>

.. _`subsec:align_common_attributes`:

Common alignment attributes
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following attributes are commonly found on ``Alignment`` objects:

-  ``sequences``: This is a list of the sequences aligned to each other.
   Depending on how the alignment was created, the sequences can have
   the following types:

   -  plain Python string;

   -  ``Seq``;

   -  ``MutableSeq``;

   -  ``SeqRecord``;

   -  ``bytes``;

   -  ``bytearray``;

   -  NumPy array with data type ``numpy.int32``;

   -  any other object with a contiguous buffer of format ``"c"``,
      ``"B"``, ``"i"``, or ``"I"``;

   -  lists or tuples of objects defined in the ``alphabet`` attribute
      of the ``PairwiseAligner`` object that created the alignment (see
      section :ref:`sec:generalized-pairwise`).

   For pairwise alignments (meaning an alignment of two sequences), the
   properties ``target`` and ``query`` are aliases for ``sequences[0]``
   and ``sequences[1]``, respectively.

-  ``coordinates``: A NumPy array of integers storing the sequence
   indices defining how the sequences are aligned to each other;

-  ``score``: The alignment score, as found by the parser in the
   alignment file, or as calculated by the ``PairwiseAligner`` (see
   section :ref:`sec:pairwise-basic`);

-  ``annotations``: A dictionary storing most other annotations
   associated with the alignment;

-  ``column_annotations``: A dictionary storing annotations that extend
   along the alignment and have the same length as the alignment, such
   as a consensus sequence (see
   section :ref:`subsec:align_clustal` for an example).

An ``Alignment`` object created by the parser in ``Bio.Align`` may have
additional attributes, depending on the alignment file format from which
the alignment was read.

.. _`subsec:slicing-indexing-alignment`:

Slicing and indexing an alignment
---------------------------------

Slices of the form ``alignment[k, i:j]``, where ``k`` is an integer and
``i`` and ``j`` are integers or are absent, return a string showing the
aligned sequence (including gaps) for the target (if ``k=0``) or the
query (if ``k=1``) that includes only the columns ``i`` through ``j`` in
the printed alignment.

To illustrate this, in the following example the printed alignment has 8
columns:

.. cont-doctest

.. code:: pycon

   >>> print(alignment)
                     1 CGGTTTTT 9
                     0 AG-TTT-- 5
                     0 AGGTTT-- 6
   <BLANKLINE>
   >>> alignment.length
   8

To get the aligned sequence strings individually, use

.. cont-doctest

.. code:: pycon

   >>> alignment[0]
   'CGGTTTTT'
   >>> alignment[1]
   'AG-TTT--'
   >>> alignment[2]
   'AGGTTT--'
   >>> alignment[0, :]
   'CGGTTTTT'
   >>> alignment[1, :]
   'AG-TTT--'
   >>> alignment[0, 1:-1]
   'GGTTTT'
   >>> alignment[1, 1:-1]
   'G-TTT-'

Columns to be included can also be selected using an iterable over
integers:

.. cont-doctest

.. code:: pycon

   >>> alignment[0, (1, 2, 4)]
   'GGT'
   >>> alignment[1, range(0, 5, 2)]
   'A-T'

To get the letter at position ``[i, j]`` of the printed alignment, use
``alignment[i, j]``; this will return ``"-"`` if a gap is found at that
position:

.. cont-doctest

.. code:: pycon

   >>> alignment[0, 2]
   'G'
   >>> alignment[2, 6]
   '-'

To get specific columns in the alignment, use

.. cont-doctest

.. code:: pycon

   >>> alignment[:, 0]
   'CAA'
   >>> alignment[:, 1]
   'GGG'
   >>> alignment[:, 2]
   'G-G'

Slices of the form ``alignment[i:j:k]`` return a new ``Alignment``
object including only sequences ``[i:j:k]`` of the alignment:

.. cont-doctest

.. code:: pycon

   >>> alignment[1:]
   <Alignment object (2 rows x 6 columns) at ...>
   >>> print(alignment[1:])
   target            0 AG-TTT 5
                     0 ||-||| 6
   query             0 AGGTTT 6
   <BLANKLINE>

Slices of the form ``alignment[:, i:j]``, where ``i`` and ``j`` are
integers or are absent, return a new ``Alignment`` object that includes
only the columns ``i`` through ``j`` in the printed alignment.

Extracting the first 4 columns for the example alignment above gives:

.. cont-doctest

.. code:: pycon

   >>> alignment[:, :4]
   <Alignment object (3 rows x 4 columns) at ...>
   >>> print(alignment[:, :4])
                     1 CGGT 5
                     0 AG-T 3
                     0 AGGT 4
   <BLANKLINE>

Similarly, extracting the last 6 columns gives:

.. cont-doctest

.. code:: pycon

   >>> alignment[:, -6:]
   <Alignment object (3 rows x 6 columns) at ...>
   >>> print(alignment[:, -6:])
                     3 GTTTTT 9
                     2 -TTT-- 5
                     2 GTTT-- 6
   <BLANKLINE>

The column index can also be an iterable of integers:

.. cont-doctest

.. code:: pycon

   >>> print(alignment[:, (1, 3, 0)])
                     0 GTC 3
                     0 GTA 3
                     0 GTA 3
   <BLANKLINE>

Calling ``alignment[:, :]`` returns a copy of the alignment.

Getting information about the alignment
---------------------------------------

Alignment shape
~~~~~~~~~~~~~~~

The number of aligned sequences is returned by ``len(alignment)``:

.. cont-doctest

.. code:: pycon

   >>> len(alignment)
   3

The alignment length is defined as the number of columns in the
alignment as printed. This is equal to the sum of the number of matches,
number of mismatches, and the total length of gaps in each sequence:

.. cont-doctest

.. code:: pycon

   >>> alignment.length
   8

The ``shape`` property returns a tuple consisting of the length of the
alignment and the number of columns in the alignment as printed:

.. cont-doctest

.. code:: pycon

   >>> alignment.shape
   (3, 8)

Comparing alignments
~~~~~~~~~~~~~~~~~~~~

Two alignments are equal to each other (meaning that
``alignment1 == alignment2`` evaluates to ``True``) if each of the
sequences in ``alignment1.sequences`` and ``alignment2.sequences`` are
equal to each other, and ``alignment1.coordinates`` and
``alignment2.coordinates`` contain the same coordinates. If either of
these conditions is not fulfilled, then ``alignment1 == alignment2``
evaluates to ``False``. Inequality of two alignments (e.g.,
``alignment1 < alignment2``) is established by first comparing
``alignment1.sequences`` and ``alignment2.sequences``, and if they are
equal, by comparing ``alignment1.coordinates`` to
``alignment2.coordinates``.

Finding the indices of aligned sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For pairwise alignments, the ``aligned`` property of an alignment
returns the start and end indices of subsequences in the target and
query sequence that were aligned to each other. If the alignment between
target (t) and query (q) consists of :math:`N` chunks, you get two
tuples of length :math:`N`:

.. code:: text

   (((t_start1, t_end1), (t_start2, t_end2), ..., (t_startN, t_endN)),
    ((q_start1, q_end1), (q_start2, q_end2), ..., (q_startN, q_endN)))

For example,

.. cont-doctest

.. code:: pycon

   >>> pairwise_alignment = alignment[:2, :]
   >>> print(pairwise_alignment)
   target            1 CGGTTTTT 9
                     0 .|-|||-- 8
   query             0 AG-TTT-- 5
   <BLANKLINE>
   >>> print(pairwise_alignment.aligned)
   [[[1 3]
     [4 7]]
   <BLANKLINE>
    [[0 2]
     [2 5]]]

Note that different alignments may have the same subsequences aligned to
each other. In particular, this may occur if alignments differ from each
other in terms of their gap placement only:

.. cont-doctest

.. code:: pycon

   >>> pairwise_alignment1 = Alignment(["AAACAAA", "AAAGAAA"],
   ...                                 np.array([[0, 3, 4, 4, 7], [0, 3, 3, 4, 7]]))  # fmt: skip
   ...
   >>> pairwise_alignment2 = Alignment(["AAACAAA", "AAAGAAA"],
   ...                                 np.array([[0, 3, 3, 4, 7], [0, 3, 4, 4, 7]]))  # fmt: skip
   ...
   >>> print(pairwise_alignment1)
   target            0 AAAC-AAA 7
                     0 |||--||| 8
   query             0 AAA-GAAA 7
   <BLANKLINE>
   >>> print(pairwise_alignment2)
   target            0 AAA-CAAA 7
                     0 |||--||| 8
   query             0 AAAG-AAA 7
   <BLANKLINE>
   >>> pairwise_alignment1.aligned
   array([[[0, 3],
           [4, 7]],
   <BLANKLINE>
          [[0, 3],
           [4, 7]]])
   >>> pairwise_alignment2.aligned
   array([[[0, 3],
           [4, 7]],
   <BLANKLINE>
          [[0, 3],
           [4, 7]]])

The property ``indices`` returns a 2D NumPy array with the sequence
index of each letter in the alignment, with gaps indicated by -1:

.. cont-doctest

.. code:: pycon

   >>> print(alignment)
                     1 CGGTTTTT 9
                     0 AG-TTT-- 5
                     0 AGGTTT-- 6
   <BLANKLINE>
   >>> alignment.indices
   array([[ 1,  2,  3,  4,  5,  6,  7,  8],
          [ 0,  1, -1,  2,  3,  4, -1, -1],
          [ 0,  1,  2,  3,  4,  5, -1, -1]])

The property ``inverse_indices`` returns a list of 1D NumPy arrays, one
for each of the aligned sequences, with the column index in the
alignment for each letter in the sequence. Letters not included in the
alignment are indicated by -1:

.. cont-doctest

.. code:: pycon

   >>> alignment.sequences
   ['CCGGTTTTT', 'AGTTTAA', 'AGGTTT']
   >>> alignment.inverse_indices  # doctest: +NORMALIZE_WHITESPACE
   [array([-1,  0,  1,  2,  3,  4,  5,  6,  7]),
    array([ 0,  1,  3,  4,  5, -1, -1]),
    array([0, 1, 2, 3, 4, 5])]

.. _`paragraph:alignment_counts`:

Counting identities, mismatches, and gaps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``counts`` method counts the number of identities, mismatches, aligned
letters, and gaps (insertions and deletions) in an alignment.  The return
value is an ``AlignmentCounts`` object:

.. cont-doctest

.. code:: pycon

   >>> print(pairwise_alignment)
   target            1 CGGTTTTT 9
                     0 .|-|||-- 8
   query             0 AG-TTT-- 5
   <BLANKLINE>
   >>> counts = pairwise_alignment.counts()
   >>> counts  # doctest: +ELLIPSIS
   <AlignmentCounts object (5 aligned letters; 4 identities; 1 mismatches; 3 gaps) at ...>

You can print the alignment to get detailed information about the alignment
statistics:

.. cont-doctest

.. code:: pycon

   >>> print(counts)
   AlignmentCounts object with
       aligned = 5:
           identities = 4,
           mismatches = 1.
       gaps = 3:
           left_gaps = 0:
               left_insertions = 0:
                   open_left_insertions = 0,
                   extend_left_insertions = 0;
               left_deletions = 0:
                   open_left_deletions = 0,
                   extend_left_deletions = 0;
           internal_gaps = 1:
               internal_insertions = 0:
                   open_internal_insertions = 0,
                   extend_internal_insertions = 0;
               internal_deletions = 1:
                   open_internal_deletions = 1,
                   extend_internal_deletions = 0;
           right_gaps = 2:
               right_insertions = 0:
                   open_right_insertions = 0,
                   extend_right_insertions = 0;
               right_deletions = 2:
                   open_right_deletions = 1,
                   extend_right_deletions = 1.
   <BLANKLINE>

These statistics can also be obtained as properties of the ``counts`` object:

.. cont-doctest

.. code:: pycon

   >>> counts.aligned
   5
   >>> counts.identities
   4
   >>> counts.mismatches
   1
   >>> counts.gaps
   3
   >>> counts.insertions
   0
   >>> counts.deletions
   3
   >>> counts.internal_deletions
   1
   >>> counts.right_deletions
   2
   >>> counts.extend_right_deletions
   1

Use a single character as the argument to ``.counts`` to specify a wildcard
character, which is ignored when counting identities, positives, and mismatches
(e.g. ``"?"`` or ``"N"`` are commonly used as wildcard characters):

.. cont-doctest

.. code:: pycon

   >>> pairwise_alignment.sequences[0] = "CCGGT?TTT"
   >>> print(pairwise_alignment)
   target            1 CGGT?TTT 9
                     0 .|-|.|-- 8
   query             0 AG-TTT-- 5
   <BLANKLINE>
   >>> counts = pairwise_alignment.counts()
   >>> counts.identities
   3
   >>> counts.mismatches
   2
   >>> counts = pairwise_alignment.counts("?")
   >>> counts.identities
   3
   >>> counts.mismatches
   1

Here, the alignment between ``?`` and ``T`` is not counted as a mismatch, as
``?`` is the wildcard character.

Use a substitution matrix (see section :ref:`sec:substitution_matrices`) as the
argument to also calculate the number of positive matches, and the total
substitution score:

.. cont-doctest

.. code:: pycon

   >>> protein_alignment = Alignment(
   ...     ["EPQSDPSVEPPLSQETFSDLWKLLPE", "EPSSETGMDPPLSQETFEDLWSLLPD"]
   ... )
   >>> print(protein_alignment)
   target            0 EPQSDPSVEPPLSQETFSDLWKLLPE 26
                     0 ||.|.....||||||||.|||.|||. 26
   query             0 EPSSETGMDPPLSQETFEDLWSLLPD 26
   <BLANKLINE>
   >>> counts = protein_alignment.counts()
   >>> print(counts.identities)
   17
   >>> print(counts.mismatches)
   9
   >>> print(counts.positives)
   None
   >>> print(counts.substitution_score)
   None

.. cont-doctest

.. code:: pycon

   >>> from Bio.Align import substitution_matrices

.. cont-doctest

.. code:: pycon

   >>> blosum62 = substitution_matrices.load("BLOSUM62")
   >>> counts = protein_alignment.counts(blosum62)
   >>> print(counts.identities)
   17
   >>> print(counts.mismatches)
   9
   >>> print(counts.positives)
   21
   >>> print(counts.substitution_score)
   101.0

.. cont-doctest

.. code:: pycon

   >>> blosum45 = substitution_matrices.load("BLOSUM45")
   >>> counts = protein_alignment.counts(blosum45)
   >>> print(counts.identities)
   17
   >>> print(counts.mismatches)
   9
   >>> print(counts.positives)
   21
   >>> print(counts.substitution_score)
   122.0

.. cont-doctest

.. code:: pycon

   >>> blosum90 = substitution_matrices.load("BLOSUM90")
   >>> counts = protein_alignment.counts(blosum90)
   >>> print(counts.identities)
   17
   >>> print(counts.mismatches)
   9
   >>> print(counts.positives)
   20
   >>> print(counts.substitution_score)
   109.0

Use a pairwise aligner object (see Chapter :ref:`chapter:pairwise`)  to
also calculate the gap score and the alignment score:

.. cont-doctest

.. code:: pycon

   >>> from Bio.Align import PairwiseAligner

.. cont-doctest

.. code:: pycon

   >>> pairwise_alignment.sequences[0] = "CCGGTTTTT"
   >>> print(pairwise_alignment)
   target            1 CGGTTTTT 9
                     0 .|-|||-- 8
   query             0 AG-TTT-- 5
   <BLANKLINE>
   >>> aligner = PairwiseAligner("blastn")
   >>> counts = pairwise_alignment.counts(aligner)
   >>> counts  # doctest: +ELLIPSIS
   <AlignmentCounts object (score = -11.0; substitution score = 5.0; gap score = -16.0; 5 aligned letters; 4 identities; 1 mismatches; 4 positives; 3 gaps) at ...>
   >>> print(counts)
   AlignmentCounts object with
       score = -11.0:
           substitution_score = 5.0,
           gap_score = -16.0.
       aligned = 5:
           identities = 4,
           positives = 4,
           mismatches = 1.
       gaps = 3:
           left_gaps = 0:
               left_insertions = 0:
                   open_left_insertions = 0,
                   extend_left_insertions = 0;
               left_deletions = 0:
                   open_left_deletions = 0,
                   extend_left_deletions = 0;
           internal_gaps = 1:
               internal_insertions = 0:
                   open_internal_insertions = 0,
                   extend_internal_insertions = 0;
               internal_deletions = 1:
                   open_internal_deletions = 1,
                   extend_internal_deletions = 0;
           right_gaps = 2:
               right_insertions = 0:
                   open_right_insertions = 0,
                   extend_right_insertions = 0;
               right_deletions = 2:
                   open_right_deletions = 1,
                   extend_right_deletions = 1.
   <BLANKLINE>

.. cont-doctest

.. code:: pycon

   >>> aligner = PairwiseAligner("blastp")
   >>> counts = protein_alignment.counts(aligner)
   >>> counts  # doctest: +ELLIPSIS
   <AlignmentCounts object (score = 101.0; substitution score = 101.0; gap score = 0.0; 26 aligned letters; 17 identities; 9 mismatches; 21 positives; 0 gaps) at ...>
   >>> print(counts)
   AlignmentCounts object with
       score = 101.0:
           substitution_score = 101.0,
           gap_score = 0.0.
       aligned = 26:
           identities = 17,
           positives = 21,
           mismatches = 9.
       gaps = 0:
           left_gaps = 0:
               left_insertions = 0:
                   open_left_insertions = 0,
                   extend_left_insertions = 0;
               left_deletions = 0:
                   open_left_deletions = 0,
                   extend_left_deletions = 0;
           internal_gaps = 0:
               internal_insertions = 0:
                   open_internal_insertions = 0,
                   extend_internal_insertions = 0;
               internal_deletions = 0:
                   open_internal_deletions = 0,
                   extend_internal_deletions = 0;
           right_gaps = 0:
               right_insertions = 0:
                   open_right_insertions = 0,
                   extend_right_insertions = 0;
               right_deletions = 0:
                   open_right_deletions = 0,
                   extend_right_deletions = 0.
   <BLANKLINE>

Note that the pairwise aligner in the argument to ``counts``  does not need to
be the same as aligner that was used the create the alignment.

For an alignment of more than two sequences, the number of identities,
mismatches, and gaps are calculated and summed for all pairs of sequences in
the alignment.

.. cont-doctest

.. code:: pycon

   >>> print(alignment)
                     1 CGGTTTTT 9
                     0 AG-TTT-- 5
                     0 AGGTTT-- 6
   <BLANKLINE>
   >>> counts = alignment.counts()
   >>> counts.aligned
   16
   >>> counts.identities
   14
   >>> counts.mismatches
   2
   >>> counts.insertions
   1
   >>> counts.deletions
   5

Here, insertions are defined as sequence insertions of a later sequence into an
earlier sequence in the alignment. In contrast to the pairwise alignment above,
the distinction between insertions and deletions may not be meaningful for a
multiple sequence alignment, and you will probably be more interested in the
number of gaps (= insertions + deletions):

.. cont-doctest

.. code:: pycon

   >>> counts.gaps
   6
   >>> counts.left_gaps
   0
   >>> counts.right_gaps
   4
   >>> counts.internal_gaps
   2

An ``AlignmentCounts`` object has the following properties:

================================ =============================================================================================================
**property**                     **description**
================================ =============================================================================================================
``score``                        Alignment score, or ``None`` if unknown
``substitution_score``           Total substitution score of letters aligned to each other, or ``None`` if unknown
``gap_score``                    Total gap score, or ``None`` if unknown
``aligned``                      The number of aligned letters in the alignment. This quantity is also calculated if some or all of the sequences are undefined. If all sequences are known, then ``aligned`` = ``identities`` + ``mismatches``. If some sequences are undefined, then ``aligned`` > ``identities`` + ``mismatches``.
``identities``                   Number of matched letters in the alignment
``mismatches``                   Number of mismatched letters in the alignment
``positives``                    Number of aligned letters with a positive substitution score
``gaps``                         Total gap length
``open_gaps``                    Number of gaps opened in the alignment
``extend_gaps``                  Number of gap extensions in the alignment
``insertions``                   Total number of letters inserted
``open_insertions``              Number of insertion gaps opened in the alignment
``extend_insertions``            Number of insertion gap extensions in the alignment
``deletions``                    Total number of letters deleted
``open_deletions``               Number of deletion gaps opened in the alignment
``extend_deletions``             Number of deletion gap extensions in the alignment
``left_gaps``                    Total gap length on the left side of the alignment
``open_left_gaps``               Number of gaps opened on the left side of the alignment
``extend_left_gaps``             Number of gap extensions on the left side of the alignment
``left_insertions``              Number of letters inserted on the left side of the alignment
``open_left_insertions``         Number of insertion gaps opened on the left side of the alignment
``extend_left_insertions``       Number of insertion gap extensions on the left side of the alignment
``left_deletions``               Number of characters deleted on the left side of the alignment
``open_left_deletions``          The number of deletion gaps opened on the left side of the alignment
``extend_left_deletions``        Number of deletion gap extensions on the left side of the alignment
``internal_gaps``                Total length of gaps within the alignment
``open_internal_gaps``           Number of gaps opened in the interior of the alignment
``extend_internal_gaps``         Number of gap extensions in the interior of the alignment
``internal_insertions``          Number of letters inserted in the interior of the alignment
``open_internal_insertions``     Number of insertion gaps opened in the interior of the alignment
``extend_internal_insertions``   Number of insertion gas extensions in the interior of the alignment
``internal_deletions``           Number of characters deleted from the alignment
``open_internal_deletions``      Number of deletion gaps opened in the interior of the alignment
``extend_internal_deletions``    Number of deletion gap extensions in the interior of the alignment
``right_gaps``                   Total gap length on the right side of the alignment
``open_right_gaps``              Number of gaps opened on the right side of the alignment
``extend_right_gaps``            Number of gap extensions on the right side of the alignment
``right_insertions``             Number of letters inserted on the right side of the alignment
``open_right_insertions``        The number of insertion gaps opened on the right side of the alignment
``extend_right_insertions``      The number of insertion gap extensions on the right side of the alignment
``right_deletions``              Number of letters deleted on the right side of the alignment
``open_right_deletions``         Number of deletion gaps opened on the right side of the alignment
``extend_right_deletions``       Number of deletion gap extensions on the right side of the alignment
================================ =============================================================================================================


Letter frequencies
~~~~~~~~~~~~~~~~~~

The ``frequencies`` method calculates how often each letter appears in
each column of the alignment:

.. cont-doctest

.. code:: pycon

   >>> alignment.frequencies  # doctest: +NORMALIZE_WHITESPACE
   {'C': array([1., 0., 0., 0., 0., 0., 0., 0.]),
    'G': array([0., 3., 2., 0., 0., 0., 0., 0.]),
    'T': array([0., 0., 0., 3., 3., 3., 1., 1.]),
    'A': array([2., 0., 0., 0., 0., 0., 0., 0.]),
    '-': array([0., 0., 1., 0., 0., 0., 2., 2.])}

Substitutions
~~~~~~~~~~~~~

Use the ``substitutions`` method to find the number of substitutions
between each pair of nucleotides:

.. cont-doctest

.. code:: pycon

   >>> m = alignment.substitutions
   >>> print(m)
       A   C   G   T
   A 1.0 0.0 0.0 0.0
   C 2.0 0.0 0.0 0.0
   G 0.0 0.0 4.0 0.0
   T 0.0 0.0 0.0 9.0
   <BLANKLINE>

Note that the matrix is not symmetric: The counts for a row letter R and
a column letter C is the number of times letter R in a sequence is
replaced by letter C in a sequence appearing below it. For example, the
number of ``C``\ ’s that are aligned to an ``A`` in a later sequence is

.. cont-doctest

.. code:: pycon

   >>> m["C", "A"]
   2.0

while the number of A’s that are aligned to a C in a later sequence is

.. cont-doctest

.. code:: pycon

   >>> m["A", "C"]
   0.0

To get a symmetric matrix, use

.. cont-doctest

.. code:: pycon

   >>> m += m.transpose()
   >>> m /= 2.0
   >>> print(m)
       A   C   G   T
   A 1.0 1.0 0.0 0.0
   C 1.0 0.0 0.0 0.0
   G 0.0 0.0 4.0 0.0
   T 0.0 0.0 0.0 9.0
   <BLANKLINE>
   >>> m["A", "C"]
   1.0
   >>> m["C", "A"]
   1.0

The total number of substitutions between ``A``\ ’s and ``T``\ ’s in the
alignment is 1.0 + 1.0 = 2.

Alignments as arrays
~~~~~~~~~~~~~~~~~~~~

Using NumPy, you can turn the ``alignment`` object into an array of
letters. In particular, this may be useful for fast calculations on the
alignment content.

.. cont-doctest

.. code:: pycon

   >>> align_array = np.array(alignment)
   >>> align_array.shape
   (3, 8)
   >>> align_array  # doctest: +NORMALIZE_WHITESPACE
   array([[b'C', b'G', b'G', b'T', b'T', b'T', b'T', b'T'],
          [b'A', b'G', b'-', b'T', b'T', b'T', b'-', b'-'],
          [b'A', b'G', b'G', b'T', b'T', b'T', b'-', b'-']], dtype='|S1')

By default, this will give you an array of ``bytes`` characters (with
data type ``dtype='|S1'``). You can create an array of Unicode (Python
string) characters by using ``dtype='U'``:

.. cont-doctest

.. code:: pycon

   >>> align_array = np.array(alignment, dtype="U")

.. code:: pycon

   >>> align_array  # doctest: +NORMALIZE_WHITESPACE
   array([['C', 'G', 'G', 'T', 'T', 'T', 'T', 'T'],
          ['A', 'G', '-', 'T', 'T', 'T', '-', '-'],
          ['A', 'G', 'G', 'T', 'T', 'T', '-', '-']], dtype='<U1')

(the printed ``dtype`` will be '<U1' or '>U1' depending on whether your system
is little-endian or big-endian, respectively).
Note that the ``alignment`` object and the NumPy array ``align_array``
are separate objects in memory - editing one will not update the other!

Operations on an alignment
--------------------------

Sorting an alignment
~~~~~~~~~~~~~~~~~~~~

The ``sort`` method sorts the alignment sequences. By default, sorting
is done based on the ``id`` attribute of each sequence if available, or
the sequence contents otherwise.

.. cont-doctest

.. code:: pycon

   >>> print(alignment)
                     1 CGGTTTTT 9
                     0 AG-TTT-- 5
                     0 AGGTTT-- 6
   <BLANKLINE>
   >>> alignment.sort()
   >>> print(alignment)
                     0 AGGTTT-- 6
                     0 AG-TTT-- 5
                     1 CGGTTTTT 9
   <BLANKLINE>

Alternatively, you can supply a ``key`` function to determine the sort
order. For example, you can sort the sequences by increasing GC content:

.. cont-doctest

.. code:: pycon

   >>> from Bio.SeqUtils import gc_fraction
   >>> alignment.sort(key=gc_fraction)
   >>> print(alignment)
                     0 AG-TTT-- 5
                     0 AGGTTT-- 6
                     1 CGGTTTTT 9
   <BLANKLINE>

Note that the ``key`` function is applied to the full sequence
(including the initial ``A`` and final ``GG`` nucleotides of ``seqB``),
not just to the aligned part.

The ``reverse`` argument lets you reverse the sort order to obtain the
sequences in decreasing GC content:

.. cont-doctest

.. code:: pycon

   >>> alignment.sort(key=gc_fraction, reverse=True)
   >>> print(alignment)
                     1 CGGTTTTT 9
                     0 AGGTTT-- 6
                     0 AG-TTT-- 5
   <BLANKLINE>

Reverse-complementing the alignment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Reverse-complementing an alignment will take the reverse complement of
each sequence, and recalculate the coordinates:

.. cont-doctest

.. code:: pycon

   >>> alignment.sequences
   ['CCGGTTTTT', 'AGGTTT', 'AGTTTAA']
   >>> rc_alignment = alignment.reverse_complement()
   >>> print(rc_alignment.sequences)
   ['AAAAACCGG', 'AAACCT', 'TTAAACT']
   >>> print(rc_alignment)
                     0 AAAAACCG 8
                     0 --AAACCT 6
                     2 --AAA-CT 7
   <BLANKLINE>
   >>> alignment[:, :4].sequences
   ['CCGGTTTTT', 'AGGTTT', 'AGTTTAA']
   >>> print(alignment[:, :4])
                     1 CGGT 5
                     0 AGGT 4
                     0 AG-T 3
   <BLANKLINE>
   >>> rc_alignment = alignment[:, :4].reverse_complement()
   >>> rc_alignment[:, :4].sequences
   ['AAAAACCGG', 'AAACCT', 'TTAAACT']
   >>> print(rc_alignment[:, :4])
                     4 ACCG 8
                     2 ACCT 6
                     4 A-CT 7
   <BLANKLINE>

Reverse-complementing an alignment preserves its column annotations (in
reverse order), but discards all other annotations.

Adding alignments
~~~~~~~~~~~~~~~~~

Alignments can be added together to form an extended alignment if they
have the same number of rows. As an example, let’s first create two
alignments:

.. cont-doctest

.. code:: pycon

   >>> from Bio.Seq import Seq
   >>> from Bio.SeqRecord import SeqRecord
   >>> a1 = SeqRecord(Seq("AAAAC"), id="Alpha")
   >>> b1 = SeqRecord(Seq("AAAC"), id="Beta")
   >>> c1 = SeqRecord(Seq("AAAAG"), id="Gamma")
   >>> a2 = SeqRecord(Seq("GTT"), id="Alpha")
   >>> b2 = SeqRecord(Seq("TT"), id="Beta")
   >>> c2 = SeqRecord(Seq("GT"), id="Gamma")
   >>> left = Alignment(
   ...     [a1, b1, c1], coordinates=np.array([[0, 3, 4, 5], [0, 3, 3, 4], [0, 3, 4, 5]])
   ... )
   >>> left.annotations = {"tool": "demo", "name": "start"}
   >>> left.column_annotations = {"stats": "CCCXC"}
   >>> right = Alignment(
   ...     [a2, b2, c2], coordinates=np.array([[0, 1, 2, 3], [0, 0, 1, 2], [0, 1, 1, 2]])
   ... )
   >>> right.annotations = {"tool": "demo", "name": "end"}
   >>> right.column_annotations = {"stats": "CXC"}

Now, let’s look at these two alignments:

.. cont-doctest

.. code:: pycon

   >>> print(left)
   Alpha             0 AAAAC 5
   Beta              0 AAA-C 4
   Gamma             0 AAAAG 5
   <BLANKLINE>
   >>> print(right)
   Alpha             0 GTT 3
   Beta              0 -TT 2
   Gamma             0 G-T 2
   <BLANKLINE>

Adding the two alignments will combine the two alignments row-wise:

.. cont-doctest

.. code:: pycon

   >>> combined = left + right
   >>> print(combined)
   Alpha             0 AAAACGTT 8
   Beta              0 AAA-C-TT 6
   Gamma             0 AAAAGG-T 7
   <BLANKLINE>

For this to work, both alignments must have the same number of sequences
(here they both have 3 rows):

.. cont-doctest

.. code:: pycon

   >>> len(left)
   3
   >>> len(right)
   3
   >>> len(combined)
   3

The sequences are ``SeqRecord`` objects, which can be added together.
Refer to Chapter :ref:`chapter:seq_annot` for
details of how the annotation is handled. This example is a special case
in that both original alignments shared the same names, meaning when the
rows are added they also get the same name.

Any common annotations are preserved, but differing annotation is lost.
This is the same behavior used in the ``SeqRecord`` annotations and is
designed to prevent accidental propagation of inappropriate values:

.. cont-doctest

.. code:: pycon

   >>> combined.annotations
   {'tool': 'demo'}

Similarly any common per-column-annotations are combined:

.. cont-doctest

.. code:: pycon

   >>> combined.column_annotations
   {'stats': 'CCCXCCXC'}

Mapping a pairwise sequence alignment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Suppose you have a pairwise alignment of a transcript to a chromosome:

.. cont-doctest

.. code:: pycon

   >>> chromosome = "AAAAAAAACCCCCCCAAAAAAAAAAAGGGGGGAAAAAAAA"
   >>> transcript = "CCCCCCCGGGGGG"
   >>> sequences1 = [chromosome, transcript]
   >>> coordinates1 = np.array([[8, 15, 26, 32], [0, 7, 7, 13]])
   >>> alignment1 = Alignment(sequences1, coordinates1)
   >>> print(alignment1)
   target            8 CCCCCCCAAAAAAAAAAAGGGGGG 32
                     0 |||||||-----------|||||| 24
   query             0 CCCCCCC-----------GGGGGG 13
   <BLANKLINE>

and a pairwise alignment between the transcript and a sequence (e.g.,
obtained by RNA-seq):

.. cont-doctest

.. code:: pycon

   >>> rnaseq = "CCCCGGGG"
   >>> sequences2 = [transcript, rnaseq]
   >>> coordinates2 = np.array([[3, 11], [0, 8]])
   >>> alignment2 = Alignment(sequences2, coordinates2)
   >>> print(alignment2)
   target            3 CCCCGGGG 11
                     0 ||||||||  8
   query             0 CCCCGGGG  8
   <BLANKLINE>

Use the ``map`` method on ``alignment1``, with ``alignment2`` as
argument, to find the alignment of the RNA-sequence to the genome:

.. cont-doctest

.. code:: pycon

   >>> alignment3 = alignment1.map(alignment2)
   >>> print(alignment3)
   target           11 CCCCAAAAAAAAAAAGGGG 30
                     0 ||||-----------|||| 19
   query             0 CCCC-----------GGGG  8
   <BLANKLINE>
   >>> print(alignment3.coordinates)
   [[11 15 26 30]
    [ 0  4  4  8]]
   >>> format(alignment3, "psl")
   '8\t0\t0\t0\t0\t0\t1\t11\t+\tquery\t8\t0\t8\ttarget\t40\t11\t30\t2\t4,4,\t0,4,\t11,26,\n'

To be able to print the sequences, in this example we constructed
``alignment1`` and ``alignment2`` using sequences with a defined
sequence contents. However, mapping the alignment does not depend on the
sequence contents; only the coordinates of ``alignment1`` and
``alignment2`` are used to construct the coordinates for ``alignment3``.

The map method can also be used to lift over an alignment between
different genome assemblies. In this case, self is a DNA alignment
between two genome assemblies, and the argument is an alignment of a
transcript against one of the genome assemblies:

.. cont-doctest

.. code:: pycon

   >>> from Bio import Align
   >>> chain = Align.read("Blat/panTro5ToPanTro6.over.chain", "chain")
   >>> chain.sequences[0].id
   'chr1'
   >>> len(chain.sequences[0].seq)
   228573443
   >>> chain.sequences[1].id
   'chr1'
   >>> len(chain.sequences[1].seq)
   224244399
   >>> import numpy as np
   >>> np.set_printoptions(threshold=5)  # print 5 array elements per row
   >>> print(chain.coordinates)
   [[122250000 122250400 122250400 ... 122909818 122909819 122909835]
    [111776384 111776784 111776785 ... 112019962 112019962 112019978]]

showing that the range 122250000:122909835 of chr1 on chimpanzee genome
assembly panTro5 aligns to range 111776384:112019978 of chr1 of
chimpanzee genome assembly panTro6. See section
:ref:`subsec:align_chain` for more information about the chain
file format.

.. cont-doctest

.. code:: pycon

   >>> transcript = Align.read("Blat/est.panTro5.psl", "psl")
   >>> transcript.sequences[0].id
   'chr1'
   >>> len(transcript.sequences[0].seq)
   228573443
   >>> transcript.sequences[1].id
   'DC525629'
   >>> len(transcript.sequences[1].seq)
   407
   >>> print(transcript.coordinates)
   [[122835789 122835847 122840993 122841145 122907212 122907314]
    [       32        90        90       242       242       344]]

This shows that nucleotide range 32:344 of expressed sequence tag
DC525629 aligns to range 122835789:122907314 of chr1 of chimpanzee
genome assembly panTro5. Note that the target sequence
``chain.sequences[0].seq`` and the target sequence
``transcript.sequences[0]`` have the same length:

.. cont-doctest

.. code:: pycon

   >>> len(chain.sequences[0].seq) == len(transcript.sequences[0].seq)
   True

We swap the target and query of the chain such that the query of
``chain`` corresponds to the target of ``transcript``:

.. cont-doctest

.. code:: pycon

   >>> chain = chain[::-1]
   >>> chain.sequences[0].id
   'chr1'
   >>> len(chain.sequences[0].seq)
   224244399
   >>> chain.sequences[1].id
   'chr1'
   >>> len(chain.sequences[1].seq)
   228573443
   >>> print(chain.coordinates)
   [[111776384 111776784 111776785 ... 112019962 112019962 112019978]
    [122250000 122250400 122250400 ... 122909818 122909819 122909835]]
   >>> np.set_printoptions(threshold=1000)  # reset the print options

Now we can get the coordinates of DC525629 against chimpanzee genome
assembly panTro6 by calling ``chain.map``, with ``transcript`` as the
argument:

.. cont-doctest

.. code:: pycon

   >>> lifted_transcript = chain.map(transcript)
   >>> lifted_transcript.sequences[0].id
   'chr1'
   >>> len(lifted_transcript.sequences[0].seq)
   224244399
   >>> lifted_transcript.sequences[1].id
   'DC525629'
   >>> len(lifted_transcript.sequences[1].seq)
   407
   >>> print(lifted_transcript.coordinates)
   [[111982717 111982775 111987921 111988073 112009200 112009302]
    [       32        90        90       242       242       344]]

This shows that nucleotide range 32:344 of expressed sequence tag
DC525629 aligns to range 111982717:112009302 of chr1 of chimpanzee
genome assembly panTro6. Note that the genome span of DC525629 on
chimpanzee genome assembly panTro5 is 122907314 - 122835789 = 71525 bp,
while on panTro6 the genome span is 112009302 - 111982717 = 26585 bp.

Mapping a multiple sequence alignment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Consider a multiple alignment of genomic sequences of chimpanzee, human,
macaque, marmoset, mouse, and rat:

.. cont-doctest

.. code:: pycon

   >>> from Bio import Align
   >>> path = "Blat/panTro5.maf"
   >>> genome_alignment = Align.read(path, "maf")
   >>> for record in genome_alignment.sequences:
   ...     print(record.id, len(record.seq))
   ...
   panTro5.chr1 228573443
   hg19.chr1 249250621
   rheMac8.chr1 225584828
   calJac3.chr18 47448759
   mm10.chr3 160039680
   rn6.chr2 266435125
   >>> print(genome_alignment.coordinates)
   [[133922962 133922962 133922970 133922970 133922972 133922972 133922995
     133922998 133923010]
    [155784573 155784573 155784581 155784581 155784583 155784583 155784606
     155784609 155784621]
    [130383910 130383910 130383918 130383918 130383920 130383920 130383943
     130383946 130383958]
    [  9790455   9790455   9790463   9790463   9790465   9790465   9790488
       9790491   9790503]
    [ 88858039  88858036  88858028  88858026  88858024  88858020  88857997
      88857997  88857985]
    [188162970 188162967 188162959 188162959 188162957 188162953 188162930
     188162930 188162918]]
   >>> print(genome_alignment)
   panTro5.c 133922962 ---ACTAGTTA--CA----GTAACAGAAAATAAAATTTAAATAGAAACTTAAAggcc
   hg19.chr1 155784573 ---ACTAGTTA--CA----GTAACAGAAAATAAAATTTAAATAGAAACTTAAAggcc
   rheMac8.c 130383910 ---ACTAGTTA--CA----GTAACAGAAAATAAAATTTAAATAGAAACTTAAAggcc
   calJac3.c   9790455 ---ACTAGTTA--CA----GTAACAGAAAATAAAATTTAAATAGAAGCTTAAAggct
   mm10.chr3  88858039 TATAATAATTGTATATGTCACAGAAAAAAATGAATTTTCAAT---GACTTAATAGCC
   rn6.chr2  188162970 TACAATAATTG--TATGTCATAGAAAAAAATGAATTTTCAAT---AACTTAATAGCC
   <BLANKLINE>
   panTro5.c 133923010
   hg19.chr1 155784621
   rheMac8.c 130383958
   calJac3.c   9790503
   mm10.chr3  88857985
   rn6.chr2  188162918
   <BLANKLINE>

Suppose we want to replace the older versions of the genome assemblies
(``panTro5``, ``hg19``, ``rheMac8``, ``calJac3``, ``mm10``, and ``rn6``)
by their current versions (``panTro6``, ``hg38``, ``rheMac10``,
``calJac4``, ``mm39``, and ``rn7``). To do so, we need the pairwise
alignment between the old and the new assembly version for each species.
These are provided by UCSC as chain files, typically used for UCSC’s
``liftOver`` tool. The ``.chain`` files in the ``Tests/Align``
subdirectory in the Biopython source distribution were extracted from
UCSC’s ``.chain`` files to only include the relevant genomic region. For
example, to lift over ``panTro5`` to ``panTro6``, we use the file
``panTro5ToPanTro6.chain`` with the following contents:

.. code:: text

   chain 1198066 chr1 228573443 + 133919957 133932620 chr1 224244399 + 130607995 130620657 1
   4990    0   2
   1362    3   0
   6308

To lift over the genome assembly for each species, we read in the
corresponding ``.chain`` file:

.. cont-doctest

.. code:: pycon

   >>> paths = [
   ...     "Blat/panTro5ToPanTro6.chain",
   ...     "Blat/hg19ToHg38.chain",
   ...     "Blat/rheMac8ToRheMac10.chain",
   ...     "Blat/calJac3ToCalJac4.chain",
   ...     "Blat/mm10ToMm39.chain",
   ...     "Blat/rn6ToRn7.chain",
   ... ]
   >>> liftover_alignments = [Align.read(path, "chain") for path in paths]
   >>> for liftover_alignment in liftover_alignments:
   ...     print(liftover_alignment.target.id, liftover_alignment.coordinates[0, :])
   ...
   chr1 [133919957 133924947 133924947 133926309 133926312 133932620]
   chr1 [155184381 156354347 156354348 157128497 157128497 157137496]
   chr1 [130382477 130383872 130383872 130384222 130384222 130388520]
   chr18 [9786631 9787941 9788508 9788508 9795062 9795065 9795737]
   chr3 [66807541 74196805 74196831 94707528 94707528 94708176 94708178 94708718]
   chr2 [188111581 188158351 188158351 188171225 188171225 188228261 188228261
    188236997]

Note that the order of species is the same in ``liftover_alignments``
and ``genome_alignment.sequences``. Now we can lift over the multiple
sequence alignment to the new genome assembly versions:

.. cont-doctest

.. code:: pycon

   >>> genome_alignment = genome_alignment.mapall(liftover_alignments)
   >>> for record in genome_alignment.sequences:
   ...     print(record.id, len(record.seq))
   ...
   chr1 224244399
   chr1 248956422
   chr1 223616942
   chr18 47031477
   chr3 159745316
   chr2 249053267
   >>> print(genome_alignment.coordinates)
   [[130611000 130611000 130611008 130611008 130611010 130611010 130611033
     130611036 130611048]
    [155814782 155814782 155814790 155814790 155814792 155814792 155814815
     155814818 155814830]
    [ 95186253  95186253  95186245  95186245  95186243  95186243  95186220
      95186217  95186205]
    [  9758318   9758318   9758326   9758326   9758328   9758328   9758351
       9758354   9758366]
    [ 88765346  88765343  88765335  88765333  88765331  88765327  88765304
      88765304  88765292]
    [174256702 174256699 174256691 174256691 174256689 174256685 174256662
     174256662 174256650]]

As the ``.chain`` files do not include the sequence contents, we cannot
print the sequence alignment directly. Instead, we read in the genomic
sequence separately (as a ``.2bit`` file, as it allows lazy loading; see
section :ref:`sec:SeqIO_directionaries`) for
each species:

.. code:: pycon

   >>> from Bio import SeqIO
   >>> names = ("panTro6", "hg38", "rheMac10", "calJac4", "mm39", "rn7")
   >>> for i, name in enumerate(names):
   ...     filename = f"{name}.2bit"
   ...     genome = SeqIO.parse(filename, "twobit")
   ...     chromosome = genome_alignment.sequences[i].id
   ...     assert len(genome_alignment.sequences[i]) == len(genome[chromosome])
   ...     genome_alignment.sequences[i] = genome[chromosome]
   ...     genome_alignment.sequences[i].id = f"{name}.{chromosome}"
   ...
   >>> print(genome_alignment)
   panTro6.c 130611000 ---ACTAGTTA--CA----GTAACAGAAAATAAAATTTAAATAGAAACTTAAAggcc
   hg38.chr1 155814782 ---ACTAGTTA--CA----GTAACAGAAAATAAAATTTAAATAGAAACTTAAAggcc
   rheMac10.  95186253 ---ACTAGTTA--CA----GTAACAGAAAATAAAATTTAAATAGAAACTTAAAggcc
   calJac4.c   9758318 ---ACTAGTTA--CA----GTAACAGAaaataaaatttaaatagaagcttaaaggct
   mm39.chr3  88765346 TATAATAATTGTATATGTCACAGAAAAAAATGAATTTTCAAT---GACTTAATAGCC
   rn7.chr2  174256702 TACAATAATTG--TATGTCATAGAAAAAAATGAATTTTCAAT---AACTTAATAGCC
   <BLANKLINE>
   panTro6.c 130611048
   hg38.chr1 155814830
   rheMac10.  95186205
   calJac4.c   9758366
   mm39.chr3  88765292
   rn7.chr2  174256650
   <BLANKLINE>

Using from_alignments_with_same_reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
The :meth:`~Bio.Align.Alignment.from_alignments_with_same_reference` method
constructs a new multiple sequence alignment from a collection of pairwise
alignments that all share the same reference sequence. This is useful when you
have aligned several sequences independently to the same reference and want to
combine those results into a single multiple alignment.

Suppose you have a pairwise alignment of a reference sequence to sequence A and
a multiple alignment of sequences B and C to the same reference. To merge
these into a single multiple alignment of the reference and sequences A, B, and C,
you can use the ``from_alignments_with_same_reference`` method as follows:

.. code:: pycon

   >>> from Bio.Seq import Seq
   >>> from Bio.Align import PairwiseAligner, Alignment
   >>> import numpy as np

   >>> reference_str = "ACGT"
   >>> seq1_str = "ACT"
   >>> seq2_str = "ACGGT"
   >>> seq3_str = "AT"

   >>> aligner = PairwiseAligner()
   >>> alignment_1 = next(aligner.align(reference_str, seq1_str))

   >>> coords = np.array([
   ...     [0, 1, 2, 3, 3, 4],
   ...     [0, 1, 2, 3, 4, 5],
   ...     [0, 1, 1, 1, 1, 2]
   ... ])
   >>> alignment_2 = Alignment([reference_str, seq2_str, seq3_str], coords)

   >>> combined_alignment = Alignment.from_alignments_with_same_reference(
   ...     [alignment_1, alignment_2]
   ... )
   >>> str(combined_alignment[0])
   'ACG-T'
   >>> str(combined_alignment[1])
   'AC--T'
   >>> str(combined_alignment[2])
   'ACGGT'
   >>> str(combined_alignment[3])
   'A---T'

The resulting alignment contains all sequences aligned to the shared reference.
This method differs from :meth:`~Bio.Align.Alignment.map` and
:meth:`~Bio.Align.Alignment.mapall` in that it *builds* a new multiple alignment
directly from a set of pairwise alignments, rather than transforming an existing
multiple alignment using mappings. Each input alignment must have the same
reference sequence (in the same orientation), otherwise an error is raised.

The ``mapall`` method can also be used to create a multiple alignment of
codon sequences from a multiple sequence alignment of the corresponding
amino acid sequences (see Section :ref:`sec:msa_codons`
for details).

.. _`sec:alignments`:

The Alignments class
--------------------

The ``Alignments`` (plural) class inherits from
``AlignmentsAbstractBaseClass`` and from ``list``, and can be used as a
list to store ``Alignment`` objects. The behavior of ``Alignments``
objects is different from that of ``list`` objects in two important
ways:

-  An ``Alignments`` object is its own iterator, consistent with iterators
   returned by ``Bio.Align.parse`` (see section :ref:`subsec:align_reading`) or
   iterators returned by the pairwise aligner (see Section
   :ref:`chapter:pairwise`). Calling ``iter`` on the iterator will
   always return the ``Alignments`` object itself. In contrast, calling
   ``iter`` on a list object creates a new iterator each time, allowing you to
   have multiple independent iterators for a given list.

   In this example, ``alignment_iterator1`` and ``alignment_iterator2`` are
   obtained from a list and act independently of each other:

   .. cont-doctest

   .. code:: pycon

      >>> alignment_list = [alignment1, alignment2, alignment3]
      >>> alignment_iterator1 = iter(alignment_list)
      >>> alignment_iterator2 = iter(alignment_list)
      >>> next(alignment_iterator1)  # doctest: +ELLIPSIS
      <Alignment object (2 rows x 24 columns) at ...>
      >>> next(alignment_iterator2)  # doctest: +ELLIPSIS
      <Alignment object (2 rows x 24 columns) at ...>
      >>> next(alignment_iterator1)  # doctest: +ELLIPSIS
      <Alignment object (2 rows x 8 columns) at ...>
      >>> next(alignment_iterator1)  # doctest: +ELLIPSIS
      <Alignment object (2 rows x 19 columns) at ...>
      >>> next(alignment_iterator2)  # doctest: +ELLIPSIS
      <Alignment object (2 rows x 8 columns) at ...>
      >>> next(alignment_iterator2)  # doctest: +ELLIPSIS
      <Alignment object (2 rows x 19 columns) at ...>

   In contrast, ``alignment_iterator1`` and ``alignment_iterator2`` obtained
   by calling ``iter`` on an ``Alignments`` object are identical to each other:

   .. cont-doctest

   .. code:: pycon

      >>> from Bio.Align import Alignments
      >>> alignments = Alignments([alignment1, alignment2, alignment3])
      >>> alignment_iterator1 = iter(alignments)
      >>> alignment_iterator2 = iter(alignments)
      >>> alignment_iterator1 is alignment_iterator2
      True
      >>> next(alignment_iterator1)  # doctest: +ELLIPSIS
      <Alignment object (2 rows x 24 columns) at ...>
      >>> next(alignment_iterator2)  # doctest: +ELLIPSIS
      <Alignment object (2 rows x 8 columns) at ...>
      >>> next(alignment_iterator1)  # doctest: +ELLIPSIS
      <Alignment object (2 rows x 19 columns) at ...>
      >>> next(alignment_iterator2)
      Traceback (most recent call last):
        File "<stdin>", line 1, in <module>
      StopIteration

   Calling ``iter`` on an ``Alignments`` object resets the iterator to its
   first item, so you can loop over it again. You can also iterate over the
   alignments multiple times using a ``for``-loop, which implicitly calls
   ``iter`` on the iterator:

   .. cont-doctest

   .. code:: pycon

      >>> for item in alignments:
      ...     print(repr(item))  # doctest: +ELLIPSIS
      ...
      <Alignment object (2 rows x 24 columns) at ...>
      <Alignment object (2 rows x 8 columns) at ...>
      <Alignment object (2 rows x 19 columns) at ...>

      >>> for item in alignments:
      ...     print(repr(item))  # doctest: +ELLIPSIS
      ...
      <Alignment object (2 rows x 24 columns) at ...>
      <Alignment object (2 rows x 8 columns) at ...>
      <Alignment object (2 rows x 19 columns) at ...>

   This behavior is consistent with regular Python lists, and with iterators
   returned by ``Bio.Align.parse`` (see section :ref:`subsec:align_reading`) or
   by the pairwise aligner (see Section :ref:`chapter:pairwise`).

-  Metadata can be stored as attributes on an ``Alignments`` object,
   whereas a plain ``list`` does not accept attributes:

   .. cont-doctest

   .. code:: pycon

      >>> alignment_list.score = 100  # doctest: +ELLIPSIS
      Traceback (most recent call last):
       ...
      AttributeError: 'list' object has no attribute 'score'...
      >>> alignments.score = 100
      >>> alignments.score
      100

.. _`sec:alignmentparsers`:

Reading and writing alignments
------------------------------

Output from sequence alignment software such as Clustal can be parsed
into ``Alignment`` objects by the ``Bio.Align.read`` and
``Bio.Align.parse`` functions. Their usage is analogous to the ``read``
and ``parse`` functions in ``Bio.SeqIO`` (see
Section :ref:`sec:Bio.SeqIO-input`): The ``read``
function is used to read an output file containing a single alignment
and returns an ``Alignment`` object, while the ``parse`` function
returns an iterator to iterate over alignments stored in an output file
containing one or more alignments. Section :ref:`sec:alignformats`
describes the alignment formats that can be parsed in ``Bio.Align``.
``Bio.Align`` also provides a ``write`` function that can write
alignments in most of these formats.

.. _`subsec:align_reading`:

Reading alignments
~~~~~~~~~~~~~~~~~~

Use ``Bio.Align.parse`` to parse a file of sequence alignments. For
example, the file ``ucsc_mm9_chr10.maf`` contains 48 multiple sequence
alignments in the MAF (Multiple Alignment Format) format (see section
:ref:`subsec:align_maf`):

.. cont-doctest

.. code:: pycon

   >>> from Bio import Align
   >>> alignments = Align.parse("MAF/ucsc_mm9_chr10.maf", "maf")
   >>> alignments  # doctest: +ELLIPSIS
   <Bio.Align.maf.AlignmentIterator object at 0x...>

where ``"maf"`` is the file format. The alignments object returned by
``Bio.Align.parse`` may contain attributes that store metadata found in
the file, such as the version number of the software that was used to
create the alignments. The specific attributes stored for each file
format are described in Section :ref:`sec:alignformats`. For MAF
files, we can obtain the file format version and the scoring scheme that
was used:

.. cont-doctest

.. code:: pycon

   >>> alignments.metadata
   {'MAF Version': '1', 'Scoring': 'autoMZ.v1'}

As alignment files can be very large, ``Align.parse`` returns an
iterator over the alignments, so you won’t have to store all alignments
in memory at the same time. You can iterate over these alignments and
print out, for example, the number of aligned sequences in each
alignment:

.. cont-doctest

.. code:: pycon

   >>> for a in alignments:
   ...     print(len(a.sequences))  # doctest: +ELLIPSIS
   ...
   2
   4
   5
   6
   ...
   15
   14
   7
   6

You can also call ``len`` on the alignments to obtain the number of
alignments.

.. cont-doctest

.. code:: pycon

   >>> len(alignments)
   48

Depending on the file format, the number of alignments may be explicitly
stored in the file (for example in the case of bigBed, bigPsl, and
bigMaf files), or otherwise the number of alignments is counted by
looping over them once (and returning the iterator to its original
position). If the file is large, it may therefore take a considerable
amount of time for ``len`` to return. However, as the number of
alignments is cached, subsequent calls to ``len`` will return quickly.

If the number of alignments is not excessively large and will fit in
memory, you can convert the alignments iterator to a list of alignments.
To do so, you could call ``list`` on the ``alignments``:

.. cont-doctest

.. code:: pycon

   >>> alignment_list = list(alignments)
   >>> len(alignment_list)
   48
   >>> alignment_list[27]  # doctest: +ELLIPSIS
   <Alignment object (3 rows x 91 columns) at 0x...>
   >>> print(alignment_list[27])
   mm9.chr10   3019377 CCCCAGCATTCTGGCAGACACAGTG-AAAAGAGACAGATGGTCACTAATAAAATCTGT-A
   felCat3.s     46845 CCCAAGTGTTCTGATAGCTAATGTGAAAAAGAAGCATGTGCCCACCAGTAAGCTTTGTGG
   canFam2.c  47545247 CCCAAGTGTTCTGATTGCCTCTGTGAAAAAGAAACATGGGCCCGCTAATAagatttgcaa
   <BLANKLINE>
   mm9.chr10   3019435 TAAATTAG-ATCTCAGAGGATGGATGGACCA  3019465
   felCat3.s     46785 TGAACTAGAATCTCAGAGGATG---GGACTC    46757
   canFam2.c  47545187 tgacctagaatctcagaggatg---ggactc 47545159
   <BLANKLINE>

But this will lose the metadata information:

.. cont-doctest

.. code:: pycon

   >>> alignment_list.metadata  # doctest: +ELLIPSIS
   Traceback (most recent call last):
     ...
   AttributeError: 'list' object has no attribute 'metadata'

Instead, you can ask for a full slice of the alignments:

.. cont-doctest

.. code:: pycon

   >>> type(alignments)
   <class 'Bio.Align.maf.AlignmentIterator'>
   >>> alignments = alignments[:]
   >>> type(alignments)
   <class 'Bio.Align.Alignments'>

This returns a ``Bio.Align.Alignments`` object, which can be used as a
list, while keeping the metadata information:

.. cont-doctest

.. code:: pycon

   >>> len(alignments)
   48
   >>> print(alignments[11])
   mm9.chr10   3014742 AAGTTCCCTCCATAATTCCTTCCTCCCACCCCCACA 3014778
   calJac1.C      6283 AAATGTA-----TGATCTCCCCATCCTGCCCTG---    6311
   otoGar1.s    175262 AGATTTC-----TGATGCCCTCACCCCCTCCGTGCA  175231
   loxAfr1.s      9317 AGGCTTA-----TG----CCACCCCCCACCCCCACA    9290
   <BLANKLINE>
   >>> alignments.metadata
   {'MAF Version': '1', 'Scoring': 'autoMZ.v1'}

.. _`subsec:align_writing`:

Writing alignments
~~~~~~~~~~~~~~~~~~

To write alignments to a file, use

.. code:: pycon

   >>> from Bio import Align
   >>> target = "myfile.txt"
   >>> Align.write(alignments, target, "clustal")

where ``alignments`` is either a single alignment or a list of
alignments, ``target`` is a file name or an open file-like object, and
``"clustal"`` is the file format to be used. As some file formats allow
or require metadata to be stored with the alignments, you may want to
use the ``Alignments`` (plural) class instead of a plain list of
alignments (see Section :ref:`sec:alignments`), allowing you to
store a metadata dictionary as an attribute on the ``alignments``
object:

.. code:: pycon

   >>> from Bio import Align
   >>> alignments = Align.Alignments(alignments)
   >>> metadata = {"Program": "Biopython", "Version": "1.81"}
   >>> alignments.metadata = metadata
   >>> target = "myfile.txt"
   >>> Align.write(alignments, target, "clustal")

.. _`subsec:align_printing`:

Printing alignments
~~~~~~~~~~~~~~~~~~~

For text (non-binary) formats, you can call Python’s built-in ``format``
function on an alignment to get a string showing the alignment in the
requested format, or use ``Alignment`` objects in formatted (f-)
strings. If called without an argument, the ``format`` function returns
the string representation of the alignment:

.. cont-doctest

.. code:: pycon

   >>> str(alignment)
   '                  1 CGGTTTTT 9\n                  0 AGGTTT-- 6\n                  0 AG-TTT-- 5\n'
   >>> format(alignment)
   '                  1 CGGTTTTT 9\n                  0 AGGTTT-- 6\n                  0 AG-TTT-- 5\n'
   >>> print(format(alignment))
                     1 CGGTTTTT 9
                     0 AGGTTT-- 6
                     0 AG-TTT-- 5
   <BLANKLINE>

As optional keyword arguments cannot be used with Python’s built-in
``format`` function or with formatted strings, the ``Alignment`` class
has a ``format`` method with optional arguments to customize the
alignment format. For example, you can use the optional ``scoring`` argument
to provide a substitution matrix (see Section
:ref:`sec:pairwise-substitution-scores`) to let the printed alignment reflect
the substitution scores as follows:

* ``|`` for identical residues,
* ``:`` for substitutions with a positive score,
* ``.`` for substitutions with a negative score,
* ``-`` for gaps.

.. cont-doctest

.. code:: pycon

   >>> M = substitution_matrices.load("NUC.4.4")
   >>> print(M[:, :])
        A    T    G    C    S    W    R    Y    K    M    B    V    H    D    N
   A  5.0 -4.0 -4.0 -4.0 -4.0  1.0  1.0 -4.0 -4.0  1.0 -4.0 -1.0 -1.0 -1.0 -2.0
   T -4.0  5.0 -4.0 -4.0 -4.0  1.0 -4.0  1.0  1.0 -4.0 -1.0 -4.0 -1.0 -1.0 -2.0
   G -4.0 -4.0  5.0 -4.0  1.0 -4.0  1.0 -4.0  1.0 -4.0 -1.0 -1.0 -4.0 -1.0 -2.0
   C -4.0 -4.0 -4.0  5.0  1.0 -4.0 -4.0  1.0 -4.0  1.0 -1.0 -1.0 -1.0 -4.0 -2.0
   S -4.0 -4.0  1.0  1.0 -1.0 -4.0 -2.0 -2.0 -2.0 -2.0 -1.0 -1.0 -3.0 -3.0 -1.0
   W  1.0  1.0 -4.0 -4.0 -4.0 -1.0 -2.0 -2.0 -2.0 -2.0 -3.0 -3.0 -1.0 -1.0 -1.0
   R  1.0 -4.0  1.0 -4.0 -2.0 -2.0 -1.0 -4.0 -2.0 -2.0 -3.0 -1.0 -3.0 -1.0 -1.0
   Y -4.0  1.0 -4.0  1.0 -2.0 -2.0 -4.0 -1.0 -2.0 -2.0 -1.0 -3.0 -1.0 -3.0 -1.0
   K -4.0  1.0  1.0 -4.0 -2.0 -2.0 -2.0 -2.0 -1.0 -4.0 -1.0 -3.0 -3.0 -1.0 -1.0
   M  1.0 -4.0 -4.0  1.0 -2.0 -2.0 -2.0 -2.0 -4.0 -1.0 -3.0 -1.0 -1.0 -3.0 -1.0
   B -4.0 -1.0 -1.0 -1.0 -1.0 -3.0 -3.0 -1.0 -1.0 -3.0 -1.0 -2.0 -2.0 -2.0 -1.0
   V -1.0 -4.0 -1.0 -1.0 -1.0 -3.0 -1.0 -3.0 -3.0 -1.0 -2.0 -1.0 -2.0 -2.0 -1.0
   H -1.0 -1.0 -4.0 -1.0 -3.0 -1.0 -3.0 -1.0 -3.0 -1.0 -2.0 -2.0 -1.0 -2.0 -1.0
   D -1.0 -1.0 -1.0 -4.0 -3.0 -1.0 -1.0 -3.0 -1.0 -3.0 -2.0 -2.0 -2.0 -1.0 -1.0
   N -2.0 -2.0 -2.0 -2.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0
   <BLANKLINE>
   >>> M["T", "Y"]
   1.0
   >>> M["T", "C"]
   -4.0
   >>> aln = Align.Alignment(["GATTACAT", "GATYACAC"])
   >>> print(aln.format(scoring=M))
   target            0 GATTACAT 8
                     0 |||:|||. 8
   query             0 GATYACAC 8
   <BLANKLINE>


Instead of the substitution matrix, you can also use a ``PairwiseAligner``
object (see Chapter :ref:`chapter:pairwise`) as the ``scoring`` argument
to use the substitution matrix associated with the aligner.

By specifying one of the formats shown in
Section :ref:`sec:alignformats`, ``format`` will create a string
showing the alignment in the requested format:

.. cont-doctest

.. code:: pycon

   >>> format(alignment, "clustal")
   'sequence_0                          CGGTTTTT\nsequence_1                          AGGTTT--\nsequence_2                          AG-TTT--\n\n\n'
   >>> print(format(alignment, "clustal"))
   sequence_0                          CGGTTTTT
   sequence_1                          AGGTTT--
   sequence_2                          AG-TTT--
   <BLANKLINE>
   <BLANKLINE>
   <BLANKLINE>
   >>> print(f"*** this is the alignment in Clustal format: ***\n{alignment:clustal}\n***")
   *** this is the alignment in Clustal format: ***
   sequence_0                          CGGTTTTT
   sequence_1                          AGGTTT--
   sequence_2                          AG-TTT--
   <BLANKLINE>
   <BLANKLINE>
   <BLANKLINE>
   ***
   >>> format(alignment, "maf")
   'a\ns sequence_0 1 8 + 9 CGGTTTTT\ns sequence_1 0 6 + 6 AGGTTT--\ns sequence_2 0 5 + 7 AG-TTT--\n\n'
   >>> print(format(alignment, "maf"))
   a
   s sequence_0 1 8 + 9 CGGTTTTT
   s sequence_1 0 6 + 6 AGGTTT--
   s sequence_2 0 5 + 7 AG-TTT--
   <BLANKLINE>
   <BLANKLINE>

As another example, we can print the alignment in BED format (see
section :ref:`subsec:align_bed`) with a specific number of columns:

.. cont-doctest

.. code:: pycon

   >>> print(pairwise_alignment)
   target            1 CGGTTTTT 9
                     0 .|-|||-- 8
   query             0 AG-TTT-- 5
   <BLANKLINE>
   >>> print(format(pairwise_alignment, "bed"))  # doctest: +NORMALIZE_WHITESPACE
   target  1   7   query   0   +   1   7   0   2   2,3,    0,3,
   <BLANKLINE>
   >>> print(pairwise_alignment.format("bed"))  # doctest: +NORMALIZE_WHITESPACE
   target  1   7   query   0   +   1   7   0   2   2,3,    0,3,
   <BLANKLINE>
   >>> print(pairwise_alignment.format("bed", bedN=3))  # doctest: +NORMALIZE_WHITESPACE
   target  1   7
   <BLANKLINE>
   >>> print(pairwise_alignment.format("bed", bedN=6))  # doctest: +NORMALIZE_WHITESPACE
   target  1   7   query   0   +
   <BLANKLINE>

.. _`sec:alignformats`:

Alignment file formats
----------------------

The table below shows the alignment formats that can be parsed in
Bio.Align. The format argument ``fmt`` used in ``Bio.Align`` functions
to specify the file format is case-insensitive. Most of these file
formats can also be written by ``Bio.Align``, as shown in the table.

.. container:: center

   +---------------+-------------+-------------+-------------+-------------+
   | File format   | Description | text /      | Supported   | Subsection  |
   | ``fmt``       |             | binary      | by          |             |
   |               |             |             | ``write``   |             |
   +---------------+-------------+-------------+-------------+-------------+
   | ``a2m``       | A2M         | text        | yes         | `1.7.11     |
   |               |             |             |             | <#subsec:al |
   |               |             |             |             | ign_a2m>`__ |
   +---------------+-------------+-------------+-------------+-------------+
   | ``bed``       | Browser     | text        | yes         | `1.7.14     |
   |               | Extensible  |             |             | <#subsec:al |
   |               | Data (BED)  |             |             | ign_bed>`__ |
   +---------------+-------------+-------------+-------------+-------------+
   | ``bigbed``    | bigBed      | binary      | yes         | `1.7.15 <#s |
   |               |             |             |             | ubsec:align |
   |               |             |             |             | _bigbed>`__ |
   +---------------+-------------+-------------+-------------+-------------+
   | ``bigmaf``    | bigMaf      | binary      | yes         | `1.7.19 <#s |
   |               |             |             |             | ubsec:align |
   |               |             |             |             | _bigmaf>`__ |
   +---------------+-------------+-------------+-------------+-------------+
   | ``bigpsl``    | bigPsl      | binary      | yes         | `1.7.17 <#s |
   |               |             |             |             | ubsec:align |
   |               |             |             |             | _bigpsl>`__ |
   +---------------+-------------+-------------+-------------+-------------+
   | ``chain``     | UCSC chain  | text        | yes         | `1.7.20 <#  |
   |               | file        |             |             | subsec:alig |
   |               |             |             |             | n_chain>`__ |
   +---------------+-------------+-------------+-------------+-------------+
   | ``clustal``   | ClustalW    | text        | yes         | `1.7.2 <#su |
   |               |             |             |             | bsec:align_ |
   |               |             |             |             | clustal>`__ |
   +---------------+-------------+-------------+-------------+-------------+
   | ``emboss``    | EMBOSS      | text        | no          | `1.7.5 <#s  |
   |               |             |             |             | ubsec:align |
   |               |             |             |             | _emboss>`__ |
   +---------------+-------------+-------------+-------------+-------------+
   | ``exonerate`` | Exonerate   | text        | yes         | `1          |
   |               |             |             |             | .7.7 <#subs |
   |               |             |             |             | ec:align_ex |
   |               |             |             |             | onerate>`__ |
   +---------------+-------------+-------------+-------------+-------------+
   | ``fasta``     | Aligned     | text        | yes         | `1.7.1 <#   |
   |               | FASTA       |             |             | subsec:alig |
   |               |             |             |             | n_fasta>`__ |
   +---------------+-------------+-------------+-------------+-------------+
   | ``hhr``       | HH-suite    | text        | no          | `1.7.10     |
   |               | output      |             |             | <#subsec:al |
   |               | files       |             |             | ign_hhr>`__ |
   +---------------+-------------+-------------+-------------+-------------+
   | ``maf``       | Multiple    | text        | yes         | `1.7.18     |
   |               | Alignment   |             |             | <#subsec:al |
   |               | Format      |             |             | ign_maf>`__ |
   |               | (MAF)       |             |             |             |
   +---------------+-------------+-------------+-------------+-------------+
   | ``mauve``     | Mauve       | text        | yes         | `1.7.12 <#  |
   |               | eXtended    |             |             | subsec:alig |
   |               | Multi-FastA |             |             | n_mauve>`__ |
   |               | (xmfa)      |             |             |             |
   |               | format      |             |             |             |
   +---------------+-------------+-------------+-------------+-------------+
   | ``msf``       | GCG         | text        | no          | `1.7.6      |
   |               | Multiple    |             |             | <#subsec:al |
   |               | Sequence    |             |             | ign_msf>`__ |
   |               | Format      |             |             |             |
   |               | (MSF)       |             |             |             |
   +---------------+-------------+-------------+-------------+-------------+
   | ``nexus``     | NEXUS       | text        | yes         | `1.7.8 <#   |
   |               |             |             |             | subsec:alig |
   |               |             |             |             | n_nexus>`__ |
   +---------------+-------------+-------------+-------------+-------------+
   | ``phylip``    | PHYLIP      | text        | yes         | `1.7.4 <#s  |
   |               | output      |             |             | ubsec:align |
   |               | files       |             |             | _phylip>`__ |
   +---------------+-------------+-------------+-------------+-------------+
   | ``psl``       | Pattern     | text        | yes         | `1.7.16     |
   |               | Space       |             |             | <#subsec:al |
   |               | Layout      |             |             | ign_psl>`__ |
   |               | (PSL)       |             |             |             |
   +---------------+-------------+-------------+-------------+-------------+
   | ``sam``       | Sequence    | text        | yes         | `1.7.13     |
   |               | Alignment/  |             |             | <#subsec:al |
   |               | Map (SAM)   |             |             | ign_sam>`__ |
   +---------------+-------------+-------------+-------------+-------------+
   | ``stockholm`` | Stockholm   | text        | yes         | `1          |
   |               |             |             |             | .7.3 <#subs |
   |               |             |             |             | ec:align_st |
   |               |             |             |             | ockholm>`__ |
   +---------------+-------------+-------------+-------------+-------------+
   | ``tabular``   | Tabular     | text        | no          | `1.7.9 <#su |
   |               | output from |             |             | bsec:align_ |
   |               | BLAST or    |             |             | tabular>`__ |
   |               | FASTA       |             |             |             |
   +---------------+-------------+-------------+-------------+-------------+

.. _`subsec:align_fasta`:

Aligned FASTA
~~~~~~~~~~~~~

Files in the aligned FASTA format store exactly one (pairwise or
multiple) sequence alignment, in which gaps in the alignment are
represented by dashes (``-``). Use ``fmt="fasta"`` to read or write
files in the aligned FASTA format. Note that this is different from
output generated by William Pearson’s FASTA alignment program (parsing
such output is described in section :ref:`subsec:align_tabular`
instead).

The file ``probcons.fa`` in Biopython’s test suite stores one multiple
alignment in the aligned FASTA format. The contents of this file is as
follows:

.. code:: text

   >plas_horvu
   D-VLLGANGGVLVFEPNDFSVKAGETITFKNNAGYPHNVVFDEDAVPSG-VD-VSKISQEEYLTAPGETFSVTLTV---PGTYGFYCEPHAGAGMVGKVTV
   >plas_chlre
   --VKLGADSGALEFVPKTLTIKSGETVNFVNNAGFPHNIVFDEDAIPSG-VN-ADAISRDDYLNAPGETYSVKLTA---AGEYGYYCEPHQGAGMVGKIIV
   >plas_anava
   --VKLGSDKGLLVFEPAKLTIKPGDTVEFLNNKVPPHNVVFDAALNPAKSADLAKSLSHKQLLMSPGQSTSTTFPADAPAGEYTFYCEPHRGAGMVGKITV
   >plas_proho
   VQIKMGTDKYAPLYEPKALSISAGDTVEFVMNKVGPHNVIFDK--VPAG-ES-APALSNTKLRIAPGSFYSVTLGT---PGTYSFYCTPHRGAGMVGTITV
   >azup_achcy
   VHMLNKGKDGAMVFEPASLKVAPGDTVTFIPTDK-GHNVETIKGMIPDG-AE-A-------FKSKINENYKVTFTA---PGVYGVKCTPHYGMGMVGVVEV

To read this file, use

.. doctest ../Tests/Clustalw lib:numpy

.. code:: pycon

   >>> from Bio import Align
   >>> alignment = Align.read("probcons.fa", "fasta")
   >>> alignment  # doctest: +ELLIPSIS
   <Alignment object (5 rows x 101 columns) at ...>

We can print the alignment to see its default representation:

.. cont-doctest

.. code:: pycon

   >>> print(alignment)
   plas_horv         0 D-VLLGANGGVLVFEPNDFSVKAGETITFKNNAGYPHNVVFDEDAVPSG-VD-VSKISQE
   plas_chlr         0 --VKLGADSGALEFVPKTLTIKSGETVNFVNNAGFPHNIVFDEDAIPSG-VN-ADAISRD
   plas_anav         0 --VKLGSDKGLLVFEPAKLTIKPGDTVEFLNNKVPPHNVVFDAALNPAKSADLAKSLSHK
   plas_proh         0 VQIKMGTDKYAPLYEPKALSISAGDTVEFVMNKVGPHNVIFDK--VPAG-ES-APALSNT
   azup_achc         0 VHMLNKGKDGAMVFEPASLKVAPGDTVTFIPTDK-GHNVETIKGMIPDG-AE-A------
   <BLANKLINE>
   plas_horv        57 EYLTAPGETFSVTLTV---PGTYGFYCEPHAGAGMVGKVTV 95
   plas_chlr        56 DYLNAPGETYSVKLTA---AGEYGYYCEPHQGAGMVGKIIV 94
   plas_anav        58 QLLMSPGQSTSTTFPADAPAGEYTFYCEPHRGAGMVGKITV 99
   plas_proh        56 KLRIAPGSFYSVTLGT---PGTYSFYCTPHRGAGMVGTITV 94
   azup_achc        51 -FKSKINENYKVTFTA---PGVYGVKCTPHYGMGMVGVVEV 88
   <BLANKLINE>

or we can print it in the aligned FASTA format:

.. cont-doctest

.. code:: pycon

   >>> print(format(alignment, "fasta"))
   >plas_horvu
   D-VLLGANGGVLVFEPNDFSVKAGETITFKNNAGYPHNVVFDEDAVPSG-VD-VSKISQEEYLTAPGETFSVTLTV---PGTYGFYCEPHAGAGMVGKVTV
   >plas_chlre
   --VKLGADSGALEFVPKTLTIKSGETVNFVNNAGFPHNIVFDEDAIPSG-VN-ADAISRDDYLNAPGETYSVKLTA---AGEYGYYCEPHQGAGMVGKIIV
   >plas_anava
   --VKLGSDKGLLVFEPAKLTIKPGDTVEFLNNKVPPHNVVFDAALNPAKSADLAKSLSHKQLLMSPGQSTSTTFPADAPAGEYTFYCEPHRGAGMVGKITV
   >plas_proho
   VQIKMGTDKYAPLYEPKALSISAGDTVEFVMNKVGPHNVIFDK--VPAG-ES-APALSNTKLRIAPGSFYSVTLGT---PGTYSFYCTPHRGAGMVGTITV
   >azup_achcy
   VHMLNKGKDGAMVFEPASLKVAPGDTVTFIPTDK-GHNVETIKGMIPDG-AE-A-------FKSKINENYKVTFTA---PGVYGVKCTPHYGMGMVGVVEV
   <BLANKLINE>

or any other available format, for example Clustal (see
section :ref:`subsec:align_clustal`):

.. cont-doctest

.. code:: pycon

   >>> print(format(alignment, "clustal"))
   plas_horvu                          D-VLLGANGGVLVFEPNDFSVKAGETITFKNNAGYPHNVVFDEDAVPSG-
   plas_chlre                          --VKLGADSGALEFVPKTLTIKSGETVNFVNNAGFPHNIVFDEDAIPSG-
   plas_anava                          --VKLGSDKGLLVFEPAKLTIKPGDTVEFLNNKVPPHNVVFDAALNPAKS
   plas_proho                          VQIKMGTDKYAPLYEPKALSISAGDTVEFVMNKVGPHNVIFDK--VPAG-
   azup_achcy                          VHMLNKGKDGAMVFEPASLKVAPGDTVTFIPTDK-GHNVETIKGMIPDG-
   <BLANKLINE>
   plas_horvu                          VD-VSKISQEEYLTAPGETFSVTLTV---PGTYGFYCEPHAGAGMVGKVT
   plas_chlre                          VN-ADAISRDDYLNAPGETYSVKLTA---AGEYGYYCEPHQGAGMVGKII
   plas_anava                          ADLAKSLSHKQLLMSPGQSTSTTFPADAPAGEYTFYCEPHRGAGMVGKIT
   plas_proho                          ES-APALSNTKLRIAPGSFYSVTLGT---PGTYSFYCTPHRGAGMVGTIT
   azup_achcy                          AE-A-------FKSKINENYKVTFTA---PGVYGVKCTPHYGMGMVGVVE
   <BLANKLINE>
   plas_horvu                          V
   plas_chlre                          V
   plas_anava                          V
   plas_proho                          V
   azup_achcy                          V
   <BLANKLINE>
   <BLANKLINE>
   <BLANKLINE>

The sequences associated with the alignment are ``SeqRecord`` objects:

.. cont-doctest

.. code:: pycon

   >>> alignment.sequences
   [SeqRecord(seq=Seq('DVLLGANGGVLVFEPNDFSVKAGETITFKNNAGYPHNVVFDEDAVPSGVDVSKI...VTV'), id='plas_horvu', name='<unknown name>', description='', dbxrefs=[]), SeqRecord(seq=Seq('VKLGADSGALEFVPKTLTIKSGETVNFVNNAGFPHNIVFDEDAIPSGVNADAIS...IIV'), id='plas_chlre', name='<unknown name>', description='', dbxrefs=[]), SeqRecord(seq=Seq('VKLGSDKGLLVFEPAKLTIKPGDTVEFLNNKVPPHNVVFDAALNPAKSADLAKS...ITV'), id='plas_anava', name='<unknown name>', description='', dbxrefs=[]), SeqRecord(seq=Seq('VQIKMGTDKYAPLYEPKALSISAGDTVEFVMNKVGPHNVIFDKVPAGESAPALS...ITV'), id='plas_proho', name='<unknown name>', description='', dbxrefs=[]), SeqRecord(seq=Seq('VHMLNKGKDGAMVFEPASLKVAPGDTVTFIPTDKGHNVETIKGMIPDGAEAFKS...VEV'), id='azup_achcy', name='<unknown name>', description='', dbxrefs=[])]

Note that these sequences do not contain gaps ("``-``" characters), as
the alignment information is stored in the ``coordinates`` attribute
instead:

.. cont-doctest

.. code:: pycon

   >>> print(alignment.coordinates)
   [[ 0  1  1 33 34 42 44 48 48 50 50 51 58 73 73 95]
    [ 0  0  0 32 33 41 43 47 47 49 49 50 57 72 72 94]
    [ 0  0  0 32 33 41 43 47 48 50 51 52 59 74 77 99]
    [ 0  1  2 34 35 43 43 47 47 49 49 50 57 72 72 94]
    [ 0  1  2 34 34 42 44 48 48 50 50 51 51 66 66 88]]

Use ``Align.write`` to write this alignment to a file (here, we’ll use a
``StringIO`` object instead of a file):

.. cont-doctest

.. code:: pycon

   >>> from io import StringIO
   >>> stream = StringIO()
   >>> Align.write(alignment, stream, "FASTA")
   1
   >>> print(stream.getvalue())
   >plas_horvu
   D-VLLGANGGVLVFEPNDFSVKAGETITFKNNAGYPHNVVFDEDAVPSG-VD-VSKISQEEYLTAPGETFSVTLTV---PGTYGFYCEPHAGAGMVGKVTV
   >plas_chlre
   --VKLGADSGALEFVPKTLTIKSGETVNFVNNAGFPHNIVFDEDAIPSG-VN-ADAISRDDYLNAPGETYSVKLTA---AGEYGYYCEPHQGAGMVGKIIV
   >plas_anava
   --VKLGSDKGLLVFEPAKLTIKPGDTVEFLNNKVPPHNVVFDAALNPAKSADLAKSLSHKQLLMSPGQSTSTTFPADAPAGEYTFYCEPHRGAGMVGKITV
   >plas_proho
   VQIKMGTDKYAPLYEPKALSISAGDTVEFVMNKVGPHNVIFDK--VPAG-ES-APALSNTKLRIAPGSFYSVTLGT---PGTYSFYCTPHRGAGMVGTITV
   >azup_achcy
   VHMLNKGKDGAMVFEPASLKVAPGDTVTFIPTDK-GHNVETIKGMIPDG-AE-A-------FKSKINENYKVTFTA---PGVYGVKCTPHYGMGMVGVVEV
   <BLANKLINE>

Note that ``Align.write`` returns the number of alignments written (1,
in this case).

.. _`subsec:align_clustal`:

ClustalW
~~~~~~~~

Clustal is a set of multiple sequence alignment programs that are
available both as standalone programs as as web servers. The file
``opuntia.aln`` (available online or in the ``Doc/examples``
subdirectory of the Biopython source code) is an output file generated
by Clustal. Its first few lines are

.. code:: text

   CLUSTAL 2.1 multiple sequence alignment


   gi|6273285|gb|AF191659.1|AF191      TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
   gi|6273284|gb|AF191658.1|AF191      TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
   gi|6273287|gb|AF191661.1|AF191      TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
   gi|6273286|gb|AF191660.1|AF191      TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
   gi|6273290|gb|AF191664.1|AF191      TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
   gi|6273289|gb|AF191663.1|AF191      TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
   gi|6273291|gb|AF191665.1|AF191      TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
                                       ******* **** *************************************

   ...

To parse this file, use

.. doctest examples lib:numpy

.. code:: pycon

   >>> from Bio import Align
   >>> alignments = Align.parse("opuntia.aln", "clustal")

The ``metadata`` attribute on ``alignments`` stores the information
shown in the file header:

.. cont-doctest

.. code:: pycon

   >>> alignments.metadata
   {'Program': 'CLUSTAL', 'Version': '2.1'}

You can call ``next`` on the ``alignments`` to pull out the first (and
only) alignment:

.. cont-doctest

.. code:: pycon

   >>> alignment = next(alignments)
   >>> print(alignment)  # doctest: +ELLIPSIS
   gi|627328         0 TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAAAAAAATGAAT
   gi|627328         0 TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAAAAAAATGAAT
   gi|627328         0 TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAAAAAAATGAAT
   gi|627328         0 TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAAAAAAATGAAT
   gi|627329         0 TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAAAAAAATGAAT
   gi|627328         0 TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAAAAAAATGAAT
   gi|627329         0 TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAAAAAAATGAAT
   <BLANKLINE>
   gi|627328        60 CTAAATGATATACGATTCCACTATGTAAGGTCTTTGAATCATATCATAAAAGACAATGTA
   gi|627328        60 CTAAATGATATACGATTCCACTATGTAAGGTCTTTGAATCATATCATAAAAGACAATGTA
   gi|627328        60 CTAAATGATATACGATTCCACTATGTAAGGTCTTTGAATCATATCATAAAAGACAATGTA
   gi|627328        60 CTAAATGATATACGATTCCACTA...

If you are not interested in the metadata, then it is more convenient to
use the ``Align.read`` function, as anyway each Clustal file contains
only one alignment:

.. cont-doctest

.. code:: pycon

   >>> from Bio import Align
   >>> alignment = Align.read("opuntia.aln", "clustal")

The consensus line below each alignment block in the Clustal output file
contains an asterisk if the sequence is conserved at each position. This
information is stored in the ``column_annotations`` attribute of the
``alignment``:

.. cont-doctest

.. code:: pycon

   >>> alignment.column_annotations  # doctest: +ELLIPSIS
   {'clustal_consensus': '******* **** **********************************...

Printing the ``alignment`` in ``clustal`` format will show the sequence
alignment, but does not include the metadata:

.. cont-doctest

.. code:: pycon

   >>> print(format(alignment, "clustal"))  # doctest: +ELLIPSIS
   gi|6273285|gb|AF191659.1|AF191      TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
   gi|6273284|gb|AF191658.1|AF191      TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
   gi|6273287|gb|AF191661.1|AF191      TATACATT...

Writing the ``alignments`` in ``clustal`` format will include both the
metadata and the sequence alignment:

.. cont-doctest

.. code:: pycon

   >>> from io import StringIO
   >>> stream = StringIO()
   >>> Align.write(alignments, stream, "clustal")
   1
   >>> print(stream.getvalue())  # doctest: +ELLIPSIS
   CLUSTAL 2.1 multiple sequence alignment
   <BLANKLINE>
   <BLANKLINE>
   gi|6273285|gb|AF191659.1|AF191      TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
   gi|6273284|gb|AF191658.1|AF191      TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
   gi|6273287|gb|AF191661.1|AF191      TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
   gi|6273286|gb|AF191660.1|AF191      TATACATAAAAGAAG...

Use an ``Alignments`` (plural) object (see
Section :ref:`sec:alignments`) if you are creating alignments by
hand, and would like to include metadata information in the output.

.. _`subsec:align_stockholm`:

Stockholm
~~~~~~~~~

This is an example of a protein sequence alignment in the Stockholm file
format used by PFAM:

.. code:: text

   # STOCKHOLM 1.0
   #=GF ID   7kD_DNA_binding
   #=GF AC   PF02294.20
   #=GF DE   7kD DNA-binding domain
   #=GF AU   Mian N;0000-0003-4284-4749
   #=GF AU   Bateman A;0000-0002-6982-4660
   #=GF SE   Pfam-B_8148 (release 5.2)
   #=GF GA   25.00 25.00;
   #=GF TC   26.60 46.20;
   #=GF NC   23.20 19.20;
   #=GF BM   hmmbuild HMM.ann SEED.ann
   #=GF SM   hmmsearch -Z 57096847 -E 1000 --cpu 4 HMM pfamseq
   #=GF TP   Domain
   #=GF CL   CL0049
   #=GF RN   [1]
   #=GF RM   3130377
   #=GF RT   Microsequence analysis of DNA-binding proteins 7a, 7b, and 7e
   #=GF RT   from the archaebacterium Sulfolobus acidocaldarius.
   #=GF RA   Choli T, Wittmann-Liebold B, Reinhardt R;
   #=GF RL   J Biol Chem 1988;263:7087-7093.
   #=GF DR   INTERPRO; IPR003212;
   #=GF DR   SCOP; 1sso; fa;
   #=GF DR   SO; 0000417; polypeptide_domain;
   #=GF CC   This family contains members of the hyper-thermophilic
   #=GF CC   archaebacterium  7kD DNA-binding/endoribonuclease P2 family.
   #=GF CC   There are five 7kD DNA-binding proteins, 7a-7e, found as
   #=GF CC   monomers in the cell. Protein 7e shows the  tightest DNA-binding
   #=GF CC   ability.
   #=GF SQ   3
   #=GS DN7_METS5/4-61   AC A4YEA2.1
   #=GS DN7A_SACS2/3-61  AC P61991.2
   #=GS DN7A_SACS2/3-61  DR PDB; 1SSO A; 2-60;
   #=GS DN7A_SACS2/3-61  DR PDB; 1JIC A; 2-60;
   #=GS DN7A_SACS2/3-61  DR PDB; 2CVR A; 2-60;
   #=GS DN7A_SACS2/3-61  DR PDB; 1B4O A; 2-60;
   #=GS DN7E_SULAC/3-60  AC P13125.2
   DN7_METS5/4-61              KIKFKYKGQDLEVDISKVKKVWKVGKMVSFTYDD.NGKTGRGAVSEKDAPKELLNMIGK
   DN7A_SACS2/3-61             TVKFKYKGEEKQVDISKIKKVWRVGKMISFTYDEGGGKTGRGAVSEKDAPKELLQMLEK
   #=GR DN7A_SACS2/3-61  SS    EEEEESSSSEEEEETTTEEEEEESSSSEEEEEE-SSSSEEEEEEETTTS-CHHHHHHTT
   DN7E_SULAC/3-60             KVRFKYKGEEKEVDTSKIKKVWRVGKMVSFTYDD.NGKTGRGAVSEKDAPKELMDMLAR
   #=GC SS_cons                EEEEESSSSEEEEETTTEEEEEESSSSEEEEEE-SSSSEEEEEEETTTS-CHHHHHHTT
   #=GC seq_cons               KVKFKYKGEEKEVDISKIKKVWRVGKMVSFTYDD.NGKTGRGAVSEKDAPKELLsMLuK
   //

This is the seed alignment for the 7kD_DNA_binding (PF02294.20) PFAM
entry, downloaded from the InterPro website
(https://www.ebi.ac.uk/interpro/). This version of the PFAM entry is
also available in the Biopython source distribution as the file
``pfam2.seed.txt`` in the subdirectory ``Tests/Stockholm/``. We can load
this file as follows:

.. doctest ../Tests/Stockholm lib:numpy

.. code:: pycon

   >>> from Bio import Align
   >>> alignment = Align.read("pfam2.seed.txt", "stockholm")
   >>> alignment  # doctest: +ELLIPSIS
   <Alignment object (3 rows x 59 columns) at ...>

We can print out a summary of the alignment:

.. cont-doctest

.. code:: pycon

   >>> print(alignment)
   DN7_METS5         0 KIKFKYKGQDLEVDISKVKKVWKVGKMVSFTYDD-NGKTGRGAVSEKDAPKELLNMIGK
   DN7A_SACS         0 TVKFKYKGEEKQVDISKIKKVWRVGKMISFTYDEGGGKTGRGAVSEKDAPKELLQMLEK
   DN7E_SULA         0 KVRFKYKGEEKEVDTSKIKKVWRVGKMVSFTYDD-NGKTGRGAVSEKDAPKELMDMLAR
   <BLANKLINE>
   DN7_METS5        58
   DN7A_SACS        59
   DN7E_SULA        58
   <BLANKLINE>

You could also call Python’s built-in ``format`` function on the
alignment object to show it in a particular file format (see
section :ref:`subsec:align_printing` for details), for example in
the Stockholm format to regenerate the file:

.. cont-doctest

.. code:: pycon

   >>> print(format(alignment, "stockholm"))
   # STOCKHOLM 1.0
   #=GF ID   7kD_DNA_binding
   #=GF AC   PF02294.20
   #=GF DE   7kD DNA-binding domain
   #=GF AU   Mian N;0000-0003-4284-4749
   #=GF AU   Bateman A;0000-0002-6982-4660
   #=GF SE   Pfam-B_8148 (release 5.2)
   #=GF GA   25.00 25.00;
   #=GF TC   26.60 46.20;
   #=GF NC   23.20 19.20;
   #=GF BM   hmmbuild HMM.ann SEED.ann
   #=GF SM   hmmsearch -Z 57096847 -E 1000 --cpu 4 HMM pfamseq
   #=GF TP   Domain
   #=GF CL   CL0049
   #=GF RN   [1]
   #=GF RM   3130377
   #=GF RT   Microsequence analysis of DNA-binding proteins 7a, 7b, and 7e from
   #=GF RT   the archaebacterium Sulfolobus acidocaldarius.
   #=GF RA   Choli T, Wittmann-Liebold B, Reinhardt R;
   #=GF RL   J Biol Chem 1988;263:7087-7093.
   #=GF DR   INTERPRO; IPR003212;
   #=GF DR   SCOP; 1sso; fa;
   #=GF DR   SO; 0000417; polypeptide_domain;
   #=GF CC   This family contains members of the hyper-thermophilic
   #=GF CC   archaebacterium  7kD DNA-binding/endoribonuclease P2 family. There
   #=GF CC   are five 7kD DNA-binding proteins, 7a-7e, found as monomers in the
   #=GF CC   cell. Protein 7e shows the  tightest DNA-binding ability.
   #=GF SQ   3
   #=GS DN7_METS5/4-61   AC A4YEA2.1
   #=GS DN7A_SACS2/3-61  AC P61991.2
   #=GS DN7A_SACS2/3-61  DR PDB; 1SSO A; 2-60;
   #=GS DN7A_SACS2/3-61  DR PDB; 1JIC A; 2-60;
   #=GS DN7A_SACS2/3-61  DR PDB; 2CVR A; 2-60;
   #=GS DN7A_SACS2/3-61  DR PDB; 1B4O A; 2-60;
   #=GS DN7E_SULAC/3-60  AC P13125.2
   DN7_METS5/4-61                  KIKFKYKGQDLEVDISKVKKVWKVGKMVSFTYDD.NGKTGRGAVSEKDAPKELLNMIGK
   DN7A_SACS2/3-61                 TVKFKYKGEEKQVDISKIKKVWRVGKMISFTYDEGGGKTGRGAVSEKDAPKELLQMLEK
   #=GR DN7A_SACS2/3-61  SS        EEEEESSSSEEEEETTTEEEEEESSSSEEEEEE-SSSSEEEEEEETTTS-CHHHHHHTT
   DN7E_SULAC/3-60                 KVRFKYKGEEKEVDTSKIKKVWRVGKMVSFTYDD.NGKTGRGAVSEKDAPKELMDMLAR
   #=GC SS_cons                    EEEEESSSSEEEEETTTEEEEEESSSSEEEEEE-SSSSEEEEEEETTTS-CHHHHHHTT
   #=GC seq_cons                   KVKFKYKGEEKEVDISKIKKVWRVGKMVSFTYDD.NGKTGRGAVSEKDAPKELLsMLuK
   //
   <BLANKLINE>

or alternatively as aligned FASTA (see section
:ref:`subsec:align_fasta`):

.. cont-doctest

.. code:: pycon

   >>> print(format(alignment, "fasta"))
   >DN7_METS5/4-61
   KIKFKYKGQDLEVDISKVKKVWKVGKMVSFTYDD-NGKTGRGAVSEKDAPKELLNMIGK
   >DN7A_SACS2/3-61
   TVKFKYKGEEKQVDISKIKKVWRVGKMISFTYDEGGGKTGRGAVSEKDAPKELLQMLEK
   >DN7E_SULAC/3-60
   KVRFKYKGEEKEVDTSKIKKVWRVGKMVSFTYDD-NGKTGRGAVSEKDAPKELMDMLAR
   <BLANKLINE>

or in the PHYLIP format (see section :ref:`subsec:align_phylip`):

.. cont-doctest

.. code:: pycon

   >>> print(format(alignment, "phylip"))
   3 59
   DN7_METS5/KIKFKYKGQDLEVDISKVKKVWKVGKMVSFTYDD-NGKTGRGAVSEKDAPKELLNMIGK
   DN7A_SACS2TVKFKYKGEEKQVDISKIKKVWRVGKMISFTYDEGGGKTGRGAVSEKDAPKELLQMLEK
   DN7E_SULACKVRFKYKGEEKEVDTSKIKKVWRVGKMVSFTYDD-NGKTGRGAVSEKDAPKELMDMLAR
   <BLANKLINE>

General information of the alignment is stored under the ``annotations``
attribute of the ``Alignment`` object, for example

.. cont-doctest

.. code:: pycon

   >>> alignment.annotations["identifier"]
   '7kD_DNA_binding'
   >>> alignment.annotations["clan"]
   'CL0049'
   >>> alignment.annotations["database references"]
   [{'reference': 'INTERPRO; IPR003212;'}, {'reference': 'SCOP; 1sso; fa;'}, {'reference': 'SO; 0000417; polypeptide_domain;'}]

The individual sequences in this alignment are stored under
``alignment.sequences`` as ``SeqRecord``\ s, including any annotations
associated with each sequence record:

.. cont-doctest

.. code:: pycon

   >>> for record in alignment.sequences:
   ...     print("%s %s %s" % (record.id, record.annotations["accession"], record.dbxrefs))
   ...
   DN7_METS5/4-61 A4YEA2.1 []
   DN7A_SACS2/3-61 P61991.2 ['PDB; 1SSO A; 2-60;', 'PDB; 1JIC A; 2-60;', 'PDB; 2CVR A; 2-60;', 'PDB; 1B4O A; 2-60;']
   DN7E_SULAC/3-60 P13125.2 []

The secondary structure of the second sequence (``DN7A_SACS2/3-61``) is
stored in the ``letter_annotations`` attribute of the ``SeqRecord``:

.. cont-doctest

.. code:: pycon

   >>> alignment.sequences[0].letter_annotations
   {}
   >>> alignment.sequences[1].letter_annotations
   {'secondary structure': 'EEEEESSSSEEEEETTTEEEEEESSSSEEEEEE-SSSSEEEEEEETTTS-CHHHHHHTT'}
   >>> alignment.sequences[2].letter_annotations
   {}

The consensus sequence and secondary structure are associated with the
sequence alignment as a whole, and are therefore stored in the
``column_annotations`` attribute of the ``Alignment`` object:

.. cont-doctest

.. code:: pycon

   >>> alignment.column_annotations  # doctest: +NORMALIZE_WHITESPACE
   {'consensus secondary structure': 'EEEEESSSSEEEEETTTEEEEEESSSSEEEEEE-SSSSEEEEEEETTTS-CHHHHHHTT',
    'consensus sequence': 'KVKFKYKGEEKEVDISKIKKVWRVGKMVSFTYDD.NGKTGRGAVSEKDAPKELLsMLuK'}

.. _`subsec:align_phylip`:

PHYLIP output files
~~~~~~~~~~~~~~~~~~~

The PHYLIP format for sequence alignments is derived from the PHYLogeny
Interference Package from Joe Felsenstein. Files in the PHYLIP format
start with two numbers for the number of rows and columns in the printed
alignment. The sequence alignment itself can be in sequential format or
in interleaved format. An example of the former is the
``sequential.phy`` file (provided in ``Tests/Phylip/`` in the Biopython
source distribution):

.. code:: text

    3 384
   CYS1_DICDI   -----MKVIL LFVLAVFTVF VSS------- --------RG IPPEEQ---- --------SQ
                FLEFQDKFNK KY-SHEEYLE RFEIFKSNLG KIEELNLIAI NHKADTKFGV NKFADLSSDE
                FKNYYLNNKE AIFTDDLPVA DYLDDEFINS IPTAFDWRTR G-AVTPVKNQ GQCGSCWSFS
                TTGNVEGQHF ISQNKLVSLS EQNLVDCDHE CMEYEGEEAC DEGCNGGLQP NAYNYIIKNG
                GIQTESSYPY TAETGTQCNF NSANIGAKIS NFTMIP-KNE TVMAGYIVST GPLAIAADAV
                E-WQFYIGGV F-DIPCN--P NSLDHGILIV GYSAKNTIFR KNMPYWIVKN SWGADWGEQG
                YIYLRRGKNT CGVSNFVSTS II--
   ALEU_HORVU   MAHARVLLLA LAVLATAAVA VASSSSFADS NPIRPVTDRA ASTLESAVLG ALGRTRHALR
                FARFAVRYGK SYESAAEVRR RFRIFSESLE EVRSTN---- RKGLPYRLGI NRFSDMSWEE
                FQATRL-GAA QTCSATLAGN HLMRDA--AA LPETKDWRED G-IVSPVKNQ AHCGSCWTFS
                TTGALEAAYT QATGKNISLS EQQLVDCAGG FNNF------ --GCNGGLPS QAFEYIKYNG
                GIDTEESYPY KGVNGV-CHY KAENAAVQVL DSVNITLNAE DELKNAVGLV RPVSVAFQVI
                DGFRQYKSGV YTSDHCGTTP DDVNHAVLAV GYGVENGV-- ---PYWLIKN SWGADWGDNG
                YFKMEMGKNM CAIATCASYP VVAA
   CATH_HUMAN   ------MWAT LPLLCAGAWL LGV------- -PVCGAAELS VNSLEK---- --------FH
                FKSWMSKHRK TY-STEEYHH RLQTFASNWR KINAHN---- NGNHTFKMAL NQFSDMSFAE
                IKHKYLWSEP QNCSAT--KS NYLRGT--GP YPPSVDWRKK GNFVSPVKNQ GACGSCWTFS
                TTGALESAIA IATGKMLSLA EQQLVDCAQD FNNY------ --GCQGGLPS QAFEYILYNK
                GIMGEDTYPY QGKDGY-CKF QPGKAIGFVK DVANITIYDE EAMVEAVALY NPVSFAFEVT
                QDFMMYRTGI YSSTSCHKTP DKVNHAVLAV GYGEKNGI-- ---PYWIVKN SWGPQWGMNG
                YFLIERGKNM CGLAACASYP IPLV

In the sequential format, the complete alignment for one sequence is
shown before proceeding to the next sequence. In the interleaved format,
the alignments for different sequences are next to each other, for
example in the file ``interlaced.phy`` (provided in ``Tests/Phylip/`` in
the Biopython source distribution):

.. code:: text

    3 384
   CYS1_DICDI   -----MKVIL LFVLAVFTVF VSS------- --------RG IPPEEQ---- --------SQ
   ALEU_HORVU   MAHARVLLLA LAVLATAAVA VASSSSFADS NPIRPVTDRA ASTLESAVLG ALGRTRHALR
   CATH_HUMAN   ------MWAT LPLLCAGAWL LGV------- -PVCGAAELS VNSLEK---- --------FH

                FLEFQDKFNK KY-SHEEYLE RFEIFKSNLG KIEELNLIAI NHKADTKFGV NKFADLSSDE
                FARFAVRYGK SYESAAEVRR RFRIFSESLE EVRSTN---- RKGLPYRLGI NRFSDMSWEE
                FKSWMSKHRK TY-STEEYHH RLQTFASNWR KINAHN---- NGNHTFKMAL NQFSDMSFAE

                FKNYYLNNKE AIFTDDLPVA DYLDDEFINS IPTAFDWRTR G-AVTPVKNQ GQCGSCWSFS
                FQATRL-GAA QTCSATLAGN HLMRDA--AA LPETKDWRED G-IVSPVKNQ AHCGSCWTFS
                IKHKYLWSEP QNCSAT--KS NYLRGT--GP YPPSVDWRKK GNFVSPVKNQ GACGSCWTFS

                TTGNVEGQHF ISQNKLVSLS EQNLVDCDHE CMEYEGEEAC DEGCNGGLQP NAYNYIIKNG
                TTGALEAAYT QATGKNISLS EQQLVDCAGG FNNF------ --GCNGGLPS QAFEYIKYNG
                TTGALESAIA IATGKMLSLA EQQLVDCAQD FNNY------ --GCQGGLPS QAFEYILYNK

                GIQTESSYPY TAETGTQCNF NSANIGAKIS NFTMIP-KNE TVMAGYIVST GPLAIAADAV
                GIDTEESYPY KGVNGV-CHY KAENAAVQVL DSVNITLNAE DELKNAVGLV RPVSVAFQVI
                GIMGEDTYPY QGKDGY-CKF QPGKAIGFVK DVANITIYDE EAMVEAVALY NPVSFAFEVT

                E-WQFYIGGV F-DIPCN--P NSLDHGILIV GYSAKNTIFR KNMPYWIVKN SWGADWGEQG
                DGFRQYKSGV YTSDHCGTTP DDVNHAVLAV GYGVENGV-- ---PYWLIKN SWGADWGDNG
                QDFMMYRTGI YSSTSCHKTP DKVNHAVLAV GYGEKNGI-- ---PYWIVKN SWGPQWGMNG

                YIYLRRGKNT CGVSNFVSTS II--
                YFKMEMGKNM CAIATCASYP VVAA
                YFLIERGKNM CGLAACASYP IPLV

The parser in ``Bio.Align`` detects from the file contents if it is in
the sequential or in the interleaved format, and then parses it
appropriately.

.. doctest ../Tests/Phylip lib:numpy

.. code:: pycon

   >>> from Bio import Align
   >>> alignment = Align.read("sequential.phy", "phylip")
   >>> alignment  # doctest: +ELLIPSIS
   <Alignment object (3 rows x 384 columns) at ...>
   >>> alignment2 = Align.read("interlaced.phy", "phylip")
   >>> alignment2  # doctest: +ELLIPSIS
   <Alignment object (3 rows x 384 columns) at ...>
   >>> alignment == alignment2
   True

Here, two alignments are considered to be equal if they have the same
sequence contents and the same alignment coordinates.

.. cont-doctest

.. code:: pycon

   >>> alignment.shape
   (3, 384)
   >>> print(alignment)
   CYS1_DICD         0 -----MKVILLFVLAVFTVFVSS---------------RGIPPEEQ------------SQ
   ALEU_HORV         0 MAHARVLLLALAVLATAAVAVASSSSFADSNPIRPVTDRAASTLESAVLGALGRTRHALR
   CATH_HUMA         0 ------MWATLPLLCAGAWLLGV--------PVCGAAELSVNSLEK------------FH
   <BLANKLINE>
   CYS1_DICD        28 FLEFQDKFNKKY-SHEEYLERFEIFKSNLGKIEELNLIAINHKADTKFGVNKFADLSSDE
   ALEU_HORV        60 FARFAVRYGKSYESAAEVRRRFRIFSESLEEVRSTN----RKGLPYRLGINRFSDMSWEE
   CATH_HUMA        34 FKSWMSKHRKTY-STEEYHHRLQTFASNWRKINAHN----NGNHTFKMALNQFSDMSFAE
   <BLANKLINE>
   CYS1_DICD        87 FKNYYLNNKEAIFTDDLPVADYLDDEFINSIPTAFDWRTRG-AVTPVKNQGQCGSCWSFS
   ALEU_HORV       116 FQATRL-GAAQTCSATLAGNHLMRDA--AALPETKDWREDG-IVSPVKNQAHCGSCWTFS
   CATH_HUMA        89 IKHKYLWSEPQNCSAT--KSNYLRGT--GPYPPSVDWRKKGNFVSPVKNQGACGSCWTFS
   <BLANKLINE>
   CYS1_DICD       146 TTGNVEGQHFISQNKLVSLSEQNLVDCDHECMEYEGEEACDEGCNGGLQPNAYNYIIKNG
   ALEU_HORV       172 TTGALEAAYTQATGKNISLSEQQLVDCAGGFNNF--------GCNGGLPSQAFEYIKYNG
   CATH_HUMA       145 TTGALESAIAIATGKMLSLAEQQLVDCAQDFNNY--------GCQGGLPSQAFEYILYNK
   <BLANKLINE>
   CYS1_DICD       206 GIQTESSYPYTAETGTQCNFNSANIGAKISNFTMIP-KNETVMAGYIVSTGPLAIAADAV
   ALEU_HORV       224 GIDTEESYPYKGVNGV-CHYKAENAAVQVLDSVNITLNAEDELKNAVGLVRPVSVAFQVI
   CATH_HUMA       197 GIMGEDTYPYQGKDGY-CKFQPGKAIGFVKDVANITIYDEEAMVEAVALYNPVSFAFEVT
   <BLANKLINE>
   CYS1_DICD       265 E-WQFYIGGVF-DIPCN--PNSLDHGILIVGYSAKNTIFRKNMPYWIVKNSWGADWGEQG
   ALEU_HORV       283 DGFRQYKSGVYTSDHCGTTPDDVNHAVLAVGYGVENGV-----PYWLIKNSWGADWGDNG
   CATH_HUMA       256 QDFMMYRTGIYSSTSCHKTPDKVNHAVLAVGYGEKNGI-----PYWIVKNSWGPQWGMNG
   <BLANKLINE>
   CYS1_DICD       321 YIYLRRGKNTCGVSNFVSTSII-- 343
   ALEU_HORV       338 YFKMEMGKNMCAIATCASYPVVAA 362
   CATH_HUMA       311 YFLIERGKNMCGLAACASYPIPLV 335
   <BLANKLINE>

When outputting the alignment in PHYLIP format, ``Bio.Align`` writes
each of the aligned sequences on one line:

.. cont-doctest

.. code:: pycon

   >>> print(format(alignment, "phylip"))
   3 384
   CYS1_DICDI-----MKVILLFVLAVFTVFVSS---------------RGIPPEEQ------------SQFLEFQDKFNKKY-SHEEYLERFEIFKSNLGKIEELNLIAINHKADTKFGVNKFADLSSDEFKNYYLNNKEAIFTDDLPVADYLDDEFINSIPTAFDWRTRG-AVTPVKNQGQCGSCWSFSTTGNVEGQHFISQNKLVSLSEQNLVDCDHECMEYEGEEACDEGCNGGLQPNAYNYIIKNGGIQTESSYPYTAETGTQCNFNSANIGAKISNFTMIP-KNETVMAGYIVSTGPLAIAADAVE-WQFYIGGVF-DIPCN--PNSLDHGILIVGYSAKNTIFRKNMPYWIVKNSWGADWGEQGYIYLRRGKNTCGVSNFVSTSII--
   ALEU_HORVUMAHARVLLLALAVLATAAVAVASSSSFADSNPIRPVTDRAASTLESAVLGALGRTRHALRFARFAVRYGKSYESAAEVRRRFRIFSESLEEVRSTN----RKGLPYRLGINRFSDMSWEEFQATRL-GAAQTCSATLAGNHLMRDA--AALPETKDWREDG-IVSPVKNQAHCGSCWTFSTTGALEAAYTQATGKNISLSEQQLVDCAGGFNNF--------GCNGGLPSQAFEYIKYNGGIDTEESYPYKGVNGV-CHYKAENAAVQVLDSVNITLNAEDELKNAVGLVRPVSVAFQVIDGFRQYKSGVYTSDHCGTTPDDVNHAVLAVGYGVENGV-----PYWLIKNSWGADWGDNGYFKMEMGKNMCAIATCASYPVVAA
   CATH_HUMAN------MWATLPLLCAGAWLLGV--------PVCGAAELSVNSLEK------------FHFKSWMSKHRKTY-STEEYHHRLQTFASNWRKINAHN----NGNHTFKMALNQFSDMSFAEIKHKYLWSEPQNCSAT--KSNYLRGT--GPYPPSVDWRKKGNFVSPVKNQGACGSCWTFSTTGALESAIAIATGKMLSLAEQQLVDCAQDFNNY--------GCQGGLPSQAFEYILYNKGIMGEDTYPYQGKDGY-CKFQPGKAIGFVKDVANITIYDEEAMVEAVALYNPVSFAFEVTQDFMMYRTGIYSSTSCHKTPDKVNHAVLAVGYGEKNGI-----PYWIVKNSWGPQWGMNGYFLIERGKNMCGLAACASYPIPLV
   <BLANKLINE>

We can write the alignment in PHYLIP format, parse the result, and
confirm it is the same as the original alignment object:

.. cont-doctest

.. code:: pycon

   >>> from io import StringIO
   >>> stream = StringIO()
   >>> Align.write(alignment, stream, "phylip")
   1
   >>> stream.seek(0)
   0
   >>> alignment3 = Align.read(stream, "phylip")
   >>> alignment == alignment3
   True
   >>> [record.id for record in alignment.sequences]
   ['CYS1_DICDI', 'ALEU_HORVU', 'CATH_HUMAN']
   >>> [record.id for record in alignment3.sequences]
   ['CYS1_DICDI', 'ALEU_HORVU', 'CATH_HUMAN']

.. _`subsec:align_emboss`:

EMBOSS
~~~~~~

EMBOSS (European Molecular Biology Open Software Suite) is a set of
open-source software tools for molecular biology and
bioinformatics [Rice2000]_. It includes software such
as ``needle`` and ``water`` for pairwise sequence alignment. This is an
example of output generated by the ``water`` program for Smith-Waterman
local pairwise sequence alignment (available as ``water.txt`` in the
``Tests/Emboss`` directory of the Biopython distribution):

.. code:: text

   ########################################
   # Program:  water
   # Rundate:  Wed Jan 16 17:23:19 2002
   # Report_file: stdout
   ########################################
   #=======================================
   #
   # Aligned_sequences: 2
   # 1: IXI_234
   # 2: IXI_235
   # Matrix: EBLOSUM62
   # Gap_penalty: 10.0
   # Extend_penalty: 0.5
   #
   # Length: 131
   # Identity:     112/131 (85.5%)
   # Similarity:   112/131 (85.5%)
   # Gaps:          19/131 (14.5%)
   # Score: 591.5
   #
   #
   #=======================================

   IXI_234            1 TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQAT     50
                        |||||||||||||||         ||||||||||||||||||||||||||
   IXI_235            1 TSPASIRPPAGPSSR---------RPSPPGPRRPTGRPCCSAAPRRPQAT     41

   IXI_234           51 GGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAG    100
                        ||||||||||||||||||||||||          ||||||||||||||||
   IXI_235           42 GGWKTCSGTCTTSTSTRHRGRSGW----------RASRKSMRAACSRSAG     81

   IXI_234          101 SRPNRFAPTLMSSCITSTTGPPAWAGDRSHE    131
                        |||||||||||||||||||||||||||||||
   IXI_235           82 SRPNRFAPTLMSSCITSTTGPPAWAGDRSHE    112


   #---------------------------------------
   #---------------------------------------

As this output file contains only one alignment, we can use
``Align.read`` to extract it directly. Here, instead we will use
``Align.parse`` so we can see the metadata of this ``water`` run:

.. doctest ../Tests/Emboss lib:numpy

.. code:: pycon

   >>> from Bio import Align
   >>> alignments = Align.parse("water.txt", "emboss")

The ``metadata`` attribute of ``alignments`` stores the information
shown in the header of the file, including the program used to generate
the output, the date and time the program was run, the output file name,
and the specific alignment file format that was used (assumed to be
``srspair`` by default):

.. cont-doctest

.. code:: pycon

   >>> alignments.metadata
   {'Align_format': 'srspair', 'Program': 'water', 'Rundate': 'Wed Jan 16 17:23:19 2002', 'Report_file': 'stdout'}

To pull out the alignment, we use

.. cont-doctest

.. code:: pycon

   >>> alignment = next(alignments)
   >>> alignment  # doctest: +ELLIPSIS
   <Alignment object (2 rows x 131 columns) at ...>
   >>> alignment.shape
   (2, 131)
   >>> print(alignment)
   IXI_234           0 TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTC
                     0 |||||||||||||||---------||||||||||||||||||||||||||||||||||||
   IXI_235           0 TSPASIRPPAGPSSR---------RPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTC
   <BLANKLINE>
   IXI_234          60 TTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTG
                    60 ||||||||||||||----------||||||||||||||||||||||||||||||||||||
   IXI_235          51 TTSTSTRHRGRSGW----------RASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTG
   <BLANKLINE>
   IXI_234         120 PPAWAGDRSHE 131
                   120 ||||||||||| 131
   IXI_235         101 PPAWAGDRSHE 112
   <BLANKLINE>
   >>> print(alignment.coordinates)
   [[  0  15  24  74  84 131]
    [  0  15  15  65  65 112]]

We can use indices to extract specific parts of the alignment:

.. cont-doctest

.. code:: pycon

   >>> alignment[0]
   'TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE'
   >>> alignment[1]
   'TSPASIRPPAGPSSR---------RPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGW----------RASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE'
   >>> alignment[1, 10:30]
   'GPSSR---------RPSPPG'

The ``annotations`` attribute of the ``alignment`` stores the
information associated with this alignment specifically:

.. cont-doctest

.. code:: pycon

   >>> alignment.annotations
   {'Matrix': 'EBLOSUM62', 'Gap_penalty': 10.0, 'Extend_penalty': 0.5, 'Identity': 112, 'Similarity': 112, 'Gaps': 19, 'Score': 591.5}

The number of gaps, identities, and mismatches can also be obtained by
calling the ``counts`` method on the ``alignment`` object:

.. cont-doctest

.. code:: pycon

   >>> counts = alignment.counts()

where the ``counts`` variable is an ``AlignmentCounts`` object collecting
information on the number of gaps, matches, and mismatches in the alignment
library (see :ref:`paragraph:alignment_counts`)):

.. cont-doctest

.. code:: pycon

   >>> counts.identities
   112
   >>> counts.mismatches
   0
   >>> counts.gaps
   19

The consensus line shown between the two sequences is stored in the
``column_annotations`` attribute:

.. cont-doctest

.. code:: pycon

   >>> alignment.column_annotations
   {'emboss_consensus': '|||||||||||||||         ||||||||||||||||||||||||||||||||||||||||||||||||||          |||||||||||||||||||||||||||||||||||||||||||||||'}

Use the ``format`` function (or the ``format`` method) to print the
alignment in other formats, for example in the PHYLIP format (see
section :ref:`subsec:align_phylip`):

.. cont-doctest

.. code:: pycon

   >>> print(format(alignment, "phylip"))
   2 131
   IXI_234   TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE
   IXI_235   TSPASIRPPAGPSSR---------RPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGW----------RASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE
   <BLANKLINE>

We can use ``alignment.sequences`` to get the individual sequences.
However, as this is a pairwise alignment, we can also use
``alignment.target`` and ``alignment.query`` to get the target and query
sequences:

.. cont-doctest

.. code:: pycon

   >>> alignment.target
   SeqRecord(seq=Seq('TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQATGGWK...SHE'), id='IXI_234', name='<unknown name>', description='<unknown description>', dbxrefs=[])
   >>> alignment.query
   SeqRecord(seq=Seq('TSPASIRPPAGPSSRRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTS...SHE'), id='IXI_235', name='<unknown name>', description='<unknown description>', dbxrefs=[])

Currently, Biopython does not support writing sequence alignments in the
output formats defined by EMBOSS.

.. _`subsec:align_msf`:

GCG Multiple Sequence Format (MSF)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Multiple Sequence Format (MSF) was created to store multiple
sequence alignments generated by the GCG (Genetics Computer Group) set
of programs. The file ``W_prot.msf`` in the ``Tests/msf`` directory of
the Biopython distribution is an example of a sequence alignment file in
the MSF format This file shows an alignment of 11 protein sequences:

.. code:: text

   !!AA_MULTIPLE_ALIGNMENT

      MSF: 99  Type: P  Oct 18, 2017  11:35  Check: 0 ..

    Name: W*01:01:01:01    Len:    99  Check: 7236  Weight:  1.00
    Name: W*01:01:01:02    Len:    99  Check: 7236  Weight:  1.00
    Name: W*01:01:01:03    Len:    99  Check: 7236  Weight:  1.00
    Name: W*01:01:01:04    Len:    99  Check: 7236  Weight:  1.00
    Name: W*01:01:01:05    Len:    99  Check: 7236  Weight:  1.00
    Name: W*01:01:01:06    Len:    99  Check: 7236  Weight:  1.00
    Name: W*02:01          Len:    93  Check: 9483  Weight:  1.00
    Name: W*03:01:01:01    Len:    93  Check: 9974  Weight:  1.00
    Name: W*03:01:01:02    Len:    93  Check: 9974  Weight:  1.00
    Name: W*04:01          Len:    93  Check: 9169  Weight:  1.00
    Name: W*05:01          Len:    99  Check: 7331  Weight:  1.00
   //

     W*01:01:01:01  GLTPFNGYTA ATWTRTAVSS VGMNIPYHGA SYLVRNQELR SWTAADKAAQ
     W*01:01:01:02  GLTPFNGYTA ATWTRTAVSS VGMNIPYHGA SYLVRNQELR SWTAADKAAQ
     W*01:01:01:03  GLTPFNGYTA ATWTRTAVSS VGMNIPYHGA SYLVRNQELR SWTAADKAAQ
     W*01:01:01:04  GLTPFNGYTA ATWTRTAVSS VGMNIPYHGA SYLVRNQELR SWTAADKAAQ
     W*01:01:01:05  GLTPFNGYTA ATWTRTAVSS VGMNIPYHGA SYLVRNQELR SWTAADKAAQ
     W*01:01:01:06  GLTPFNGYTA ATWTRTAVSS VGMNIPYHGA SYLVRNQELR SWTAADKAAQ
           W*02:01  GLTPSNGYTA ATWTRTAASS VGMNIPYDGA SYLVRNQELR SWTAADKAAQ
     W*03:01:01:01  GLTPSSGYTA ATWTRTAVSS VGMNIPYHGA SYLVRNQELR SWTAADKAAQ
     W*03:01:01:02  GLTPSSGYTA ATWTRTAVSS VGMNIPYHGA SYLVRNQELR SWTAADKAAQ
           W*04:01  GLTPSNGYTA ATWTRTAASS VGMNIPYDGA SYLVRNQELR SWTAADKAAQ
           W*05:01  GLTPSSGYTA ATWTRTAVSS VGMNIPYHGA SYLVRNQELR SWTAADKAAQ

     W*01:01:01:01  MPWRRNRQSC SKPTCREGGR SGSAKSLRMG RRGCSAQNPK DSHDPPPHL
     W*01:01:01:02  MPWRRNRQSC SKPTCREGGR SGSAKSLRMG RRGCSAQNPK DSHDPPPHL
     W*01:01:01:03  MPWRRNRQSC SKPTCREGGR SGSAKSLRMG RRGCSAQNPK DSHDPPPHL
     W*01:01:01:04  MPWRRNRQSC SKPTCREGGR SGSAKSLRMG RRGCSAQNPK DSHDPPPHL
     W*01:01:01:05  MPWRRNRQSC SKPTCREGGR SGSAKSLRMG RRGCSAQNPK DSHDPPPHL
     W*01:01:01:06  MPWRRNRQSC SKPTCREGGR SGSAKSLRMG RRGCSAQNPK DSHDPPPHL
           W*02:01  MPWRRNMQSC SKPTCREGGR SGSAKSLRMG RRRCTAQNPK RLT
     W*03:01:01:01  MPWRRNRQSC SKPTCREGGR SGSAKSLRMG RRGCSAQNPK RLT
     W*03:01:01:02  MPWRRNRQSC SKPTCREGGR SGSAKSLRMG RRGCSAQNPK RLT
           W*04:01  MPWRRNMQSC SKPTCREGGR SGSAKSLRMG RRGCSAQNPK RLT
           W*05:01  MPWRRNRQSC SKPTCREGGR SGSAKSLRMG RRGCSAQNPK DSHDPPPHL

To parse this file with Biopython, use

.. doctest ../Tests/msf lib:numpy

.. code:: pycon

   >>> from Bio import Align
   >>> alignment = Align.read("W_prot.msf", "msf")

The parser skips all lines up to and including the line starting with
"``MSF:``". The following lines (until the "``//``" demarcation) are
read by the parser to verify the length of each sequence. The alignment
section (after the "``//``" demarcation) is read by the parser and
stored as an ``Alignment`` object:

.. cont-doctest

.. code:: pycon

   >>> alignment  # doctest: +ELLIPSIS
   <Alignment object (11 rows x 99 columns) at ...>
   >>> print(alignment)
   W*01:01:0         0 GLTPFNGYTAATWTRTAVSSVGMNIPYHGASYLVRNQELRSWTAADKAAQMPWRRNRQSC
   W*01:01:0         0 GLTPFNGYTAATWTRTAVSSVGMNIPYHGASYLVRNQELRSWTAADKAAQMPWRRNRQSC
   W*01:01:0         0 GLTPFNGYTAATWTRTAVSSVGMNIPYHGASYLVRNQELRSWTAADKAAQMPWRRNRQSC
   W*01:01:0         0 GLTPFNGYTAATWTRTAVSSVGMNIPYHGASYLVRNQELRSWTAADKAAQMPWRRNRQSC
   W*01:01:0         0 GLTPFNGYTAATWTRTAVSSVGMNIPYHGASYLVRNQELRSWTAADKAAQMPWRRNRQSC
   W*01:01:0         0 GLTPFNGYTAATWTRTAVSSVGMNIPYHGASYLVRNQELRSWTAADKAAQMPWRRNRQSC
   W*02:01           0 GLTPSNGYTAATWTRTAASSVGMNIPYDGASYLVRNQELRSWTAADKAAQMPWRRNMQSC
   W*03:01:0         0 GLTPSSGYTAATWTRTAVSSVGMNIPYHGASYLVRNQELRSWTAADKAAQMPWRRNRQSC
   W*03:01:0         0 GLTPSSGYTAATWTRTAVSSVGMNIPYHGASYLVRNQELRSWTAADKAAQMPWRRNRQSC
   W*04:01           0 GLTPSNGYTAATWTRTAASSVGMNIPYDGASYLVRNQELRSWTAADKAAQMPWRRNMQSC
   W*05:01           0 GLTPSSGYTAATWTRTAVSSVGMNIPYHGASYLVRNQELRSWTAADKAAQMPWRRNRQSC
   <BLANKLINE>
   W*01:01:0        60 SKPTCREGGRSGSAKSLRMGRRGCSAQNPKDSHDPPPHL 99
   W*01:01:0        60 SKPTCREGGRSGSAKSLRMGRRGCSAQNPKDSHDPPPHL 99
   W*01:01:0        60 SKPTCREGGRSGSAKSLRMGRRGCSAQNPKDSHDPPPHL 99
   W*01:01:0        60 SKPTCREGGRSGSAKSLRMGRRGCSAQNPKDSHDPPPHL 99
   W*01:01:0        60 SKPTCREGGRSGSAKSLRMGRRGCSAQNPKDSHDPPPHL 99
   W*01:01:0        60 SKPTCREGGRSGSAKSLRMGRRGCSAQNPKDSHDPPPHL 99
   W*02:01          60 SKPTCREGGRSGSAKSLRMGRRRCTAQNPKRLT------ 93
   W*03:01:0        60 SKPTCREGGRSGSAKSLRMGRRGCSAQNPKRLT------ 93
   W*03:01:0        60 SKPTCREGGRSGSAKSLRMGRRGCSAQNPKRLT------ 93
   W*04:01          60 SKPTCREGGRSGSAKSLRMGRRGCSAQNPKRLT------ 93
   W*05:01          60 SKPTCREGGRSGSAKSLRMGRRGCSAQNPKDSHDPPPHL 99
   <BLANKLINE>

The sequences and their names are stored in the ``alignment.sequences``
attribute:

.. cont-doctest

.. code:: pycon

   >>> len(alignment.sequences)
   11
   >>> alignment.sequences[0].id
   'W*01:01:01:01'
   >>> alignment.sequences[0].seq
   Seq('GLTPFNGYTAATWTRTAVSSVGMNIPYHGASYLVRNQELRSWTAADKAAQMPWR...PHL')

The alignment coordinates are stored in the ``alignment.coordinates``
attribute:

.. cont-doctest

.. code:: pycon

   >>> print(alignment.coordinates)
   [[ 0 93 99]
    [ 0 93 99]
    [ 0 93 99]
    [ 0 93 99]
    [ 0 93 99]
    [ 0 93 99]
    [ 0 93 93]
    [ 0 93 93]
    [ 0 93 93]
    [ 0 93 93]
    [ 0 93 99]]

Currently, Biopython does not support writing sequence alignments in the
MSF format.

.. _`subsec:align_exonerate`:

Exonerate
~~~~~~~~~

Exonerate is a generic program for pairwise sequence
alignments [Slater2005]_. The sequence alignments found
by Exonerate can be output in a human-readable form, in the "cigar"
(Compact Idiosyncratic Gapped Alignment Report) format, or in the
"vulgar" (Verbose Useful Labelled Gapped Alignment Report) format. The
user can request to include one or more of these formats in the output.
The parser in ``Bio.Align`` can only parse alignments in the cigar or
vulgar formats, and will not parse output that includes alignments in
human-readable format.

The file ``exn_22_m_cdna2genome_vulgar.exn`` in the Biopython test suite
is an example of an Exonerate output file showing the alignments in
vulgar format:

.. code:: text

   Command line: [exonerate -m cdna2genome ../scer_cad1.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showcigar no --showvulgar yes]
   Hostname: [blackbriar]
   vulgar: gi|296143771|ref|NM_001180731.1| 0 1230 + gi|330443520|ref|NC_001136.10| 1319275 1318045 - 6146 M 1 1 C 3 3 M 1226 1226
   vulgar: gi|296143771|ref|NM_001180731.1| 1230 0 - gi|330443520|ref|NC_001136.10| 1318045 1319275 + 6146 M 129 129 C 3 3 M 1098 1098
   vulgar: gi|296143771|ref|NM_001180731.1| 0 516 + gi|330443688|ref|NC_001145.3| 85010 667216 + 518 M 11 11 G 1 0 M 15 15 G 2 0 M 4 4 G 1 0 M 1 1 G 1 0 M 8 8 G 4 0 M 17 17 5 0 2 I 0 168904 3 0 2 M 4 4 G 0 1 M 8 8 G 2 0 M 3 3 G 1 0 M 33 33 G 0 2 M 7 7 G 0 1 M 102 102 5 0 2 I 0 96820 3 0 2 M 14 14 G 0 2 M 10 10 G 2 0 M 5 5 G 0 2 M 10 10 G 2 0 M 4 4 G 0 1 M 20 20 G 1 0 M 15 15 G 0 1 M 5 5 G 3 0 M 4 4 5 0 2 I 0 122114 3 0 2 M 20 20 G 0 5 M 6 6 5 0 2 I 0 193835 3 0 2 M 12 12 G 0 2 M 5 5 G 1 0 M 7 7 G 0 2 M 1 1 G 0 1 M 12 12 C 75 75 M 6 6 G 1 0 M 4 4 G 0 1 M 2 2 G 0 1 M 3 3 G 0 1 M 41 41
   -- completed exonerate analysis

This file includes three alignments. To parse this file, use

.. doctest ../Tests/Exonerate lib:numpy

.. code:: pycon

   >>> from Bio import Align
   >>> alignments = Align.parse("exn_22_m_cdna2genome_vulgar.exn", "exonerate")

The dictionary ``alignments.metadata`` stores general information about
these alignments, shown at the top of the output file:

.. cont-doctest

.. code:: pycon

   >>> alignments.metadata  # doctest: +NORMALIZE_WHITESPACE
   {'Program': 'exonerate',
    'Command line': 'exonerate -m cdna2genome ../scer_cad1.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showcigar no --showvulgar yes',
    'Hostname': 'blackbriar'}

Now we can iterate over the alignments. The first alignment, with
alignment score 6146.0, has no gaps:

.. cont-doctest

.. code:: pycon

   >>> alignment = next(alignments)
   >>> alignment.score
   6146.0
   >>> print(alignment.coordinates)
   [[1319275 1319274 1319271 1318045]
    [      0       1       4    1230]]
   >>> print(alignment)  # doctest: +ELLIPSIS
   gi|330443   1319275 ????????????????????????????????????????????????????????????
                     0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
   gi|296143         0 ????????????????????????????????????????????????????????????
   ...
   gi|330443   1318075 ?????????????????????????????? 1318045
                  1200 ||||||||||||||||||||||||||||||    1230
   gi|296143      1200 ??????????????????????????????    1230
   <BLANKLINE>

Note that the target (the first sequence) in the printed alignment is on
the reverse strand while the query (the second sequence) is on the
forward strand, with the target coordinate decreasing and the query
coordinate increasing. Printing this alignment in ``exonerate`` format
using Python’s built-in ``format`` function writes a vulgar line:

.. cont-doctest

.. code:: pycon

   >>> print(format(alignment, "exonerate"))
   vulgar: gi|296143771|ref|NM_001180731.1| 0 1230 + gi|330443520|ref|NC_001136.10| 1319275 1318045 - 6146 M 1 1 C 3 3 M 1226 1226
   <BLANKLINE>

Using the ``format`` method allows us to request either a vulgar line
(default) or a cigar line:

.. cont-doctest

.. code:: pycon

   >>> print(alignment.format("exonerate", "vulgar"))
   vulgar: gi|296143771|ref|NM_001180731.1| 0 1230 + gi|330443520|ref|NC_001136.10| 1319275 1318045 - 6146 M 1 1 C 3 3 M 1226 1226
   <BLANKLINE>
   >>> print(alignment.format("exonerate", "cigar"))
   cigar: gi|296143771|ref|NM_001180731.1| 0 1230 + gi|330443520|ref|NC_001136.10| 1319275 1318045 - 6146 M 1 M 3 M 1226
   <BLANKLINE>

The vulgar line contains information about the alignment (in the section
``M 1 1 C 3 3 M 1226``) that is missing from the cigar line
``M 1 M 3 M 1226``. The vulgar line specifies that the alignment starts
with a single aligned nucleotides, followed by three aligned nucleotides
that form a codon (``C``), followed by 1226 aligned nucleotides. In the
cigar line, we see a single aligned nucleotide, followed by three
aligned nucleotides, followed by 1226 aligned nucleotides; it does not
specify that the three aligned nucleotides form a codon. This
information from the vulgar line is stored in the ``operations``
attribute:

.. cont-doctest

.. code:: pycon

   >>> alignment.operations
   bytearray(b'MCM')

See the Exonerate documentation for the definition of other operation
codes.

Similarly, the ``"vulgar"`` or ``"cigar"`` argument can be used when
calling ``Bio.Align.write`` to write a file with vulgar or cigar
alignment lines.

We can print the alignment in BED and PSL format:

.. cont-doctest

.. code:: pycon

   >>> print(format(alignment, "bed"))  # doctest: +NORMALIZE_WHITESPACE
   gi|330443520|ref|NC_001136.10|  1318045 1319275 gi|296143771|ref|NM_001180731.1| 6146   -   1318045 1319275 0   3   1226,3,1,   0,1226,1229,
   <BLANKLINE>
   >>> print(format(alignment, "psl"))  # doctest: +NORMALIZE_WHITESPACE
   1230    0   0   0   0   0   0   0   -   gi|296143771|ref|NM_001180731.1|    1230    0   1230    gi|330443520|ref|NC_001136.10|  1319275 1318045 1319275 3   1226,3,1,   0,1226,1229,    1318045,1319271,1319274,
   <BLANKLINE>

The SAM format parser defines its own (optional) ``operations``
attribute (section :ref:`subsec:align_sam`), which is not quite
consistent with the ``operations`` attribute defined in the Exonerate
format parser. As the ``operations`` attribute is optional, we delete it
before printing the alignment in SAM format:

.. cont-doctest

.. code:: pycon

   >>> del alignment.operations
   >>> print(format(alignment, "sam"))  # doctest: +NORMALIZE_WHITESPACE
   gi|296143771|ref|NM_001180731.1|    16  gi|330443520|ref|NC_001136.10|  1318046 255 1226M3M1M   *   0   0   *   *   AS:i:6146
   <BLANKLINE>

The third alignment contains four long gaps:

.. cont-doctest

.. code:: pycon

   >>> alignment = next(alignments)  # second alignment
   >>> alignment = next(alignments)  # third alignment
   >>> print(alignment)  # doctest: +ELLIPSIS
   gi|330443     85010 ???????????-???????????????--????-?-????????----????????????
                     0 |||||||||||-|||||||||||||||--||||-|-||||||||----||||||||||||
   gi|296143         0 ????????????????????????????????????????????????????????????
   <BLANKLINE>
   gi|330443     85061 ????????????????????????????????????????????????????????????
                    60 |||||-------------------------------------------------------
   gi|296143        60 ?????-------------------------------------------------------
   ...
   gi|330443    666990 ????????????????????????????????????????????????????????????
                582000 --------------------------------------------------||||||||||
   gi|296143       346 --------------------------------------------------??????????
   <BLANKLINE>
   gi|330443    667050 ?????????-??????????????????????????????????????????????????
                582060 ||--|||||-|||||||--|-|||||||||||||||||||||||||||||||||||||||
   gi|296143       356 ??--?????????????--?-???????????????????????????????????????
   <BLANKLINE>
   gi|330443    667109 ??????????????????????????????????????????????????????-?????
                582120 ||||||||||||||||||||||||||||||||||||||||||||||||||||||-||||-
   gi|296143       411 ???????????????????????????????????????????????????????????-
   <BLANKLINE>
   gi|330443    667168 ???????????????????????????????????????????????? 667216
                582180 ||-|||-||||||||||||||||||||||||||||||||||||||||| 582228
   gi|296143       470 ??-???-?????????????????????????????????????????    516
   <BLANKLINE>
   >>> print(format(alignment, "exonerate"))  # doctest: +NORMALIZE_WHITESPACE
   vulgar: gi|296143771|ref|NM_001180731.1| 0 516 + gi|330443688|ref|NC_001145.3|
   85010 667216 + 518 M 11 11 G 1 0 M 15 15 G 2 0 M 4 4 G 1 0 M 1 1 G 1 0 M 8 8
    G 4 0 M 17 17 5 0 2 I 0 168904 3 0 2 M 4 4 G 0 1 M 8 8 G 2 0 M 3 3 G 1 0
    M 33 33 G 0 2 M 7 7 G 0 1 M 102 102 5 0 2 I 0 96820 3 0 2 M 14 14 G 0 2 M 10 10
    G 2 0 M 5 5 G 0 2 M 10 10 G 2 0 M 4 4 G 0 1 M 20 20 G 1 0 M 15 15 G 0 1 M 5 5
    G 3 0 M 4 4 5 0 2 I 0 122114 3 0 2 M 20 20 G 0 5 M 6 6 5 0 2 I 0 193835 3 0 2
    M 12 12 G 0 2 M 5 5 G 1 0 M 7 7 G 0 2 M 1 1 G 0 1 M 12 12 C 75 75 M 6 6 G 1 0
    M 4 4 G 0 1 M 2 2 G 0 1 M 3 3 G 0 1 M 41 41
   <BLANKLINE>

.. _`subsec:align_nexus`:

NEXUS
~~~~~

The NEXUS file format [Maddison1997]_ is used by
several programs to store phylogenetic information. This is an example
of a file in the NEXUS format (available as ``codonposset.nex`` in the
``Tests/Nexus`` subdirectory in the Biopython distribution):

.. code:: text

   #NEXUS
   [MacClade 4.05 registered to Computational Biologist, University]


   BEGIN DATA;
          DIMENSIONS  NTAX=2 NCHAR=22;
          FORMAT DATATYPE=DNA  MISSING=? GAP=- ;
   MATRIX
   [                           10        20 ]
   [                           .         .  ]

   Aegotheles         AAAAAGGCATTGTGGTGGGAAT   [22]
   Aerodramus         ?????????TTGTGGTGGGAAT   [13]
   ;
   END;


   BEGIN CODONS;
          CODONPOSSET * CodonPositions =
                  N: 1-10,
                  1: 11-22\3,
                  2: 12-22\3,
                  3: 13-22\3;
          CODESET  * UNTITLED = Universal: all ;
   END;

In general, files in the NEXUS format can be much more complex.
``Bio.Align`` relies heavily on NEXUS parser in ``Bio.Nexus`` (see
Chapter :ref:`chapter:phylo`) to extract ``Alignment``
objects from NEXUS files.

To read the alignment in this NEXUS file, use

.. doctest ../Tests/Nexus lib:numpy

.. code:: pycon

   >>> from Bio import Align
   >>> alignment = Align.read("codonposset.nex", "nexus")
   >>> print(alignment)
   Aegothele         0 AAAAAGGCATTGTGGTGGGAAT 22
                     0 .........||||||||||||| 22
   Aerodramu         0 ?????????TTGTGGTGGGAAT 22
   <BLANKLINE>
   >>> alignment.shape
   (2, 22)

The sequences are stored under the ``sequences`` attribute:

.. cont-doctest

.. code:: pycon

   >>> alignment.sequences[0].id
   'Aegotheles'
   >>> alignment.sequences[0].seq
   Seq('AAAAAGGCATTGTGGTGGGAAT')
   >>> alignment.sequences[0].annotations
   {'molecule_type': 'DNA'}
   >>> alignment.sequences[1].id
   'Aerodramus'
   >>> alignment.sequences[1].seq
   Seq('?????????TTGTGGTGGGAAT')
   >>> alignment.sequences[1].annotations
   {'molecule_type': 'DNA'}

To print this alignment in the NEXUS format, use

.. cont-doctest

.. code:: pycon

   >>> print(format(alignment, "nexus"))
   #NEXUS
   begin data;
   dimensions ntax=2 nchar=22;
   format datatype=dna missing=? gap=-;
   matrix
   Aegotheles AAAAAGGCATTGTGGTGGGAAT
   Aerodramus ?????????TTGTGGTGGGAAT
   ;
   end;
   <BLANKLINE>

Similarly, you can use
``Align.write(alignment, "myfilename.nex", "nexus")`` to write the
alignment in the NEXUS format to the file ``myfilename.nex``.

.. _`subsec:align_tabular`:

Tabular output from BLAST or FASTA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Alignment output in tabular output is generated by the FASTA
aligner [Pearson1988]_ run with the ``-m 8CB`` or
``-m 8CC`` argument, or by BLAST [Altschul1990]_ run
with the ``-outfmt 7`` argument.

The file ``nucleotide_m8CC.txt`` in the ``Tests/Fasta`` subdirectory of
the Biopython source distribution is an example of an output file
generated by FASTA with the ``-m 8CC`` argument:

.. code:: text

   # fasta36 -m 8CC seq/mgstm1.nt seq/gst.nlib
   # FASTA 36.3.8h May, 2020
   # Query: pGT875  - 657 nt
   # Database: seq/gst.nlib
   # Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score, aln_code
   # 12 hits found
   pGT875  pGT875  100.00  657 0   0   1   657 38  694 4.6e-191    655.6   657M
   pGT875  RABGLTR 79.10   646 135 0   1   646 34  679 1.6e-116    408.0   646M
   pGT875  BTGST   59.56   413 167 21  176 594 228 655 1.9e-07 45.7    149M1D7M1I17M3D60M5I6M1I13M2I13M4I30M2I6M2D112M
   pGT875  RABGSTB 66.93   127 42  8   159 289 157 287 3.2e-07 45.0    15M2I17M2D11M1I58M1I11M1D7M1D8M
   pGT875  OCDHPR  91.30   23  2   1   266 289 2303    2325    0.012   29.7    17M1D6M
   ...
   # FASTA processed 1 queries

To parse this file, use

.. doctest ../Tests/Fasta lib:numpy

.. code:: pycon

   >>> from Bio import Align
   >>> alignments = Align.parse("nucleotide_m8CC.txt", "tabular")

Information shown in the file header is stored in the ``metadata``
attribute of ``alignments``:

.. cont-doctest

.. code:: pycon

   >>> alignments.metadata  # doctest: +NORMALIZE_WHITESPACE
   {'Command line': 'fasta36 -m 8CC seq/mgstm1.nt seq/gst.nlib',
    'Program': 'FASTA',
    'Version': '36.3.8h May, 2020',
    'Database': 'seq/gst.nlib'}

Extract a specific alignment by iterating over the ``alignments``. As an
example, let’s go to the fourth alignment:

.. cont-doctest

.. code:: pycon

   >>> alignment = next(alignments)
   >>> alignment = next(alignments)
   >>> alignment = next(alignments)
   >>> alignment = next(alignments)
   >>> print(alignment)
   RABGSTB         156 ??????????????????????????????????--????????????????????????
                     0 |||||||||||||||--|||||||||||||||||--|||||||||||-||||||||||||
   pGT875          158 ???????????????--??????????????????????????????-????????????
   <BLANKLINE>
   RABGSTB         214 ??????????????????????????????????????????????????????????-?
                    60 ||||||||||||||||||||||||||||||||||||||||||||||-|||||||||||-|
   pGT875          215 ??????????????????????????????????????????????-?????????????
   <BLANKLINE>
   RABGSTB         273 ??????-???????? 287
                   120 ||||||-|||||||| 135
   pGT875          274 ??????????????? 289
   <BLANKLINE>
   >>> print(alignment.coordinates)
   [[156 171 173 190 190 201 202 260 261 272 272 279 279 287]
    [158 173 173 190 192 203 203 261 261 272 273 280 281 289]]
   >>> alignment.aligned
   array([[[156, 171],
           [173, 190],
           [190, 201],
           [202, 260],
           [261, 272],
           [272, 279],
           [279, 287]],
   <BLANKLINE>
          [[158, 173],
           [173, 190],
           [192, 203],
           [203, 261],
           [261, 272],
           [273, 280],
           [281, 289]]])

The sequence information of the target and query sequences is stored in
the ``target`` and ``query`` attributes (as well as under
``alignment.sequences``):

.. cont-doctest

.. code:: pycon

   >>> alignment.target
   SeqRecord(seq=Seq(None, length=287), id='RABGSTB', name='<unknown name>', description='<unknown description>', dbxrefs=[])
   >>> alignment.query
   SeqRecord(seq=Seq(None, length=657), id='pGT875', name='<unknown name>', description='<unknown description>', dbxrefs=[])

Information of the alignment is stored under the ``annotations``
attribute of the ``alignment``:

.. cont-doctest

.. code:: pycon

   >>> alignment.annotations  # doctest: +NORMALIZE_WHITESPACE
   {'% identity': 66.93,
    'mismatches': 42,
    'gap opens': 8,
    'evalue': 3.2e-07,
    'bit score': 45.0}

BLAST in particular offers many options to customize tabular output by
including or excluding specific columns; see the BLAST documentation for
details. This information is stored in the dictionaries
``alignment.annotations``, ``alignment.target.annotations``, or
``alignment.query.annotations``, as appropriate.

.. _`subsec:align_hhr`:

HH-suite output files
~~~~~~~~~~~~~~~~~~~~~

Alignment files in the ``hhr`` format are generated by ``hhsearch`` or
``hhblits`` in HH-suite [Steinegger2019]_. As an
example, see the file ``2uvo_hhblits.hhr`` in Biopython’s test suite:

.. code:: text

   Query         2UVO:A|PDBID|CHAIN|SEQUENCE
   Match_columns 171
   No_of_seqs    1560 out of 4005
   Neff          8.3
   Searched_HMMs 34
   Date          Fri Feb 15 16:34:13 2019
   Command       hhblits -i 2uvoAh.fasta -d /pdb70

    No Hit                             Prob E-value P-value  Score    SS Cols Query HMM  Template HMM
     1 2uvo_A Agglutinin isolectin 1; 100.0 3.7E-34 4.8E-38  210.3   0.0  171    1-171     1-171 (171)
     2 2wga   ; lectin (agglutinin);   99.9 1.1E-30 1.4E-34  190.4   0.0  162    2-169     2-163 (164)
     3 1ulk_A Lectin-C; chitin-bindin  99.8 5.2E-24 6.6E-28  148.2   0.0  120    1-124     2-121 (126)
   ...
    31 4z8i_A BBTPGRP3, peptidoglycan  79.6    0.12 1.5E-05   36.1   0.0   37    1-37      9-54  (236)
    32 1wga   ; lectin (agglutinin);   40.4     2.6 0.00029   25.9   0.0  106   54-163    11-116 (164)

   No 1
   >2uvo_A Agglutinin isolectin 1; carbohydrate-binding protein, hevein domain, chitin-binding, GERM agglutinin, chitin-binding protein; HET: NDG NAG GOL; 1.40A {Triticum aestivum} PDB: 1wgc_A* 2cwg_A* 2x3t_A* 4aml_A* 7wga_A 9wga_A 2wgc_A 1wgt_A 1k7t_A* 1k7v_A* 1k7u_A 2x52_A* 1t0w_A*
   Probab=99.95  E-value=3.7e-34  Score=210.31  Aligned_cols=171  Identities=100%  Similarity=2.050  Sum_probs=166.9

   Q 2UVO:A|PDBID|C    1 ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQG   80 (171)
   Q Consensus         1 ~~cg~~~~~~~c~~~~CCs~~g~CG~~~~~c~~~c~~~~c~~~~~Cg~~~~~~~c~~~~CCs~~g~CG~~~~~c~~~c~~   80 (171)
                         ||||++.++..||++.|||+|+|||.+.+||+++||.+.|++..+|+++++.++|....|||.++||+.+.+||+.+||.
   T Consensus         1 ~~cg~~~~~~~c~~~~CCS~~g~Cg~~~~~Cg~gC~~~~c~~~~~cg~~~~~~~c~~~~CCs~~g~Cg~~~~~c~~~c~~   80 (171)
   T 2uvo_A            1 ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQG   80 (171)
   T ss_dssp             CBCBGGGTTBBCGGGCEECTTSBEEBSHHHHSTTCCBSSCSSCCBCBGGGTTBCCSTTCEECTTSBEEBSHHHHSTTCCB
   T ss_pred             CCCCCCCCCcCCCCCCeeCCCCeECCCcccccCCccccccccccccCcccCCcccCCccccCCCceeCCCccccCCCccc
   Confidence            79999999999999999999999999999999999999999999999999999999999999999999999999999999


   Q 2UVO:A|PDBID|C   81 GPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYC  160 (171)
   Q Consensus        81 ~~~~~~~~Cg~~~~~~~c~~~~CCS~~G~CG~~~~~C~~~Cq~~~c~~~~~Cg~~~~~~~c~~~~CCS~~G~CG~~~~~C  160 (171)
                         +++++|+.|+...+++.||++.|||.|||||...+||+.+||+++|++|.+|++.+++++|..+.|||+++-||+...||
   T Consensus        81 ~~~~~~~~cg~~~~~~~c~~~~CCs~~g~CG~~~~~C~~gCq~~~c~~~~~cg~~~~~~~c~~~~ccs~~g~Cg~~~~~C  160 (171)
   T 2uvo_A           81 GPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYC  160 (171)
   T ss_dssp             SSCSSCCBCBGGGTTBCCGGGCEECTTSBEEBSHHHHSTTCCBSSCSSCCCCBTTTTTBCCSTTCEECTTSCEEBSHHHH
   T ss_pred             ccccccccccccccCCCCCCCcccCCCCccCCCcccccCCCcCCccccccccccccccccCCCCCCcCCCCEecCchhhc
   Confidence            99999999999988999999999999999999999999999999999999999999999999999999999999999999


   Q 2UVO:A|PDBID|C  161 GAGCQSGGCDG  171 (171)
   Q Consensus       161 ~~gCq~~~c~~  171 (171)
                         +++||++.|||
   T Consensus       161 ~~~cq~~~~~~  171 (171)
   T 2uvo_A          161 GAGCQSGGCDG  171 (171)
   T ss_dssp             STTCCBSSCC-
   T ss_pred             ccccccCCCCC
   Confidence            99999999986


   No 2
   ...


   No 32
   >1wga   ; lectin (agglutinin); NMR {}
   Probab=40.43  E-value=2.6  Score=25.90  Aligned_cols=106  Identities=20%  Similarity=0.652  Sum_probs=54.7

   Q 2UVO:A|PDBID|C   54 TCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCG  133 (171)
   Q Consensus        54 ~c~~~~CCs~~g~CG~~~~~c~~~c~~~~~~~~~~Cg~~~~~~~c~~~~CCS~~G~CG~~~~~C~~~Cq~~~c~~~~~Cg  133 (171)
                         .|....||.....|......|...|....|.....|...  ...|....||.....|......|...|....+.....|.
   T Consensus        11 ~c~~~~cc~~~~~c~~~~~~c~~~c~~~~c~~~~~c~~~--~~~c~~~~cc~~~~~c~~~~~~c~~~c~~~~c~~~~~c~   88 (164)
   T 1wga             11 XCXXXXCCXXXXXCXXXXXXCXXXCXXXXCXXXXXCXXX--XXXCXXXXCCXXXXXCXXXXXXCXXXCXXXXCXXXXXCX   88 (164)
   T ss_pred             ccccccccccccccccccccccccccccccccccccccc--ccccccccccccccccccccccccccccccccccccccc
   Confidence            344556666666666666566555543333223333321  234666677777777777766666655544332223333


   Q 2UVO:A|PDBID|C  134 KDAGGRVCTNNYCCSKWGSCGIGPGYCGAG  163 (171)
   Q Consensus       134 ~~~~~~~c~~~~CCS~~G~CG~~~~~C~~g  163 (171)
                         ..  ...|....||.....|......|...
   T Consensus        89 ~~--~~~c~~~~cc~~~~~c~~~~~~c~~~  116 (164)
   T 1wga             89 XX--XXXCXXXXCCXXXXXCXXXXXXCXXX  116 (164)
   T ss_pred             cc--cccccccccccccccccccccccccc
   Confidence            22  23344455555555555555544433


   Done!

The file contains three sections:

-  A header with general information about the alignments;

-  A summary with one line for each of the alignments obtained;

-  The alignments shown consecutively in detail.

To parse this file, use

.. doctest ../Tests/HHsuite lib:numpy

.. code:: pycon

   >>> from Bio import Align
   >>> alignments = Align.parse("2uvo_hhblits.hhr", "hhr")

Most of the header information is stored in the ``metadata`` attribute
of ``alignments``:

.. cont-doctest

.. code:: pycon

   >>> alignments.metadata  # doctest: +NORMALIZE_WHITESPACE
   {'Match_columns': 171,
    'No_of_seqs': (1560, 4005),
    'Neff': 8.3,
    'Searched_HMMs': 34,
    'Rundate': 'Fri Feb 15 16:34:13 2019',
    'Command line': 'hhblits -i 2uvoAh.fasta -d /pdb70'}

except the query name, which is stored as an attribute:

.. cont-doctest

.. code:: pycon

   >>> alignments.query_name
   '2UVO:A|PDBID|CHAIN|SEQUENCE'

as it will reappear in each of the alignments.

Iterate over the alignments:

.. cont-doctest

.. code:: pycon

   >>> for alignment in alignments:
   ...     print(alignment.target.id)  # doctest: +ELLIPSIS
   ...
   2uvo_A
   2wga
   1ulk_A
   ...
   4z8i_A
   1wga

Let’s look at the first alignment in more detail:

.. cont-doctest

.. code:: pycon

   >>> alignments = iter(alignments)
   >>> alignment = next(alignments)
   >>> alignment  # doctest: +ELLIPSIS
   <Alignment object (2 rows x 171 columns) at ...>
   >>> print(alignment)
   2uvo_A            0 ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQC
                     0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
   2UVO:A|PD         0 ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQC
   <BLANKLINE>
   2uvo_A           60 CSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGG
                    60 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
   2UVO:A|PD        60 CSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGG
   <BLANKLINE>
   2uvo_A          120 CQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQSGGCDG 171
                   120 ||||||||||||||||||||||||||||||||||||||||||||||||||| 171
   2UVO:A|PD       120 CQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQSGGCDG 171
   <BLANKLINE>

The target and query sequences are stored in ``alignment.sequences``. As
these are pairwise alignments, we can also access them through
``alignment.target`` and ``alignment.query``:

.. cont-doctest

.. code:: pycon

   >>> alignment.target is alignment.sequences[0]
   True
   >>> alignment.query is alignment.sequences[1]
   True

The ID of the query is set from the ``alignments.query_name`` (note that
the query ID printed in the alignment in the ``hhr`` file is
abbreviated):

.. cont-doctest

.. code:: pycon

   >>> alignment.query.id
   '2UVO:A|PDBID|CHAIN|SEQUENCE'

The ID of the target is taken from the sequence alignment block (the
line starting with ``T 2uvo_A``):

.. code:: pycon

   >>> alignment.target.id
   '2uvo_A'

The sequence contents of the target and query are filled in from the
information available in this alignment:

.. cont-doctest

.. code:: pycon

   >>> alignment.target.seq
   Seq('ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGAT...CDG')
   >>> alignment.query.seq
   Seq('ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGAT...CDG')

The sequence contents will be incomplete (a partially defined sequence;
see Section :ref:`sec:partial-seq`) if the alignment
does not extend over the full sequence.

The second line of this alignment block, starting with "``>``", shows
the name and description of the Hidden Markov Model from which the
target sequence was taken. These are stored under the keys
``"hmm_name"`` and ``"hmm_description"`` in the
``alignment.target.annotations`` dictionary:

.. cont-doctest

.. code:: pycon

   >>> alignment.target.annotations  # doctest: +NORMALIZE_WHITESPACE
   {'hmm_name': '2uvo_A',
    'hmm_description': 'Agglutinin isolectin 1; carbohydrate-binding protein, hevein domain, chitin-binding, GERM agglutinin, chitin-binding protein; HET: NDG NAG GOL; 1.40A {Triticum aestivum} PDB: 1wgc_A* 2cwg_A* 2x3t_A* 4aml_A* 7wga_A 9wga_A 2wgc_A 1wgt_A 1k7t_A* 1k7v_A* 1k7u_A 2x52_A* 1t0w_A*'}

The dictionary ``alignment.target.letter_annotations`` stores the target
alignent consensus sequence, the secondary structure as predicted by
PSIPRED, and the target secondary structure as determined by DSSP:

.. cont-doctest

.. code:: pycon

   >>> alignment.target.letter_annotations  # doctest: +NORMALIZE_WHITESPACE
   {'Consensus': '~~cg~~~~~~~c~~~~CCS~~g~Cg~~~~~Cg~gC~~~~c~~~~~cg~~~~~~~c~~~~CCs~~g~Cg~~~~~c~~~c~~~~~~~~~~cg~~~~~~~c~~~~CCs~~g~CG~~~~~C~~gCq~~~c~~~~~cg~~~~~~~c~~~~ccs~~g~Cg~~~~~C~~~cq~~~~~~',
    'ss_pred': 'CCCCCCCCCcCCCCCCeeCCCCeECCCcccccCCccccccccccccCcccCCcccCCccccCCCceeCCCccccCCCcccccccccccccccccCCCCCCCcccCCCCccCCCcccccCCCcCCccccccccccccccccCCCCCCcCCCCEecCchhhcccccccCCCCC',
    'ss_dssp': 'CBCBGGGTTBBCGGGCEECTTSBEEBSHHHHSTTCCBSSCSSCCBCBGGGTTBCCSTTCEECTTSBEEBSHHHHSTTCCBSSCSSCCBCBGGGTTBCCGGGCEECTTSBEEBSHHHHSTTCCBSSCSSCCCCBTTTTTBCCSTTCEECTTSCEEBSHHHHSTTCCBSSCC '}

In this example, for the query sequence only the consensus sequence is
available:

.. cont-doctest

.. code:: pycon

   >>> alignment.query.letter_annotations
   {'Consensus': '~~cg~~~~~~~c~~~~CCs~~g~CG~~~~~c~~~c~~~~c~~~~~Cg~~~~~~~c~~~~CCs~~g~CG~~~~~c~~~c~~~~~~~~~~Cg~~~~~~~c~~~~CCS~~G~CG~~~~~C~~~Cq~~~c~~~~~Cg~~~~~~~c~~~~CCS~~G~CG~~~~~C~~gCq~~~c~~'}

The ``alignment.annotations`` dictionary stores information about the
alignment shown on the third line of the alignment block:

.. cont-doctest

.. code:: pycon

   >>> alignment.annotations  # doctest: +NORMALIZE_WHITESPACE
   {'Probab': 99.95,
    'E-value': 3.7e-34,
    'Score': 210.31,
    'Identities': 100.0,
    'Similarity': 2.05,
    'Sum_probs': 166.9}

Confidence values for the pairwise alignment are stored under the
``"Confidence"`` key in the ``alignment.column_annotations`` dictionary.
This dictionary also stores the score for each column, shown between the
query and the target section of each alignment block:

.. cont-doctest

.. code:: pycon

   >>> alignment.column_annotations  # doctest: +NORMALIZE_WHITESPACE
   {'column score': '||||++.++..||++.|||+|+|||.+.+||+++||.+.|++..+|+++++.++|....|||.++||+.+.+||+.+||.+++++|+.|+...+++.||++.|||.|||||...+||+.+||+++|++|.+|++.+++++|..+.|||+++-||+...||+++||++.|||',
    'Confidence': '799999999999999999999999999999999999999999999999999999999999999999999999999999999999999999998899999999999999999999999999999999999999999999999999999999999999999999999999986'}

.. _`subsec:align_a2m`:

A2M
~~~

A2M files are alignment files created by ``align2model`` or ``hmmscore``
in the SAM Sequence Alignment and Modeling Software
System [Krogh1994]_, [Hughey1996]_. An A2M file contains
one multiple alignment. The A2M file format is similar to aligned FASTA
(see section :ref:`subsec:align_fasta`). However, to distinguish
insertions from deletions, A2M uses both dashes and periods to represent
gaps, and both upper and lower case characters in the aligned sequences.
Matches are represented by upper case letters and deletions by dashes in
alignment columns containing matches or deletions only. Insertions are
represented by lower case letters, with gaps aligned to the insertion
shown as periods. Header lines start with "``>``" followed by the name
of the sequence, and optionally a description.

The file ``probcons.a2m`` in Biopython’s test suite is an example of an
A2M file (see section :ref:`subsec:align_fasta` for the same
alignment in aligned FASTA format):

.. code:: text

   >plas_horvu
   D.VLLGANGGVLVFEPNDFSVKAGETITFKNNAGYPHNVVFDEDAVPSG.VD.VSKISQEEYLTAPGETFSVTLTV...PGTYGFYCEPHAGAGMVGKVT
   V
   >plas_chlre
   -.VKLGADSGALEFVPKTLTIKSGETVNFVNNAGFPHNIVFDEDAIPSG.VN.ADAISRDDYLNAPGETYSVKLTA...AGEYGYYCEPHQGAGMVGKII
   V
   >plas_anava
   -.VKLGSDKGLLVFEPAKLTIKPGDTVEFLNNKVPPHNVVFDAALNPAKsADlAKSLSHKQLLMSPGQSTSTTFPAdapAGEYTFYCEPHRGAGMVGKIT
   V
   >plas_proho
   VqIKMGTDKYAPLYEPKALSISAGDTVEFVMNKVGPHNVIFDK--VPAG.ES.APALSNTKLRIAPGSFYSVTLGT...PGTYSFYCTPHRGAGMVGTIT
   V
   >azup_achcy
   VhMLNKGKDGAMVFEPASLKVAPGDTVTFIPTDK-GHNVETIKGMIPDG.AE.A-------FKSKINENYKVTFTA...PGVYGVKCTPHYGMGMVGVVE
   V

To parse this alignment, use

.. doctest ../Tests/Clustalw lib:numpy

.. code:: pycon

   >>> from Bio import Align
   >>> alignment = Align.read("probcons.a2m", "a2m")
   >>> alignment  # doctest: +ELLIPSIS
   <Alignment object (5 rows x 101 columns) at ...>
   >>> print(alignment)
   plas_horv         0 D-VLLGANGGVLVFEPNDFSVKAGETITFKNNAGYPHNVVFDEDAVPSG-VD-VSKISQE
   plas_chlr         0 --VKLGADSGALEFVPKTLTIKSGETVNFVNNAGFPHNIVFDEDAIPSG-VN-ADAISRD
   plas_anav         0 --VKLGSDKGLLVFEPAKLTIKPGDTVEFLNNKVPPHNVVFDAALNPAKSADLAKSLSHK
   plas_proh         0 VQIKMGTDKYAPLYEPKALSISAGDTVEFVMNKVGPHNVIFDK--VPAG-ES-APALSNT
   azup_achc         0 VHMLNKGKDGAMVFEPASLKVAPGDTVTFIPTDK-GHNVETIKGMIPDG-AE-A------
   <BLANKLINE>
   plas_horv        57 EYLTAPGETFSVTLTV---PGTYGFYCEPHAGAGMVGKVTV 95
   plas_chlr        56 DYLNAPGETYSVKLTA---AGEYGYYCEPHQGAGMVGKIIV 94
   plas_anav        58 QLLMSPGQSTSTTFPADAPAGEYTFYCEPHRGAGMVGKITV 99
   plas_proh        56 KLRIAPGSFYSVTLGT---PGTYSFYCTPHRGAGMVGTITV 94
   azup_achc        51 -FKSKINENYKVTFTA---PGVYGVKCTPHYGMGMVGVVEV 88
   <BLANKLINE>

The parser analyzes the pattern of dashes, periods, and lower and upper
case letters in the A2M file to determine if a column is an
match/mismatch/deletion ("``D``") or an insertion ("``I``"). This
information is stored under the ``match`` key of the
``alignment.column_annotations`` dictionary:

.. cont-doctest

.. code:: pycon

   >>> alignment.column_annotations
   {'state': 'DIDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDIDDIDDDDDDDDDDDDDDDDDDDDDDDIIIDDDDDDDDDDDDDDDDDDDDDD'}

As the state information is stored in the ``alignment``, we can print
the alignment in the A2M format:

.. cont-doctest

.. code:: pycon

   >>> print(format(alignment, "a2m"))
   >plas_horvu
   D.VLLGANGGVLVFEPNDFSVKAGETITFKNNAGYPHNVVFDEDAVPSG.VD.VSKISQEEYLTAPGETFSVTLTV...PGTYGFYCEPHAGAGMVGKVTV
   >plas_chlre
   -.VKLGADSGALEFVPKTLTIKSGETVNFVNNAGFPHNIVFDEDAIPSG.VN.ADAISRDDYLNAPGETYSVKLTA...AGEYGYYCEPHQGAGMVGKIIV
   >plas_anava
   -.VKLGSDKGLLVFEPAKLTIKPGDTVEFLNNKVPPHNVVFDAALNPAKsADlAKSLSHKQLLMSPGQSTSTTFPAdapAGEYTFYCEPHRGAGMVGKITV
   >plas_proho
   VqIKMGTDKYAPLYEPKALSISAGDTVEFVMNKVGPHNVIFDK--VPAG.ES.APALSNTKLRIAPGSFYSVTLGT...PGTYSFYCTPHRGAGMVGTITV
   >azup_achcy
   VhMLNKGKDGAMVFEPASLKVAPGDTVTFIPTDK-GHNVETIKGMIPDG.AE.A-------FKSKINENYKVTFTA...PGVYGVKCTPHYGMGMVGVVEV
   <BLANKLINE>

Similarly, the alignment can be written in the A2M format to an output
file using ``Align.write`` (see
section :ref:`subsec:align_writing`).

.. _`subsec:align_mauve`:

Mauve eXtended Multi-FastA (xmfa) format
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Mauve [Darling2004]_ is a software package for
constructing multiple genome alignments. These alignments are stored in
the eXtended Multi-FastA (xmfa) format. Depending on how exactly
``progressiveMauve`` (the aligner program in Mauve) was called, the xmfa
format is slightly different.

If ``progressiveMauve`` is called with a single sequence input file, as
in

.. code:: text

   progressiveMauve combined.fasta  --output=combined.xmfa ...

where ``combined.fasta`` contains the genome sequences:

.. code:: text

   >equCab1
   GAAAAGGAAAGTACGGCCCGGCCACTCCGGGTGTGTGCTAGGAGGGCTTA
   >mm9
   GAAGAGGAAAAGTAGATCCCTGGCGTCCGGAGCTGGGACGT
   >canFam2
   CAAGCCCTGCGCGCTCAGCCGGAGTGTCCCGGGCCCTGCTTTCCTTTTC

then the output file ``combined.xmfa`` is as follows:

.. code:: text

   #FormatVersion Mauve1
   #Sequence1File  combined.fa
   #Sequence1Entry 1
   #Sequence1Format    FastA
   #Sequence2File  combined.fa
   #Sequence2Entry 2
   #Sequence2Format    FastA
   #Sequence3File  combined.fa
   #Sequence3Entry 3
   #Sequence3Format    FastA
   #BackboneFile   combined.xmfa.bbcols
   > 1:2-49 - combined.fa
   AAGCCCTCCTAGCACACACCCGGAGTGG-CCGGGCCGTACTTTCCTTTT
   > 2:0-0 + combined.fa
   -------------------------------------------------
   > 3:2-48 + combined.fa
   AAGCCCTGC--GCGCTCAGCCGGAGTGTCCCGGGCCCTGCTTTCCTTTT
   =
   > 1:1-1 + combined.fa
   G
   =
   > 1:50-50 + combined.fa
   A
   =
   > 2:1-41 + combined.fa
   GAAGAGGAAAAGTAGATCCCTGGCGTCCGGAGCTGGGACGT
   =
   > 3:1-1 + combined.fa
   C
   =
   > 3:49-49 + combined.fa
   C
   =

with numbers (1, 2, 3) referring to the input genome sequences for horse
(``equCab1``), mouse (``mm9``), and dog (``canFam2``), respectively.
This xmfa file consists of six alignment blocks, separated by ``=``
characters. Use ``Align.parse`` to extract these alignments:

.. doctest ../Tests/Mauve lib:numpy

.. code:: pycon

   >>> from Bio import Align
   >>> alignments = Align.parse("combined.xmfa", "mauve")

The file header data are stored in the ``metadata`` attribute:

.. cont-doctest

.. code:: pycon

   >>> alignments.metadata  # doctest: +NORMALIZE_WHITESPACE
   {'FormatVersion': 'Mauve1',
    'BackboneFile': 'combined.xmfa.bbcols',
    'File': 'combined.fa'}

The ``identifiers`` attribute stores the sequence identifiers for the
three sequences, which in this case is the three numbers:

.. cont-doctest

.. code:: pycon

   >>> alignments.identifiers
   ['0', '1', '2']

These identifiers are used in the individual alignments:

.. cont-doctest

.. code:: pycon

   >>> for alignment in alignments:
   ...     print([record.id for record in alignment.sequences])
   ...     print(alignment)
   ...     print("******")
   ...
   ['0', '1', '2']
   0                49 AAGCCCTCCTAGCACACACCCGGAGTGG-CCGGGCCGTACTTTCCTTTT  1
   1                 0 -------------------------------------------------  0
   2                 1 AAGCCCTGC--GCGCTCAGCCGGAGTGTCCCGGGCCCTGCTTTCCTTTT 48
   <BLANKLINE>
   ******
   ['0']
   0                 0 G 1
   <BLANKLINE>
   ******
   ['0']
   0                49 A 50
   <BLANKLINE>
   ******
   ['1']
   1                 0 GAAGAGGAAAAGTAGATCCCTGGCGTCCGGAGCTGGGACGT 41
   <BLANKLINE>
   ******
   ['2']
   2                 0 C 1
   <BLANKLINE>
   ******
   ['2']
   2                48 C 49
   <BLANKLINE>
   ******

Note that only the first block is a real alignment; the other blocks
contain only a single sequence. By including these blocks, the xmfa file
contains the full sequence that was provided in the ``combined.fa``
input file.

If ``progressiveMauve`` is called with a separate input file for each
genome, as in

.. code:: text

   progressiveMauve equCab1.fa canFam2.fa mm9.fa --output=separate.xmfa ...

where each Fasta file contains the genome sequence for one species only,
then the output file ``separate.xmfa`` is as follows:

.. code:: text

   #FormatVersion Mauve1
   #Sequence1File  equCab1.fa
   #Sequence1Format    FastA
   #Sequence2File  canFam2.fa
   #Sequence2Format    FastA
   #Sequence3File  mm9.fa
   #Sequence3Format    FastA
   #BackboneFile   separate.xmfa.bbcols
   > 1:1-50 - equCab1.fa
   TAAGCCCTCCTAGCACACACCCGGAGTGGCC-GGGCCGTAC-TTTCCTTTTC
   > 2:1-49 + canFam2.fa
   CAAGCCCTGC--GCGCTCAGCCGGAGTGTCCCGGGCCCTGC-TTTCCTTTTC
   > 3:1-19 - mm9.fa
   ---------------------------------GGATCTACTTTTCCTCTTC
   =
   > 3:20-41 + mm9.fa
   CTGGCGTCCGGAGCTGGGACGT
   =

The identifiers ``equCab1`` for horse, ``mm9`` for mouse, and
``canFam2`` for dog are now shown explicitly in the output file. This
xmfa file consists of two alignment blocks, separated by ``=``
characters. Use ``Align.parse`` to extract these alignments:

.. doctest ../Tests/Mauve lib:numpy

.. code:: pycon

   >>> from Bio import Align
   >>> alignments = Align.parse("separate.xmfa", "mauve")

The file header data now does not include the input file name:

.. cont-doctest

.. code:: pycon

   >>> alignments.metadata  # doctest: +NORMALIZE_WHITESPACE
   {'FormatVersion': 'Mauve1',
    'BackboneFile': 'separate.xmfa.bbcols'}

The ``identifiers`` attribute stores the sequence identifiers for the
three sequences:

.. cont-doctest

.. code:: pycon

   >>> alignments.identifiers
   ['equCab1.fa', 'canFam2.fa', 'mm9.fa']

These identifiers are used in the individual alignments:

.. cont-doctest

.. code:: pycon

   >>> for alignment in alignments:
   ...     print([record.id for record in alignment.sequences])
   ...     print(alignment)
   ...     print("******")
   ...
   ['equCab1.fa', 'canFam2.fa', 'mm9.fa']
   equCab1.f        50 TAAGCCCTCCTAGCACACACCCGGAGTGGCC-GGGCCGTAC-TTTCCTTTTC  0
   canFam2.f         0 CAAGCCCTGC--GCGCTCAGCCGGAGTGTCCCGGGCCCTGC-TTTCCTTTTC 49
   mm9.fa           19 ---------------------------------GGATCTACTTTTCCTCTTC  0
   <BLANKLINE>
   ******
   ['mm9.fa']
   mm9.fa           19 CTGGCGTCCGGAGCTGGGACGT 41
   <BLANKLINE>
   ******

To output the alignments in Mauve format, use ``Align.write``:

.. cont-doctest

.. code:: pycon

   >>> from io import StringIO
   >>> stream = StringIO()
   >>> alignments = Align.parse("separate.xmfa", "mauve")
   >>> Align.write(alignments, stream, "mauve")
   2
   >>> print(stream.getvalue())  # doctest: +NORMALIZE_WHITESPACE
   #FormatVersion Mauve1
   #Sequence1File  equCab1.fa
   #Sequence1Format    FastA
   #Sequence2File  canFam2.fa
   #Sequence2Format    FastA
   #Sequence3File  mm9.fa
   #Sequence3Format    FastA
   #BackboneFile   separate.xmfa.bbcols
   > 1:1-50 - equCab1.fa
   TAAGCCCTCCTAGCACACACCCGGAGTGGCC-GGGCCGTAC-TTTCCTTTTC
   > 2:1-49 + canFam2.fa
   CAAGCCCTGC--GCGCTCAGCCGGAGTGTCCCGGGCCCTGC-TTTCCTTTTC
   > 3:1-19 - mm9.fa
   ---------------------------------GGATCTACTTTTCCTCTTC
   =
   > 3:20-41 + mm9.fa
   CTGGCGTCCGGAGCTGGGACGT
   =
   <BLANKLINE>

Here, the writer makes use of the information stored in
``alignments.metadata`` and ``alignments.identifiers`` to create this
format. If your ``alignments`` object does not have these attributes,
you can provide them as keyword arguments to ``Align.write``:

.. cont-doctest

.. code:: pycon

   >>> stream = StringIO()
   >>> alignments = Align.parse("separate.xmfa", "mauve")
   >>> metadata = alignments.metadata
   >>> identifiers = alignments.identifiers
   >>> alignments = list(alignments)  # this drops the attributes
   >>> alignments.metadata  # doctest: +ELLIPSIS
   Traceback (most recent call last):
    ...
   AttributeError: 'list' object has no attribute 'metadata'
   >>> alignments.identifiers  # doctest: +ELLIPSIS
   Traceback (most recent call last):
    ...
   AttributeError: 'list' object has no attribute 'identifiers'
   >>> Align.write(alignments, stream, "mauve", metadata=metadata, identifiers=identifiers)
   2
   >>> print(stream.getvalue())  # doctest: +NORMALIZE_WHITESPACE
   #FormatVersion Mauve1
   #Sequence1File  equCab1.fa
   #Sequence1Format    FastA
   #Sequence2File  canFam2.fa
   #Sequence2Format    FastA
   #Sequence3File  mm9.fa
   #Sequence3Format    FastA
   #BackboneFile   separate.xmfa.bbcols
   > 1:1-50 - equCab1.fa
   TAAGCCCTCCTAGCACACACCCGGAGTGGCC-GGGCCGTAC-TTTCCTTTTC
   > 2:1-49 + canFam2.fa
   CAAGCCCTGC--GCGCTCAGCCGGAGTGTCCCGGGCCCTGC-TTTCCTTTTC
   > 3:1-19 - mm9.fa
   ---------------------------------GGATCTACTTTTCCTCTTC
   =
   > 3:20-41 + mm9.fa
   CTGGCGTCCGGAGCTGGGACGT
   =
   <BLANKLINE>

Python does not allow you to add these attributes to the ``alignments``
object directly, as in this example it was converted to a plain list.
However, you can construct an ``Alignments`` object and add attributes
to it (see Section :ref:`sec:alignments`):

.. cont-doctest

.. code:: pycon

   >>> alignments = Align.Alignments(alignments)
   >>> alignments.metadata = metadata
   >>> alignments.identifiers = identifiers
   >>> stream = StringIO()
   >>> Align.write(alignments, stream, "mauve", metadata=metadata, identifiers=identifiers)
   2
   >>> print(stream.getvalue())  # doctest: +NORMALIZE_WHITESPACE
   #FormatVersion Mauve1
   #Sequence1File  equCab1.fa
   #Sequence1Format    FastA
   #Sequence2File  canFam2.fa
   #Sequence2Format    FastA
   #Sequence3File  mm9.fa
   #Sequence3Format    FastA
   #BackboneFile   separate.xmfa.bbcols
   > 1:1-50 - equCab1.fa
   TAAGCCCTCCTAGCACACACCCGGAGTGGCC-GGGCCGTAC-TTTCCTTTTC
   > 2:1-49 + canFam2.fa
   CAAGCCCTGC--GCGCTCAGCCGGAGTGTCCCGGGCCCTGC-TTTCCTTTTC
   > 3:1-19 - mm9.fa
   ---------------------------------GGATCTACTTTTCCTCTTC
   =
   > 3:20-41 + mm9.fa
   CTGGCGTCCGGAGCTGGGACGT
   =
   <BLANKLINE>

When printing a single alignment in ``Mauve`` format, use keyword
arguments to provide the metadata and identifiers:

.. cont-doctest

.. code:: pycon

   >>> alignment = alignments[0]
   >>> print(alignment.format("mauve", metadata=metadata, identifiers=identifiers))
   > 1:1-50 - equCab1.fa
   TAAGCCCTCCTAGCACACACCCGGAGTGGCC-GGGCCGTAC-TTTCCTTTTC
   > 2:1-49 + canFam2.fa
   CAAGCCCTGC--GCGCTCAGCCGGAGTGTCCCGGGCCCTGC-TTTCCTTTTC
   > 3:1-19 - mm9.fa
   ---------------------------------GGATCTACTTTTCCTCTTC
   =
   <BLANKLINE>

.. _`subsec:align_sam`:

Sequence Alignment/Map (SAM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Files in the Sequence Alignment/Map (SAM) format
[Li2009]_ store pairwise sequence alignments, usually
of next-generation sequencing data against a reference genome. The file
``ex1.sam`` in Biopython’s test suite is an example of a minimal file in
the SAM format. Its first few lines are as follows:

.. code:: text

   EAS56_57:6:190:289:82   69      chr1    100     0       *       =       100     0       CTCAAGGTTGTTGCAAGGGGGTCTATGTGAACAAA     <<<7<<<;<<<<<<<<8;;<7;4<;<;;;;;94<;     MF:i:192
   EAS56_57:6:190:289:82   137     chr1    100     73      35M     =       100     0       AGGGGTGCAGAGCCGAGTCACGGGGTTGCCAGCAC     <<<<<<;<<<<<<<<<<;<<;<<<<;8<6;9;;2;     MF:i:64 Aq:i:0  NM:i:0  UQ:i:0  H0:i:1  H1:i:0
   EAS51_64:3:190:727:308  99      chr1    103     99      35M     =       263     195     GGTGCAGAGCCGAGTCACGGGGTTGCCAGCACAGG     <<<<<<<<<<<<<<<<<<<<<<<<<<<::<<<844     MF:i:18 Aq:i:73 NM:i:0  UQ:i:0  H0:i:1  H1:i:0
   ...

To parse this file, use

.. doctest ../Tests/SamBam lib:numpy

.. code:: pycon

   >>> from Bio import Align
   >>> alignments = Align.parse("ex1.sam", "sam")
   >>> alignment = next(alignments)

The ``flag`` of the first line is 69. According to the SAM/BAM file
format specification, lines for which the flag contains the bitwise flag
4 are unmapped. As 69 has the bit corresponding to this position set to
True, this sequence is unmapped and was not aligned to the genome (in
spite of the first line showing ``chr1``). The target of this alignment
(or the first item in ``alignment.sequences``) is therefore ``None``:

.. cont-doctest

.. code:: pycon

   >>> alignment.flag
   69
   >>> bin(69)
   '0b1000101'
   >>> bin(4)
   '0b100'
   >>> if alignment.flag & 4:
   ...     print("unmapped")
   ... else:
   ...     print("mapped")
   ...
   unmapped
   >>> alignment.sequences
   [None, SeqRecord(seq=Seq('CTCAAGGTTGTTGCAAGGGGGTCTATGTGAACAAA'), id='EAS56_57:6:190:289:82', name='<unknown name>', description='', dbxrefs=[])]
   >>> alignment.target is None
   True

The second line represents an alignment to chromosome 1:

.. cont-doctest

.. code:: pycon

   >>> alignment = next(alignments)
   >>> if alignment.flag & 4:
   ...     print("unmapped")
   ... else:
   ...     print("mapped")
   ...
   mapped
   >>> alignment.target
   SeqRecord(seq=None, id='chr1', name='<unknown name>', description='', dbxrefs=[])

As this SAM file does not store the genome sequence information for each
alignment, we cannot print the alignment. However, we can print the
alignment information in SAM format or any other format (such as BED,
see section :ref:`subsec:align_bed`) that does not require the
target sequence information:

.. cont-doctest

.. code:: pycon

   >>> format(alignment, "sam")
   'EAS56_57:6:190:289:82\t137\tchr1\t100\t73\t35M\t=\t100\t0\tAGGGGTGCAGAGCCGAGTCACGGGGTTGCCAGCAC\t<<<<<<;<<<<<<<<<<;<<;<<<<;8<6;9;;2;\tMF:i:64\tAq:i:0\tNM:i:0\tUQ:i:0\tH0:i:1\tH1:i:0\n'
   >>> format(alignment, "bed")
   'chr1\t99\t134\tEAS56_57:6:190:289:82\t0\t+\t99\t134\t0\t1\t35,\t0,\n'

However, we cannot print the alignment in PSL format (see
section :ref:`subsec:align_psl`) as that would require knowing
the size of the target sequence chr1:

.. cont-doctest

.. code:: pycon

   >>> format(alignment, "psl")  # doctest: +ELLIPSIS
   Traceback (most recent call last):
    ...
   TypeError: ...

If you know the size of the target sequences, you can set them by hand:

.. cont-doctest

.. code:: pycon

   >>> from Bio.Seq import Seq
   >>> from Bio.SeqRecord import SeqRecord
   >>> target = SeqRecord(Seq(None, length=1575), id="chr1")
   >>> alignment.target = target
   >>> format(alignment, "psl")  # doctest: +ELLIPSIS
   '35\t0\t0\t0\t0\t0\t0\t0\t+\tEAS56_57:6:190:289:82\t35\t0\t35\tchr1\t1575\t99\t134\t1\t35,\t0,\t99,\n'

The file ``ex1_header.sam`` in Biopython’s test suite contains the same
alignments, but now also includes a header. Its first few lines are as
follows:

.. code:: text

   @HD\tVN:1.3\tSO:coordinate
   @SQ\tSN:chr1\tLN:1575
   @SQ\tSN:chr2\tLN:1584
   EAS56_57:6:190:289:82   69      chr1    100     0       *       =       100     0       CTCAAGGTTGTTGCAAGGGGGTCTATGTGAACAAA     <<<7<<<;<<<<<<<<8;;<7;4<;<;;;;;94<;     MF:i:192
   ...

The header stores general information about the alignments, including
the size of the target chromosomes. The target information is stored in
the ``targets`` attribute of the ``alignments`` object:

.. doctest ../Tests/SamBam lib:numpy

.. code:: pycon

   >>> from Bio import Align
   >>> alignments = Align.parse("ex1_header.sam", "sam")
   >>> len(alignments.targets)
   2
   >>> alignments.targets[0]
   SeqRecord(seq=Seq(None, length=1575), id='chr1', name='<unknown name>', description='', dbxrefs=[])
   >>> alignments.targets[1]
   SeqRecord(seq=Seq(None, length=1584), id='chr2', name='<unknown name>', description='', dbxrefs=[])

Other information provided in the header is stored in the ``metadata``
attribute:

.. cont-doctest

.. code:: pycon

   >>> alignments.metadata
   {'HD': {'VN': '1.3', 'SO': 'coordinate'}}

With the target information, we can now also print the alignment in PSL
format:

.. cont-doctest

.. code:: pycon

   >>> alignment = next(alignments)  # the unmapped sequence; skip it
   >>> alignment = next(alignments)
   >>> format(alignment, "psl")
   '35\t0\t0\t0\t0\t0\t0\t0\t+\tEAS56_57:6:190:289:82\t35\t0\t35\tchr1\t1575\t99\t134\t1\t35,\t0,\t99,\n'

We can now also print the alignment in human-readable form, but note
that the target sequence contents is not available from this file:

.. cont-doctest

.. code:: pycon

   >>> print(alignment)
   chr1             99 ??????????????????????????????????? 134
                     0 ...................................  35
   EAS56_57:         0 AGGGGTGCAGAGCCGAGTCACGGGGTTGCCAGCAC  35
   <BLANKLINE>

Alignments in the file ``sam1.sam`` in the Biopython test suite contain
an additional ``MD`` tag that shows how the query sequence differs from
the target sequence:

.. code:: text

   @SQ     SN:1    LN:239940
   @PG     ID:bwa  PN:bwa  VN:0.6.2-r126
   HWI-1KL120:88:D0LRBACXX:1:1101:1780:2146        77      *       0       0       *       *       0       0       GATGGGAAACCCATGGCCGAGTGGGAAGAAACCAGCTGAGGTCACATCACCAGAGGAGGGAGAGTGTGGCCCCTGACTCAGTCCATCAGCTTGTGGAGCTG   @=?DDDDBFFFF7A;E?GGEGE8BB?FF?F>G@F=GIIDEIBCFF<FEFEC@EEEE2?8B8/=@((-;?@2<B9@##########################
   ...
   HWI-1KL120:88:D0LRBACXX:1:1101:2852:2134        137     1       136186  25      101M    =       136186  0       TCACGGTGGCCTGTTGAGGCAGGGGCTCACGCTGACCTCTCTCGGCGTGGGAGGGGCCGGTGTGAGGCAAGGGCTCACGCTGACCTCTCTCGGCGTGGGAG   @C@FFFDFHGHHHJJJIJJJJIJJJGEDHHGGHGBGIIGIIAB@GEE=BDBBCCDD@D@B7@;@DDD?<A?DD728:>8()009>:>>C@>5??B######   XT:A:U  NM:i:5  SM:i:25 AM:i:0  X0:i:1  X1:i:0  XM:i:5  XO:i:0  XG:i:0  MD:Z:25G14G2C34A12A9

The parser reconstructs the local genome sequence from the ``MD`` tag,
allowing us to see the target sequence explicitly when printing the
alignment:

.. doctest ../Tests/SamBam lib:numpy

.. code:: pycon

   >>> from Bio import Align
   >>> alignments = Align.parse("sam1.sam", "sam")
   >>> for alignment in alignments:
   ...     if not alignment.flag & 4:  # Skip the unmapped lines
   ...         break
   ...
   >>> alignment  # doctest: +ELLIPSIS
   <Alignment object (2 rows x 101 columns) at ...>
   >>> print(alignment)
   1            136185 TCACGGTGGCCTGTTGAGGCAGGGGGTCACGCTGACCTCTGTCCGCGTGGGAGGGGCCGG
                     0 |||||||||||||||||||||||||.||||||||||||||.||.||||||||||||||||
   HWI-1KL12         0 TCACGGTGGCCTGTTGAGGCAGGGGCTCACGCTGACCTCTCTCGGCGTGGGAGGGGCCGG
   <BLANKLINE>
   1            136245 TGTGAGGCAAGGGCTCACACTGACCTCTCTCAGCGTGGGAG 136286
                    60 ||||||||||||||||||.||||||||||||.|||||||||    101
   HWI-1KL12        60 TGTGAGGCAAGGGCTCACGCTGACCTCTCTCGGCGTGGGAG    101
   <BLANKLINE>

SAM files may include additional information to distinguish simple
sequence insertions and deletions from skipped regions of the genome
(e.g. introns), hard and soft clipping, and padded sequence regions. As
this information cannot be stored in the ``coordinates`` attribute of an
``Alignment`` object, and is stored in a dedicated ``operations``
attribute instead. Let’s use the third alignment in this SAM file as an
example:

.. doctest ../Tests/Blat lib:numpy

.. code:: pycon

   >>> from Bio import Align
   >>> alignments = Align.parse("dna_rna.sam", "sam")
   >>> alignment = next(alignments)
   >>> alignment = next(alignments)
   >>> alignment = next(alignments)
   >>> print(format(alignment, "SAM"))  # doctest: +NORMALIZE_WHITESPACE
   NR_111921.1 0   chr3    48663768    0   46M1827N82M3376N76M12H  *   0   0   CACGAGAGGAGCGGAGGCGAGGGGTGAACGCGGAGCACTCCAATCGCTCCCAACTAGAGGTCCACCCAGGACCCAGAGACCTGGATTTGAGGCTGCTGGGCGGCAGATGGAGCGATCAGAAGACCAGGAGACGGGAGCTGGAGTGCAGTGGCTGTTCACAAGCGTGAAAGCAAAGATTAAAAAATTTGTTTTTATATTAAAAAA    *   AS:i:1000   NM:i:0
   <BLANKLINE>
   >>> print(alignment.coordinates)
   [[48663767 48663813 48665640 48665722 48669098 48669174]
    [       0       46       46      128      128      204]]
   >>> alignment.operations
   bytearray(b'MNMNM')
   >>> alignment.query.annotations["hard_clip_right"]
   12

In this alignment, the cigar string ``63M1062N75M468N43M`` defines 46
aligned nucleotides, an intron of 1827 nucleotides, 82 aligned
nucleotides, an intron of 3376 nucleotides, 76 aligned nucleotides, and
12 hard-clipped nucleotides. These operations are shown in the
``operations`` attribute, except for hard-clipping, which is stored in
``alignment.query.annotations["hard_clip_right"]`` (or
``alignment.query.annotations["hard_clip_left"]``, if applicable)
instead.

To write a SAM file with alignments created from scratch, use an
``Alignments`` (plural) object (see Section :ref:`sec:alignments`)
to store the alignments as well as the metadata and targets:

.. doctest . lib:numpy

.. code:: pycon

   >>> from io import StringIO
   >>> import numpy as np

   >>> from Bio import Align
   >>> from Bio.Seq import Seq
   >>> from Bio.SeqRecord import SeqRecord

   >>> alignments = Align.Alignments()

   >>> seq1 = Seq(None, length=10000)
   >>> target1 = SeqRecord(seq1, id="chr1")
   >>> seq2 = Seq(None, length=15000)
   >>> target2 = SeqRecord(seq2, id="chr2")
   >>> alignments.targets = [target1, target2]
   >>> alignments.metadata = {"HD": {"VN": "1.3", "SO": "coordinate"}}

   >>> seqA = Seq(None, length=20)
   >>> queryA = SeqRecord(seqA, id="readA")
   >>> sequences = [target1, queryA]
   >>> coordinates = np.array([[4300, 4320], [0, 20]])
   >>> alignment = Align.Alignment(sequences, coordinates)
   >>> alignments.append(alignment)

   >>> seqB = Seq(None, length=25)
   >>> queryB = SeqRecord(seqB, id="readB")
   >>> sequences = [target1, queryB]
   >>> coordinates = np.array([[5900, 5925], [25, 0]])
   >>> alignment = Align.Alignment(sequences, coordinates)
   >>> alignments.append(alignment)

   >>> seqC = Seq(None, length=40)
   >>> queryC = SeqRecord(seqC, id="readC")
   >>> sequences = [target2, queryC]
   >>> coordinates = np.array([[12300, 12318], [0, 18]])
   >>> alignment = Align.Alignment(sequences, coordinates)
   >>> alignments.append(alignment)

   >>> stream = StringIO()
   >>> Align.write(alignments, stream, "sam")
   3
   >>> print(stream.getvalue())  # doctest: +NORMALIZE_WHITESPACE
   @HD VN:1.3  SO:coordinate
   @SQ SN:chr1 LN:10000
   @SQ SN:chr2 LN:15000
   readA   0   chr1    4301    255 20M *   0   0   *   *
   readB   16  chr1    5901    255 25M *   0   0   *   *
   readC   0   chr2    12301   255 18M22S  *   0   0   *       *
   <BLANKLINE>

.. _`subsec:align_bed`:

Browser Extensible Data (BED)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BED (Browser Extensible Data) files are typically used to store the
alignments of gene transcripts to the genome. See the `description from
UCSC <http://genome.cse.ucsc.edu/FAQ/FAQformat.html#format1>`__ for a
full explanation of the BED format.

BED files have three required fields and nine optional fields. The file
``bed12.bed`` in subdirectory ``Tests/Blat`` is an example of a BED file
with 12 fields:

.. code:: text

   chr22   1000    5000    mRNA1   960 +   1200    4900    255,0,0 2   567,488,    0,3512,
   chr22   2000    6000    mRNA2   900 -   2300    5960    0,255,0 2   433,399,    0,3601,

To parse this file, use

.. doctest ../Tests/Blat lib:numpy

.. code:: pycon

   >>> from Bio import Align
   >>> alignments = Align.parse("bed12.bed", "bed")
   >>> len(alignments)
   2
   >>> for alignment in alignments:
   ...     print(alignment.coordinates)
   ...
   [[1000 1567 4512 5000]
    [   0  567  567 1055]]
   [[2000 2433 5601 6000]
    [ 832  399  399    0]]

Note that the first sequence ("``mRNA1``") was mapped to the forward
strand, while the second sequence ("``mRNA2``") was mapped to the
reverse strand.

As a BED file does not store the length of each chromosome, the length
of the target sequence is set to its maximum:

.. code:: pycon

   >>> alignment.target
   SeqRecord(seq=Seq(None, length=9223372036854775807), id='chr22', name='<unknown name>', description='', dbxrefs=[])

The length of the query sequence can be inferred from its alignment
information:

.. cont-doctest

.. code:: pycon

   >>> alignment.query
   SeqRecord(seq=Seq(None, length=832), id='mRNA2', name='<unknown name>', description='', dbxrefs=[])

The alignment score (field 5) and information stored in fields 7-9
(referred to as ``thickStart``, ``thickEnd``, and ``itemRgb`` in the BED
format specification) are stored as attributes on the ``alignment``
object:

.. cont-doctest

.. code:: pycon

   >>> alignment.score
   900.0
   >>> alignment.thickStart
   2300
   >>> alignment.thickEnd
   5960
   >>> alignment.itemRgb
   '0,255,0'

To print an alignment in the BED format, you can use Python’s built-in
``format`` function:

.. cont-doctest

.. code:: pycon

   >>> print(format(alignment, "bed"))  # doctest: +NORMALIZE_WHITESPACE
   chr22   2000    6000    mRNA2   900 -   2300    5960    0,255,0 2   433,399,    0,3601,
   <BLANKLINE>

or you can use the ``format`` method of the ``alignment`` object. This
allows you to specify the number of fields to be written as the ``bedN``
keyword argument:

.. cont-doctest

.. code:: pycon

   >>> print(alignment.format("bed"))  # doctest: +NORMALIZE_WHITESPACE
   chr22   2000    6000    mRNA2   900 -   2300    5960    0,255,0 2   433,399,    0,3601,
   <BLANKLINE>
   >>> print(alignment.format("bed", 3))  # doctest: +NORMALIZE_WHITESPACE
   chr22   2000    6000
   <BLANKLINE>
   >>> print(alignment.format("bed", 6))  # doctest: +NORMALIZE_WHITESPACE
   chr22   2000    6000    mRNA2   900 -
   <BLANKLINE>

The same keyword argument can be used with ``Align.write``:

.. code:: pycon

   >>> Align.write(alignments, "mybed3file.bed", "bed", bedN=3)
   2
   >>> Align.write(alignments, "mybed6file.bed", "bed", bedN=6)
   2
   >>> Align.write(alignments, "mybed12file.bed", "bed")
   2

.. _`subsec:align_bigbed`:

bigBed
~~~~~~

The bigBed file format is an indexed binary version of a BED
file :ref:`subsec:align_bed`. To create a bigBed file, you can
either use the ``bedToBigBed`` program from UCSC
(`) <https://genome.ucsc.edu/goldenPath/help/bigBed.html>`__. or you can
use Biopython for it by calling the ``Bio.Align.write`` function with
``fmt="bigbed"``. While the two methods should result in identical
bigBed files, using ``bedToBigBed`` is much faster and may be more
reliable, as it is the gold standard. As bigBed files come with a
built-in index, it allows you to quickly search a specific genomic
region.

As an example, let’s parse the bigBed file ``dna_rna.bb``, available in
the ``Tests/Blat`` subdirectory in the Biopython distribution:

.. doctest ../Tests/Blat lib:numpy

.. code:: pycon

   >>> from Bio import Align
   >>> alignments = Align.parse("dna_rna.bb", "bigbed")
   >>> len(alignments)
   4
   >>> print(alignments.declaration)  # doctest: +NORMALIZE_WHITESPACE
   table bed
   "Browser Extensible Data"
   (
      string          chrom;          "Reference sequence chromosome or scaffold"
      uint            chromStart;     "Start position in chromosome"
      uint            chromEnd;       "End position in chromosome"
      string          name;           "Name of item."
      uint            score;          "Score (0-1000)"
      char[1]         strand;         "+ or - for strand"
      uint            thickStart;     "Start of where display should be thick (start codon)"
      uint            thickEnd;       "End of where display should be thick (stop codon)"
      uint            reserved;       "Used as itemRgb as of 2004-11-22"
      int             blockCount;     "Number of blocks"
      int[blockCount] blockSizes;     "Comma separated list of block sizes"
      int[blockCount] chromStarts;    "Start positions relative to chromStart"
   )
   <BLANKLINE>

The ``declaration`` contains the specification of the columns, in
AutoSql format, that was used to create the bigBed file. Target
sequences (typically, the chromosomes against which the sequences were
aligned) are stored in the ``targets`` attribute. In the bigBed format,
only the identifier and the size of each target is stored. In this
example, there is only a single chromosome:

.. cont-doctest

.. code:: pycon

   >>> alignments.targets
   [SeqRecord(seq=Seq(None, length=198295559), id='chr3', name='<unknown name>', description='<unknown description>', dbxrefs=[])]

Let’s look at the individual alignments. The alignment information is
stored in the same way as for a BED file (see section
:ref:`subsec:align_bed`):

.. cont-doctest

.. code:: pycon

   >>> alignment = next(alignments)
   >>> alignment.target.id
   'chr3'
   >>> alignment.query.id
   'NR_046654.1'
   >>> alignment.coordinates
   array([[42530895, 42530958, 42532020, 42532095, 42532563, 42532606],
          [     181,      118,      118,       43,       43,        0]])
   >>> alignment.thickStart
   42530895
   >>> alignment.thickEnd
   42532606
   >>> print(alignment)  # doctest: +ELLIPSIS
   chr3       42530895 ????????????????????????????????????????????????????????????
                     0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
   NR_046654       181 ????????????????????????????????????????????????????????????
   <BLANKLINE>
   chr3       42530955 ????????????????????????????????????????????????????????????
                    60 |||---------------------------------------------------------
   NR_046654       121 ???---------------------------------------------------------
   ...
   chr3       42532515 ????????????????????????????????????????????????????????????
                  1620 ------------------------------------------------||||||||||||
   NR_046654        43 ------------------------------------------------????????????
   <BLANKLINE>
   chr3       42532575 ??????????????????????????????? 42532606
                  1680 |||||||||||||||||||||||||||||||     1711
   NR_046654        31 ???????????????????????????????        0
   <BLANKLINE>

The default bigBed format does not store the sequence contents of the
target and query. If these are available elsewhere (for example, a Fasta
file), you can set ``alignment.target.seq`` and ``alignment.query.seq``
to show the sequence contents when printing the alignment, or to write
the alignment in formats that require the sequence contents (such as
Clustal, see section :ref:`subsec:align_clustal`). The test script
``test_Align_bigbed.py`` in the ``Tests`` subdirectory in the Biopython
distribution gives some examples on how to do that.

Now let’s see how to search for a sequence region. These are the
sequences stored in the bigBed file, printed in BED format (see section
:ref:`subsec:align_bed`):

.. cont-doctest

.. code:: pycon

   >>> for alignment in alignments:
   ...     print(format(alignment, "bed"))  # doctest: +NORMALIZE_WHITESPACE
   ...
   chr3    42530895    42532606    NR_046654.1 1000    -   42530895    42532606    0   3   63,75,43,   0,1125,1668,
   <BLANKLINE>
   chr3    42530895    42532606    NR_046654.1_modified    978 -   42530895    42532606    0   5   27,36,17,56,43, 0,27,1125,1144,1668,
   <BLANKLINE>
   chr3    48663767    48669174    NR_111921.1 1000    +   48663767    48669174    0   3   46,82,76,   0,1873,5331,
   <BLANKLINE>
   chr3    48663767    48669174    NR_111921.1_modified    972 +   48663767    48669174    0   5   28,17,76,6,76,  0,29,1873,1949,5331,
   <BLANKLINE>

Use the ``search`` method on the ``alignments`` object to find regions
on chr3 between positions 48000000 and 49000000. This method returns an
iterator:

.. cont-doctest

.. code:: pycon

   >>> selected_alignments = alignments.search("chr3", 48000000, 49000000)
   >>> for alignment in selected_alignments:
   ...     print(alignment.query.id)
   ...
   NR_111921.1
   NR_111921.1_modified

The chromosome name may be ``None`` to include all chromosomes, and the
start and end positions may be ``None`` to start searching from position
0 or to continue searching until the end of the chromosome,
respectively.

Writing alignments in the bigBed format is as easy as calling
``Bio.Align.write``:

.. code:: pycon

   >>> Align.write(alignments, "output.bb", "bigbed")

You can specify the number of BED fields to be included in the bigBed
file. For example, to write a BED6 file, use

.. code:: pycon

   >>> Align.write(alignments, "output.bb", "bigbed", bedN=6)

Same as for ``bedToBigBed``, you can include additional columns in the
bigBed output. Suppose the file ``bedExample2.as`` (available in the
``Tests/Blat`` subdirectory of the Biopython distribution) stores the
declaration of the included BED fields in AutoSql format. We can read
this declaration as follows:

.. cont-doctest

.. code:: pycon

   >>> from Bio.Align import bigbed
   >>> with open("bedExample2.as") as stream:
   ...     autosql_data = stream.read()
   ...
   >>> declaration = bigbed.AutoSQLTable.from_string(autosql_data)
   >>> type(declaration)
   <class 'Bio.Align.bigbed.AutoSQLTable'>
   >>> print(declaration)
   table hg18KGchr7
   "UCSC Genes for chr7 with color plus GeneSymbol and SwissProtID"
   (
      string  chrom;         "Reference sequence chromosome or scaffold"
      uint    chromStart;    "Start position of feature on chromosome"
      uint    chromEnd;      "End position of feature on chromosome"
      string  name;          "Name of gene"
      uint    score;         "Score"
      char[1] strand;        "+ or - for strand"
      uint    thickStart;    "Coding region start"
      uint    thickEnd;      "Coding region end"
      uint    reserved;      "Green on + strand, Red on - strand"
      string  geneSymbol;    "Gene Symbol"
      string  spID;          "SWISS-PROT protein Accession number"
   )
   <BLANKLINE>

Now we can write a bigBed file with the 9 BED fields plus the additional
fields ``geneSymbol`` and ``spID`` by calling

.. code:: pycon

   >>> Align.write(
   ...     alignments,
   ...     "output.bb",
   ...     "bigbed",
   ...     bedN=9,
   ...     declaration=declaration,
   ...     extraIndex=["name", "geneSymbol"],
   ... )

Here, we also requested to include additional indices on the ``name``
and ``geneSymbol`` in the bigBed file. ``Align.write`` expects to find
the keys ``geneSymbol`` and ``spID`` in the ``alignment.annotations``
dictionary. Please refer to the test script ``test_Align_bigbed.py`` in
the ``Tests`` subdirectory in the Biopython distribution for more
examples of writing alignment files in the bigBed format.

Optional arguments are ``compress`` (default value is ``True``), ``blockSize``
(default value is 256), and ``itemsPerSlot`` (default value is 512). See the
documentation of UCSC's ``bedToBigBed`` program for a description of these
arguments.  Searching a ``bigBed`` file can be faster by using
``compress=False`` and ``itemsPerSlot=1`` when creating the bigBed file.

.. _`subsec:align_psl`:

Pattern Space Layout (PSL)
~~~~~~~~~~~~~~~~~~~~~~~~~~

PSL (Pattern Space Layout) files are are generated by the BLAST-Like
Alignment Tool BLAT [Kent2002]_. Like BED files (see
section :ref:`subsec:align_bed`), PSL files are typically used to
store alignments of transcripts to genomes. This is an example of a
short BLAT file (available as ``dna_rna.psl`` in the ``Tests/Blat``
subdirectory of the Biopython distribution), with the standard PSL
header consisting of 5 lines:

.. code:: text

   psLayout version 3

   match   mis-    rep.    N's Q gap   Q gap   T gap   T gap   strand  Q           Q       Q       Q   T           T       T       T   block   blockSizes  qStarts  tStarts
           match   match       count   bases   count   bases           name        size    start   end name        size    start   end count
   ---------------------------------------------------------------------------------------------------------------------------------------------------------------
   165 0   39  0   0   0   2   5203    +   NR_111921.1 216 0   204 chr3    198295559   48663767    48669174    3   46,82,76,   0,46,128,   48663767,48665640,48669098,
   175 0   6   0   0   0   2   1530    -   NR_046654.1 181 0   181 chr3    198295559   42530895    42532606    3   63,75,43,   0,63,138,   42530895,42532020,42532563,
   162 2   39  0   1   2   3   5204    +   NR_111921.1_modified    220 3   208 chr3    198295559   48663767    48669174    5   28,17,76,6,76,  3,31,48,126,132,    48663767,48663796,48665640,48665716,48669098,
   172 1   6   0   1   3   3   1532    -   NR_046654.1_modified    190 3   185 chr3    198295559   42530895    42532606    5   27,36,17,56,43, 5,35,71,88,144, 42530895,42530922,42532020,42532039,42532563,

To parse this file, use

.. doctest ../Tests/Blat lib:numpy

.. code:: pycon

   >>> from Bio import Align
   >>> alignments = Align.parse("dna_rna.psl", "psl")
   >>> alignments.metadata
   {'psLayout version': '3'}

Iterate over the alignments to get one ``Alignment`` object for each
line:

.. cont-doctest

.. code:: pycon

   >>> for alignment in alignments:
   ...     print(alignment.target.id, alignment.query.id)
   ...
   chr3 NR_046654.1
   chr3 NR_046654.1_modified
   chr3 NR_111921.1
   chr3 NR_111921.1_modified

Let’s look at the last alignment in more detail. The first four columns
in the PSL file show the number of matches, the number of mismatches,
the number of nucleotides aligned to repeat regions, and the number of
nucleotides aligned to N (unknown) characters. These values are stored
as attributes to the ``Alignment`` object:

.. cont-doctest

.. code:: pycon

   >>> alignment.matches
   162
   >>> alignment.misMatches
   2
   >>> alignment.repMatches
   39
   >>> alignment.nCount
   0

As the sequence data of the target and query are not stored explicitly
in the PSL file, the sequence content of ``alignment.target`` and
``alignment.query`` is undefined. However, their sequence lengths are
known:

.. cont-doctest

.. code:: pycon

   >>> alignment.target  # doctest: +ELLIPSIS
   SeqRecord(seq=Seq(None, length=198295559), id='chr3', ...)
   >>> alignment.query  # doctest: +ELLIPSIS
   SeqRecord(seq=Seq(None, length=220), id='NR_111921.1_modified', ...)

We can print the alignment in BED or PSL format:

.. cont-doctest

.. code:: pycon

   >>> print(format(alignment, "bed"))  # doctest: +NORMALIZE_WHITESPACE
   chr3    48663767    48669174    NR_111921.1_modified    0   +   48663767    48669174    0   5   28,17,76,6,76,  0,29,1873,1949,5331,
   <BLANKLINE>
   >>> print(format(alignment, "psl"))  # doctest: +NORMALIZE_WHITESPACE
   162 2   39  0   1   2   3   5204    +   NR_111921.1_modified    220 3   208 chr3    198295559   48663767    48669174    5   28,17,76,6,76,  3,31,48,126,132,    48663767,48663796,48665640,48665716,48669098,
   <BLANKLINE>

Here, the number of matches, mismatches, repeat region matches, and
matches to unknown nucleotides were taken from the corresponding
attributes of the ``Alignment`` object. If these attributes are not
available, for example if the alignment did not come from a PSL file,
then these numbers are calculated using the sequence contents, if
available. Repeat regions in the target sequence are indicated by
masking the sequence as lower-case or upper-case characters, as defined
by the following values for the ``mask`` keyword argument:

-  ``False`` (default): Do not count matches to masked sequences
   separately;

-  ``"lower"``: Count and report matches to lower-case characters as
   matches to repeat regions;

-  ``"upper"``: Count and report matches to upper-case characters as
   matches to repeat regions;

The character used for unknown nucleotides is defined by the
``wildcard`` argument. For consistency with BLAT, the wildcard character
is ``"N"`` by default. Use ``wildcard=None`` if you don’t want to count
matches to any unknown nucleotides separately.

.. doctest . lib:numpy

.. code:: pycon

   >>> import numpy as np
   >>> from Bio import Align
   >>> query = "GGTGGGGG"
   >>> target = "AAAAAAAggggGGNGAAAAA"
   >>> coordinates = np.array([[0, 7, 15, 20], [0, 0, 8, 8]])
   >>> alignment = Align.Alignment([target, query], coordinates)
   >>> print(alignment)
   target            0 AAAAAAAggggGGNGAAAAA 20
                     0 -------....||.|----- 20
   query             0 -------GGTGGGGG-----  8
   <BLANKLINE>
   >>> line = alignment.format("psl")
   >>> print(line)  # doctest: +NORMALIZE_WHITESPACE
   6   1   0   1   0   0   0   0   +   query   8   0   8   target   20   7   15   1   8,   0,   7,
   >>> line = alignment.format("psl", mask="lower")
   >>> print(line)  # doctest: +NORMALIZE_WHITESPACE
   3   1   3   1   0   0   0   0   +   query   8   0   8   target   20   7   15   1   8,   0,   7,
   >>> line = alignment.format("psl", mask="lower", wildcard=None)
   >>> print(line)  # doctest: +NORMALIZE_WHITESPACE
   3   2   3   0   0   0   0   0   +   query   8   0   8   target   20   7   15   1   8,   0,   7,

The same arguments can be used when writing alignments to an output file
in PSL format using ``Bio.Align.write``. This function has an additional
keyword ``header`` (``True`` by default) specifying if the PSL header
should be written.

In addition to the ``format`` method, you can use Python’s built-in
``format`` function:

.. cont-doctest

.. code:: pycon

   >>> print(format(alignment, "psl"))  # doctest: +NORMALIZE_WHITESPACE
   6   1   0   1   0   0   0   0   +   query   8   0   8   target   20   7   15   1   8,   0,   7,

allowing ``Alignment`` objects to be used in formatted (f-) strings in
Python:

.. code:: pycon

   >>> line = f"The alignment in PSL format is '{alignment:psl}'."
   >>> print(line)  # doctest: +NORMALIZE_WHITESPACE
   The alignment in PSL format is '6   1   0   1   0   0   0   0   +   query   8   0   8   target   20   7   15   1   8,   0,   7,
   '

Note that optional keyword arguments cannot be used with the ``format``
function or with formatted strings.

.. _`subsec:align_bigpsl`:

bigPsl
~~~~~~

A bigPsl file is a bigBed file with a BED12+13 format consisting of the
12 predefined BED fields and 13 custom fields defined in the AutoSql
file `bigPsl.as <https://genome.ucsc.edu/goldenPath/help/bigPsl.html>`__
provided by UCSC, creating an indexed binary version of a PSL file (see
section :ref:`subsec:align_psl`). To create a bigPsl file, you
can either use the ``pslToBigPsl`` and ``bedToBigBed`` programs from
UCSC. or you can use Biopython by calling the ``Bio.Align.write``
function with ``fmt="bigpsl"``. While the two methods should result in
identical bigPsl files, the UCSC tools are much faster and may be more
reliable, as it is the gold standard. As bigPsl files are bigBed files,
they come with a built-in index, allowing you to quickly search a
specific genomic region.

As an example, let’s parse the bigBed file ``dna_rna.psl.bb``, available
in the ``Tests/Blat`` subdirectory in the Biopython distribution. This
file is the bigPsl equivalent of the bigBed file ``dna_rna.bb`` (see
section :ref:`subsec:align_bigbed`) and of the PSL file
``dna_rna.psl`` (see section :ref:`subsec:align_psl`).

.. doctest ../Tests/Blat lib:numpy

.. code:: pycon

   >>> from Bio import Align
   >>> alignments = Align.parse("dna_rna.psl.bb", "bigpsl")
   >>> len(alignments)
   4
   >>> print(alignments.declaration)  # doctest: +NORMALIZE_WHITESPACE
   table bigPsl
   "bigPsl pairwise alignment"
   (
      string          chrom;           "Reference sequence chromosome or scaffold"
      uint            chromStart;      "Start position in chromosome"
      uint            chromEnd;        "End position in chromosome"
      string          name;            "Name or ID of item, ideally both human readable and unique"
      uint            score;           "Score (0-1000)"
      char[1]         strand;          "+ or - indicates whether the query aligns to the + or - strand on the reference"
      uint            thickStart;      "Start of where display should be thick (start codon)"
      uint            thickEnd;        "End of where display should be thick (stop codon)"
      uint            reserved;        "RGB value (use R,G,B string in input file)"
      int             blockCount;      "Number of blocks"
      int[blockCount] blockSizes;      "Comma separated list of block sizes"
      int[blockCount] chromStarts;     "Start positions relative to chromStart"
      uint            oChromStart;     "Start position in other chromosome"
      uint            oChromEnd;       "End position in other chromosome"
      char[1]         oStrand;         "+ or -, - means that psl was reversed into BED-compatible coordinates"
      uint            oChromSize;      "Size of other chromosome."
      int[blockCount] oChromStarts;    "Start positions relative to oChromStart or from oChromStart+oChromSize depending on strand"
      lstring         oSequence;       "Sequence on other chrom (or edit list, or empty)"
      string          oCDS;            "CDS in NCBI format"
      uint            chromSize;       "Size of target chromosome"
      uint            match;           "Number of bases matched."
      uint            misMatch;        "Number of bases that don't match"
      uint            repMatch;        "Number of bases that match but are part of repeats"
      uint            nCount;          "Number of 'N' bases"
      uint            seqType;         "0=empty, 1=nucleotide, 2=amino_acid"
   )
   <BLANKLINE>

The declaration contains the specification of the columns as defined by
the ``bigPsl.as`` AutoSql file from UCSC. Target sequences (typically,
the chromosomes against which the sequences were aligned) are stored in
the ``targets`` attribute. In the bigBed format, only the identifier and
the size of each target is stored. In this example, there is only a
single chromosome:

.. cont-doctest

.. code:: pycon

   >>> alignments.targets
   [SeqRecord(seq=Seq(None, length=198295559), id='chr3', name='<unknown name>', description='<unknown description>', dbxrefs=[])]

Iterating over the alignments gives one Alignment object for each line:

.. cont-doctest

.. code:: pycon

   >>> for alignment in alignments:
   ...     print(alignment.target.id, alignment.query.id)
   ...
   chr3 NR_046654.1
   chr3 NR_046654.1_modified
   chr3 NR_111921.1
   chr3 NR_111921.1_modified

Let’s look at the individual alignments. The alignment information is
stored in the same way as for the corresponding PSL file (see
section :ref:`subsec:align_psl`):

.. cont-doctest

.. code:: pycon

   >>> alignment.coordinates
   array([[48663767, 48663795, 48663796, 48663813, 48665640, 48665716,
           48665716, 48665722, 48669098, 48669174],
          [       3,       31,       31,       48,       48,      124,
                126,      132,      132,      208]])
   >>> alignment.thickStart
   48663767
   >>> alignment.thickEnd
   48669174
   >>> alignment.matches
   162
   >>> alignment.misMatches
   2
   >>> alignment.repMatches
   39
   >>> alignment.nCount
   0

We can print the alignment in BED or PSL format:

.. cont-doctest

.. code:: pycon

   >>> print(format(alignment, "bed"))  # doctest: +NORMALIZE_WHITESPACE
   chr3    48663767    48669174    NR_111921.1_modified    1000    +   48663767    48669174    0   5   28,17,76,6,76,  0,29,1873,1949,5331,
   <BLANKLINE>
   >>> print(format(alignment, "psl"))  # doctest: +NORMALIZE_WHITESPACE
   162 2   39  0   1   2   3   5204    +   NR_111921.1_modified    220 3   208 chr3    198295559   48663767    48669174    5   28,17,76,6,76,  3,31,48,126,132,    48663767,48663796,48665640,48665716,48669098,
   <BLANKLINE>

As a bigPsl file is a special case of a bigBed file, you can use the
``search`` method on the alignments object to find alignments to
specific genomic regions. For example, we can look for regions on chr3
between positions 48000000 and 49000000:

.. cont-doctest

.. code:: pycon

   >>> selected_alignments = alignments.search("chr3", 48000000, 49000000)
   >>> for alignment in selected_alignments:
   ...     print(alignment.query.id)
   ...
   NR_111921.1
   NR_111921.1_modified

The chromosome name may be ``None`` to include all chromosomes, and the
start and end positions may be ``None`` to start searching from position
0 or to continue searching until the end of the chromosome,
respectively.

To write a bigPsl file with Biopython, use
``Bio.Align.write(alignments, "myfilename.bb", fmt="bigpsl")``, where
``myfilename.bb`` is the name of the output bigPsl file. Alternatively,
you can use a (binary) stream for output. Additional options are

-  ``compress``: If ``True`` (default), compress data using zlib; if
   ``False``, do not compress data.

-  ``extraIndex``: List of strings with the names of extra columns to be
   indexed.

-  ``cds``: If ``True``, look for a query feature of type CDS and write
   it in NCBI style in the PSL file (default: ``False``).

-  ``fa``: If ``True``, include the query sequence in the PSL file
   (default: ``False``).

-  ``mask``: Specify if repeat regions in the target sequence are masked
   and should be reported in the ``repMatches`` field instead of in the
   ``matches`` field. Acceptable values are

   -  ``None``: no masking (default);

   -  ``"lower"``: masking by lower-case characters;

   -  ``"upper"``: masking by upper-case characters.

-  ``wildcard``: Report alignments to the wildcard character
   (representing unknown nucleotides) in the target or query sequence in
   the ``nCount`` field instead of in the ``matches``, ``misMatches``,
   or ``repMatches`` fields. Default value is ``"N"``.

See section :ref:`subsec:align_psl` for an explanation on how the
number of matches, mismatches, repeat region matches, and matches to
unknown nucleotides are obtained.

Further optional arguments are ``blockSize`` (default value is 256), and
``itemsPerSlot`` (default value is 512). See the documentation of UCSC's
``bedToBigBed`` program for a description of these arguments.  Searching a
``bigPsl`` file can be faster by using ``compress=False`` and
``itemsPerSlot=1`` when creating the bigPsl file.

.. _`subsec:align_maf`:

Multiple Alignment Format (MAF)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MAF (Multiple Alignment Format) files store a series of multiple
sequence alignments in a human-readable format. MAF files are typically
used to store alignment s of genomes to each other. The file
``ucsc_test.maf`` in the ``Tests/MAF`` subdirectory of the Biopython
distribution is an example of a simple MAF file:

.. code:: text

   track name=euArc visibility=pack mafDot=off frames="multiz28wayFrames" speciesOrder="hg16 panTro1 baboon mm4 rn3" description="A sample alignment"
   ##maf version=1 scoring=tba.v8
   # tba.v8 (((human chimp) baboon) (mouse rat))
   # multiz.v7
   # maf_project.v5 _tba_right.maf3 mouse _tba_C
   # single_cov2.v4 single_cov2 /dev/stdin

   a score=23262.0
   s hg16.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
   s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
   s baboon         116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
   s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
   s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG

   a score=5062.0
   s hg16.chr7    27699739 6 + 158545518 TAAAGA
   s panTro1.chr6 28862317 6 + 161576975 TAAAGA
   s baboon         241163 6 +   4622798 TAAAGA
   s mm4.chr6     53303881 6 + 151104725 TAAAGA
   s rn3.chr4     81444246 6 + 187371129 taagga

   a score=6636.0
   s hg16.chr7    27707221 13 + 158545518 gcagctgaaaaca
   s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
   s baboon         249182 13 +   4622798 gcagctgaaaaca
   s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA

To parse this file, use

.. doctest ../Tests/MAF lib:numpy

.. code:: pycon

   >>> from Bio import Align
   >>> alignments = Align.parse("ucsc_test.maf", "maf")

Information shown in the file header (the track line and subsequent
lines starting with "``#``")) is stored in the ``metadata`` attribute of
the ``alignments`` object:

.. cont-doctest

.. code:: pycon

   >>> alignments.metadata  # doctest: +NORMALIZE_WHITESPACE
   {'name': 'euArc',
    'visibility': 'pack',
    'mafDot': 'off',
    'frames': 'multiz28wayFrames',
    'speciesOrder': ['hg16', 'panTro1', 'baboon', 'mm4', 'rn3'],
    'description': 'A sample alignment',
    'MAF Version': '1',
    'Scoring': 'tba.v8',
    'Comments': ['tba.v8 (((human chimp) baboon) (mouse rat))',
                 'multiz.v7',
                 'maf_project.v5 _tba_right.maf3 mouse _tba_C',
                 'single_cov2.v4 single_cov2 /dev/stdin']}

By iterating over the ``alignments`` we obtain one ``Alignment`` object
for each alignment block in the MAF file:

.. cont-doctest

.. code:: pycon

   >>> alignment = next(alignments)
   >>> alignment.score
   23262.0
   >>> {seq.id: len(seq) for seq in alignment.sequences}  # doctest: +NORMALIZE_WHITESPACE
   {'hg16.chr7': 158545518,
    'panTro1.chr6': 161576975,
    'baboon': 4622798,
    'mm4.chr6': 151104725,
    'rn3.chr4': 187371129}
   >>> print(alignment.coordinates)
   [[27578828 27578829 27578831 27578831 27578850 27578850 27578866]
    [28741140 28741141 28741143 28741143 28741162 28741162 28741178]
    [  116834   116835   116837   116837   116856   116856   116872]
    [53215344 53215344 53215346 53215347 53215366 53215366 53215382]
    [81344243 81344243 81344245 81344245 81344264 81344267 81344283]]
   >>> print(alignment)
   hg16.chr7  27578828 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG 27578866
   panTro1.c  28741140 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG 28741178
   baboon       116834 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG   116872
   mm4.chr6   53215344 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG 53215382
   rn3.chr4   81344243 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG 81344283
   <BLANKLINE>
   >>> print(format(alignment, "phylip"))
   5 42
   hg16.chr7 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
   panTro1.chAAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
   baboon    AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
   mm4.chr6  -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
   rn3.chr4  -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
   <BLANKLINE>

In addition to the "``a``" (alignment block) and "``s``" (sequence)
lines, MAF files may contain "``i``" lines with information about the
genome sequence before and after this block, "``e``" lines with
information about empty parts of the alignment, and "``q``" lines
showing the quality of each aligned base. This is an example of an
alignment block including such lines:

.. code:: text

   a score=19159.000000
   s mm9.chr10                         3014644 45 + 129993255 CCTGTACC---CTTTGGTGAGAATTTTTGTTTCAGTGTTAAAAGTTTG
   s hg18.chr6                        15870786 46 - 170899992 CCTATACCTTTCTTTTATGAGAA-TTTTGTTTTAATCCTAAAC-TTTT
   i hg18.chr6                        I 9085 C 0
   s panTro2.chr6                     16389355 46 - 173908612 CCTATACCTTTCTTTTATGAGAA-TTTTGTTTTAATCCTAAAC-TTTT
   q panTro2.chr6                                             99999999999999999999999-9999999999999999999-9999
   i panTro2.chr6                     I 9106 C 0
   s calJac1.Contig6394                   6182 46 +    133105 CCTATACCTTTCTTTCATGAGAA-TTTTGTTTGAATCCTAAAC-TTTT
   i calJac1.Contig6394               N 0 C 0
   s loxAfr1.scaffold_75566               1167 34 -     10574 ------------TTTGGTTAGAA-TTATGCTTTAATTCAAAAC-TTCC
   q loxAfr1.scaffold_75566                                   ------------99999699899-9999999999999869998-9997
   i loxAfr1.scaffold_75566           N 0 C 0
   e tupBel1.scaffold_114895.1-498454   167376 4145 -    498454 I
   e echTel1.scaffold_288249             87661 7564 +    100002 I
   e otoGar1.scaffold_334.1-359464      181217 2931 -    359464 I
   e ponAbe2.chr6                     16161448 8044 - 174210431 I

This is the 10th alignment block in the file ``ucsc_mm9_chr10.maf``
(available in the ``Tests/MAF`` subdirectory of the Biopython
distribution):

.. doctest ../Tests/MAF lib:numpy

.. code:: pycon

   >>> from Bio import Align
   >>> alignments = Align.parse("ucsc_mm9_chr10.maf", "maf")
   >>> for i in range(10):
   ...     alignment = next(alignments)
   ...
   >>> alignment.score
   19159.0
   >>> print(alignment)
   mm9.chr10   3014644 CCTGTACC---CTTTGGTGAGAATTTTTGTTTCAGTGTTAAAAGTTTG   3014689
   hg18.chr6 155029206 CCTATACCTTTCTTTTATGAGAA-TTTTGTTTTAATCCTAAAC-TTTT 155029160
   panTro2.c 157519257 CCTATACCTTTCTTTTATGAGAA-TTTTGTTTTAATCCTAAAC-TTTT 157519211
   calJac1.C      6182 CCTATACCTTTCTTTCATGAGAA-TTTTGTTTGAATCCTAAAC-TTTT      6228
   loxAfr1.s      9407 ------------TTTGGTTAGAA-TTATGCTTTAATTCAAAAC-TTCC      9373
   <BLANKLINE>

The "``i``" lines show the relationship between the sequence in the
current alignment block to the ones in the preceding and subsequent
alignment block. This information is stored in the ``annotations``
attribute of the corresponding sequence:

.. cont-doctest

.. code:: pycon

   >>> alignment.sequences[0].annotations
   {}
   >>> alignment.sequences[1].annotations
   {'leftStatus': 'I', 'leftCount': 9085, 'rightStatus': 'C', 'rightCount': 0}

showing that there are 9085 bases inserted ("``I``") between this block
and the preceding one, while the block is contiguous ("``C``") with the
subsequent one. See the `UCSC
documentation <https://genome.ucsc.edu/FAQ/FAQformat.html#format5>`__
for the full description of these fields and status characters.

The "``q``" lines show the sequence quality, which is stored under the
"``quality``" dictionary key of the\ ``annotations`` attribute of the
corresponding sequence:

.. cont-doctest

.. code:: pycon

   >>> alignment.sequences[2].annotations["quality"]
   '9999999999999999999999999999999999999999999999'
   >>> alignment.sequences[4].annotations["quality"]
   '9999969989999999999999998699989997'

The "``e``" lines show information about species with a contiguous
sequence before and after this alignment bloack, but with no aligning
nucleotides in this alignment block. This is stored under the
"``empty``" key of the ``alignment.annotations`` dictionary:

.. cont-doctest

.. code:: pycon

   >>> alignment.annotations["empty"]  # doctest: +NORMALIZE_WHITESPACE
   [(SeqRecord(seq=Seq(None, length=498454), id='tupBel1.scaffold_114895.1-498454', name='', description='', dbxrefs=[]), (331078, 326933), 'I'),
    (SeqRecord(seq=Seq(None, length=100002), id='echTel1.scaffold_288249', name='', description='', dbxrefs=[]), (87661, 95225), 'I'),
    (SeqRecord(seq=Seq(None, length=359464), id='otoGar1.scaffold_334.1-359464', name='', description='', dbxrefs=[]), (178247, 175316), 'I'),
    (SeqRecord(seq=Seq(None, length=174210431), id='ponAbe2.chr6', name='', description='', dbxrefs=[]), (158048983, 158040939), 'I')]

This shows for example that there were non-aligning bases inserted
("``I``") from position 158040939 to 158048983 on the opposite strand of
the ``ponAbe2.chr6`` genomic sequence. Again, see the `UCSC
documentation <https://genome.ucsc.edu/FAQ/FAQformat.html#format5>`__
for the full definition of "``e``" lines.

To print an alignment in MAF format, you can use Python’s built-in
``format`` function:

.. cont-doctest

.. code:: pycon

   >>> print(format(alignment, "MAF"))
   a score=19159.000000
   s mm9.chr10                         3014644   45 + 129993255 CCTGTACC---CTTTGGTGAGAATTTTTGTTTCAGTGTTAAAAGTTTG
   s hg18.chr6                        15870786   46 - 170899992 CCTATACCTTTCTTTTATGAGAA-TTTTGTTTTAATCCTAAAC-TTTT
   i hg18.chr6                        I 9085 C 0
   s panTro2.chr6                     16389355   46 - 173908612 CCTATACCTTTCTTTTATGAGAA-TTTTGTTTTAATCCTAAAC-TTTT
   q panTro2.chr6                                               99999999999999999999999-9999999999999999999-9999
   i panTro2.chr6                     I 9106 C 0
   s calJac1.Contig6394                   6182   46 +    133105 CCTATACCTTTCTTTCATGAGAA-TTTTGTTTGAATCCTAAAC-TTTT
   i calJac1.Contig6394               N 0 C 0
   s loxAfr1.scaffold_75566               1167   34 -     10574 ------------TTTGGTTAGAA-TTATGCTTTAATTCAAAAC-TTCC
   q loxAfr1.scaffold_75566                                     ------------99999699899-9999999999999869998-9997
   i loxAfr1.scaffold_75566           N 0 C 0
   e tupBel1.scaffold_114895.1-498454   167376 4145 -    498454 I
   e echTel1.scaffold_288249             87661 7564 +    100002 I
   e otoGar1.scaffold_334.1-359464      181217 2931 -    359464 I
   e ponAbe2.chr6                     16161448 8044 - 174210431 I
   <BLANKLINE>
   <BLANKLINE>

To write a complete MAF file, use
``Bio.Align.write(alignments, "myfilename.maf", fmt="maf")``, where
``myfilename.maf`` is the name of the output MAF file. Alternatively,
you can use a (text) stream for output. File header information will be
taken from the ``metadata`` attribute of the ``alignments`` object. If
you are creating the alignments from scratch, you can use the
``Alignments`` (plural) class to create a list-like ``alignments``
object (see Section :ref:`sec:alignments`) and give it a
``metadata`` attribute.

.. _`subsec:align_bigmaf`:

bigMaf
~~~~~~

A bigMaf file is a bigBed file with a BED3+1 format consisting of the 3
required BED fields plus a custom field that stores a MAF alignment
block as a string, creating an indexed binary version of a MAF file (see
section :ref:`subsec:align_maf`). The associated AutoSql file
`bigMaf.as <https://genome.ucsc.edu/goldenPath/help/examples/bigMaf.as>`__
is provided by UCSC. To create a bigMaf file, you can either use the
``mafToBigMaf`` and ``bedToBigBed`` programs from UCSC. or you can use
Biopython by calling the Bio.Align.write function with ``fmt="bigmaf"``.
While the two methods should result in identical bigMaf files, the UCSC
tools are much faster and may be more reliable, as it is the gold
standard. As bigMaf files are bigBed files, they come with a built-in
index, allowing you to quickly search a specific region of the reference
genome.

The file ``ucsc_test.bb`` in the ``Tests/MAF`` subdirectory of the
Biopython distribution is an example of a bigMaf file. This file is
equivalent to the MAF file ``ucsc_test.maf`` (see
section :ref:`subsec:align_maf`). To parse this file, use

.. doctest ../Tests/MAF lib:numpy

.. code:: pycon

   >>> from Bio import Align
   >>> alignments = Align.parse("ucsc_test.bb", "bigmaf")
   >>> len(alignments)
   3
   >>> print(alignments.declaration)  # doctest: +NORMALIZE_WHITESPACE
   table bedMaf
   "Bed3 with MAF block"
   (
      string  chrom;         "Reference sequence chromosome or scaffold"
      uint    chromStart;    "Start position in chromosome"
      uint    chromEnd;      "End position in chromosome"
      lstring mafBlock;      "MAF block"
   )
   <BLANKLINE>

The declaration contains the specification of the columns as defined by
the bigMaf.as AutoSql file from UCSC.

The bigMaf file does not store the header information found in the MAF
file, but it does define a reference genome. The corresponding
``SeqRecord`` is stored in the ``targets`` attribute of the
``alignments`` object:

.. cont-doctest

.. code:: pycon

   >>> alignments.reference
   'hg16'
   >>> alignments.targets  # doctest: +ELLIPSIS
   [SeqRecord(seq=Seq(None, length=158545518), id='hg16.chr7', ...)]

By iterating over the ``alignments`` we obtain one ``Alignment`` object
for each alignment block in the bigMaf file:

.. cont-doctest

.. code:: pycon

   >>> alignment = next(alignments)
   >>> alignment.score
   23262.0
   >>> {seq.id: len(seq) for seq in alignment.sequences}  # doctest: +NORMALIZE_WHITESPACE
   {'hg16.chr7': 158545518,
    'panTro1.chr6': 161576975,
    'baboon': 4622798,
    'mm4.chr6': 151104725,
    'rn3.chr4': 187371129}
   >>> print(alignment.coordinates)
   [[27578828 27578829 27578831 27578831 27578850 27578850 27578866]
    [28741140 28741141 28741143 28741143 28741162 28741162 28741178]
    [  116834   116835   116837   116837   116856   116856   116872]
    [53215344 53215344 53215346 53215347 53215366 53215366 53215382]
    [81344243 81344243 81344245 81344245 81344264 81344267 81344283]]
   >>> print(alignment)
   hg16.chr7  27578828 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG 27578866
   panTro1.c  28741140 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG 28741178
   baboon       116834 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG   116872
   mm4.chr6   53215344 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG 53215382
   rn3.chr4   81344243 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG 81344283
   <BLANKLINE>
   >>> print(format(alignment, "phylip"))
   5 42
   hg16.chr7 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
   panTro1.chAAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
   baboon    AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
   mm4.chr6  -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
   rn3.chr4  -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
   <BLANKLINE>

Information in the "``i``", "``e``", and "``q``" lines is stored in the
same way as in the corresponding MAF file (see
section :ref:`subsec:align_maf`):

.. doctest ../Tests/MAF lib:numpy

.. code:: pycon

   >>> from Bio import Align
   >>> alignments = Align.parse("ucsc_mm9_chr10.bb", "bigmaf")
   >>> for i in range(10):
   ...     alignment = next(alignments)
   ...
   >>> alignment.score
   19159.0
   >>> print(alignment)
   mm9.chr10   3014644 CCTGTACC---CTTTGGTGAGAATTTTTGTTTCAGTGTTAAAAGTTTG   3014689
   hg18.chr6 155029206 CCTATACCTTTCTTTTATGAGAA-TTTTGTTTTAATCCTAAAC-TTTT 155029160
   panTro2.c 157519257 CCTATACCTTTCTTTTATGAGAA-TTTTGTTTTAATCCTAAAC-TTTT 157519211
   calJac1.C      6182 CCTATACCTTTCTTTCATGAGAA-TTTTGTTTGAATCCTAAAC-TTTT      6228
   loxAfr1.s      9407 ------------TTTGGTTAGAA-TTATGCTTTAATTCAAAAC-TTCC      9373
   <BLANKLINE>
   >>> print(format(alignment, "MAF"))
   a score=19159.000000
   s mm9.chr10                         3014644   45 + 129993255 CCTGTACC---CTTTGGTGAGAATTTTTGTTTCAGTGTTAAAAGTTTG
   s hg18.chr6                        15870786   46 - 170899992 CCTATACCTTTCTTTTATGAGAA-TTTTGTTTTAATCCTAAAC-TTTT
   i hg18.chr6                        I 9085 C 0
   s panTro2.chr6                     16389355   46 - 173908612 CCTATACCTTTCTTTTATGAGAA-TTTTGTTTTAATCCTAAAC-TTTT
   q panTro2.chr6                                               99999999999999999999999-9999999999999999999-9999
   i panTro2.chr6                     I 9106 C 0
   s calJac1.Contig6394                   6182   46 +    133105 CCTATACCTTTCTTTCATGAGAA-TTTTGTTTGAATCCTAAAC-TTTT
   i calJac1.Contig6394               N 0 C 0
   s loxAfr1.scaffold_75566               1167   34 -     10574 ------------TTTGGTTAGAA-TTATGCTTTAATTCAAAAC-TTCC
   q loxAfr1.scaffold_75566                                     ------------99999699899-9999999999999869998-9997
   i loxAfr1.scaffold_75566           N 0 C 0
   e tupBel1.scaffold_114895.1-498454   167376 4145 -    498454 I
   e echTel1.scaffold_288249             87661 7564 +    100002 I
   e otoGar1.scaffold_334.1-359464      181217 2931 -    359464 I
   e ponAbe2.chr6                     16161448 8044 - 174210431 I
   <BLANKLINE>
   <BLANKLINE>
   >>> alignment.sequences[1].annotations
   {'leftStatus': 'I', 'leftCount': 9085, 'rightStatus': 'C', 'rightCount': 0}
   >>> alignment.sequences[2].annotations["quality"]
   '9999999999999999999999999999999999999999999999'
   >>> alignment.sequences[4].annotations["quality"]
   '9999969989999999999999998699989997'
   >>> alignment.annotations["empty"]  # doctest: +NORMALIZE_WHITESPACE
   [(SeqRecord(seq=Seq(None, length=498454), id='tupBel1.scaffold_114895.1-498454', name='', description='', dbxrefs=[]), (331078, 326933), 'I'),
    (SeqRecord(seq=Seq(None, length=100002), id='echTel1.scaffold_288249', name='', description='', dbxrefs=[]), (87661, 95225), 'I'),
    (SeqRecord(seq=Seq(None, length=359464), id='otoGar1.scaffold_334.1-359464', name='', description='', dbxrefs=[]), (178247, 175316), 'I'),
    (SeqRecord(seq=Seq(None, length=174210431), id='ponAbe2.chr6', name='', description='', dbxrefs=[]), (158048983, 158040939), 'I')]

To write a complete bigMaf file, use
``Bio.Align.write(alignments, "myfilename.bb", fmt="bigMaf")``, where
``myfilename.bb`` is the name of the output bigMaf file. Alternatively,
you can use a (binary) stream for output. If you are creating the
alignments from scratch, you can use the ``Alignments`` (plural) class
to create a list-like ``alignments`` object (see
Section :ref:`sec:alignments`) and give it a ``targets`` attribute.
The latter must be a list of ``SeqRecord`` objects for the chromosomes
for the reference species in the order in which they appear in the
alignments. Alternatively, you can use the ``targets`` keyword argument
when calling ``Bio.Align.write``. The ``id`` of each ``SeqRecord`` must
be of the form ``reference.chromosome``, where ``reference`` refers to
the reference species. ``Bio.Align.write`` has the additional keyword
argument ``compress`` (``True`` by default) specifying whether the data
should be compressed using zlib.
Further optional arguments are ``blockSize`` (default value is 256), and
``itemsPerSlot`` (default value is 512). See the documentation of UCSC's
``bedToBigBed`` program for a description of these arguments.

As a bigMaf file is a special case of a bigBed file, you can use the
``search`` method on the ``alignments`` object to find alignments to
specific regions of the reference species. For example, we can look for
regions on chr10 between positions 3018000 and 3019000 on chromosome 10:

.. cont-doctest

.. code:: pycon

   >>> selected_alignments = alignments.search("mm9.chr10", 3018000, 3019000)
   >>> for alignment in selected_alignments:
   ...     start, end = alignment.coordinates[0, 0], alignment.coordinates[0, -1]
   ...     print(start, end)
   ...
   3017743 3018161
   3018161 3018230
   3018230 3018359
   3018359 3018482
   3018482 3018644
   3018644 3018822
   3018822 3018932
   3018932 3019271

The chromosome name may be ``None`` to include all chromosomes, and the
start and end positions may be ``None`` to start searching from position
0 or to continue searching until the end of the chromosome,
respectively. Note that we can search on genomic position for the
reference species only.

Searching a ``bigMaf`` file can be faster by using ``compress=False`` and
``itemsPerSlot=1`` when creating the bigMaf file.

.. _`subsec:align_chain`:

UCSC chain file format
~~~~~~~~~~~~~~~~~~~~~~

Chain files describe a pairwise alignment between two nucleotide
sequences, allowing gaps in both sequences. Only the length of each
aligned subsequences and the gap lengths are stored in a chain file; the
sequences themselves are not stored. Chain files are typically used to
store alignments between two genome assembly versions, allowing
alignments to one genome assembly version to be lifted over to the other
genome assembly. This is an example of a chain file (available as
``psl_34_001.chain`` in the ``Tests/Blat`` subdirectory of the Biopython
distribution):

.. code:: text

   chain 16 chr4 191154276 + 61646095 61646111 hg18_dna 33 + 11 27 1
   16
   chain 33 chr1 249250621 + 10271783 10271816 hg18_dna 33 + 0 33 2
   33
   chain 17 chr2 243199373 + 53575980 53575997 hg18_dna 33 - 8 25 3
   17
   chain 35 chr9 141213431 + 85737865 85737906 hg19_dna 50 + 9 50 4
   41
   chain 41 chr8 146364022 + 95160479 95160520 hg19_dna 50 + 8 49 5
   41
   chain 30 chr22 51304566 + 42144400 42144436 hg19_dna 50 + 11 47 6
   36
   chain 41 chr2 243199373 + 183925984 183926028 hg19_dna 50 + 1 49 7
   6       0       4
   38
   chain 31 chr19 59128983 + 35483340 35483510 hg19_dna 50 + 10 46 8
   25      134     0
   11
   chain 39 chr18 78077248 + 23891310 23891349 hg19_dna 50 + 10 49 9
   39
   ...

This file was generated by running UCSC’s ``pslToChain`` program on the
PSL file ``psl_34_001.psl``. According to the chain file format
specification, there should be a blank line after each chain block, but
some tools (including ``pslToChain``) apparently do not follow this
rule.

To parse this file, use

.. doctest ../Tests/Blat lib:numpy

.. code:: pycon

   >>> from Bio import Align
   >>> alignments = Align.parse("psl_34_001.chain", "chain")

Iterate over alignments to get one ``Alignment`` object for each chain:

.. cont-doctest

.. code:: pycon

   >>> for alignment in alignments:
   ...     print(alignment.target.id, alignment.query.id)  # doctest: +ELLIPSIS
   ...
   chr4 hg18_dna
   chr1 hg18_dna
   chr2 hg18_dna
   chr9 hg19_dna
   chr8 hg19_dna
   chr22 hg19_dna
   chr2 hg19_dna
   ...
   chr1 hg19_dna

Iterate from the start until we reach the seventh alignment:

.. cont-doctest

.. code:: pycon

   >>> alignments = iter(alignments)
   >>> for i in range(7):
   ...     alignment = next(alignments)
   ...

Check the alignment score and chain ID (the first and last number,
respectively, in the header line of each chain block) to confirm that we
got the seventh alignment:

.. cont-doctest

.. code:: pycon

   >>> alignment.score
   41.0
   >>> alignment.annotations["id"]
   '7'

We can print the alignment in the chain file format. The alignment
coordinates are consistent with the information in the chain block, with
an aligned section of 6 nucleotides, a gap of 4 nucleotides, and an
aligned section of 38 nucleotides:

.. cont-doctest

.. code:: pycon

   >>> print(format(alignment, "chain"))  # doctest: +NORMALIZE_WHITESPACE
   chain 41 chr2 243199373 + 183925984 183926028 hg19_dna 50 + 1 49 7
   6   0   4
   38
   <BLANKLINE>
   <BLANKLINE>
   >>> alignment.coordinates
   array([[183925984, 183925990, 183925990, 183926028],
          [        1,         7,        11,        49]])
   >>> print(alignment)
   chr2      183925984 ??????----?????????????????????????????????????? 183926028
                     0 ||||||----||||||||||||||||||||||||||||||||||||||        48
   hg19_dna          1 ????????????????????????????????????????????????        49
   <BLANKLINE>

We can also print the alignment in a few other alignment fite formats:

.. cont-doctest

.. code:: pycon

   >>> print(format(alignment, "BED"))  # doctest: +NORMALIZE_WHITESPACE
   chr2    183925984   183926028   hg19_dna    41  +   183925984   183926028   0   2   6,38,   0,6,
   <BLANKLINE>
   >>> print(format(alignment, "PSL"))  # doctest: +NORMALIZE_WHITESPACE
   44  0   0   0   1   4   0   0   +   hg19_dna    50  1   49  chr2    243199373   183925984   183926028   2   6,38,   1,11,   183925984,183925990,
   <BLANKLINE>
   >>> print(format(alignment, "exonerate"))
   vulgar: hg19_dna 1 49 + chr2 183925984 183926028 + 41 M 6 6 G 4 0 M 38 38
   <BLANKLINE>
   >>> print(alignment.format("exonerate", "cigar"))
   cigar: hg19_dna 1 49 + chr2 183925984 183926028 + 41 M 6 I 4 M 38
   <BLANKLINE>
   >>> print(format(alignment, "sam"))  # doctest: +NORMALIZE_WHITESPACE
   hg19_dna    0   chr2    183925985   255 1S6M4I38M1S *   0   0   *   *   AS:i:41 id:A:7
   <BLANKLINE>
