.. _`chapter:pairwise`:

Pairwise sequence alignment
===========================

Pairwise sequence alignment is the process of aligning two sequences to
each other by optimizing the similarity score between them. The
``Bio.Align`` module contains the ``PairwiseAligner`` class for global
and local alignments using the Needleman-Wunsch, Smith-Waterman, Gotoh
(three-state), and Waterman-Smith-Beyer global and local pairwise
alignment algorithms, and the Fast Optimal Global Alignment Algorithm (FOGSAA),
with numerous options to change the alignment parameters. We refer to Durbin
*et al.* [Durbin1998]_ for in-depth information on sequence alignment
algorithms.

.. _`sec:pairwise-basic`:

Basic usage
-----------

To generate pairwise alignments, first create a ``PairwiseAligner``
object:

.. doctest examples

.. code:: pycon

   >>> from Bio import Align
   >>> aligner = Align.PairwiseAligner()

The ``PairwiseAligner`` object ``aligner`` (see
Section :ref:`sec:pairwise-aligner`) stores the alignment parameters
to be used for the pairwise alignments. These attributes can be set in
the constructor of the object:

.. cont-doctest

.. code:: pycon

   >>> aligner = Align.PairwiseAligner(match_score=1.0)

or after the object is made:

.. cont-doctest

.. code:: pycon

   >>> aligner.match_score = 1.0

Use the ``aligner.score`` method to calculate the alignment score
between two sequences:

.. cont-doctest

.. code:: pycon

   >>> target = "GAACT"
   >>> query = "GAT"
   >>> score = aligner.score(target, query)
   >>> score
   3.0

The ``aligner.align`` method returns ``Alignment`` objects, each
representing one alignment between the two sequences:

.. cont-doctest

.. code:: pycon

   >>> alignments = aligner.align(target, query)
   >>> alignment = alignments[0]
   >>> alignment  # doctest: +ELLIPSIS
   <Alignment object (2 rows x 5 columns) at ...>

Iterate over the ``Alignment`` objects and print them to see the
alignments:

.. cont-doctest

.. code:: pycon

   >>> for alignment in alignments:
   ...     print(alignment)
   ...
   target            0 GAACT 5
                     0 ||--| 5
   query             0 GA--T 3
   <BLANKLINE>
   target            0 GAACT 5
                     0 |-|-| 5
   query             0 G-A-T 3
   <BLANKLINE>

Use indices to get the aligned sequence (see :ref:`subsec:slicing-indexing-alignment`):

.. cont-doctest

.. code:: pycon

   >>> alignment[0]
   'GAACT'
   >>> alignment[1]
   'G-A-T'

Each alignment stores the alignment score:

.. cont-doctest

.. code:: pycon

   >>> alignment.score
   3.0

as well as pointers to the sequences that were aligned:

.. cont-doctest

.. code:: pycon

   >>> alignment.target
   'GAACT'
   >>> alignment.query
   'GAT'

Internally, the alignment is stored in terms of the sequence
coordinates:

.. cont-doctest

.. code:: pycon

   >>> alignment = alignments[0]
   >>> alignment.coordinates
   array([[0, 2, 4, 5],
          [0, 2, 2, 3]])

Here, the two rows refer to the target and query sequence. These
coordinates show that the alignment consists of the following three
blocks:

-  ``target[0:2]`` aligned to ``query[0:2]``;

-  ``target[2:4]`` aligned to a gap, since ``query[2:2]`` is an empty
   string (i.e., a deletion);

-  ``target[4:5]`` aligned to ``query[2:3]``.

The number of aligned sequences is always 2 for a pairwise alignment:

.. cont-doctest

.. code:: pycon

   >>> len(alignment)
   2

The alignment length is defined as the number of columns in the
alignment as printed. This is equal to the sum of the number of matches,
number of mismatches, and the total length of gaps in the target and
query:

.. cont-doctest

.. code:: pycon

   >>> alignment.length
   5

The ``aligned`` property, which returns the start and end indices of
aligned subsequences, returns two tuples of length 2 for the first
alignment:

.. cont-doctest

.. code:: pycon

   >>> alignment.aligned
   array([[[0, 2],
           [4, 5]],
   <BLANKLINE>
          [[0, 2],
           [2, 3]]])

while for the alternative alignment, two tuples of length 3 are
returned:

.. cont-doctest

.. code:: pycon

   >>> alignment = alignments[1]
   >>> print(alignment)
   target            0 GAACT 5
                     0 |-|-| 5
   query             0 G-A-T 3
   <BLANKLINE>
   >>> alignment.aligned
   array([[[0, 1],
           [2, 3],
           [4, 5]],
   <BLANKLINE>
          [[0, 1],
           [1, 2],
           [2, 3]]])

Note that different alignments may have the same subsequences aligned to
each other. In particular, this may occur if alignments differ from each
other in terms of their gap placement only:

.. cont-doctest

.. code:: pycon

   >>> aligner.mode = "global"
   >>> aligner.mismatch_score = -10
   >>> alignments = aligner.align("AAACAAA", "AAAGAAA")
   >>> len(alignments)
   2
   >>> print(alignments[0])
   target            0 AAAC-AAA 7
                     0 |||--||| 8
   query             0 AAA-GAAA 7
   <BLANKLINE>
   >>> alignments[0].aligned
   array([[[0, 3],
           [4, 7]],
   <BLANKLINE>
          [[0, 3],
           [4, 7]]])
   >>> print(alignments[1])
   target            0 AAA-CAAA 7
                     0 |||--||| 8
   query             0 AAAG-AAA 7
   <BLANKLINE>
   >>> alignments[1].aligned
   array([[[0, 3],
           [4, 7]],
   <BLANKLINE>
          [[0, 3],
           [4, 7]]])

The ``map`` method can be applied on a pairwise alignment ``alignment1``
to find the pairwise alignment of the query of ``alignment2`` to the
target of ``alignment1``, where the target of ``alignment2`` and the
query of ``alignment1`` are identical. A typical example is where
``alignment1`` is the pairwise alignment between a chromosome and a
transcript, ``alignment2`` is the pairwise alignment between the
transcript and a sequence (e.g., an RNA-seq read), and we want to find
the alignment of the sequence to the chromosome:

.. cont-doctest

.. code:: pycon

   >>> aligner.mode = "local"
   >>> aligner.open_gap_score = -1
   >>> aligner.extend_gap_score = 0
   >>> chromosome = "AAAAAAAACCCCCCCAAAAAAAAAAAGGGGGGAAAAAAAA"
   >>> transcript = "CCCCCCCGGGGGG"
   >>> alignments1 = aligner.align(chromosome, transcript)
   >>> len(alignments1)
   1
   >>> alignment1 = alignments1[0]
   >>> print(alignment1)
   target            8 CCCCCCCAAAAAAAAAAAGGGGGG 32
                     0 |||||||-----------|||||| 24
   query             0 CCCCCCC-----------GGGGGG 13
   <BLANKLINE>
   >>> sequence = "CCCCGGGG"
   >>> alignments2 = aligner.align(transcript, sequence)
   >>> len(alignments2)
   1
   >>> alignment2 = alignments2[0]
   >>> print(alignment2)
   target            3 CCCCGGGG 11
                     0 ||||||||  8
   query             0 CCCCGGGG  8
   <BLANKLINE>
   >>> mapped_alignment = alignment1.map(alignment2)
   >>> print(mapped_alignment)
   target           11 CCCCAAAAAAAAAAAGGGG 30
                     0 ||||-----------|||| 19
   query             0 CCCC-----------GGGG  8
   <BLANKLINE>
   >>> format(mapped_alignment, "psl")
   '8\t0\t0\t0\t0\t0\t1\t11\t+\tquery\t8\t0\t8\ttarget\t40\t11\t30\t2\t4,4,\t0,4,\t11,26,\n'

Mapping the alignment does not depend on the sequence contents. If we
delete the sequence contents, the same alignment is found in PSL format
(though we obviously lose the ability to print the sequence alignment):

.. cont-doctest

.. code:: pycon

   >>> from Bio.Seq import Seq
   >>> alignment1.target = Seq(None, len(alignment1.target))
   >>> alignment1.query = Seq(None, len(alignment1.query))
   >>> alignment2.target = Seq(None, len(alignment2.target))
   >>> alignment2.query = Seq(None, len(alignment2.query))
   >>> mapped_alignment = alignment1.map(alignment2)
   >>> format(mapped_alignment, "psl")
   '8\t0\t0\t0\t0\t0\t1\t11\t+\tquery\t8\t0\t8\ttarget\t40\t11\t30\t2\t4,4,\t0,4,\t11,26,\n'

By default, a global pairwise alignment is performed, which finds the
optimal alignment over the whole length of ``target`` and ``query``.
Instead, a local alignment will find the subsequence of ``target`` and
``query`` with the highest alignment score. Local alignments can be
generated by setting ``aligner.mode`` to ``"local"``:

.. cont-doctest

.. code:: pycon

   >>> aligner.mode = "local"
   >>> target = "AGAACTC"
   >>> query = "GAACT"
   >>> score = aligner.score(target, query)
   >>> score
   5.0
   >>> alignments = aligner.align(target, query)
   >>> for alignment in alignments:
   ...     print(alignment)
   ...
   target            1 GAACT 6
                     0 ||||| 5
   query             0 GAACT 5
   <BLANKLINE>

Note that there is some ambiguity in the definition of the best local
alignments if segments with a score 0 can be added to the alignment. We
follow the suggestion by Waterman & Eggert
[Waterman1987]_ and disallow such extensions.

If `aligner.mode` is set to `"fogsaa"`, then the Fast Optimal Global Alignment
Algorithm [Chakraborty2013]_ with some modifications is used. This mode
calculates a global alignment, but it is not like the regular `"global"` mode.
It is best suited for long alignments between similar sequences. Rather than
calculating all possible alignments like other algorithms do, FOGSAA uses a
heuristic to detect steps in an alignment that cannot lead to an optimal
alignment. This can speed up alignment, however, the heuristic makes
assumptions about your match, mismatch, and gap scores. If the match score is
less than the mismatch score or any gap score, or if any gap score is greater
than the mismatch score, then a warning is raised and the algorithm may return
incorrect results. Unlike other modes that may return more than one alignment,
FOGSAA always returns only one alignment.

.. cont-doctest

.. code:: pycon

   >>> aligner.mode = "fogsaa"
   >>> aligner.mismatch_score = -10
   >>> alignments = aligner.align("AAACAAA", "AAAGAAA")
   >>> len(alignments)
   1
   >>> print(alignments[0])
   target            0 AAAC-AAA 7
                     0 |||--||| 8
   query             0 AAA-GAAA 7
   <BLANKLINE>

.. _`sec:pairwise-aligner`:

The pairwise aligner object
---------------------------

The ``PairwiseAligner`` object stores all alignment parameters to be
used for the pairwise alignments. To see an overview of the values for
all parameters, use

.. doctest

.. code:: pycon

   >>> from Bio import Align
   >>> aligner = Align.PairwiseAligner(match_score=1.0, mode="local")
   >>> print(aligner)
   Pairwise sequence aligner with parameters
     wildcard: None
     match_score: 1.000000
     mismatch_score: 0.000000
     open_internal_insertion_score: 0.000000
     extend_internal_insertion_score: 0.000000
     open_left_insertion_score: 0.000000
     extend_left_insertion_score: 0.000000
     open_right_insertion_score: 0.000000
     extend_right_insertion_score: 0.000000
     open_internal_deletion_score: 0.000000
     extend_internal_deletion_score: 0.000000
     open_left_deletion_score: 0.000000
     extend_left_deletion_score: 0.000000
     open_right_deletion_score: 0.000000
     extend_right_deletion_score: 0.000000
     mode: local
   <BLANKLINE>

See Sections :ref:`sec:pairwise-substitution-scores`,
:ref:`sec:pairwise-affine-gapscores`, and
:ref:`sec:pairwise-general-gapscores` below for the definition of
these parameters. The attribute ``mode`` (described above in
Section :ref:`sec:pairwise-basic`) can be set equal to ``"global"``
or ``"local"`` to specify global or local pairwise alignment,
respectively.

Depending on the gap scoring parameters (see
Sections :ref:`sec:pairwise-affine-gapscores` and
:ref:`sec:pairwise-general-gapscores`) and mode, a
``PairwiseAligner`` object automatically chooses the appropriate
algorithm to use for pairwise sequence alignment. To verify the selected
algorithm, use

.. cont-doctest

.. code:: pycon

   >>> aligner.algorithm
   'Smith-Waterman'

This attribute is read-only.

A ``PairwiseAligner`` object also stores the precision :math:`\epsilon`
to be used during alignment. The value of :math:`\epsilon` is stored in
the attribute ``aligner.epsilon``, and by default is equal to
:math:`10^{-6}`:

.. cont-doctest

.. code:: pycon

   >>> aligner.epsilon
   1e-06

Two scores will be considered equal to each other for the purpose of the
alignment if the absolute difference between them is less than
:math:`\epsilon`.

.. _`sec:pairwise-substitution-scores`:

Substitution scores
-------------------

Substitution scores define the value to be added to the total score when
two letters (nucleotides or amino acids) are aligned to each other. The
substitution scores to be used by the ``PairwiseAligner`` can be
specified in two ways:

-  By specifying a match score for identical letters, and a mismatch
   scores for mismatched letters. Nucleotide sequence alignments are
   typically based on match and mismatch scores. For example, by default
   BLAST [Altschul1990]_ uses a match score of
   :math:`+1` and a mismatch score of :math:`-2` for nucleotide
   alignments by ``megablast``, with a gap penalty of 2.5 (see section
   :ref:`sec:pairwise-affine-gapscores` for more information on gap
   scores). Match and mismatch scores can be specified by setting the
   ``match`` and ``mismatch`` attributes of the ``PairwiseAligner``
   object:

   .. doctest examples lib:numpy

   .. code:: pycon

      >>> from Bio import Align
      >>> aligner = Align.PairwiseAligner()
      >>> aligner.match_score
      1.0
      >>> aligner.mismatch_score
      0.0
      >>> score = aligner.score("ACGT", "ACAT")
      >>> print(score)
      3.0
      >>> aligner.match_score = 1.0
      >>> aligner.mismatch_score = -2.0
      >>> aligner.gap_score = -2.5
      >>> score = aligner.score("ACGT", "ACAT")
      >>> print(score)
      1.0

   When using match and mismatch scores, you can specify a wildcard
   character (``None`` by default) for unknown letters. These will get a
   zero score in alignments, irrespective of the value of the match or
   mismatch score:

   .. cont-doctest

   .. code:: pycon

      >>> aligner.wildcard = "?"
      >>> score = aligner.score("ACGT", "AC?T")
      >>> print(score)
      3.0

-  Alternatively, you can use the ``substitution_matrix`` attribute of
   the ``PairwiseAligner`` object to specify a substitution matrix. This
   allows you to apply different scores for different pairs of matched
   and mismatched letters. This is typically used for amino acid
   sequence alignments. For example, by default BLAST
   [Altschul1990]_ uses the BLOSUM62 substitution
   matrix for protein alignments by ``blastp``. This substitution matrix
   is available from Biopython:

   .. cont-doctest

   .. code:: pycon

      >>> from Bio.Align import substitution_matrices
      >>> substitution_matrices.load()  # doctest: +ELLIPSIS
      ['BENNER22', 'BENNER6', 'BENNER74', 'BLASTN', 'BLASTP', 'BLOSUM45', 'BLOSUM50', 'BLOSUM62', ..., 'TRANS']
      >>> matrix = substitution_matrices.load("BLOSUM62")
      >>> print(matrix)  # doctest: +ELLIPSIS
      #  Matrix made by matblas from blosum62.iij
      ...
           A    R    N    D    C    Q ...
      A  4.0 -1.0 -2.0 -2.0  0.0 -1.0 ...
      R -1.0  5.0  0.0 -2.0 -3.0  1.0 ...
      N -2.0  0.0  6.0  1.0 -3.0  0.0 ...
      D -2.0 -2.0  1.0  6.0 -3.0  0.0 ...
      C  0.0 -3.0 -3.0 -3.0  9.0 -3.0 ...
      Q -1.0  1.0  0.0  0.0 -3.0  5.0 ...
      ...
      >>> aligner.substitution_matrix = matrix
      >>> score = aligner.score("ACDQ", "ACDQ")
      >>> score
      24.0
      >>> score = aligner.score("ACDQ", "ACNQ")
      >>> score
      19.0

   When using a substitution matrix, ``X`` is *not* interpreted as an
   unknown character. Instead, the score provided by the substitution
   matrix will be used:

   .. cont-doctest

   .. code:: pycon

      >>> matrix["D", "X"]
      -1.0
      >>> score = aligner.score("ACDQ", "ACXQ")
      >>> score
      17.0

By default, ``aligner.substitution_matrix`` is ``None``. The attributes
``aligner.match_score`` and ``aligner.mismatch_score`` are ignored if
``aligner.substitution_matrix`` is not ``None``. Setting
``aligner.match_score`` or ``aligner.mismatch_score`` to valid values
will reset ``aligner.substitution_matrix`` to ``None``.

.. _`sec:pairwise-affine-gapscores`:

Affine gap scores
-----------------

Affine gap scores are defined by a score to open a gap, and a score to
extend an existing gap:

:math:`\textrm{gap score} = \textrm{open gap score} + (n-1) \times \textrm{extend gap score}`,

where :math:`n` is the length of the gap. Biopython’s pairwise sequence
aligner allows fine-grained control over the gap scoring scheme by
specifying the following twelve attributes of a ``PairwiseAligner``
object:

================================== ====================================
**Opening scores**                 **Extending scores**
================================== ====================================
``open_left_deletion_score``       ``extend_left_deletion_score``
``open_internal_deletion_score``   ``extend_internal_deletion_score``
``open_right_deletion_score``      ``extend_right_deletion_score``
``open_left_insertion_score``      ``extend_left_insertion_score``
``open_internal_insertion_score``  ``extend_internal_insertion_score``
``open_right_insertion_score``     ``extend_right_insertion_score``
================================== ====================================

These attributes allow for different gap scores for internal gaps and on
either end of the sequence, as shown in this example:

========== ========= ===============================
**target** **query** **score**
========== ========= ===============================
A          -         open left deletion score
C          -         extend left deletion score
C          -         extend left deletion score
G          G         match score
G          T         mismatch score
G          -         open internal deletion score
A          -         extend internal deletion score
A          -         extend internal deletion score
T          T         match score
A          A         match score
G          -         open internal deletion score
C          C         match score
-          C         open internal insertion score
-          C         extend internal insertion score
C          C         match score
T          G         mismatch score
C          C         match score
-          C         open internal insertion score
A          A         match score
-          T         open right insertion score
-          A         extend right insertion score
-          A         extend right insertion score
========== ========= ===============================

For convenience, ``PairwiseAligner`` objects have additional attributes
that refer to a number of these values collectively, as shown
(hierarchically) in Table :ref:`table:align-meta-attributes`.

.. table:: Meta-attributes of the pairwise aligner objects.
   :name: table:align-meta-attributes

   +--------------------------------+--------------------------------------+
   | Meta-attribute                 | Attributes it maps to                |
   +================================+======================================+
   | ``gap_score``                  | ``insertion_score``,                 |
   |                                | ``deletion_score``                   |
   +--------------------------------+--------------------------------------+
   | ``open_gap_score``             | ``open_insertion_score``,            |
   |                                | ``open_deletion_score``              |
   +--------------------------------+--------------------------------------+
   | ``extend_gap_score``           | ``extend_insertion_score``,          |
   |                                | ``extend_deletion_score``            |
   +--------------------------------+--------------------------------------+
   | ``internal_gap_score``         | ``internal_insertion_score``,        |
   |                                | ``internal_deletion_score``          |
   +--------------------------------+--------------------------------------+
   | ``open_internal_gap_score``    | ``open_internal_insertion_score``,   |
   |                                | ``open_internal_deletion_score``     |
   +--------------------------------+--------------------------------------+
   | ``extend_internal_gap_score``  | ``extend_internal_insertion_score``, |
   |                                | ``extend_internal_deletion_score``   |
   +--------------------------------+--------------------------------------+
   | ``end_gap_score``              | ``end_insertion_score``,             |
   |                                | ``end_deletion_score``               |
   +--------------------------------+--------------------------------------+
   | ``open_end_gap_score``         | ``open_end_insertion_score``,        |
   |                                | ``open_end_deletion_score``          |
   +--------------------------------+--------------------------------------+
   | ``extend_end_gap_score``       | ``extend_end_insertion_score``,      |
   |                                | ``extend_end_deletion_score``        |
   +--------------------------------+--------------------------------------+
   | ``left_gap_score``             | ``left_insertion_score``,            |
   |                                | ``left_deletion_score``              |
   +--------------------------------+--------------------------------------+
   | ``right_gap_score``            | ``right_insertion_score``,           |
   |                                | ``right_deletion_score``             |
   +--------------------------------+--------------------------------------+
   | ``open_left_gap_score``        | ``open_left_insertion_score``,       |
   |                                | ``open_left_deletion_score``         |
   +--------------------------------+--------------------------------------+
   | ``extend_left_gap_score``      | ``extend_left_insertion_score``,     |
   |                                | ``extend_left_deletion_score``       |
   +--------------------------------+--------------------------------------+
   | ``open_right_gap_score``       | ``open_right_insertion_score``,      |
   |                                | ``open_right_deletion_score``        |
   +--------------------------------+--------------------------------------+
   | ``extend_right_gap_score``     | ``extend_right_insertion_score``,    |
   |                                | ``extend_right_deletion_score``      |
   +--------------------------------+--------------------------------------+
   | ``open_insertion_score``       | ``open_internal_insertion_score``,   |
   |                                | ``open_left_insertion_score``,       |
   |                                | ``open_right_insertion_score``       |
   +--------------------------------+--------------------------------------+
   | ``extend_insertion_score``     | ``extend_internal_insertion_score``, |
   |                                | ``extend_left_insertion_score``,     |
   |                                | ``extend_right_insertion_score``     |
   +--------------------------------+--------------------------------------+
   | ``insertion_score``            | ``open_insertion_score``,            |
   |                                | ``extend_insertion_score``           |
   +--------------------------------+--------------------------------------+
   | ``open_deletion_score``        | ``open_internal_deletion_score``,    |
   |                                | ``open_left_deletion_score``,        |
   |                                | ``open_right_deletion_score``        |
   +--------------------------------+--------------------------------------+
   | ``extend_deletion_score``      | ``extend_internal_deletion_score``,  |
   |                                | ``extend_left_deletion_score``,      |
   |                                | ``extend_right_deletion_score``      |
   +--------------------------------+--------------------------------------+
   | ``deletion_score``             | ``open_deletion_score``,             |
   |                                | ``extend_deletion_score``            |
   +--------------------------------+--------------------------------------+
   | ``internal_insertion_score``   | ``open_internal_insertion_score``,   |
   |                                | ``extend_internal_insertion_score``  |
   +--------------------------------+--------------------------------------+
   | ``end_insertion_score``        | ``open_end_insertion_score``,        |
   |                                | ``extend_end_insertion_score``       |
   +--------------------------------+--------------------------------------+
   | ``open_end_insertion_score``   | ``open_left_insertion_score``,       |
   |                                | ``open_right_insertion_score``       |
   +--------------------------------+--------------------------------------+
   | ``extend_end_insertion_score`` | ``extend_left_insertion_score``,     |
   |                                | ``extend_right_insertion_score``     |
   +--------------------------------+--------------------------------------+
   | ``left_insertion_score``       | ``open_left_insertion_score``,       |
   |                                | ``extend_left_insertion_score``      |
   +--------------------------------+--------------------------------------+
   | ``right_insertion_score``      | ``open_right_insertion_score``,      |
   |                                | ``extend_right_insertion_score``     |
   +--------------------------------+--------------------------------------+
   | ``end_deletion_score``         | ``open_end_deletion_score``,         |
   |                                | ``extend_end_deletion_score``        |
   +--------------------------------+--------------------------------------+
   | ``open_end_deletion_score``    | ``open_left_deletion_score``,        |
   |                                | ``open_right_deletionp_score``       |
   +--------------------------------+--------------------------------------+
   | ``extend_end_deletion_score``  | ``extend_left_deletion_score``,      |
   |                                | ``extend_right_deletion_score``      |
   +--------------------------------+--------------------------------------+
   | ``internal_deletion_score``    | ``open_internal_deletion_score``,    |
   |                                | ``extend_internal_deletion_score``   |
   +--------------------------------+--------------------------------------+
   | ``left_deletion_score``        | ``open_left_deletion_score``,        |
   |                                | ``extend_left_deletion_score``       |
   +--------------------------------+--------------------------------------+
   | ``right_deletion_score``       | ``open_right_deletion_score``,       |
   |                                | ``extend_right_deletion_score``      |
   +--------------------------------+--------------------------------------+

.. _`sec:pairwise-general-gapscores`:

General gap scores
------------------

For even more fine-grained control over the gap scores, you can specify
a gap scoring function. For example, the gap scoring function below
disallows a deletion after two nucleotides in the query sequence:

.. doctest

.. code:: pycon

   >>> from Bio import Align
   >>> aligner = Align.PairwiseAligner()
   >>> def my_gap_score_function(start, length):
   ...     if start == 2:
   ...         return -1000
   ...     else:
   ...         return -1 * length
   ...
   >>> aligner.deletion_score = my_gap_score_function
   >>> alignments = aligner.align("AACTT", "AATT")
   >>> for alignment in alignments:
   ...     print(alignment)
   ...
   target            0 AACTT 5
                     0 -|.|| 5
   query             0 -AATT 4
   <BLANKLINE>
   target            0 AACTT 5
                     0 |-.|| 5
   query             0 A-ATT 4
   <BLANKLINE>
   target            0 AACTT 5
                     0 ||.-| 5
   query             0 AAT-T 4
   <BLANKLINE>
   target            0 AACTT 5
                     0 ||.|- 5
   query             0 AATT- 4
   <BLANKLINE>

.. _`sec:pairwise-predefined-scoring`:

Using a pre-defined substitution matrix and gap scores
------------------------------------------------------

By default, a ``PairwiseAligner`` object is initialized with a match
score of +1.0, a mismatch score of 0.0, and all gap scores equal to 0.0,
While this has the benefit of being a simple scoring scheme, in general
it does not give the best performance. Instead, you can use the argument
``scoring`` to select a predefined scoring scheme when initializing a
``PairwiseAligner`` object. Currently, the provided scoring schemes are
``blastn`` and ``megablast``, which are suitable for nucleotide
alignments, and ``blastp``, which is suitable for protein alignments.
Selecting these scoring schemes will initialize the ``PairwiseAligner``
object to the default scoring parameters used by BLASTN, MegaBLAST, and
BLASTP, respectively.

.. doctest

.. code:: pycon

   >>> from Bio import Align
   >>> aligner = Align.PairwiseAligner(scoring="blastn")
   >>> print(aligner)
   Pairwise sequence aligner with parameters
     substitution_matrix: <Array object at ...>
     open_internal_insertion_score: -7.000000
     extend_internal_insertion_score: -2.000000
     open_left_insertion_score: -7.000000
     extend_left_insertion_score: -2.000000
     open_right_insertion_score: -7.000000
     extend_right_insertion_score: -2.000000
     open_internal_deletion_score: -7.000000
     extend_internal_deletion_score: -2.000000
     open_left_deletion_score: -7.000000
     extend_left_deletion_score: -2.000000
     open_right_deletion_score: -7.000000
     extend_right_deletion_score: -2.000000
     mode: global
   <BLANKLINE>
   >>> print(aligner.substitution_matrix[:, :])
        A    T    G    C    S    W    R    Y    K    M    B    V    H    D    N
   A  2.0 -3.0 -3.0 -3.0 -3.0 -1.0 -1.0 -3.0 -3.0 -1.0 -3.0 -1.0 -1.0 -1.0 -2.0
   T -3.0  2.0 -3.0 -3.0 -3.0 -1.0 -3.0 -1.0 -1.0 -3.0 -1.0 -3.0 -1.0 -1.0 -2.0
   G -3.0 -3.0  2.0 -3.0 -1.0 -3.0 -1.0 -3.0 -1.0 -3.0 -1.0 -1.0 -3.0 -1.0 -2.0
   C -3.0 -3.0 -3.0  2.0 -1.0 -3.0 -3.0 -1.0 -3.0 -1.0 -1.0 -1.0 -1.0 -3.0 -2.0
   S -3.0 -3.0 -1.0 -1.0 -1.0 -3.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -2.0
   W -1.0 -1.0 -3.0 -3.0 -3.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -2.0
   R -1.0 -3.0 -1.0 -3.0 -1.0 -1.0 -1.0 -3.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -2.0
   Y -3.0 -1.0 -3.0 -1.0 -1.0 -1.0 -3.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -2.0
   K -3.0 -1.0 -1.0 -3.0 -1.0 -1.0 -1.0 -1.0 -1.0 -3.0 -1.0 -1.0 -1.0 -1.0 -2.0
   M -1.0 -3.0 -3.0 -1.0 -1.0 -1.0 -1.0 -1.0 -3.0 -1.0 -1.0 -1.0 -1.0 -1.0 -2.0
   B -3.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -2.0
   V -1.0 -3.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -2.0
   H -1.0 -1.0 -3.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -2.0
   D -1.0 -1.0 -1.0 -3.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -2.0
   N -2.0 -2.0 -2.0 -2.0 -2.0 -2.0 -2.0 -2.0 -2.0 -2.0 -2.0 -2.0 -2.0 -2.0 -2.0
   <BLANKLINE>

Iterating over alignments
-------------------------

The ``alignments`` returned by ``aligner.align`` are a kind of immutable
iterable objects (similar to ``range``). While they appear similar to a
``tuple`` or ``list`` of ``Alignment`` objects, they are different in
the sense that each ``Alignment`` object is created dynamically when it
is needed. This approach was chosen because the number of alignments can
be extremely large, in particular for poor alignments (see
Section :ref:`sec:pairwise-examples` for an example).

You can perform the following operations on ``alignments``:

-  ``len(alignments)`` returns the number of alignments stored. This
   function returns quickly, even if the number of alignments is huge.
   If the number of alignments is extremely large (typically, larger
   than 9,223,372,036,854,775,807, which is the largest integer that can
   be stored as a ``long int`` on 64 bit machines), ``len(alignments)``
   will raise an ``OverflowError``. A large number of alignments
   suggests that the alignment quality is low.

   .. doctest

   .. code:: pycon

      >>> from Bio import Align
      >>> aligner = Align.PairwiseAligner()
      >>> alignments = aligner.align("AAA", "AA")
      >>> len(alignments)
      3

-  You can extract a specific alignment by index:

   .. doctest

   .. code:: pycon

      >>> from Bio import Align
      >>> aligner = Align.PairwiseAligner()
      >>> alignments = aligner.align("AAA", "AA")
      >>> print(alignments[2])
      target            0 AAA 3
                        0 -|| 3
      query             0 -AA 2
      <BLANKLINE>
      >>> print(alignments[0])
      target            0 AAA 3
                        0 ||- 3
      query             0 AA- 2
      <BLANKLINE>

-  You can iterate over alignments, for example as in

   .. code:: pycon

      >>> for alignment in alignments:
      ...     print(alignment)
      ...

   The ``alignments`` iterator can be converted into a ``list`` or ``tuple``:

   .. code:: pycon

      >>> alignments = list(alignments)

   It is wise to check the number of alignments by calling
   ``len(alignments)`` before attempting to call ``list(alignments)`` to
   save all alignments as a list.

-  The alignment score (which has the same value for each alignment in
   ``alignments``) is stored as an attribute. This allows you to check
   the alignment score before proceeding to extract individual
   alignments:

   .. cont-doctest

   .. code:: pycon

      >>> print(alignments.score)
      2.0

Aligning to the reverse strand
------------------------------

By default, the pairwise aligner aligns the forward strand of the query
to the forward strand of the target. To calculate the alignment score
for ``query`` to the reverse strand of ``target``, use ``strand="-"``:

.. doctest

.. code:: pycon

   >>> from Bio import Align
   >>> from Bio.Seq import reverse_complement
   >>> target = "AAAACCC"
   >>> query = "AACC"
   >>> aligner = Align.PairwiseAligner()
   >>> aligner.mismatch_score = -1
   >>> aligner.internal_gap_score = -1
   >>> aligner.score(target, query)  # strand is "+" by default
   4.0
   >>> aligner.score(target, reverse_complement(query), strand="-")
   4.0
   >>> aligner.score(target, query, strand="-")
   0.0
   >>> aligner.score(target, reverse_complement(query))
   0.0

The alignments against the reverse strand can be obtained by specifying
``strand="-"`` when calling ``aligner.align``:

.. cont-doctest

.. code:: pycon

   >>> alignments = aligner.align(target, query)
   >>> len(alignments)
   1
   >>> print(alignments[0])
   target            0 AAAACCC 7
                     0 --||||- 7
   query             0 --AACC- 4
   <BLANKLINE>
   >>> print(alignments[0].format("bed"))  # doctest: +NORMALIZE_WHITESPACE
   target   2   6   query   4   +   2   6   0   1   4,   0,
   <BLANKLINE>
   >>> alignments = aligner.align(target, reverse_complement(query), strand="-")
   >>> len(alignments)
   1
   >>> print(alignments[0])
   target            0 AAAACCC 7
                     0 --||||- 7
   query             4 --AACC- 0
   <BLANKLINE>
   >>> print(alignments[0].format("bed"))  # doctest: +NORMALIZE_WHITESPACE
   target   2   6   query   4   -   2   6   0   1   4,   0,
   <BLANKLINE>
   >>> alignments = aligner.align(target, query, strand="-")
   >>> len(alignments)
   2
   >>> print(alignments[0])
   target            0 AAAACCC----  7
                     0 ----------- 11
   query             4 -------GGTT  0
   <BLANKLINE>
   >>> print(alignments[1])
   target            0 ----AAAACCC  7
                     0 ----------- 11
   query             4 GGTT-------  0
   <BLANKLINE>

Note that the score for aligning ``query`` to the reverse strand of
``target`` may be different from the score for aligning the reverse
complement of ``query`` to the forward strand of ``target`` if the left
and right gap scores are different:

.. cont-doctest

.. code:: pycon

   >>> aligner.left_gap_score = -0.5
   >>> aligner.right_gap_score = -0.2
   >>> aligner.score(target, query)
   2.8
   >>> alignments = aligner.align(target, query)
   >>> len(alignments)
   1
   >>> print(alignments[0])
   target            0 AAAACCC 7
                     0 --||||- 7
   query             0 --AACC- 4
   <BLANKLINE>
   >>> aligner.score(target, reverse_complement(query), strand="-")
   3.1
   >>> alignments = aligner.align(target, reverse_complement(query), strand="-")
   >>> len(alignments)
   1
   >>> print(alignments[0])
   target            0 AAAACCC 7
                     0 --||||- 7
   query             4 --AACC- 0
   <BLANKLINE>

.. _`sec:substitution_matrices`:

Substitution matrices
---------------------

Substitution matrices [Durbin1998]_ provide the scoring
terms for classifying how likely two different residues are to
substitute for each other. This is essential in doing sequence
comparisons. Biopython provides a ton of common substitution matrices,
including the famous PAM and BLOSUM series of matrices, and also
provides functionality for creating your own substitution matrices.

Array objects
~~~~~~~~~~~~~

You can think of substitutions matrices as 2D arrays in which the
indices are letters (nucleotides or amino acids) rather than integers.
The ``Array`` class in ``Bio.Align.substitution_matrices`` is a subclass
of numpy arrays that supports indexing both by integers and by specific
strings. An ``Array`` instance can either be a one-dimensional array or
a square two-dimensional arrays. A one-dimensional ``Array`` object can
for example be used to store the nucleotide frequency of a DNA sequence,
while a two-dimensional ``Array`` object can be used to represent a
scoring matrix for sequence alignments.

To create a one-dimensional ``Array``, only the alphabet of allowed
letters needs to be specified:

.. doctest . lib:numpy

.. code:: pycon

   >>> from Bio.Align.substitution_matrices import Array
   >>> counts = Array("ACGT")
   >>> print(counts)
   A 0.0
   C 0.0
   G 0.0
   T 0.0
   <BLANKLINE>

The allowed letters are stored in the ``alphabet`` property:

.. cont-doctest

.. code:: pycon

   >>> counts.alphabet
   'ACGT'

This property is read-only; modifying the underlying ``_alphabet``
attribute may lead to unexpected results. Elements can be accessed both
by letter and by integer index:

.. cont-doctest

.. code:: pycon

   >>> counts["C"] = -3
   >>> counts[2] = 7
   >>> print(counts)
   A  0.0
   C -3.0
   G  7.0
   T  0.0
   <BLANKLINE>
   >>> counts[1]
   -3.0

Using a letter that is not in the alphabet, or an index that is out of
bounds, will cause a ``IndexError``:

.. cont-doctest

.. code:: pycon

   >>> counts["U"]
   Traceback (most recent call last):
       ...
   IndexError: 'U'
   >>> counts["X"] = 6
   Traceback (most recent call last):
       ...
   IndexError: 'X'
   >>> counts[7]
   Traceback (most recent call last):
       ...
   IndexError: index 7 is out of bounds for axis 0 with size 4

A two-dimensional ``Array`` can be created by specifying ``dims=2``:

.. doctest . lib:numpy

.. code:: pycon

   >>> from Bio.Align.substitution_matrices import Array
   >>> counts = Array("ACGT", dims=2)
   >>> print(counts)
       A   C   G   T
   A 0.0 0.0 0.0 0.0
   C 0.0 0.0 0.0 0.0
   G 0.0 0.0 0.0 0.0
   T 0.0 0.0 0.0 0.0
   <BLANKLINE>

Again, both letters and integers can be used for indexing, and
specifying a letter that is not in the alphabet will cause an
``IndexError``:

.. cont-doctest

.. code:: pycon

   >>> counts["A", "C"] = 12.0
   >>> counts[2, 1] = 5.0
   >>> counts[3, "T"] = -2
   >>> print(counts)
       A    C   G    T
   A 0.0 12.0 0.0  0.0
   C 0.0  0.0 0.0  0.0
   G 0.0  5.0 0.0  0.0
   T 0.0  0.0 0.0 -2.0
   <BLANKLINE>
   >>> counts["X", 1]
   Traceback (most recent call last):
       ...
   IndexError: 'X'
   >>> counts["A", 5]
   Traceback (most recent call last):
       ...
   IndexError: index 5 is out of bounds for axis 1 with size 4

Selecting a row or column from the two-dimensional array will return a
one-dimensional ``Array``:

.. cont-doctest

.. code:: pycon

   >>> counts = Array("ACGT", dims=2)
   >>> counts["A", "C"] = 12.0
   >>> counts[2, 1] = 5.0
   >>> counts[3, "T"] = -2

.. code:: pycon

   >>> counts["G"]
   Array([0., 5., 0., 0.],
         alphabet='ACGT')
   >>> counts[:, "C"]
   Array([12.,  0.,  5.,  0.],
         alphabet='ACGT')

``Array`` objects can thus be used as an array and as a dictionary. They
can be converted to plain numpy arrays or plain dictionary objects:

.. cont-doctest

.. code:: pycon

   >>> import numpy as np
   >>> x = Array("ACGT")
   >>> x["C"] = 5

.. code:: pycon

   >>> x
   Array([0., 5., 0., 0.],
         alphabet='ACGT')
   >>> a = np.array(x)  # create a plain numpy array
   >>> a
   array([0., 5., 0., 0.])
   >>> d = dict(x)  # create a plain dictionary
   >>> d
   {'A': 0.0, 'C': 5.0, 'G': 0.0, 'T': 0.0}

While the alphabet of an ``Array`` is usually a string, you may also use
a tuple of (immutable) objects. This is used for example for a codon
substitution matrix (as in the
``substitution_matrices.load("SCHNEIDER")`` example shown later), where
the keys are not individual nucleotides or amino acids but instead
three-nucleotide codons.

While the ``alphabet`` property of an ``Array`` is immutable, you can
create a new ``Array`` object by selecting the letters you are
interested in from the alphabet. For example,

.. cont-doctest

.. code:: pycon

   >>> a = Array("ABCD", dims=2, data=np.arange(16).reshape(4, 4))
   >>> print(a)
        A    B    C    D
   A  0.0  1.0  2.0  3.0
   B  4.0  5.0  6.0  7.0
   C  8.0  9.0 10.0 11.0
   D 12.0 13.0 14.0 15.0
   <BLANKLINE>
   >>> b = a.select("CAD")
   >>> print(b)
        C    A    D
   C 10.0  8.0 11.0
   A  2.0  0.0  3.0
   D 14.0 12.0 15.0
   <BLANKLINE>

Note that this also allows you to reorder the alphabet.

Data for letters that are not found in the alphabet are set to zero:

.. cont-doctest

.. code:: pycon

   >>> c = a.select("DEC")
   >>> print(c)
        D   E    C
   D 15.0 0.0 14.0
   E  0.0 0.0  0.0
   C 11.0 0.0 10.0
   <BLANKLINE>

As the ``Array`` class is a subclass of numpy array, it can be used as
such. A ``ValueError`` is triggered if the ``Array`` objects appearing
in a mathematical operation have different alphabets, for example

.. doctest . lib:numpy

.. code:: pycon

   >>> from Bio.Align.substitution_matrices import Array
   >>> d = Array("ACGT")
   >>> r = Array("ACGU")
   >>> d + r
   Traceback (most recent call last):
       ...
   ValueError: alphabets are inconsistent

Calculating a substitution matrix from a pairwise sequence alignment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As ``Array`` is a subclass of a numpy array, you can apply mathematical
operations on an ``Array`` object in much the same way. Here, we
illustrate this by calculating a scoring matrix from the alignment of
the 16S ribosomal RNA gene sequences of *Escherichia coli* and *Bacillus
subtilis*. First, we create a ``PairwiseAligner`` object (see
Chapter :ref:`chapter:pairwise`) and initialize it with the default
scores used by ``blastn``:

.. doctest ../Tests/Align lib:numpy

.. code:: pycon

   >>> from Bio.Align import PairwiseAligner
   >>> aligner = PairwiseAligner(scoring="blastn")
   >>> aligner.mode = "local"

Next, we read in the 16S ribosomal RNA gene sequence of *Escherichia
coli* and *Bacillus subtilis* (provided in ``Tests/Align/ecoli.fa`` and
``Tests/Align/bsubtilis.fa``), and align them to each other:

.. cont-doctest

.. code:: pycon

   >>> from Bio import SeqIO
   >>> sequence1 = SeqIO.read("ecoli.fa", "fasta")
   >>> sequence2 = SeqIO.read("bsubtilis.fa", "fasta")
   >>> alignments = aligner.align(sequence1, sequence2)

The number of alignments generated is very large:

.. cont-doctest

.. code:: pycon

   >>> len(alignments)
   1990656

However, as they only differ trivially from each other, we arbitrarily
choose the first alignment, and count the number of each substitution:

.. cont-doctest

.. code:: pycon

   >>> alignment = alignments[0]
   >>> substitutions = alignment.substitutions
   >>> print(substitutions)
         A     C     G     T
   A 307.0  19.0  34.0  19.0
   C  15.0 280.0  25.0  29.0
   G  34.0  24.0 401.0  20.0
   T  24.0  36.0  20.0 228.0
   <BLANKLINE>

We normalize against the total number to find the probability of each
substitution, and create a symmetric matrix of observed frequencies:

.. cont-doctest

.. code:: pycon

   >>> observed_frequencies = substitutions / substitutions.sum()
   >>> observed_frequencies = (observed_frequencies + observed_frequencies.transpose()) / 2.0
   >>> print(format(observed_frequencies, "%.4f"))
          A      C      G      T
   A 0.2026 0.0112 0.0224 0.0142
   C 0.0112 0.1848 0.0162 0.0215
   G 0.0224 0.0162 0.2647 0.0132
   T 0.0142 0.0215 0.0132 0.1505
   <BLANKLINE>

The background probability is the probability of finding an A, C, G, or
T nucleotide in each sequence separately. This can be calculated as the
sum of each row or column:

.. cont-doctest

.. code:: pycon

   >>> background = observed_frequencies.sum(0)
   >>> print(format(background, "%.4f"))
   A 0.2505
   C 0.2337
   G 0.3165
   T 0.1993
   <BLANKLINE>

The number of substitutions expected at random is simply the product of
the background distribution with itself:

.. cont-doctest

.. code:: pycon

   >>> expected_frequencies = background[:, None].dot(background[None, :])
   >>> print(format(expected_frequencies, "%.4f"))
          A      C      G      T
   A 0.0627 0.0585 0.0793 0.0499
   C 0.0585 0.0546 0.0740 0.0466
   G 0.0793 0.0740 0.1002 0.0631
   T 0.0499 0.0466 0.0631 0.0397
   <BLANKLINE>

The scoring matrix can then be calculated as the logarithm of the
odds-ratio of the observed and the expected probabilities:

.. cont-doctest

.. code:: pycon

   >>> oddsratios = observed_frequencies / expected_frequencies
   >>> import numpy as np
   >>> scoring_matrix = np.log2(oddsratios)
   >>> print(scoring_matrix)
        A    C    G    T
   A  1.7 -2.4 -1.8 -1.8
   C -2.4  1.8 -2.2 -1.1
   G -1.8 -2.2  1.4 -2.3
   T -1.8 -1.1 -2.3  1.9
   <BLANKLINE>

The matrix can be used to set the substitution matrix for the pairwise
aligner (see Chapter :ref:`chapter:pairwise`):

.. cont-doctest

.. code:: pycon

   >>> aligner.substitution_matrix = scoring_matrix

.. _`subsec:subs_mat_ex`:

Calculating a substitution matrix from a multiple sequence alignment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this example, we’ll first read a protein sequence alignment from the
Clustalw file `protein.aln <examples/protein.aln>`__ (also available
online
`here <https://raw.githubusercontent.com/biopython/biopython/master/Tests/Clustalw/protein.aln>`__)

.. doctest ../Tests/Clustalw lib:numpy

.. code:: pycon

   >>> from Bio import Align
   >>> filename = "protein.aln"
   >>> alignment = Align.read(filename, "clustal")

Section :ref:`subsec:align_clustal` contains more
information on doing this.

The ``substitutions`` property of the alignment stores the number of
times different residues substitute for each other:

.. cont-doctest

.. code:: pycon

   >>> substitutions = alignment.substitutions

To make the example more readable, we’ll select only amino acids with
polar charged side chains:

.. cont-doctest

.. code:: pycon

   >>> substitutions = substitutions.select("DEHKR")
   >>> print(substitutions)
          D      E      H      K      R
   D 2360.0  270.0   15.0    1.0   48.0
   E  241.0 3305.0   15.0   45.0    2.0
   H    0.0   18.0 1235.0    8.0    0.0
   K    0.0    9.0   24.0 3218.0  130.0
   R    2.0    2.0   17.0  103.0 2079.0
   <BLANKLINE>

Rows and columns for other amino acids were removed from the matrix.

Next, we normalize the matrix and make it symmetric.

.. cont-doctest

.. code:: pycon

   >>> observed_frequencies = substitutions / substitutions.sum()
   >>> observed_frequencies = (observed_frequencies + observed_frequencies.transpose()) / 2.0
   >>> print(format(observed_frequencies, "%.4f"))
          D      E      H      K      R
   D 0.1795 0.0194 0.0006 0.0000 0.0019
   E 0.0194 0.2514 0.0013 0.0021 0.0002
   H 0.0006 0.0013 0.0939 0.0012 0.0006
   K 0.0000 0.0021 0.0012 0.2448 0.0089
   R 0.0019 0.0002 0.0006 0.0089 0.1581
   <BLANKLINE>

Summing over rows or columns gives the relative frequency of occurrence
of each residue:

.. cont-doctest

.. code:: pycon

   >>> background = observed_frequencies.sum(0)
   >>> print(format(background, "%.4f"))
   D 0.2015
   E 0.2743
   H 0.0976
   K 0.2569
   R 0.1697
   <BLANKLINE>
   >>> sum(background) == 1.0
   True

The expected frequency of residue pairs is then

.. cont-doctest

.. code:: pycon

   >>> expected_frequencies = background[:, None].dot(background[None, :])
   >>> print(format(expected_frequencies, "%.4f"))
          D      E      H      K      R
   D 0.0406 0.0553 0.0197 0.0518 0.0342
   E 0.0553 0.0752 0.0268 0.0705 0.0465
   H 0.0197 0.0268 0.0095 0.0251 0.0166
   K 0.0518 0.0705 0.0251 0.0660 0.0436
   R 0.0342 0.0465 0.0166 0.0436 0.0288
   <BLANKLINE>

Here, ``background[:, None]`` creates a 2D array consisting of a single
column with the values of ``expected_frequencies``, and
``expected_frequencies[None, :]`` a 2D array with these values as a
single row. Taking their dot product (inner product) creates a matrix of
expected frequencies where each entry consists of two
``expected_frequencies`` values multiplied with each other. For example,
``expected_frequencies['D', 'E']`` is equal to
``residue_frequencies['D'] * residue_frequencies['E']``.

We can now calculate the log-odds matrix by dividing the observed
frequencies by the expected frequencies and taking the logarithm:

.. cont-doctest

.. code:: pycon

   >>> import numpy as np
   >>> scoring_matrix = np.log2(observed_frequencies / expected_frequencies)
   >>> print(scoring_matrix)
         D    E    H     K    R
   D   2.1 -1.5 -5.1 -10.4 -4.2
   E  -1.5  1.7 -4.4  -5.1 -8.3
   H  -5.1 -4.4  3.3  -4.4 -4.7
   K -10.4 -5.1 -4.4   1.9 -2.3
   R  -4.2 -8.3 -4.7  -2.3  2.5
   <BLANKLINE>

This matrix can be used as the substitution matrix when performing
alignments. For example,

.. cont-doctest

.. code:: pycon

   >>> from Bio.Align import PairwiseAligner
   >>> aligner = PairwiseAligner()
   >>> aligner.substitution_matrix = scoring_matrix
   >>> aligner.gap_score = -3.0
   >>> alignments = aligner.align("DEHEK", "DHHKK")
   >>> print(alignments[0])
   target            0 DEHEK 5
                     0 |.|.| 5
   query             0 DHHKK 5
   <BLANKLINE>
   >>> print("%.2f" % alignments.score)
   -2.18
   >>> score = (
   ...     scoring_matrix["D", "D"]
   ...     + scoring_matrix["E", "H"]
   ...     + scoring_matrix["H", "H"]
   ...     + scoring_matrix["E", "K"]
   ...     + scoring_matrix["K", "K"]
   ... )
   >>> print("%.2f" % score)
   -2.18

(see Chapter :ref:`chapter:pairwise` for details).

Reading ``Array`` objects from file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``Bio.Align.substitution_matrices`` includes a parser to read one- and
two-dimensional ``Array`` objects from file. One-dimensional arrays are
represented by a simple two-column format, with the first column
containing the key and the second column the corresponding value. For
example, the file ``hg38.chrom.sizes`` (obtained from UCSC), available
in the ``Tests/Align`` subdirectory of the Biopython distribution,
contains the size in nucleotides of each chromosome in human genome
assembly hg38:

.. code:: text

   chr1    248956422
   chr2    242193529
   chr3    198295559
   chr4    190214555
   ...
   chrUn_KI270385v1    990
   chrUn_KI270423v1    981
   chrUn_KI270392v1    971
   chrUn_KI270394v1    970

To parse this file, use

.. doctest ../Tests/Align lib:numpy

.. code:: pycon

   >>> from Bio.Align import substitution_matrices
   >>> with open("hg38.chrom.sizes") as handle:
   ...     table = substitution_matrices.read(handle)
   ...
   >>> print(table)  # doctest: +ELLIPSIS
   chr1 248956422.0
   chr2 242193529.0
   chr3 198295559.0
   chr4 190214555.0
   ...
   chrUn_KI270423v1       981.0
   chrUn_KI270392v1       971.0
   chrUn_KI270394v1       970.0
   <BLANKLINE>

Use ``dtype=int`` to read the values as integers:

.. cont-doctest

.. code:: pycon

   >>> with open("hg38.chrom.sizes") as handle:
   ...     table = substitution_matrices.read(handle, int)
   ...
   >>> print(table)  # doctest: +ELLIPSIS
   chr1 248956422
   chr2 242193529
   chr3 198295559
   chr4 190214555
   ...
   chrUn_KI270423v1       981
   chrUn_KI270392v1       971
   chrUn_KI270394v1       970
   <BLANKLINE>

For two-dimensional arrays, we follow the file format of substitution
matrices provided by NCBI. For example, the BLOSUM62 matrix, which is
the default substitution matrix for NCBI’s protein-protein BLAST
[Altschul1990]_ program ``blastp``, is stored as
follows:

.. code:: text

   #  Matrix made by matblas from blosum62.iij
   #  * column uses minimum score
   #  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
   #  Blocks Database = /data/blocks_5.0/blocks.dat
   #  Cluster Percentage: >= 62
   #  Entropy =   0.6979, Expected =  -0.5209
      A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
   A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
   R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
   N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
   D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
   C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
   Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
   E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
   G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
   H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
   ...

This file is included in the Biopython distribution under
``Bio/Align/substitution_matrices/data``. To parse this file, use

.. doctest ../Bio/Align/substitution_matrices/data lib:numpy

.. code:: pycon

   >>> from Bio.Align import substitution_matrices
   >>> with open("BLOSUM62") as handle:
   ...     matrix = substitution_matrices.read(handle)
   ...
   >>> print(matrix.alphabet)
   ARNDCQEGHILKMFPSTWYVBZX*
   >>> print(matrix["A", "D"])
   -2.0

The header lines starting with ``#`` are stored in the attribute
``header``:

.. cont-doctest

.. code:: pycon

   >>> matrix.header[0]
   'Matrix made by matblas from blosum62.iij'

We can now use this matrix as the substitution matrix on an aligner
object:

.. cont-doctest

.. code:: pycon

   >>> from Bio.Align import PairwiseAligner
   >>> aligner = PairwiseAligner()
   >>> aligner.substitution_matrix = matrix

To save an Array object, create a string first:

.. cont-doctest

.. code:: pycon

   >>> text = str(matrix)
   >>> print(text)  # doctest: +ELLIPSIS
   #  Matrix made by matblas from blosum62.iij
   #  * column uses minimum score
   #  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
   #  Blocks Database = /data/blocks_5.0/blocks.dat
   #  Cluster Percentage: >= 62
   #  Entropy =   0.6979, Expected =  -0.5209
        A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S ...
   A  4.0 -1.0 -2.0 -2.0  0.0 -1.0 -1.0  0.0 -2.0 -1.0 -1.0 -1.0 -1.0 -2.0 -1.0  1.0 ...
   R -1.0  5.0  0.0 -2.0 -3.0  1.0  0.0 -2.0  0.0 -3.0 -2.0  2.0 -1.0 -3.0 -2.0 -1.0 ...
   N -2.0  0.0  6.0  1.0 -3.0  0.0  0.0  0.0  1.0 -3.0 -3.0  0.0 -2.0 -3.0 -2.0  1.0 ...
   D -2.0 -2.0  1.0  6.0 -3.0  0.0  2.0 -1.0 -1.0 -3.0 -4.0 -1.0 -3.0 -3.0 -1.0  0.0 ...
   C  0.0 -3.0 -3.0 -3.0  9.0 -3.0 -4.0 -3.0 -3.0 -1.0 -1.0 -3.0 -1.0 -2.0 -3.0 -1.0 ...
   ...

and write the ``text`` to a file.

Loading predefined substitution matrices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Biopython contains a large set of substitution matrices defined in the
literature, including BLOSUM (Blocks Substitution Matrix)
[Henikoff1992]_ and PAM (Point Accepted Mutation)
matrices [Dayhoff1978]_. These matrices are available
as flat files in the ``Bio/Align/substitution_matrices/data`` directory,
and can be loaded into Python using the ``load`` function in the
``substitution_matrices`` submodule. For example, the BLOSUM62 matrix
can be loaded by running

.. doctest . lib:numpy

.. code:: pycon

   >>> from Bio.Align import substitution_matrices
   >>> m = substitution_matrices.load("BLOSUM62")

This substitution matrix has an alphabet consisting of the 20 amino
acids used in the genetic code, the three ambiguous amino acids B
(asparagine or aspartic acid), Z (glutamine or glutamic acid), and X
(representing any amino acid), and the stop codon represented by an
asterisk:

.. cont-doctest

.. code:: pycon

   >>> m.alphabet
   'ARNDCQEGHILKMFPSTWYVBZX*'

To get a full list of available substitution matrices, use ``load``
without an argument:

.. cont-doctest

.. code:: pycon

   >>> substitution_matrices.load()  # doctest: +ELLIPSIS
   ['BENNER22', 'BENNER6', 'BENNER74', 'BLASTN', 'BLASTP', 'BLOSUM45', 'BLOSUM50', ..., 'TRANS']

Note that the substitution matrix provided by Schneider *et al.*
[Schneider2005]_ uses an alphabet consisting of
three-nucleotide codons:

.. cont-doctest

.. code:: pycon

   >>> m = substitution_matrices.load("SCHNEIDER")
   >>> m.alphabet  # doctest: +ELLIPSIS
   ('AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', ..., 'TTG', 'TTT')

.. _`sec:pairwise-examples`:

Examples
--------

Suppose you want to do a global pairwise alignment between the same two
hemoglobin sequences from above (``HBA_HUMAN``, ``HBB_HUMAN``) stored in
``alpha.faa`` and ``beta.faa``:

.. doctest examples

.. code:: pycon

   >>> from Bio import Align
   >>> from Bio import SeqIO
   >>> seq1 = SeqIO.read("alpha.faa", "fasta")
   >>> seq2 = SeqIO.read("beta.faa", "fasta")
   >>> aligner = Align.PairwiseAligner()
   >>> score = aligner.score(seq1.seq, seq2.seq)
   >>> print(score)
   72.0

showing an alignment score of 72.0. To see the individual alignments, do

.. cont-doctest

.. code:: pycon

   >>> alignments = aligner.align(seq1.seq, seq2.seq)

In this example, the total number of optimal alignments is huge (more
than :math:`4 \times 10^{37}`), and calling ``len(alignments)`` will
raise an ``OverflowError``:

.. code:: pycon

   >>> len(alignments)
   Traceback (most recent call last):
   ...
   OverflowError: number of optimal alignments is larger than 9223372036854775807

Let’s have a look at the first alignment:

.. cont-doctest

.. code:: pycon

   >>> alignment = alignments[0]

The alignment object stores the alignment score, as well as the
alignment itself:

.. cont-doctest

.. code:: pycon

   >>> print(alignment.score)
   72.0
   >>> print(alignment)
   target            0 MV-LS-PAD--KTN--VK-AA-WGKV-----GAHAGEYGAEALE-RMFLSF----P-TTK
                     0 ||-|--|----|----|--|--||||-----|---||--|--|--|--|------|-|--
   query             0 MVHL-TP--EEK--SAV-TA-LWGKVNVDEVG---GE--A--L-GR--L--LVVYPWT--
   <BLANKLINE>
   target           41 TY--FPHF----DLSHGS---AQVK-G------HGKKV--A--DA-LTNAVAHV-DDMPN
                    60 ----|--|----|||------|-|--|------|||||--|--|--|--|--|--|---|
   query            39 --QRF--FESFGDLS---TPDA-V-MGNPKVKAHGKKVLGAFSD-GL--A--H-LD---N
   <BLANKLINE>
   target           79 ALS----A-LSD-LHAH--KLR-VDPV-NFK-LLSHC---LLVT--LAAHLPA----EFT
                   120 -|-----|-||--||----||--|||--||--||------|-|---||-|-------|||
   query            81 -L-KGTFATLS-ELH--CDKL-HVDP-ENF-RLL---GNVL-V-CVLA-H---HFGKEFT
   <BLANKLINE>
   target          119 PA-VH-ASLDKFLAS---VSTV------LTS--KYR- 142
                   180 |--|--|------|----|--|------|----||-- 217
   query           124 P-PV-QA------A-YQKV--VAGVANAL--AHKY-H 147
   <BLANKLINE>

Better alignments are usually obtained by penalizing gaps: higher costs
for opening a gap and lower costs for extending an existing gap. For
amino acid sequences match scores are usually encoded in matrices like
``PAM`` or ``BLOSUM``. Thus, a more meaningful alignment for our example
can be obtained by using the BLOSUM62 matrix, together with a gap open
penalty of 10 and a gap extension penalty of 0.5:

.. doctest examples lib:numpy

.. code:: pycon

   >>> from Bio import Align
   >>> from Bio import SeqIO
   >>> from Bio.Align import substitution_matrices
   >>> seq1 = SeqIO.read("alpha.faa", "fasta")
   >>> seq2 = SeqIO.read("beta.faa", "fasta")
   >>> aligner = Align.PairwiseAligner()
   >>> aligner.open_gap_score = -10
   >>> aligner.extend_gap_score = -0.5
   >>> aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
   >>> score = aligner.score(seq1.seq, seq2.seq)
   >>> print(score)
   292.5
   >>> alignments = aligner.align(seq1.seq, seq2.seq)
   >>> len(alignments)
   2
   >>> print(alignments[0].score)
   292.5
   >>> print(alignments[0])
   target            0 MV-LSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF-DLS-----HGS
                     0 ||-|.|..|..|.|.||||--...|.|.|||.|.....|.|...|..|-|||-----.|.
   query             0 MVHLTPEEKSAVTALWGKV--NVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGN
   <BLANKLINE>
   target           53 AQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAH
                    60 ..||.|||||..|.....||.|........||.||..||.|||.||.||...|...||.|
   query            58 PKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHH
   <BLANKLINE>
   target          113 LPAEFTPAVHASLDKFLASVSTVLTSKYR 142
                   120 ...||||.|.|...|..|.|...|..||. 149
   query           118 FGKEFTPPVQAAYQKVVAGVANALAHKYH 147
   <BLANKLINE>

This alignment has the same score that we obtained earlier with EMBOSS
needle using the same sequences and the same parameters.

To perform a local alignment, set ``aligner.mode`` to ``'local'``:

.. cont-doctest

.. code:: pycon

   >>> aligner.mode = "local"
   >>> aligner.open_gap_score = -10
   >>> aligner.extend_gap_score = -1
   >>> alignments = aligner.align("LSPADKTNVKAA", "PEEKSAV")
   >>> print(len(alignments))
   1
   >>> alignment = alignments[0]
   >>> print(alignment)
   target            2 PADKTNV 9
                     0 |..|..| 7
   query             0 PEEKSAV 7
   <BLANKLINE>
   >>> print(alignment.score)
   16.0

.. _`sec:generalized-pairwise`:

Generalized pairwise alignments
-------------------------------

In most cases, ``PairwiseAligner`` is used to perform alignments of
sequences (strings or ``Seq`` objects) consisting of single-letter
nucleotides or amino acids. More generally, ``PairwiseAligner`` can also
be applied to lists or tuples of arbitrary objects. This section will
describe some examples of such generalized pairwise alignments.

Generalized pairwise alignments using a substitution matrix and alphabet
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Schneider *et al.* [Schneider2005]_ created a
substitution matrix for aligning three-nucleotide codons (see
`below <#codonmatrix>`__ in section :ref:`sec:substitution_matrices`
for more information). This substitution matrix is associated with an
alphabet consisting of all three-letter codons:

.. doctest . lib:numpy

.. code:: pycon

   >>> from Bio.Align import substitution_matrices
   >>> m = substitution_matrices.load("SCHNEIDER")
   >>> m.alphabet  # doctest: +ELLIPSIS
   ('AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', ..., 'TTG', 'TTT')

We can use this matrix to align codon sequences to each other:

.. cont-doctest

.. code:: pycon

   >>> from Bio import Align
   >>> aligner = Align.PairwiseAligner()
   >>> aligner.substitution_matrix = m
   >>> aligner.gap_score = -1.0
   >>> s1 = ("AAT", "CTG", "TTT", "TTT")
   >>> s2 = ("AAT", "TTA", "TTT")
   >>> alignments = aligner.align(s1, s2)
   >>> len(alignments)
   2
   >>> print(alignments[0])
   AAT CTG TTT TTT
   ||| ... ||| ---
   AAT TTA TTT ---
   <BLANKLINE>
   >>> print(alignments[1])
   AAT CTG TTT TTT
   ||| ... --- |||
   AAT TTA --- TTT
   <BLANKLINE>

Note that aligning ``TTT`` to ``TTA``, as in this example:

.. code:: pycon

   AAT CTG TTT TTT
   ||| --- ... |||
   AAT --- TTA TTT

would get a much lower score:

.. cont-doctest

.. code:: pycon

   >>> print(m["CTG", "TTA"])
   7.6
   >>> print(m["TTT", "TTA"])
   -0.3

presumably because ``CTG`` and ``TTA`` both code for leucine, while
``TTT`` codes for phenylalanine. The three-letter codon substitution
matrix also reveals a preference among codons representing the same
amino acid. For example, ``TTA`` has a preference for ``CTG`` preferred
compared to ``CTC``, though all three code for leucine:

.. cont-doctest

.. code:: pycon

   >>> s1 = ("AAT", "CTG", "CTC", "TTT")
   >>> s2 = ("AAT", "TTA", "TTT")
   >>> alignments = aligner.align(s1, s2)
   >>> len(alignments)
   1
   >>> print(alignments[0])
   AAT CTG CTC TTT
   ||| ... --- |||
   AAT TTA --- TTT
   <BLANKLINE>
   >>> print(m["CTC", "TTA"])
   6.5

Generalized pairwise alignments using match/mismatch scores and an alphabet
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Using the three-letter amino acid symbols, the sequences above translate
to

.. doctest

.. code:: pycon

   >>> s1 = ("Asn", "Leu", "Leu", "Phe")
   >>> s2 = ("Asn", "Leu", "Phe")

We can align these sequences directly to each other by using a
three-letter amino acid alphabet:

.. cont-doctest

.. code:: pycon

   >>> from Bio import Align
   >>> aligner = Align.PairwiseAligner()
   >>> aligner.alphabet = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys',
   ...                     'Gln', 'Glu', 'Gly', 'His', 'Ile',
   ...                     'Leu', 'Lys', 'Met', 'Phe', 'Pro',
   ...                     'Ser', 'Thr', 'Trp', 'Tyr', 'Val']  # fmt: skip
   ...

We use +6/-1 match and mismatch scores as an approximation of the
BLOSUM62 matrix, and align these sequences to each other:

.. cont-doctest

.. code:: pycon

   >>> aligner.match = +6
   >>> aligner.mismatch = -1
   >>> alignments = aligner.align(s1, s2)
   >>> print(len(alignments))
   2
   >>> print(alignments[0])
   Asn Leu Leu Phe
   ||| ||| --- |||
   Asn Leu --- Phe
   <BLANKLINE>
   >>> print(alignments[1])
   Asn Leu Leu Phe
   ||| --- ||| |||
   Asn --- Leu Phe
   <BLANKLINE>
   >>> print(alignments.score)
   18.0

Generalized pairwise alignments using match/mismatch scores and integer sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Internally, the first step when performing an alignment is to replace
the two sequences by integer arrays consisting of the indices of each
letter in each sequence in the alphabet associated with the aligner.
This step can be bypassed by passing integer arrays directly:

.. doctest . lib:numpy

.. code:: pycon

   >>> import numpy as np
   >>> from Bio import Align
   >>> aligner = Align.PairwiseAligner()
   >>> s1 = np.array([2, 10, 10, 13], np.int32)
   >>> s2 = np.array([2, 10, 13], np.int32)
   >>> aligner.match = +6
   >>> aligner.mismatch = -1
   >>> alignments = aligner.align(s1, s2)
   >>> print(len(alignments))
   2
   >>> print(alignments[0])
   2 10 10 13
   | || -- ||
   2 10 -- 13
   <BLANKLINE>
   >>> print(alignments[1])
   2 10 10 13
   | -- || ||
   2 -- 10 13
   <BLANKLINE>
   >>> print(alignments.score)
   18.0

Note that the indices should consist of 32-bit integers, as specified in
this example by ``numpy.int32``.

Unknown letters can again be included by defining a wildcard character,
and using the corresponding Unicode code point number as the index:

.. cont-doctest

.. code:: pycon

   >>> aligner.wildcard = "?"
   >>> ord(aligner.wildcard)
   63
   >>> s2 = np.array([2, 63, 13], np.int32)
   >>> aligner.gap_score = -3
   >>> alignments = aligner.align(s1, s2)
   >>> print(len(alignments))
   2
   >>> print(alignments[0])
   2 10 10 13
   | .. -- ||
   2 63 -- 13
   <BLANKLINE>
   >>> print(alignments[1])
   2 10 10 13
   | -- .. ||
   2 -- 63 13
   <BLANKLINE>
   >>> print(alignments.score)
   9.0

Generalized pairwise alignments using a substitution matrix and integer sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Integer sequences can also be aligned using a substitution matrix, in
this case a numpy square array without an alphabet associated with it.
In this case, all index values must be non-negative, and smaller than
the size of the substitution matrix:

.. doctest . lib:numpy

.. code:: pycon

   >>> from Bio import Align
   >>> import numpy as np
   >>> aligner = Align.PairwiseAligner()
   >>> m = np.eye(5)
   >>> m[0, 1:] = m[1:, 0] = -2
   >>> m[2, 2] = 3
   >>> print(m)
   [[ 1. -2. -2. -2. -2.]
    [-2.  1.  0.  0.  0.]
    [-2.  0.  3.  0.  0.]
    [-2.  0.  0.  1.  0.]
    [-2.  0.  0.  0.  1.]]
   >>> aligner.substitution_matrix = m
   >>> aligner.gap_score = -1
   >>> s1 = np.array([0, 2, 3, 4], np.int32)
   >>> s2 = np.array([0, 3, 2, 1], np.int32)
   >>> alignments = aligner.align(s1, s2)
   >>> print(len(alignments))
   2
   >>> print(alignments[0])
   0 - 2 3 4
   | - | . -
   0 3 2 1 -
   <BLANKLINE>
   >>> print(alignments[1])
   0 - 2 3 4
   | - | - .
   0 3 2 - 1
   <BLANKLINE>
   >>> print(alignments.score)
   2.0

.. _`sec:codon_alignments`:

Codon alignments
----------------

The ``CodonAligner`` class in the ``Bio.Align`` module implements a
specialized aligner for aligning a nucleotide sequence to the amino acid
sequence it encodes. Such alignments are non-trivial if frameshifts
occur during translation.

Aligning a nucleotide sequence to an amino acid sequence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To align a nucleotide sequence to an amino acid sequence, first create a
``CodonAligner`` object:

.. doctest

.. code:: pycon

   >>> from Bio import Align
   >>> aligner = Align.CodonAligner()

The ``CodonAligner`` object ``aligner`` stores the alignment parameters
to be used for the alignments:

.. cont-doctest

.. code:: pycon

   >>> print(aligner)
   Codon aligner with parameters
     wildcard: 'X'
     match_score: 1.0
     mismatch_score: 0.0
     frameshift_minus_two_score: -3.0
     frameshift_minus_one_score: -3.0
     frameshift_plus_one_score: -3.0
     frameshift_plus_two_score: -3.0
   <BLANKLINE>

The ``wildcard``, ``match_score``, and ``mismatch_score`` parameters are
defined in the same was as for the ``PairwiseAligner`` class described
above (see Section :ref:`sec:pairwise-aligner`). The values
specified by the ``frameshift_minus_two_score``,
``frameshift_minus_one_score``, ``frameshift_plus_one_score``, and
``frameshift_plus_two_score`` parameters are added to the alignment
score whenever a -2, -1, +1, or +2 frame shift, respectively, occurs in
the alignment. By default, the frame shift scores are set to -3.0.
Similar to the ``PairwiseAligner`` class
(Table :ref:`table:align-meta-attributes`), the ``CodonAligner``
class defines additional attributes that refer to a number of these
values collectively, as shown in
Table :ref:`table:codonalign-meta-attributes`.

.. table:: Meta-attributes of CodonAligner objects.
   :name: table:codonalign-meta-attributes

   +----------------------------+---------------------------------+
   | Meta-attribute             | Attributes it maps to           |
   +============================+=================================+
   | ``frameshift_minus_score`` | ``frameshift_minus_two_score``, |
   |                            | ``frameshift_minus_one_score``  |
   +----------------------------+---------------------------------+
   | ``frameshift_plus_score``  | ``frameshift_plus_two_score``,  |
   |                            | ``frameshift_plus_one_score``   |
   +----------------------------+---------------------------------+
   | ``frameshift_two_score``   | ``frameshift_minus_two_score``, |
   |                            | ``frameshift_plus_two_score``   |
   +----------------------------+---------------------------------+
   | ``frameshift_one_score``   | ``frameshift_minus_one_score``, |
   |                            | ``frameshift_plus_one_score``   |
   +----------------------------+---------------------------------+
   | ``frameshift_score``       | ``frameshift_minus_two_score``, |
   |                            | ``frameshift_minus_one_score``, |
   |                            | ``frameshift_plus_one_score``,  |
   |                            | ``frameshift_plus_two_score``   |
   +----------------------------+---------------------------------+

Now let’s consider two nucleotide sequences and the amino acid sequences
they encode:

.. cont-doctest

.. code:: pycon

   >>> from Bio.Seq import Seq
   >>> from Bio.SeqRecord import SeqRecord
   >>> nuc1 = Seq("TCAGGGACTGCGAGAACCAAGCTACTGCTGCTGCTGGCTGCGCTCTGCGCCGCAGGTGGGGCGCTGGAG")
   >>> rna1 = SeqRecord(nuc1, id="rna1")
   >>> nuc2 = Seq("TCAGGGACTTCGAGAACCAAGCGCTCCTGCTGCTGGCTGCGCTCGGCGCCGCAGGTGGAGCACTGGAG")
   >>> rna2 = SeqRecord(nuc2, id="rna2")
   >>> aa1 = Seq("SGTARTKLLLLLAALCAAGGALE")
   >>> aa2 = Seq("SGTSRTKRLLLLAALGAAGGALE")
   >>> pro1 = SeqRecord(aa1, id="pro1")
   >>> pro2 = SeqRecord(aa2, id="pro2")

While the two protein sequences both consist of 23 amino acids, the
first nucleotide sequence consists of :math:`3 \times 23 = 69`
nucleotides while the second nucleotide sequence tonsists of only 68
nucleotides:

.. cont-doctest

.. code:: pycon

   >>> len(pro1)
   23
   >>> len(pro2)
   23
   >>> len(rna1)
   69
   >>> len(rna2)
   68

This is due to a -1 frame shift event during translation of the second
nucleotide sequence. Use ``CodonAligner.align`` to align ``rna1`` to
``pro1``, and ``rna2`` to ``pro2``, returning an iterator of
``Alignment`` objects:

.. cont-doctest

.. code:: pycon

   >>> alignments1 = aligner.align(pro1, rna1)
   >>> len(alignments1)
   1
   >>> alignment1 = next(alignments1)
   >>> print(alignment1)
   pro1              0 S  G  T  A  R  T  K  L  L  L  L  L  A  A  L  C  A  A  G  G  
   rna1              0 TCAGGGACTGCGAGAACCAAGCTACTGCTGCTGCTGGCTGCGCTCTGCGCCGCAGGTGGG
   <BLANKLINE>
   pro1             20 A  L  E   23
   rna1             60 GCGCTGGAG 69
   <BLANKLINE>
   >>> alignment1.coordinates
   array([[ 0, 23],
          [ 0, 69]])
   >>> alignment1[0]
   'SGTARTKLLLLLAALCAAGGALE'
   >>> alignment1[1]
   'TCAGGGACTGCGAGAACCAAGCTACTGCTGCTGCTGGCTGCGCTCTGCGCCGCAGGTGGGGCGCTGGAG'
   >>> alignments2 = aligner.align(pro2, rna2)
   >>> len(alignments2)
   1
   >>> alignment2 = next(alignments2)
   >>> print(alignment2)
   pro2              0 S  G  T  S  R  T  K  R   8
   rna2              0 TCAGGGACTTCGAGAACCAAGCGC 24
   <BLANKLINE>
   pro2              8 L  L  L  L  A  A  L  G  A  A  G  G  A  L  E   23
   rna2             23 CTCCTGCTGCTGGCTGCGCTCGGCGCCGCAGGTGGAGCACTGGAG 68
   <BLANKLINE>
   >>> alignment2[0]
   'SGTSRTKRLLLLAALGAAGGALE'
   >>> alignment2[1]
   'TCAGGGACTTCGAGAACCAAGCGCCTCCTGCTGCTGGCTGCGCTCGGCGCCGCAGGTGGAGCACTGGAG'
   >>> alignment2.coordinates
   array([[ 0,  8,  8, 23],
          [ 0, 24, 23, 68]])

While ``alignment1`` is a continuous alignment of the 69 nucleotides to
the 23 amino acids, in ``alignment2`` we find a -1 frame shift after 24
nucleotides. As ``alignment2[1]`` contains the nucleotide sequence after
applying the -1 frame shift, it is one nucleotide longer than ``nuc2``
and can be translated directly, resulting in the amino acid sequence
``aa2``:

.. cont-doctest

.. code:: pycon

   >>> from Bio.Seq import translate
   >>> len(nuc2)
   68
   >>> len(alignment2[1])
   69
   >>> translate(alignment2[1])
   'SGTSRTKRLLLLAALGAAGGALE'
   >>> _ == aa2
   True

The alignment score is stored as an attribute on the ``alignments1`` and
``alignments2`` iterators, and on the individual alignments
``alignment1`` and ``alignment2``:

.. cont-doctest

.. code:: pycon

   >>> alignments1.score
   23.0
   >>> alignment1.score
   23.0
   >>> alignments2.score
   20.0
   >>> alignment2.score
   20.0

where the score of the ``rna1``-``pro1`` alignment is equal to the
number of aligned amino acids, and the score of the ``rna2``-``pro2``
alignment is 3 less due to the penalty for the frame shift. To calculate
the alignment score without calculating the alignment itself, the
``score`` method can be used:

.. cont-doctest

.. code:: pycon

   >>> score = aligner.score(pro1, rna1)
   >>> print(score)
   23.0
   >>> score = aligner.score(pro2, rna2)
   >>> print(score)
   20.0

.. _`sec:msa_codons`:

Generating a multiple sequence alignment of codon sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Suppose we have a third related amino acid sequence and its associated
nucleotide sequence:

.. cont-doctest

.. code:: pycon

   >>> aa3 = Seq("MGTALLLLLAALCAAGGALE")
   >>> pro3 = SeqRecord(aa3, id="pro3")
   >>> nuc3 = Seq("ATGGGAACCGCGCTGCTTTTGCTACTGGCCGCGCTCTGCGCCGCAGGTGGGGCCCTGGAG")
   >>> rna3 = SeqRecord(nuc3, id="rna3")
   >>> nuc3.translate() == aa3
   True

As above, we use the ``CodonAligner`` to align the nucleotide sequence
to the amino acid sequence:

.. cont-doctest

.. code:: pycon

   >>> alignments3 = aligner.align(pro3, rna3)
   >>> len(alignments3)
   1
   >>> alignment3 = next(alignments3)
   >>> print(alignment3)
   pro3              0 M  G  T  A  L  L  L  L  L  A  A  L  C  A  A  G  G  A  L  E  
   rna3              0 ATGGGAACCGCGCTGCTTTTGCTACTGGCCGCGCTCTGCGCCGCAGGTGGGGCCCTGGAG
   <BLANKLINE>
   pro3             20 
   rna3             60 
   <BLANKLINE>

The three amino acid sequences can be aligned to each other, for example
using ClustalW. Here, we create the alignment by hand:

.. cont-doctest

.. code:: pycon

   >>> import numpy as np
   >>> from Bio.Align import Alignment
   >>> sequences = [pro1, pro2, pro3]
   >>> protein_alignment = Alignment(
   ...     sequences, coordinates=np.array([[0, 4, 7, 23], [0, 4, 7, 23], [0, 4, 4, 20]])
   ... )
   >>> print(protein_alignment)
   pro1              0 SGTARTKLLLLLAALCAAGGALE 23
   pro2              0 SGTSRTKRLLLLAALGAAGGALE 23
   pro3              0 MGTA---LLLLLAALCAAGGALE 20
   <BLANKLINE>

Now we can use the ``mapall`` method on the protein alignment, with the
nucleotide-to-protein pairwise alignments as the argument, to obtain the
corresponding codon alignment:

.. cont-doctest

.. code:: pycon

   >>> codon_alignment = protein_alignment.mapall([alignment1, alignment2, alignment3])
   >>> print(codon_alignment)
   rna1              0 TCAGGGACTGCGAGAACCAAGCTA 24
   rna2              0 TCAGGGACTTCGAGAACCAAGCGC 24
   rna3              0 ATGGGAACCGCG---------CTG 15
   <BLANKLINE>
   rna1             24 CTGCTGCTGCTGGCTGCGCTCTGCGCCGCAGGTGGGGCGCTGGAG 69
   rna2             23 CTCCTGCTGCTGGCTGCGCTCGGCGCCGCAGGTGGAGCACTGGAG 68
   rna3             15 CTTTTGCTACTGGCCGCGCTCTGCGCCGCAGGTGGGGCCCTGGAG 60
   <BLANKLINE>

Analyzing a codon alignment
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Calculating the number of nonsynonymous and synonymous substitutions per site
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The most important application of a codon alignment is to estimate the
number of nonsynonymous substitutions per site (dN) and synonymous
substitutions per site (dS). These can be calculated by the
``calculate_dn_ds`` function in ``Bio.Align.analysis``. This function
takes a pairwise codon alignment and input, as well as the optional
arguments ``method`` specifying the calculation method, ``codon_table``
(defaulting to the Standard Code), the ratio ``k`` of the transition and
transversion rates, and ``cfreq`` to specify the equilibrium codon
frequency. Biopython currently supports three counting based methods
(``NG86``, ``LWL85``, ``YN00``) as well as the maximum likelihood method
(``ML``) to estimate dN and dS:

-  ``NG86``: Nei and Gojobori (1986) [Nei1986]_
   (default). With this method, you can also specify the ratio of the
   transition and transversion rates via the argument ``k``, defaulting
   to ``1.0``.

-  ``LWL85``: Li *et al.* (1985) [Li1985]_.

-  ``YN00``: Yang and Nielsen (2000) [Yang2000]_.

-  ``ML``: Goldman and Yang (1994) [Goldman1994]_. With
   this method, you can also specify the equilibrium codon frequency via
   the ``cfreq`` argument, with the following options:

   -  ``F1x4``: count the nucleotide frequency in the provided codon
      sequences, and use it to calculate the background codon frequency;

   -  ``F3x4``: (default) count the nucleotide frequency separately for
      the first, second, and third position in the provided codons, and
      use it to calculate the background codon frequency;

   -  ``F61``: count the frequency of codons from the provided codon
      sequences, with a pseudocount of 0.1.

The ``calculate_dN_dS`` method can be applied to a pairwise codon
alignment. In general, the different calculation methods will result in
slightly different estimates for dN and dS:

.. cont-doctest

.. code:: pycon

   >>> from Bio.Align import analysis
   >>> pairwise_codon_alignment = codon_alignment[:2]
   >>> print(pairwise_codon_alignment)
   rna1              0 TCAGGGACTGCGAGAACCAAGCTA 24
                     0 |||||||||.||||||||||||..
   rna2              0 TCAGGGACTTCGAGAACCAAGCGC 24
   <BLANKLINE>
   rna1             24 CTGCTGCTGCTGGCTGCGCTCTGCGCCGCAGGTGGGGCGCTGGAG 69
                    24 ||.||||||||||||||||||.|||||||||||||.||.|||||| 69
   rna2             23 CTCCTGCTGCTGGCTGCGCTCGGCGCCGCAGGTGGAGCACTGGAG 68
   <BLANKLINE>
   >>> dN, dS = analysis.calculate_dn_ds(pairwise_codon_alignment, method="NG86")
   >>> print(dN, dS)  # doctest: +ELLIPSIS
   0.067715... 0.201197...
   >>> dN, dS = analysis.calculate_dn_ds(pairwise_codon_alignment, method="LWL85")
   >>> print(dN, dS)  # doctest: +ELLIPSIS
   0.068728... 0.207551...

.. code:: pycon

   >>> dN, dS = analysis.calculate_dn_ds(pairwise_codon_alignment, method="YN00")
   >>> print(dN, dS)  # doctest: +ELLIPSIS
   0.081468... 0.127706...
   >>> dN, dS = analysis.calculate_dn_ds(pairwise_codon_alignment, method="ML")
   >>> print(dN, dS)  # doctest: +ELLIPSIS
   0.069475... 0.205754...

For a multiple alignment of codon sequences, you can calculate a matrix
of dN and dS values:

.. cont-doctest

.. code:: pycon

   >>> dN, dS = analysis.calculate_dn_ds_matrix(codon_alignment, method="NG86")
   >>> print(dN)
   rna1    0.000000
   rna2    0.067715    0.000000
   rna3    0.060204    0.145469    0.000000
       rna1    rna2    rna3
   >>> print(dS)
   rna1    0.000000
   rna2    0.201198    0.000000
   rna3    0.664268    0.798957    0.000000
       rna1    rna2    rna3

The objects ``dN`` and ``dS`` returned by ``calculate_dn_ds_matrix`` are
instances of the ``DistanceMatrix`` class in
``Bio.Phylo.TreeConstruction``. This function only takes ``codon_table``
as an optional argument.

From these two sequences, you can create a dN tree and a dS tree using
``Bio.Phylo.TreeConstruction``:

.. cont-doctest

.. code:: pycon

   >>> from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
   >>> dn_constructor = DistanceTreeConstructor()
   >>> ds_constructor = DistanceTreeConstructor()
   >>> dn_tree = dn_constructor.upgma(dN)
   >>> ds_tree = ds_constructor.upgma(dS)
   >>> print(type(dn_tree))
   <class 'Bio.Phylo.BaseTree.Tree'>
   >>> print(dn_tree)  # doctest: +ELLIPSIS
   Tree(rooted=True)
       Clade(branch_length=0, name='Inner2')
           Clade(branch_length=0.053296..., name='rna2')
           Clade(branch_length=0.023194..., name='Inner1')
               Clade(branch_length=0.0301021..., name='rna3')
               Clade(branch_length=0.0301021..., name='rna1')
   >>> print(ds_tree)  # doctest: +ELLIPSIS
   Tree(rooted=True)
       Clade(branch_length=0, name='Inner2')
           Clade(branch_length=0.365806..., name='rna3')
           Clade(branch_length=0.265207..., name='Inner1')
               Clade(branch_length=0.100598..., name='rna2')
               Clade(branch_length=0.100598..., name='rna1')

Performing the McDonald-Kreitman test
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The McDonald-Kreitman test assesses the amount of adaptive evolution by
comparing the within species synonymous substitutions and nonsynonymous
substitutions to the between species synonymous substitutions and
nonsynonymous substitutions to see if they are from the same
evolutionary process. The test requires gene sequences sampled from
different individuals of the same species. In the following example, we
will use Adh gene from fruit fly. The data includes 11 individuals from
*Drosophila melanogaster*, 4 individuals from *Drosophila simulans*, and
12 individuals from *Drosophila yakuba*. The protein alignment data and
the nucleotide sequences are available in the ``Tests/codonalign``
directory as the files ``adh.aln`` and ``drosophila.fasta``,
respectively, in the Biopython distribution. The function ``mktest`` in
``Bio.Align.analysis`` implements the Mcdonald-Kreitman test.

.. doctest ../Tests/codonalign lib:numpy

.. code:: pycon

   >>> from Bio import SeqIO
   >>> from Bio import Align
   >>> from Bio.Align import CodonAligner
   >>> from Bio.Align.analysis import mktest
   >>> aligner = CodonAligner()
   >>> nucleotide_records = SeqIO.index("drosophila.fasta", "fasta")
   >>> for nucleotide_record in nucleotide_records.values():
   ...     print(nucleotide_record.description)  # doctest: +ELLIPSIS
   ...
   gi|9097|emb|X57361.1| Drosophila simulans (individual c) ...
   gi|9099|emb|X57362.1| Drosophila simulans (individual d) ...
   gi|9101|emb|X57363.1| Drosophila simulans (individual e) ...
   gi|9103|emb|X57364.1| Drosophila simulans (individual f) ...
   gi|9217|emb|X57365.1| Drosophila yakuba (individual a) ...
   gi|9219|emb|X57366.1| Drosophila yakuba (individual b) ...
   gi|9221|emb|X57367.1| Drosophila yakuba (individual c) ...
   gi|9223|emb|X57368.1| Drosophila yakuba (individual d) ...
   gi|9225|emb|X57369.1| Drosophila yakuba (individual e) ...
   gi|9227|emb|X57370.1| Drosophila yakuba (individual f) ...
   gi|9229|emb|X57371.1| Drosophila yakuba (individual g) ...
   gi|9231|emb|X57372.1| Drosophila yakuba (individual h) ...
   gi|9233|emb|X57373.1| Drosophila yakuba (individual i) ...
   gi|9235|emb|X57374.1| Drosophila yakuba (individual j) ...
   gi|9237|emb|X57375.1| Drosophila yakuba (individual k) ...
   gi|9239|emb|X57376.1| Drosophila yakuba (individual l) ...
   gi|156879|gb|M17837.1|DROADHCK D.melanogaster (strain Ja-F) ...
   gi|156863|gb|M19547.1|DROADHCC D.melanogaster (strain Af-S) ...
   gi|156877|gb|M17836.1|DROADHCJ D.melanogaster (strain Af-F) ...
   gi|156875|gb|M17835.1|DROADHCI D.melanogaster (strain Wa-F) ...
   gi|156873|gb|M17834.1|DROADHCH D.melanogaster (strain Fr-F) ...
   gi|156871|gb|M17833.1|DROADHCG D.melanogaster (strain Fl-F) ...
   gi|156869|gb|M17832.1|DROADHCF D.melanogaster (strain Ja-S) ...
   gi|156867|gb|M17831.1|DROADHCE D.melanogaster (strain Fl-2S) ...
   gi|156865|gb|M17830.1|DROADHCD D.melanogaster (strain Fr-S) ...
   gi|156861|gb|M17828.1|DROADHCB D.melanogaster (strain Fl-1S) ...
   gi|156859|gb|M17827.1|DROADHCA D.melanogaster (strain Wa-S) ...
   >>> protein_alignment = Align.read("adh.aln", "clustal")
   >>> len(protein_alignment)
   27
   >>> print(protein_alignment)  # doctest: +ELLIPSIS
   gi|9217|e         0 MAFTLTNKNVVFVAGLGGIGLDTSKELVKRDLKNLVILDRIENPAAIAELKAINPKVTVT
   gi|9219|e         0 MAFTLTNKNVVFVAGLGGIGLDTSKELVKRDLKNLVILDRIENPAAIAELKAINPKVTVT
   gi|9221|e         0 MAFTLTNKNVVFVAGLGGIGLDTSKELVKRDLKNLVILDRIENPAAIAELKAINPKVTVT
   ...
   gi|156859         0 MSFTLTNKNVIFVAGLGGIGLDTSKELLKRDLKNLVILDRIENPAAIAELKAINPKVTVT
   <BLANKLINE>
   ...
   <BLANKLINE>
   gi|9217|e       240 GTLEAIQWSKHWDSGI 256
   gi|9219|e       240 GTLEAIQWSKHWDSGI 256
   gi|9221|e       240 GTLEAIQWSKHWDSGI 256
   ...
   gi|156859       240 GTLEAIQWTKHWDSGI 256
   <BLANKLINE>
   >>> codon_alignments = []
   >>> for protein_record in protein_alignment.sequences:
   ...     nucleotide_record = nucleotide_records[protein_record.id]
   ...     alignments = aligner.align(protein_record, nucleotide_record)
   ...     assert len(alignments) == 1
   ...     codon_alignment = next(alignments)
   ...     codon_alignments.append(codon_alignment)
   ...
   >>> print(codon_alignment)  # doctest: +ELLIPSIS
   gi|156859         0 M  S  F  T  L  T  N  K  N  V  I  F  V  A  G  L  G  G  I  G  
   gi|156859         0 ATGTCGTTTACTTTGACCAACAAGAACGTGATTTTCGTTGCCGGTCTGGGAGGCATTGGT
   <BLANKLINE>
   gi|156859        20 L  D  T  S  K  E  L  L  K  R  D  L  K  N  L  V  I  L  D  R  
   gi|156859        60 CTGGACACCAGCAAGGAGCTGCTCAAGCGCGATCTGAAGAACCTGGTGATCCTCGACCGC
   <BLANKLINE>
   ...
   <BLANKLINE>
   gi|156859       240 G  T  L  E  A  I  Q  W  T  K  H  W  D  S  G  I   256
   gi|156859       720 GGCACCCTGGAGGCCATCCAGTGGACCAAGCACTGGGACTCCGGCATC 768
   <BLANKLINE>
   >>> nucleotide_records.close()  # Close indexed FASTA file
   >>> alignment = protein_alignment.mapall(codon_alignments)
   >>> print(alignment)  # doctest: +ELLIPSIS
   gi|9217|e         0 ATGGCGTTTACCTTGACCAACAAGAACGTGGTTTTCGTGGCCGGTCTGGGAGGCATTGGT
   gi|9219|e         0 ATGGCGTTTACCTTGACCAACAAGAACGTGGTTTTCGTGGCCGGTCTGGGAGGCATTGGT
   gi|9221|e         0 ATGGCGTTTACCTTGACCAACAAGAACGTGGTTTTCGTGGCCGGTCTGGGAGGCATTGGT
   ...
   gi|156859         0 ATGTCGTTTACTTTGACCAACAAGAACGTGATTTTCGTTGCCGGTCTGGGAGGCATTGGT
   <BLANKLINE>
   ...
   <BLANKLINE>
   gi|9217|e       720 GGCACCCTGGAGGCCATCCAGTGGTCCAAGCACTGGGACTCCGGCATC 768
   gi|9219|e       720 GGCACCCTGGAGGCCATCCAGTGGTCCAAGCACTGGGACTCCGGCATC 768
   gi|9221|e       720 GGTACCCTGGAGGCCATCCAGTGGTCCAAGCACTGGGACTCCGGCATC 768
   ...
   gi|156859       720 GGCACCCTGGAGGCCATCCAGTGGACCAAGCACTGGGACTCCGGCATC 768
   <BLANKLINE>
   >>> unique_species = ["Drosophila simulans", "Drosophila yakuba", "D.melanogaster"]
   >>> species = []
   >>> for record in alignment.sequences:
   ...     description = record.description
   ...     for s in unique_species:
   ...         if s in description:
   ...             break
   ...     else:
   ...         raise Exception(f"Failed to find species for {description}")
   ...     species.append(s)
   ...
   >>> print(species)
   ['Drosophila yakuba', 'Drosophila yakuba', 'Drosophila yakuba', 'Drosophila yakuba', 'Drosophila yakuba', 'Drosophila yakuba', 'Drosophila yakuba', 'Drosophila yakuba', 'Drosophila yakuba', 'Drosophila yakuba', 'Drosophila yakuba', 'Drosophila yakuba', 'Drosophila simulans', 'Drosophila simulans', 'Drosophila simulans', 'Drosophila simulans', 'D.melanogaster', 'D.melanogaster', 'D.melanogaster', 'D.melanogaster', 'D.melanogaster', 'D.melanogaster', 'D.melanogaster', 'D.melanogaster', 'D.melanogaster', 'D.melanogaster', 'D.melanogaster']
   >>> pvalue = mktest(alignment, species)
   >>> print(pvalue)  # doctest: +ELLIPSIS
   0.00206457...

In addition to the multiple codon alignment, the function ``mktest``
takes as input the species to which each sequence in the alignment
belongs to. The codon table can be provided as an optional argument
``codon_table``.
