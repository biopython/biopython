.. _`chapter:motifs`:

Sequence motif analysis using Bio.motifs
========================================

This chapter gives an overview of the functionality of the
``Bio.motifs`` package included in Biopython. It is intended for people
who are involved in the analysis of sequence motifs, so I’ll assume that
you are familiar with basic notions of motif analysis. In case something
is unclear, please look at Section :ref:`sec:links` for some
relevant links.

Most of this chapter describes the new ``Bio.motifs`` package included
in Biopython 1.61 onwards, which is replacing the older ``Bio.Motif``
package introduced with Biopython 1.50, which was in turn based on two
older former Biopython modules, ``Bio.AlignAce`` and ``Bio.MEME``. It
provides most of their functionality with a unified motif object
implementation.

Speaking of other libraries, if you are reading this you might be
interested in `TAMO <http://fraenkel-nsf.csbi.mit.edu/TAMO/>`__, another
python library designed to deal with sequence motifs. It supports more
*de-novo* motif finders, but it is not a part of Biopython and has some
restrictions on commercial use.

.. _`sec:motif_object`:

Motif objects
-------------

Since we are interested in motif analysis, we need to take a look at
``Motif`` objects in the first place. For that we need to import the
Bio.motifs library:

.. doctest ../Tests/motifs

.. code:: pycon

   >>> from Bio import motifs

and we can start creating our first motif objects. We can either create
a ``Motif`` object from a list of instances of the motif, or we can
obtain a ``Motif`` object by parsing a file from a motif database or
motif finding software.

.. _`subsec:creating_motif`:

Creating a motif from instances
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Suppose we have these instances of a DNA motif:

.. cont-doctest

.. code:: pycon

   >>> from Bio.Seq import Seq
   >>> instances = [
   ...     Seq("TACAA"),
   ...     Seq("TACGC"),
   ...     Seq("TACAC"),
   ...     Seq("TACCC"),
   ...     Seq("AACCC"),
   ...     Seq("AATGC"),
   ...     Seq("AATGC"),
   ... ]

then we can create a Motif object as follows:

.. cont-doctest

.. code:: pycon

   >>> m = motifs.create(instances)

The instances from which this motif was created is stored in the
``.alignment`` property:

.. cont-doctest

.. code:: pycon

   >>> print(m.alignment.sequences)
   [Seq('TACAA'), Seq('TACGC'), Seq('TACAC'), Seq('TACCC'), Seq('AACCC'), Seq('AATGC'), Seq('AATGC')]

Printing the Motif object shows the instances from which it was
constructed:

.. cont-doctest

.. code:: pycon

   >>> print(m)
   TACAA
   TACGC
   TACAC
   TACCC
   AACCC
   AATGC
   AATGC

The length of the motif is defined as the sequence length, which should
be the same for all instances:

.. cont-doctest

.. code:: pycon

   >>> len(m)
   5

The Motif object has an attribute ``.counts`` containing the counts of
each nucleotide at each position. Printing this counts matrix shows it
in an easily readable format:

.. cont-doctest

.. code:: pycon

   >>> print(m.counts)
           0      1      2      3      4
   A:   3.00   7.00   0.00   2.00   1.00
   C:   0.00   0.00   5.00   2.00   6.00
   G:   0.00   0.00   0.00   3.00   0.00
   T:   4.00   0.00   2.00   0.00   0.00
   <BLANKLINE>

You can access these counts as a dictionary:

.. cont-doctest

.. code:: pycon

   >>> m.counts["A"]
   [3.0, 7.0, 0.0, 2.0, 1.0]

but you can also think of it as a 2D array with the nucleotide as the
first dimension and the position as the second dimension:

.. cont-doctest

.. code:: pycon

   >>> m.counts["T", 0]
   4.0
   >>> m.counts["T", 2]
   2.0
   >>> m.counts["T", 3]
   0.0

You can also directly access columns of the counts matrix

.. code:: pycon

   >>> m.counts[:, 3]
   {'A': 2.0, 'C': 2.0, 'T': 0.0, 'G': 3.0}

Instead of the nucleotide itself, you can also use the index of the
nucleotide in the alphabet of the motif:

.. cont-doctest

.. code:: pycon

   >>> m.alphabet
   'ACGT'
   >>> m.counts["A", :]
   (3.0, 7.0, 0.0, 2.0, 1.0)
   >>> m.counts[0, :]
   (3.0, 7.0, 0.0, 2.0, 1.0)

.. _`sec:motif_consensus`:

Obtaining a consensus sequence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The consensus sequence of a motif is defined as the sequence of letters
along the positions of the motif for which the largest value in the
corresponding columns of the ``.counts`` matrix is obtained:

.. cont-doctest

.. code:: pycon

   >>> m.consensus
   Seq('TACGC')

Conversely, the anticonsensus sequence corresponds to the smallest
values in the columns of the ``.counts`` matrix:

.. cont-doctest

.. code:: pycon

   >>> m.anticonsensus
   Seq('CCATG')

Note that there is some ambiguity in the definition of the consensus and
anticonsensus sequence if in some columns multiple nucleotides have the
maximum or minimum count.

For DNA sequences, you can also ask for a degenerate consensus sequence,
in which ambiguous nucleotides are used for positions where there are
multiple nucleotides with high counts:

.. cont-doctest

.. code:: pycon

   >>> m.degenerate_consensus
   Seq('WACVC')

Here, W and R follow the IUPAC nucleotide ambiguity codes: W is either A
or T, and V is A, C, or G [Cornish1985]_. The
degenerate consensus sequence is constructed following the rules
specified by Cavener [Cavener1987]_.

The ``motif.counts.calculate_consensus`` method lets you specify in
detail how the consensus sequence should be calculated. This method
largely follows the conventions of the EMBOSS program ``cons``, and
takes the following arguments:

substitution_matrix
   The scoring matrix used when comparing sequences. By default, it is
   ``None``, in which case we simply count the frequency of each letter.
   Instead of the default value, you can use the substitution matrices
   available in ``Bio.Align.substitution\_matrices``. Common choices are
   BLOSUM62 (also known as EBLOSUM62) for protein, and NUC.4.4 (also
   known as EDNAFULL) for nucleotides. NOTE: Currently, this method has
   not yet been implemented for values other than the default value
   ``None``.

plurality
   Threshold value for the number of positive matches, divided by the
   total count in a column, required to reach consensus. If
   ``substitution_matrix`` is ``None``, then this argument must also be
   ``None``, and is ignored; a ``ValueError`` is raised otherwise. If
   ``substitution_matrix`` is not ``None``, then the default value of
   the plurality is 0.5.

identity
   Number of identities, divided by the total count in a column,
   required to define a consensus value. If the number of identities is
   less than identity multiplied by the total count in a column, then
   the undefined character (``N`` for nucleotides and ``X`` for amino
   acid sequences) is used in the consensus sequence. If ``identity`` is
   1.0, then only columns of identical letters contribute to the
   consensus. Default value is zero.

setcase
   threshold for the positive matches, divided by the total count in a
   column, above which the consensus is is upper-case and below which
   the consensus is in lower-case. By default, this is equal to 0.5.

This is an example:

.. cont-doctest

.. code:: pycon

   >>> m.counts.calculate_consensus(identity=0.5, setcase=0.7)
   'tACNC'

Reverse-complementing a motif
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We can get the reverse complement of a motif by calling the
``reverse_complement`` method on it:

.. cont-doctest

.. code:: pycon

   >>> r = m.reverse_complement()
   >>> r.consensus
   Seq('GCGTA')
   >>> r.degenerate_consensus
   Seq('GBGTW')
   >>> print(r)
   TTGTA
   GCGTA
   GTGTA
   GGGTA
   GGGTT
   GCATT
   GCATT

The reverse complement is only defined for DNA motifs.

Slicing a motif
~~~~~~~~~~~~~~~

You can slice the motif to obtain a new ``Motif`` object for the
selected positions:

.. cont-doctest

.. code:: pycon

   >>> m_sub = m[2:-1]
   >>> print(m_sub)
   CA
   CG
   CA
   CC
   CC
   TG
   TG
   >>> m_sub.consensus
   Seq('CG')
   >>> m_sub.degenerate_consensus
   Seq('CV')

.. _`subsec:relative_entropy`:

Relative entropy
~~~~~~~~~~~~~~~~

The relative entropy (or Kullback-Leibler distance) :math:`H_j` of
column :math:`j` of the motif is defined as in
[Schneider1986]_ [Durbin1998]_:

.. math:: H_{j} = \sum_{i=1}^{M} p_{ij} \log\left(\frac{p_{ij}}{b_{i}}\right)

where:

-  :math:`M` – The number of letters in the alphabet (given by
   ``len(m.alphabet)``);

-  :math:`p_{ij}` – The observed frequency of letter :math:`i`,
   normalized, in the :math:`j`-th column (see below);

-  :math:`b_{i}` – The background probability of letter :math:`i` (given
   by ``m.background[i]``).

The observed frequency :math:`p_{ij}` is computed as follows:

.. math:: p_{ij} = \frac{c_{ij} + k_i}{C_{j} + k}

where:

-  :math:`c_{ij}` – the number of times letter :math:`i` appears in
   column :math:`j` of the alignment (given by ``m.counts[i, j]``);

-  :math:`C_{j}` – The total number of letters in column :math:`j`:
   :math:`C_{j} = \sum_{i=1}^{M} c_{ij}` (given by
   ``sum(m.counts[:, j])``).

-  :math:`k_i` – the pseudocount of letter :math:`i` (given by
   ``m.pseudocounts[i]``).

-  :math:`k` – the total pseudocount: :math:`k = \sum_{i=1}^{M} k_i`
   (given by ``sum(m.pseudocounts.values())``).

With these definitions, both :math:`p_{ij}` and :math:`b_{i}` are
normalized to 1:

.. math:: \sum_{i=1}^{M} p_{ij} = 1

.. math:: \sum_{i=1}^{M} b_i = 1

The relative entropy is the same as the information content if the
background distribution is uniform.

The relative entropy for each column of motif ``m`` can be obtained
using the ``relative_entropy`` property:

.. code:: pycon

   >>> m.relative_entropy
   array([1.01477186, 2.        , 1.13687943, 0.44334329, 1.40832722])

These values are calculated using the base-2 logarithm, and are
therefore in units of bits. The second column (which consists of ``A``
nucleotides only) has the highest relative entropy; the fourth column
(which consists of ``A``, ``C``, or ``G`` nucleotides) has the lowest
relative entropy). The relative entropy of the motif can be calculated
by summing over the columns:

.. cont-doctest

.. code:: pycon

   >>> print(f"Relative entropy is {sum(m.relative_entropy):0.5f}")
   Relative entropy is 6.00332

Creating a sequence logo
~~~~~~~~~~~~~~~~~~~~~~~~

If we have internet access, we can create a
`weblogo <https://weblogo.berkeley.edu>`__:

.. code:: pycon

   >>> m.weblogo("mymotif.png")

We should get our logo saved as a PNG in the specified file.

.. _`sec:io`:

Reading motifs
--------------

Creating motifs from instances by hand is a bit boring, so it’s useful
to have some I/O functions for reading and writing motifs. There are not
any really well established standards for storing motifs, but there are
a couple of formats that are more used than others.

JASPAR
~~~~~~

One of the most popular motif databases is
`JASPAR <http://jaspar.genereg.net>`__. In addition to the motif
sequence information, the JASPAR database stores a lot of
meta-information for each motif. The module ``Bio.motifs`` contains a
specialized class ``jaspar.Motif`` in which this meta-information is
represented as attributes:

-  ``matrix_id`` - the unique JASPAR motif ID, e.g. ’MA0004.1’

-  ``name`` - the name of the TF, e.g. ’Arnt’

-  ``collection`` - the JASPAR collection to which the motif belongs,
   e.g. ’CORE’

-  ``tf_class`` - the structural class of this TF, e.g. ’Zipper-Type’

-  ``tf_family`` - the family to which this TF belongs, e.g.
   ’Helix-Loop-Helix’

-  ``species`` - the species to which this TF belongs, may have multiple
   values, these are specified as taxonomy IDs, e.g. 10090

-  ``tax_group`` - the taxonomic supergroup to which this motif belongs,
   e.g. ’vertebrates’

-  ``acc`` - the accession number of the TF protein, e.g. ’P53762’

-  ``data_type`` - the type of data used to construct this motif, e.g.
   ’SELEX’

-  ``medline`` - the Pubmed ID of literature supporting this motif, may
   be multiple values, e.g. 7592839

-  ``pazar_id`` - external reference to the TF in the PAZAR database,
   e.g. ’TF0000003’

-  ``comment`` - free form text containing notes about the construction
   of the motif

The ``jaspar.Motif`` class inherits from the generic ``Motif`` class and
therefore provides all the facilities of any of the motif formats —
reading motifs, writing motifs, scanning sequences for motif instances
etc.

JASPAR stores motifs in several different ways including three different
flat file formats and as an SQL database. All of these formats
facilitate the construction of a counts matrix. However, the amount of
meta information described above that is available varies with the
format.

The JASPAR ``sites`` format
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The first of the three flat file formats contains a list of instances.
As an example, these are the beginning and ending lines of the JASPAR
``Arnt.sites`` file showing known binding sites of the mouse
helix-loop-helix transcription factor Arnt.

.. code:: text

   >MA0004 ARNT 1
   CACGTGatgtcctc
   >MA0004 ARNT 2
   CACGTGggaggtac
   >MA0004 ARNT 3
   CACGTGccgcgcgc
   ...
   >MA0004 ARNT 18
   AACGTGacagccctcc
   >MA0004 ARNT 19
   AACGTGcacatcgtcc
   >MA0004 ARNT 20
   aggaatCGCGTGc

The parts of the sequence in capital letters are the motif instances
that were found to align to each other.

We can create a ``Motif`` object from these instances as follows:

.. cont-doctest

.. code:: pycon

   >>> from Bio import motifs
   >>> with open("Arnt.sites") as handle:
   ...     arnt = motifs.read(handle, "sites")
   ...

The instances from which this motif was created is stored in the
``.alignment`` property:

.. cont-doctest

.. code:: pycon

   >>> print(arnt.alignment.sequences[:3])
   [Seq('CACGTG'), Seq('CACGTG'), Seq('CACGTG')]
   >>> for sequence in arnt.alignment.sequences:
   ...     print(sequence)
   ...
   CACGTG
   CACGTG
   CACGTG
   CACGTG
   CACGTG
   CACGTG
   CACGTG
   CACGTG
   CACGTG
   CACGTG
   CACGTG
   CACGTG
   CACGTG
   CACGTG
   CACGTG
   AACGTG
   AACGTG
   AACGTG
   AACGTG
   CGCGTG

The counts matrix of this motif is automatically calculated from the
instances:

.. cont-doctest

.. code:: pycon

   >>> print(arnt.counts)
           0      1      2      3      4      5
   A:   4.00  19.00   0.00   0.00   0.00   0.00
   C:  16.00   0.00  20.00   0.00   0.00   0.00
   G:   0.00   1.00   0.00  20.00   0.00  20.00
   T:   0.00   0.00   0.00   0.00  20.00   0.00
   <BLANKLINE>

This format does not store any meta information.

The JASPAR ``pfm`` format
^^^^^^^^^^^^^^^^^^^^^^^^^

JASPAR also makes motifs available directly as a count matrix, without
the instances from which it was created. This ``pfm`` format only stores
the counts matrix for a single motif. For example, this is the JASPAR
file ``SRF.pfm`` containing the counts matrix for the human SRF
transcription factor:

.. code:: text

    2 9 0 1 32 3 46 1 43 15 2 2
    1 33 45 45 1 1 0 0 0 1 0 1
   39 2 1 0 0 0 0 0 0 0 44 43
    4 2 0 0 13 42 0 45 3 30 0 0

We can create a motif for this count matrix as follows:

.. cont-doctest

.. code:: pycon

   >>> with open("SRF.pfm") as handle:
   ...     srf = motifs.read(handle, "pfm")
   ...
   >>> print(srf.counts)
           0      1      2      3      4      5      6      7      8      9     10     11
   A:   2.00   9.00   0.00   1.00  32.00   3.00  46.00   1.00  43.00  15.00   2.00   2.00
   C:   1.00  33.00  45.00  45.00   1.00   1.00   0.00   0.00   0.00   1.00   0.00   1.00
   G:  39.00   2.00   1.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00  44.00  43.00
   T:   4.00   2.00   0.00   0.00  13.00  42.00   0.00  45.00   3.00  30.00   0.00   0.00
   <BLANKLINE>

As this motif was created from the counts matrix directly, it has no
instances associated with it:

.. cont-doctest

.. code:: pycon

   >>> print(srf.alignment)
   None

We can now ask for the consensus sequence of these two motifs:

.. cont-doctest

.. code:: pycon

   >>> print(arnt.counts.consensus)
   CACGTG
   >>> print(srf.counts.consensus)
   GCCCATATATGG

As with the instances file, no meta information is stored in this
format.

The JASPAR format ``jaspar``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``jaspar`` file format allows multiple motifs to be specified in a
single file. In this format each of the motif records consist of a
header line followed by four lines defining the counts matrix. The
header line begins with a ``>`` character (similar to the Fasta file
format) and is followed by the unique JASPAR matrix ID and the TF name.
The following example shows a ``jaspar`` formatted file containing the
three motifs Arnt, RUNX1 and MEF2A:

.. code:: text

   >MA0004.1 Arnt
   A  [ 4 19  0  0  0  0 ]
   C  [16  0 20  0  0  0 ]
   G  [ 0  1  0 20  0 20 ]
   T  [ 0  0  0  0 20  0 ]
   >MA0002.1 RUNX1
   A  [10 12  4  1  2  2  0  0  0  8 13 ]
   C  [ 2  2  7  1  0  8  0  0  1  2  2 ]
   G  [ 3  1  1  0 23  0 26 26  0  0  4 ]
   T  [11 11 14 24  1 16  0  0 25 16  7 ]
   >MA0052.1 MEF2A
   A  [ 1  0 57  2  9  6 37  2 56  6 ]
   C  [50  0  1  1  0  0  0  0  0  0 ]
   G  [ 0  0  0  0  0  0  0  0  2 50 ]
   T  [ 7 58  0 55 49 52 21 56  0  2 ]

The motifs are read as follows:

.. code:: pycon

   >>> fh = open("jaspar_motifs.txt")
   >>> for m in motifs.parse(fh, "jaspar"):
   ...     print(m)
   ...
   TF name  Arnt
   Matrix ID   MA0004.1
   Matrix:
           0      1      2      3      4      5
   A:   4.00  19.00   0.00   0.00   0.00   0.00
   C:  16.00   0.00  20.00   0.00   0.00   0.00
   G:   0.00   1.00   0.00  20.00   0.00  20.00
   T:   0.00   0.00   0.00   0.00  20.00   0.00



   TF name  RUNX1
   Matrix ID   MA0002.1
   Matrix:
           0      1      2      3      4      5      6      7      8      9     10
   A:  10.00  12.00   4.00   1.00   2.00   2.00   0.00   0.00   0.00   8.00  13.00
   C:   2.00   2.00   7.00   1.00   0.00   8.00   0.00   0.00   1.00   2.00   2.00
   G:   3.00   1.00   1.00   0.00  23.00   0.00  26.00  26.00   0.00   0.00   4.00
   T:  11.00  11.00  14.00  24.00   1.00  16.00   0.00   0.00  25.00  16.00   7.00



   TF name  MEF2A
   Matrix ID   MA0052.1
   Matrix:
           0      1      2      3      4      5      6      7      8      9
   A:   1.00   0.00  57.00   2.00   9.00   6.00  37.00   2.00  56.00   6.00
   C:  50.00   0.00   1.00   1.00   0.00   0.00   0.00   0.00   0.00   0.00
   G:   0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   2.00  50.00
   T:   7.00  58.00   0.00  55.00  49.00  52.00  21.00  56.00   0.00   2.00

Note that printing a JASPAR motif yields both the counts data and the
available meta-information.

Accessing the JASPAR database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In addition to parsing these flat file formats, we can also retrieve
motifs from a JASPAR SQL database. Unlike the flat file formats, a
JASPAR database allows storing of all possible meta information defined
in the JASPAR ``Motif`` class. It is beyond the scope of this document
to describe how to set up a JASPAR database (please see the main
`JASPAR <http://jaspar.genereg.net>`__ website). Motifs are read from a
JASPAR database using the ``Bio.motifs.jaspar.db`` module. First connect
to the JASPAR database using the JASPAR5 class which models the the
latest JASPAR schema:

.. code:: pycon

   >>> from Bio.motifs.jaspar.db import JASPAR5
   >>>
   >>> JASPAR_DB_HOST = "yourhostname"  # fill in these values
   >>> JASPAR_DB_NAME = "yourdatabase"
   >>> JASPAR_DB_USER = "yourusername"
   >>> JASPAR_DB_PASS = "yourpassword"
   >>>
   >>> jdb = JASPAR5(
   ...     host=JASPAR_DB_HOST,
   ...     name=JASPAR_DB_NAME,
   ...     user=JASPAR_DB_USER,
   ...     password=JASPAR_DB_PASS,
   ... )

Now we can fetch a single motif by its unique JASPAR ID with the
``fetch_motif_by_id`` method. Note that a JASPAR ID consists of a base
ID and a version number separated by a decimal point, e.g. ’MA0004.1’.
The ``fetch_motif_by_id`` method allows you to use either the fully
specified ID or just the base ID. If only the base ID is provided, the
latest version of the motif is returned.

.. code:: pycon

   >>> arnt = jdb.fetch_motif_by_id("MA0004")

Printing the motif reveals that the JASPAR SQL database stores much more
meta-information than the flat files:

.. code:: pycon

   >>> print(arnt)
   TF name Arnt
   Matrix ID   MA0004.1
   Collection  CORE
   TF class    Zipper-Type
   TF family   Helix-Loop-Helix
   Species 10090
   Taxonomic group vertebrates
   Accession   ['P53762']
   Data type used  SELEX
   Medline 7592839
   PAZAR ID    TF0000003
   Comments    -
   Matrix:
       0      1      2      3      4      5
   A:   4.00  19.00   0.00   0.00   0.00   0.00
   C:  16.00   0.00  20.00   0.00   0.00   0.00
   G:   0.00   1.00   0.00  20.00   0.00  20.00
   T:   0.00   0.00   0.00   0.00  20.00   0.00

We can also fetch motifs by name. The name must be an exact match
(partial matches or database wildcards are not currently supported).
Note that as the name is not guaranteed to be unique, the
``fetch_motifs_by_name`` method actually returns a list.

.. code:: pycon

   >>> motifs = jdb.fetch_motifs_by_name("Arnt")
   >>> print(motifs[0])
   TF name Arnt
   Matrix ID   MA0004.1
   Collection  CORE
   TF class    Zipper-Type
   TF family   Helix-Loop-Helix
   Species 10090
   Taxonomic group vertebrates
   Accession   ['P53762']
   Data type used  SELEX
   Medline 7592839
   PAZAR ID    TF0000003
   Comments    -
   Matrix:
       0      1      2      3      4      5
   A:   4.00  19.00   0.00   0.00   0.00   0.00
   C:  16.00   0.00  20.00   0.00   0.00   0.00
   G:   0.00   1.00   0.00  20.00   0.00  20.00
   T:   0.00   0.00   0.00   0.00  20.00   0.00

The ``fetch_motifs`` method allows you to fetch motifs which match a
specified set of criteria. These criteria include any of the above
described meta information as well as certain matrix properties such as
the minimum information content (``min_ic`` in the example below), the
minimum length of the matrix or the minimum number of sites used to
construct the matrix. Only motifs which pass ALL the specified criteria
are returned. Note that selection criteria which correspond to meta
information which allow for multiple values may be specified as either a
single value or a list of values, e.g. ``tax_group`` and ``tf_family``
in the example below.

.. code:: pycon

   >>> motifs = jdb.fetch_motifs(
   ...     collection="CORE",
   ...     tax_group=["vertebrates", "insects"],
   ...     tf_class="Winged Helix-Turn-Helix",
   ...     tf_family=["Forkhead", "Ets"],
   ...     min_ic=12,
   ... )
   >>> for motif in motifs:
   ...     pass  # do something with the motif
   ...

Compatibility with Perl TFBS modules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An important thing to note is that the JASPAR ``Motif`` class was
designed to be compatible with the popular `Perl TFBS
modules <http://tfbs.genereg.net/>`__. Therefore some specifics about
the choice of defaults for background and pseudocounts as well as how
information content is computed and sequences searched for instances is
based on this compatibility criteria. These choices are noted in the
specific subsections below.

-  | **Choice of background:**
   | The Perl ``TFBS`` modules appear to allow a choice of custom
     background probabilities (although the documentation states that
     uniform background is assumed). However the default is to use a
     uniform background. Therefore it is recommended that you use a
     uniform background for computing the position-specific scoring
     matrix (PSSM). This is the default when using the Biopython
     ``motifs`` module.

-  | **Choice of pseudocounts:**
   | By default, the Perl ``TFBS`` modules use a pseudocount equal to
     :math:`\sqrt{N} * \textrm{bg}[\textrm{nucleotide}]`, where
     :math:`N` represents the total number of sequences used to
     construct the matrix. To apply this same pseudocount formula, set
     the motif ``pseudocounts`` attribute using the
     ``jaspar.calculate\_pseudcounts()`` function:

   .. code:: pycon

      >>> motif.pseudocounts = motifs.jaspar.calculate_pseudocounts(motif)

   Note that it is possible for the counts matrix to have an unequal
   number of sequences making up the columns. The pseudocount
   computation uses the average number of sequences making up the
   matrix. However, when ``normalize`` is called on the counts matrix,
   each count value in a column is divided by the total number of
   sequences making up that specific column, not by the average number
   of sequences. This differs from the Perl ``TFBS`` modules because the
   normalization is not done as a separate step and so the average
   number of sequences is used throughout the computation of the pssm.
   Therefore, for matrices with unequal column counts, the PSSM computed
   by the ``motifs`` module will differ somewhat from the pssm computed
   by the Perl ``TFBS`` modules.

-  | **Computation of matrix information content:**
   | The information content (IC) or specificity of a matrix is computed
     using the ``mean`` method of the ``PositionSpecificScoringMatrix``
     class. However of note, in the Perl ``TFBS`` modules the default
     behavior is to compute the IC without first applying pseudocounts,
     even though by default the PSSMs are computed using pseudocounts as
     described above.

-  | **Searching for instances:**
   | Searching for instances with the Perl ``TFBS`` motifs was usually
     performed using a relative score threshold, i.e. a score in the
     range 0 to 1. In order to compute the absolute PSSM score
     corresponding to a relative score one can use the equation:

   .. code:: pycon

      >>> abs_score = (pssm.max - pssm.min) * rel_score + pssm.min

   To convert the absolute score of an instance back to a relative
   score, one can use the equation:

   .. code:: pycon

      >>> rel_score = (abs_score - pssm.min) / (pssm.max - pssm.min)

   For example, using the Arnt motif before, let’s search a sequence
   with a relative score threshold of 0.8.

   .. code:: pycon

      >>> test_seq = Seq("TAAGCGTGCACGCGCAACACGTGCATTA")
      >>> arnt.pseudocounts = motifs.jaspar.calculate_pseudocounts(arnt)
      >>> pssm = arnt.pssm
      >>> max_score = pssm.max
      >>> min_score = pssm.min
      >>> abs_score_threshold = (max_score - min_score) * 0.8 + min_score
      >>> for pos, score in pssm.search(test_seq, threshold=abs_score_threshold):
      ...     rel_score = (score - min_score) / (max_score - min_score)
      ...     print(f"Position {pos}: score = {score:5.3f}, rel. score = {rel_score:5.3f}")
      ...
      Position 2: score = 5.362, rel. score = 0.801
      Position 8: score = 6.112, rel. score = 0.831
      Position -20: score = 7.103, rel. score = 0.870
      Position 17: score = 10.351, rel. score = 1.000
      Position -11: score = 10.351, rel. score = 1.000

MEME
~~~~

MEME [Bailey1994]_ is a tool for discovering motifs in
a group of related DNA or protein sequences. It takes as input a group
of DNA or protein sequences and outputs as many motifs as requested.
Therefore, in contrast to JASPAR files, MEME output files typically
contain multiple motifs. This is an example.

At the top of an output file generated by MEME shows some background
information about the MEME and the version of MEME used:

.. code:: text

   ********************************************************************************
   MEME - Motif discovery tool
   ********************************************************************************
   MEME version 3.0 (Release date: 2004/08/18 09:07:01)
   ...

Further down, the input set of training sequences is recapitulated:

.. code:: text

   ********************************************************************************
   TRAINING SET
   ********************************************************************************
   DATAFILE= INO_up800.s
   ALPHABET= ACGT
   Sequence name            Weight Length  Sequence name            Weight Length
   -------------            ------ ------  -------------            ------ ------
   CHO1                     1.0000    800  CHO2                     1.0000    800
   FAS1                     1.0000    800  FAS2                     1.0000    800
   ACC1                     1.0000    800  INO1                     1.0000    800
   OPI3                     1.0000    800
   ********************************************************************************

and the exact command line that was used:

.. code:: text

   ********************************************************************************
   COMMAND LINE SUMMARY
   ********************************************************************************
   This information can also be useful in the event you wish to report a
   problem with the MEME software.

   command: meme -mod oops -dna -revcomp -nmotifs 2 -bfile yeast.nc.6.freq INO_up800.s
   ...

Next is detailed information on each motif that was found:

.. code:: text

   ********************************************************************************
   MOTIF  1        width =   12   sites =   7   llr = 95   E-value = 2.0e-001
   ********************************************************************************
   --------------------------------------------------------------------------------
           Motif 1 Description
   --------------------------------------------------------------------------------
   Simplified        A  :::9:a::::3:
   pos.-specific     C  ::a:9:11691a
   probability       G  ::::1::94:4:
   matrix            T  aa:1::9::11:

To parse this file (stored as ``meme.dna.oops.txt``), use

.. cont-doctest

.. code:: pycon

   >>> with open("meme.INO_up800.classic.oops.xml") as handle:
   ...     record = motifs.parse(handle, "meme")
   ...

The ``motifs.parse`` command reads the complete file directly, so you
can close the file after calling ``motifs.parse``. The header
information is stored in attributes:

.. cont-doctest

.. code:: pycon

   >>> record.version
   '5.0.1'
   >>> record.datafile
   'common/INO_up800.s'
   >>> record.command
   'meme common/INO_up800.s -oc results/meme10 -mod oops -dna -revcomp -bfile common/yeast.nc.6.freq -nmotifs 2 -objfun classic -minw 8 -nostatus '
   >>> record.alphabet
   'ACGT'
   >>> record.sequences
   ['sequence_0', 'sequence_1', 'sequence_2', 'sequence_3', 'sequence_4', 'sequence_5', 'sequence_6']

The record is an object of the ``Bio.motifs.meme.Record`` class. The
class inherits from list, and you can think of ``record`` as a list of
Motif objects:

.. cont-doctest

.. code:: pycon

   >>> len(record)
   2
   >>> motif = record[0]
   >>> print(motif.consensus)
   GCGGCATGTGAAA
   >>> print(motif.degenerate_consensus)
   GSKGCATGTGAAA

In addition to these generic motif attributes, each motif also stores
its specific information as calculated by MEME. For example,

.. cont-doctest

.. code:: pycon

   >>> motif.num_occurrences
   7
   >>> motif.length
   13
   >>> evalue = motif.evalue
   >>> print("%3.1g" % evalue)
   0.2
   >>> motif.name
   'GSKGCATGTGAAA'
   >>> motif.id
   'motif_1'

In addition to using an index into the record, as we did above, you can
also find it by its name:

.. cont-doctest

.. code:: pycon

   >>> motif = record["GSKGCATGTGAAA"]

Each motif has an attribute ``.alignment`` with the sequence alignment
in which the motif was found, providing some information on each of the
sequences:

.. cont-doctest

.. code:: pycon

   >>> len(motif.alignment)
   7
   >>> motif.alignment.sequences[0]
   Seq('GCGGCATGTGAAA')
   >>> motif.alignment.sequences[0].motif_name
   'GSKGCATGTGAAA'
   >>> motif.alignment.sequences[0].sequence_name
   'INO1'
   >>> motif.alignment.sequences[0].sequence_id
   'sequence_5'
   >>> motif.alignment.sequences[0].start
   620
   >>> motif.alignment.sequences[0].strand
   '+'
   >>> motif.alignment.sequences[0].length
   13
   >>> pvalue = motif.alignment.sequences[0].pvalue
   >>> print("%5.3g" % pvalue)
   1.21e-08

MAST
^^^^

TRANSFAC
~~~~~~~~

TRANSFAC is a manually curated database of transcription factors,
together with their genomic binding sites and DNA binding profiles
[Matys2003]_. While the file format used in the
TRANSFAC database is nowadays also used by others, we will refer to it
as the TRANSFAC file format.

A minimal file in the TRANSFAC format looks as follows:

.. code:: text

   ID  motif1
   P0      A      C      G      T
   01      1      2      2      0      S
   02      2      1      2      0      R
   03      3      0      1      1      A
   04      0      5      0      0      C
   05      5      0      0      0      A
   06      0      0      4      1      G
   07      0      1      4      0      G
   08      0      0      0      5      T
   09      0      0      5      0      G
   10      0      1      2      2      K
   11      0      2      0      3      Y
   12      1      0      3      1      G
   //

This file shows the frequency matrix of motif ``motif1`` of 12
nucleotides. In general, one file in the TRANSFAC format can contain
multiple motifs. For example, this is the contents of the example
TRANSFAC file ``transfac.dat``:

.. code:: text

   VV  EXAMPLE January 15, 2013
   XX
   //
   ID  motif1
   P0      A      C      G      T
   01      1      2      2      0      S
   02      2      1      2      0      R
   03      3      0      1      1      A
   ...
   11      0      2      0      3      Y
   12      1      0      3      1      G
   //
   ID  motif2
   P0      A      C      G      T
   01      2      1      2      0      R
   02      1      2      2      0      S
   ...
   09      0      0      0      5      T
   10      0      2      0      3      Y
   //

To parse a TRANSFAC file, use

.. cont-doctest

.. code:: pycon

   >>> with open("transfac.dat") as handle:
   ...     record = motifs.parse(handle, "TRANSFAC")
   ...

If any discrepancies between the file contents and the TRANSFAC file
format are detected, a ``ValueError`` is raised. Note that you may
encounter files that do not follow the TRANSFAC format strictly. For
example, the number of spaces between columns may be different, or a tab
may be used instead of spaces. Use ``strict=False`` to enable parsing
such files without raising a ``ValueError``:

.. code:: pycon

   >>> record = motifs.parse(handle, "TRANSFAC", strict=False)

When parsing a non-compliant file, we recommend to check the record
returned by ``motif.parse`` to ensure that it is consistent with the
file contents.

The overall version number, if available, is stored as
``record.version``:

.. cont-doctest

.. code:: pycon

   >>> record.version
   'EXAMPLE January 15, 2013'

Each motif in ``record`` is in instance of the
``Bio.motifs.transfac.Motif`` class, which inherits both from the
``Bio.motifs.Motif`` class and from a Python dictionary. The dictionary
uses the two-letter keys to store any additional information about the
motif:

.. cont-doctest

.. code:: pycon

   >>> motif = record[0]
   >>> motif.degenerate_consensus  # Using the Bio.motifs.Motif property
   Seq('SRACAGGTGKYG')
   >>> motif["ID"]  # Using motif as a dictionary
   'motif1'

TRANSFAC files are typically much more elaborate than this example,
containing lots of additional information about the motif. Table
:ref:`table:transfaccodes` lists the two-letter field codes that are
commonly found in TRANSFAC files:

.. table:: Fields commonly found in TRANSFAC files
   :name: table:transfaccodes

   ====== ===============================================
   ``AC`` Accession number
   ``AS`` Accession numbers, secondary
   ``BA`` Statistical basis
   ``BF`` Binding factors
   ``BS`` Factor binding sites underlying the matrix
   ``CC`` Comments
   ``CO`` Copyright notice
   ``DE`` Short factor description
   ``DR`` External databases
   ``DT`` Date created/updated
   ``HC`` Subfamilies
   ``HP`` Superfamilies
   ``ID`` Identifier
   ``NA`` Name of the binding factor
   ``OC`` Taxonomic classification
   ``OS`` Species/Taxon
   ``OV`` Older version
   ``PV`` Preferred version
   ``TY`` Type
   ``XX`` Empty line; these are not stored in the Record.
   ====== ===============================================

Each motif also has an attribute ``.references`` containing the
references associated with the motif, using these two-letter keys:

.. container:: center

   .. table:: Fields used to store references in TRANSFAC files

      ====== =================
      ``RN`` Reference number
      ``RA`` Reference authors
      ``RL`` Reference data
      ``RT`` Reference title
      ``RX`` PubMed ID
      ====== =================

Printing the motifs writes them out in their native TRANSFAC format:

.. cont-doctest

.. code:: pycon

   >>> print(record)
   VV  EXAMPLE January 15, 2013
   XX
   //
   ID  motif1
   XX
   P0      A      C      G      T
   01      1      2      2      0      S
   02      2      1      2      0      R
   03      3      0      1      1      A
   04      0      5      0      0      C
   05      5      0      0      0      A
   06      0      0      4      1      G
   07      0      1      4      0      G
   08      0      0      0      5      T
   09      0      0      5      0      G
   10      0      1      2      2      K
   11      0      2      0      3      Y
   12      1      0      3      1      G
   XX
   //
   ID  motif2
   XX
   P0      A      C      G      T
   01      2      1      2      0      R
   02      1      2      2      0      S
   03      0      5      0      0      C
   04      3      0      1      1      A
   05      0      0      4      1      G
   06      5      0      0      0      A
   07      0      1      4      0      G
   08      0      0      5      0      G
   09      0      0      0      5      T
   10      0      2      0      3      Y
   XX
   //
   <BLANKLINE>

You can export the motifs in the TRANSFAC format by capturing this
output in a string and saving it in a file:

.. code:: pycon

   >>> text = str(record)
   >>> with open("mytransfacfile.dat", "w") as out_handle:
   ...     out_handle.write(text)
   ...

The generic ``pfm-four-columns`` format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If none of the tool-specific motif formats work for your PFM file
and your PFM file has the values organized in a 4 columns format,
you can try the generic ``pfm-four-columns`` motif parser:

.. code:: text

    # CIS-BP
    Pos A   C   G   T
    1   0.00961538461538462 0.00961538461538462 0.00961538461538462 0.971153846153846
    2   0.00961538461538462 0.00961538461538462 0.00961538461538462 0.971153846153846
    3   0.971153846153846   0.00961538461538462 0.00961538461538462 0.00961538461538462
    4   0.00961538461538462 0.00961538461538462 0.00961538461538462 0.971153846153846
    5   0.00961538461538462 0.971153846153846   0.00961538461538462 0.00961538461538462
    6   0.971153846153846   0.00961538461538462 0.00961538461538462 0.00961538461538462
    7   0.00961538461538462 0.971153846153846   0.00961538461538462 0.00961538461538462
    8   0.00961538461538462 0.00961538461538462 0.00961538461538462 0.971153846153846

    # C2H2-ZFs
    Gene    ENSG00000197372
    Pos A   C   G   T
    1   0.341303    0.132427    0.117054    0.409215
    2   0.283785    0.077066    0.364552    0.274597
    3   0.491055    0.078208    0.310520    0.120217
    4   0.492621    0.076117    0.131007    0.300256
    5   0.250645    0.361464    0.176504    0.211387
    6   0.276694    0.498070    0.197793    0.027444
    7   0.056317    0.014631    0.926202    0.002850
    8   0.004470    0.007769    0.983797    0.003964
    9   0.936213    0.058787    0.002387    0.002613
    10  0.004352    0.004030    0.002418    0.989200
    11  0.013277    0.008165    0.001991    0.976567
    12  0.968132    0.002263    0.002868    0.026737
    13  0.397623    0.052017    0.350783    0.199577
    14  0.000000    0.000000    1.000000    0.000000
    15  1.000000    0.000000    0.000000    0.000000
    16  0.000000    0.000000    1.000000    0.000000
    17  0.000000    0.000000    1.000000    0.000000
    18  1.000000    0.000000    0.000000    0.000000
    19  0.000000    1.000000    0.000000    0.000000
    20  1.000000    0.000000    0.000000    0.000000

    # C2H2-ZFs
    Gene    FBgn0000210
    Motif   M1734_0.90
    Pos A   C   G   T
    1   0.25    0.0833333   0.0833333   0.583333
    2   0.75    0.166667    0.0833333   0
    3   0.833333    0   0   0.166667
    4   1   0   0   0
    5   0   0.833333    0.0833333   0.0833333
    6   0.333333    0   0   0.666667
    7   0.833333    0   0   0.166667
    8   0.5 0   0.333333    0.166667
    9   0.5 0.0833333   0.166667    0.25
    10  0.333333    0.25    0.166667    0.25
    11  0.166667    0.25    0.416667    0.166667

    # FlyFactorSurvey (Cluster Buster)
    >AbdA_Cell_FBgn0000014
    1   3   0   14
    0   0   0   18
    16  0   0   2
    18  0   0   0
    1   0   0   17
    0   0   6   12
    15  1   2   0

    # HOMER
    >ATGACTCATC AP-1(bZIP)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer    6.049537    -1.782996e+03   0   9805.3,5781.0,3085.1,2715.0,0.00e+00
    0.419   0.275   0.277   0.028
    0.001   0.001   0.001   0.997
    0.010   0.002   0.965   0.023
    0.984   0.003   0.001   0.012
    0.062   0.579   0.305   0.054
    0.026   0.001   0.001   0.972
    0.043   0.943   0.001   0.012
    0.980   0.005   0.001   0.014
    0.050   0.172   0.307   0.471
    0.149   0.444   0.211   0.195

    # HOCOMOCO
    > AHR_si
    40.51343240527031  18.259112547756697  56.41253757072521  38.77363485291994
    10.877470982533044  11.870876719950774  34.66312982331297  96.54723985087516
    21.7165707818416  43.883079837598544  20.706746561638717  67.6523201955933
    2.5465132509466635  1.3171620263517245  145.8637051322628  4.231336967110781
    0.0  150.35847450464382  1.4927836298652875  2.1074592421627525
    3.441039751299748  0.7902972158110341  149.37613720253387  0.3512432070271259
    0.0  3.441039751299748  0.7024864140542533  149.81519121131782
    0.0  0.0  153.95871737667187  0.0
    43.07922333291745  66.87558226865211  16.159862546986584  27.844049228115868

    # Neph
    UW.Motif.0001   atgactca
    0.772949    0.089579    0.098612    0.038860
    0.026652    0.004653    0.025056    0.943639
    0.017663    0.023344    0.918728    0.040264
    0.919596    0.025414    0.029759    0.025231
    0.060312    0.772259    0.104968    0.062462
    0.037406    0.020643    0.006667    0.935284
    0.047316    0.899024    0.026928    0.026732
    0.948639    0.019497    0.005737    0.026128

    # Tiffin
    T   A   G   C
    30  0   28  40
    0   0   0   99
    0   55  14  29
    0   99  0   0
    20  78  0   0
    0   52  7   39
    19  46  11  22
    0   60  38  0
    0   33  0   66
    73  0   25  0
    99  0   0   0

The motifs are read as follows:

.. code:: pycon

   >>> with open("fourcolumns.pfm") as fh:
   ...     for m in motifs.parse(fh, "pfm-four-columns"):
   ...         print(m.name, m.counts, sep="\n")
   ...

           0      1      2      3      4      5      6      7
   G:   0.01   0.01   0.01   0.01   0.01   0.01   0.01   0.01
   A:   0.01   0.01   0.97   0.01   0.01   0.97   0.01   0.01
   T:   0.97   0.97   0.01   0.97   0.01   0.01   0.01   0.97
   C:   0.01   0.01   0.01   0.01   0.97   0.01   0.97   0.01

   ENSG00000197372
           0      1      2      3      4      5      6      7      8      9     10     11     12     13     14     15     16     17     18     19
   G:   0.12   0.36   0.31   0.13   0.18   0.20   0.93   0.98   0.00   0.00   0.00   0.00   0.35   1.00   0.00   1.00   1.00   0.00   0.00   0.00
   A:   0.34   0.28   0.49   0.49   0.25   0.28   0.06   0.00   0.94   0.00   0.01   0.97   0.40   0.00   1.00   0.00   0.00   1.00   0.00   1.00
   T:   0.41   0.27   0.12   0.30   0.21   0.03   0.00   0.00   0.00   0.99   0.98   0.03   0.20   0.00   0.00   0.00   0.00   0.00   0.00   0.00
   C:   0.13   0.08   0.08   0.08   0.36   0.50   0.01   0.01   0.06   0.00   0.01   0.00   0.05   0.00   0.00   0.00   0.00   0.00   1.00   0.00

   M1734_0.90
           0      1      2      3      4      5      6      7      8      9     10
   G:   0.08   0.08   0.00   0.00   0.08   0.00   0.00   0.33   0.17   0.17   0.42
   A:   0.25   0.75   0.83   1.00   0.00   0.33   0.83   0.50   0.50   0.33   0.17
   T:   0.58   0.00   0.17   0.00   0.08   0.67   0.17   0.17   0.25   0.25   0.17
   C:   0.08   0.17   0.00   0.00   0.83   0.00   0.00   0.00   0.08   0.25   0.25

   AbdA_Cell_FBgn0000014
           0      1      2      3      4      5      6
   G:   0.00   0.00   0.00   0.00   0.00   6.00   2.00
   A:   1.00   0.00  16.00  18.00   1.00   0.00  15.00
   T:  14.00  18.00   2.00   0.00  17.00  12.00   0.00
   C:   3.00   0.00   0.00   0.00   0.00   0.00   1.00

   ATGACTCATC AP-1(bZIP)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer    6.049537    -1.782996e+03   0   9805.3,5781.0,3085.1,2715.0,0.00e+00
           0      1      2      3      4      5      6      7      8      9
   G:   0.28   0.00   0.96   0.00   0.30   0.00   0.00   0.00   0.31   0.21
   A:   0.42   0.00   0.01   0.98   0.06   0.03   0.04   0.98   0.05   0.15
   T:   0.03   1.00   0.02   0.01   0.05   0.97   0.01   0.01   0.47   0.20
   C:   0.28   0.00   0.00   0.00   0.58   0.00   0.94   0.01   0.17   0.44

   AHR_si
           0      1      2      3      4      5      6      7      8
   G:  56.41  34.66  20.71 145.86   1.49 149.38   0.70 153.96  16.16
   A:  40.51  10.88  21.72   2.55   0.00   3.44   0.00   0.00  43.08
   T:  38.77  96.55  67.65   4.23   2.11   0.35 149.82   0.00  27.84
   C:  18.26  11.87  43.88   1.32 150.36   0.79   3.44   0.00  66.88


           0      1      2      3      4      5      6      7
   G:   0.10   0.03   0.92   0.03   0.10   0.01   0.03   0.01
   A:   0.77   0.03   0.02   0.92   0.06   0.04   0.05   0.95
   T:   0.04   0.94   0.04   0.03   0.06   0.94   0.03   0.03
   C:   0.09   0.00   0.02   0.03   0.77   0.02   0.90   0.02


           0      1      2      3      4      5      6      7      8      9     10
   G:  28.00   0.00  14.00   0.00   0.00   7.00  11.00  38.00   0.00  25.00   0.00
   A:   0.00   0.00  55.00  99.00  78.00  52.00  46.00  60.00  33.00   0.00   0.00
   T:  30.00   0.00   0.00   0.00  20.00   0.00  19.00   0.00   0.00  73.00  99.00
   C:  40.00  99.00  29.00   0.00   0.00  39.00  22.00   0.00  66.00   0.00   0.00

The generic ``pfm-four-rows`` format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If none of the tool-specific motif formats work for your PFM file
and your PFM file has the values organized in a 4 rows format,
you can try the generic ``pfm-four-rows`` motif parser:

.. code:: text

    # hDPI
    A   0   5   6   5   1   0
    C   1   1   0   0   0   4
    G   5   0   0   0   3   0
    T   0   0   0   1   2   2

    # YeTFaSCo
    A   0.5 0.0 0.0 0.25    0.25    0.25    0.25    0.25    0.25    0.25    0.25    0.25    0.5 0.0 0.0833333334583333
    T   0.0 0.0 0.0 0.25    0.25    0.25    0.25    0.25    0.25    0.25    0.25    0.25    0.0 0.0 0.0833333334583333
    G   0.0 1.0 0.0 0.25    0.25    0.25    0.25    0.25    0.25    0.25    0.25    0.25    0.0 1.0 0.249999999875
    C   0.5 0.0 1.0 0.25    0.25    0.25    0.25    0.25    0.25    0.25    0.25    0.25    0.5 0.0 0.583333333208333

    # FlyFactorSurvey ZFP finger
    A |     92    106    231    135      0      1    780     28      0    700    739     94     60    127    130
    C |    138     82    129     81    774      1      3      1      0      6     17     49    193    122    148
    G |    270    398     54    164      7    659      1    750    755     65      1     41    202    234    205
    T |    290    204    375    411      9    127      6     11     36     20     31    605    335    307    308

    # ScerTF pcm
    A | 9 1 1 97 1 94
    T | 80 1 97 1 1 2
    C | 9 97 1 1 1 2
    G | 2 1 1 1 97 2

    # ScerTF pfm
    A | 0.090 0.010 0.010 0.970 0.010 0.940
    C | 0.090 0.970 0.010 0.010 0.010 0.020
    G | 0.020 0.010 0.010 0.010 0.970 0.020
    T | 0.800 0.010 0.970 0.010 0.010 0.020

    # iDMMPMM
    > abd-A
    0.218451749734889 0.0230646871686108 0.656680805938494 0.898197242841994 0.040694591728526 0.132953340402969 0.74907211028632 0.628313891834571
    0.0896076352067868 0.317338282078473 0.321580063626723 0.0461293743372216 0.0502386002120891 0.040694591728526 0.0284994697773065 0.0339342523860021
    0.455991516436904 0.0691940615058324 0.0108695652173913 0.0217391304347826 0.0284994697773065 0.0284994697773065 0.016304347826087 0.160127253446448
    0.235949098621421 0.590402969247084 0.0108695652173913 0.0339342523860021 0.880567338282079 0.797852598091198 0.206124072110286 0.17762460233298

    # JASPAR
    >MA0001.1 AGL3
    A  [ 0  3 79 40 66 48 65 11 65  0 ]
    C  [94 75  4  3  1  2  5  2  3  3 ]
    G  [ 1  0  3  4  1  0  5  3 28 88 ]
    T  [ 2 19 11 50 29 47 22 81  1  6 ]

    # JASPAR
    >MA0001.1 AGL3
    0  3 79 40 66 48 65 11 65  0
    94 75  4  3  1  2  5  2  3  3
    1  0  3  4  1  0  5  3 28 88
    2 19 11 50 29 47 22 81  1  6

    # Cys2His2 Zinc Finger Proteins PWM Predictor
    base       1       2       3       4       5       6       7       8       9
       a   0.116   0.974   0.444   0.116   0.974   0.444   0.667   0.939   0.068  # noqa: RST301
       c   0.718   0.006   0.214   0.718   0.006   0.214   0.143   0.006   0.107  # noqa: RST301
       g   0.016   0.020   0.028   0.016   0.020   0.028   0.047   0.045   0.216  # noqa: RST301
       t   0.150   0.001   0.314   0.150   0.001   0.314   0.143   0.009   0.609  # noqa: RST301

    Ent=   1.210   0.202   1.665   1.210   0.202   1.665   1.399   0.396   1.521

The motifs are read as follows:

.. code:: pycon

   >>> with open("fourrows.pfm") as fh:
   ...     for m in motifs.parse(fh, "pfm-four-rows"):
   ...         print(m.name, m.counts, sep="\n")
   ...

           0      1      2      3      4      5
   G:   5.00   0.00   0.00   0.00   3.00   0.00
   A:   0.00   5.00   6.00   5.00   1.00   0.00
   T:   0.00   0.00   0.00   1.00   2.00   2.00
   C:   1.00   1.00   0.00   0.00   0.00   4.00


           0      1      2      3      4      5      6      7      8      9     10     11     12     13     14
   G:   0.00   1.00   0.00   0.25   0.25   0.25   0.25   0.25   0.25   0.25   0.25   0.25   0.00   1.00   0.25
   A:   0.50   0.00   0.00   0.25   0.25   0.25   0.25   0.25   0.25   0.25   0.25   0.25   0.50   0.00   0.08
   T:   0.00   0.00   0.00   0.25   0.25   0.25   0.25   0.25   0.25   0.25   0.25   0.25   0.00   0.00   0.08
   C:   0.50   0.00   1.00   0.25   0.25   0.25   0.25   0.25   0.25   0.25   0.25   0.25   0.50   0.00   0.58


           0      1      2      3      4      5      6      7      8      9     10     11     12     13     14
   G: 270.00 398.00  54.00 164.00   7.00 659.00   1.00 750.00 755.00  65.00   1.00  41.00 202.00 234.00 205.00
   A:  92.00 106.00 231.00 135.00   0.00   1.00 780.00  28.00   0.00 700.00 739.00  94.00  60.00 127.00 130.00
   T: 290.00 204.00 375.00 411.00   9.00 127.00   6.00  11.00  36.00  20.00  31.00 605.00 335.00 307.00 308.00
   C: 138.00  82.00 129.00  81.00 774.00   1.00   3.00   1.00   0.00   6.00  17.00  49.00 193.00 122.00 148.00


           0      1      2      3      4      5
   G:   2.00   1.00   1.00   1.00  97.00   2.00
   A:   9.00   1.00   1.00  97.00   1.00  94.00
   T:  80.00   1.00  97.00   1.00   1.00   2.00
   C:   9.00  97.00   1.00   1.00   1.00   2.00


           0      1      2      3      4      5
   G:   0.02   0.01   0.01   0.01   0.97   0.02
   A:   0.09   0.01   0.01   0.97   0.01   0.94
   T:   0.80   0.01   0.97   0.01   0.01   0.02
   C:   0.09   0.97   0.01   0.01   0.01   0.02

   abd-A
           0      1      2      3      4      5      6      7
   G:   0.46   0.07   0.01   0.02   0.03   0.03   0.02   0.16
   A:   0.22   0.02   0.66   0.90   0.04   0.13   0.75   0.63
   T:   0.24   0.59   0.01   0.03   0.88   0.80   0.21   0.18
   C:   0.09   0.32   0.32   0.05   0.05   0.04   0.03   0.03

   MA0001.1 AGL3
           0      1      2      3      4      5      6      7      8      9
   G:   1.00   0.00   3.00   4.00   1.00   0.00   5.00   3.00  28.00  88.00
   A:   0.00   3.00  79.00  40.00  66.00  48.00  65.00  11.00  65.00   0.00
   T:   2.00  19.00  11.00  50.00  29.00  47.00  22.00  81.00   1.00   6.00
   C:  94.00  75.00   4.00   3.00   1.00   2.00   5.00   2.00   3.00   3.00

   MA0001.1 AGL3
           0      1      2      3      4      5      6      7      8      9
   G:   1.00   0.00   3.00   4.00   1.00   0.00   5.00   3.00  28.00  88.00
   A:   0.00   3.00  79.00  40.00  66.00  48.00  65.00  11.00  65.00   0.00
   T:   2.00  19.00  11.00  50.00  29.00  47.00  22.00  81.00   1.00   6.00
   C:  94.00  75.00   4.00   3.00   1.00   2.00   5.00   2.00   3.00   3.00


           0      1      2      3      4      5      6      7      8
   G:   0.02   0.02   0.03   0.02   0.02   0.03   0.05   0.04   0.22
   A:   0.12   0.97   0.44   0.12   0.97   0.44   0.67   0.94   0.07
   T:   0.15   0.00   0.31   0.15   0.00   0.31   0.14   0.01   0.61
   C:   0.72   0.01   0.21   0.72   0.01   0.21   0.14   0.01   0.11

Writing motifs
--------------

Speaking of exporting, let’s look at export functions in general. We can
use the ``format`` built-in function to write the motif in the simple
JASPAR ``pfm`` format:

.. code:: pycon

   >>> print(format(arnt, "pfm"))
     4.00  19.00   0.00   0.00   0.00   0.00
    16.00   0.00  20.00   0.00   0.00   0.00
     0.00   1.00   0.00  20.00   0.00  20.00
     0.00   0.00   0.00   0.00  20.00   0.00

Similarly, we can use ``format`` to write the motif in the JASPAR
``jaspar`` format:

.. code:: pycon

   >>> print(format(arnt, "jaspar"))
   >MA0004.1  Arnt
   A [  4.00  19.00   0.00   0.00   0.00   0.00]
   C [ 16.00   0.00  20.00   0.00   0.00   0.00]
   G [  0.00   1.00   0.00  20.00   0.00  20.00]
   T [  0.00   0.00   0.00   0.00  20.00   0.00]

To write the motif in a TRANSFAC-like matrix format, use

.. cont-doctest

.. code:: pycon

   >>> print(format(m, "transfac"))
   P0      A      C      G      T
   01      3      0      0      4      W
   02      7      0      0      0      A
   03      0      5      0      2      C
   04      2      2      3      0      V
   05      1      6      0      0      C
   XX
   //
   <BLANKLINE>

To write out multiple motifs, you can use ``motifs.write``. This
function can be used regardless of whether the motifs originated from a
TRANSFAC file. For example,

.. cont-doctest

.. code:: pycon

   >>> two_motifs = [arnt, srf]
   >>> print(motifs.write(two_motifs, "transfac"))
   P0      A      C      G      T
   01      4     16      0      0      C
   02     19      0      1      0      A
   03      0     20      0      0      C
   04      0      0     20      0      G
   05      0      0      0     20      T
   06      0      0     20      0      G
   XX
   //
   P0      A      C      G      T
   01      2      1     39      4      G
   02      9     33      2      2      C
   03      0     45      1      0      C
   04      1     45      0      0      C
   05     32      1      0     13      A
   06      3      1      0     42      T
   07     46      0      0      0      A
   08      1      0      0     45      T
   09     43      0      0      3      A
   10     15      1      0     30      W
   11      2      0     44      0      G
   12      2      1     43      0      G
   XX
   //
   <BLANKLINE>

Or, to write multiple motifs in the ``jaspar`` format:

.. code:: pycon

   >>> two_motifs = [arnt, mef2a]
   >>> print(motifs.write(two_motifs, "jaspar"))
   >MA0004.1  Arnt
   A [  4.00  19.00   0.00   0.00   0.00   0.00]
   C [ 16.00   0.00  20.00   0.00   0.00   0.00]
   G [  0.00   1.00   0.00  20.00   0.00  20.00]
   T [  0.00   0.00   0.00   0.00  20.00   0.00]
   >MA0052.1  MEF2A
   A [  1.00   0.00  57.00   2.00   9.00   6.00  37.00   2.00  56.00   6.00]
   C [ 50.00   0.00   1.00   1.00   0.00   0.00   0.00   0.00   0.00   0.00]
   G [  0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   2.00  50.00]
   T [  7.00  58.00   0.00  55.00  49.00  52.00  21.00  56.00   0.00   2.00]

Position-Weight Matrices
------------------------

The ``.counts`` attribute of a Motif object shows how often each
nucleotide appeared at each position along the alignment. We can
normalize this matrix by dividing by the number of instances in the
alignment, resulting in the probability of each nucleotide at each
position along the alignment. We refer to these probabilities as the
position-weight matrix. However, beware that in the literature this term
may also be used to refer to the position-specific scoring matrix, which
we discuss below.

Usually, pseudocounts are added to each position before normalizing.
This avoids overfitting of the position-weight matrix to the limited
number of motif instances in the alignment, and can also prevent
probabilities from becoming zero. To add a fixed pseudocount to all
nucleotides at all positions, specify a number for the ``pseudocounts``
argument:

.. cont-doctest

.. code:: pycon

   >>> pwm = m.counts.normalize(pseudocounts=0.5)
   >>> print(pwm)
           0      1      2      3      4
   A:   0.39   0.83   0.06   0.28   0.17
   C:   0.06   0.06   0.61   0.28   0.72
   G:   0.06   0.06   0.06   0.39   0.06
   T:   0.50   0.06   0.28   0.06   0.06
   <BLANKLINE>

Alternatively, ``pseudocounts`` can be a dictionary specifying the
pseudocounts for each nucleotide. For example, as the GC content of the
human genome is about 40%, you may want to choose the pseudocounts
accordingly:

.. cont-doctest

.. code:: pycon

   >>> pwm = m.counts.normalize(pseudocounts={"A": 0.6, "C": 0.4, "G": 0.4, "T": 0.6})
   >>> print(pwm)
           0      1      2      3      4
   A:   0.40   0.84   0.07   0.29   0.18
   C:   0.04   0.04   0.60   0.27   0.71
   G:   0.04   0.04   0.04   0.38   0.04
   T:   0.51   0.07   0.29   0.07   0.07
   <BLANKLINE>

The position-weight matrix has its own methods to calculate the
consensus, anticonsensus, and degenerate consensus sequences:

.. cont-doctest

.. code:: pycon

   >>> pwm.consensus
   Seq('TACGC')
   >>> pwm.anticonsensus
   Seq('CCGTG')
   >>> pwm.degenerate_consensus
   Seq('WACNC')

Note that due to the pseudocounts, the degenerate consensus sequence
calculated from the position-weight matrix is slightly different from
the degenerate consensus sequence calculated from the instances in the
motif:

.. cont-doctest

.. code:: pycon

   >>> m.degenerate_consensus
   Seq('WACVC')

The reverse complement of the position-weight matrix can be calculated
directly from the ``pwm``:

.. cont-doctest

.. code:: pycon

   >>> rpwm = pwm.reverse_complement()
   >>> print(rpwm)
           0      1      2      3      4
   A:   0.07   0.07   0.29   0.07   0.51
   C:   0.04   0.38   0.04   0.04   0.04
   G:   0.71   0.27   0.60   0.04   0.04
   T:   0.18   0.29   0.07   0.84   0.40
   <BLANKLINE>

Position-Specific Scoring Matrices
----------------------------------

Using the background distribution and PWM with pseudo-counts added, it’s
easy to compute the log-odds ratios, telling us what are the log odds of
a particular symbol to be coming from a motif against the background. We
can use the ``.log_odds()`` method on the position-weight matrix:

.. cont-doctest

.. code:: pycon

   >>> pssm = pwm.log_odds()
   >>> print(pssm)
           0      1      2      3      4
   A:   0.68   1.76  -1.91   0.21  -0.49
   C:  -2.49  -2.49   1.26   0.09   1.51
   G:  -2.49  -2.49  -2.49   0.60  -2.49
   T:   1.03  -1.91   0.21  -1.91  -1.91
   <BLANKLINE>

Here we can see positive values for symbols more frequent in the motif
than in the background and negative for symbols more frequent in the
background. :math:`0.0` means that it’s equally likely to see a symbol
in the background and in the motif.

This assumes that A, C, G, and T are equally likely in the background.
To calculate the position-specific scoring matrix against a background
with unequal probabilities for A, C, G, T, use the ``background``
argument. For example, against a background with a 40% GC content, use

.. cont-doctest

.. code:: pycon

   >>> background = {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}
   >>> pssm = pwm.log_odds(background)
   >>> print(pssm)
           0      1      2      3      4
   A:   0.42   1.49  -2.17  -0.05  -0.75
   C:  -2.17  -2.17   1.58   0.42   1.83
   G:  -2.17  -2.17  -2.17   0.92  -2.17
   T:   0.77  -2.17  -0.05  -2.17  -2.17
   <BLANKLINE>

The maximum and minimum score obtainable from the PSSM are stored in the
``.max`` and ``.min`` properties:

.. cont-doctest

.. code:: pycon

   >>> print("%4.2f" % pssm.max)
   6.59
   >>> print("%4.2f" % pssm.min)
   -10.85

The mean and standard deviation of the PSSM scores with respect to a
specific background are calculated by the ``.mean`` and ``.std``
methods.

.. cont-doctest

.. code:: pycon

   >>> mean = pssm.mean(background)
   >>> std = pssm.std(background)
   >>> print("mean = %0.2f, standard deviation = %0.2f" % (mean, std))
   mean = 3.21, standard deviation = 2.59

A uniform background is used if ``background`` is not specified. The
mean is equal to the Kullback-Leibler divergence or relative entropy
described in Section :ref:`subsec:relative_entropy`.

The ``.reverse_complement``, ``.consensus``, ``.anticonsensus``, and
``.degenerate_consensus`` methods can be applied directly to PSSM
objects.

.. _`sec:search`:

Searching for instances
-----------------------

The most frequent use for a motif is to find its instances in some
sequence. For the sake of this section, we will use an artificial
sequence like this:

.. cont-doctest

.. code:: pycon

   >>> test_seq = Seq("TACACTGCATTACAACCCAAGCATTA")
   >>> len(test_seq)
   26

Searching for exact matches
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The simplest way to find instances, is to look for exact matches of the
true instances of the motif:

.. cont-doctest

.. code:: pycon

   >>> for pos, seq in test_seq.search(m.alignment):
   ...     print("%i %s" % (pos, seq))
   ...
   0 TACAC
   10 TACAA
   13 AACCC

We can do the same with the reverse complement (to find instances on the
complementary strand):

.. cont-doctest

.. code:: pycon

   >>> for pos, seq in test_seq.search(r.alignment):
   ...     print("%i %s" % (pos, seq))
   ...
   6 GCATT
   20 GCATT

Searching for matches using the PSSM score
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It’s just as easy to look for positions, giving rise to high log-odds
scores against our motif:

.. cont-doctest

.. code:: pycon

   >>> for position, score in pssm.search(test_seq, threshold=3.0):
   ...     print("Position %d: score = %5.3f" % (position, score))
   ...
   Position 0: score = 5.622
   Position -20: score = 4.601
   Position 10: score = 3.037
   Position 13: score = 5.738
   Position -6: score = 4.601

The negative positions refer to instances of the motif found on the
reverse strand of the test sequence, and follow the Python convention on
negative indices. Therefore, the instance of the motif at ``pos`` is
located at ``test_seq[pos:pos+len(m)]`` both for positive and for
negative values of ``pos``.

You may notice the threshold parameter, here set arbitrarily to
:math:`3.0`. This is in :math:`log_2`, so we are now looking only for
words, which are eight times more likely to occur under the motif model
than in the background. The default threshold is :math:`0.0`, which
selects everything that looks more like the motif than the background.

You can also calculate the scores at all positions along the sequence:

.. code:: pycon

   >>> pssm.calculate(test_seq)
   array([  5.62230396,  -5.6796999 ,  -3.43177247,   0.93827754,
           -6.84962511,  -2.04066086, -10.84962463,  -3.65614533,
           -0.03370807,  -3.91102552,   3.03734159,  -2.14918518,
           -0.6016975 ,   5.7381525 ,  -0.50977498,  -3.56422281,
           -8.73414803,  -0.09919716,  -0.6016975 ,  -2.39429784,
          -10.84962463,  -3.65614533], dtype=float32)

In general, this is the fastest way to calculate PSSM scores. The scores
returned by ``pssm.calculate`` are for the forward strand only. To
obtain the scores on the reverse strand, you can take the reverse
complement of the PSSM:

.. code:: pycon

   >>> rpssm = pssm.reverse_complement()
   >>> rpssm.calculate(test_seq)
   array([ -9.43458748,  -3.06172252,  -7.18665981,  -7.76216221,
           -2.04066086,  -4.26466274,   4.60124254,  -4.2480607 ,
           -8.73414803,  -2.26503372,  -6.49598789,  -5.64668512,
           -8.73414803, -10.84962463,  -4.82356262,  -4.82356262,
           -5.64668512,  -8.73414803,  -4.15613794,  -5.6796999 ,
            4.60124254,  -4.2480607 ], dtype=float32)

Selecting a score threshold
~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you want to use a less arbitrary way of selecting thresholds, you can
explore the distribution of PSSM scores. Since the space for a score
distribution grows exponentially with motif length, we are using an
approximation with a given precision to keep computation cost
manageable:

.. cont-doctest

.. code:: pycon

   >>> distribution = pssm.distribution(background=background, precision=10**4)

The ``distribution`` object can be used to determine a number of
different thresholds. We can specify the requested false-positive rate
(probability of “finding” a motif instance in background generated
sequence):

.. cont-doctest

.. code:: pycon

   >>> threshold = distribution.threshold_fpr(0.01)
   >>> print("%5.3f" % threshold)
   4.009

or the false-negative rate (probability of “not finding” an instance
generated from the motif):

.. cont-doctest

.. code:: pycon

   >>> threshold = distribution.threshold_fnr(0.1)
   >>> print("%5.3f" % threshold)
   -0.510

or a threshold (approximately) satisfying some relation between the
false-positive rate and the false-negative rate
(:math:`\frac{\textrm{fnr}}{\textrm{fpr}}\simeq t`):

.. cont-doctest

.. code:: pycon

   >>> threshold = distribution.threshold_balanced(1000)
   >>> print("%5.3f" % threshold)
   6.241

or a threshold satisfying (roughly) the equality between the
:math:`-log` of the false-positive rate and the information content (as
used in patser software by Hertz and Stormo):

.. cont-doctest

.. code:: pycon

   >>> threshold = distribution.threshold_patser()
   >>> print("%5.3f" % threshold)
   0.346

For example, in case of our motif, you can get the threshold giving you
exactly the same results (for this sequence) as searching for instances
with balanced threshold with rate of :math:`1000`.

.. cont-doctest

.. code:: pycon

   >>> threshold = distribution.threshold_fpr(0.01)
   >>> print("%5.3f" % threshold)
   4.009
   >>> for position, score in pssm.search(test_seq, threshold=threshold):
   ...     print("Position %d: score = %5.3f" % (position, score))
   ...
   Position 0: score = 5.622
   Position -20: score = 4.601
   Position 13: score = 5.738
   Position -6: score = 4.601

Each motif object has an associated Position-Specific Scoring Matrix
--------------------------------------------------------------------

To facilitate searching for potential TFBSs using PSSMs, both the
position-weight matrix and the position-specific scoring matrix are
associated with each motif. Using the Arnt motif as an example:

.. cont-doctest

.. code:: pycon

   >>> from Bio import motifs
   >>> with open("Arnt.sites") as handle:
   ...     motif = motifs.read(handle, "sites")
   ...
   >>> print(motif.counts)
           0      1      2      3      4      5
   A:   4.00  19.00   0.00   0.00   0.00   0.00
   C:  16.00   0.00  20.00   0.00   0.00   0.00
   G:   0.00   1.00   0.00  20.00   0.00  20.00
   T:   0.00   0.00   0.00   0.00  20.00   0.00
   <BLANKLINE>
   >>> print(motif.pwm)
           0      1      2      3      4      5
   A:   0.20   0.95   0.00   0.00   0.00   0.00
   C:   0.80   0.00   1.00   0.00   0.00   0.00
   G:   0.00   0.05   0.00   1.00   0.00   1.00
   T:   0.00   0.00   0.00   0.00   1.00   0.00
   <BLANKLINE>
   >>> print(motif.pssm)
           0      1      2      3      4      5
   A:  -0.32   1.93   -inf   -inf   -inf   -inf
   C:   1.68   -inf   2.00   -inf   -inf   -inf
   G:   -inf  -2.32   -inf   2.00   -inf   2.00
   T:   -inf   -inf   -inf   -inf   2.00   -inf
   <BLANKLINE>

The negative infinities appear here because the corresponding entry in
the frequency matrix is 0, and we are using zero pseudocounts by
default:

.. cont-doctest

.. code:: pycon

   >>> for letter in "ACGT":
   ...     print("%s: %4.2f" % (letter, motif.pseudocounts[letter]))
   ...
   A: 0.00
   C: 0.00
   G: 0.00
   T: 0.00

If you change the ``.pseudocounts`` attribute, the position-frequency
matrix and the position-specific scoring matrix are recalculated
automatically:

.. cont-doctest

.. code:: pycon

   >>> motif.pseudocounts = 3.0
   >>> for letter in "ACGT":
   ...     print("%s: %4.2f" % (letter, motif.pseudocounts[letter]))
   ...
   A: 3.00
   C: 3.00
   G: 3.00
   T: 3.00
   >>> print(motif.pwm)
           0      1      2      3      4      5
   A:   0.22   0.69   0.09   0.09   0.09   0.09
   C:   0.59   0.09   0.72   0.09   0.09   0.09
   G:   0.09   0.12   0.09   0.72   0.09   0.72
   T:   0.09   0.09   0.09   0.09   0.72   0.09
   <BLANKLINE>

.. cont-doctest

.. code:: pycon

   >>> print(motif.pssm)
           0      1      2      3      4      5
   A:  -0.19   1.46  -1.42  -1.42  -1.42  -1.42
   C:   1.25  -1.42   1.52  -1.42  -1.42  -1.42
   G:  -1.42  -1.00  -1.42   1.52  -1.42   1.52
   T:  -1.42  -1.42  -1.42  -1.42   1.52  -1.42
   <BLANKLINE>

You can also set the ``.pseudocounts`` to a dictionary over the four
nucleotides if you want to use different pseudocounts for them. Setting
``motif.pseudocounts`` to ``None`` resets it to its default value of
zero.

The position-specific scoring matrix depends on the background
distribution, which is uniform by default:

.. cont-doctest

.. code:: pycon

   >>> for letter in "ACGT":
   ...     print("%s: %4.2f" % (letter, motif.background[letter]))
   ...
   A: 0.25
   C: 0.25
   G: 0.25
   T: 0.25

Again, if you modify the background distribution, the position-specific
scoring matrix is recalculated:

.. cont-doctest

.. code:: pycon

   >>> motif.background = {"A": 0.2, "C": 0.3, "G": 0.3, "T": 0.2}
   >>> print(motif.pssm)
           0      1      2      3      4      5
   A:   0.13   1.78  -1.09  -1.09  -1.09  -1.09
   C:   0.98  -1.68   1.26  -1.68  -1.68  -1.68
   G:  -1.68  -1.26  -1.68   1.26  -1.68   1.26
   T:  -1.09  -1.09  -1.09  -1.09   1.85  -1.09
   <BLANKLINE>

Setting ``motif.background`` to ``None`` resets it to a uniform
distribution:

.. cont-doctest

.. code:: pycon

   >>> motif.background = None
   >>> for letter in "ACGT":
   ...     print("%s: %4.2f" % (letter, motif.background[letter]))
   ...
   A: 0.25
   C: 0.25
   G: 0.25
   T: 0.25

If you set ``motif.background`` equal to a single value, it will be
interpreted as the GC content:

.. cont-doctest

.. code:: pycon

   >>> motif.background = 0.8
   >>> for letter in "ACGT":
   ...     print("%s: %4.2f" % (letter, motif.background[letter]))
   ...
   A: 0.10
   C: 0.40
   G: 0.40
   T: 0.10

Note that you can now calculate the mean of the PSSM scores over the
background against which it was computed:

.. cont-doctest

.. code:: pycon

   >>> print("%f" % motif.pssm.mean(motif.background))
   4.703928

as well as its standard deviation:

.. cont-doctest

.. code:: pycon

   >>> print("%f" % motif.pssm.std(motif.background))
   3.290900

and its distribution:

.. cont-doctest

.. code:: pycon

   >>> distribution = motif.pssm.distribution(background=motif.background)
   >>> threshold = distribution.threshold_fpr(0.01)
   >>> print("%f" % threshold)
   3.854375

Note that the position-weight matrix and the position-specific scoring
matrix are recalculated each time you call ``motif.pwm`` or
``motif.pssm``, respectively. If speed is an issue and you want to use
the PWM or PSSM repeatedly, you can save them as a variable, as in

.. code:: pycon

   >>> pssm = motif.pssm

.. _`sec:comp`:

Comparing motifs
----------------

Once we have more than one motif, we might want to compare them.

Before we start comparing motifs, I should point out that motif
boundaries are usually quite arbitrary. This means we often need to
compare motifs of different lengths, so comparison needs to involve some
kind of alignment. This means we have to take into account two things:

-  alignment of motifs

-  some function to compare aligned motifs

To align the motifs, we use ungapped alignment of PSSMs and substitute
zeros for any missing columns at the beginning and end of the matrices.
This means that effectively we are using the background distribution for
columns missing from the PSSM. The distance function then returns the
minimal distance between motifs, as well as the corresponding offset in
their alignment.

To give an example, let us first load another motif, which is similar to
our test motif ``m``:

.. cont-doctest

.. code:: pycon

   >>> with open("REB1.pfm") as handle:
   ...     m_reb1 = motifs.read(handle, "pfm")
   ...
   >>> m_reb1.consensus
   Seq('GTTACCCGG')
   >>> print(m_reb1.counts)
           0      1      2      3      4      5      6      7      8
   A:  30.00   0.00   0.00 100.00   0.00   0.00   0.00   0.00  15.00
   C:  10.00   0.00   0.00   0.00 100.00 100.00 100.00   0.00  15.00
   G:  50.00   0.00   0.00   0.00   0.00   0.00   0.00  60.00  55.00
   T:  10.00 100.00 100.00   0.00   0.00   0.00   0.00  40.00  15.00
   <BLANKLINE>

To make the motifs comparable, we choose the same values for the
pseudocounts and the background distribution as our motif ``m``:

.. cont-doctest

.. code:: pycon

   >>> m_reb1.pseudocounts = {"A": 0.6, "C": 0.4, "G": 0.4, "T": 0.6}
   >>> m_reb1.background = {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}
   >>> pssm_reb1 = m_reb1.pssm
   >>> print(pssm_reb1)
           0      1      2      3      4      5      6      7      8
   A:   0.00  -5.67  -5.67   1.72  -5.67  -5.67  -5.67  -5.67  -0.97
   C:  -0.97  -5.67  -5.67  -5.67   2.30   2.30   2.30  -5.67  -0.41
   G:   1.30  -5.67  -5.67  -5.67  -5.67  -5.67  -5.67   1.57   1.44
   T:  -1.53   1.72   1.72  -5.67  -5.67  -5.67  -5.67   0.41  -0.97
   <BLANKLINE>

We’ll compare these motifs using the Pearson correlation. Since we want
it to resemble a distance measure, we actually take :math:`1-r`, where
:math:`r` is the Pearson correlation coefficient (PCC):

.. cont-doctest

.. code:: pycon

   >>> distance, offset = pssm.dist_pearson(pssm_reb1)
   >>> print("distance = %5.3g" % distance)
   distance = 0.239
   >>> print(offset)
   -2

This means that the best PCC between motif ``m`` and ``m_reb1`` is
obtained with the following alignment:

.. code:: text

   m:      bbTACGCbb
   m_reb1: GTTACCCGG

where ``b`` stands for background distribution. The PCC itself is
roughly :math:`1-0.239=0.761`.

.. _`sec:find`:

*De novo* motif finding
-----------------------

Currently, Biopython has only limited support for *de novo* motif
finding. Namely, we support running ``xxmotif`` and also parsing of
MEME. Since the number of motif finding tools is growing rapidly,
contributions of new parsers are welcome.

.. _`sec:meme`:

MEME
~~~~

Let’s assume, you have run MEME on sequences of your choice with your
favorite parameters and saved the output in the file ``meme.out``. You
can retrieve the motifs reported by MEME by running the following piece
of code:

.. doctest ../Tests/motifs

.. code:: pycon

   >>> from Bio import motifs
   >>> with open("meme.psp_test.classic.zoops.xml") as handle:
   ...     motifsM = motifs.parse(handle, "meme")
   ...

.. code:: pycon

   >>> motifsM
   [<Bio.motifs.meme.Motif object at 0xc356b0>]

Besides the most wanted list of motifs, the result object contains more
useful information, accessible through properties with self-explanatory
names:

-  ``.alphabet``

-  ``.datafile``

-  ``.sequences``

-  ``.version``

-  ``.command``

The motifs returned by the MEME Parser can be treated exactly like
regular Motif objects (with instances), they also provide some extra
functionality, by adding additional information about the instances.

.. cont-doctest

.. code:: pycon

   >>> motifsM[0].consensus
   Seq('GCTTATGTAA')
   >>> motifsM[0].alignment.sequences[0].sequence_name
   'iYFL005W'
   >>> motifsM[0].alignment.sequences[0].sequence_id
   'sequence_15'
   >>> motifsM[0].alignment.sequences[0].start
   480
   >>> motifsM[0].alignment.sequences[0].strand
   '+'

.. code:: pycon

   >>> motifsM[0].alignment.sequences[0].pvalue
   1.97e-06

.. _`sec:links`:

Useful links
------------

-  `Sequence motif <https://en.wikipedia.org/wiki/Sequence_motif>`__ in
   wikipedia

-  `PWM <https://en.wikipedia.org/wiki/Position_weight_matrix>`__ in
   wikipedia

-  `Consensus
   sequence <https://en.wikipedia.org/wiki/Consensus_sequence>`__ in
   wikipedia

-  `Comparison of different motif finding
   programs <http://bio.cs.washington.edu/assessment/>`__
