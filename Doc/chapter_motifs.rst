.. _chapter:motifs:

Sequence motif analysis using Bio.motifs
========================================

This chapter gives an overview of the functionality of the
``Bio.motifs`` package included in Biopython. It is intended for people
who are involved in the analysis of sequence motifs, so I’ll assume that
you are familiar with basic notions of motif analysis. In case something
is unclear, please look at Section :ref:`sec:links` for some
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

.. _sec:object:

Motif objects
-------------

Since we are interested in motif analysis, we need to take a look at
``Motif`` objects in the first place. For that we need to import the
Bio.motifs library:

.. code:: pycon

   >>> from Bio import motifs

and we can start creating our first motif objects. We can either create
a ``Motif`` object from a list of instances of the motif, or we can
obtain a ``Motif`` object by parsing a file from a motif database or
motif finding software.

Creating a motif from instances
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Suppose we have these instances of a DNA motif:

.. code:: pycon

   >>> from Bio.Seq import Seq
   >>> instances = [Seq("TACAA"),
   ...              Seq("TACGC"),
   ...              Seq("TACAC"),
   ...              Seq("TACCC"),
   ...              Seq("AACCC"),
   ...              Seq("AATGC"),
   ...              Seq("AATGC"),
   ...             ]

then we can create a Motif object as follows:

.. code:: pycon

   >>> m = motifs.create(instances)

The instances are saved in an attribute ``m.instances``, which is
essentially a Python list with some added functionality, as described
below. Printing out the Motif object shows the instances from which it
was constructed:

.. code:: pycon

   >>> print(m)
   TACAA
   TACGC
   TACAC
   TACCC
   AACCC
   AATGC
   AATGC
   <BLANKLINE>

The length of the motif is defined as the sequence length, which should
be the same for all instances:

.. code:: pycon

   >>> len(m)
   5

The Motif object has an attribute ``.counts`` containing the counts of
each nucleotide at each position. Printing this counts matrix shows it
in an easily readable format:

.. code:: pycon

   >>> print(m.counts)
           0      1      2      3      4
   A:   3.00   7.00   0.00   2.00   1.00
   C:   0.00   0.00   5.00   2.00   6.00
   G:   0.00   0.00   0.00   3.00   0.00
   T:   4.00   0.00   2.00   0.00   0.00
   <BLANKLINE>

You can access these counts as a dictionary:

.. code:: pycon

   >>> m.counts["A"]
   [3, 7, 0, 2, 1]

but you can also think of it as a 2D array with the nucleotide as the
first dimension and the position as the second dimension:

.. code:: pycon

   >>> m.counts["T", 0]
   4
   >>> m.counts["T", 2]
   2
   >>> m.counts["T", 3]
   0

You can also directly access columns of the counts matrix

.. code:: pycon

   >>> m.counts[:, 3]
   {'A': 2, 'C': 2, 'T': 0, 'G': 3}

Instead of the nucleotide itself, you can also use the index of the
nucleotide in the alphabet of the motif:

.. code:: pycon

   >>> m.alphabet
   'ACGT'
   >>> m.counts["A",:]
   (3, 7, 0, 2, 1)
   >>> m.counts[0,:]
   (3, 7, 0, 2, 1)

The motif has an associated consensus sequence, defined as the sequence
of letters along the positions of the motif for which the largest value
in the corresponding columns of the ``.counts`` matrix is obtained:

.. code:: pycon

   >>> m.consensus
   Seq('TACGC')

as well as an anticonsensus sequence, corresponding to the smallest
values in the columns of the ``.counts`` matrix:

.. code:: pycon

   >>> m.anticonsensus
   Seq('CCATG')

Note that there is some ambiguity in the definition of the consensus and
anticonsensus sequence if in some columns multiple nucleotides have the
maximum or minimum count.

You can also ask for a degenerate consensus sequence, in which ambiguous
nucleotides are used for positions where there are multiple nucleotides
with high counts:

.. code:: pycon

   >>> m.degenerate_consensus
   Seq('WACVC')

Here, W and R follow the IUPAC nucleotide ambiguity codes: W is either A
or T, and V is A, C, or G :raw-latex:`\cite{cornish1985}`. The
degenerate consensus sequence is constructed following the rules
specified by Cavener :raw-latex:`\cite{cavener1987}`.

We can also get the reverse complement of a motif:

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
   <BLANKLINE>

The reverse complement and the degenerate consensus sequence are only
defined for DNA motifs.

Creating a sequence logo
~~~~~~~~~~~~~~~~~~~~~~~~

If we have internet access, we can create a
`weblogo <https://weblogo.berkeley.edu>`__:

.. code:: pycon

   >>> m.weblogo("mymotif.png")

We should get our logo saved as a PNG in the specified file.

.. _sec:io:

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

-  ``tf_class`` - the structual class of this TF, e.g. ’Zipper-Type’

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

.. code:: pycon

   >>> from Bio import motifs
   >>> with open("Arnt.sites") as handle:
   ...     arnt = motifs.read(handle, "sites")
   ...

The instances from which this motif was created is stored in the
``.instances`` property:

.. code:: pycon

   >>> print(arnt.instances[:3])
   [Seq('CACGTG'), Seq('CACGTG'), Seq('CACGTG')]
   >>> for instance in arnt.instances:
   ...     print(instance)
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

.. code:: pycon

   >>> print(srf.instances)
   None

We can now ask for the consensus sequence of these two motifs:

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
   >>> for m in motifs.parse(fh, "jaspar"))
   ...     print(m)
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
   >>> JASPAR_DB_HOST = <hostname>
   >>> JASPAR_DB_NAME = <db_name>
   >>> JASPAR_DB_USER = <user>
   >>> JASPAR_DB_PASS = <passord>
   >>>
   >>> jdb = JASPAR5(
   ...     host=JASPAR_DB_HOST,
   ...     name=JASPAR_DB_NAME,
   ...     user=JASPAR_DB_USER,
   ...     password=JASPAR_DB_PASS
   ... )

Now we can fetch a single motif by its unique JASPAR ID with the
``fetch_motif_by_id`` method. Note that a JASPAR ID conists of a base ID
and a version number seperated by a decimal point, e.g. ’MA0004.1’. The
``fetch_motif_by_id`` method allows you to use either the fully
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
   ...     collection = "CORE",
   ...     tax_group = ["vertebrates", "insects"],
   ...     tf_class = "Winged Helix-Turn-Helix",
   ...     tf_family = ["Forkhead", "Ets"],
   ...     min_ic = 12
   ... )
   >>> for motif in motifs:
   ...     pass # do something with the motif

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
     behaviour is to compute the IC without first applying pseudocounts,
     even though by default the PSSMs are computed using pseudocounts as
     described above.

-  | **Searching for instances:**
   | Searching for instances with the Perl ``TFBS`` motifs was usually
     performed using a relative score threshold, i.e. a score in the
     range 0 to 1. In order to compute the absolute PSSM score
     corresponding to a relative score one can use the equation:

   .. code:: pycon

      >>> abs_score =  (pssm.max - pssm.min) * rel_score + pssm.min

   To convert the absolute score of an instance back to a relative
   score, one can use the equation:

   .. code:: pycon

      >>> rel_score = (abs_score - pssm.min) / (pssm.max - pssm.min)

   For example, using the Arnt motif before, let’s search a sequence
   with a relative score threshold of 0.8.

   .. code:: pycon

      >>> test_seq=Seq("TAAGCGTGCACGCGCAACACGTGCATTA")
      >>> arnt.pseudocounts = motifs.jaspar.calculate_pseudocounts(arnt)
      >>> pssm = arnt.pssm
      >>> max_score = pssm.max
      >>> min_score = pssm.min
      >>> abs_score_threshold = (max_score - min_score) * 0.8 + min_score
      >>> for position, score in pssm.search(test_seq,
                                             threshold=abs_score_threshold):
      ...     rel_score = (score - min_score) / (max_score - min_score)
      ...     print("Position %d: score = %5.3f, rel. score = %5.3f" % (
                  position, score, rel_score))
      ...
      Position 2: score = 5.362, rel. score = 0.801
      Position 8: score = 6.112, rel. score = 0.831
      Position -20: score = 7.103, rel. score = 0.870
      Position 17: score = 10.351, rel. score = 1.000
      Position -11: score = 10.351, rel. score = 1.000

MEME
~~~~

MEME :raw-latex:`\cite{bailey1994}` is a tool for discovering motifs in
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

.. code:: pycon

   >>> with open("meme.dna.oops.txt") as handle:
   ...     record = motifs.parse(handle, "meme")
   ...

The ``motifs.parse`` command reads the complete file directly, so you
can close the file after calling ``motifs.parse``. The header
information is stored in attributes:

.. code:: pycon

   >>> record.version
   '3.0'
   >>> record.datafile
   'INO_up800.s'
   >>> record.command
   'meme -mod oops -dna -revcomp -nmotifs 2 -bfile yeast.nc.6.freq INO_up800.s'
   >>> record.alphabet
   'ACGT'
   >>> record.sequences
   ['CHO1', 'CHO2', 'FAS1', 'FAS2', 'ACC1', 'INO1', 'OPI3']

The record is an object of the ``Bio.motifs.meme.Record`` class. The
class inherits from list, and you can think of ``record`` as a list of
Motif objects:

.. code:: pycon

   >>> len(record)
   2
   >>> motif = record[0]
   >>> print(motif.consensus)
   TTCACATGCCGC
   >>> print(motif.degenerate_consensus)
   TTCACATGSCNC

In addition to these generic motif attributes, each motif also stores
its specific information as calculated by MEME. For example,

.. code:: pycon

   >>> motif.num_occurrences
   7
   >>> motif.length
   12
   >>> evalue = motif.evalue
   >>> print("%3.1g" % evalue)
   0.2
   >>> motif.name
   'Motif 1'

In addition to using an index into the record, as we did above, you can
also find it by its name:

.. code:: pycon

   >>> motif = record["Motif 1"]

Each motif has an attribute ``.instances`` with the sequence instances
in which the motif was found, providing some information on each
instance:

.. code:: pycon

   >>> len(motif.instances)
   7
   >>> motif.instances[0]
   Instance('TTCACATGCCGC', 'ACGT')
   >>> motif.instances[0].motif_name
   'Motif 1'
   >>> motif.instances[0].sequence_name
   'INO1'
   >>> motif.instances[0].start
   620
   >>> motif.instances[0].strand
   '-'
   >>> motif.instances[0].length
   12
   >>> pvalue = motif.instances[0].pvalue
   >>> print("%5.3g" % pvalue)
   1.85e-08

MAST
^^^^

TRANSFAC
~~~~~~~~

TRANSFAC is a manually curated database of transcription factors,
together with their genomic binding sites and DNA binding profiles
:raw-latex:`\cite{matys2003}`. While the file format used in the
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

.. code:: pycon

   >>> record.version
   'EXAMPLE January 15, 2013'

Each motif in ``record`` is in instance of the
``Bio.motifs.transfac.Motif`` class, which inherits both from the
``Bio.motifs.Motif`` class and from a Python dictionary. The dictionary
uses the two-letter keys to store any additional information about the
motif:

.. code:: pycon

   >>> motif = record[0]
   >>> motif.degenerate_consensus # Using the Bio.motifs.Motif method
   Seq('SRACAGGTGKYG')
   >>> motif["ID"] # Using motif as a dictionary
   'motif1'

TRANSFAC files are typically much more elaborate than this example,
containing lots of additional information about the motif. Table
:ref:`table:transfaccodes` lists the two-letter
field codes that are commonly found in TRANSFAC files:

[table:transfaccodes]

.. table:: Fields commonly found in TRANSFAC files

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

.. table:: Fields used to store references in TRANSFAC files

   ====== =================
   ``RN`` Reference number
   ``RA`` Reference authors
   ``RL`` Reference data
   ``RT`` Reference title
   ``RX`` PubMed ID
   ====== =================

Printing the motifs writes them out in their native TRANSFAC format:

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

Writing motifs
--------------

Speaking of exporting, let’s look at export functions in general. We can
use the ``format`` method to write the motif in the simple JASPAR
``pfm`` format:

.. code:: pycon

   >>> print(arnt.format("pfm"))
     4.00  19.00   0.00   0.00   0.00   0.00
    16.00   0.00  20.00   0.00   0.00   0.00
     0.00   1.00   0.00  20.00   0.00  20.00
     0.00   0.00   0.00   0.00  20.00   0.00

Similarly, we can use ``format`` to write the motif in the JASPAR
``jaspar`` format:

.. code:: pycon

   >>> print(arnt.format("jaspar"))
   >MA0004.1  Arnt
   A [  4.00  19.00   0.00   0.00   0.00   0.00]
   C [ 16.00   0.00  20.00   0.00   0.00   0.00]
   G [  0.00   1.00   0.00  20.00   0.00  20.00]
   T [  0.00   0.00   0.00   0.00  20.00   0.00]

To write the motif in a TRANSFAC-like matrix format, use

.. code:: pycon

   >>> print(m.format("transfac"))
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

.. code:: pycon

   >>> pwm = m.counts.normalize(pseudocounts={"A":0.6, "C": 0.4, "G": 0.4, "T": 0.6})
   >>> print(pwm)
           0      1      2      3      4
   A:   0.40   0.84   0.07   0.29   0.18
   C:   0.04   0.04   0.60   0.27   0.71
   G:   0.04   0.04   0.04   0.38   0.04
   T:   0.51   0.07   0.29   0.07   0.07
   <BLANKLINE>

The position-weight matrix has its own methods to calculate the
consensus, anticonsensus, and degenerate consensus sequences:

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

.. code:: pycon

   >>> m.degenerate_consensus
   Seq('WACVC')

The reverse complement of the position-weight matrix can be calculated
directly from the ``pwm``:

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

.. code:: pycon

   >>> background = {"A":0.3,"C":0.2,"G":0.2,"T":0.3}
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

.. code:: pycon

   >>> print("%4.2f" % pssm.max)
   6.59
   >>> print("%4.2f" % pssm.min)
   -10.85

The mean and standard deviation of the PSSM scores with respect to a
specific background are calculated by the ``.mean`` and ``.std``
methods.

.. code:: pycon

   >>> mean = pssm.mean(background)
   >>> std = pssm.std(background)
   >>> print("mean = %0.2f, standard deviation = %0.2f" % (mean, std))
   mean = 3.21, standard deviation = 2.59

A uniform background is used if ``background`` is not specified. The
mean is particularly important, as its value is equal to the
Kullback-Leibler divergence or relative entropy, and is a measure for
the information content of the motif compared to the background. As in
Biopython the base-2 logarithm is used in the calculation of the
log-odds scores, the information content has units of bits.

The ``.reverse_complement``, ``.consensus``, ``.anticonsensus``, and
``.degenerate_consensus`` methods can be applied directly to PSSM
objects.

.. _sec:search:

Searching for instances
-----------------------

The most frequent use for a motif is to find its instances in some
sequence. For the sake of this section, we will use an artificial
sequence like this:

.. code:: pycon

   >>> test_seq=Seq("TACACTGCATTACAACCCAAGCATTA")
   >>> len(test_seq)
   26

Searching for exact matches
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The simplest way to find instances, is to look for exact matches of the
true instances of the motif:

.. code:: pycon

   >>> for pos, seq in m.instances.search(test_seq):
   ...     print("%i %s" % (pos, seq))
   ...
   0 TACAC
   10 TACAA
   13 AACCC

We can do the same with the reverse complement (to find instances on the
complementary strand):

.. code:: pycon

   >>> for pos, seq in r.instances.search(test_seq):
   ...     print("%i %s" % (pos, seq))
   ...
   6 GCATT
   20 GCATT

Searching for matches using the PSSM score
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It’s just as easy to look for positions, giving rise to high log-odds
scores against our motif:

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

.. code:: pycon

   >>> distribution = pssm.distribution(background=background, precision=10**4)

The ``distribution`` object can be used to determine a number of
different thresholds. We can specify the requested false-positive rate
(probability of “finding” a motif instance in background generated
sequence):

.. code:: pycon

   >>> threshold = distribution.threshold_fpr(0.01)
   >>> print("%5.3f" % threshold)
   4.009

or the false-negative rate (probability of “not finding” an instance
generated from the motif):

.. code:: pycon

   >>> threshold = distribution.threshold_fnr(0.1)
   >>> print("%5.3f" % threshold)
   -0.510

or a threshold (approximately) satisfying some relation between the
false-positive rate and the false-negative rate
(:math:`\frac{\textrm{fnr}}{\textrm{fpr}}\simeq t`):

.. code:: pycon

   >>> threshold = distribution.threshold_balanced(1000)
   >>> print("%5.3f" % threshold)
   6.241

or a threshold satisfying (roughly) the equality between the
:math:`-log` of the false-positive rate and the information content (as
used in patser software by Hertz and Stormo):

.. code:: pycon

   >>> threshold = distribution.threshold_patser()
   >>> print("%5.3f" % threshold)
   0.346

For example, in case of our motif, you can get the threshold giving you
exactly the same results (for this sequence) as searching for instances
with balanced threshold with rate of :math:`1000`.

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

.. code:: pycon

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

.. code:: pycon

   >>> motif.pseudocounts = 3.0
   >>> for letter in "ACGT":
   ...     print("%s: %4.2f" % (letter, motif.pseudocounts[letter]))
   ...
   A: 3.00
   C: 3.00
   G: 3.00
   T: 3.00

.. code:: pycon

   >>> print(motif.pwm)
           0      1      2      3      4      5
   A:   0.22   0.69   0.09   0.09   0.09   0.09
   C:   0.59   0.09   0.72   0.09   0.09   0.09
   G:   0.09   0.12   0.09   0.72   0.09   0.72
   T:   0.09   0.09   0.09   0.09   0.72   0.09
   <BLANKLINE>

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

.. code:: pycon

   >>> print("%f" % motif.pssm.mean(motif.background))
   4.703928

as well as its standard deviation:

.. code:: pycon

   >>> print("%f" % motif.pssm.std(motif.background))
   3.290900

and its distribution:

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

.. _sec:comp:

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

.. code:: pycon

   >>> with open("REB1.pfm") as handle:
   ...    m_reb1 = motifs.read(handle, "pfm")
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

.. code:: pycon

   >>> m_reb1.pseudocounts = {"A":0.6, "C": 0.4, "G": 0.4, "T": 0.6}
   >>> m_reb1.background = {"A":0.3,"C":0.2,"G":0.2,"T":0.3}
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

.. _sec:find:

*De novo* motif finding
-----------------------

Currently, Biopython has only limited support for *de novo* motif
finding. Namely, we support running ``xxmotif`` and also parsing of
MEME. Since the number of motif finding tools is growing rapidly,
contributions of new parsers are welcome.

.. _sec:meme:

MEME
~~~~

Let’s assume, you have run MEME on sequences of your choice with your
favorite parameters and saved the output in the file ``meme.out``. You
can retrieve the motifs reported by MEME by running the following piece
of code:

.. code:: pycon

   >>> from Bio import motifs
   >>> with open("meme.out") as handle:
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

-  ``.sequence_names``

-  ``.version``

-  ``.command``

The motifs returned by the MEME Parser can be treated exactly like
regular Motif objects (with instances), they also provide some extra
functionality, by adding additional information about the instances.

.. code:: pycon

   >>> motifsM[0].consensus
   Seq('CTCAATCGTA')
   >>> motifsM[0].instances[0].sequence_name
   'SEQ10;'
   >>> motifsM[0].instances[0].start
   3
   >>> motifsM[0].instances[0].strand
   '+'

.. code:: pycon

   >>> motifsM[0].instances[0].pvalue
   8.71e-07

.. _sec:links:

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
