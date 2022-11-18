News for the Biopython Project
==============================

This file contains release notes and general news about the Biopython project.
See also the DEPRECATED file which tracks the removal of obsolete modules or
functions, and online https://biopython.org/wiki/News and
https://www.open-bio.org/category/obf-projects/biopython/

The latest news is at the top of this file.

18 November 2022: Biopython 1.80
================================

This release of Biopython supports Python 3.7, 3.8, 3.9, 3.10, 3.11. It has
also been tested on PyPy3.7 v7.3.5.

Functions ``read``, ``parse``, and ``write`` were added to ``Bio.Align`` to
read and write ``Alignment`` objects.

Because dict retains the item order by default since Python3.6, all instances
of ``collections.OrderedDict`` have been replaced by either standard ``dict``
or where appropriate by ``collections.defaultsdict``.

The ``Bio.motifs.jaspar.db`` now returns ``tf_family`` and ``tf_class`` as a
string array since the JASPAR 2018 release.

The Local Composition Complexity functions from ``Bio.SeqUtils`` now uses
base 4 log instead of 2 as stated in the original reference Konopka (2005),
Sequence Complexity and Composition. https://doi.org/10.1038/npg.els.0005260

Append mode is now supported in ``Bio.bgzf`` (and a bug parsing blocked GZIP
files with an internal empty block fixed).

The experimental warning was dropped from ``Bio.phenotype`` (which was new in
Biopython 1.67).

Sequences now have a ``defined`` attribute that returns a boolean indicating
if the underlying data is defined or not.

The ``Bio.PDB`` module now includes a structural alignment module, using the
combinatorial extension algorithm of Shindyalov and Bourne, commonly known as
CEAlign. The module allows for two structures to be aligned based solely on
their 3D conformation, ie. in a sequence-independent manner. The method is
particularly powerful when the structures shared a very low degree of sequence
similarity. The new module is available in ``Bio.PDB.CEAligner`` with an
interface similar to other 3D superimposition modules.

A new module ``Bio.PDB.qcprot`` implements the QCP superposition algorithm in
pure Python, deprecating the existing C implementation. This leads to a slight
performance improvement and to much better maintainability. The refactored
``qcprot.QCPSuperimposer`` class has small changes to its API, to better mirror
that of ``Bio.PDB.Superimposer``.

The ``Bio.PDB.PDBList`` module now allows downloading biological assemblies,
for one or more entries of the wwPDB.

In the ``Bio.Restriction`` module, each restriction enzyme now includes an `id`
property giving the numerical identifier for the REBASE database identifier
from which the enzyme object was created, and a `uri` property with a canonical
`identifiers.org` link to the database, for use in linked-data representations.

Add new ``gc_fraction`` function in ``SeqUtils`` and marks ``GC`` for future
deprecation.

Support for the old format (dating back to 2004) of the GN line in SwissProt
files was dropped in ``Bio.SwissProt``.

Additionally, a number of small bugs and typos have been fixed with additions
to the test suite.

Many thanks to the Biopython developers and community for making this release
possible, especially the following contributors:

- Andrius Merkys
- Arup Ghosh (first contribution)
- Alessio Quercia
- Aziz Khan
- Alex Morehead
- Caio Fontes
- Chenghao Zhu
- Christian Brueffer
- Damien Goutte-Gattat
- Erik Weßels (first contribution)
- Erik  Whiting
- Fabian Egli
- Fredric Johansson
- Hongbo Zhu
- Hussein Faara (first contribution)
- Manuel Lera Ramirez
- Jacob Beal (first contribution)
- Jarrod Millman
- João Rodrigues
- Josha Inglis
- Markus Piotrowski
- Michiel de Hoon
- Neil P. (first contribution)
- Peter Cock
- Robert Sawicki (first contribution)
- Sebastian Bassi
- Sean Aubin
- Sean Workman (first contribution)
- Soroush Saffari (first contribution)
- Tim Burke
- Valentin Vareškić (first contribution)

3 June 2021: Biopython 1.79
===========================

This is intended to be our final release supporting Python 3.6. It also
supports Python 3.7, 3.8 and 3.9, and has also been tested on PyPy3.6.1 v7.1.1.

The ``Seq`` and ``MutableSeq`` classes in ``Bio.Seq`` now store their sequence
contents as ``bytes`` and ``bytearray`` objects, respectively. Previously, for
``Seq`` objects a string object was used, and a Unicode array object for
``MutableSeq`` objects. This was maintained during the transition from Python2
to Python3. However, a Python2 string object corresponds to a ``bytes`` object
in Python3, storing the string as a series of 256-bit characters. While
non-ASCII characters could be stored in Python2 strings, they were not treated
as such. For example:

In Python2::

    >>> s = "Генетика"
    >>> type(s)
    <class 'str'>
    >>> len(s)
    16

In Python3::

    >>> s = "Генетика"
    >>> type(s)
    <class 'str'>
    >>> len(s)
    8

In Python3, storing the sequence contents as ``bytes`` and ``bytearray``
objects has the further advantage that both support the buffer protocol.

Taking advantage of the similarity between ``bytes`` and ``bytearray``, the
``Seq`` and ``MutableSeq`` classes now inherit from an abstract base class
``_SeqAbstractBaseClass`` in ``Bio.Seq`` that implements most of the ``Seq``
and ``MutableSeq`` methods, ensuring their consistency with each other. For
methods that modify the sequence contents, an optional ``inplace`` argument to
specify if a new sequence object should be returned with the new sequence
contents (if ``inplace`` is ``False``, the default) or if the sequence object
itself should be modified (if ``inplace`` is ``True``). For ``Seq`` objects,
which are immutable, using ``inplace=True`` raises an exception. For
``inplace=False``, the default, ``Seq`` objects and ``MutableSeq`` behave
consistently.

As before, ``Seq`` and ``MutableSeq`` objects can be initialized using a string
object, which will be converted to a ``bytes`` or ``bytearray`` object assuming
an ASCII encoding. Alternatively, a ``bytes`` or ``bytearray`` object can be
used, or an instance of any class inheriting from the new
``SequenceDataAbstractBaseClass`` abstract base class in ``Bio.Seq``. This
requires that the class implements the ``__len__`` and ``__getitem`` methods
that return the sequence length and sequence contents on demand. Initializing a
``Seq`` instance using an instance of a class inheriting from
``SequenceDataAbstractBaseClass`` allows the ``Seq`` object to be lazy, meaning
that its sequence is provided on demand only, without requiring to initialize
the full sequence. This feature is now used in ``BioSQL``, providing on-demand
sequence loading from an SQL database, as well as in a new parser for twoBit
(.2bit) sequence data added to ``Bio.SeqIO``. This is a lazy parser that allows
fast access to genome-size DNA sequence files by not having to read the full
genome sequence. The new ``_UndefinedSequenceData`` class in ``Bio.Seq``  also
inherits from ``SequenceDataAbstractBaseClass`` to represent sequences of known
length but unknown sequence contents. This provides an alternative to
``UnknownSeq``, which is now deprecated as its definition was ambiguous. For
example, in these examples the ``UnknownSeq`` is interpreted as a sequence with
a well-defined sequence contents::

    >>> s = UnknownSeq(3, character="A")
    >>> s.translate()
    UnknownSeq(1, character='K')
    >>> s + "A"
    Seq("AAAA")

A sequence object with an undefined sequence contents can now be created by
using ``None`` when creating the ``Seq`` object, together with the sequence
length. Trying to access its sequence contents raises an
``UndefinedSequenceError``::

    >>> s = Seq(None, length=6)
    >>> s
    Seq(None, length=6)
    >>> len(s)
    6
    >>> "A" in s
    Traceback (most recent call last):
    ...
    Bio.Seq.UndefinedSequenceError: Sequence content is undefined
    >>> print(s)
    Traceback (most recent call last):
    ....
    Bio.Seq.UndefinedSequenceError: Sequence content is undefined

Element assignment in Bio.PDB.Atom now returns "X" when the element cannot be
unambiguously guessed from the atom name, in accordance with PDB structures.

Bio.PDB entities now have a ``center_of_mass()`` method that calculates either
centers of gravity or geometry.

New method ``disordered_remove()`` implemented in Bio.PDB DisorderedAtom and
DisorderedResidue to remove children.

New module Bio.PDB.SASA implements the Shrake-Rupley algorithm to calculate
atomic solvent accessible areas without third-party tools.

Expected ``TypeError`` behaviour has been restored to the ``Seq`` object's
string like methods (fixing a regression in Biopython 1.78).

The KEGG ``KGML_Pathway`` KGML output was fixed to produce output that complies
with KGML v0.7.2.

Parsing motifs in ``pfm-four-rows`` format can now handle motifs with values
in scientific notation.

Parsing motifs in ``minimal`` MEME format will use ``nsites`` when making
the count matrix from the frequency matrix, instead of multiply the frequency
matrix by 1000000.

Bio.UniProt.GOA now parses Gene Product Information (GPI) files version 1.2,
files can be downloaded from the EBI ftp site:
ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/

Many thanks to the Biopython developers and community for making this release
possible, especially the following contributors:

- Damien Goutte-Gattat
- Gert Hulselmans
- João Rodrigues
- Markus Piotrowski
- Pascal Schläpfer (first contribution)
- Leighton Pritchard
- Sergio Valqui
- Suyash Gupta
- Vini Salazar (first contribution)


4 September 2020: Biopython 1.78
================================

This release of Biopython supports Python 3.6, 3.7 and 3.8. It has also been
tested on PyPy3.6.1 v7.1.1.

The main change is that ``Bio.Alphabet`` is no longer used. In some cases you
will now have to specify expected letters, molecule type (DNA, RNA, protein),
or gap character explicitly. Please consult the updated Tutorial and API
documentation for guidance. This simplification has sped up many ``Seq``
object methods. See https://biopython.org/wiki/Alphabet for more information.

``Bio.SeqIO.parse()`` is faster with "fastq" format due to small improvements
in the ``Bio.SeqIO.QualityIO`` module.

The ``SeqFeature`` object's ``.extract()`` method can now be used for
trans-spliced locations via an optional dictionary of references.

As in recent releases, more of our code is now explicitly available under
either our original "Biopython License Agreement", or the very similar but
more commonly used "3-Clause BSD License".  See the ``LICENSE.rst`` file for
more details.

Additionally, a number of small bugs and typos have been fixed with additions
to the test suite. There has been further work to follow the Python PEP8,
PEP257 and best practice standard coding style, and all of the tests have
been reformatted with the ``black`` tool to match the main code base.

Many thanks to the Biopython developers and community for making this release
possible, especially the following contributors:

- Adam Sjøgren (first contribution)
- Carlos Pena
- Chris Daley
- Chris Rands
- Christian Brueffer
- Damien Goutte-Gattat
- João Rodrigues
- João Vitor F Cavalcante (first contribution)
- Marie Crane
- Markus Piotrowski
- Michiel de Hoon
- Peter Cock
- Sergio Valqui
- Yogesh Kulkarni (first contribution)
- Zheng Ruan

25 May 2020: Biopython 1.77
===========================

This release of Biopython supports Python 3.6, 3.7 and 3.8 It has also been
tested on PyPy3.6.1 v7.1.1-beta0.

**We have dropped support for Python 2 now.**

``pairwise2`` now allows the input of parameters with keywords and returns the
alignments as a list of ``namedtuples``.

The codon tables have been updated to NCBI genetic code table version 4.5,
which adds Cephalodiscidae mitochondrial as table 33.

Updated ``Bio.Restriction`` to the January 2020 release of REBASE.

A major contribution by Rob Miller to ``Bio.PDB`` provides new methods to
handle protein structure transformations using dihedral angles (internal
coordinates). The new framework supports lossless interconversion between
internal and cartesian coordinates, which, among other uses, simplifies the
analysis and manipulation of coordinates of proteins structures.

As in recent releases, more of our code is now explicitly available under
either our original "Biopython License Agreement", or the very similar but
more commonly used "3-Clause BSD License".  See the ``LICENSE.rst`` file for
more details.

Additionally, a number of small bugs and typos have been fixed with further
additions to the test suite. There has been further work to follow the Python
PEP8, PEP257 and best practice standard coding style, and all the main code
base has been reformatted with the ``black`` tool.

Many thanks to the Biopython developers and community for making this release
possible, especially the following contributors:

- Alexander Decurnou (first contribution)
- Andrei Istrate (first contribution)
- Andrey Raspopov
- Artemi Bendandi (first contribution)
- Austin Varela (first contribution)
- Chris Daley
- Chris Rands
- Deepak Khatri
- Hielke Walinga (first contribution)
- Kai Blin
- Karthikeyan Singaravelan (first contribution)
- Konstantinos Zisis (first contribution)
- Markus Piotrowski
- Michiel de Hoon
- Peter Cock
- Rob Miller
- Sergio Valqui
- Steve Bond
- Sujan Dulal (first contribution)
- Tianyi Shi (first contribution)

20 December 2019: Biopython 1.76
================================

This release of Biopython supports Python 2.7, 3.5, 3.6, 3.7 and 3.8. It has
also been tested on PyPy2.7.13 v7.1.1 and PyPy3.6.1 v7.1.1-beta0.

We intend this to be our final release supporting Python 2.7 and 3.5.

As in recent releases, more of our code is now explicitly available under
either our original "Biopython License Agreement", or the very similar but
more commonly used "3-Clause BSD License".  See the ``LICENSE.rst`` file for
more details.


``PDBParser`` and ``PDBIO`` now support PQR format file parsing and input/
output.

In addition to the mainstream ``x86_64`` aka ``AMD64`` CPU architecture, we
now also test every contribution on the ``ARM64``, ``ppc64le``, and ``s390x``
CPUs under Linux thanks to Travis CI. Further post-release testing done by
Debian and other packagers and distributors of Biopython also covers these
CPUs.

``Bio.motifs.PositionSpecificScoringMatrix.search()`` method has been
re-written: it now applies ``.calculate()`` to chunks of the sequence
to maintain a low memory footprint for long sequences.

Additionally, a number of small bugs and typos have been fixed with further
additions to the test suite. There has been further work to follow the Python
PEP8, PEP257 and best practice standard coding style, and more of the code
style has been reformatted with the ``black`` tool.

Many thanks to the Biopython developers and community for making this release
possible, especially the following contributors:

- Chris Daley (first contribution)
- Chris Rands
- Christian Brueffer
- Ilya Flyamer (first contribution)
- Jakub Lipinski (first contribution)
- Michael R. Crusoe (first contribution)
- Michiel de Hoon
- Peter Cock
- Sergio Valqui

6 November 2019: Biopython 1.75
===============================

This release of Biopython supports Python 2.7, 3.5, 3.6, 3.7 and is expected
to work on the soon to be released Python 3.8. It has also been tested on
PyPy2.7.13 v7.1.1 and PyPy3.6.1 v7.1.1-beta0.

Note we intend to drop Python 2.7 support in early 2020.

The restriction enzyme list in ``Bio.Restriction`` has been updated to the
August 2019 release of REBASE.

``Bio.SeqIO`` now supports reading and writing files in the native format of
Christian Marck's DNA Strider program ("xdna" format, also used by Serial
Cloner), as well as reading files in the native formats of GSL Biotech's
SnapGene ("snapgene") and Textco Biosoftware's Gene Construction Kit ("gck").

``Bio.AlignIO`` now supports GCG MSF multiple sequence alignments as the "msf"
format (work funded by the National Marrow Donor Program).

The main ``Seq`` object now has string-like ``.index()`` and ``.rindex()``
methods, matching the existing ``.find()`` and ``.rfind()`` implementations.
The ``MutableSeq`` object retains its more list-like ``.index()`` behaviour.

The ``MMTFIO`` class has been added that allows writing of MMTF file format
files from a Biopython structure object. ``MMTFIO`` has a similar interface to
``PDBIO`` and ``MMCIFIO``, including the use of a ``Select`` class to write
out a specified selection. This final addition to read/write support for
PDB/mmCIF/MMTF in Biopython allows conversion between all three file formats.

Values from mmCIF files are now read in as a list even when they consist of a
single value. This change improves consistency and reduces the likelihood of
making an error, but will require user code to be updated accordingly.

`Bio.motifs.meme` has been updated to parse XML output files from MEME over
the plain-text output file. The goal of this change is to parse a more
structured data source with minimal loss of functionality upon future MEME
releases.

``Bio.PDB`` has been updated to support parsing REMARK 99 header entries from
PDB-style Astral files.

A new keyword parameter ``full_sequences`` was added to ``Bio.pairwise2``'s
pretty print method ``format_alignment`` to restore the output of local
alignments to the 'old' format (showing the whole sequences including the
un-aligned parts instead of only showing the aligned parts).

A new function ``charge_at_pH(pH)`` has been added to ``ProtParam`` and
``IsoelectricPoint`` in ``Bio.SeqUtils``.

The ``PairwiseAligner`` in ``Bio.Align`` was extended to allow generalized
pairwise alignments, i.e. alignments of any Python object, for example
three-letter amino acid sequences, three-nucleotide codons, and arrays of
integers.

A new module ``substitution_matrices`` was added to ``Bio.Align``, which
includes an ``Array`` class that can be used as a substitution matrix. As
the ``Array`` class is a subclass of a numpy array, mathematical operations
can be applied to it directly, and C code that makes use of substitution
matrices can directly access the numerical values stored in the substitution
matrices. This module is intended as a replacement of ``Bio.SubsMat``,
which is currently unmaintained.

As in recent releases, more of our code is now explicitly available under
either our original "Biopython License Agreement", or the very similar but
more commonly used "3-Clause BSD License".  See the ``LICENSE.rst`` file for
more details.

Additionally, a number of small bugs and typos have been fixed with further
additions to the test suite, and there has been further work to follow the
Python PEP8, PEP257 and best practice standard coding style. We have also
started to use the ``black`` Python code formatting tool.

Many thanks to the Biopython developers and community for making this release
possible, especially the following contributors:

- Chris MacRaild
- Chris Rands
- Damien Goutte-Gattat (first contribution)
- Devang Thakkar
- Harry Jubb
- Joe Greener
- Kiran Mukhyala (first contribution)
- Konstantin Vdovkin
- Mark Amery
- Markus Piotrowski
- Michiel de Hoon
- Mike Moritz (first contribution)
- Mustafa Anil Tuncel
- Nick Negretti
- Osvaldo Zagordi (first contribution)
- Peter Cock
- Peter Kerpedjiev
- Sergio Valqui
- Spencer Bliven
- Victor Lin


16 July 2019: Biopython 1.74
============================

This release of Biopython supports Python 2.7, 3.4, 3.5, 3.6 and 3.7. However,
it will be the last release to support Python 3.4 which is now at end-of-life.
It has also been tested on PyPy2.7 v6.0.0 and PyPy3.5 v6.0.0.

As in recent releases, more of our code is now explicitly available under
either our original "Biopython License Agreement", or the very similar but
more commonly used "3-Clause BSD License".  See the ``LICENSE.rst`` file for
more details.

Our core sequence objects (``Seq``, ``UnknownSeq``, and ``MutableSeq``) now
have a string-like ``.join()`` method.

The NCBI now allows longer accessions in the GenBank file LOCUS line, meaning
the fields may not always follow the historical column based positions. We
no longer give a warning when parsing these. We now allow writing such files
(although with a warning as support for reading them is not yet widespread).

Support for the ``mysqlclient`` package, a fork of MySQLdb, has been added.

We now capture the IDcode field from PDB Header records.

``Bio.pairwise2``'s pretty-print output from ``format_alignment`` has been
optimized for local alignments: If they do not consist of the whole sequences,
only the aligned section of the sequences are shown, together with the start
positions of the sequences (in 1-based notation). Alignments of lists will now
also be prettily printed.

``Bio.SearchIO`` now supports parsing the text output of the HHsuite protein
sequence search tool. The format name is ``hhsuite2-text`` and
``hhsuite3-text``, for versions 2 and 3 of HHsuite, respectively.

``Bio.SearchIO`` HSP objects has a new attribute called ``output_index``. This
attribute is meant for capturing the order by which the HSP were output in the
parsed file and is set with a default value of -1 for all HSP objects. It is
also used for sorting the output of ``QueryResult.hsps``.

``Bio.SeqIO.AbiIO`` has been updated to preserve bytes value when parsing. The
goal of this change is make the parser more robust by being able to extract
string-values that are not utf-8-encoded. This affects all tag values, except
for ID and description values, where they need to be extracted as strings
to conform to the ``SeqRecord`` interface. In this case, the parser will
attempt to decode using ``utf-8`` and fall back to the system encoding if that
fails. This change affects Python 3 only.

``Bio.motifs.mast`` has been updated to parse XML output files from MAST over
the plain-text output file. The goal of this change is to parse a more
structured data source with minimal loss of functionality upon future MAST
releases. Class structure remains the same plus an additional attribute
``Record.strand_handling`` required for diagram parsing.

``Bio.Entrez`` now automatically retries HTTP requests on failure. The
maximum number of tries and the sleep between them can be configured by
changing ``Bio.Entrez.max_tries`` and ``Bio.Entrez.sleep_between_tries``.
(The defaults are 3 tries and 15 seconds, respectively.)

The restriction enzyme list in ``Bio.Restriction`` has been updated to the May
2019 release of REBASE.

All tests using the older print-and-compare approach have been replaced by
unittests following Python's standard testing framework.

On the documentation side, all the public modules, classes, methods and
functions now have docstrings (built in help strings). Furthermore, the PDF
version of the *Biopython Tutorial and Cookbook* now uses syntax coloring
for code snippets.

Additionally, a number of small bugs and typos have been fixed with further
additions to the test suite, and there has been further work to follow the
Python PEP8, PEP257 and best practice standard coding style.

Many thanks to the Biopython developers and community for making this release
possible, especially the following contributors:

- Andrey Raspopov (first contribution)
- Antony Lee
- Benjamin Rowell (first contribution)
- Bernhard Thiel
- Brandon Invergo
- Catherine Lesuisse
- Chris Rands
- Deepak Khatri (first contribution)
- Gert Hulselmans
- Jared Andrews
- Jens Thomas (first contribution)
- Konstantin Vdovkin
- Lenna Peterson
- Mark Amery
- Markus Piotrowski
- Micky Yun Chan (first contribution)
- Nick Negretti
- Peter Cock
- Peter Kerpedjiev
- Ralf Stephan
- Rob Miller (first contribution)
- Sergio Valqui
- Victor Lin
- Wibowo 'Bow' Arindrarto
- Zheng Ruan


18 December 2018: Biopython 1.73
================================

This release of Biopython supports Python 2.7, 3.4, 3.5, 3.6 and 3.7.
It has also been tested on PyPy2.7 v6.0.0 and PyPy3.5 v6.0.0.

As in recent releases, more of our code is now explicitly available under
either our original "Biopython License Agreement", or the very similar but
more commonly used "3-Clause BSD License".  See the ``LICENSE.rst`` file for
more details.

The dictionary-like indexing in SeqIO and SearchIO will now explicitly preserve
record order to match a behaviour change in the Python standard dict object.
This means looping over the index will load the records in the on-disk order,
which will be much faster (previously it would be effectively at random, based
on the key hash sorting).

The "grant" matrix in Bio.SubsMat.MatrixInfo has been replaced as our original
values taken from Gerhard Vogt's old webpages at EMBL Heidelberg were
discovered to be in error. The new values have been transformed following
Vogt's approach, taking the global maximum 215 minus the similarity scores
from the original paper Grantham (1974), to give a distance measure.

Additionally, a number of small bugs and typos have been fixed with further
additions to the test suite, and there has been further work to follow the
Python PEP8, PEP257 and best practice standard coding style.

Double-quote characters in GenBank feature qualifier values in ``Bio.SeqIO``
are now escaped as per the NCBI standard. Improperly escaped values trigger a
warning on parsing.

There is a new command line wrapper for the BWA-MEM sequence mapper.

The string-based FASTA parsers in ``Bio.SeqIO.FastaIO`` have been optimised,
which also speeds up parsing FASTA files using ``Bio.SeqIO.parse()``.

Many thanks to the Biopython developers and community for making this release
possible, especially the following contributors:

- Alona Levy-Jurgenson (first contribution)
- Ariel Aptekmann
- Brandon Invergo
- Catherine Lesuisse
- Chris Rands
- Darcy Mason (first contribution)
- Devang Thakkar (first contribution)
- Ivan Antonov (first contribution)
- Jeremy LaBarage (first contribution)
- Juraj Szász (first contribution)
- Kai Blin
- Konstantin Vdovkin (first contribution)
- Manuel Nuno Melo (first contribution)
- Maximilian Greil
- Nick Negretti (first contribution)
- Peter Cock
- Rona Costello (first contribution)
- Spencer Bliven
- Wibowo 'Bow' Arindrarto
- Yi Hsiao (first contribution)


21 June 2018: Biopython 1.72
============================

This release of Biopython supports Python 2.7, 3.4, 3.5 and 3.6.
It has also been tested on PyPy2.7 v6.0.0 and PyPy3.5 v6.0.0.

Internal changes to Bio.SeqIO have sped up the SeqRecord .format method and
SeqIO.write (especially when used in a for loop).

The MAF alignment indexing in Bio.AlignIO.MafIO has been updated to use
inclusive end coordinates to better handle searches at end points. This
will require you to rebuild any existing MAF index files.

In this release more of our code is now explicitly available under either our
original "Biopython License Agreement", or the very similar but more commonly
used "3-Clause BSD License".  See the ``LICENSE.rst`` file for more details.

The Entrez module now supports the NCBI API key. Also you can now set a custom
directory for DTD and XSD files. This allows Entrez to be used in environments
like AWS Lambda, which restricts write access to specific directories.
Improved support for parsing NCBI Entrez XML files that use XSD schemas.

Internal changes to our C code mean that NumPy is no longer required at
compile time - only at run time (and only for those modules which use NumPy).

Seq, UnknownSeq, MutableSeq and derived classes now support integer
multiplication methods, matching native Python string methods.

A translate method has been added to Bio.SeqFeature that will extract a
feature and translate it using the codon_start and transl_table qualifiers
of the feature if they are present.

Bio.SearchIO is no longer considered experimental, and so it does not raise
warnings anymore when imported.

A new pairwise sequence aligner is available in Bio.Align, as an alternative
to the existing pairwise sequence aligner in Bio.pairwise2.

Many thanks to the Biopython developers and community for making this release
possible, especially the following contributors:

- Benjamin Vaisvil (first contribution)
- Blaise Li
- Chad Parmet
- Chris Rands
- Connor T. Skennerton
- Francesco Gastaldello
- Michiel de Hoon
- Pamela Russell (first contribution)
- Peter Cock
- Spencer Bliven
- Stefans Mezulis
- Wibowo 'Bow' Arindrarto


3 April 2018: Biopython 1.71
============================

This release of Biopython supports Python 2.7, 3.4, 3.5 and 3.6.
It has also been tested on PyPy2.7 v5.10.0 and PyPy3.5 v5.10.1.

Python 3 is the primary development platform for Biopython. We will drop
support for Python 2.7 no later than 2020, in line with the end-of-life or
sunset date for Python 2.7 itself.

Encoding issues have been fixed in several parsers when reading data files
with non-ASCII characters, like accented letters in people's names. This would
raise ``UnicodeDecodeError: 'ascii' codec can't decode byte ...`` under some
system locale settings.

Bio.KEGG can now parse Gene files.

The multiple-sequence-alignment object used by Bio.AlignIO etc now supports
a per-column annotation dictionary, useful for richly annotated alignments
in the Stockholm/PFAM format.

The SeqRecord object now has a translate method, following the approach used
for its existing reverse_complement method etc.

The output of function ``format_alignment`` in ``Bio.pairwise2`` for displaying
a pairwise sequence alignment as text now indicates gaps and mis-matches.

Bio.SeqIO now supports reading and writing two-line-per-record FASTA files
under the format name "fasta-2line", useful if you wish to work without
line-wrapped sequences.

Bio.PDB now contains a writer for the mmCIF file format, which has been the
standard PDB archive format since 2014. This allows structural objects to be
written out and facilitates conversion between the PDB and mmCIF file formats.

Bio.Emboss.Applications has been updated to fix a wrong parameter in fuzznuc
wrapper and include a new wrapper for fuzzpro.

The restriction enzyme list in ``Bio.Restriction`` has been updated to the
November 2017 release of REBASE.

New codon tables 27-31 from NCBI (NCBI genetic code table version 4.2)
were added to Bio.Data.CodonTable. Note that tables 27, 28 and 31 contain
no dedicated stop codons; the stop codons in these codes have a context
dependent encoding as either STOP or as amino acid.

IO functions such as ``SeqIO.parse`` now accept any objects which can be passed
to the builtin ``open`` function. Specifically, this allows using
``pathlib.Path`` objects under Python 3.6 and newer, as per `PEP 519
<https://www.python.org/dev/peps/pep-0519/>`_.

Bio.SearchIO can now parse InterProScan XML files.

For Python 3 compatibility, comparison operators for the entities within a
Bio.PDB Structure object were implemented. These allow the comparison of
models, chains, residues, and atoms with the common operators  (==, !=, >, ...)
Comparisons are based on IDs and take the parents of the entity up to the
model level into account. For consistent behaviour of all entities the
operators for atoms were modified to also consider the parent IDs. NOTE: this
represents a change in behaviour in respect to v1.70 for Atom comparisons. In
order to mimic the behaviour of previous versions, comparison will have to be
done for Atom IDs and alternative locations specifically.

In this release more of our code is now explicitly available under either our
original "Biopython License Agreement", or the very similar but more commonly
used "3-Clause BSD License".  See the ``LICENSE.rst`` file for more details.

Additionally, a number of small bugs and typos have been fixed with further
additions to the test suite, and there has been further work to follow the
Python PEP8, PEP257 and best practice standard coding style.

Many thanks to the Biopython developers and community for making this release
possible, especially the following contributors:

- Adhemar Zerlotini
- Ariel Aptekmann
- Chris Rands
- Christian Brueffer
- Connor T. Skennerton
- Erik Cederstrand (first contribution)
- Fei Qi (first contribution)
- Francesco Gastaldello
- James Jeffryes (first contribution)
- Jerven Bolleman (first contribution)
- Joe Greener (first contribution)
- Joerg Schaarschmidt (first contribution)
- João Rodrigues
- Jeroen Van Goey
- Jun Aruga (first contribution)
- Kai Blin
- Kozo Nishida
- Lewis A. Marshall (first contribution)
- Markus Piotrowski
- Michiel de Hoon
- Nicolas Fontrodona (first contribution)
- Peter Cock
- Philip Bergstrom (first contribution)
- rht (first contribution)
- Saket Choudhary
- Shuichiro MAKIGAKI (first contribution)
- Shyam Saladi (first contribution)
- Siong Kong
- Spencer Bliven
- Stefans Mezulis
- Steve Bond
- Yasar L. Ahmed (first contribution)
- Zachary Sailer (first contribution)
- Zaid Ur-Rehman (first contribution)


10 July 2017: Biopython 1.70
============================

This release of Biopython supports Python 2.7, 3.4, 3.5 and 3.6 (we have now
dropped support for Python 3.3). It has also been tested on PyPy v5.7,
PyPy3.5 v5.8 beta, and Jython 2.7 (although support for Jython is deprecated).

Biopython now has a new logo, contributed by Patrick Kunzmann. Drawing on our
original logo and the current Python logo, this shows a yellow and blue snake
forming a double helix.

For installation Biopython now assumes ``setuptools`` is present, and takes
advantage of this to declare we require NumPy at install time (except under
Jython). This should help ensure ``pip install biopython`` works smoothly.

Bio.AlignIO now supports Mauve's eXtended Multi-FastA (XMFA) file format
under the format name "mauve" (contributed by Eric Rasche).

Bio.ExPASy was updated to fix fetching PROSITE and PRODOC records, and return
text-mode handles for use under Python 3.

Two new arguments for reading and writing blast-xml files have been added
to the Bio.SearchIO functions (read/parse and write, respectively). They
are 'use_raw_hit_ids' and 'use_raw_query_ids'. Check out the relevant
SearchIO.BlastIO documentation for a complete description of what these
arguments do.

Bio.motifs was updated to support changes in MEME v4.11.4 output.

The Bio.Seq sequence objects now have a ``.count_overlap()`` method to
supplement the Python string like non-overlap based ``.count()`` method.

The Bio.SeqFeature location objects can now be compared for equality.

Bio.Phylo.draw_graphviz is now deprecated. We recommend using Bio.Phylo.draw
instead, or another library or program if more advanced plotting functionality
is needed.

In Bio.Phylo.TreeConstruction, the DistanceMatrix class (previously
_DistanceMatrix) has a new method 'format_phylip' to write Phylip-compatible
distance matrix files (contributed by Jordan Willis).

Additionally, a number of small bugs have been fixed with further additions
to the test suite, and there has been further work to follow the Python PEP8,
PEP257 and best practice standard coding style.

Many thanks to the Biopython developers and community for making this release
possible, especially the following contributors:

- Aaron Kitzmiller (first contribution)
- Adil Iqbal (first contribution)
- Allis Tauri
- Andrew Guy
- Ariel Aptekmann (first contribution)
- Ben Fulton
- Bertrand Caron (first contribution)
- Chris Rands (first contribution)
- Connor T. Skennerton
- Eric Rasche
- Eric Talevich
- Francesco Gastaldello
- François Coste (first contribution)
- Frederic Sapet (first contribution)
- Jimmy O'Donnell (first contribution)
- Jared Andrews (first contribution)
- John Kern (first contribution)
- Jordan Willis (first contribution)
- João Rodrigues
- Kai Blin
- Markus Piotrowski
- Mateusz Korycinski (first contribution)
- Maximilian Greil
- Michiel de Hoon
- morrme (first contribution)
- Noam Kremen (first contribution)
- Patrick Kunzmann (first contribution)
- Peter Cock
- Rasmus Fonseca (first contribution)
- Rodrigo Dorantes-Gilardi (first contribution)
- Sacha Laurent (first contribution)
- Sourav Singh
- Ted Cybulski (first contribution)
- Tiago Antao
- Wibowo 'Bow' Arindrarto
- Zheng Ruan


6 April 2017: Biopython 1.69
============================

This release of Biopython supports Python 2.7, 3.3, 3.4, 3.5 and 3.6 (we have
now dropped support for Python 2.6). It has also been tested on PyPy v5.7,
PyPy3.5 v5.7 beta, and Jython 2.7.

We have started to dual-license Biopython under both our original liberal
"Biopython License Agreement", and the very similar but more commonly used
"3-Clause BSD License". In this release a small number of the Python files
are explicitly available under either license, but most of the code remains
under the "Biopython License Agreement" only. See the ``LICENSE.rst`` file
for more details.

We now expect and take advantage of NumPy under PyPy, and compile most of the
Biopython C code modules as well.

Bio.AlignIO now supports the UCSC Multiple Alignment Format (MAF) under the
format name "maf", using new module Bio.AlignIO.MafIO which also offers
indexed access to these potentially large files using SQLite3 (contributed by
Andrew Sczesnak, with additional refinements from Adam Novak).

Bio.SearchIO.AbiIO has been extended to support parsing FSA files. The
underlying format (ABIF) remains the same as AB1 files and so the string
'abif' is the expected format argument in the main SeqIO functions. AbiIO
determines whether the file is AB1 or FSA based on the presence of specific
tags.

The Uniprot parser is now able to parse "submittedName" elements in XML files.

The NEXUS parser handling of internal node comments has been improved, which
should help if working with tools like the BEAST TreeAnnotator. Slashes are
now also allowed in identifiers.

New parser for ExPASy Cellosaurus, a cell line database, cell line catalogue,
and cell line ontology (contributed by Steve Marshall).

For consistency the Bio.Seq module now offers a complement function (already
available as a method on the Seq and MutableSeq objects).

The SeqFeature object's qualifiers is now an explicitly ordered dictionary
(note that as of Python 3.6 the Python dict is ordered by default anyway).
This helps reproduce GenBank/EMBL files on input/output.

The Bio.SeqIO UniProt-XML parser was updated to cope with features with
unknown locations which can be found in mass spec data.

The Bio.SeqIO GenBank, EMBL, and IMGT parsers now record the molecule type
from the LOCUS/ID line explicitly in the record.annotations dictionary.
The Bio.SeqIO EMBL parser was updated to cope with more variants seen in
patent data files, and the related IMGT parser was updated to cope with
IPD-IMGT/HLA database files after release v3.16.0 when their ID line changed.
The GenBank output now uses colon space to match current NCBI DBLINK lines.

The Bio.Affy package supports Affymetrix version 4 of the CEL file format,
in addition to version 3.

The restriction enzyme list in ``Bio.Restriction`` has been updated to the
February 2017 release of REBASE.

Bio.PDB.PDBList now can download PDBx/mmCif (new default), PDB (old default),
PDBML/XML and mmtf format protein structures.  This is inline with the RCSB
recommendation to use PDBx/mmCif and deprecate the PDB file format. Biopython
already has support for parsing mmCif files.

Additionally, a number of small bugs have been fixed with further additions
to the test suite, and there has been further work to follow the Python PEP8,
PEP257 and best practice standard coding style.

Many thanks to the Biopython developers and community for making this release
possible, especially the following contributors:

- Aaron Rosenfeld
- Adam Kurkiewicz (first contribution)
- Adam Novak (first contribution)
- Adrian Altenhoff (first contribution)
- Allis Tauri (first contribution)
- Andrew Dalke
- Andrew Guy (first contribution)
- Andrew Sczesnak (first contribution)
- Ben Fulton
- Bernhard Thiel (first contribution)
- Bertrand Néron
- Blaise Li (first contribution)
- Brandon Carter (first contribution)
- Brandon Invergo
- Carlos Pena
- Carlos Ríos
- Chris Warth
- Emmanuel Noutahi
- Foen Peng (first contribution)
- Francesco Gastaldello (first contribution)
- Francisco Pina-Martins (first contribution)
- Hector Martinez (first contribution)
- Jacek Śmietański
- Jack Twilley (first contribution)
- Jeroen Van Goey (first contribution)
- Joshua Meyers (first contribution)
- Kurt Graff (first contribution)
- Lenna Peterson
- Leonhard Heizinger (first contribution)
- Marcin Magnus (first contribution)
- Markus Piotrowski
- Maximilian Greil (first contribution)
- Michał J. Gajda (first contribution)
- Michiel de Hoon
- Milind Luthra (first contribution)
- Oscar G. Garcia (first contribution)
- Owen Solberg
- Peter Cock
- Richard Neher (first contribution)
- Sebastian Bassi
- Sourav Singh (first contribution)
- Spencer Bliven (first contribution)
- Stefans Mezulis
- Steve Bond
- Steve Marshall (first contribution)
- Uri Laserson
- Veronika Berman (first contribution)
- Vincent Davis
- Wibowo 'Bow' Arindrarto


25 August 2016: Biopython 1.68
==============================

This release of Biopython supports Python 2.6, 2.7, 3.3, 3.4 and 3.5, but
this will be our final release to run on Python 2.6. It has also been tested
on PyPy 5.0, PyPy3 version 2.4, and Jython 2.7.

Bio.PDB has been extended to parse the RSSB's new binary Macromolecular
Transmission Format (MMTF, see http://mmtf.rcsb.org), in addition to the
mmCIF and PDB file formats (contributed by Anthony Bradley). This requires
an optional external dependency on the mmtf-python library.

Module Bio.pairwise2 has been re-written (contributed by Markus Piotrowski).
It is now faster, addresses some problems with local alignments, and also
now allows gap insertions after deletions, and vice versa, inspired by the
https://doi.org/10.1101/031500 preprint from Flouri et al.

The two sample graphical tools SeqGui (Sequence Graphical User Interface)
and xbbtools were rewritten (SeqGui) or updated (xbbtools) using the tkinter
library (contributed by Markus Piotrowski). SeqGui allows simple nucleotide
transcription, back-transcription and translation into amino acids using
Bio.Seq internally, offering of the NCBI genetic codes supported in Biopython.
xbbtools is able to open Fasta formatted files, does simple nucleotide
operations and translations in any reading frame using one of the NCBI genetic
codes. In addition, it supports standalone Blast installations to do local
Blast searches.

New NCBI genetic code table 26 (Pachysolen tannophilus Nuclear Code) has been
added to Bio.Data (and the translation functionality), and table 11 is now
also available under the alias Archaeal.

In line with NCBI website changes, Biopython now uses HTTPS rather than HTTP
to connect to the NCBI Entrez and QBLAST API.

Additionally, a number of small bugs have been fixed with further additions
to the test suite, and there has been further work to follow the Python PEP8
and best practice standard coding style.

Many thanks to the Biopython developers and community for making this release
possible, especially the following contributors:

- Anthony Bradley (first contribution)
- Ben Fulton
- Carlos Pena
- Connor T. Skennerton
- Iddo Friedberg
- Kai Blin
- Kristian Davidsen (first contribution)
- Markus Piotrowski
- Olivier Morelle (first contribution)
- Peter Cock
- Stefans Mezulis (first contribution)
- Tiago Antao
- Travis Wrightsman
- Uwe Schmitt (first contribution)
- Xiaoyu Zhuo (first contribution)


8 June 2016: Biopython 1.67
===========================

This release of Biopython supports Python 2.6, 2.7, 3.3, 3.4 and 3.5, but
support for Python 2.6 is considered to be deprecated. It has also been
tested on PyPy 5.0, PyPy3 version 2.4, and Jython 2.7.

Comparison of SeqRecord objects until now has used the default Python object
comparison (are they the same instance in memory?). This can be surprising, but
comparing all of the attributes would be too complex. As of this release
attempting to compare SeqRecord objects should raise an exception instead. If
you want the old behaviour, use id(record1) == id(record2) instead.

New experimental module Bio.phenotype is for working with Phenotype Microarray
plates in JSON and the machine vendor's CSV format (contributed by Marco
Galardini).

Following the convention used elsewhere in Biopython, there is a new function
Bio.KEGG.read(...) for parsing KEGG files expected to contain a single record
only - the existing function Bio.KEGG.parse(...) is intended to be used to
iterate over multi-record files.

When a gap character is defined, Bio.Seq will now translate gap codons
(e.g. "---") into a single gap ("-") in the protein sequence. The gap character
is inferred from the Seq object's alphabet, but it can also be passed as an
argument to the translate method.

The new NCBI genetic code table 25, covering Candidate Division SR1 and
Gracilibacteria, has been added to Bio.Data (and the translation
functionality).

The Bio.Entrez interface will automatically use an HTTP POST rather than
HTTP GET if the URL would exceed 1000 characters. This is based on NCBI
guidelines and the fact that very long queries like complex searches can
otherwise trigger an HTTP Error 414 Request URI too long.

Foreign keys are now used when creating BioSQL databases with SQLite3 (this
was not possible until SQLite version 3.6.19). The BioSQL taxonomy code now
updates the taxon table left/right keys when updating the taxonomy.

There have been some fixes to the MMCIF structure parser which now uses
identifiers which better match results from the PDB structure parse.

The restriction enzyme list in ``Bio.Restriction`` has been updated to the
May 2016 release of REBASE.

The mmCIF parser in Bio.PDB.MMCIFParser has been joined by a second version
which only looks at the ATOM and HETATM lines and can be much faster.

The Bio.KEGG.REST will now return unicode text-based handles, except for
images which remain as binary bytes-based handles, making it easier to use
with the mostly text-based parsers in Biopython.

Note that the BioSQL test configuration information is now in a new file
Tests/biosql.ini rather than directly in Tests/test_BioSQL_*.py as before.
You can make a copy of the provided example file Tests/biosql.ini.sample
as Tests/biosql.ini and edit this if you wish to run the BioSQL tests.

Additionally, a number of small bugs have been fixed with further additions
to the test suite, and there has been further work to follow the Python PEP8
standard coding style, and in converting our docstring documentation to use
the reStructuredText markup style.

Many thanks to the Biopython developers and community for making this release
possible, especially the following contributors:

- Aaron Rosenfeld (first contribution)
- Anders Pitman (first contribution)
- Barbara Mühlemann (first contribution)
- Ben Fulton
- Ben Woodcroft (first contribution)
- Brandon Invergo
- Brian Osborne (first contribution)
- Carlos Pena
- Chaitanya Gupta (first contribution)
- Chris Warth (first contribution)
- Christiam Camacho (first contribution)
- Connor T. Skennerton
- David Koppstein (first contribution)
- Eric Talevich
- Jacek Śmietański (first contribution)
- João D Ferreira (first contribution)
- João Rodrigues
- Joe Cora (first contribution)
- Kai Blin
- Leighton Pritchard
- Lenna Peterson
- Marco Galardini (first contribution)
- Markus Piotrowski
- Matt Ruffalo (first contribution)
- Matteo Sticco (first contribution)
- Nader Morshed (first contribution)
- Owen Solberg (first contribution)
- Peter Cock
- Steve Bond (first contribution)
- Terry Jones (first contribution)
- Vincent Davis
- Zheng Ruan


21 October 2015: Biopython 1.66
===============================

This release of Biopython supports Python 2.6, 2.7, 3.3, 3.4 and 3.5, but
support for Python 2.6 is considered to be deprecated. It has also been
tested on PyPy 2.4 to 2.6, PyPy3 version 2.4, and Jython 2.7.

Further work on the Bio.KEGG and Bio.Graphics modules now allows drawing KGML
pathways with transparency.

The Bio.SeqIO "abi" parser now decodes almost all the documented fields used
by the ABIF instruments - including the individual color channels.

Bio.PDB now has a QCPSuperimposer module using the Quaternion Characteristic
Polynomial algorithm for superimposing structures. This is a fast alternative
to the existing SVDSuperimposer code using singular value decomposition.

Bio.Entrez now implements the NCBI Entrez Citation Matching function
(ECitMatch), which retrieves PubMed IDs (PMIDs) that correspond to a set of
input citation strings.

Bio.Entrez.parse(...) now supports NCBI XML files using XSD schemas, which
will be downloaded and cached like NCBI DTD files.

A subtle bug in how multi-part GenBank/EMBL locations on the reverse strand
were parsed into CompoundLocations was fixed: complement(join(...)) as used
by NCBI worked, but join(complement(...),complement(...),...) as used by
EMBL/ENSEMBL gave the CompoundLocation parts in the wrong order. A related
bug when taking the reverse complement of a SeqRecord containing features
with CompoundLocations was also fixed.

Additionally, a number of small bugs have been fixed with further additions
to the test suite, and there has been further work on conforming to the
Python PEP8 standard coding style.

Many thanks to the Biopython developers and community for making this release
possible, especially the following contributors:

- Alan Medlar (first contribution)
- Anthony Mathelier (first contribution)
- Antony Lee (first contribution)
- Anuj Sharma (first contribution)
- Ben Fulton (first contribution)
- Bertrand Néron
- Brandon Invergo
- Carlos Pena
- Christian Brueffer
- Connor T. Skennerton (first contribution)
- David Arenillas (first contribution)
- David Nicholson (first contribution)
- Emmanuel Noutahi (first contribution)
- Eric Rasche (first contribution)
- Fabio Madeira (first contribution)
- Franco Caramia (first contribution)
- Gert Hulselmans (first contribution)
- Gleb Kuznetsov (first contribution)
- João Rodrigues
- John Bradley (first contribution)
- Kai Blin
- Kian Ho (first contribution)
- Kozo Nishida (first contribution)
- Kuan-Yi Li (first contribution)
- Leighton Pritchard
- Lucas Sinclair
- Michiel de Hoon
- Peter Cock
- Saket Choudhary
- Sunhwan Jo (first contribution)
- Tarcisio Fedrizzi (first contribution)
- Tiago Antao
- Vincent Davis


17 December 2014: Biopython 1.65 released.
==========================================

The Biopython sequence objects now use string comparison, rather than Python's
object comparison. This has been planned for a long time with warning messages
in place (under Python 2, the warnings were sadly missing under Python 3).

The Bio.KEGG and Bio.Graphics modules have been expanded with support for
the online KEGG REST API, and parsing, representing and drawing KGML pathways.

The Pterobranchia Mitochondrial genetic code has been added to Bio.Data (and
the translation functionality), which is the new NCBI genetic code table 24.

The Bio.SeqIO parser for the ABI capillary file format now exposes all the raw
data in the SeqRecord's annotation as a dictionary. This allows further
in-depth analysis by advanced users.

Bio.SearchIO QueryResult objects now allow Hit retrieval using its alternative
IDs (any IDs listed after the first one, for example as used with the NCBI
BLAST NR database).

We have also done some more work applying PEP8 coding styles to Biopython.

Bio.SeqUtils.MeltingTemp has been rewritten with new functionality.

The new experimental module Bio.CodonAlign has been renamed Bio.codonalign
(and similar lower case PEP8 style module names have been used for the
sub-modules within this).

Bio.SeqIO.index_db(...) and Bio.SearchIO.index_db(...) now store any relative
filenames relative to the index file, rather than (as before) relative to the
current directory at the time the index was built. This makes the indexes
less fragile, so that they can be used from other working directories. NOTE:
This change is backward compatible (old index files work as before), however
relative paths in new indexes will not work on older versions of Biopython!

Biopython also seems to work fine under PyPy3 2.4 which implements Python 3.2
plus unicode string literals.

Many thanks to the Biopython developers and community for making this release
possible, especially the following contributors:

- Alan Du (first contribution)
- Carlos Pena (first contribution)
- Colin Lappala (first contribution)
- Christian Brueffer
- David Bulger (first contribution)
- Eric Talevich
- Evan Parker (first contribution)
- Hongbo Zhu
- Kai Blin
- Kevin Wu (first contribution)
- Leighton Pritchard
- Leszek Pryszcz (first contribution)
- Markus Piotrowski
- Matt Shirley (first contribution)
- Mike Cariaso (first contribution)
- Peter Cock
- Seth Sims (first contribution)
- Tiago Antao
- Travis Wrightsman (first contribution)
- Tyghe Vallard (first contribution)
- Vincent Davis
- Wibowo 'Bow' Arindrarto
- Zheng Ruan


29 May 2014: Biopython 1.64 released.
=====================================

This release of Biopython supports Python 2.6 and 2.7, 3.3 and also the
new 3.4 version. It is also tested on PyPy 2.0 to 2.3, and Jython 2.7b2.

The new experimental module Bio.CodonAlign facilitates building codon
alignment and further analysis upon it. This work is from the Google
Summer of Code (GSoC) project by Zheng Ruan.

Bio.Phylo now has tree construction and consensus modules, from the
GSoC work by Yanbo Ye.

Bio.Entrez will now automatically download and cache new NCBI DTD files for
XML parsing under the user's home directory (using ``~/.biopython`` on
Unix like systems, and ``$APPDATA/biopython`` on Windows).

Bio.Sequencing.Applications now includes a wrapper for the samtools command
line tool.

Bio.PopGen.SimCoal now also supports fastsimcoal.

SearchIO hmmer3-text, hmmer3-tab, and hmmer3-domtab now support output from
hmmer3.1b1.

The ``accession`` of QueryResult and Hit objects created when using the
'hmmer3-tab' format are now properly named as ``accession`` (previously they
were ``acc``, deviating from the documentation).

The ``homology` key in the ``aln_annotation`` attribute of an HSP object in
Bio.SearchIO has been renamed to ``similarity``.

The Bio.SeqUtils masses and molecular_weight function have been updated.

BioSQL can now use the mysql-connector package (available for Python 2, 3
and PyPy) as an alternative to MySQLdb (Python 2 only) to connect to a MySQL
database.

Many thanks to the Biopython developers and community for making this release
possible, especially the following contributors:

- Chunlei Wu (first contribution)
- Edward Liaw (first contribution)
- Eric Talevich
- Leighton Pritchard
- Manlio Calvi (first contribution)
- Markus Piotrowski (first contribution)
- Melissa Gymrek (first contribution)
- Michiel de Hoon
- Nigel Delaney
- Peter Cock
- Saket Choudhary
- Tiago Antao
- Vincent Davis (first contribution)
- Wibowo 'Bow' Arindrarto
- Yanbo Ye (first contribution)
- Zheng Ruan (first contribution)


4 December 2013: Biopython 1.63 released.
=========================================

This release supports Python 3.3 onwards without conversion via the 2to3
library. See the Biopython 1.63 beta release notes below for details. Since
the beta release we have made some minor bug fixes and test improvements.

The restriction enzyme list in Bio.Restriction has been updated to the
December 2013 release of REBASE.

Additional contributors since the beta:

- Gokcen Eraslan (first contribution)


12 November 2013: Biopython 1.63 beta released.
===============================================

This is a beta release for testing purposes, the main reason for a
beta version is the large amount of changes imposed by the removal of
the 2to3 library previously required for the support of Python 3.X.
This was made possible by dropping Python 2.5 (and Jython 2.5).

This release of Biopython supports Python 2.6 and 2.7, and also Python
3.3.

The Biopython Tutorial & Cookbook, and the docstring examples in the source
code, now use the Python 3 style print function in place of the Python 2
style print statement. This language feature is available under Python 2.6
and 2.7 via::

    from __future__ import print_function

Similarly we now use the Python 3 style built-in next function in place of
the Python 2 style iterators' .next() method. This language feature is also
available under Python 2.6 and 2.7.

Many thanks to the Biopython developers and community for making this release
possible, especially the following contributors:

- Chris Mitchell (first contribution)
- Christian Brueffer
- Eric Talevich
- Josha Inglis (first contribution)
- Konstantin Tretyakov (first contribution)
- Lenna Peterson
- Martin Mokrejs
- Nigel Delaney (first contribution)
- Peter Cock
- Sergei Lebedev (first contribution)
- Tiago Antao
- Wayne Decatur (first contribution)
- Wibowo 'Bow' Arindrarto


28 August 2013: Biopython 1.62 released.
========================================

This is our first release to officially support Python 3, however it is
also our final release supporting Python 2.5. Specifically this release
is supported and tested on standard Python 2.5, 2.6, 2.7 and 3.3.
It was also tested under Jython 2.5, 2.7 and PyPy 1.9, 2.0.

See the Biopython 1.62 beta release notes below for most changes. Since the
beta release we have added several minor bug fixes and test improvements.
Additional contributors since the beta:

- Bertrand Néron (first contribution)
- Lenna Peterson
- Martin Mokrejs
- Matsuyuki Shirota (first contribution)


15 July 2013: Biopython 1.62 beta released.
===========================================

This is a beta release for testing purposes, both for new features added,
and changes to location parsing, but more importantly Biopython 1.62 will
be our first release to officially support Python 3.

Specifically we intend Biopython 1.62 to support standard Python 2.5, 2.6, 2.7
and 3.3, but the release will also be tested under Jython 2.5, 2.7 and PyPy
1.9, 2.0 as well. It will be our final release supporting Python 2.5.

The translation functions will give a warning on any partial codons (and this
will probably become an error in a future release). If you know you are dealing
with partial sequences, either pad with N to extend the sequence length to a
multiple of three, or explicitly trim the sequence.

The handling of joins and related complex features in Genbank/EMBL files has
been changed with the introduction of a CompoundLocation object. Previously
a SeqFeature for something like a multi-exon CDS would have a child SeqFeature
(under the sub_features attribute) for each exon. The sub_features property
will still be populated for now, but is deprecated and will in future be
removed. Please consult the examples in the help (docstrings) and Tutorial.

Thanks to the efforts of Ben Morris, the Phylo module now supports the file
formats NeXML and CDAO. The Newick parser is also significantly faster, and can
now optionally extract bootstrap values from the Newick comment field (like
Molphy and Archaeopteryx do). Nate Sutton added a wrapper for FastTree to
Bio.Phylo.Applications.

New module Bio.UniProt adds parsers for the GAF, GPA and GPI formats from
UniProt-GOA.

The BioSQL module is now supported in Jython. MySQL and PostgreSQL databases
can be used. The relevant JDBC driver should be available in the CLASSPATH.

Feature labels on circular GenomeDiagram figures now support the label_position
argument (start, middle or end) in addition to the current default placement,
and in a change to prior releases these labels are outside the features which
is now consistent with the linear diagrams.

The code for parsing 3D structures in mmCIF files was updated to use the
Python standard library's shlex module instead of C code using flex.

The Bio.Sequencing.Applications module now includes a BWA command line wrapper.

Bio.motifs supports JASPAR format files with multiple position-frequence
matrices.

Additionally there have been other minor bug fixes and more unit tests.

Many thanks to the Biopython developers and community for making this release
possible, especially the following contributors:

- Alexander Campbell (first contribution)
- Andrea Rizzi (first contribution)
- Anthony Mathelier (first contribution)
- Ben Morris (first contribution)
- Brad Chapman
- Christian Brueffer
- David Arenillas (first contribution)
- David Martin (first contribution)
- Eric Talevich
- Iddo Friedberg
- Jian-Long Huang (first contribution)
- Joao Rodrigues
- Kai Blin
- Michiel de Hoon
- Nate Sutton (first contribution)
- Peter Cock
- Petra Kubincová (first contribution)
- Phillip Garland
- Saket Choudhary (first contribution)
- Tiago Antao
- Wibowo 'Bow' Arindrarto
- Xabier Bello (first contribution)


5 February 2013: Biopython 1.61 released.
=========================================

GenomeDiagram has three new sigils (shapes to illustrate features). OCTO shows
an octagonal shape, like the existing BOX sigil but with the corners cut off.
JAGGY shows a box with jagged edges at the start and end, intended for things
like NNNNN regions in draft genomes. Finally BIGARROW is like the existing
ARROW sigil but is drawn straddling the axis. This is useful for drawing
vertically compact figures where you do not have overlapping genes.

New module Bio.Graphics.ColorSpiral can generate colors along a spiral path
through HSV color space. This can be used to make arbitrary 'rainbow' scales,
for example to color features or cross-links on a GenomeDiagram figure.

The Bio.SeqIO module now supports reading sequences from PDB files in two
different ways. The "pdb-atom" format determines the sequence as it appears in
the structure based on the atom coordinate section of the file (via Bio.PDB,
so NumPy is currently required for this). Alternatively, you can use the
"pdb-seqres" format to read the complete protein sequence as it is listed in
the PDB header, if available.

The Bio.SeqUtils module how has a seq1 function to turn a sequence using three
letter amino acid codes into one using the more common one letter codes. This
acts as the inverse of the existing seq3 function.

The multiple-sequence-alignment object used by Bio.AlignIO etc now supports
an annotation dictionary. Additional support for per-column annotation is
planned, with addition and splicing to work like that for the SeqRecord
per-letter annotation.

A new warning, Bio.BiopythonExperimentalWarning, has been introduced. This
marks any experimental code included in the otherwise stable release. Such
'beta' level code is ready for wider testing, but still likely to change and
should only be tried by early adopters to give feedback via the biopython-dev
mailing list. We'd expect such experimental code to reach stable status in
one or two releases time, at which point our normal policies about trying to
preserve backwards compatibility would apply. See also the README file.

This release also includes Bow's Google Summer of Code work writing a unified
parsing framework for NCBI BLAST (assorted formats including tabular and XML),
HMMER, BLAT, and other sequence searching tools. This is currently available
with the new BiopythonExperimentalWarning to indicate that this is still
somewhat experimental. We're bundling it with the main release to get more
public feedback, but with the big warning that the API is likely to change.
In fact, even the current name of Bio.SearchIO may change since unless you
are familiar with BioPerl its purpose isn't immediately clear.

The Bio.Motif module has been updated and reorganized. To allow for a clean
deprecation of the old code, the new motif code is stored in a new module
Bio.motifs, and a PendingDeprecationWarning was added to Bio.Motif.

A faster low level string FASTA based parser SimpleFastaParser has been added
to Bio.SeqIO.FastaIO which like its sister function for FASTQ files does not
have the overhead of constructing SeqRecord objects.

Additionally there have been other minor bug fixes and more unit tests.

Finally, we are phasing out support for Python 2.5. We will continue support
for at least one further release (Biopython 1.62). This could be extended
given feedback from our users (or if the Jython 2.7 release is delayed, since
the current stable release Jython 2.5 implemented Python 2.5 only). Focusing
on Python 2.6 and 2.7 only will make writing Python 3 compatible code easier.

Many thanks to the Biopython developers and community for making this release
possible, especially the following contributors:

- Brandon Invergo
- Bryan Lunt (first contribution)
- Christian Brueffer (first contribution)
- David Cain
- Eric Talevich
- Grace Yeo (first contribution)
- Jeffrey Chang
- Jingping Li (first contribution)
- Kai Blin (first contribution)
- Leighton Pritchard
- Lenna Peterson
- Lucas Sinclair (first contribution)
- Michiel de Hoon
- Nick Semenkovich (first contribution)
- Peter Cock
- Robert Ernst (first contribution)
- Tiago Antao
- Wibowo 'Bow' Arindrarto


25 June 2012: Biopython 1.60 released.
======================================

New module Bio.bgzf supports reading and writing BGZF files (Blocked GNU
Zip Format), a variant of GZIP with efficient random access, most commonly
used as part of the BAM file format. This uses Python's zlib library
internally, and provides a simple interface like Python's gzip library.
Using this the Bio.SeqIO indexing functions now support BGZF compressed
sequence files.

The GenBank/EMBL parser will now give a warning on unrecognised feature
locations and continue parsing (leaving the feature's location as None).
Previously it would abort with an exception, which was often unhelpful.

The Bio.PDB.MMCIFParser is now compiled by default (but is still not
available under Jython, PyPy or Python 3).

The SFF parser in Bio.SeqIO now decodes Roche 454 'universal accession
number' 14 character read names, which encode the timestamp of the run,
the region the read came from, and the location of the well.

In the Phylo module, the "draw" function for plotting tree objects has become
much more flexible, with improved support for matplotlib conventions and new
parameters for specifying branch and taxon labels. Writing in the PhyloXML
format has been updated to more closely match the output of other programs. A
wrapper for the program RAxML has been added under Bio.Phylo.Applications,
alongside the existing wrapper for PhyML.

Additionally there have been other minor bug fixes and more unit tests.

Many thanks to the Biopython developers and community for making this release
possible, especially the following contributors:

- Brandon Invergo
- Eric Talevich
- Jeff Hussmann (first contribution)
- John Comeau (first contribution)
- Kamil Slowikowski (first contribution)
- Kevin Jacobs
- Lenna Peterson (first contribution)
- Matt Fenwick (first contribution)
- Peter Cock
- Paul T. Bathen
- Wibowo Arindrarto


24 February 2012: Biopython 1.59 released.
==========================================

Please note that this release will *not* work on Python 2.4 (while the recent
releases have worked despite us not officially supporting this).

The position objects used in Bio.SeqFeature now act almost like integers,
making dealing with fuzzy locations in EMBL/GenBank files much easier. Note as
part of this work, the arguments to create fuzzy positions OneOfPosition and
WithinPosition have changed in a non-backwards compatible way.

The SeqFeature's strand and any database reference are now properties of the
FeatureLocation object (a more logical placement), with proxy methods for
backwards compatibility. As part of this change, if you print a location
object it will now display any strand and database reference information.

The installation setup.py now supports 'install_requires' when setuptools
is installed. This avoids the manual dialog when installing Biopython via
easy_install or pip and numpy is not installed. It also allows user libraries
that require Biopython to include it in their install_requires and get
automatic installation of dependencies.

Bio.Graphics.BasicChromosome has been extended to allow simple sub-features to
be drawn on chromosome segments, suitable to show the position of genes, SNPs
or other loci. Note Bio.Graphics requires the ReportLab library.

Bio.Graphics.GenomeDiagram has been extended to allow cross-links between
tracks, and track specific start/end positions for showing regions. This can
be used to imitate the output from the Artemis Comparison Tool (ACT).
Also, a new attribute circle_core makes it easier to have an empty space in
the middle of a circular diagram (see tutorial).

Bio.Align.Applications now includes a wrapper for command line tool Clustal
Omega for protein multiple sequence alignment.

Bio.AlignIO now supports sequential PHYLIP files (as well as interlaced
PHYLIP files) as a separate format variant.

New module Bio.TogoWS offers a wrapper for the TogoWS REST API, a web service
based in Japan offering access to KEGG, DDBJ, PDBj, CBRC plus access to some
NCBI, EBI resources including PubMed, GenBank and UniProt. This is much easier
to use than the NCBI Entrez API, but should be especially useful for Biopython
users based in Asia.

Bio.Entrez function efetch has been updated to handle the NCBI's stricter
handling of multiple ID arguments in EFetch 2.0, however the NCBI have also
changed the retmode default argument so you may need to make this explicit.
e.g. retmode="text"

Additionally there have been other minor bug fixes and more unit tests.

Many thanks to the Biopython developers and community for making this release
possible, especially the following contributors:

- Andreas Wilm (first contribution)
- Alessio Papini (first contribution)
- Brad Chapman
- Brandon Invergo
- Connor McCoy
- Eric Talevich
- João Rodrigues
- Konrad Förstner (first contribution)
- Michiel de Hoon
- Matej Repič (first contribution)
- Leighton Pritchard
- Peter Cock


18 August 2011: Biopython 1.58 released.
========================================

A new interface and parsers for the PAML (Phylogenetic Analysis by Maximum
Likelihood) package of programs, supporting codeml, baseml and yn00 as well
as a Python re-implementation of chi2 was added as the Bio.Phylo.PAML module.

Bio.SeqIO now includes read and write support for the SeqXML, a simple XML
format offering basic annotation support. See Schmitt et al (2011) in
Briefings in Bioinformatics, https://doi.org/10.1093/bib/bbr025

Bio.SeqIO now includes read support for ABI files ("Sanger" capillary
sequencing trace files, containing called sequence with PHRED qualities).

The Bio.AlignIO "fasta-m10" parser was updated to cope with the >>><<< lines
as used in Bill Pearson's FASTA version 3.36, without this fix the parser
would only return alignments for the first query sequence.

The Bio.AlignIO "phylip" parser and writer now treat a dot/period in the
sequence as an error, in line with the official PHYLIP specification. Older
versions of our code didn't do anything special with this character. Also,
support for "phylip-relaxed" has been added which allows longer record names
as used in RAxML and PHYML.

Of potential interest to anyone subclassing Biopython objects, any remaining
"old style" Python classes have been switched to "new style" classes. This
allows things like defining properties.

Bio.HMM's Viterbi algorithm now expects the initial probabilities explicitly.

Many thanks to the Biopython developers and community for making this release
possible, especially the following contributors:

- Aaron Gallagher (first contribution)
- Bartek Wilczynski
- Bogdan T. (first contribution)
- Brandon Invergo (first contribution)
- Connor McCoy (first contribution)
- David Cain (first contribution)
- Eric Talevich
- Fábio Madeira (first contribution)
- Hongbo Zhu
- Joao Rodrigues
- Michiel de Hoon
- Peter Cock
- Thomas Schmitt (first contribution)
- Tiago Antao
- Walter Gillett
- Wibowo Arindrarto (first contribution)


2 April 2011: Biopython 1.57 released.
======================================

Bio.SeqIO now includes an index_db() function which extends the existing
indexing functionality to allow indexing many files, and more importantly
this keeps the index on disk in a simple SQLite3 database rather than in
memory in a Python dictionary.

Bio.Blast.Applications now includes a wrapper for the BLAST+ blast_formatter
tool from NCBI BLAST 2.2.24+ or later. This release of BLAST+ added the
ability to run the BLAST tools and save the output as ASN.1 format, and then
convert this to any other supported BLAST output format (plain text, tabular,
XML, or HTML) with the blast_formatter tool. The wrappers were also updated
to include new arguments added in BLAST 2.2.25+ such as -db_hard_mask.

The SeqRecord object now has a reverse_complement method (similar to that of
the Seq object). This is most useful to reversing per-letter-annotation (such
as quality scores from FASTQ) or features (such as annotation from GenBank).

Bio.SeqIO.write's QUAL output has been sped up, and Bio.SeqIO.convert now
uses an optimised routine for FASTQ to QUAL making this much faster.

Biopython can now be installed with pip. Thanks to David Koppstein and
James Casbon for reporting the problem.

Bio.SeqIO.write now uses lower case for the sequence for GenBank, EMBL and
IMGT output.

The Bio.PDB module received several fixes and improvements, including starting
to merge João's work from GSoC 2010; consequently Atom objects now know
their element type and IUPAC mass. (The new features that use these
attributes won't be included in Biopython until the next release, though, so
stay tuned.)

The nodetype hierarchy in the Bio.SCOP.Cla.Record class is now a dictionary
(previously it was a list of key,value tuples) to better match the standard.

Many thanks to the Biopython developers and community for making this release
possible, especially the following contributors:

- Brad Chapman
- Eric Talevich
- Erick Matsen (first contribution)
- Hongbo Zhu
- Jeffrey Finkelstein (first contribution)
- Joanna & Dominik Kasprzak (first contribution)
- Joao Rodrigues
- Kristian Rother
- Leighton Pritchard
- Michiel de Hoon
- Peter Cock
- Peter Thorpe (first contribution)
- Phillip Garland
- Walter Gillett (first contribution)


26 November 2010: Biopython 1.56 released.
==========================================

This is planned to be our last release to support Python 2.4, however this
could be delayed given immediate feedback from our users (e.g. if this proves
to be a problem in combination with other libraries or a popular Linux
distribution).

Bio.SeqIO can now read and index UniProt XML files (under format name
"uniprot-xml", which was agreed with EMBOSS and BioPerl for when/if they
support it too).

Bio.SeqIO can now read, write and index IMGT files. These are a variant of
the EMBL sequence text file format with longer feature indentation.

Bio.SeqIO now supports protein EMBL files (used in the EMBL patents database
file epo_prt.dat) - previously we only expected nucleotide EMBL files.

The Bio.Seq translation methods and function will now accept an arbitrary
CodonTable object (for those of you working on very unusual organisms).

The SeqFeature object now supports len(feature) giving the length consistent
with the existing extract method. Also, it now supports iteration giving the
coordinate (with respect to the parent sequence) of each letter within the
feature (in frame aware order), and "in" which allows you to check if a
(parent based) coordinate is within the feature location.

Bio.Entrez will now try to download any missing NCBI DTD files and cache them
in the user's home directory.

The provisional database schema for BioSQL support on SQLite which Biopython
has been using since Release 1.53 has now been added to BioSQL, and updated
slightly.

Bio.PopGen.FDist now supports the DFDist command line tool as well as FDist2.

Bio.Motif now has a chapter in the Tutorial.

(At least) 13 people have contributed to this release, including 6 new people:

- Andrea Pierleoni (first contribution)
- Bart de Koning (first contribution)
- Bartek Wilczynski
- Bartosz Telenczuk (first contribution)
- Cymon Cox
- Eric Talevich
- Frank Kauff
- Michiel de Hoon
- Peter Cock
- Phillip Garland (first contribution)
- Siong Kong (first contribution)
- Tiago Antao
- Uri Laserson (first contribution)


31 August 2010: Biopython 1.55 released.
========================================

See the notes below for the Biopython 1.55 beta release for changes since
Biopython 1.54 was released. Since the beta release we have marked a few
modules as obsolete or deprecated, and removed some deprecated code. There
have also been a few bug fixes, extra unit tests, and documentation
improvements.

(At least) 12 people have contributed to this release, including 6 new people:

- Andres Colubri (first contribution)
- Carlos Ríos (first contribution)
- Claude Paroz (first contribution)
- Cymon Cox
- Eric Talevich
- Frank Kauff
- Joao Rodrigues (first contribution)
- Konstantin Okonechnikov (first contribution)
- Michiel de Hoon
- Nathan Edwards (first contribution)
- Peter Cock
- Tiago Antao


18 August 2010: Biopython 1.55 beta released.
=============================================

This is a beta release for testing purposes, both for new features added,
and more importantly updates to avoid code deprecated in Python 2.7 or in
Python 3. This is an important step towards Python 3 support.

We are phasing out support for Python 2.4. We will continue to support it
for at least one further release (Biopython 1.56). This could be delayed
given feedback from our users (e.g. if this proves to be a problem in
combination with other libraries or a popular Linux distribution).

The SeqRecord object now has upper and lower methods (like the Seq object and
Python strings), which return a new SeqRecord with the sequence in upper or
lower case and a copy of all the annotation unchanged.

Several small issues with Bio.PDB have been resolved, which includes better
handling of model numbers, and files missing the element column.

Feature location parsing for GenBank and EMBL files has been rewritten,
making the parser much faster.

Ace parsing by SeqIO now uses zero rather than None for the quality score of
any gaps (insertions) in the contig sequence.

The BioSQL classes DBServer and BioSeqDatabase now act more like Python
dictionaries, making it easier to count, delete, iterate over, or check for
membership of namespaces and records.

The command line tool application wrapper classes are now executable, so you
can use them to call the tool (using the subprocess module internally) and
capture the output and any error messages as strings (stdout and stderr).
This avoids having to worry about the details of how best to use subprocess.

(At least) 10 people have contributed to this release, including 5 new people:

- Andres Colubri (first contribution)
- Carlos Ríos (first contribution)
- Claude Paroz (first contribution)
- Eric Talevich
- Frank Kauff
- Joao Rodrigues (first contribution)
- Konstantin Okonechnikov (first contribution)
- Michiel de Hoon
- Peter Cock
- Tiago Antao


May 20, 2010: Biopython 1.54 released.
======================================

See the notes below for the Biopython 1.54 beta release for changes since
Biopython 1.53 was released. Since then there have been some changes to
the new Bio.Phylo module, more documentation, and a number of smaller
bug fixes.


April 2, 2010: Biopython 1.54 beta released.
============================================

We are phasing out support for Python 2.4. We will continue to support it
for at least two further releases, and at least one year (whichever takes
longer), before dropping support for Python 2.4. This could be delayed
given feedback from our users (e.g. if this proves to be a problem in
combination with other libraries or a popular Linux distribution).

New module Bio.Phylo includes support for reading, writing and working with
phylogenetic trees from Newick, Nexus and phyloXML files. This was work by
Eric Talevich on a Google Summer of Code 2009 project, under The National
Evolutionary Synthesis Center (NESCent), mentored by Brad Chapman and
Christian Zmasek.

Bio.Entrez includes some more DTD files, in particular eLink_090910.dtd,
needed for our NCBI Entrez Utilities XML parser.

The parse, read and write functions in Bio.SeqIO and Bio.AlignIO will now
accept filenames as well as handles. This follows a general shift from
other Python libraries, and does make usage a little simpler. Also
the write functions will now accept a single SeqRecord or alignment.

Bio.SeqIO now supports writing EMBL files (DNA and RNA sequences only).

The dictionary-like objects from Bio.SeqIO.index() now support a get_raw
method for most file formats, giving you the original unparsed data from the
file as a string. This is useful for selecting a subset of records from a
file where Bio.SeqIO.write() does not support the file format (e.g. the
"swiss" format) or where you need to exactly preserve the original layout.

Based on code from Jose Blanca (author of sff_extract), Bio.SeqIO now
supports reading, indexing and writing Standard Flowgram Format (SFF)
files which are used by 454 Life Sciences (Roche) sequencers. This means
you can use SeqIO to convert from SFF to FASTQ, FASTA and QUAL (as
trimmed or untrimmed reads).

An improved multiple sequence alignment object has been introduced,
and is used by Bio.AlignIO for input. This is a little stricter than the
old class but should otherwise be backwards compatible.

(At least) 11 people contributed to this release, including 5 new people:

- Anne Pajon (first contribution)
- Brad Chapman
- Christian Zmasek
- Diana Jaunzeikare (first contribution)
- Eric Talevich
- Jose Blanca (first contribution)
- Kevin Jacobs (first contribution)
- Leighton Pritchard
- Michiel de Hoon
- Peter Cock
- Thomas Holder (first contribution)


December 15, 2009: Biopython 1.53 released.
===========================================

Biopython is now using git for source code control, currently on github. Our
old CVS repository will remain on the OBF servers in the short/medium term
as a backup, but will not be updated in future.

The Bio.Blast.Applications wrappers now covers the new NCBI BLAST C++ tools
(where blastall is replaced by blastp, blastn, etc, and the command line
switches have all been renamed). These will be replacing the old wrappers in
Bio.Blast.NCBIStandalone which are now obsolete, and will be deprecated in
our next release.

The plain text BLAST parser has been updated, and should cope with recent
versions of NCBI BLAST, including the new C++ based version. Nevertheless,
we (and the NCBI) still recommend using the XML output for parsing.

The Seq (and related UnknownSeq) objects gained upper and lower methods,
like the string methods of the same name but alphabet aware. The Seq object
also gained a new ungap method for removing gap characters in an alphabet
aware manner.

The SeqFeature object now has an extract method, used with the parent
sequence (as a string or Seq object) to get the region of that sequence
described by the feature's location information (including the strand and
any sub-features for a join). As an example, this is useful to get the
nucleotide sequence for features in GenBank or EMBL files.

SeqRecord objects now support addition, giving a new SeqRecord with the
combined sequence, all the SeqFeatures, and any common annotation.

Bio.Entrez includes the new (Jan 2010) DTD files from the NCBI for parsing
MedLine/PubMed data.

The NCBI codon tables have been updated from version 3.4 to 3.9, which adds
a few extra start codons, and a few new tables (Tables 16, 21, 22 and 23).
Note that Table 14 which used to be called "Flatworm Mitochondrial" is now
called "Alternative Flatworm Mitochondrial", and "Flatworm Mitochondrial" is
now an alias for Table 9 ("Echinoderm Mitochondrial").

The restriction enzyme list in Bio.Restriction has been updated to the
Nov 2009 release of REBASE.

The Bio.PDB parser and output code has been updated to understand the
element column in ATOM and HETATM lines (based on patches contributed by
Hongbo Zhu and Frederik Gwinner). Bio.PDB.PDBList has also been updated
for recent changes to the PDB FTP site (Paul T. Bathen).

SQLite support was added for BioSQL databases (Brad Chapman), allowing access
to BioSQL through a lightweight embedded SQL engine. Python 2.5+ includes
support for SQLite built in, but on Python 2.4 the optional sqlite3 library
must be installed to use this. We currently use a draft BioSQL on SQLite
schema, which will be merged with the main BioSQL release for use in other
projects.

Support for running Biopython under Jython (using the Java Virtual Machine)
has been much improved thanks to input from Kyle Ellrott. Note that Jython
does not support C code - this means NumPy isn't available, and nor are a
selection of Biopython modules (including Bio.Cluster, Bio.PDB and BioSQL).
Also, currently Jython does not parse DTD files, which means the XML parser
in Bio.Entrez won't work. However, most of the Biopython modules seem fine
from testing Jython 2.5.0 and 2.5.1.

(At least) 12 people contributed to this release, including 3 first timers:

- Bartek Wilczynski
- Brad Chapman
- Chris Lasher
- Cymon Cox
- Frank Kauff
- Frederik Gwinner (first contribution)
- Hongbo Zhu (first contribution)
- Kyle Ellrott
- Leighton Pritchard
- Michiel de Hoon
- Paul Bathen (first contribution)
- Peter Cock


September 22, 2009: Biopython 1.52 released.
============================================

The Population Genetics module now allows the calculation of several tests,
and statistical estimators via a wrapper to GenePop. Supported are tests for
Hardy-Weinberg equilibrium, linkage disequilibrium and estimates for various
F statistics (Cockerham and Wier Fst and Fis, Robertson and Hill Fis, etc),
null allele frequencies and number of migrants among many others. Isolation
By Distance (IBD) functionality is also supported.

New helper functions Bio.SeqIO.convert() and Bio.AlignIO.convert() allow an
easier way to use Biopython for simple file format conversions. Additionally,
these new functions allow Biopython to offer important file format specific
optimisations (e.g. FASTQ to FASTA, and interconverting FASTQ variants).

New function Bio.SeqIO.index() allows indexing of most sequence file formats
(but not alignment file formats), allowing dictionary like random access to
all the entries in the file as SeqRecord objects, keyed on the record id.
This is especially useful for very large sequencing files, where all the
records cannot be held in memory at once. This supplements the more flexible
but memory demanding Bio.SeqIO.to_dict() function.

Bio.SeqIO can now write "phd" format files (used by PHRED, PHRAD and CONSED),
allowing interconversion with FASTQ files, or FASTA+QUAL files.

Bio.Emboss.Applications now includes wrappers for the "new" PHYLIP EMBASSY
package (e.g. fneighbor) which replace the "old" PHYLIP EMBASSY package (e.g.
eneighbor) whose Biopython wrappers are now obsolete.

See also the DEPRECATED file, as several old deprecated modules have finally
been removed (e.g. Bio.EUtils which had been replaced by Bio.Entrez).

On a technical note, this will be the last release using CVS for source code
control. Biopython is moving from CVS to git.


August 17, 2009: Biopython 1.51 released.
=========================================

FASTQ support in Bio.SeqIO has been improved, extended and sped up since
Biopython 1.50. Support for Illumina 1.3+ style FASTQ files was added in the
1.51 beta release. Furthermore, we now follow the interpretation agreed on
the OBF mailing lists with EMBOSS, BioPerl, BioJava and BioRuby for inter-
conversion and the valid score range for each FASTQ variant. This means
Solexa FASTQ scores can be from -5 to 62 (format name "fastq-solexa" in
Bio.SeqIO), Illumina 1.3+ FASTQ files have PHRED scores from 0 to 62 (format
name "fastq-illumina"), and Sanger FASTQ files have PHRED scores from 0 to
93 (format name "fastq" or "fastq-sanger").

Bio.Sequencing.Phd has been updated, for example to cope with missing peak
positions. The "phd" support in Bio.SeqIO has also been updated to record
the PHRED qualities (and peak positions) in the SeqRecord's per-letter
annotation. This allows conversion of PHD files into FASTQ or QUAL which may
be useful for meta-assembly.

See the notes below for the Biopython 1.50 beta release for changes since
Biopython 1.49 was released. This includes dropping support for Python 2.3,
removing our deprecated parsing infrastructure (Martel and Bio.Mindy), and
hence removing any dependence on mxTextTools.

Additionally, since the beta, a number of small bugs have been fixed, and
there have been further additions to the test suite and documentation.


June 23, 2009: Biopython 1.51 beta released.
============================================

Biopython no longer supports Python 2.3.  Currently we support Python 2.4,
2.5 and 2.6.

Our deprecated parsing infrastructure (Martel and Bio.Mindy) has been
removed.  This means Biopython no longer has any dependence on mxTextTools.

A few cosmetic issues in GenomeDiagram with arrow sigils and labels on
circular diagrams have been fixed.

Bio.SeqIO will now write GenBank files with the feature table (previously
omitted), and a couple of obscure errors parsing ambiguous locations have
been fixed.

Bio.SeqIO can now read and write Illumina 1.3+ style FASTQ files (which use
PHRED quality scores with an ASCII offset of 64) under the format name
"fastq-illumina". Biopython 1.50 supported just "fastq" (the original Sanger
style FASTQ files using PHRED scores with an ASCII offset of 33), and
"fastq-solexa" (the original Solexa/Illumina FASTQ format variant holding
Solexa scores with an ASCII offset of 64) .

For parsing the "swiss" format, Bio.SeqIO now uses the new Bio.SwissProt
parser, making it about twice as fast as in Biopython 1.50, where the older
now deprecated Bio.SwissProt.SProt was used. There should be no functional
differences as a result of this change.

Our command line wrapper objects have been updated to support accessing
parameters via python properties, and setting of parameters at initiation
with keyword arguments.  Additionally Cymon Cox has contributed several new
multiple alignment wrappers under Bio.Align.Applications.

A few more issues with Biopython's BioSQL support have been fixed (mostly by
Cymon Cox). In particular, the default PostgreSQL schema includes some rules
intended for BioPerl support only, which were causing problems in Biopython
(see BioSQL bug 2839).

There have also been additions to the tutorial, such as the new alignment
wrappers, with a whole chapter for the SeqRecord object. We have also added
to the unit test coverage.


April 20, 2009: Biopython 1.50 released.
========================================

See the notes below for the Biopython 1.50 beta release for more details,
but the highlights are:

* The SeqRecord supports slicing and per-letter-annotation
* Bio.SeqIO can read and write FASTQ and QUAL files
* Bio.Seq now has an UnknownSeq object
* GenomeDiagram has been integrated into Biopython
* New module Bio.Motif will later replace Bio.AlignAce and Bio.MEME
* This will be the final release to support Python 2.3
* This will be the final release with Martel and Bio.Mindy

Since the 1.50 beta release:

* The NCBI's Entrez EFetch no longer supports rettype="genbank"
  and "gb" (or "gp") should be used instead.
* Bio.SeqIO now supports "gb" as an alias for "genbank".
* The Seq object now has string-like startswith and endswith methods
* Bio.Blast.NCBIXML now has a read function for single record files
* A few more unit tests were added
* More documentation


April 3, 2009: Biopython 1.50 beta released.
============================================

The SeqRecord object has a new dictionary attribute, letter_annotations,
which is for holding per-letter-annotation information like sequence
quality scores or secondary structure predictions.  As part of this work,
the SeqRecord object can now be sliced to give a new SeqRecord covering
just part of the sequence.  This will slice the per-letter-annotation to
match, and will also include any SeqFeature objects as appropriate.

Bio.SeqIO can now read and write FASTQ and QUAL quality files using PHRED
quality scores (Sanger style, also used for Roche 454 sequencing), and FASTQ
files using Solexa/Illumina quality scores.

The Bio.Seq module now has an UnknownSeq object, used for when we have a
sequence of known length, but unknown content.  This is used in parsing
GenBank and EMBL files where the sequence may not be present (e.g. for a
contig record) and when parsing QUAL files (which don't have the sequence)

GenomeDiagram by Leighton Pritchard has been integrated into Biopython as
the Bio.Graphics.GenomeDiagram module  If you use this code, please cite the
publication Pritchard et al. (2006), Bioinformatics 22 616-617.  Note that
like Bio.Graphics, this requires the ReportLab python library.

A new module Bio.Motif has been added, which is intended to replace the
existing Bio.AlignAce and Bio.MEME modules.

The set of NCBI DTD files included with Bio.Entrez has been updated with the
revised files the NCBI introduced on 1 Jan 2009.

Minor fix to BioSQL for retrieving references and comments.

Bio.SwissProt has a new faster parser which will be replacing the older
slower code in Bio.SwissProt.SProt (which we expect to deprecate in the next
release).

We've also made some changes to our test framework, which is now given a
whole chapter in the tutorial.  This intended to help new developers or
contributors wanting to improve our unit test coverage.


November 21, 2008: Biopython 1.49 released.
===========================================

See the notes below for the Biopython 1.49 beta release for more details,
but the highlights are:

* Biopython has transitioned from Numeric to NumPy
* Martel and Bio.Mindy are now deprecated

Since the 1.49 beta release:

* A couple of NumPy issues have been resolved
* Further small improvements to BioSQL
* Bio.PopGen.SimCoal should now work on Windows
* A few more unit tests were added


November 7, 2008: Biopython 1.49 beta released.
===============================================

Biopython has transitioned from Numeric to NumPy.  Please move to NumPy.

A number of small changes have been made to support Python 2.6 (mostly
avoiding deprecated functionality), and further small changes have been
made for better compatibility with Python 3 (this work is still ongoing).
However, we intend to support Python 2.3 for only a couple more releases.

As part of the Numeric to NumPy migration, Bio.KDTree has been rewritten in
C instead of C++ which therefore simplifies building Biopython from source.

Martel and Bio.Mindy are now considered to be deprecated, meaning mxTextTools
is no longer required to use Biopython.  See the DEPRECATED file for details
of other deprecations.

The Seq object now supports more string like methods (gaining find, rfind,
split, rsplit, strip, lstrip and rstrip in addition to previously supported
methods like count).  Also, biological methods transcribe, back_transcribe
and translate have been added, joining the pre-existing reverse_complement
and complement methods.  Together these changes allow a more object
orientated programming style using the Seq object.

The behaviour of the Bio.Seq module's translate function has changed so that
ambiguous codons which could be a stop codon like "TAN" or "NNN" are now
translated as "X" (consistent with EMBOSS and BioPerl - Biopython previously
raised an exception), and a bug was fixed so that invalid codons (like "A-T")
now raise an exception (previously these were translated as stop codons).

BioSQL had a few bugs fixed, and can now optionally fetch the NCBI taxonomy
on demand when loading sequences (via Bio.Entrez) allowing you to populate
the taxon/taxon_name tables gradually.  This has been tested in combination
with the BioSQL load_ncbi_taxonomy.pl script used to populate or update the
taxon/taxon_name tables.  BioSQL should also now work with the psycopg2
driver for PostgreSQL as well as the older psycopg driver.

The PDB and PopGen sections of the Tutorial have been promoted to full
chapters, and a new chapter has been added on supervised learning methods
like logistic regression.  The "Cookbook" section now has a few graphical
examples using Biopython to calculate sequence properties, and matplotlib
(pylab) to plot them.

The input functions in Bio.SeqIO and Bio.AlignIO now accept an optional
argument to specify the expected sequence alphabet.

The somewhat quirky unit test GUI has been removed, the unit tests are now
run via the command line by default.


September 8, 2008: Biopython 1.48 released.
===========================================

The SeqRecord and Alignment objects have a new method to format the object as
a string in a requested file format (handled via Bio.SeqIO and Bio.AlignIO).

Additional file formats supported in Bio.SeqIO and Bio.AlignIO:

- reading and writing "tab" format (simple tab separated)
- writing "nexus" files.
- reading "pir" files (NBRF/PIR)
- basic support for writing "genbank" files (GenBank plain text)

Fixed some problems reading Clustal alignments (introduced in Biopython 1.46
when consolidating Bio.AlignIO and Bio.Clustalw).

Updates to the Bio.Sequencing parsers.

Bio.PubMed and the online code in Bio.GenBank are now considered obsolete,
and we intend to deprecate them after the next release. For accessing PubMed
and GenBank, please use Bio.Entrez instead.

Bio.Fasta is now considered to be obsolete, please use Bio.SeqIO instead. We
do intend to deprecate this module eventually, however, for several years
this was the primary FASTA parsing module in Biopython and is likely to be in
use in many existing scripts.

Martel and Bio.Mindy are now considered to be obsolete, and are likely to be
deprecated and removed in a future release.

In addition a number of other modules have been deprecated, including:
Bio.MetaTool, Bio.EUtils, Bio.Saf, Bio.NBRF, and Bio.IntelliGenetics
See the DEPRECATED file for full details.


July 5, 2008: Biopython 1.47 released.
======================================

Improved handling of ambiguous nucleotides in Bio.Seq.Translate().
Better handling of stop codons in the alphabet from a translation.
Fixed some codon tables (problem introduced in Biopython 1.46).

Updated Nexus file handling.

Fixed a bug in Bio.Cluster potentially causing segfaults in the
single-linkage hierarchical clustering library.

Added some DTDs to be able to parse EFetch results from the
nucleotide database.

Added IntelliGenetics/MASE parsing to Bio.SeqIO (as the "ig" format).


June 29, 2008: Biopython 1.46 released.
=======================================

Bio.Entrez now has several Entrez format XML parsers, and a chapter
in the tutorial.

Addition of new Bio.AlignIO module for working with sequence alignments
in the style introduced with Bio.SeqIO in recent releases, with a whole
chapter in the tutorial.

A problem parsing certain EMBL files was fixed.

Several minor fixes were made to the NCBI BLAST XML parser, including
support for the online version 2.2.18+ introduced in May 2008.

The NCBIWWW.qblast() function now allows other programs (blastx, tblastn,
tblastx) in addition to just blastn and blastp.

Bio.EUtils has been updated to explicitly enforce the NCBI's rule of at
most one query every 3 seconds, rather than assuming the user would obey
this.

Iterators in Bio.Medline, Bio.SCOP, Bio.Prosite, Bio.Prosite.Prodoc,
Bio.SwissProt, and others to make them more generally usable.

Phylip export added to Bio.Nexus.

Improved handling of ambiguous nucleotides and stop codons in
Bio.Seq.Translate (plus introduced a regression fixed in Biopython 1.47).


March 22, 2008: Biopython 1.45 released.
========================================

The Seq and MutableSeq objects act more like python strings, in particular
str(object) now returns the full sequence as a plain string.  The existing
tostring() method is preserved for backwards compatibility.

BioSQL has had some bugs fixed, and has an additional unit test which loads
records into a database using Bio.SeqIO and then checks the records can be
retrieved correctly.  The DBSeq and DBSeqRecord classes now subclass the
Seq and SeqRecord classes, which provides more functionality.

The modules under Bio.WWW are being deprecated.
Functionality in Bio.WWW.NCBI, Bio.WWW.SCOP, Bio.WWW.InterPro and
Bio.WWW.ExPASy is now available from Bio.Entrez, Bio.SCOP, Bio.InterPro and
Bio.ExPASy instead. Bio.Entrez was used to fix a nasty bug in Bio.GenBank.

Tiago Antao has included more functionality in the Population Genetics
module, Bio.PopGen.

The Bio.Cluster module has been updated to be more consistent with other
Biopython code.

The tutorial has been updated, including devoting a whole chapter to
Swiss-Prot, Prosite, Prodoc, and ExPASy. There is also a new chapter on
Bio.Entrez.

Bio.biblio was deprecated.


October 28, 2007: Biopython 1.44 released.
==========================================

NOTE: This release includes some rather drastic code changes, which were
necessary to get Biopython to work with the new release of mxTextTools.

The (reverse)complement functions in Bio.Seq support ambiguous nucleotides.

Bio.Kabat, which was previously deprecated, is now removed from Biopython.

Bio.MarkupEditor was deprecated, as it does not appear to have any users.

Bio.Blast.NCBI.qblast() updated with more URL options, thanks to a patch
from Chang Soon Ong.

Several fixes to the Blast parser.

The deprecated Bio.Blast.NCBIWWW functions blast and blasturl were removed.

The standalone Blast functions blastall, blastpgp now create XML output by
default.

Bio.SeqIO.FASTA and Bio.SeqIO.generic have been deprecated in favour of
the new Bio.SeqIO module.

Bio.FormatIO has been removed (a gradual deprecation was not possible).
Please look at Bio.SeqIO for sequence input/output instead.

Fix for a bug in Bio.Cluster, which caused kcluster() to hang on some
platforms.

Bio.expressions has been deprecated.

Bio.SeqUtils.CheckSum created, including new methods from Sebastian Bassi,
and functions crc32 and crc64 which were moved from Bio/crc.py.
Bio.crc is now deprecated. Bio.lcc was updated and moved to Bio.SeqUtils.lcc.

Bio.SwissProt parser updated to cope with recent file format updates.

Bio.Fasta, Bio.KEGG and Bio.Geo updated to pure python parsers which
don't rely on Martel.

Numerous fixes in the Genbank parser.

Several fixes in Bio.Nexus.

Bio.MultiProc and Bio.Medline.NLMMedlineXML were deprecating, as they failed
on some platforms, and seemed to have no users. Deprecated concurrent
behavior in Bio.config.DBRegistry and timeouts in Bio.dbdefs.swissprot,
which relies on Bio.MultiProc.

Tiago Antao has started work on a Population Genetics module, Bio.PopGen

Updates to the tutorial, including giving Bio.Seq and Bio.SeqIO a whole
chapter each.


March 17, 2007: Biopython 1.43 released.
========================================

New Bio.SeqIO module for reading and writing biological sequence files
in various formats, based on SeqRecord objects.  This includes a new fasta
parser which is much faster than Bio.Fasta, particularly for larger files.
Easier to use, too.

Various improvements in Bio.SeqRecord.

Running Blast using Bio.Blast.NCBIStandalone now generates output in XML
format by default.
The new function Bio.Blast.NCBIXML.parse can parse multiple Blast records
in XML format.

Bio.Cluster no longer uses ranlib, but uses its own random number generator
instead. Some modifications to make Bio.Cluster more compatible with the new
NumPy (we're not quite there yet though).

New Bio.UniGene parser.

Numerous improvements in Bio.PDB.

Bug fixes in Bio.SwissProt, BioSQL, Bio.Nexus, and other modules.

Faster parsing of large GenBank files.

New EMBL parser under Bio.GenBank and also integrated into (new) Bio.SeqIO

Compilation of KDTree (C++ code) is optional (setup.py asks the user if it
should be compiled). For the Windows installer, C++ code is now included.

Nominating Bio.Kabat for removal.

Believe it or not, even the documentation was updated.


July 16, 2006: Biopython 1.42 released.
=======================================

Bio.GenBank: New parser by Peter, which doesn't rely on Martel.

Numerous updates in Bio.Nexus and Bio.Geo.

Bio.Cluster became (somewhat) object-oriented.

Lots of bug fixes, and updates to the documentation.


October 28, 2005: Biopython 1.41 released.
==========================================

Major changes:

NEW: Bio.MEME -- thanks to Jason Hackney

Added transcribe, translate, and reverse_complement functions to Bio.Seq that
work both on Seq objects and plain strings.

Major code optimization in cpairwise2module.

CompareACE support added to AlignAce.

Updates to Blast parsers in Bio.Blast, in particular use of the XML parser
in NCBIXML contributed by Bertrand Frottier, and the BLAT parser by Yair
Benita.

Pairwise single-linkage hierarchical clustering in Bio.Cluster became much
faster and memory-efficient, allowing clustering of large data sets.

Bio.Emboss: Added command lines for einverted and palindrome.

Bio.Nexus: Added support for StringIO objects.

Numerous updates in Bio.PDB.

Lots of fixes in the documentation.

March 29, 2005: MEME parser added. Thanks to Jason Hackney


Feb 18, 2005: Biopython 1.40 beta
=================================
Major Changes since v1.30. For a full list of changes please see the CVS

IMPORTANT: Biopython now works with Python version >= 2.3

NEW: Bio.Nexus -- thanks to Frank Kauff
Bio.Nexus is a Nexus file parser. Nexus is a common format for phylogenetic
trees.

NEW: CAPS module -- Thanks to Jonathan Taylor.

NEW: Restriction enzyme package contributed by Frederic Sohm. This includes
classes for manipulating enzymes, updating from Rebase, as well as
documentation and Tests.

CHANGED: Bio.PDB -- thanks to Thomas Hamelryck.

- Added atom serial number.
- Epydoc style documentation.
- Added secondary structure support (through DSSP).
- Added Accessible Surface Area support (through DSSP).
- Added Residue Depth support (through MSMS).
- Added Half Sphere Exposure.
- Added Fragment classification of the protein backbone (see Kolodny et al.,
- JMB, 2002).
- Corrected problem on Windows with PDBList (thanks to Matt Dimmic)
- Added StructureAlignment module to superimpose structures based on a FASTA
  sequence alignment.
- Various additions to Polypeptide.
- Various bug corrections in Vector.
- Lots of smaller bug corrections and additional features

CHANGED: MutableSeq -- thanks to Michiel De Hoon
Added the functions 'complement' and 'reverse_complement' to Bio.Seq's Seq and
MutableSeq objects. Similar functions previously existed in various locations
in BioPython:

- forward_complement, reverse_complement in Bio.GFF.easy
- complement, antiparallel in Bio.SeqUtils

These functions have now been deprecated, and will issue a DeprecationWarning
when used. The functions complement and reverse_complement, when applied to a
Seq object, will return a new Seq object. The same function applied to a
MutableSeq object will modify the MutableSeq object itself, and don't return
anything.


May 14, 2004: Biopython 1.30
============================

- Affy package added for dealing with Affymetrix cel files -- thanks to Harry
  Zuzan.
- Added code for parsing Blast XML output -- thanks to Bertrand Frottier.
- Added code for parsing Compass output -- thanks to James Casbon.
- New melting temperature calculation module -- thanks to Sebastian Bassi.
- Added lowess function for non-parameteric regression -- thanks to Michiel.
- Reduced protein alphabet supported added -- thanks to Iddo.

- Added documentation for Logistic Regression and Bio.PDB -- thanks to Michiel
  and Thomas.
- Documentation added for converting between file formats.
- Updates to install documentation for non-root users -- thanks to Jakob
  Fredslund.
- epydoc now used for automatic generation of documentation.

- Fasta parser updated to use Martel for parsing and indexing, allowing better
  speed and dealing with large data files.
- Updated to Registry code. Now 'from Bio import db' gives you a number of new
  retrieval options, including embl, fasta, genbak, interpro, prodoc and
  swissprot.
- GenBank parser uses new Martel format. GenBank retrieval now uses EUtils
  instead of the old non-working entrez scripts. GenBank indexing uses standard
  Mindy indexing. Fix for valueless qualifiers in feature keys -- thanks to
  Leighton Pritchard.
- Numerous updated to Bio.PDB modules -- thanks to Thomas. PDB can now parse
  headers -- thanks to Kristian Rother.
- Updates to the Ace parser -- thanks to Frank Kauff and Leighton Pritchard.

- Added pgdb (PyGreSQL) support to BioSQL -- thanks to Marc Colosimo.
- Fix problems with using py2exe and Biopython -- thanks to Michael Cariaso.
- PSIBlast parser fixes -- thanks to Jer-Yee John Chuang and James Casbon.
- Fix to NCBIWWW retrieval so that HTML results are returned correctly.
- Fix to Clustalw to handle question marks in title names -- thanks to Ashleigh
  Smythe.
- Fix to NBRF parsing to it accepts files produced by Clustalw -- thanks to
  Ashleigh Smythe.
- Fixes to the Enyzme module -- thanks to Marc Colosimo.
- Fix for bugs in SeqUtils -- thanks to Frank Kauff.
- Fix for optional hsps in ncbiblast Martel format -- thanks to Heiko.
- Fix to Fasta parsing to allow # comment lines -- thanks to Karl Diedrich.
- Updates to the C clustering library -- thanks to Michiel.
- Fixes for breakage in the SCOP module and addition of regression tests to
  framework -- thanks to Gavin.
- Various fixes to Bio.Wise -- thanks to Michael.
- Fix for bug in FastaReader -- thanks to Micheal.
- Fix EUtils bug where efetch would only return 500 sequences.
- Updates for Emboss commandlines, water and tranalign.
- Fixes to the FormatIO system of file conversion.

- C++ code (KDTree, Affy) now compiled by default on most platforms -- thanks
  to Michael for some nice distutils hacks and many people for testing.
- Deprecated Bio.sequtils -- use Bio.SeqUtils instead.
- Deprecated Bio.SVM -- use libsvm instead.
- Deprecated Bio.kMeans and Bio.xkMeans -- use Bio.cluster instead.
- Deprecated RecordFile -- doesn't appear to be finished code.


Feb 16, 2004: Biopython 1.24
============================

- New parsers for Phred and Ace format files -- thanks to Frank Kauff
- New Code for dealing with NMR data -- thanks to Bob Bussell
- New SeqUtils modules for codon usage, isoelectric points and other
  protein properties -- thanks to Yair Benita
- New code for dealing with Wise contributed by Michael
- EZ-Retrieve sequence retrieval now supported thanks to Jeff
- Bio.Cluster updated along with documentation by Michiel
- BioSQL fixed so it now works with the current SQL schema -- thanks to Yves
  Bastide for patches
- Patches to Bio/__init__ to make it compatible with py2exe -- thanks to
  Leighton Pritchard
- Added __iter__ to all Biopython Iterators to make them Python 2.2 compatible
- Fixes to NCBIWWW for retrieving from NCBI -- thanks to Chris Wroe
- Retrieval of multiple alignment objects from BLAST records -- thanks to
  James Casbon
- Fixes to GenBank format for new tags by Peter
- Parsing fixes in clustalw parsed -- thanks to Greg Singer and Iddo
- Fasta Indexes can have a specified filename -- thanks to Chunlei Wu
- Fix to Prosite parser -- thanks to Mike Liang
- Fix in GenBank parsing -- mRNAs now get strand information


Oct 18, 2003: Biopython 1.23
============================

- Fixed distribution of files in Bio/Cluster
- Now distributing Bio/KDTree/_KDTree.swig.C
- minor updates in installation code
- added mmCIF support for PDB files


Oct 9, 2003: Biopython 1.22
===========================

- Added Peter Slicker's patches for speeding up modules under Python 2.3
- Fixed Martel installation.
- Does not install Bio.Cluster without Numeric.
- Distribute EUtils DTDs.
- Yves Bastide patched NCBIStandalone.Iterator to be Python 2.0 iterator
- Ashleigh's string coercion fixes in Clustalw.
- Yair Benita added precision to the protein molecular weights.
- Bartek updated AlignAce.Parser and added Motif.sim method
- bug fixes in Michiel De Hoon's clustering library
- Iddo's bug fixes to Bio.Enzyme and new RecordConsumer
- Guido Draheim added patches for fixing import path to xbb scripts
- regression tests updated to be Python 2.3 compatible
- GenBank.NCBIDictionary is smarter about guessing the format


Jul 28, 2003: Biopython 1.21
============================

- Martel added back into the released package
- new AlignACE module by Bartek Wilczynski
- Andreas Kuntzagk fix for GenBank Iterator on empty files


Jul 27, 2003: Biopython 1.20
============================

- added Andrew Dalke's EUtils library
- added Michiel de Hoon's gene expression analysis package
- updates to setup code, now smarter about dependencies
- updates to test suite, now smarter about code that is imported
- Michael Hoffman's fixes to DocSQL
- syntax fixes in triemodule.c to compile on SGI, Python 2.1 compatible
- updates in NCBIStandalone, short query error
- Sebastian Bassi submitted code to calculate LCC complexity
- Greg Kettler's NCBIStandalone fix for long query lengths
- slew of miscellaneous fixes from George Paci
- miscellaneous cleanups and updates from Andreas Kuntzagk
- Peter Bienstman's fixes to Genbank code -- now parses whole database
- Kayte Lindner's LocusLink package
- miscellaneous speedups and code cleanup in ParserSupport by Brad Chapman
- miscellaneous BLAST fixes and updates
- Iddo added new code to parse BLAST table output format
- Karl Diedrich's patch to read T_Coffee files
- Larry Heisler's fix for primer3 output
- Bio.Medline now uses proper iterator objects
- copen now handles SIGTERM correctly
- small bugfixes and updates in Thomas Hamelryck's PDB package
- bugfixes and updates to SeqIO.FASTA reader
- updates to Registry system, conforms to 2003 hackathon OBDA spec
- Yu Huang patch to support tblastn in wublast expression


Dec 17, 2002: Biopython 1.10
============================

- Python requirement bumped up to 2.2
- hierarchy reorg, many things moved upwards into Bio namespace
- pairwise2 replaces fastpairwise and pairwise
- removed deprecated Sequence.py package
- minor bug fix in File.SGMLStripper
- added Scripts/debug/debug_blast_parser.py to diagnose blast parsing errors
- IPI supported by SwissProt/SProt.py parser
- large speedup for kmeans
- new registry framework for generic access to databases and parsers
- small bug fix in stringfns.split
- scripts that access NCBI moved over to new EUtils system
- new crc module
- biblio.py supports the EBI Bibliographic database
- new CDD parser
- new Ndb parser
- new ECell parser
- new Geo parser
- access to GFF databases
- new KDTree data structure
- new LocusLink parser
- new MarkovModel algorithm
- new Saf parser
- miscellaneous sequence handling functions in sequtils
- new SVDSuperimpose algorithm


Dec 18, 2001: Biopython1.00a4
=============================

- minor bug fix in NCBIStandalone.blastall
- optimization in dynamic programming code
- new modules for logistic regression and maximum entropy
- minor bug fix in ParserSupport
- minor bug fixes in SCOP package
- minor updates in the kMeans cluster selection code
- minor bug fixes in SubsMat code
- support for XML-formatted MEDLINE files
- added MultiProc.run to simplify splitting code across processors
- listfns.items now supports lists with unhashable items
- new data type for pathways
- new support for intelligenetics format
- new support for metatool format
- new support for NBRF format
- new support for generalized launching of applications
- new support for genetic algorithms
- minor bug fixes in GenBank parsing
- new support for Primer in the Emboss package
- new support for chromosome graphics
- new support for HMMs
- new support for NeuralNetwork
- slew of Martel fixes (see Martel docs)


Sept 3, 2001: Biopython1.00a3
=============================

- added package to support KEGG
- added sequtils module for computations on sequences
- added pairwise sequence alignment algorithm
- major bug fixes in UndoHandle
- format updates in PubMed
- Tk interface to kMeans clustering


July 5, 2001: Biopython1.00a2
=============================

- deprecated old regression testing frameworks
- deprecated Sequence.py
- Swiss-Prot parser bug fixes
- GenBank parser bug fixes
- Can now output GenBank format
- can now download many sequences at a time from GenBank
- kMeans clustering algorithm
- Kabat format now supported
- FSSP format now supported
- more functionality for alignment code
- SubsMat bug fixes and updates
- fixed memory leak in listfns bug fixes
- Martel bundled and part of the install procedure
- Medline.Parser bug fixes
- PubMed.download_many handles broken IDs better


Mar 3, 2001: Biopython 1.00a1
=============================

- Refactoring of modules.  X/X.py moved to X/__init__.py.
- Can search sequences for Prosite patterns at ExPASy
- Can do BLAST searches against stable URL at NCBI
- Prosite Pattern bug fixes
- GenBank parser
- Complete Seq and SeqFeatures framework
- distutils cleanup
- compile warning cleanups
- support for UniGene
- code for working with substitution matrices
- Tools.MultiProc package for rudimentary multiprocessing stuff


Nov 10, 2000: Biopython 0.90d04
===============================

- Added support for multiple alignments, ClustalW
- BLAST updates, bug fixes, and BlastErrorParser
- Fixes for PSI-BLAST in master-slave mode
- Minor update in stringfns, split separators can be negated
- Added download_many function to PubMed
- xbbtools updates
- Prodoc parser now accepts a copyright at the end of a record
- Swiss-Prot parser now handles taxonomy ID tag


Sept 6, 2000: Biopython 0.90d03
===============================

- Blast updates:

  - bug fixes in NCBIStandalone, NCBIWWW
  - some __str__ methods in Record.py implemented (incomplete)

- Tests:

  - new BLAST regression tests
  - prosite tests fixed

- New parsers for Rebase, Gobase
- pure python implementation of C-based tools
- Thomas Sicheritz-Ponten's xbbtools
- can now generate documentation from docstrings using HappyDoc


Aug17-18, 2000: Bioinformatics Open Source Conference 2000
==========================================================

We had a very good Birds-of-a-Feather meeting:
http://mailman.open-bio.org/pipermail/biopython/2000-August/000360.html


Aug 2, 2000: Biopython 0.90d02 is released.
===========================================

- Blast updates:
  - now works with v2.0.14
  - HSP.identities and HSP.positives now tuples
  - HSP.gaps added
- SCOP updates:
  - Lin.Iterator now works with release 50
- Starting a tutorial
- New regression tests for Prodoc


July 6, 2000: Biopython 0.90d01 is released.
============================================


February 8, 2000: Anonymous CVS made available.
===============================================


August 1999: Biopython project founded.
=======================================

Call for Participation sent out to relevant mailing lists, news
groups.

The Biopython Project (https://www.biopython.org/) is a new open
collaborative effort to develop freely available Python libraries and
applications that address the needs of current and future work in
bioinformatics, including sequence analysis, structural biology,
pathways, expression data, etc.  When available, the source code will
be released as open source (https://github.com/biopython/biopython/blob/9c4785fc9eaf8a3bc436c6c0b16e7a05019cade1/LICENSE)
under terms similar to Python.

This is a Call for Participation for interested people to join the
project.  We are hoping to attract people from a diverse set of
backgrounds to help with code development, site maintenance,
scientific discussion, etc.  This project is open to everyone.  If
you're interested, please visit the web page, join the biopython
mailing list, and let us know what you think!

Jeffrey Chang <jchang@smi.stanford.edu>
Andrew Dalke <dalke@bioreason.com>
