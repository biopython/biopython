This file provides documentation for modules in Biopython that have been moved
or deprecated in favor of other modules. This provides some quick and easy
to find documentation about how to update your code to work again.
Python releases go first, then code (modules, methods, functions).

Python
======

Python 3.9
----------
First supported in release 1.79, although it was mostly working in 1.78.

Python 3.8
----------
First supported in release 1.75. Support deprecated as of Release 1.83.

Python 3.7
----------
No longer supported as of Release 1.82. First supported in release 1.73.

Python 3.6
----------
No longer supported as of Release 1.80, having triggered a deprecation
warning as of release 1.79. First supported in release 1.69.

Python 3.5
----------
No longer supported as of Release 1.77. First supported in release 1.66.

Python 3.4
----------
No longer supported as of Release 1.75, having triggered a deprecation
warning in release 1.74. First supported in release 1.64.

Python 3.3
----------
No longer supported as of Release 1.70, having triggered a warning with
release 1.67 onwards.

Python 3.0, 3.1, 3.2
--------------------
Never officially supported, these triggered a warning in Release 1.62
recommending Python 3.3 or later. As of Biopython Release 1.63 onwards,
installation simply aborts with a error message.

Python 2.7
----------
No longer supported as of Release 1.77 (2020, in line with end-of-life or
sunset date for Python 2.7 itself), having triggered a warning in prior
releases.

Jython
------
No longer supported as of Release 1.77 with the end of Python 2 support.
Biopython was mostly working under Jython 2.7.0, but support for Jython
was deprecated as of Release 1.70.

Biopython modules, methods, functions
=====================================

Bio.AlignInfo
-------------
The ``pos_specific_score_matrix`` method of the ``SummaryInfo`` class and the
``PSSM`` class were deprecated in release 1.82. As an alternative, please use
the ``alignment`` property of a ``MultipleSeqAlignment`` object to obtains a
new-style ``Alignment`` object, and use it to create a ``Bio.motifs.Motif``
object. For example,

>>> alignment = msa.alignment
>>> from Bio.motifs import Motif
>>> motif = Motif('ACGT', alignment)
>>> counts = motif.counts

The ``counts`` object contains the same information as the PSSM returned by
``pos_specific_score_matrix``, but note that the indices are reversed:

>>> counts[letter][i] == pssm[index][letter]
True

The ``information_content`` method and the ``ic_vector`` attribute of the
``SummaryInfo`` class were deprecated in release 1.82. As an alternative,
please use the ``relative_entropy`` attribute of the ``motif`` instance (see
above); it contains the same values as the ``ic_vector`` attribute, while
``sum(relative_entropy)`` is equal to the value returned by
``information_content``.

The ``replacement_dictionary`` method of the ``SummaryInfo`` class was
deprecated in release 1.82. As an alternative, please use the ``alignment``
property of the ``MultipleSeqAlignment`` object to obtain a new-style
``Alignment`` object, and use its ``substitutions`` attribute to obtain the
replacement dictionary:

>>> alignment = msa.alignment
>>> dictionary = alignment.substitutions

If the multiple sequence alignment object ``msa`` was obtained using
``Bio.AlignIO``, then you can obtain a new-style ``Alignment`` object directly
by using ``Bio.Align.read`` instead of ``Bio.AlignIO.read``, or
``Bio.Align.parse`` instead of ``Bio.AlignIO.parse``.

The ``dumb_consensus`` and ``gap_consensus`` methods of the ``SummaryInfo``
class were deprecated in Release 1.82.

The ``print_info_content`` function in ``Bio.Align.AlignInfo`` was deprecated
in Release 1.82.

Bio.kNN
-------
Deprecated in release 1.82, consider using scikit-learn instead.

Bio.LogisticRegression
----------------------
Deprecated in release 1.82, consider using scikit-learn instead.

Bio.NaiveBayes
--------------
Deprecated in release 1.82, consider using skikit-learn instead.

Bio.MaxEntropy
--------------
Deprecated in release 1.82, consider using scikit-learn instead.

Bio.MarkovModel
---------------
Deprecated in release 1.82, consider using hmmlearn
(https://pypi.org/project/hmmlearn/) instead.

Bio.HMM
-------
The `Bio.HMM.DynamicProgramming`, `Bio.HMM.Trainer`, `Bio.HMM.MarkovModel`, and
`Bio.HMM.Utilities` modules were deprecated in release 1.82. Consider using
hmmlearn (https://pypi.org/project/hmmlearn/) instead.


Bio.Data.SCOPData
-----------------
Declared obsolete in release 1.80, and removed in release 1.82. Please use
Bio.Data.PDBData instead.

Bio.Application and the command line wrappers using it
------------------------------------------------------
Declared obsolete in release 1.79, and deprecated in release 1.82. Please use
the standard library subprocess module directly instead.

Bio.Index
---------
Deprecated in release 1.75, removed in release 1.77. Was not used anywhere in
Biopython.

Bio.Crystal
-----------
Declared obsolete in release 1.75, deprecated in release 1.76, removed in
release 1.79. PDB NDB files can be opened with Bio.PDB.

Bio.motifs
----------
``Bio.motifs.mast`` plain-text parsing deprecated in favor of XML parsing as of
release 1.74. Also affects ``Bio.motifs.read`` and ``Bio.motifs.parse`` for the
``mast`` format.
The ``format`` method of the ``Motif`` class in ``Bio.motifs`` was deprecated
in release 1.77, in favor of a ``__format__`` method that can be used from the
``format`` built-in function. This decision was reversed in release 1.79.
The ``search`` method of the ``Instances`` class in ``Bio.motifs`` was
deprecated in release 1.82. Instead of ``instances.search(sequence)``,
``sequence.search(instances)`` can be used, where sequence is a Seq object.
This allows instances to have different lengths.

The ``Instances`` class and the ``instances`` argument of the ``Motif`` class
initializer in ``Bio.motifs`` were deprecated in release 1.82. Instead of

>>> from Bio.motifs import Instances
>>> instances = Instances([Seq('ACGT'), Seq('ACCT'), Seq('AAGT')])
>>> motif = Motif(alphabet='ACGT', instances=instances)

please use

>>> from Bio.Align import Alignment
>>> alignment = Alignment([Seq('ACGT'), Seq('ACCT'), Seq('AAGT')])
>>> motif = Motif(alphabet='ACGT', alignment=alignment)

The ``instances`` attribute of the ``Motif`` class  in ``Bio.motifs`` was
deprecated in release 1.82. Instead of ``mymotif.instances``, please use
``mymotif.alignment.sequences``.


Bio.Restriction.RanaConfig
--------------------------
Removed in Biopython 1.74 without explicit deprecation period. RanaConfig was
a configuration file containing some constants for Bio.Restriction.PrintFormat
and ranacompiler.py, a script to update Bio.Restriction.Restriction_Dictionary,
and which is not part of the Biopython installation. The constants were
implemented in the respective modules.

Bio.Alphabet
------------
Declared obsolete in Biopython release 1.74, and removed from Biopython in
release 1.78. This module defined an ``Alphabet`` class and various subclasses,
which were used as attributes to ``Seq`` and objects to describe how the
individual characters in the sequence string should be interpreted. For
example, a string "AGTACACTGGT" could be a DNA sequence or a protein sequence
that happens to be rich in Alanines, Glycines, Cysteines and Threonines.
However, as the exact definition of the alphabet and its purpose remained
unclear, this class was removed from Biopython.
Starting with Biopython 1.78, the molecule type, if specified in the input
file, is stored by the ``SeqIO`` parser as ``molecule_type`` in the annotations
of each ``SeqRecord``. We urge users to use this attribute with caution, as the
molecule type means different things in different sequence file formats, and in
a sense the interpretation of ``molecule_type`` can still be ambiguous.


Bio.ExPASy.sprot_search_ful and ExPASy.sprot_search_de
------------------------------------------------------
These two functions were labelled as broken in Release 1.70, and removed in
Release 1.73, since the underlying web-server API no longer exists.

Bio.GA
------
This was deprecated in Biopython 1.70, and removed in Release 1.73.
Please consider using a dedicated genetic algorithm library like DEAP
instead.

Bio.NeuralNetwork
-----------------
This was deprecated in Biopython 1.70, and removed in Release 1.73.
Please consider using a dedicated machine learning library like
scikit-learn or TensorFlow instead.

Bio.Phylo.CDAOIO.CDAOError
--------------------------
This exception was deprecated as of Release 1.70 as it was no longer used
within Biopython, and removed in Release 1.75.

Bio.DocSQL
----------
This was deprecated in Biopython 1.69, and removed in Release 1.71.

Bio.CodonAlign
--------------
This new experimental module included in Biopython 1.64 was renamed to
Bio.codonalign in Biopython 1.65 to follow PEP8 module naming rules.

Bio.SeqRecord equality
----------------------
As of Release 1.67, the SeqRecord objects (and their subclasses) no longer use
the default Python object comparison. Instead they will raise an exception if
you try to compare them.

For backward compatibility and/or to explicitly use object comparison, please
use id(record1) == id(record2) instead.

Otherwise please test whichever specific attributes you are interested in
explicitly, for example record1.id == record2.id or record1.seq == record.seq
(see also the note below about sequence equality).

Bio.Seq sequence equality
-------------------------
As of Release 1.65, the Seq and MutableSeq objects (and their subclasses)
use string-like equality testing and hashing (ignoring any difference in
alphabet except to issue warnings).

Prior releases used Python's object comparison. Warnings of this change
were first added in Release 1.54 (May 2010), with hash warnings present
from Release 1.62 (August 2013) to Release 1.76 (December 2019).

For backward compatibility and/or to silence warnings about this, please use
explicit string comparison, str(seq1) == str(seq2), or object comparison,
id(seq1) == id(seq2), as required.

Bio.Seq.Seq.tostring() and Bio.Seq.MutableSeq.tostring()
--------------------------------------------------------
Deprecated in release 1.64, and removed in release 1.73.
You should now use str(Bio.Seq.Seq) or str(Bio.Seq.MutableSeq) instead of
the tostring() methods.

Bio.Seq.Seq.tomutable() and Bio.Seq.MutableSeq.toseq()
------------------------------------------------------
Deprecated in release 1.79, removed in release 1.81.
Instead of myseq.tomutable() or mymutableseq.toseq(), you should now use
Bio.Seq.MutableSeq(myseq) or Bio.Seq.Seq(mymutableseq), respectively.

Bio.Seq.Seq.ungap()
-------------------
Declared obsolete in release 1.79, deprecated in release 1.80, and removed in
release 1.82.  Instead of myseq.ungap(), please use myseq.replace("-", "").

Bio.Seq.UnknownSeq
------------------
Deprecated in release 1.79, and removed in release 1.81.
Instead of ``UnknownSeq(length)``, please use ``Seq(None, length=length)``.
Note that the sequence contents of a ``Seq`` object constructed in this way
is considered to be unknown, and any attempt to access the sequence contents
(for example, by calling ``print`` on the object) will result in an
``UndefinedSequenceError``.

Bio.Seq: Functions and methods ``complement`` and ``reverse_complement``
------------------------------------------------------------------------
Starting from release 1.82, the ``inplace`` argument of ``complement`` and
``reverse_complement`` in ``Bio.Seq`` always default to ``False`` both for
``Seq`` and ``MutableSeq`` objects.
To modify a ``MutableSeq`` in-place, use ``inplace=True``.

Iterator .next() methods
------------------------
The .next() method defined for any Biopython iterator is deprecated as of
Biopython 1.63 under Python 2 (and not present on Python 3). Please replace
my_iterator.next() with next(my_iterator) using the new built-in function
next() instead. Python 2 support and the remaining next methods were removed
in release 1.77.

Bio.SVDSuperimposer
-------------------
As of Release 1.63, the main class (confusingly also called) SVDSuperimposer
is best imported as follows:

>>> from Bio.SVDSuperimposer import SVDSuperimposer
>>> super_imposer = SVDSuperimposer()

This short form also works on older releases. The longer even more
confusing historical alternatives dependent on the double module name
no longer work, e.g. you can no longer do this:

>>> from Bio.SVDSuperimposer.SVDSuperimposer import SVDSuperimposer
>>> super_imposer = SVDSuperimposer()

Bio.PDB.Vector (the module)
---------------------------
Due to a long standing name shadowing problem, ``Bio.PDB.Vector`` was
both a class and a module, which defined the class and various other
functions imported to the ``Bio.PDB`` namespace.

As of Release 1.70, the module has been renamed ``Bio.PDB.vectors``, leaving
``Bio.PDB.Vector`` to unambiguously mean the class. This is in line with the
PEP8 naming conventions. A deprecated compatibility stub was left in place
so that any imports via the old module name will work but raise a warning.
This compatibility stub was removed in Release 1.74.

We expect this to have no impact for the majority of users, unless you do
something like ``from Bio.PDB.Vector import calc_dihedral`` in which case
use ``from Bio.PDB import calc_dihedral`` (which will work on older versions
of Biopython as well).

Bio.PDB.Residue
---------------
The ``get_atom`` and ``sort`` methods of the ``Residue`` class were deprecated
in Release 1.71 and 1.70 respectively, and removed in Release 1.79.

Bio.PDB.ResidueDepth
--------------------
Use of the ``PDB_TO_XYZR`` bash script was removed from ``get_surface`` in
Release 1.79.

Bio.PDB.QCPSuperimposer
-----------------------
The ``Bio.PDB.QCPSuperimposer`` module was deprecated in release 1.80, and
removed in release 1.82. Please use the ``Bio.PDB.qcprot`` module instead.

Bio.SeqFeature
--------------

Release 1.82 unfortunately removed the ``.strand``, ``.ref``, and ``.ref_db``
attributes of the ``SeqFeature`` without a deprecation period. Release 1.83
restored but deprecated them. Please use ``.location.strand`` etc instead.

With the introduction of the CompoundLocation in Release 1.62, the SeqFeature
attribute sub_features was deprecated. It was removed in Release 1.68.

Note that in Release 1.80 the location_operator argument can no longer be
used, instead do this via the CompoundLocation object. The location_operator
argument was removed from the SeqFeature initializer in Release 1.82.

There were multiple deprecations in Release 1.80, listed below. The
deprecated code was removed in Release 1.82.

* Class ``FeatureLocation`` renamed to ``SimpleLocation``, with the old
  name preserved for now solely for backward compatibility.
* Arguments ``strand``, ``ref`` and ``ref_db`` to the ``SeqFeature``
  class - set them via the location object
* Unused class ``PositionGap`` - originally for very old GenBank files.
* Location attributes ``location.nofuzzy_start`` and ``location.nofuzzy_end`` -
  use the location directly or if required ``int(location.start)`` and
  ``int(location.end)``. This will fail for the ``UnknownPosition``
  where the nofuzzy aliases returned ``None``.
* Position attribute ``.position`` returned the (left) position as an
  integer - use the location directly or if required ``int(position)``,
  however for ``OneOfPosition``, ``BetweenPosition``, and
  ``WithinPosition`` that will give the default position rather than
  the left-most (minimum) value.
* Position attribute ``.extension`` returned the "width", typically
  zero except for ``OneOfPosition``, ``BetweenPosition``, and
  ``WithinPosition`` where this must be handled explicitly now.
* Base class ``AbstractPosition`` was renamed to ``Position``.

Bio.Motif
---------
Declared obsolete with a PendingDeprecationWarning in Release 1.61, formally
deprecated in Release 1.62, removed in Release 1.67. Please use the newer
Bio.motifs module instead.

AlignAceCommandline and CompareAceCommandline
---------------------------------------------
Deprecated in release 1.62, removed in Release 1.67. An up to date version of
the software cannot be obtained anymore (affects Bio.Motif and its replacement
Bio.motifs).

Bio.SeqIO.Interfaces
--------------------
Unused class InterlacedSequenceIterator was deprecated in Release 1.61, and
removed in Release 1.64.

Class SequentialSequenceWriter was declared obsolete in Release 1.77,
deprecated in Release 1.78, and removed in Release 1.80.

Bio.HotRand
-----------
Obsolete file Bio/HotRand.py was deprecated in Release 1.61, and removed in
Release 1.64. Consider using an alternative RNG, or the Python module
"randomdotorg".

Bio.Search
----------
Long obsolete file Bio/Search.py was deprecated in Release 1.61, and removed
in Release 1.64.

Bio.Blast.NCBIStandalone
------------------------
The three functions for calling the "legacy" NCBI BLAST command line tools
blastall, blastpgp and rpsblast were declared obsolete in Biopython Release
1.53, deprecated in Release 1.61, and removed in Release 1.64. Please use
the BLAST+ wrappers in Bio.Blast.Applications instead.

The remainder of this module is a parser for the plain text BLAST output,
which was declared obsolete in Release 1.54, and deprecated in Release 1.63.
The module was removed in Release 1.72 from the public API. It lives now
in maintenance mode in Bio.SearchIO._legacy to preserve existing functionality.
A BiopythonDeprecationWarning was added to this module in Release 1.80.
The Bio.SearchIO._legacy module was removed from Biopython in Release 1.82.

For some time now, both the NCBI and Biopython have encouraged people to
parse the XML output instead.

Bio.Blast.Applications
----------------------
NCBI "legacy" BLAST tool wrappers FastacmdCommandline, BlastallCommandline,
BlastpgpCommandline and RpsBlastCommandline were declared obsolete in Release
1.53, deprecated in Release 1.61, and removed in Release 1.64, having been
replaced with wrappers for the new NCBI BLAST+ tools (e.g.
NcbiblastpCommandline and NcbipsiblastCommandline).

Bio.Blast.ParseBlastTable
-------------------------
The parser in ``Bio.Blast.ParseBlastTable`` for tabular output generated by
NCBI blastpgp was deprecated in Biopython release 1.80, and removed in release
1.82. To parse tabular output generated by BLAST programs, please use the
``parse`` function in ``Bio.Align``.

BioSQL.BioSeqDatabase
---------------------
The ``remove_database`` and ``get_all_primary_ids`` methods were removed from
the ``DBServer`` class in Release 1.79.
The ``get_Seq_by_primary_id`` method was removed from the ``BioSeqDatabase``
class in Release 1.79.

Bio.Graphics.GenomeDiagram and colour/color, centre/center
----------------------------------------------------------
GenomeDiagram originally used colour and centre (UK spelling of color and
center) for argument names.  As part of its integration into Biopython 1.50,
this will support both colour and color, and both centre and center, to help
people port existing scripts written for the standalone version of
GenomeDiagram.  However, these were deprecated in Release 1.55 final.
Support for centre was removed in Release 1.62, and we intend to eventually
remove support for colour in later releases of Biopython.

Bio.Seq, Bio.MutableSeq and the data property
---------------------------------------------
Direct use of the Seq object (and MutableSeq object) .data property is
deprecated.  As of Release 1.49, writing to the Seq object's .data property
triggered a warning, and this property was made read only in Release 1.53. In
Release 1.55 final, accessing the .data property of a Seq object gives a
DeprecationWarning. The Seq object's .data property was removed in Release
1.61.  Starting from Release 1.78, accessing the .data property of a MutableSeq
object similarly gives a deprecation warning.

Bio.SeqUtils
------------
Function quick_FASTA_reader was declared obsolete in Release 1.61,
deprecated in Release 1.64, and removed in Release 1.67. Use function
list(SimpleFastaParser(handle)) from Bio.SeqIO.FastaIO instead (but
ideally convert your code to using an iterator approach).

The 'title2ids' argument to FastaIterator in Bio.SeqIO.FastaIO and
FastqPhredIterator in Bio.SeqIO.QualityIO was deprecated in Release 1.80, and
removed in Release 1.82.
Please use a generator function to modify the records returned by the parser.

Function Tm_staluc in Bio.SeqUtils.MeltingTemp was deprecated in Release 1.78,
and removed in Release 1.80.

The modules Bio.SeqUtils.CodonUsage and Bio.SeqUtils.CodonUsageIndices were
deprecated in Release 1.80, and removed in Release 1.82. Please use the new
CodonAdaptationIndex class in Bio.SeqUtils instead. Note that this class has
been updated to use modern Python, and may give slightly different results from
the CodonAdaptationIndex class in Bio.SeqUtils.CodonUsage, as the calculation
was updated to be consistent with the calculated values by Sharp & Li.

Function 'GC' in Bio.SeqUtils was deprecated in Release 1.80, and removed in
Release 1.82. Instead use function 'gc_fraction'.

Bio.PopGen.Async
----------------
``Bio.PopGen.Async`` was deprecated in Release 1.68, removed in Release 1.70.

Bio.PopGen.FDist
----------------
``Bio.PopGen.FDist`` was deprecated in Release 1.68, removed in Release 1.70.

Bio.PopGen.SimCoal
------------------
``Bio.PopGen.SimCoal`` was deprecated in Release 1.68, and removed in Release
1.70.

Bio.UniGene
-----------
Submodule Bio.UniGene.UniGene which was an HTML parser was declared obsolete
in Release 1.59, deprecated in Release 1.61, and removed in Release 1.64.

Bio.SubsMat
-----------
The methods ``print_full_mat`` and ``print_mat`` were removed from the
`SeqMat`` class in Bio.SubsMat in Release 1.79.
The Bio.SubsMat module was deprecated in Release 1.78, and removed in Release
1.80. As an alternative, please consider using Bio.Align.substitution_matrices.

Bio.Align
---------
The ``get_column`` method of the MultipleSeqAlignment was deprecated in
Release 1.57 and removed in Release 1.69.

The ``add_sequence`` method of the MultipleSeqAlignment was deprecated in
Release 1.57 and should have been removed in Release 1.69. It was actually
removed in Release 1.79.

The ``format`` method of the MultipleSeqAlignment class and the
PairwiseAlignment class were deprecated in Release 1.76. This decision was
reversed in Release 1.79.

The ``__format__`` method of the Array class in Bio.Align.substitution_matrices
was deprecated in Release 1.79.

The PairwiseAlignment class was deprecated in Release 1.80, and removed in
Release 1.82. Please use the new Alignment class instead.

Bio.Align.Generic
-----------------
This module which defined to original (Multiple-Sequence) Alignment class was
deprecated in Release 1.57 and removed in Release 1.69.

Bio.ParserSupport
-----------------
``Bio.ParserSupport`` was declared obsolete in Release 1.59, and deprecated in
Release 1.63. The Martel specific ``EventGenerator`` was removed in Release
1.67, and the entire module was removed in Release 1.72.

Bio.KDTree
----------
This module was declared obsolete in Release 1.72, deprecated in Release 1.74,
and removed in Release 1.77. As of Release 1.72, KDTree data structures and
the functionality previously available in ``Bio.KDTree`` are provided in a new
module ``Bio.PDB.kdtrees``.

Bio.trie, Bio.triefind
----------------------
These modules were declared obsolete in Release 1.72, deprecated in Release
1.73, and removed in Release 1.77. We suggest pygtrie as an alternative library
implementing a trie data structure.

Bio.Statistics
--------------
This module was declared obsolete in Release 1.74, deprecated in Release 1.76,
and removed in Release 1.79.

Bio.File
--------
The UndoHandle class was deprecated in Release 1.77, and moved to
Bio/SearchIO/_legacy/ParserSupport.py, which was the only module in
Biopython still using this class. The UndoHandle class in Bio.File was removed
in Release 1.79.

Bio.FSSP
-----------
Deprecated in Release 1.77, and removed in Release 1.79.

Bio.Phylo._utils
----------------
The ``draw_graphviz`` function was removed in Release 1.79.

Bio.pairwise2
-------------
The ``Bio.pairwise2`` module was deprecated in Release 1.80.

Bio.Wise
--------
The ``Bio.Wise`` module was deprecated in Release 1.80, and removed in Release
1.82.

Scripts/Restriction/ranacompiler.py
-----------------------------------
The ``is_palindrom`` function was removed in Release 1.79.
