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
First supported in release 1.75.

Python 3.7
----------
First supported in release 1.73.

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

Python 2.6
----------
No longer supported as of Release 1.69, having triggered a warning with
release 1.66 onwards.

Python 2.5
----------
No longer supported as of Release 1.63, having triggered a warning with
Release 1.62, with advance notice in the release notes for Release 1.61.

Python 2.4
----------
No longer supported as of Release 1.59, having triggered a warning since
Release 1.55, with advance notice in the release notes for Release 1.54.

Jython
------
No longer supported as of Release 1.77 with the end of Python 2 support.
Biopython was mostly working under Jython 2.7.0, but support for Jython
was deprecated as of Release 1.70.

Biopython modules, methods, functions
=====================================

Bio.Data.SCOPData
-----------------
Declared obsolete in release 1.80. Please use Bio.Data.PDBData instead.

Bio.Application and the command line wrappers using it
------------------------------------------------------
Declared obsolete in release 1.79. Please use the standard library subprocess
module directly instead.

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
Deprecated in release 1.79.
Instead of myseq.tomutable() or mymutableseq.toseq(), you should now use
Bio.Seq.MutableSeq(myseq) or Bio.Seq.Seq(mymutableseq), respectively.

Bio.Seq.Seq.ungap()
-------------------
Declared obsolete in release 1.79, and deprecated in release 1.80.
Instead of myseq.ungap(), please use myseq.replace("-", "").

Bio.Seq.UnknownSeq
------------------
Deprecated in release 1.79.
Instead of ``UnknownSeq(length)``, please use ``Seq(None, length=length)``.
Note that the sequence contents of a ``Seq`` object constructed in this way
is considered to be unknown, and any attempt to access the sequence contents
(for example, by calling ``print`` on the object) will result in an
``UndefinedSequenceError``.

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

Bio.PDB.mmCIF
-------------
This was removed in Release 1.62, when MMCIF2Dict was updated to use shlex
from the standard library instead. This had required manual intervention to
include when installing Biopython from source due to a dependency on flex.

Bio.PDB.Residue
---------------
The ``get_atom`` and ``sort`` methods of the ``Residue`` class were deprecated
in Release 1.71 and 1.70 respectively, and removed in Release 1.79.

Bio.PDB.ResidueDepth
--------------------
Use of the ``PDB_TO_XYZR`` bash script was removed from ``get_surface`` in
Release 1.79.

Bio.SeqFeature
--------------
With the introduction of the CompoundLocation in Release 1.62, the SeqFeature
attribute sub_features was deprecated. It was removed in Release 1.68.

Note that in Release 1.80 the location_operator argument can no longer be
used, instead do this via the CompoundLocation object.

There were multiple deprecations in Release 1.80:

* Class ``FeatureLocation`` renamed to ``SimpleLocation``, with the old
  name preserved for now solely for backard compatibility.
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

Bio.Motif
---------
Declared obsolete with a PendingDeprecationWarning in Release 1.61, formally
deprecated in Release 1.62, removed in Release 1.67. Please use the newer
Bio.motifs module instead.

Before this, ``CompareAceParser`` and ``CompareAceConsumer`` from
``Bio.Motif.Parsers.AlignAce`` were declared obsolete in Release 1.53,
deprecated in Release 1.55 final, and removed in Release 1.57.

``AlignAceConsumer``, ``AlignAceParser``, and ``AlignAceScanner`` were
declared obsolete in Release 1.53 and deprecated in Release 1.55 final;
their functionality is now available through a read() function in
``Bio.Motif.Parsers.AlignAce``.

``MEMEParser``, ``_MEMEScanner``, ``_MEMEConsumer``, ``_MASTConsumer``,
``MASTParser``, ``_MASTScanner``, and ``MASTRecord`` were declared obsolete in
Release 1.54 and deprecated in Release 1.55 final; their functionality is now
available through a ``read()`` function in ``Bio.Motif.Parsers.MEME`` and
``Bio.Motif.Parsers.MAST``, respectively.

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

Bio.Pathway.Rep.HashSet
-----------------------
Deprecated in Release 1.59, removed in Release 1.62. Use Python's built in
set object.

Bio.SeqFeature.WithinPosition and OneOfPosition
-----------------------------------------------
The arguments to create these fuzzy positions changed in Release 1.59.

Bio.Encodings
-------------
Explicitly declared obsolete in Release 1.55, deprecated in Release 1.56, and
removed in Release 1.57.

Bio.PropertyManager
-------------------
Explicitly declared obsolete in Release 1.55, deprecated in Release 1.56, and
removed in Release 1.57.

Bio.InterPro
------------
This module was a parser for the EBI InterPro webpages, but no longer worked
with their current website. Deprecated in Release 1.55, and removed in
Release 1.58.

Bio.GenBank.LocationParser
--------------------------
This module used to be used for parsing GenBank and EMBL feature locations.
It has been replaced with faster code using regular expressions, and is no
longer needed. Declared obsolete in Release 1.55, deprecated in Release 1.56,
and removed in Release 1.59.

Bio.Parsers and Bio.Parsers.spark
---------------------------------
This module was a copy of John Aycock's SPARK parser included with Biopython
solely for use in Bio.GenBank.LocationParser. Declared obsolete in Release
1.55, deprecated in Release 1.56, and removed in Release 1.59.

Bio.Restriction.DNAUtils and check_bases
----------------------------------------
This module (originally in C) offered complement and antiparallel functions
(duplicating functionality in Bio.Seq) and a rather odd function called
check_bases (also available as Bio.Restriction.Restriction.check_bases).
Deprecated in Release 1.53, removed in Release 1.57.

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

For some time now, both the NCBI and Biopython have encouraged people to
parse the XML output instead, however Bio.SearchIO will initially attempt
to support plain text BLAST output.

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
NCBI blastpgp was deprecated in Biopython 1.80. To parse tabular output
generated by BLAST programs, please use the ``parse`` function in
``Bio.Align``.

Bio.Clustalw
------------
Declared obsolete in Release 1.52, deprecated in Release 1.55 final, and
removed in Release 1.58. Replaced with Bio.AlignIO for parsing and writing
clustal format alignments (since Release 1.46), and Bio.Align.Applications
for calling the ClustalW command line tool (since Release 1.51). See the
Tutorial for examples.

BioSQL and psycopg
------------------
Support for psycopg (version one) in Biopython's BioSQL code was deprecated
in Release 1.51, and removed in Release 1.55. Please use psycopg2 instead.

BioSQL.BioSeqDatabase
---------------------
The ``remove_database`` and ``get_all_primary_ids`` methods were removed from
the ``DBServer`` class in Release 1.79.
The ``get_Seq_by_primary_id`` method was removed from the ``BioSeqDatabase``
class in Release 1.79.

Bio.Application.generic_run and ApplicationResult
-------------------------------------------------
Declared obsolete in Release 1.51, deprecated in Release 1.53, and removed in
Release 1.57. Please use the Python subprocess module instead, or as of
Release 1.55 the application wrappers can be used directly to execute the
command.

Bio.Entrez.efetch and rettype="genbank"
---------------------------------------
As of Easter 2009, the NCBI have stopped supporting the unofficial return type
of "genbank" in EFetch.  Instead we should be using "gb" (GenBank) or "gp"
(GenPept).  As of Biopython 1.50, Bio.Entrez.efetch will correct this
automatically, but issues a deprecation warning. The code to check and correct
for "genbank" was removed in Biopython 1.55 final.

Bio.SwissProt.SProt
-------------------
Declared obsolete in Release 1.50, deprecated in Release 1.51, and removed in
Release 1.56. Most of the functionality in Bio.SwissProt.SProt is available
from Bio.SwissProt.

Bio.Prosite and Bio.Enzyme
--------------------------
Declared obsolete in Release 1.50, deprecated in Release 1.53, and removed in
Release 1.57. Most of the functionality has moved to Bio.ExPASy.Prosite and
Bio.ExPASy.Enzyme, respectively.

Bio.EZRetrieve, Bio.NetCatch, Bio.FilteredReader
------------------------------------------------
Declared obsolete in Release 1.50, deprecated in Release 1.52, and removed in
Release 1.56.

Bio.File
--------
Bio.File.SGMLHandle was declared obsolete in Release 1.50, deprecated in
Release 1.52, and removed in Release 1.56. Bio.File.SGMLStripper was deprecated
in Release 1.57, removed in Release 1.61. Bio.File.StringHandle was deprecated
in Release 1.59, removed in Release 1.61.

Bio.Graphics.GenomeDiagram and colour/color, centre/center
----------------------------------------------------------
GenomeDiagram originally used colour and centre (UK spelling of color and
center) for argument names.  As part of its integration into Biopython 1.50,
this will support both colour and color, and both centre and center, to help
people port existing scripts written for the standalone version of
GenomeDiagram.  However, these were deprecated in Release 1.55 final.
Support for centre was removed in Release 1.62, and we intend to eventually
remove support for colour in later releases of Biopython.

Bio.AlignAce and Bio.MEME
-------------------------
Declared obsolete in Release 1.50, deprecated in Release 1.52, and removed
in Release 1.56. Please use Bio.Motif instead.

Bio.Seq, Bio.MutableSeq and the data property
---------------------------------------------
Direct use of the Seq object (and MutableSeq object) .data property is
deprecated.  As of Release 1.49, writing to the Seq object's .data property
triggered a warning, and this property was made read only in Release 1.53. In
Release 1.55 final, accessing the .data property of a Seq object gives a
DeprecationWarning. The Seq object's .data property was removed in Release
1.61.  Starting from Release 1.78, accessing the .data property of a MutableSeq
object similarly gives a deprecation warning.

Bio.Transcribe and Bio.Translate
--------------------------------
Declared obsolete in Release 1.49, deprecated in Release 1.51, and removed
in Release 1.57. Please use the methods or functions in Bio.Seq instead.

Bio.mathfns, Bio.stringfns and Bio.listfns (and their C code variants)
----------------------------------------------------------------------
Declared obsolete in Release 1.49. Bio.mathfns and Bio.stringfns were
deprecated in Release 1.50, Bio.listfns was deprecated in Release 1.53.
The three C implementations were all removed in Release 1.53. Bio.mathfns
and Bio.stringfns were removed in Release 1.55. Bio.listfns was removed in
Release 1.57.

Bio.Fasta (including Bio.Fasta.FastaAlign)
------------------------------------------
Declared obsolete in Release 1.48, deprecated in Release 1.51, and removed
in Release 1.55 final. Please use the "fasta" support in Bio.SeqIO or
Bio.AlignIO instead.

Note that ``Bio.Fasta`` could be used with a ``RecordParser`` which gave
``FastaRecord`` objects, for example::

    # Old code which won't work	any more
    from Bio import Fasta
    handle = open("example.fas")
    for record in Fasta.Iterator(handle, Fasta.RecordParser()) :
        # Here record was a Bio.Fasta.Record object
        print record.title # The full title line as a string
        print record.sequence # The sequence as a string
    handle.close()

Alternatively using the old ``SequenceParser`` would give ``SeqRecord``
objects like those from the new ``Bio.SeqIO`` code, for example::

    # Old code which won't work any more
    from Bio import Fasta
    handle = open("example.fas")
    for seq_record in Fasta.Iterator(handle, Fasta.SequenceParser()) :
        print seq_record.description # The full title line as a string
        print str(seq_record.seq) # The sequence as a string
    handle.close()

Either of those examples using ``Bio.SeqIO`` becomes just::

    # Updated versions of above examples using Bio.SeqIO instead
    from Bio import SeqIO
    for seq_record in SeqIO.parse("example.fas", "fasta") :
        print seq_record.description # The full title line as a string
        print str(seq_record.seq) # The sequence as a string

You can also continue to use handles with ``Bio.SeqIO`` if you want to.

Bio.Align.FormatConvert
-----------------------
Declared obsolete in Release 1.48, deprecated in Release 1.51, and
removed in Release 1.55 final. Instead, please use Bio.AlignIO or call the
format built-in function on the Alignment object.

Bio.Emboss.Applications
-----------------------
The wrappers for the "old" EMBOSS PHYLIP tools (e.g. eneighbor) were declared
obsolete in Biopython 1.52, deprecated in Release 1.55 final, and removed in
release 1.58. please use the wrappers for the "new" EMBOSS PHYLIP tools (e.g.
fneighbor) instead. Specifically, EProtDistCommandline, ENeighborCommandline,
EProtParsCommandline, EConsenseCommandline, and ESeqBootCommandline are
replaced by FProtDistCommandline, FNeighborCommandline, FProtParsCommandline,
FConsenseCommandline, and FSeqBootCommandline, respectively.

Bio.SeqIO.to_alignment()
------------------------
This function was made obsolete with the introduction of Bio.AlignIO,
deprecated in Release 1.54, and removed in Release 1.58. Use either the
Bio.AlignIO functions, or the Bio.Align.MultipleSeqAlignment class
directly instead.

Bio.SeqUtils
------------
Function makeTableX and classes ProteinX and MissingTable were deprecated
in Release 1.54, and removed in Release 1.58. These were remnants of the
removed translate function, and no longer served any useful purpose.

Function 'reverse' in Bio.SeqUtils was deprecated in Release 1.54, and
removed in Release 1.58. Instead just use the string's slice method with
a step of minus one.

Functions GC_Frame, fasta_uniqids, apply_on_multi_fasta, and
quicker_apply_on_multi_fasta were deprecated in Release 1.55, and removed
in Release 1.58.

Function quick_FASTA_reader was declared obsolete in Release 1.61,
deprecated in Release 1.64, and removed in Release 1.67. Use function
list(SimpleFastaParser(handle)) from Bio.SeqIO.FastaIO instead (but
ideally convert your code to using an iterator approach).

The 'title2ids' argument to FastaIterator in Bio.SeqIO.FastaIO and
FastqPhredIterator in Bio.SeqIO.QualityIO was deprecated in Release 1.80.
Please use a generator function to modify the records returned by the parser.

Function Tm_staluc in Bio.SeqUtils.MeltingTemp was deprecated in Release 1.78,
and removed in Release 1.80.

The method 'print_index' of the CodonAdaptationIndex class in
Bio.SeqUtils.CodonUsage was deprecated in Release 1.80. Instead of
self.print_index(), please use print(self).

The modules Bio.SeqUtils.CodonUsage and Bio.SeqUtils.CodonUsageIndices were
deprecated in Release 1.80. Please use the new CodonAdaptationIndex class in
Bio.SeqUtils instead. Note that this class has been updated to use modern
Python, and may give slightly different results from the CodonAdaptationIndex
class in Bio.SeqUtils.CodonUsage, as the calculation was updated to be
consistent with the calculated values by Sharp & Li.

Function 'GC' in Bio.SeqUtils was deprecated in Release 1.80. Instead use
function 'gc_fraction'.

Bio.GFF (for accessing a MySQL database created with BioPerl, etc)
------------------------------------------------------------------
The whole of the old ``Bio.GFF`` module was deprecated in Release 1.53, and
removed in Release 1.57 (with the intention of reusing this name space for a
GFF parser).

Bio.utils
---------
Function 'ungap' was deprecated in Release 1.53. Use Bio.Seq instead.
The whole of Bio.utils was declared obsolete in Release 1.55, deprecated in
Release 1.56, and removed in Release 1.57.

Bio.Compass
-----------
The RecordParser and Iterator classes were declared obsolete in Release 1.54,
deprecated in Release 1.55, removed in Release 1.59. Their functionality is
now available through a read() and a parse() function, respectively.

Bio.Affy.CelFile
----------------
The CelScanner, CelConsumer, CelRecord, and CelParser were declared obsolete
in Release 1.54, deprecated in Release 1.55 and removed in Release 1.59.
Their functionality is now available through a read() function.

Bio.PopGen.Async
----------------
``Bio.PopGen.Async`` was deprecated in Release 1.68, removed in Release 1.70.

Bio.PopGen.FDist
----------------
``Bio.PopGen.FDist`` was deprecated in Release 1.68, removed in Release 1.70.

Prior to this, the ``RecordParser``, ``_Scanner``, and ``_RecordConsumer``
classes were declared obsolete in Release 1.54, deprecated in Release 1.55,
and removed in Release 1.58. Their functionality is now available through
a ``read()`` function.

Bio.PopGen.SimCoal
------------------
``Bio.PopGen.SimCoal`` was deprecated in Release 1.68, and removed in Release
1.70.

Bio.UniGene
-----------
The classes UnigeneSequenceRecord, UnigeneProtsimRecord, UnigeneSTSRecord,
UnigeneRecord, _RecordConsumer, _Scanner, RecordParser, and Iterator in
Bio.UniGene were declared obsolete in Release 1.54, deprecated in Release 1.55,
and removed in Release 1.59. Their functionality is now available through a
read() and a parse() function in Bio.UniGene.

Submodule Bio.UniGene.UniGene which was an HTML parser was declared obsolete
in Release 1.59, deprecated in Release 1.61, and removed in Release 1.64.

Bio.SubsMat
-----------
The methods ``letter_sum`` and ``all_letters_sum`` were removed from the
``SeqMat`` class in Bio.SubsMat in Release 1.57.
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

The PairwiseAlignment class was deprecated in Release 1.80; please use the new
Alignment class instead.

Bio.Align.Generic
-----------------
This module which defined to original (Multiple-Sequence) Alignment class was
deprecated in Release 1.57 and removed in Release 1.69.

Bio.ParserSupport
-----------------
``Bio.ParserSupport`` was declared obsolete in Release 1.59, and deprecated in
Release 1.63. The Martel specific ``EventGenerator`` was removed in Release
1.67, and the entire module was removed in Release 1.72.

``Bio.ParserSupport.SGMLStrippingConsumer`` was deprecated in Release 1.59, and
removed in Release 1.61.

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
The ``Bio.Wise`` module was deprecated in Release 1.80.

Scripts/Restriction/ranacompiler.py
-----------------------------------
The ``is_palindrom`` function was removed in Release 1.79.
