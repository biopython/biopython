This file provides documentation for modules in Biopython that have been moved
or deprecated in favor of other modules. This provides some quick and easy
to find documentation about how to update your code to work again.
Python releases go first, then code (modules, methods, functions).

Python
======

Python 2.3
----------
No longer supported as of Release 1.51, having triggered a warning with
Release 1.50, with a mention in the release notes for Release 1.49.

Python 2.4
----------
No longer supported as of Release 1.59, having triggered a warning since
Release 1.55, with advance notice in the release notes for Release 1.54.

Python 2.5
----------
No longer supported as of Release 1.63, having triggered a warning with
Release 1.62, with advance notice in the release notes for Release 1.61.

Python 2.6
----------
No longer supported as of Release 1.69, having triggered a warning with
release 1.66 onwards.

Python 2.7
----------
No longer supported as of Release 1.77 (2020, in line with end-of-life or
sunset date for Python 2.7 itself), having triggered a warning in prior
releases.

Python 3.0, 3.1, 3.2
--------------------
Never officially supported, these triggered a warning in Release 1.62
recommending Python 3.3 or later. As of Biopython Release 1.63 onwards,
installation simply aborts with a error message.

Python 3.3
----------
No longer supported as of Release 1.70, having triggered a warning with
release 1.67 onwards.

Python 3.4
----------
No longer supported as of Release 1.75, having triggered a deprecation
warning in release 1.74. First supported in release 1.64.

Python 3.5
----------
No longer supported as of Release 1.77. First supported in release 1.66.

Python 3.6
----------
First supported in release 1.69.

Python 3.7
----------
First supported in release 1.73.

Python 3.8
----------
First supported in release 1.75.

Jython
------
No longer supported as of Release 1.77 with the end of Python 2 support.
Biopython was mostly working under Jython 2.7.0, but support for Jython
was deprecated as of Release 1.70.

Biopython modules, methods, functions
=====================================

Bio.Index
---------
Deprecated in release 1.75, removed in release 1.77. Was not used anywhere in
Biopython.

Bio.Crystal
-----------
Declared obsolete in release 1.75, deprecated in release 1.76. PDB NDB files
can be opened with Bio.PDB.

Bio.motifs
----------
``Bio.motifs.mast`` plain-text parsing deprecated in favor of XML parsing as of
release 1.74. Also affects ``Bio.motifs.read`` and ``Bio.motifs.parse`` for the
``mast`` format.
The ``format`` method of the ``Motif`` class in ``Bio.motifs`` was deprecated
in release 1.77, in favor of a ``__format__`` method that can be used from the
``format`` built-in function.

Bio.Restriction.RanaConfig
--------------------------
Removed in Biopython 1.74 without explicit depreciation period. RanaConfig was
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
use string-like equality testing and hashing (ingoring any difference in
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

Bio.SeqFeature
--------------
With the introduction of the CompoundLocation in Release 1.62, the SeqFeature
attribute sub_features was deprecated. It was removed in Release 1.68.

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

Class SequentialSequenceWriter was declared obsolete in Release 1.77, and
deprecated in Release 1.78.

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

Earlier the function ``pairlist_to_dict`` was deprecated in Release 1.45, and
removed in Release 1.53.

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

For some time now, both the NCBI and Biopython have encouraged people to
parse the XML output instead, however Bio.SearchIO will initially attempt
to support plain text BLAST output.

The module was removed in Release 1.72 from the public API. It lives now
in maintenance mode in Bio.SearchIO._legacy to preserve existing functionality.

Bio.Blast.Applications
----------------------
NCBI "legacy" BLAST tool wrappers FastacmdCommandline, BlastallCommandline,
BlastpgpCommandline and RpsBlastCommandline were declared obsolete in Release
1.53, deprecated in Release 1.61, and removed in Release 1.64, having been
replaced with wrappers for the new NCBI BLAST+ tools (e.g.
NcbiblastpCommandline and NcbipsiblastCommandline).

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

Bio.Entrez.query function
-------------------------
Deprecated in Release 1.47, removed in Release 1.52.

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

Numeric support
---------------
Following the Release of 1.48, Numeric support in Biopython is discontinued.
Please move to NumPy for Biopython 1.49 or later.

Bio.Seq and the data property
-----------------------------
Direct use of the Seq object (and MutableSeq object) .data property is
deprecated.  As of Release 1.49, writing to the Seq object's .data property
triggered a warning, and this property was made read only in Release 1.53. In
Release 1.55 final, accessing the .data property gives a DeprecationWarning.
The Seq object's .data property was removed in Release 1.61.

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

Bio.distance (and Bio.cdistance)
--------------------------------
Bio.distance was deprecated in Release 1.49, at which point its C code
implementation Bio.cdistance was removed (this was not intended as a public
API). Removed in Release 1.53.

Bio.Ndb
-------
Deprecated in Release 1.49, as the website this parsed has been redesigned.
Removed in Release 1.53.

Martel
------
Declared obsolete in Release 1.48, deprecated in Release 1.49, and removed
in Release 1.51.  The source code for Martel is still in our repository if
anyone wanted to develop this outside of Biopython.

Bio.Mindy and associated modules.
---------------------------------
Declared obsolete in Release 1.48, deprecated in Release 1.49, removed in
Release 1.51.  This includes the Bio.Writer, Bio.writers, Bio.builders,
Bio.Std, Bio.StdHandler, Bio.Decode and Bio.DBXRef modules

Bio.Fasta index_file and Dictionary
-----------------------------------
Deprecated in Release 1.44, removed in Biopython 1.46. For small to medium
sized files, use Bio.SeqIO.to_dict() to make an in memory dictionary of
SeqRecord objects. Biopython 1.52 onwards provides Bio.SeqIO.index()
which is suitable even for very large files.

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

Bio.Emboss.Primer
-----------------
Deprecated in Release 1.48, and removed in Release 1.51, this parser was
replaced by Bio.Emboss.Primer3 and Bio.Emboss.PrimerSearch instead.

Bio.Emboss.Applications
-----------------------
The wrappers for the "old" EMBOSS PHYLIP tools (e.g. eneighbor) were declared
obsolete in Biopython 1.52, deprecated in Release 1.55 final, and removed in
release 1.58. please use the wrappers for the "new" EMBOSS PHYLIP tools (e.g.
fneighbor) instead. Specifically, EProtDistCommandline, ENeighborCommandline,
EProtParsCommandline, EConsenseCommandline, and ESeqBootCommandline are
replaced by FProtDistCommandline, FNeighborCommandline, FProtParsCommandline,
FConsenseCommandline, and FSeqBootCommandline, respectively.

Bio.MetaTool
------------
Deprecated in Release 1.48, and removed in Release 1.51, this was a parser
for the output of MetaTool 3.5 which is now obsolete.

Bio.GenBank
-----------
The online functionality (search_for, download_many, and NCBIDictionary) was
declared obsolete in Release 1.48, deprecated in Release 1.50, and removed
in Release 1.54. Please use Bio.Entrez instead.

Bio.PubMed
----------
Declared obsolete in Release 1.48, deprecated in Release 1.49, and
removed in Release 1.53. Please use Bio.Entrez instead.

Bio.EUtils
----------
Deprecated in favor of Bio.Entrez in Release 1.48, removed in Release 1.52.

Bio.Sequencing & Bio.Medline
----------------------------
A revised API was added and the old one deprecated in Release 1.48,
and removed in Biopython 1.52:

* Bio.Sequencing.Ace.RecordParser --> Bio.Sequencing.Ace.read(handle)
* Bio.Sequencing.Ace.Iterator --> Bio.Sequencing.Ace.parse(handle)
* Bio.Sequencing.Phd.RecordParser --> Bio.Sequencing.Phd.read(handle)
* Bio.Sequencing.Phd.Iterator --> Bio.Sequencing.Phd.parse(handle)
* Bio.Medline.RecordParser --> Bio.Medline.read(handle)
* Bio.Medline.Iterator --> Bio.Medline.parse(handle)

Bio.Blast.NCBIWWW
-----------------
The HTML BLAST parser was deprecated in Release 1.48, and removed in 1.52.
The deprecated functions blast and blasturl were removed in Release 1.44.

Bio.Saf
-------
Deprecated as of Release 1.48, removed in Release 1.51.  If useful, a parser
for this "simple alignment format" could be developed for Bio.AlignIO instead.

Bio.NBRF
--------
Deprecated as of Release 1.48 in favor of the "pir" format in Bio.SeqIO,
removed in Release 1.51.

Bio.IntelliGenetics
-------------------
Deprecated as of Release 1.48 in favor of the "ig" format in Bio.SeqIO,
removed in Release 1.51.

Bio.SeqIO submodules PhylipIO, ClustalIO, NexusIO and StockholmIO
-----------------------------------------------------------------
You can still use the "phylip", "clustal", "nexus" and "stockholm" formats
in Bio.SeqIO, however these are now supported via Bio.AlignIO, with the
old code deprecated in Releases 1.46 or 1.47, and removed in Release 1.49.

Bio.SeqIO.to_alignment()
------------------------
This function was made obsolete with the introduction of Bio.AlignIO,
deprecated in Release 1.54, and removed in Release 1.58. Use either the
Bio.AlignIO functions, or the Bio.Align.MultipleSeqAlignment class
directly instead.

Bio.ECell
---------
Deprecated as of Release 1.47, as it appears to have no users, and the code
does not seem relevant for ECell 3.  Removed in Release 1.49.

Bio.Ais
-------
Deprecated as of Release 1.45, removed in Release 1.49.

Bio.LocusLink
-------------
Deprecated as of Release 1.45, removed in Release 1.49.
The NCBI's LocusLink was superseded by Entrez Gene.

Bio.SGMLExtractor
-----------------
Deprecated as of Release 1.46, removed in Release 1.49.

Bio.Rebase
----------
Deprecated as of Release 1.46, removed in Release 1.49.

Bio.Gobase
----------
Deprecated as of Release 1.46, removed in Release 1.49.

Bio.CDD
-------
Deprecated as of Release 1.46, removed in Release 1.49.

Bio.biblio
----------
Deprecated as of Release 1.45, removed in Release 1.48

Bio.WWW
-------
The modules under Bio.WWW were deprecated in Release 1.45, and removed in
Release 1.48.  The remaining stub Bio.WWW was deprecated in Release 1.48,
and removed in Release 1.53.

The functionality in Bio.WWW.SCOP, Bio.WWW.InterPro, Bio.WWW.ExPASy and
Bio.WWW.NCBI is now available from Bio.SCOP, Bio.InterPro, Bio.ExPASy and
Bio.Entrez instead.

Bio.SeqIO
---------
The old Bio.SeqIO.FASTA and Bio.SeqIO.generic were deprecated in favour of
the new Bio.SeqIO module as of Release 1.44, removed in Release 1.47.

Bio.Medline.NLMMedlineXML
-------------------------
Deprecated in Release 1.44, removed in 1.46.

Bio.MultiProc
-------------
Deprecated in Release 1.44, removed in 1.46.

Bio.MarkupEditor
----------------
Deprecated in Release 1.44, removed in 1.46.

Bio.lcc
-------
Deprecated in favor of Bio.SeqUtils.lcc in Release 1.44, removed in 1.46.

Bio.crc
-------
Deprecated in favor of Bio.SeqUtils.CheckSum in Release 1.44, removed in 1.46.

Bio.FormatIO
------------
This was removed in Release 1.44 (a deprecation was not possible).

Bio.expressions, Bio.config, Bio.dbdefs, Bio.formatdefs and Bio.dbdefs
----------------------------------------------------------------------
These were deprecated in Release 1.44, and removed in Release 1.49.

Bio.Kabat
---------
This was deprecated in Release 1.43 and removed in Release 1.44.

Bio.SeqUtils
------------
Functions 'complement' and 'antiparallel' in Bio.SeqUtils were deprecated
in Release 1.31, and removed in Release 1.43.  Function 'translate' was
deprecated in Release 1.49, and removed in Release 1.53. Use the functions
and methods in Bio.Seq instead.

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

Bio.GFF (for accessing a MySQL database created with BioPerl, etc)
------------------------------------------------------------------
The functions ``forward_complement`` and ``antiparallel`` in ``Bio.GFF.easy``
have been deprecated as of Release 1.31, and removed in Release 1.43.
Use the functions ``complement`` and ``reverse_complement`` in ``Bio.Seq``
instead.

The whole of the old ``Bio.GFF`` module was deprecated in Release 1.53, and
removed in Release 1.57 (with the intention of reusing this name space for a
GFF parser).

Bio.sequtils
------------
Deprecated as of Release 1.30, removed in Release 1.42. Use ``Bio.SeqUtils``
instead.

Bio.SVM
-------
Deprecated as of Release 1.30, removed in Release 1.42.
The Support Vector Machine code in Biopython has been superseded by a
more robust (and maintained) SVM library, which includes a python
interface. We recommend using LIBSVM:

http://www.csie.ntu.edu.tw/~cjlin/libsvm/

Bio.RecordFile
--------------
Deprecated as of Release 1.30, removed in Release 1.42.  RecordFile wasn't
completely implemented and duplicates the work of most standard parsers.

Bio.kMeans and Bio.xkMeans
--------------------------
Deprecated as of Release 1.30, removed in Release 1.42.  Instead, please use
the function kcluster in Bio.Cluster which performs k-means or k-medians
clustering.

Bio.SCOP
--------
The module Bio.SCOP.FileIndex was deprecated in Release 1.46, and removed in
Release 1.53. The class Parser in Bio.SCOP.Dom was removed in Release 1.55
final. The class Iterator in Bio.SCOP.Dom was removed in Release 1.56.

Dictionary to_one_letter_code in module Bio.SCOP.three_to_one_dict was moved
to protein_letters_3to1 in module Bio.Data.SCOPData in Release 1.62. The old
alias was preserved with a deprecation warning, until it was removed in
Release 1.66.

Bio.utils
---------
Functions 'translate', 'translate_to_stop', 'back_translate', 'transcribe',
and 'back_transcribe' were deprecated in Release 1.49, and removed in Release
1.53. Function 'ungap' was deprecated in Release 1.53. Use Bio.Seq instead.
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
The methods letter_sum and all_letters_sum were removed from the SeqMat class
in Bio.SubsMat in Release 1.57.
The Bio.SubsMat module was deprecated in Release 1.78. As an alternative,
please consider using Bio.Align.substitution_matrices.

Bio.Align
---------
The methods get_column and add_sequence of the MultipleSeqAlignment class were
deprecated in Release 1.57 and removed in Release 1.69.
The format method of the MultipleSeqAlignment class and the PairwiseAlignment
class were deprecated in Release 1.76.

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
This module was declared obsolete in Release 1.74, and deprecated in Release
1.76.

Bio.File
--------
The UndoHandle class was deprecated in Release 1.77, and moved to
Bio/SearchIO/_legacy/ParserSupport.py, which was the only module in
Biopython still using this class.

Bio.FSSP
-----------
Deprecated in release 1.77.
