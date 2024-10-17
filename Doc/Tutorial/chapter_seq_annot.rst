.. _`chapter:seq_annot`:

Sequence annotation objects
===========================

Chapter :ref:`chapter:seq_objects` introduced the
sequence classes. Immediately “above” the ``Seq`` class is the Sequence
Record or ``SeqRecord`` class, defined in the ``Bio.SeqRecord`` module.
This class allows higher level features such as identifiers and features
(as ``SeqFeature`` objects) to be associated with the sequence, and is
used throughout the sequence input/output interface ``Bio.SeqIO``
described fully in Chapter :ref:`chapter:seqio`.

If you are only going to be working with simple data like FASTA files,
you can probably skip this chapter for now. If on the other hand you are
going to be using richly annotated sequence data, say from GenBank or
EMBL files, this information is quite important.

While this chapter should cover most things to do with the ``SeqRecord``
and ``SeqFeature`` objects in this chapter, you may also want to read
the ``SeqRecord`` wiki page (http://biopython.org/wiki/SeqRecord), and
the built-in documentation (:py:mod:`Bio.SeqRecord` and
:py:mod:`Bio.SeqFeature`):

.. code:: pycon

   >>> from Bio.SeqRecord import SeqRecord
   >>> help(SeqRecord)

.. _`sec:SeqRecord`:

The SeqRecord object
--------------------

The ``SeqRecord`` (Sequence Record) class is defined in the
``Bio.SeqRecord`` module. This class allows higher level features such
as identifiers and features to be associated with a sequence (see
Chapter :ref:`chapter:seq_objects`), and is the
basic data type for the ``Bio.SeqIO`` sequence input/output interface
(see Chapter :ref:`chapter:seqio`).

The ``SeqRecord`` class itself is quite simple, and offers the following
information as attributes:

.seq
   The sequence itself, typically a ``Seq`` object.

.id
   The primary ID used to identify the sequence – a string. In most
   cases this is something like an accession number.

.name
   A “common” name/id for the sequence – a string. In some cases this
   will be the same as the accession number, but it could also be a
   clone name. I think of this as being analogous to the LOCUS id in a
   GenBank record.

.description
   A human readable description or expressive name for the sequence –
   a string.

.letter_annotations
   Holds per-letter-annotations using a (restricted) dictionary of
   additional information about the letters in the sequence. The keys
   are the name of the information, and the information is contained in
   the value as a Python sequence (i.e. a list, tuple or string) with
   the same length as the sequence itself. This is often used for
   quality scores (e.g.
   Section :ref:`sec:FASTQ-filtering-example`)
   or secondary structure information (e.g. from Stockholm/PFAM
   alignment files).

.annotations
   A dictionary of additional information about the sequence. The keys
   are the name of the information, and the information is contained in
   the value. This allows the addition of more “unstructured”
   information to the sequence.

.features
   A list of ``SeqFeature`` objects with more structured information
   about the features on a sequence (e.g. position of genes on a genome,
   or domains on a protein sequence). The structure of sequence features
   is described below in Section :ref:`sec:seq_features`.

.dbxrefs
   A list of database cross-references as strings.

Creating a SeqRecord
--------------------

Using a ``SeqRecord`` object is not very complicated, since all of the
information is presented as attributes of the class. Usually you won’t
create a ``SeqRecord`` “by hand”, but instead use ``Bio.SeqIO`` to read
in a sequence file for you (see
Chapter :ref:`chapter:seqio` and the examples below).
However, creating ``SeqRecord`` can be quite simple.

SeqRecord objects from scratch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To create a ``SeqRecord`` at a minimum you just need a ``Seq`` object:

.. doctest

.. code:: pycon

   >>> from Bio.Seq import Seq
   >>> simple_seq = Seq("GATC")
   >>> from Bio.SeqRecord import SeqRecord
   >>> simple_seq_r = SeqRecord(simple_seq)

Additionally, you can also pass the id, name and description to the
initialization function, but if not they will be set as strings
indicating they are unknown, and can be modified subsequently:

.. cont-doctest

.. code:: pycon

   >>> simple_seq_r.id
   '<unknown id>'
   >>> simple_seq_r.id = "AC12345"
   >>> simple_seq_r.description = "Made up sequence I wish I could write a paper about"
   >>> print(simple_seq_r.description)
   Made up sequence I wish I could write a paper about
   >>> simple_seq_r.seq
   Seq('GATC')

Including an identifier is very important if you want to output your
``SeqRecord`` to a file. You would normally include this when creating
the object:

.. doctest

.. code:: pycon

   >>> from Bio.Seq import Seq
   >>> simple_seq = Seq("GATC")
   >>> from Bio.SeqRecord import SeqRecord
   >>> simple_seq_r = SeqRecord(simple_seq, id="AC12345")

As mentioned above, the ``SeqRecord`` has an dictionary attribute
``annotations``. This is used for any miscellaneous annotations that
doesn’t fit under one of the other more specific attributes. Adding
annotations is easy, and just involves dealing directly with the
annotation dictionary:

.. cont-doctest

.. code:: pycon

   >>> simple_seq_r.annotations["evidence"] = "None. I just made it up."
   >>> print(simple_seq_r.annotations)
   {'evidence': 'None. I just made it up.'}
   >>> print(simple_seq_r.annotations["evidence"])
   None. I just made it up.

Working with per-letter-annotations is similar, ``letter_annotations``
is a dictionary like attribute which will let you assign any Python
sequence (i.e. a string, list or tuple) which has the same length as the
sequence:

.. cont-doctest

.. code:: pycon

   >>> simple_seq_r.letter_annotations["phred_quality"] = [40, 40, 38, 30]
   >>> print(simple_seq_r.letter_annotations)
   {'phred_quality': [40, 40, 38, 30]}
   >>> print(simple_seq_r.letter_annotations["phred_quality"])
   [40, 40, 38, 30]

The ``dbxrefs`` and ``features`` attributes are just Python lists, and
should be used to store strings and ``SeqFeature`` objects (discussed
later in this chapter) respectively.

SeqRecord objects from FASTA files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example uses a fairly large FASTA file containing the whole
sequence for *Yersinia pestis biovar Microtus* str. 91001 plasmid pPCP1,
originally downloaded from the NCBI. This file is included with the
Biopython unit tests under the GenBank folder, or online
`NC_005816.fna <https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.fna>`__
from our website.

The file starts like this - and you can check there is only one record
present (i.e. only one line starting with a greater than symbol):

.. code:: text

   >gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus ... pPCP1, complete sequence
   TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGGGGGTAATCTGCTCTCC
   ...

Back in Chapter :ref:`chapter:quick_start` you
will have seen the function ``Bio.SeqIO.parse(...)`` used to loop over
all the records in a file as ``SeqRecord`` objects. The ``Bio.SeqIO``
module has a sister function for use on files which contain just one
record which we’ll use here (see
Chapter :ref:`chapter:seqio` for details):

.. doctest ../Tests/GenBank

.. code:: pycon

   >>> from Bio import SeqIO
   >>> record = SeqIO.read("NC_005816.fna", "fasta")
   >>> record
   SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='gi|45478711|ref|NC_005816.1|', name='gi|45478711|ref|NC_005816.1|', description='gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])

Now, let’s have a look at the key attributes of this ``SeqRecord``
individually – starting with the ``seq`` attribute which gives you a
``Seq`` object:

.. cont-doctest

.. code:: pycon

   >>> record.seq
   Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')

Next, the identifiers and description:

.. cont-doctest

.. code:: pycon

   >>> record.id
   'gi|45478711|ref|NC_005816.1|'
   >>> record.name
   'gi|45478711|ref|NC_005816.1|'
   >>> record.description
   'gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'

As you can see above, the first word of the FASTA record’s title line
(after removing the greater than symbol) is used for both the ``id`` and
``name`` attributes. The whole title line (after removing the greater
than symbol) is used for the record description. This is deliberate,
partly for backwards compatibility reasons, but it also makes sense if
you have a FASTA file like this:

.. code:: text

   >Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1
   TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGGGGGTAATCTGCTCTCC
   ...

Note that none of the other annotation attributes get populated when
reading a FASTA file:

.. cont-doctest

.. code:: pycon

   >>> record.dbxrefs
   []
   >>> record.annotations
   {}
   >>> record.letter_annotations
   {}
   >>> record.features
   []

In this case our example FASTA file was from the NCBI, and they have a
fairly well defined set of conventions for formatting their FASTA lines.
This means it would be possible to parse this information and extract
the GI number and accession for example. However, FASTA files from other
sources vary, so this isn’t possible in general.

SeqRecord objects from GenBank files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As in the previous example, we’re going to look at the whole sequence
for *Yersinia pestis biovar Microtus* str. 91001 plasmid pPCP1,
originally downloaded from the NCBI, but this time as a GenBank file.
Again, this file is included with the Biopython unit tests under the
GenBank folder, or online
`NC_005816.gb <https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.gb>`__
from our website.

This file contains a single record (i.e. only one LOCUS line) and
starts:

.. code:: text

   LOCUS       NC_005816               9609 bp    DNA     circular BCT 21-JUL-2008
   DEFINITION  Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete
               sequence.
   ACCESSION   NC_005816
   VERSION     NC_005816.1  GI:45478711
   PROJECT     GenomeProject:10638
   ...

Again, we’ll use ``Bio.SeqIO`` to read this file in, and the code is
almost identical to that for used above for the FASTA file (see
Chapter :ref:`chapter:seqio` for details):

.. doctest ../Tests/GenBank

.. code:: pycon

   >>> from Bio import SeqIO
   >>> record = SeqIO.read("NC_005816.gb", "genbank")
   >>> record
   SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])

.. cont-doctest

.. code:: pycon

   >>> record.seq
   Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')

The ``name`` comes from the LOCUS line, while the ``id`` includes the
version suffix. The description comes from the DEFINITION line:

.. cont-doctest

.. code:: pycon

   >>> record.id
   'NC_005816.1'
   >>> record.name
   'NC_005816'
   >>> record.description
   'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'

GenBank files don’t have any per-letter annotations:

.. cont-doctest

.. code:: pycon

   >>> record.letter_annotations
   {}

Most of the annotations information gets recorded in the ``annotations``
dictionary, for example:

.. cont-doctest

.. code:: pycon

   >>> len(record.annotations)
   13
   >>> record.annotations["source"]
   'Yersinia pestis biovar Microtus str. 91001'

The ``dbxrefs`` list gets populated from any PROJECT or DBLINK lines:

.. cont-doctest

.. code:: pycon

   >>> record.dbxrefs
   ['Project:58037']

Finally, and perhaps most interestingly, all the entries in the features
table (e.g. the genes or CDS features) get recorded as ``SeqFeature``
objects in the ``features`` list.

.. cont-doctest

.. code:: pycon

   >>> len(record.features)
   41

We’ll talk about ``SeqFeature`` objects next, in
Section :ref:`sec:seq_features`.

.. _`sec:seq_features`:

Feature, location and position objects
--------------------------------------

SeqFeature objects
~~~~~~~~~~~~~~~~~~

Sequence features are an essential part of describing a sequence. Once
you get beyond the sequence itself, you need some way to organize and
easily get at the more “abstract” information that is known about the
sequence. While it is probably impossible to develop a general sequence
feature class that will cover everything, the Biopython ``SeqFeature``
class attempts to encapsulate as much of the information about the
sequence as possible. The design is heavily based on the GenBank/EMBL
feature tables, so if you understand how they look, you’ll probably have
an easier time grasping the structure of the Biopython classes.

The key idea about each ``SeqFeature`` object is to describe a region on
a parent sequence, typically a ``SeqRecord`` object. That region is
described with a location object, typically a range between two
positions (see Section :ref:`sec:locations` below).

The ``SeqFeature`` class has a number of attributes, so first we’ll list
them and their general features, and then later in the chapter work
through examples to show how this applies to a real life example. The
attributes of a SeqFeature are:

.type
   This is a textual description of the type of feature (for instance,
   this will be something like ‘CDS’ or ‘gene’).

.location
   The location of the ``SeqFeature`` on the sequence that you are
   dealing with, see Section :ref:`sec:locations` below. The
   ``SeqFeature`` delegates much of its functionality to the location
   object, and includes a number of shortcut attributes for properties
   of the location:

   .ref
      shorthand for ``.location.ref`` – any (different) reference
      sequence the location is referring to. Usually just None.

   .ref_db
      shorthand for ``.location.ref_db`` – specifies the database any
      identifier in ``.ref`` refers to. Usually just None.

   .strand
      shorthand for ``.location.strand`` – the strand on the sequence
      that the feature is located on. For double stranded nucleotide
      sequence this may either be :math:`1` for the top strand,
      :math:`-1` for the bottom strand, :math:`0` if the strand is
      important but is unknown, or ``None`` if it doesn’t matter. This
      is None for proteins, or single stranded sequences.

.qualifiers
   This is a Python dictionary of additional information about the
   feature. The key is some kind of terse one-word description of what
   the information contained in the value is about, and the value is the
   actual information. For example, a common key for a qualifier might
   be “evidence” and the value might be “computational
   (non-experimental).” This is just a way to let the person who is
   looking at the feature know that it has not be experimentally
   (i. e. in a wet lab) confirmed. Note that other the value will be a
   list of strings (even when there is only one string). This is a
   reflection of the feature tables in GenBank/EMBL files.

.sub_features
   This used to be used to represent features with complicated
   locations like ‘joins’ in GenBank/EMBL files. This has been
   deprecated with the introduction of the ``CompoundLocation`` object,
   and should now be ignored.

.. _`sec:locations`:

Positions and locations
~~~~~~~~~~~~~~~~~~~~~~~

The key idea about each ``SeqFeature`` object is to describe a region on
a parent sequence, for which we use a location object, typically
describing a range between two positions. Two try to clarify the
terminology we’re using:

position
   This refers to a single position on a sequence, which may be fuzzy
   or not. For instance, 5, 20, ``<100`` and ``>200`` are all positions.

location
   A location is region of sequence bounded by some positions. For
   instance ``5..20`` (i. e. 5 to 20) is a location.

I just mention this because sometimes I get confused between the two.

SimpleLocation object
^^^^^^^^^^^^^^^^^^^^^

Unless you work with eukaryotic genes, most ``SeqFeature`` locations are
extremely simple - you just need start and end coordinates and a strand.
That’s essentially all the basic ``SimpleLocation`` object does.

In practice of course, things can be more complicated. First of all we
have to handle compound locations made up of several regions. Secondly,
the positions themselves may be fuzzy (inexact).

CompoundLocation object
^^^^^^^^^^^^^^^^^^^^^^^

Biopython 1.62 introduced the ``CompoundLocation`` as part of a
restructuring of how complex locations made up of multiple regions are
represented. The main usage is for handling ‘join’ locations in
EMBL/GenBank files.

Fuzzy Positions
^^^^^^^^^^^^^^^

So far we’ve only used simple positions. One complication in dealing
with feature locations comes in the positions themselves. In biology
many times things aren’t entirely certain (as much as us wet lab
biologists try to make them certain!). For instance, you might do a
dinucleotide priming experiment and discover that the start of mRNA
transcript starts at one of two sites. This is very useful information,
but the complication comes in how to represent this as a position. To
help us deal with this, we have the concept of fuzzy positions.
Basically there are several types of fuzzy positions, so we have five
classes to deal with them:

ExactPosition
   As its name suggests, this class represents a position which is
   specified as exact along the sequence. This is represented as just a
   number, and you can get the position by looking at the ``position``
   attribute of the object.

BeforePosition
   This class represents a fuzzy position that occurs prior to some
   specified site. In GenBank/EMBL notation, this is represented as
   something like ``<13``, signifying that the real position
   is located somewhere less than 13. To get the specified upper
   boundary, look at the ``position`` attribute of the object.

AfterPosition
   Contrary to ``BeforePosition``, this class represents a position
   that occurs after some specified site. This is represented in GenBank
   as ``>13``, and like ``BeforePosition``, you get the
   boundary number by looking at the ``position`` attribute of the
   object.

WithinPosition
   Occasionally used for GenBank/EMBL locations, this class models a
   position which occurs somewhere between two specified nucleotides. In
   GenBank/EMBL notation, this would be represented as ``(1.5)``, to
   represent that the position is somewhere within the range 1 to 5.

OneOfPosition
   Occasionally used for GenBank/EMBL locations, this class deals with
   a position where several possible values exist, for instance you
   could use this if the start codon was unclear and there where two
   candidates for the start of the gene. Alternatively, that might be
   handled explicitly as two related gene features.

UnknownPosition
   This class deals with a position of unknown location. This is not
   used in GenBank/EMBL, but corresponds to the ‘?’ feature coordinate
   used in UniProt.

Here’s an example where we create a location with fuzzy end points:

.. doctest

.. code:: pycon

   >>> from Bio import SeqFeature
   >>> start_pos = SeqFeature.AfterPosition(5)
   >>> end_pos = SeqFeature.BetweenPosition(9, left=8, right=9)
   >>> my_location = SeqFeature.SimpleLocation(start_pos, end_pos)

Note that the details of some of the fuzzy-locations changed in
Biopython 1.59, in particular for BetweenPosition and WithinPosition you
must now make it explicit which integer position should be used for
slicing etc. For a start position this is generally the lower (left)
value, while for an end position this would generally be the higher
(right) value.

If you print out a ``SimpleLocation`` object, you can get a nice
representation of the information:

.. cont-doctest

.. code:: pycon

   >>> print(my_location)
   [>5:(8^9)]

We can access the fuzzy start and end positions using the start and end
attributes of the location:

.. cont-doctest

.. code:: pycon

   >>> my_location.start
   AfterPosition(5)
   >>> print(my_location.start)
   >5
   >>> my_location.end
   BetweenPosition(9, left=8, right=9)
   >>> print(my_location.end)
   (8^9)

If you don’t want to deal with fuzzy positions and just want numbers,
they are actually subclasses of integers so should work like integers:

.. cont-doctest

.. code:: pycon

   >>> int(my_location.start)
   5
   >>> int(my_location.end)
   9

Similarly, to make it easy to create a position without worrying about
fuzzy positions, you can just pass in numbers to the ``FeaturePosition``
constructors, and you’ll get back out ``ExactPosition`` objects:

.. cont-doctest

.. code:: pycon

   >>> exact_location = SeqFeature.SimpleLocation(5, 9)
   >>> print(exact_location)
   [5:9]
   >>> exact_location.start
   ExactPosition(5)
   >>> int(exact_location.start)
   5

That is most of the nitty gritty about dealing with fuzzy positions in
Biopython. It has been designed so that dealing with fuzziness is not
that much more complicated than dealing with exact positions, and
hopefully you find that true!

Location testing
^^^^^^^^^^^^^^^^

You can use the Python keyword ``in`` with a ``SeqFeature`` or location
object to see if the base/residue for a parent coordinate is within the
feature/location or not.

For example, suppose you have a SNP of interest and you want to know
which features this SNP is within, and lets suppose this SNP is at index
4350 (Python counting!). Here is a simple brute force solution where we
just check all the features one by one in a loop:

.. doctest ../Tests/GenBank

.. code:: pycon

   >>> from Bio import SeqIO
   >>> my_snp = 4350
   >>> record = SeqIO.read("NC_005816.gb", "genbank")
   >>> for feature in record.features:
   ...     if my_snp in feature:
   ...         print("%s %s" % (feature.type, feature.qualifiers.get("db_xref")))
   ...
   source ['taxon:229193']
   gene ['GeneID:2767712']
   CDS ['GI:45478716', 'GeneID:2767712']

Note that gene and CDS features from GenBank or EMBL files defined with
joins are the union of the exons – they do not cover any introns.

Sequence described by a feature or location
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A ``SeqFeature`` or location object doesn’t directly contain a sequence,
instead the location (see Section :ref:`sec:locations`) describes
how to get this from the parent sequence. For example consider a (short)
gene sequence with location ``5:18`` on the reverse strand, which in
GenBank/EMBL notation using 1-based counting would be
``complement(6..18)``, like this:

.. doctest

.. code:: pycon

   >>> from Bio.Seq import Seq
   >>> from Bio.SeqFeature import SeqFeature, SimpleLocation
   >>> seq = Seq("ACCGAGACGGCAAAGGCTAGCATAGGTATGAGACTTCCTTCCTGCCAGTGCTGAGGAACTGGGAGCCTAC")
   >>> feature = SeqFeature(SimpleLocation(5, 18, strand=-1), type="gene")

You could take the parent sequence, slice it to extract ``5:18``, and then
take the reverse complement. The feature location’s start and end are
integer-like so this works:

.. cont-doctest

.. code:: pycon

   >>> feature_seq = seq[feature.location.start : feature.location.end].reverse_complement()
   >>> print(feature_seq)
   AGCCTTTGCCGTC

This is a simple example so this isn’t too bad – however once you have
to deal with compound features (joins) this is rather messy. Instead,
the ``SeqFeature`` object has an ``extract`` method to take care of all
this (and since Biopython 1.78 can handle trans-splicing by supplying a
dictionary of referenced sequences):

.. cont-doctest

.. code:: pycon

   >>> feature_seq = feature.extract(seq)
   >>> print(feature_seq)
   AGCCTTTGCCGTC

The length of a ``SeqFeature`` or location matches that of the region of
sequence it describes.

.. cont-doctest

.. code:: pycon

   >>> print(len(feature_seq))
   13
   >>> print(len(feature))
   13
   >>> print(len(feature.location))
   13

For ``SimpleLocation`` objects the length is just the difference between
the start and end positions. However, for a ``CompoundLocation`` the
length is the sum of the constituent regions.

Comparison
----------

The ``SeqRecord`` objects can be very complex, but here’s a simple
example:

.. doctest

.. code:: pycon

   >>> from Bio.Seq import Seq
   >>> from Bio.SeqRecord import SeqRecord
   >>> record1 = SeqRecord(Seq("ACGT"), id="test")
   >>> record2 = SeqRecord(Seq("ACGT"), id="test")

What happens when you try to compare these “identical” records?

.. code:: pycon

   >>> record1 == record2

Perhaps surprisingly older versions of Biopython would use Python’s
default object comparison for the ``SeqRecord``, meaning
``record1 == record2`` would only return ``True`` if these variables
pointed at the same object in memory. In this example,
``record1 == record2`` would have returned ``False`` here!

.. code:: pycon

   >>> record1 == record2  # on old versions of Biopython!
   False

As of Biopython 1.67, ``SeqRecord`` comparison like
``record1 == record2`` will instead raise an explicit error to avoid
people being caught out by this:

.. cont-doctest

.. code:: pycon

   >>> record1 == record2
   Traceback (most recent call last):
   ...
   NotImplementedError: SeqRecord comparison is deliberately not implemented. Explicitly compare the attributes of interest.

Instead you should check the attributes you are interested in, for
example the identifier and the sequence:

.. cont-doctest

.. code:: pycon

   >>> record1.id == record2.id
   True
   >>> record1.seq == record2.seq
   True

Beware that comparing complex objects quickly gets complicated (see also
Section :ref:`sec:seq-comparison`).

References
----------

Another common annotation related to a sequence is a reference to a
journal or other published work dealing with the sequence. We have a
fairly simple way of representing a Reference in Biopython – we have a
``Bio.SeqFeature.Reference`` class that stores the relevant information
about a reference as attributes of an object.

The attributes include things that you would expect to see in a
reference like ``journal``, ``title`` and ``authors``. Additionally, it
also can hold the ``medline_id`` and ``pubmed_id`` and a ``comment``
about the reference. These are all accessed simply as attributes of the
object.

A reference also has a ``location`` object so that it can specify a
particular location on the sequence that the reference refers to. For
instance, you might have a journal that is dealing with a particular
gene located on a BAC, and want to specify that it only refers to this
position exactly. The ``location`` is a potentially fuzzy location, as
described in section :ref:`sec:locations`.

Any reference objects are stored as a list in the ``SeqRecord`` object’s
``annotations`` dictionary under the key “references”. That’s all there
is too it. References are meant to be easy to deal with, and hopefully
general enough to cover lots of usage cases.

.. _`sec:SeqRecord-format`:

The format method
-----------------

The ``format()`` method of the ``SeqRecord`` class gives a string
containing your record formatted using one of the output file formats
supported by ``Bio.SeqIO``, such as FASTA:

.. doctest

.. code:: pycon

   >>> from Bio.Seq import Seq
   >>> from Bio.SeqRecord import SeqRecord
   >>> record = SeqRecord(
   ...     Seq(
   ...         "MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVGQALFGD"
   ...         "GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK"
   ...         "NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM"
   ...         "SSAC"
   ...     ),
   ...     id="gi|14150838|gb|AAK54648.1|AF376133_1",
   ...     description="chalcone synthase [Cucumis sativus]",
   ... )
   >>> print(record.format("fasta"))

which should give:

.. cont-doctest

.. code:: pycon

   >gi|14150838|gb|AAK54648.1|AF376133_1 chalcone synthase [Cucumis sativus]
   MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVGQALFGD
   GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK
   NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM
   SSAC
   <BLANKLINE>

This ``format`` method takes a single mandatory argument, a lower case
string which is supported by ``Bio.SeqIO`` as an output format (see
Chapter :ref:`chapter:seqio`). However, some of the file
formats ``Bio.SeqIO`` can write to *require* more than one record
(typically the case for multiple sequence alignment formats), and thus
won’t work via this ``format()`` method. See also
Section :ref:`sec:Bio.SeqIO-and-StringIO`.

.. _`sec:SeqRecord-slicing`:

Slicing a SeqRecord
-------------------

You can slice a ``SeqRecord``, to give you a new ``SeqRecord`` covering
just part of the sequence. What is important here is that any per-letter
annotations are also sliced, and any features which fall completely
within the new sequence are preserved (with their locations adjusted).

For example, taking the same GenBank file used earlier:

.. doctest ../Tests/GenBank

.. code:: pycon

   >>> from Bio import SeqIO
   >>> record = SeqIO.read("NC_005816.gb", "genbank")
   >>> record
   SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])
   >>> len(record)
   9609
   >>> len(record.features)
   41

For this example we’re going to focus in on the ``pim`` gene,
``YP_pPCP05``. If you have a look at the GenBank file directly you’ll
find this gene/CDS has location string ``4343..4780``, or in Python
counting ``4342:4780``. From looking at the file you can work out that
these are the twelfth and thirteenth entries in the file, so in Python
zero-based counting they are entries :math:`11` and :math:`12` in the
``features`` list:

.. cont-doctest

.. code:: pycon

   >>> print(record.features[20])
   type: gene
   location: [4342:4780](+)
   qualifiers:
       Key: db_xref, Value: ['GeneID:2767712']
       Key: gene, Value: ['pim']
       Key: locus_tag, Value: ['YP_pPCP05']
   <BLANKLINE>
   >>> print(record.features[21])
   type: CDS
   location: [4342:4780](+)
   qualifiers:
       Key: codon_start, Value: ['1']
       Key: db_xref, Value: ['GI:45478716', 'GeneID:2767712']
       Key: gene, Value: ['pim']
       Key: locus_tag, Value: ['YP_pPCP05']
       Key: note, Value: ['similar to many previously sequenced pesticin immunity protein entries of Yersinia pestis plasmid pPCP, e.g. gi| 16082683|,ref|NP_395230.1| (NC_003132) , gi|1200166|emb|CAA90861.1| (Z54145 ) , gi|1488655| emb|CAA63439.1| (X92856) , gi|2996219|gb|AAC62543.1| (AF053945) , and gi|5763814|emb|CAB531 67.1| (AL109969)']
       Key: product, Value: ['pesticin immunity protein']
       Key: protein_id, Value: ['NP_995571.1']
       Key: transl_table, Value: ['11']
       Key: translation, Value: ['MGGGMISKLFCLALIFLSSSGLAEKNTYTAKDILQNLELNTFGNSLSHGIYGKQTTFKQTEFTNIKSNTKKHIALINKDNSWMISLKILGIKRDEYTVCFEDFSLIRPPTYVAIHPLLIKKVKSGNFIVVKEIKKSIPGCTVYYH']
   <BLANKLINE>

Let’s slice this parent record from 4300 to 4800 (enough to include the
``pim`` gene/CDS), and see how many features we get:

.. cont-doctest

.. code:: pycon

   >>> sub_record = record[4300:4800]
   >>> sub_record
   SeqRecord(seq=Seq('ATAAATAGATTATTCCAAATAATTTATTTATGTAAGAACAGGATGGGAGGGGGA...TTA'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])
   >>> len(sub_record)
   500
   >>> len(sub_record.features)
   2

Our sub-record just has two features, the gene and CDS entries for
``YP_pPCP05``:

.. cont-doctest

.. code:: pycon

   >>> print(sub_record.features[0])
   type: gene
   location: [42:480](+)
   qualifiers:
       Key: db_xref, Value: ['GeneID:2767712']
       Key: gene, Value: ['pim']
       Key: locus_tag, Value: ['YP_pPCP05']
   <BLANKLINE>
   >>> print(sub_record.features[1])
   type: CDS
   location: [42:480](+)
   qualifiers:
       Key: codon_start, Value: ['1']
       Key: db_xref, Value: ['GI:45478716', 'GeneID:2767712']
       Key: gene, Value: ['pim']
       Key: locus_tag, Value: ['YP_pPCP05']
       Key: note, Value: ['similar to many previously sequenced pesticin immunity protein entries of Yersinia pestis plasmid pPCP, e.g. gi| 16082683|,ref|NP_395230.1| (NC_003132) , gi|1200166|emb|CAA90861.1| (Z54145 ) , gi|1488655| emb|CAA63439.1| (X92856) , gi|2996219|gb|AAC62543.1| (AF053945) , and gi|5763814|emb|CAB531 67.1| (AL109969)']
       Key: product, Value: ['pesticin immunity protein']
       Key: protein_id, Value: ['NP_995571.1']
       Key: transl_table, Value: ['11']
       Key: translation, Value: ['MGGGMISKLFCLALIFLSSSGLAEKNTYTAKDILQNLELNTFGNSLSHGIYGKQTTFKQTEFTNIKSNTKKHIALINKDNSWMISLKILGIKRDEYTVCFEDFSLIRPPTYVAIHPLLIKKVKSGNFIVVKEIKKSIPGCTVYYH']
   <BLANKLINE>

Notice that their locations have been adjusted to reflect the new parent
sequence!

While Biopython has done something sensible and hopefully intuitive with
the features (and any per-letter annotation), for the other annotation
it is impossible to know if this still applies to the sub-sequence or
not. To avoid guessing, with the exception of the molecule type, the
``.annotations`` and ``.dbxrefs`` are omitted from the sub-record, and
it is up to you to transfer any relevant information as appropriate.

.. cont-doctest

.. code:: pycon

   >>> sub_record.annotations
   {'molecule_type': 'DNA'}
   >>> sub_record.dbxrefs
   []

You may wish to preserve other entries like the organism? Beware of
copying the entire annotations dictionary as in this case your partial
sequence is no longer circular DNA - it is now linear:

.. cont-doctest

.. code:: pycon

   >>> sub_record.annotations["topology"] = "linear"

The same point could be made about the record ``id``, ``name`` and
``description``, but for practicality these are preserved:

.. cont-doctest

.. code:: pycon

   >>> sub_record.id
   'NC_005816.1'
   >>> sub_record.name
   'NC_005816'
   >>> sub_record.description
   'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'

This illustrates the problem nicely though, our new sub-record is *not*
the complete sequence of the plasmid, so the description is wrong! Let’s
fix this and then view the sub-record as a reduced GenBank file using
the ``format`` method described above in
Section :ref:`sec:SeqRecord-format`:

.. cont-doctest

.. code:: pycon

   >>> sub_record.description = (
   ...     "Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, partial"
   ... )
   >>> print(sub_record.format("genbank")[:200] + "...")
   LOCUS       NC_005816                500 bp    DNA     linear   UNK 01-JAN-1980
   DEFINITION  Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, partial.
   ACCESSION   NC_005816
   VERSION     NC_0058...

See
Sections :ref:`sec:FASTQ-slicing-off-primer`
and :ref:`sec:FASTQ-slicing-off-adaptor`
for some FASTQ examples where the per-letter annotations (the read
quality scores) are also sliced.

.. _`sec:SeqRecord-addition`:

Adding SeqRecord objects
------------------------

You can add ``SeqRecord`` objects together, giving a new ``SeqRecord``.
What is important here is that any common per-letter annotations are
also added, all the features are preserved (with their locations
adjusted), and any other common annotation is also kept (like the id,
name and description).

For an example with per-letter annotation, we’ll use the first record in
a FASTQ file. Chapter :ref:`chapter:seqio` will explain
the ``SeqIO`` functions:

.. doctest ../Tests/Quality

.. code:: pycon

   >>> from Bio import SeqIO
   >>> record = next(SeqIO.parse("example.fastq", "fastq"))
   >>> len(record)
   25
   >>> print(record.seq)
   CCCTTCTTGTCTTCAGCGTTTCTCC
   >>> print(record.letter_annotations["phred_quality"])
   [26, 26, 18, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 22, 26, 26, 26, 26, 26, 26, 26, 23, 23]

Let’s suppose this was Roche 454 data, and that from other information
you think the ``TTT`` should be only ``TT``. We can make a new edited
record by first slicing the ``SeqRecord`` before and after the “extra”
third ``T``:

.. cont-doctest

.. code:: pycon

   >>> left = record[:20]
   >>> print(left.seq)
   CCCTTCTTGTCTTCAGCGTT
   >>> print(left.letter_annotations["phred_quality"])
   [26, 26, 18, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 22, 26, 26, 26, 26]
   >>> right = record[21:]
   >>> print(right.seq)
   CTCC
   >>> print(right.letter_annotations["phred_quality"])
   [26, 26, 23, 23]

Now add the two parts together:

.. cont-doctest

.. code:: pycon

   >>> edited = left + right
   >>> len(edited)
   24
   >>> print(edited.seq)
   CCCTTCTTGTCTTCAGCGTTCTCC
   >>> print(edited.letter_annotations["phred_quality"])
   [26, 26, 18, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 22, 26, 26, 26, 26, 26, 26, 23, 23]

Easy and intuitive? We hope so! You can make this shorter with just:

.. cont-doctest

.. code:: pycon

   >>> edited = record[:20] + record[21:]

Now, for an example with features, we’ll use a GenBank file. Suppose you
have a circular genome:

.. doctest ../Tests/GenBank

.. code:: pycon

   >>> from Bio import SeqIO
   >>> record = SeqIO.read("NC_005816.gb", "genbank")
   >>> record
   SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])
   >>> len(record)
   9609
   >>> len(record.features)
   41
   >>> record.dbxrefs
   ['Project:58037']
   >>> record.annotations.keys()
   dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references', 'comment'])

You can shift the origin like this:

.. cont-doctest

.. code:: pycon

   >>> shifted = record[2000:] + record[:2000]
   >>> shifted
   SeqRecord(seq=Seq('GATACGCAGTCATATTTTTTACACAATTCTCTAATCCCGACAAGGTCGTAGGTC...GGA'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])
   >>> len(shifted)
   9609

Note that this isn’t perfect in that some annotation like the database
cross references, all the annotations except molecule type, and one of
the features (the source feature) have been lost:

.. cont-doctest

.. code:: pycon

   >>> len(shifted.features)
   40
   >>> shifted.dbxrefs
   []
   >>> shifted.annotations.keys()
   dict_keys(['molecule_type'])

This is because the ``SeqRecord`` slicing step is cautious in what
annotation it preserves (erroneously propagating annotation can cause
major problems). If you want to keep the database cross references or
the annotations dictionary, this must be done explicitly:

.. cont-doctest

.. code:: pycon

   >>> shifted.dbxrefs = record.dbxrefs[:]
   >>> shifted.annotations = record.annotations.copy()
   >>> shifted.dbxrefs
   ['Project:58037']
   >>> shifted.annotations.keys()
   dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references', 'comment'])

Also note that in an example like this, you should probably change the
record identifiers since the NCBI references refer to the *original*
unmodified sequence.

.. _`sec:SeqRecord-reverse-complement`:

Reverse-complementing SeqRecord objects
---------------------------------------

One of the new features in Biopython 1.57 was the ``SeqRecord`` object’s
``reverse_complement`` method. This tries to balance easy of use with
worries about what to do with the annotation in the reverse complemented
record.

For the sequence, this uses the Seq object’s reverse complement method.
Any features are transferred with the location and strand recalculated.
Likewise any per-letter-annotation is also copied but reversed (which
makes sense for typical examples like quality scores). However, transfer
of most annotation is problematical.

For instance, if the record ID was an accession, that accession should
not really apply to the reverse complemented sequence, and transferring
the identifier by default could easily cause subtle data corruption in
downstream analysis. Therefore by default, the ``SeqRecord``’s id,
name, description, annotations and database cross references are all
*not* transferred by default.

The ``SeqRecord`` object’s ``reverse_complement`` method takes a number
of optional arguments corresponding to properties of the record. Setting
these arguments to ``True`` means copy the old values, while ``False``
means drop the old values and use the default value. You can
alternatively provide the new desired value instead.

Consider this example record:

.. doctest ../Tests/GenBank

.. code:: pycon

   >>> from Bio import SeqIO
   >>> rec = SeqIO.read("NC_005816.gb", "genbank")
   >>> print(rec.id, len(rec), len(rec.features), len(rec.dbxrefs), len(rec.annotations))
   NC_005816.1 9609 41 1 13

Here we take the reverse complement and specify a new identifier – but
notice how most of the annotation is dropped (but not the features):

.. cont-doctest

.. code:: pycon

   >>> rc = rec.reverse_complement(id="TESTING")
   >>> print(rc.id, len(rc), len(rc.features), len(rc.dbxrefs), len(rc.annotations))
   TESTING 9609 41 0 0
