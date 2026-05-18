.. _`chapter:msa`:

Multiple Sequence Alignment objects
===================================

This chapter describes the older ``MultipleSeqAlignment`` class and the
parsers in ``Bio.AlignIO`` that parse the output of sequence alignment
software, generating ``MultipleSeqAlignment`` objects. By Multiple
Sequence Alignments we mean a collection of multiple sequences which
have been aligned together – usually with the insertion of gap
characters, and addition of leading or trailing gaps – such that all the
sequence strings are the same length. Such an alignment can be regarded
as a matrix of letters, where each row is held as a ``SeqRecord`` object
internally.

We will introduce the ``MultipleSeqAlignment`` object which holds this
kind of data, and the ``Bio.AlignIO`` module for reading and writing
them as various file formats (following the design of the ``Bio.SeqIO``
module from the previous chapter). Note that both ``Bio.SeqIO`` and
``Bio.AlignIO`` can read and write sequence alignment files. The
appropriate choice will depend largely on what you want to do with the
data.

The final part of this chapter is about using common multiple sequence
alignment tools like ClustalW and MUSCLE from Python, and parsing the
results with Biopython.

Parsing or Reading Sequence Alignments
--------------------------------------

We have two functions for reading in sequence alignments,
``Bio.AlignIO.read()`` and ``Bio.AlignIO.parse()`` which following the
convention introduced in ``Bio.SeqIO`` are for files containing one or
multiple alignments respectively.

Using ``Bio.AlignIO.parse()`` will return an *iterator* which gives
``MultipleSeqAlignment`` objects. Iterators are typically used in a for
loop. Examples of situations where you will have multiple different
alignments include resampled alignments from the PHYLIP tool
``seqboot``, or multiple pairwise alignments from the EMBOSS tools
``water`` or ``needle``, or Bill Pearson’s FASTA tools.

However, in many situations you will be dealing with files which contain
only a single alignment. In this case, you should use the
``Bio.AlignIO.read()`` function which returns a single
``MultipleSeqAlignment`` object.

Both functions expect two mandatory arguments:

#. The first argument is a *handle* to read the data from, typically an
   open file (see
   Section :ref:`sec:appendix-handles`), or a
   filename.

#. The second argument is a lower case string specifying the alignment
   format. As in ``Bio.SeqIO`` we don’t try and guess the file format
   for you! See http://biopython.org/wiki/AlignIO for a full listing of
   supported formats.

There is also an optional ``seq_count`` argument which is discussed in
Section :ref:`sec:AlignIO-count-argument` below for dealing with
ambiguous file formats which may contain more than one alignment.

Single Alignments
~~~~~~~~~~~~~~~~~

As an example, consider the following annotation rich protein alignment
in the PFAM or Stockholm file format:

.. code:: text

   # STOCKHOLM 1.0
   #=GS COATB_BPIKE/30-81  AC P03620.1
   #=GS COATB_BPIKE/30-81  DR PDB; 1ifl ; 1-52;
   #=GS Q9T0Q8_BPIKE/1-52  AC Q9T0Q8.1
   #=GS COATB_BPI22/32-83  AC P15416.1
   #=GS COATB_BPM13/24-72  AC P69541.1
   #=GS COATB_BPM13/24-72  DR PDB; 2cpb ; 1-49;
   #=GS COATB_BPM13/24-72  DR PDB; 2cps ; 1-49;
   #=GS COATB_BPZJ2/1-49   AC P03618.1
   #=GS Q9T0Q9_BPFD/1-49   AC Q9T0Q9.1
   #=GS Q9T0Q9_BPFD/1-49   DR PDB; 1nh4 A; 1-49;
   #=GS COATB_BPIF1/22-73  AC P03619.2
   #=GS COATB_BPIF1/22-73  DR PDB; 1ifk ; 1-50;
   COATB_BPIKE/30-81             AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA
   #=GR COATB_BPIKE/30-81  SS    -HHHHHHHHHHHHHH--HHHHHHHH--HHHHHHHHHHHHHHHHHHHHH----
   Q9T0Q8_BPIKE/1-52             AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA
   COATB_BPI22/32-83             DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA
   COATB_BPM13/24-72             AEGDDP...AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA
   #=GR COATB_BPM13/24-72  SS    ---S-T...CHCHHHHCCCCTCCCTTCHHHHHHHHHHHHHHHHHHHHCTT--
   COATB_BPZJ2/1-49              AEGDDP...AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA
   Q9T0Q9_BPFD/1-49              AEGDDP...AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA
   #=GR Q9T0Q9_BPFD/1-49   SS    ------...-HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH--
   COATB_BPIF1/22-73             FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA
   #=GR COATB_BPIF1/22-73  SS    XX-HHHH--HHHHHH--HHHHHHH--HHHHHHHHHHHHHHHHHHHHHHH---
   #=GC SS_cons                  XHHHHHHHHHHHHHHHCHHHHHHHHCHHHHHHHHHHHHHHHHHHHHHHHC--
   #=GC seq_cons                 AEssss...AptAhDSLpspAT-hIu.sWshVsslVsAsluIKLFKKFsSKA
   //

This is the seed alignment for the Phage_Coat_Gp8 (PF05371) PFAM entry,
downloaded from a now out of date release of PFAM from
https://pfam.xfam.org/. We can load this file as follows (assuming it
has been saved to disk as “PF05371_seed.sth” in the current working
directory):

.. doctest examples

.. code:: pycon

   >>> from Bio import AlignIO
   >>> alignment = AlignIO.read("PF05371_seed.sth", "stockholm")

This code will print out a summary of the alignment:

.. cont-doctest

.. code:: pycon

   >>> print(alignment)
   Alignment with 7 rows and 52 columns
   AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRL...SKA COATB_BPIKE/30-81
   AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKL...SRA Q9T0Q8_BPIKE/1-52
   DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRL...SKA COATB_BPI22/32-83
   AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPM13/24-72
   AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPZJ2/1-49
   AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA Q9T0Q9_BPFD/1-49
   FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA COATB_BPIF1/22-73

You’ll notice in the above output the sequences have been truncated. We
could instead write our own code to format this as we please by
iterating over the rows as ``SeqRecord`` objects:

.. doctest examples

.. code:: pycon

   >>> from Bio import AlignIO
   >>> alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
   >>> print("Alignment length %i" % alignment.get_alignment_length())
   Alignment length 52
   >>> for record in alignment:
   ...     print("%s - %s" % (record.seq, record.id))
   ...
   AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA - COATB_BPIKE/30-81
   AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA - Q9T0Q8_BPIKE/1-52
   DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA - COATB_BPI22/32-83
   AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - COATB_BPM13/24-72
   AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA - COATB_BPZJ2/1-49
   AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - Q9T0Q9_BPFD/1-49
   FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA - COATB_BPIF1/22-73

You could also call Python’s built-in ``format`` function on the
alignment object to show it in a particular file format – see
Section :ref:`sec:alignment-format` for details.

Did you notice in the raw file above that several of the sequences
include database cross-references to the PDB and the associated known
secondary structure? Try this:

.. cont-doctest

.. code:: pycon

   >>> for record in alignment:
   ...     if record.dbxrefs:
   ...         print("%s %s" % (record.id, record.dbxrefs))
   ...
   COATB_BPIKE/30-81 ['PDB; 1ifl ; 1-52;']
   COATB_BPM13/24-72 ['PDB; 2cpb ; 1-49;', 'PDB; 2cps ; 1-49;']
   Q9T0Q9_BPFD/1-49 ['PDB; 1nh4 A; 1-49;']
   COATB_BPIF1/22-73 ['PDB; 1ifk ; 1-50;']

To have a look at all the sequence annotation, try this:

.. code:: pycon

   >>> for record in alignment:
   ...     print(record)
   ...

PFAM provide a nice web interface at http://pfam.xfam.org/family/PF05371
which will actually let you download this alignment in several other
formats. This is what the file looks like in the FASTA file format:

.. code:: text

   >COATB_BPIKE/30-81
   AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA
   >Q9T0Q8_BPIKE/1-52
   AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA
   >COATB_BPI22/32-83
   DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA
   >COATB_BPM13/24-72
   AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA
   >COATB_BPZJ2/1-49
   AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA
   >Q9T0Q9_BPFD/1-49
   AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA
   >COATB_BPIF1/22-73
   FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA

Note the website should have an option about showing gaps as periods
(dots) or dashes, we’ve shown dashes above. Assuming you download and
save this as file “PF05371_seed.faa” then you can load it with almost
exactly the same code:

.. code:: pycon

   >>> from Bio import AlignIO
   >>> alignment = AlignIO.read("PF05371_seed.faa", "fasta")
   >>> print(alignment)

All that has changed in this code is the filename and the format string.
You’ll get the same output as before, the sequences and record
identifiers are the same. However, as you should expect, if you check
each ``SeqRecord`` there is no annotation nor database cross-references
because these are not included in the FASTA file format.

Note that rather than using the Sanger website, you could have used
``Bio.AlignIO`` to convert the original Stockholm format file into a
FASTA file yourself (see below).

With any supported file format, you can load an alignment in exactly the
same way just by changing the format string. For example, use “phylip”
for PHYLIP files, “nexus” for NEXUS files or “emboss” for the alignments
output by the EMBOSS tools. There is a full listing on the wiki page
(http://biopython.org/wiki/AlignIO) and in the built-in documentation,
:py:mod:`Bio.AlignIO`:

.. code:: pycon

   >>> from Bio import AlignIO
   >>> help(AlignIO)

Multiple Alignments
~~~~~~~~~~~~~~~~~~~

The previous section focused on reading files containing a single
alignment. In general however, files can contain more than one
alignment, and to read these files we must use the
``Bio.AlignIO.parse()`` function.

Suppose you have a small alignment in PHYLIP format:

.. code:: text

       5    6
   Alpha     AACAAC
   Beta      AACCCC
   Gamma     ACCAAC
   Delta     CCACCA
   Epsilon   CCAAAC

If you wanted to bootstrap a phylogenetic tree using the PHYLIP tools,
one of the steps would be to create a set of many resampled alignments
using the tool ``bootseq``. This would give output something like this,
which has been abbreviated for conciseness:

.. code:: text

       5     6
   Alpha     AAACCA
   Beta      AAACCC
   Gamma     ACCCCA
   Delta     CCCAAC
   Epsilon   CCCAAA
       5     6
   Alpha     AAACAA
   Beta      AAACCC
   Gamma     ACCCAA
   Delta     CCCACC
   Epsilon   CCCAAA
       5     6
   Alpha     AAAAAC
   Beta      AAACCC
   Gamma     AACAAC
   Delta     CCCCCA
   Epsilon   CCCAAC
   ...
       5     6
   Alpha     AAAACC
   Beta      ACCCCC
   Gamma     AAAACC
   Delta     CCCCAA
   Epsilon   CAAACC

If you wanted to read this in using ``Bio.AlignIO`` you could use:

.. code:: pycon

   >>> from Bio import AlignIO
   >>> alignments = AlignIO.parse("resampled.phy", "phylip")
   >>> for alignment in alignments:
   ...     print(alignment)
   ...     print()
   ...

This would give the following output, again abbreviated for display:

.. code:: text

   Alignment with 5 rows and 6 columns
   AAACCA Alpha
   AAACCC Beta
   ACCCCA Gamma
   CCCAAC Delta
   CCCAAA Epsilon

   Alignment with 5 rows and 6 columns
   AAACAA Alpha
   AAACCC Beta
   ACCCAA Gamma
   CCCACC Delta
   CCCAAA Epsilon

   Alignment with 5 rows and 6 columns
   AAAAAC Alpha
   AAACCC Beta
   AACAAC Gamma
   CCCCCA Delta
   CCCAAC Epsilon

   ...

   Alignment with 5 rows and 6 columns
   AAAACC Alpha
   ACCCCC Beta
   AAAACC Gamma
   CCCCAA Delta
   CAAACC Epsilon

As with the function ``Bio.SeqIO.parse()``, using
``Bio.AlignIO.parse()`` returns an iterator. If you want to keep all the
alignments in memory at once, which will allow you to access them in any
order, then turn the iterator into a list:

.. code:: pycon

   >>> from Bio import AlignIO
   >>> alignments = list(AlignIO.parse("resampled.phy", "phylip"))
   >>> last_align = alignments[-1]
   >>> first_align = alignments[0]

.. _`sec:AlignIO-count-argument`:

Ambiguous Alignments
~~~~~~~~~~~~~~~~~~~~

Many alignment file formats can explicitly store more than one
alignment, and the division between each alignment is clear. However,
when a general sequence file format has been used there is no such block
structure. The most common such situation is when alignments have been
saved in the FASTA file format. For example consider the following:

.. code:: text

   >Alpha
   ACTACGACTAGCTCAG--G
   >Beta
   ACTACCGCTAGCTCAGAAG
   >Gamma
   ACTACGGCTAGCACAGAAG
   >Alpha
   ACTACGACTAGCTCAGG--
   >Beta
   ACTACCGCTAGCTCAGAAG
   >Gamma
   ACTACGGCTAGCACAGAAG

This could be a single alignment containing six sequences (with repeated
identifiers). Or, judging from the identifiers, this is probably two
different alignments each with three sequences, which happen to all have
the same length.

What about this next example?

.. code:: text

   >Alpha
   ACTACGACTAGCTCAG--G
   >Beta
   ACTACCGCTAGCTCAGAAG
   >Alpha
   ACTACGACTAGCTCAGG--
   >Gamma
   ACTACGGCTAGCACAGAAG
   >Alpha
   ACTACGACTAGCTCAGG--
   >Delta
   ACTACGGCTAGCACAGAAG

Again, this could be a single alignment with six sequences. However this
time based on the identifiers we might guess this is three pairwise
alignments which by chance have all got the same lengths.

This final example is similar:

.. code:: text

   >Alpha
   ACTACGACTAGCTCAG--G
   >XXX
   ACTACCGCTAGCTCAGAAG
   >Alpha
   ACTACGACTAGCTCAGG
   >YYY
   ACTACGGCAAGCACAGG
   >Alpha
   --ACTACGAC--TAGCTCAGG
   >ZZZ
   GGACTACGACAATAGCTCAGG

In this third example, because of the differing lengths, this cannot be
treated as a single alignment containing all six records. However, it
could be three pairwise alignments.

Clearly trying to store more than one alignment in a FASTA file is not
ideal. However, if you are forced to deal with these as input files
``Bio.AlignIO`` can cope with the most common situation where all the
alignments have the same number of records. One example of this is a
collection of pairwise alignments, which can be produced by the EMBOSS
tools ``needle`` and ``water`` – although in this situation,
``Bio.AlignIO`` should be able to understand their native output using
“emboss” as the format string.

To interpret these FASTA examples as several separate alignments, we can
use ``Bio.AlignIO.parse()`` with the optional ``seq_count`` argument
which specifies how many sequences are expected in each alignment (in
these examples, 3, 2 and 2 respectively). For example, using the third
example as the input data:

.. code:: pycon

   >>> for alignment in AlignIO.parse(handle, "fasta", seq_count=2):
   ...     print("Alignment length %i" % alignment.get_alignment_length())
   ...     for record in alignment:
   ...         print("%s - %s" % (record.seq, record.id))
   ...     print()
   ...

giving:

.. code:: text

   Alignment length 19
   ACTACGACTAGCTCAG--G - Alpha
   ACTACCGCTAGCTCAGAAG - XXX

   Alignment length 17
   ACTACGACTAGCTCAGG - Alpha
   ACTACGGCAAGCACAGG - YYY

   Alignment length 21
   --ACTACGAC--TAGCTCAGG - Alpha
   GGACTACGACAATAGCTCAGG - ZZZ

Using ``Bio.AlignIO.read()`` or ``Bio.AlignIO.parse()`` without the
``seq_count`` argument would give a single alignment containing all six
records for the first two examples. For the third example, an exception
would be raised because the lengths differ preventing them being turned
into a single alignment.

If the file format itself has a block structure allowing ``Bio.AlignIO``
to determine the number of sequences in each alignment directly, then
the ``seq_count`` argument is not needed. If it is supplied, and doesn’t
agree with the file contents, an error is raised.

Note that this optional ``seq_count`` argument assumes each alignment in
the file has the same number of sequences. Hypothetically you may come
across stranger situations, for example a FASTA file containing several
alignments each with a different number of sequences – although I would
love to hear of a real world example of this. Assuming you cannot get
the data in a nicer file format, there is no straight forward way to
deal with this using ``Bio.AlignIO``. In this case, you could consider
reading in the sequences themselves using ``Bio.SeqIO`` and batching
them together to create the alignments as appropriate.

Writing Alignments
------------------

We’ve talked about using ``Bio.AlignIO.read()`` and
``Bio.AlignIO.parse()`` for alignment input (reading files), and now
we’ll look at ``Bio.AlignIO.write()`` which is for alignment output
(writing files). This is a function taking three arguments: some
``MultipleSeqAlignment`` objects (or for backwards compatibility the
obsolete ``Alignment`` objects), a handle or filename to write to, and a
sequence format.

Here is an example, where we start by creating a few
``MultipleSeqAlignment`` objects the hard way (by hand, rather than by
loading them from a file). Note we create some ``SeqRecord`` objects to
construct the alignment from.

.. doctest

.. code:: pycon

   >>> from Bio.Seq import Seq
   >>> from Bio.SeqRecord import SeqRecord
   >>> from Bio.Align import MultipleSeqAlignment
   >>> align1 = MultipleSeqAlignment(
   ...     [
   ...         SeqRecord(Seq("ACTGCTAGCTAG"), id="Alpha"),
   ...         SeqRecord(Seq("ACT-CTAGCTAG"), id="Beta"),
   ...         SeqRecord(Seq("ACTGCTAGDTAG"), id="Gamma"),
   ...     ]
   ... )
   >>> align2 = MultipleSeqAlignment(
   ...     [
   ...         SeqRecord(Seq("GTCAGC-AG"), id="Delta"),
   ...         SeqRecord(Seq("GACAGCTAG"), id="Epsilon"),
   ...         SeqRecord(Seq("GTCAGCTAG"), id="Zeta"),
   ...     ]
   ... )
   >>> align3 = MultipleSeqAlignment(
   ...     [
   ...         SeqRecord(Seq("ACTAGTACAGCTG"), id="Eta"),
   ...         SeqRecord(Seq("ACTAGTACAGCT-"), id="Theta"),
   ...         SeqRecord(Seq("-CTACTACAGGTG"), id="Iota"),
   ...     ]
   ... )
   >>> my_alignments = [align1, align2, align3]

Now we have a list of ``Alignment`` objects, we’ll write them to a
PHYLIP format file:

.. code:: pycon

   >>> from Bio import AlignIO
   >>> AlignIO.write(my_alignments, "my_example.phy", "phylip")

And if you open this file in your favorite text editor it should look
like this:

.. code:: text

    3 12
   Alpha      ACTGCTAGCT AG
   Beta       ACT-CTAGCT AG
   Gamma      ACTGCTAGDT AG
    3 9
   Delta      GTCAGC-AG
   Epislon    GACAGCTAG
   Zeta       GTCAGCTAG
    3 13
   Eta        ACTAGTACAG CTG
   Theta      ACTAGTACAG CT-
   Iota       -CTACTACAG GTG

Its more common to want to load an existing alignment, and save that,
perhaps after some simple manipulation like removing certain rows or
columns.

Suppose you wanted to know how many alignments the
``Bio.AlignIO.write()`` function wrote to the handle? If your alignments
were in a list like the example above, you could just use
``len(my_alignments)``, however you can’t do that when your records come
from a generator/iterator. Therefore the ``Bio.AlignIO.write()``
function returns the number of alignments written to the file.

*Note* - If you tell the ``Bio.AlignIO.write()`` function to write to a
file that already exists, the old file will be overwritten without any
warning.

.. _`sec:converting-alignments`:

Converting between sequence alignment file formats
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Converting between sequence alignment file formats with ``Bio.AlignIO``
works in the same way as converting between sequence file formats with
``Bio.SeqIO``
(Section :ref:`sec:SeqIO-conversion`). We load
generally the alignment(s) using ``Bio.AlignIO.parse()`` and then save
them using the ``Bio.AlignIO.write()`` – or just use the
``Bio.AlignIO.convert()`` helper function.

For this example, we’ll load the PFAM/Stockholm format file used earlier
and save it as a Clustal W format file:

.. code:: pycon

   >>> from Bio import AlignIO
   >>> count = AlignIO.convert("PF05371_seed.sth", "stockholm", "PF05371_seed.aln", "clustal")
   >>> print("Converted %i alignments" % count)
   Converted 1 alignments

Or, using ``Bio.AlignIO.parse()`` and ``Bio.AlignIO.write()``:

.. code:: pycon

   >>> from Bio import AlignIO
   >>> alignments = AlignIO.parse("PF05371_seed.sth", "stockholm")
   >>> count = AlignIO.write(alignments, "PF05371_seed.aln", "clustal")
   >>> print("Converted %i alignments" % count)
   Converted 1 alignments

The ``Bio.AlignIO.write()`` function expects to be given multiple
alignment objects. In the example above we gave it the alignment
iterator returned by ``Bio.AlignIO.parse()``.

In this case, we know there is only one alignment in the file so we
could have used ``Bio.AlignIO.read()`` instead, but notice we have to
pass this alignment to ``Bio.AlignIO.write()`` as a single element list:

.. code:: pycon

   >>> from Bio import AlignIO
   >>> alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
   >>> AlignIO.write([alignment], "PF05371_seed.aln", "clustal")

Either way, you should end up with the same new Clustal W format file
“PF05371_seed.aln” with the following content:

.. code:: text

   CLUSTAL X (1.81) multiple sequence alignment


   COATB_BPIKE/30-81                   AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSS
   Q9T0Q8_BPIKE/1-52                   AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVS
   COATB_BPI22/32-83                   DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSS
   COATB_BPM13/24-72                   AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTS
   COATB_BPZJ2/1-49                    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFAS
   Q9T0Q9_BPFD/1-49                    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTS
   COATB_BPIF1/22-73                   FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVS

   COATB_BPIKE/30-81                   KA
   Q9T0Q8_BPIKE/1-52                   RA
   COATB_BPI22/32-83                   KA
   COATB_BPM13/24-72                   KA
   COATB_BPZJ2/1-49                    KA
   Q9T0Q9_BPFD/1-49                    KA
   COATB_BPIF1/22-73                   RA

Alternatively, you could make a PHYLIP format file which we’ll name
“PF05371_seed.phy”:

.. code:: pycon

   >>> from Bio import AlignIO
   >>> AlignIO.convert("PF05371_seed.sth", "stockholm", "PF05371_seed.phy", "phylip")

This time the output looks like this:

.. code:: text

    7 52
   COATB_BPIK AEPNAATNYA TEAMDSLKTQ AIDLISQTWP VVTTVVVAGL VIRLFKKFSS
   Q9T0Q8_BPI AEPNAATNYA TEAMDSLKTQ AIDLISQTWP VVTTVVVAGL VIKLFKKFVS
   COATB_BPI2 DGTSTATSYA TEAMNSLKTQ ATDLIDQTWP VVTSVAVAGL AIRLFKKFSS
   COATB_BPM1 AEGDDP---A KAAFNSLQAS ATEYIGYAWA MVVVIVGATI GIKLFKKFTS
   COATB_BPZJ AEGDDP---A KAAFDSLQAS ATEYIGYAWA MVVVIVGATI GIKLFKKFAS
   Q9T0Q9_BPF AEGDDP---A KAAFDSLQAS ATEYIGYAWA MVVVIVGATI GIKLFKKFTS
   COATB_BPIF FAADDATSQA KAAFDSLTAQ ATEMSGYAWA LVVLVVGATV GIKLFKKFVS

              KA
              RA
              KA
              KA
              KA
              KA
              RA

One of the big handicaps of the original PHYLIP alignment file format is
that the sequence identifiers are strictly truncated at ten characters.
In this example, as you can see the resulting names are still unique -
but they are not very readable. As a result, a more relaxed variant of
the original PHYLIP format is now quite widely used:

.. code:: pycon

   >>> from Bio import AlignIO
   >>> AlignIO.convert("PF05371_seed.sth", "stockholm", "PF05371_seed.phy", "phylip-relaxed")

This time the output looks like this, using a longer indentation to
allow all the identifiers to be given in full:

.. code:: text

    7 52
   COATB_BPIKE/30-81  AEPNAATNYA TEAMDSLKTQ AIDLISQTWP VVTTVVVAGL VIRLFKKFSS
   Q9T0Q8_BPIKE/1-52  AEPNAATNYA TEAMDSLKTQ AIDLISQTWP VVTTVVVAGL VIKLFKKFVS
   COATB_BPI22/32-83  DGTSTATSYA TEAMNSLKTQ ATDLIDQTWP VVTSVAVAGL AIRLFKKFSS
   COATB_BPM13/24-72  AEGDDP---A KAAFNSLQAS ATEYIGYAWA MVVVIVGATI GIKLFKKFTS
   COATB_BPZJ2/1-49   AEGDDP---A KAAFDSLQAS ATEYIGYAWA MVVVIVGATI GIKLFKKFAS
   Q9T0Q9_BPFD/1-49   AEGDDP---A KAAFDSLQAS ATEYIGYAWA MVVVIVGATI GIKLFKKFTS
   COATB_BPIF1/22-73  FAADDATSQA KAAFDSLTAQ ATEMSGYAWA LVVLVVGATV GIKLFKKFVS

                      KA
                      RA
                      KA
                      KA
                      KA
                      KA
                      RA

If you have to work with the original strict PHYLIP format, then you may
need to compress the identifiers somehow – or assign your own names or
numbering system. This following bit of code manipulates the record
identifiers before saving the output:

.. code:: pycon

   >>> from Bio import AlignIO
   >>> alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
   >>> name_mapping = {}
   >>> for i, record in enumerate(alignment):
   ...     name_mapping[i] = record.id
   ...     record.id = "seq%i" % i
   ...
   >>> print(name_mapping)
   {0: 'COATB_BPIKE/30-81', 1: 'Q9T0Q8_BPIKE/1-52', 2: 'COATB_BPI22/32-83', 3: 'COATB_BPM13/24-72', 4: 'COATB_BPZJ2/1-49', 5: 'Q9T0Q9_BPFD/1-49', 6: 'COATB_BPIF1/22-73'}
   >>> AlignIO.write([alignment], "PF05371_seed.phy", "phylip")

This code used a Python dictionary to record a simple mapping from the
new sequence system to the original identifier:

.. code:: python

   {
       0: "COATB_BPIKE/30-81",
       1: "Q9T0Q8_BPIKE/1-52",
       2: "COATB_BPI22/32-83",
       # ...
   }

Here is the new (strict) PHYLIP format output:

.. code:: text

    7 52
   seq0       AEPNAATNYA TEAMDSLKTQ AIDLISQTWP VVTTVVVAGL VIRLFKKFSS
   seq1       AEPNAATNYA TEAMDSLKTQ AIDLISQTWP VVTTVVVAGL VIKLFKKFVS
   seq2       DGTSTATSYA TEAMNSLKTQ ATDLIDQTWP VVTSVAVAGL AIRLFKKFSS
   seq3       AEGDDP---A KAAFNSLQAS ATEYIGYAWA MVVVIVGATI GIKLFKKFTS
   seq4       AEGDDP---A KAAFDSLQAS ATEYIGYAWA MVVVIVGATI GIKLFKKFAS
   seq5       AEGDDP---A KAAFDSLQAS ATEYIGYAWA MVVVIVGATI GIKLFKKFTS
   seq6       FAADDATSQA KAAFDSLTAQ ATEMSGYAWA LVVLVVGATV GIKLFKKFVS

              KA
              RA
              KA
              KA
              KA
              KA
              RA

In general, because of the identifier limitation, working with *strict*
PHYLIP file formats shouldn’t be your first choice. Using the
PFAM/Stockholm format on the other hand allows you to record a lot of
additional annotation too.

.. _`sec:alignment-format`:

Getting your alignment objects as formatted strings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``Bio.AlignIO`` interface is based on handles, which means if you
want to get your alignment(s) into a string in a particular file format
you need to do a little bit more work (see below). However, you will
probably prefer to call Python’s built-in ``format`` function on the
alignment object. This takes an output format specification as a single
argument, a lower case string which is supported by ``Bio.AlignIO`` as
an output format. For example:

.. code:: pycon

   >>> from Bio import AlignIO
   >>> alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
   >>> print(format(alignment, "clustal"))
   CLUSTAL X (1.81) multiple sequence alignment


   COATB_BPIKE/30-81                   AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSS
   Q9T0Q8_BPIKE/1-52                   AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVS
   COATB_BPI22/32-83                   DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSS
   ...

Without an output format specification, ``format`` returns the same
output as ``str``.

As described in
Section :ref:`sec:SeqRecord-format`, the
``SeqRecord`` object has a similar method using output formats supported
by ``Bio.SeqIO``.

Internally ``format`` is calling ``Bio.AlignIO.write()`` with a
``StringIO`` handle. You can do this in your own code if for example you
are using an older version of Biopython:

.. code:: pycon

   >>> from io import StringIO
   >>> from Bio import AlignIO
   >>> alignments = AlignIO.parse("PF05371_seed.sth", "stockholm")
   >>> out_handle = StringIO()
   >>> AlignIO.write(alignments, out_handle, "clustal")
   1
   >>> clustal_data = out_handle.getvalue()
   >>> print(clustal_data)
   CLUSTAL X (1.81) multiple sequence alignment


   COATB_BPIKE/30-81                   AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSS
   Q9T0Q8_BPIKE/1-52                   AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVS
   COATB_BPI22/32-83                   DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSS
   COATB_BPM13/24-72                   AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTS
   ...

.. _`sec:manipulating-alignments`:

Manipulating Alignments
-----------------------

Now that we’ve covered loading and saving alignments, we’ll look at what
else you can do with them.

Slicing alignments
~~~~~~~~~~~~~~~~~~

First of all, in some senses the alignment objects act like a Python
``list`` of ``SeqRecord`` objects (the rows). With this model in mind
hopefully the actions of ``len()`` (the number of rows) and iteration
(each row as a ``SeqRecord``) make sense:

.. doctest examples

.. code:: pycon

   >>> from Bio import AlignIO
   >>> alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
   >>> print("Number of rows: %i" % len(alignment))
   Number of rows: 7
   >>> for record in alignment:
   ...     print("%s - %s" % (record.seq, record.id))
   ...
   AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA - COATB_BPIKE/30-81
   AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA - Q9T0Q8_BPIKE/1-52
   DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA - COATB_BPI22/32-83
   AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - COATB_BPM13/24-72
   AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA - COATB_BPZJ2/1-49
   AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - Q9T0Q9_BPFD/1-49
   FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA - COATB_BPIF1/22-73

You can also use the list-like ``append`` and ``extend`` methods to add
more rows to the alignment (as ``SeqRecord`` objects). Keeping the list
metaphor in mind, simple slicing of the alignment should also make sense
- it selects some of the rows giving back another alignment object:

.. cont-doctest

.. code:: pycon

   >>> print(alignment)
   Alignment with 7 rows and 52 columns
   AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRL...SKA COATB_BPIKE/30-81
   AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKL...SRA Q9T0Q8_BPIKE/1-52
   DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRL...SKA COATB_BPI22/32-83
   AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPM13/24-72
   AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPZJ2/1-49
   AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA Q9T0Q9_BPFD/1-49
   FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA COATB_BPIF1/22-73
   >>> print(alignment[3:7])
   Alignment with 4 rows and 52 columns
   AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPM13/24-72
   AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPZJ2/1-49
   AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA Q9T0Q9_BPFD/1-49
   FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA COATB_BPIF1/22-73

What if you wanted to select by column? Those of you who have used the
NumPy matrix or array objects won’t be surprised at this - you use a
double index.

.. cont-doctest

.. code:: pycon

   >>> print(alignment[2, 6])
   T

Using two integer indices pulls out a single letter, short hand for
this:

.. cont-doctest

.. code:: pycon

   >>> print(alignment[2].seq[6])
   T

You can pull out a single column as a string like this:

.. cont-doctest

.. code:: pycon

   >>> print(alignment[:, 6])
   TTT---T

You can also select a range of columns. For example, to pick out those
same three rows we extracted earlier, but take just their first six
columns:

.. cont-doctest

.. code:: pycon

   >>> print(alignment[3:6, :6])
   Alignment with 3 rows and 6 columns
   AEGDDP COATB_BPM13/24-72
   AEGDDP COATB_BPZJ2/1-49
   AEGDDP Q9T0Q9_BPFD/1-49

Leaving the first index as ``:`` means take all the rows:

.. cont-doctest

.. code:: pycon

   >>> print(alignment[:, :6])
   Alignment with 7 rows and 6 columns
   AEPNAA COATB_BPIKE/30-81
   AEPNAA Q9T0Q8_BPIKE/1-52
   DGTSTA COATB_BPI22/32-83
   AEGDDP COATB_BPM13/24-72
   AEGDDP COATB_BPZJ2/1-49
   AEGDDP Q9T0Q9_BPFD/1-49
   FAADDA COATB_BPIF1/22-73

This brings us to a neat way to remove a section. Notice columns 7, 8
and 9 which are gaps in three of the seven sequences:

.. cont-doctest

.. code:: pycon

   >>> print(alignment[:, 6:9])
   Alignment with 7 rows and 3 columns
   TNY COATB_BPIKE/30-81
   TNY Q9T0Q8_BPIKE/1-52
   TSY COATB_BPI22/32-83
   --- COATB_BPM13/24-72
   --- COATB_BPZJ2/1-49
   --- Q9T0Q9_BPFD/1-49
   TSQ COATB_BPIF1/22-73

Again, you can slice to get everything after the ninth column:

.. cont-doctest

.. code:: pycon

   >>> print(alignment[:, 9:])
   Alignment with 7 rows and 43 columns
   ATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
   ATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
   ATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
   AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
   AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
   AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49
   AKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73

Now, the interesting thing is that addition of alignment objects works
by column. This lets you do this as a way to remove a block of columns:

.. cont-doctest

.. code:: pycon

   >>> edited = alignment[:, :6] + alignment[:, 9:]
   >>> print(edited)
   Alignment with 7 rows and 49 columns
   AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
   AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
   DGTSTAATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
   AEGDDPAKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
   AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
   AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49
   FAADDAAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73

Another common use of alignment addition would be to combine alignments
for several different genes into a meta-alignment. Watch out though -
the identifiers need to match up (see
Section :ref:`sec:SeqRecord-addition` for how
adding ``SeqRecord`` objects works). You may find it helpful to first
sort the alignment rows alphabetically by id:

.. cont-doctest

.. code:: pycon

   >>> edited.sort()
   >>> print(edited)
   Alignment with 7 rows and 49 columns
   DGTSTAATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
   FAADDAAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73
   AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
   AEGDDPAKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
   AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
   AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
   AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49

Note that you can only add two alignments together if they have the same
number of rows.

Alignments as arrays
~~~~~~~~~~~~~~~~~~~~

Depending on what you are doing, it can be more useful to turn the
alignment object into an array of letters – and you can do this with
NumPy:

.. doctest examples lib:numpy

.. code:: pycon

   >>> import numpy as np
   >>> from Bio import AlignIO
   >>> alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
   >>> align_array = np.array(alignment)
   >>> print("Array shape %i by %i" % align_array.shape)
   Array shape 7 by 52
   >>> align_array[:, :10]
   array([['A', 'E', 'P', 'N', 'A', 'A', 'T', 'N', 'Y', 'A'],
          ['A', 'E', 'P', 'N', 'A', 'A', 'T', 'N', 'Y', 'A'],
          ['D', 'G', 'T', 'S', 'T', 'A', 'T', 'S', 'Y', 'A'],
          ['A', 'E', 'G', 'D', 'D', 'P', '-', '-', '-', 'A'],
          ['A', 'E', 'G', 'D', 'D', 'P', '-', '-', '-', 'A'],
          ['A', 'E', 'G', 'D', 'D', 'P', '-', '-', '-', 'A'],
          ['F', 'A', 'A', 'D', 'D', 'A', 'T', 'S', 'Q', 'A']],...

Note that this leaves the original Biopython alignment object and the
NumPy array in memory as separate objects - editing one will not update
the other!

Counting substitutions
~~~~~~~~~~~~~~~~~~~~~~

The ``substitutions`` property of an alignment reports how often letters
in the alignment are substituted for each other. This is calculated by
taking all pairs of rows in the alignment, counting the number of times
two letters are aligned to each other, and summing this over all pairs.
For example,

.. doctest

.. code:: pycon

   >>> from Bio.Seq import Seq
   >>> from Bio.SeqRecord import SeqRecord
   >>> from Bio.Align import MultipleSeqAlignment
   >>> msa = MultipleSeqAlignment(
   ...     [
   ...         SeqRecord(Seq("ACTCCTA"), id="seq1"),
   ...         SeqRecord(Seq("AAT-CTA"), id="seq2"),
   ...         SeqRecord(Seq("CCTACT-"), id="seq3"),
   ...         SeqRecord(Seq("TCTCCTC"), id="seq4"),
   ...     ]
   ... )
   >>> print(msa)
   Alignment with 4 rows and 7 columns
   ACTCCTA seq1
   AAT-CTA seq2
   CCTACT- seq3
   TCTCCTC seq4
   >>> substitutions = msa.substitutions
   >>> print(substitutions)
       A    C    T
   A 2.0  4.5  1.0
   C 4.5 10.0  0.5
   T 1.0  0.5 12.0
   <BLANKLINE>

As the ordering of pairs is arbitrary, counts are divided equally above
and below the diagonal. For example, the 9 alignments of ``A`` to ``C``
are stored as 4.5 at position ``['A', 'C']`` and 4.5 at position
``['C', 'A']``. This arrangement helps to make the math easier when
calculating a substitution matrix from these counts, as described in
Section :ref:`sec:substitution_matrices`.

Note that ``msa.substitutions`` contains entries for the letters
appearing in the alignment only. You can use the ``select`` method to
add entries for missing letters, for example

.. cont-doctest

.. code:: pycon

   >>> m = substitutions.select("ATCG")
   >>> print(m)
       A    T    C   G
   A 2.0  1.0  4.5 0.0
   T 1.0 12.0  0.5 0.0
   C 4.5  0.5 10.0 0.0
   G 0.0  0.0  0.0 0.0
   <BLANKLINE>

This also allows you to change the order of letters in the alphabet.

.. _`sec:alignment_newstyle`:

Getting a new-style Alignment object
------------------------------------

Use the ``alignment`` property to create a new-style ``Alignment``
object (see section :ref:`sec:alignmentobject`)
from an old-style ``MultipleSeqAlignment`` object:

.. cont-doctest

.. code:: pycon

   >>> type(msa)
   <class 'Bio.Align.MultipleSeqAlignment'>
   >>> print(msa)
   Alignment with 4 rows and 7 columns
   ACTCCTA seq1
   AAT-CTA seq2
   CCTACT- seq3
   TCTCCTC seq4
   >>> alignment = msa.alignment
   >>> type(alignment)
   <class 'Bio.Align.Alignment'>
   >>> print(alignment)
   seq1              0 ACTCCTA 7
   seq2              0 AAT-CTA 6
   seq3              0 CCTACT- 6
   seq4              0 TCTCCTC 7
   <BLANKLINE>

Note that the ``alignment`` property creates and returns a new
``Alignment`` object that is consistent with the information stored in
the ``MultipleSeqAlignment`` object at the time the ``Alignment`` object
is created. Any changes to the ``MultipleSeqAlignment`` after calling
the ``alignment`` property will not propagate to the ``Alignment``
object. However, you can of course call the ``alignment`` property again
to create a new ``Alignment`` object consistent with the updated
``MultipleSeqAlignment`` object.

.. _`sec:subs_mat_ex`:

Calculating a substitution matrix from a multiple sequence alignment
--------------------------------------------------------------------

You can create your own substitution matrix from an alignment. In this
example, we’ll first read a protein sequence alignment from the Clustalw
file `protein.aln <examples/protein.aln>`__ (also available online
`here <https://raw.githubusercontent.com/biopython/biopython/master/Tests/Clustalw/protein.aln>`__)

.. doctest ../Tests/Clustalw

.. code:: pycon

   >>> from Bio import AlignIO
   >>> filename = "protein.aln"
   >>> msa = AlignIO.read(filename, "clustal")

Section :ref:`sec:alignio_clustal` contains more information on
doing this.

The ``substitutions`` property of the alignment stores the number of
times different residues substitute for each other:

.. cont-doctest

.. code:: pycon

   >>> observed_frequencies = msa.substitutions

To make the example more readable, we’ll select only amino acids with
polar charged side chains:

.. cont-doctest

.. code:: pycon

   >>> observed_frequencies = observed_frequencies.select("DEHKR")
   >>> print(observed_frequencies)
          D      E      H      K      R
   D 2360.0  255.5    7.5    0.5   25.0
   E  255.5 3305.0   16.5   27.0    2.0
   H    7.5   16.5 1235.0   16.0    8.5
   K    0.5   27.0   16.0 3218.0  116.5
   R   25.0    2.0    8.5  116.5 2079.0
   <BLANKLINE>

Rows and columns for other amino acids were removed from the matrix.

Next, we normalize the matrix:

.. cont-doctest

.. code:: pycon

   >>> import numpy as np
   >>> observed_frequencies /= np.sum(observed_frequencies)

Summing over rows or columns gives the relative frequency of occurrence
of each residue:

.. cont-doctest

.. code:: pycon

   >>> residue_frequencies = np.sum(observed_frequencies, 0)
   >>> print(residue_frequencies.format("%.4f"))
   D 0.2015
   E 0.2743
   H 0.0976
   K 0.2569
   R 0.1697
   <BLANKLINE>
   >>> sum(residue_frequencies) == 1.0
   True

The expected frequency of residue pairs is then

.. cont-doctest

.. code:: pycon

   >>> expected_frequencies = np.dot(
   ...     residue_frequencies[:, None], residue_frequencies[None, :]
   ... )
   >>> print(expected_frequencies.format("%.4f"))
          D      E      H      K      R
   D 0.0406 0.0553 0.0197 0.0518 0.0342
   E 0.0553 0.0752 0.0268 0.0705 0.0465
   H 0.0197 0.0268 0.0095 0.0251 0.0166
   K 0.0518 0.0705 0.0251 0.0660 0.0436
   R 0.0342 0.0465 0.0166 0.0436 0.0288
   <BLANKLINE>

Here, ``residue_frequencies[:, None]`` creates a 2D array consisting of
a single column with the values of ``residue_frequencies``, and
``residue_frequencies[None, :]`` a 2D array with these values as a
single row. Taking their dot product (inner product) creates a matrix of
expected frequencies where each entry consists of two
``residue_frequencies`` values multiplied with each other. For example,
``expected_frequencies['D', 'E']`` is equal to
``residue_frequencies['D'] * residue_frequencies['E']``.

We can now calculate the log-odds matrix by dividing the observed
frequencies by the expected frequencies and taking the logarithm:

.. cont-doctest

.. code:: pycon

   >>> m = np.log2(observed_frequencies / expected_frequencies)
   >>> print(m)
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
   >>> aligner.substitution_matrix = m
   >>> aligner.gap_score = -3.0
   >>> alignments = aligner.align("DEHEK", "DHHKK")
   >>> print(alignments[0])
   target            0 DEHEK 5
                     0 |.|.| 5
   query             0 DHHKK 5
   <BLANKLINE>
   >>> print("%.2f" % alignments.score)
   -2.18
   >>> score = m["D", "D"] + m["E", "H"] + m["H", "H"] + m["E", "K"] + m["K", "K"]
   >>> print("%.2f" % score)
   -2.18

.. _`sec:alignment-tools`:

Alignment Tools
---------------

There are *lots* of algorithms out there for aligning sequences, both
pairwise alignments and multiple sequence alignments. These calculations
are relatively slow, and you generally wouldn’t want to write such an
algorithm in Python. For pairwise alignments, you can use Biopython’s
``PairwiseAligner`` (see
Chapter :ref:`chapter:pairwise`), which is
implemented in C and therefore fast. Alternatively, you can run an
external alignment program by invoking it from Python. Normally you
would:

#. Prepare an input file of your unaligned sequences, typically this
   will be a FASTA file which you might create using ``Bio.SeqIO`` (see
   Chapter :ref:`chapter:seqio`).

#. Run the alignment program by running its command using Python’s
   ``subprocess`` module.

#. Read the output from the tool, i.e. your aligned sequences, typically
   using ``Bio.AlignIO`` (see earlier in this chapter).

Here, we will show a few examples of this workflow.

.. _`sec:alignio_clustal`:

ClustalW
~~~~~~~~

ClustalW is a popular command line tool for multiple sequence alignment
(there is also a graphical interface called ClustalX). Before trying to
use ClustalW from within Python, you should first try running the
ClustalW tool yourself by hand at the command line, to familiarize
yourself the other options.

For the most basic usage, all you need is to have a FASTA input file,
such as
`opuntia.fasta <https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/opuntia.fasta>`__
(available online or in the Doc/examples subdirectory of the Biopython
source code). This is a small FASTA file containing seven prickly-pear
DNA sequences (from the cactus family *Opuntia*). By default ClustalW
will generate an alignment and guide tree file with names based on the
input FASTA file, in this case ``opuntia.aln`` and ``opuntia.dnd``, but
you can override this or make it explicit:

.. code:: pycon

   >>> import subprocess
   >>> cmd = "clustalw2 -infile=opuntia.fasta"
   >>> results = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, text=True)

Notice here we have given the executable name as ``clustalw2``,
indicating we have version two installed, which has a different filename
to version one (``clustalw``, the default). Fortunately both versions
support the same set of arguments at the command line (and indeed,
should be functionally identical).

You may find that even though you have ClustalW installed, the above
command doesn’t work – you may get a message about “command not found”
(especially on Windows). This indicated that the ClustalW executable is
not on your PATH (an environment variable, a list of directories to be
searched). You can either update your PATH setting to include the
location of your copy of ClustalW tools (how you do this will depend on
your OS), or simply type in the full path of the tool. Remember, in
Python strings ``\n`` and ``\t`` are by default interpreted as a new
line and a tab – which is why we’re put a letter “r” at the start for a
raw string that isn’t translated in this way. This is generally good
practice when specifying a Windows style file name.

.. code:: pycon

   >>> import os
   >>> clustalw_exe = r"C:\Program Files\new clustal\clustalw2.exe"
   >>> assert os.path.isfile(clustalw_exe), "Clustal W executable missing"
   >>> cmd = clustalw_exe + " -infile=opuntia.fasta"
   >>> results = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, text=True)

Now, at this point it helps to know about how command line tools “work”.
When you run a tool at the command line, it will often print text output
directly to screen. This text can be captured or redirected, via two
“pipes”, called standard output (the normal results) and standard error
(for error messages and debug messages). There is also standard input,
which is any text fed into the tool. These names get shortened to stdin,
stdout and stderr. When the tool finishes, it has a return code (an
integer), which by convention is zero for success, while a non-zero
return code indicates that an error has occurred.

In the example of ClustalW above, when run at the command line all the
important output is written directly to the output files. Everything
normally printed to screen while you wait is captured in
``results.stdout`` and ``results.stderr``, while the return code is
stored in ``results.returncode``.

What we care about are the two output files, the alignment and the guide
tree. We didn’t tell ClustalW what filenames to use, but it defaults to
picking names based on the input file. In this case the output should be
in the file ``opuntia.aln``. You should be able to work out how to read
in the alignment using ``Bio.AlignIO`` by now:

.. doctest examples

.. code:: pycon

   >>> from Bio import AlignIO
   >>> align = AlignIO.read("opuntia.aln", "clustal")
   >>> print(align)
   Alignment with 7 rows and 906 columns
   TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273285|gb|AF191659.1|AF191
   TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273284|gb|AF191658.1|AF191
   TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273287|gb|AF191661.1|AF191
   TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273286|gb|AF191660.1|AF191
   TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273290|gb|AF191664.1|AF191
   TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273289|gb|AF191663.1|AF191
   TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273291|gb|AF191665.1|AF191

In case you are interested (and this is an aside from the main thrust of
this chapter), the ``opuntia.dnd`` file ClustalW creates is just a
standard Newick tree file, and ``Bio.Phylo`` can parse these:

.. doctest examples

.. code:: pycon

   >>> from Bio import Phylo
   >>> tree = Phylo.read("opuntia.dnd", "newick")
   >>> Phylo.draw_ascii(tree)
                                _______________ gi|6273291|gb|AF191665.1|AF191665
     __________________________|
    |                          |   ______ gi|6273290|gb|AF191664.1|AF191664
    |                          |__|
    |                             |_____ gi|6273289|gb|AF191663.1|AF191663
    |
   _|_________________ gi|6273287|gb|AF191661.1|AF191661
    |
    |__________ gi|6273286|gb|AF191660.1|AF191660
    |
    |    __ gi|6273285|gb|AF191659.1|AF191659
    |___|
        | gi|6273284|gb|AF191658.1|AF191658
   <BLANKLINE>

Chapter :ref:`chapter:phylo` covers Biopython’s support
for phylogenetic trees in more depth.

MUSCLE
~~~~~~

MUSCLE is a more recent multiple sequence alignment tool than ClustalW.
As before, we recommend you try using MUSCLE from the command line
before trying to run it from Python.

For the most basic usage, all you need is to have a FASTA input file,
such as
`opuntia.fasta <https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/opuntia.fasta>`__
(available online or in the Doc/examples subdirectory of the Biopython
source code). You can then tell MUSCLE to read in this FASTA file, and
write the alignment to an output file named ``opuntia.txt``:

.. code:: pycon

   >>> import subprocess
   >>> cmd = "muscle -align opuntia.fasta -output opuntia.txt"
   >>> results = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, text=True)

MUSCLE will output the alignment as a FASTA file (using gapped
sequences). The ``Bio.AlignIO`` module is able to read this alignment
using ``format="fasta"``:

.. code:: pycon

   >>> from Bio import AlignIO
   >>> align = AlignIO.read("opuntia.txt", "fasta")
   >>> print(align)
   Alignment with 7 rows and 906 columns
   TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273289|gb|AF191663.1|AF191663
   TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273291|gb|AF191665.1|AF191665
   TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273290|gb|AF191664.1|AF191664
   TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273287|gb|AF191661.1|AF191661
   TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273286|gb|AF191660.1|AF191660
   TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273285|gb|AF191659.1|AF191659
   TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273284|gb|AF191658.1|AF191658

You can also set the other optional parameters; see MUSCLE’s built-in
help for details.

.. _`sec:emboss-needle-water`:

EMBOSS needle and water
~~~~~~~~~~~~~~~~~~~~~~~

The `EMBOSS <http://emboss.sourceforge.net/>`__ suite includes the
``water`` and ``needle`` tools for Smith-Waterman algorithm local
alignment, and Needleman-Wunsch global alignment. The tools share the
same style interface, so switching between the two is trivial – we’ll
just use ``needle`` here.

Suppose you want to do a global pairwise alignment between two
sequences, prepared in FASTA format as follows:

.. code:: text

   >HBA_HUMAN
   MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG
   KKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTP
   AVHASLDKFLASVSTVLTSKYR

in a file ``alpha.faa``, and secondly in a file ``beta.faa``:

.. code:: text

   >HBB_HUMAN
   MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK
   VKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFG
   KEFTPPVQAAYQKVVAGVANALAHKYH

You can find copies of these example files with the Biopython source
code under the ``Doc/examples/`` directory.

The command to align these two sequences against each other using
``needle`` is as follows:

.. code:: text

   needle -outfile=needle.txt -asequence=alpha.faa -bsequence=beta.faa -gapopen=10 -gapextend=0.5

Why not try running this by hand at the command prompt? You should see
it does a pairwise comparison and records the output in the file
``needle.txt`` (in the default EMBOSS alignment file format).

Even if you have EMBOSS installed, running this command may not work –
you might get a message about “command not found” (especially on
Windows). This probably means that the EMBOSS tools are not on your PATH
environment variable. You can either update your PATH setting, or simply
use the full path to the tool, for example:

.. code:: text

   C:\EMBOSS\needle.exe -outfile=needle.txt -asequence=alpha.faa -bsequence=beta.faa -gapopen=10 -gapextend=0.5

Next we want to use Python to run this command for us. As explained
above, for full control, we recommend you use Python’s built-in
``subprocess`` module:

.. code:: pycon

   >>> import sys
   >>> import subprocess
   >>> cmd = "needle -outfile=needle.txt -asequence=alpha.faa -bsequence=beta.faa -gapopen=10 -gapextend=0.5"
   >>> results = subprocess.run(
   ...     cmd,
   ...     stdout=subprocess.PIPE,
   ...     stderr=subprocess.PIPE,
   ...     text=True,
   ...     shell=(sys, platform != "win32"),
   ... )
   >>> print(results.stdout)

   >>> print(results.stderr)
   Needleman-Wunsch global alignment of two sequences

Next we can load the output file with ``Bio.AlignIO`` as discussed
earlier in this chapter, as the ``emboss`` format:

.. code:: pycon

   >>> from Bio import AlignIO
   >>> align = AlignIO.read("needle.txt", "emboss")
   >>> print(align)
   Alignment with 2 rows and 149 columns
   MV-LSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTY...KYR HBA_HUMAN
   MVHLTPEEKSAVTALWGKV--NVDEVGGEALGRLLVVYPWTQRF...KYH HBB_HUMAN

In this example, we told EMBOSS to write the output to a file, but you
*can* tell it to write the output to stdout instead (useful if you don’t
want a temporary output file to get rid of – use ``outfile=stdout``
argument):

.. code:: pycon

   >>> cmd = "needle -outfile=stdout -asequence=alpha.faa -bsequence=beta.faa -gapopen=10 -gapextend=0.5"
   >>> child = subprocess.Popen(
   ...     cmd,
   ...     stdout=subprocess.PIPE,
   ...     stderr=subprocess.PIPE,
   ...     text=True,
   ...     shell=(sys.platform != "win32"),
   ... )
   >>> align = AlignIO.read(child.stdout, "emboss")
   >>> print(align)
   Alignment with 2 rows and 149 columns
   MV-LSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTY...KYR HBA_HUMAN
   MVHLTPEEKSAVTALWGKV--NVDEVGGEALGRLLVVYPWTQRF...KYH HBB_HUMAN

Similarly, it is possible to read *one* of the inputs from stdin (e.g.
``asequence="stdin"``).

This has only scratched the surface of what you can do with ``needle``
and ``water``. One useful trick is that the second file can contain
multiple sequences (say five), and then EMBOSS will do five pairwise
alignments.
