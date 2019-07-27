Multiple Sequence Alignment objects
===================================

This chapter is about Multiple Sequence Alignments, by which we mean a
collection of multiple sequences which have been aligned together –
usually with the insertion of gap characters, and addition of leading or
trailing gaps – such that all the sequence strings are the same length.
Such an alignment can be regarded as a matrix of letters, where each row
is held as a ``SeqRecord`` object internally.

We will introduce the ``MultipleSeqAlignment`` object which holds this
kind of data, and the ``Bio.AlignIO`` module for reading and writing
them as various file formats (following the design of the ``Bio.SeqIO``
module from the previous chapter). Note that both ``Bio.SeqIO`` and
``Bio.AlignIO`` can read and write sequence alignment files. The
appropriate choice will depend largely on what you want to do with the
data.

The final part of this chapter is about our command line wrappers for
common multiple sequence alignment tools like ClustalW and MUSCLE.

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
   open file (see Section [sec:appendix-handles]), or a filename.

#. The second argument is a lower case string specifying the alignment
   format. As in ``Bio.SeqIO`` we don’t try and guess the file format
   for you! See http://biopython.org/wiki/AlignIO for a full listing of
   supported formats.

There is also an optional ``seq_count`` argument which is discussed in
Section [sec:AlignIO-count-argument] below for dealing with ambiguous
file formats which may contain more than one alignment.

A further optional ``alphabet`` argument allowing you to specify the
expected alphabet. This can be useful as many alignment file formats do
not explicitly label the sequences as RNA, DNA or protein – which means
``Bio.AlignIO`` will default to using a generic alphabet.

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

This is the seed alignment for the Phage\_Coat\_Gp8 (PF05371) PFAM
entry, downloaded from a now out of date release of PFAM from
https://pfam.xfam.org/. We can load this file as follows (assuming it
has been saved to disk as “PF05371\_seed.sth” in the current working
directory):

.. doctest examples

.. code:: pycon

    >>> from Bio import AlignIO
    >>> alignment = AlignIO.read("PF05371_seed.sth", "stockholm")

This code will print out a summary of the alignment:

.. cont-doctest

.. code:: pycon

    >>> print(alignment)
    SingleLetterAlphabet() alignment with 7 rows and 52 columns
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

You could also use the alignment object’s ``format`` method to show it
in a particular file format – see Section [sec:alignment-format-method]
for details.

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

PFAM provide a nice web interface at
` http://pfam.xfam.org/family/PF05371 < http://pfam.xfam.org/family/PF05371>`__
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
save this as file “PF05371\_seed.faa” then you can load it with almost
exactly the same code:

.. code:: python

    from Bio import AlignIO
    alignment = AlignIO.read("PF05371_seed.faa", "fasta")
    print(alignment)

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
(http://biopython.org/wiki/AlignIO) and in the built in documentation
(also
`online <http://biopython.org/DIST/docs/api/Bio.AlignIO-module.html>`__):

.. code:: pycon

    >>> from Bio import AlignIO
    >>> help(AlignIO)
    ...

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

.. code:: python

    from Bio import AlignIO
    alignments = AlignIO.parse("resampled.phy", "phylip")
    for alignment in alignments:
        print(alignment)
        print("")

This would give the following output, again abbreviated for display:

.. code:: text

    SingleLetterAlphabet() alignment with 5 rows and 6 columns
    AAACCA Alpha
    AAACCC Beta
    ACCCCA Gamma
    CCCAAC Delta
    CCCAAA Epsilon

    SingleLetterAlphabet() alignment with 5 rows and 6 columns
    AAACAA Alpha
    AAACCC Beta
    ACCCAA Gamma
    CCCACC Delta
    CCCAAA Epsilon

    SingleLetterAlphabet() alignment with 5 rows and 6 columns
    AAAAAC Alpha
    AAACCC Beta
    AACAAC Gamma
    CCCCCA Delta
    CCCAAC Epsilon

    ...

    SingleLetterAlphabet() alignment with 5 rows and 6 columns
    AAAACC Alpha
    ACCCCC Beta
    AAAACC Gamma
    CCCCAA Delta
    CAAACC Epsilon

As with the function ``Bio.SeqIO.parse()``, using
``Bio.AlignIO.parse()`` returns an iterator. If you want to keep all the
alignments in memory at once, which will allow you to access them in any
order, then turn the iterator into a list:

.. code:: python

    from Bio import AlignIO
    alignments = list(AlignIO.parse("resampled.phy", "phylip"))
    last_align = alignments[-1]
    first_align = alignments[0]

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

.. code:: python

    for alignment in AlignIO.parse(handle, "fasta", seq_count=2):
        print("Alignment length %i" % alignment.get_alignment_length())
        for record in alignment:
            print("%s - %s" % (record.seq, record.id))
        print("")

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

.. code:: python

    from Bio.Alphabet import generic_dna
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import MultipleSeqAlignment

    align1 = MultipleSeqAlignment([
                 SeqRecord(Seq("ACTGCTAGCTAG", generic_dna), id="Alpha"),
                 SeqRecord(Seq("ACT-CTAGCTAG", generic_dna), id="Beta"),
                 SeqRecord(Seq("ACTGCTAGDTAG", generic_dna), id="Gamma"),
             ])

    align2 = MultipleSeqAlignment([
                 SeqRecord(Seq("GTCAGC-AG", generic_dna), id="Delta"),
                 SeqRecord(Seq("GACAGCTAG", generic_dna), id="Epsilon"),
                 SeqRecord(Seq("GTCAGCTAG", generic_dna), id="Zeta"),
             ])

    align3 = MultipleSeqAlignment([
                 SeqRecord(Seq("ACTAGTACAGCTG", generic_dna), id="Eta"),
                 SeqRecord(Seq("ACTAGTACAGCT-", generic_dna), id="Theta"),
                 SeqRecord(Seq("-CTACTACAGGTG", generic_dna), id="Iota"),
             ])

    my_alignments = [align1, align2, align3]

Now we have a list of ``Alignment`` objects, we’ll write them to a
PHYLIP format file:

.. code:: python

    from Bio import AlignIO
    AlignIO.write(my_alignments, "my_example.phy", "phylip")

And if you open this file in your favourite text editor it should look
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

Converting between sequence alignment file formats
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Converting between sequence alignment file formats with ``Bio.AlignIO``
works in the same way as converting between sequence file formats with
``Bio.SeqIO`` (Section [sec:SeqIO-conversion]). We load generally the
alignment(s) using ``Bio.AlignIO.parse()`` and then save them using the
``Bio.AlignIO.write()`` – or just use the ``Bio.AlignIO.convert()``
helper function.

For this example, we’ll load the PFAM/Stockholm format file used earlier
and save it as a Clustal W format file:

.. code:: python

    from Bio import AlignIO
    count = AlignIO.convert("PF05371_seed.sth", "stockholm", "PF05371_seed.aln", "clustal")
    print("Converted %i alignments" % count)

Or, using ``Bio.AlignIO.parse()`` and ``Bio.AlignIO.write()``:

.. code:: python

    from Bio import AlignIO
    alignments = AlignIO.parse("PF05371_seed.sth", "stockholm")
    count = AlignIO.write(alignments, "PF05371_seed.aln", "clustal")
    print("Converted %i alignments" % count)

The ``Bio.AlignIO.write()`` function expects to be given multiple
alignment objects. In the example above we gave it the alignment
iterator returned by ``Bio.AlignIO.parse()``.

In this case, we know there is only one alignment in the file so we
could have used ``Bio.AlignIO.read()`` instead, but notice we have to
pass this alignment to ``Bio.AlignIO.write()`` as a single element list:

.. code:: python

    from Bio import AlignIO
    alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
    AlignIO.write([alignment], "PF05371_seed.aln", "clustal")

Either way, you should end up with the same new Clustal W format file
“PF05371\_seed.aln” with the following content:

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
“PF05371\_seed.phy”:

.. code:: python

    from Bio import AlignIO
    AlignIO.convert("PF05371_seed.sth", "stockholm", "PF05371_seed.phy", "phylip")

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

.. code:: python

    from Bio import AlignIO
    AlignIO.convert("PF05371_seed.sth", "stockholm", "PF05371_seed.phy", "phylip-relaxed")

This time the output looks like this, using a longer indentation to
allow all the identifers to be given in full:

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
need to compress the identifers somehow – or assign your own names or
numbering system. This following bit of code manipulates the record
identifiers before saving the output:

.. code:: python

    from Bio import AlignIO
    alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
    name_mapping = {}
    for i, record in enumerate(alignment):
        name_mapping[i] = record.id
        record.id = "seq%i" % i
    print(name_mapping)

    AlignIO.write([alignment], "PF05371_seed.phy", "phylip")

This code used a Python dictionary to record a simple mapping from the
new sequence system to the original identifier:

.. code:: python

    {0: 'COATB_BPIKE/30-81', 1: 'Q9T0Q8_BPIKE/1-52', 2: 'COATB_BPI22/32-83', ...}

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

Getting your alignment objects as formatted strings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``Bio.AlignIO`` interface is based on handles, which means if you
want to get your alignment(s) into a string in a particular file format
you need to do a little bit more work (see below). However, you will
probably prefer to take advantage of the alignment object’s ``format()``
method. This takes a single mandatory argument, a lower case string
which is supported by ``Bio.AlignIO`` as an output format. For example:

.. code:: python

    from Bio import AlignIO
    alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
    print(alignment.format("clustal"))

As described in Section [sec:SeqRecord-format], the ``SeqRecord`` object
has a similar method using output formats supported by ``Bio.SeqIO``.

Internally the ``format()`` method is using the ``StringIO`` string
based handle and calling ``Bio.AlignIO.write()``. You can do this in
your own code if for example you are using an older version of
Biopython:

.. code:: python

    from Bio import AlignIO
    from StringIO import StringIO

    alignments = AlignIO.parse("PF05371_seed.sth", "stockholm")

    out_handle = StringIO()
    AlignIO.write(alignments, out_handle, "clustal")
    clustal_data = out_handle.getvalue()

    print(clustal_data)

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
    SingleLetterAlphabet() alignment with 7 rows and 52 columns
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRL...SKA COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKL...SRA Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRL...SKA COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA COATB_BPIF1/22-73
    >>> print(alignment[3:7])
    SingleLetterAlphabet() alignment with 4 rows and 52 columns
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
    SingleLetterAlphabet() alignment with 3 rows and 6 columns
    AEGDDP COATB_BPM13/24-72
    AEGDDP COATB_BPZJ2/1-49
    AEGDDP Q9T0Q9_BPFD/1-49

Leaving the first index as ``:`` means take all the rows:

.. cont-doctest

.. code:: pycon

    >>> print(alignment[:, :6])
    SingleLetterAlphabet() alignment with 7 rows and 6 columns
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
    SingleLetterAlphabet() alignment with 7 rows and 3 columns
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
    SingleLetterAlphabet() alignment with 7 rows and 43 columns
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
    SingleLetterAlphabet() alignment with 7 rows and 49 columns
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
    DGTSTAATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
    AEGDDPAKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49
    FAADDAAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73

Another common use of alignment addition would be to combine alignments
for several different genes into a meta-alignment. Watch out though -
the identifiers need to match up (see Section [sec:SeqRecord-addition]
for how adding ``SeqRecord`` objects works). You may find it helpful to
first sort the alignment rows alphabetically by id:

.. cont-doctest

.. code:: pycon

    >>> edited.sort()
    >>> print(edited)
    SingleLetterAlphabet() alignment with 7 rows and 49 columns
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
    >>> align_array = np.array([list(rec) for rec in alignment], np.character)
    >>> print("Array shape %i by %i" % align_array.shape)
    Array shape 7 by 52

If you will be working heavily with the columns, you can tell NumPy to
store the array by column (as in Fortran) rather then its default of by
row (as in C):

.. code:: pycon

    >>> align_array = np.array([list(rec) for rec in alignment], np.character, order="F")

Note that this leaves the original Biopython alignment object and the
NumPy array in memory as separate objects - editing one will not update
the other!

Alignment Tools
---------------

There are *lots* of algorithms out there for aligning sequences, both
pairwise alignments and multiple sequence alignments. These calculations
are relatively slow, and you generally wouldn’t want to write such an
algorithm in Python. For pairwise alignments Biopython contains the
``Bio.pairwise2`` module , which is supplemented by functions written in
C for speed enhancements and the new ``PairwiseAligner`` (see
Section [sec:pairwise]). In addition, you can use Biopython to invoke a
command line tool on your behalf. Normally you would:

#. Prepare an input file of your unaligned sequences, typically this
   will be a FASTA file which you might create using ``Bio.SeqIO`` (see
   Chapter [chapter:seqio]).

#. Call the command line tool to process this input file, typically via
   one of Biopython’s command line wrappers (which we’ll discuss here).

#. Read the output from the tool, i.e. your aligned sequences, typically
   using ``Bio.AlignIO`` (see earlier in this chapter).

All the command line wrappers we’re going to talk about in this chapter
follow the same style. You create a command line object specifying the
options (e.g. the input filename and the output filename), then invoke
this command line via a Python operating system call (e.g. using the
``subprocess`` module).

Most of these wrappers are defined in the ``Bio.Align.Applications``
module:

.. code:: pycon

    >>> import Bio.Align.Applications
    >>> dir(Bio.Align.Applications)
    ...
    ['ClustalwCommandline', 'DialignCommandline', 'MafftCommandline', 'MuscleCommandline',
    'PrankCommandline', 'ProbconsCommandline', 'TCoffeeCommandline' ...]

(Ignore the entries starting with an underscore – these have special
meaning in Python.) The module ``Bio.Emboss.Applications`` has wrappers
for some of the `EMBOSS suite <http://emboss.sourceforge.net/>`__,
including ``needle`` and ``water``, which are described below in
Section [sec:emboss-needle-water], and wrappers for the EMBOSS packaged
versions of the PHYLIP tools (which EMBOSS refer to as one of their
EMBASSY packages - third party tools with an EMBOSS style interface). We
won’t explore all these alignment tools here in the section, just a
sample, but the same principles apply.

ClustalW
~~~~~~~~

ClustalW is a popular command line tool for multiple sequence alignment
(there is also a graphical interface called ClustalX). Biopython’s
``Bio.Align.Applications`` module has a wrapper for this alignment tool
(and several others).

Before trying to use ClustalW from within Python, you should first try
running the ClustalW tool yourself by hand at the command line, to
familiarise yourself the other options. You’ll find the Biopython
wrapper is very faithful to the actual command line API:

.. code:: pycon

    >>> from Bio.Align.Applications import ClustalwCommandline
    >>> help(ClustalwCommandline)
    ...

For the most basic usage, all you need is to have a FASTA input file,
such as
`opuntia.fasta <https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/opuntia.fasta>`__
(available online or in the Doc/examples subdirectory of the Biopython
source code). This is a small FASTA file containing seven prickly-pear
DNA sequences (from the cactus family *Opuntia*).

By default ClustalW will generate an alignment and guide tree file with
names based on the input FASTA file, in this case ``opuntia.aln`` and
``opuntia.dnd``, but you can override this or make it explicit:

.. doctest

.. code:: pycon

    >>> from Bio.Align.Applications import ClustalwCommandline
    >>> cline = ClustalwCommandline("clustalw2", infile="opuntia.fasta")
    >>> print(cline)
    clustalw2 -infile=opuntia.fasta

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
your OS), or simply type in the full path of the tool. For example:

.. doctest

.. code:: pycon

    >>> import os
    >>> from Bio.Align.Applications import ClustalwCommandline
    >>> clustalw_exe = r"C:\Program Files\new clustal\clustalw2.exe"
    >>> clustalw_cline = ClustalwCommandline(clustalw_exe, infile="opuntia.fasta")

.. code:: pycon

    >>> assert os.path.isfile(clustalw_exe), "Clustal W executable missing"
    >>> stdout, stderr = clustalw_cline()

Remember, in Python strings ``\n`` and ``\t`` are by default interpreted
as a new line and a tab – which is why we’re put a letter “r” at the
start for a raw string that isn’t translated in this way. This is
generally good practice when specifying a Windows style file name.

Internally this uses the ``subprocess`` module which is now the
recommended way to run another program in Python. This replaces older
options like the ``os.system()`` and the ``os.popen*`` functions.

Now, at this point it helps to know about how command line tools “work”.
When you run a tool at the command line, it will often print text output
directly to screen. This text can be captured or redirected, via two
“pipes”, called standard output (the normal results) and standard error
(for error messages and debug messages). There is also standard input,
which is any text fed into the tool. These names get shortened to stdin,
stdout and stderr. When the tool finishes, it has a return code (an
integer), which by convention is zero for success.

When you run the command line tool like this via the Biopython wrapper,
it will wait for it to finish, and check the return code. If this is non
zero (indicating an error), an exception is raised. The wrapper then
returns two strings, stdout and stderr.

In the case of ClustalW, when run at the command line all the important
output is written directly to the output files. Everything normally
printed to screen while you wait (via stdout or stderr) is boring and
can be ignored (assuming it worked).

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
    SingleLetterAlphabet() alignment with 7 rows and 906 columns
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

Chapter [chapter:phylo] covers Biopython’s support for phylogenetic
trees in more depth.

MUSCLE
~~~~~~

MUSCLE is a more recent multiple sequence alignment tool than ClustalW,
and Biopython also has a wrapper for it under the
``Bio.Align.Applications`` module. As before, we recommend you try using
MUSCLE from the command line before trying it from within Python, as the
Biopython wrapper is very faithful to the actual command line API:

.. code:: pycon

    >>> from Bio.Align.Applications import MuscleCommandline
    >>> help(MuscleCommandline)
    ...

For the most basic usage, all you need is to have a FASTA input file,
such as
`opuntia.fasta <https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/opuntia.fasta>`__
(available online or in the Doc/examples subdirectory of the Biopython
source code). You can then tell MUSCLE to read in this FASTA file, and
write the alignment to an output file:

.. doctest

.. code:: pycon

    >>> from Bio.Align.Applications import MuscleCommandline
    >>> cline = MuscleCommandline(input="opuntia.fasta", out="opuntia.txt")
    >>> print(cline)
    muscle -in opuntia.fasta -out opuntia.txt

Note that MUSCLE uses “-in” and “-out” but in Biopython we have to use
“input” and “out” as the keyword arguments or property names. This is
because “in” is a reserved word in Python.

By default MUSCLE will output the alignment as a FASTA file (using
gapped sequences). The ``Bio.AlignIO`` module should be able to read
this alignment using ``format=fasta``. You can also ask for
ClustalW-like output:

.. doctest

.. code:: pycon

    >>> from Bio.Align.Applications import MuscleCommandline
    >>> cline = MuscleCommandline(input="opuntia.fasta", out="opuntia.aln", clw=True)
    >>> print(cline)
    muscle -in opuntia.fasta -out opuntia.aln -clw

Or, strict ClustalW output where the original ClustalW header line is
used for maximum compatibility:

.. doctest

.. code:: pycon

    >>> from Bio.Align.Applications import MuscleCommandline
    >>> cline = MuscleCommandline(input="opuntia.fasta", out="opuntia.aln", clwstrict=True)
    >>> print(cline)
    muscle -in opuntia.fasta -out opuntia.aln -clwstrict

The ``Bio.AlignIO`` module should be able to read these alignments using
``format=clustal``.

MUSCLE can also output in GCG MSF format (using the ``msf`` argument),
but Biopython can’t currently parse that, or using HTML which would give
a human readable web page (not suitable for parsing).

You can also set the other optional parameters, for example the maximum
number of iterations. See the built in help for details.

You would then run MUSCLE command line string as described above for
ClustalW, and parse the output using ``Bio.AlignIO`` to get an alignment
object.

MUSCLE using stdout
~~~~~~~~~~~~~~~~~~~

Using a MUSCLE command line as in the examples above will write the
alignment to a file. This means there will be no important information
written to the standard out (stdout) or standard error (stderr) handles.
However, by default MUSCLE will write the alignment to standard output
(stdout). We can take advantage of this to avoid having a temporary
output file! For example:

.. doctest

.. code:: pycon

    >>> from Bio.Align.Applications import MuscleCommandline
    >>> muscle_cline = MuscleCommandline(input="opuntia.fasta")
    >>> print(muscle_cline)
    muscle -in opuntia.fasta

If we run this via the wrapper, we get back the output as a string. In
order to parse this we can use ``StringIO`` to turn it into a handle.
Remember that MUSCLE defaults to using FASTA as the output format:

.. code:: pycon

    >>> from Bio.Align.Applications import MuscleCommandline
    >>> muscle_cline = MuscleCommandline(input="opuntia.fasta")
    >>> stdout, stderr = muscle_cline()
    >>> from StringIO import StringIO
    >>> from Bio import AlignIO
    >>> align = AlignIO.read(StringIO(stdout), "fasta")
    >>> print(align)
    SingleLetterAlphabet() alignment with 7 rows and 906 columns
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273289|gb|AF191663.1|AF191663
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273291|gb|AF191665.1|AF191665
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273290|gb|AF191664.1|AF191664
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273287|gb|AF191661.1|AF191661
    TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273286|gb|AF191660.1|AF191660
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273285|gb|AF191659.1|AF191659
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273284|gb|AF191658.1|AF191658

The above approach is fairly simple, but if you are dealing with very
large output text the fact that all of stdout and stderr is loaded into
memory as a string can be a potential drawback. Using the ``subprocess``
module we can work directly with handles instead:

.. code:: pycon

    >>> import subprocess
    >>> from Bio.Align.Applications import MuscleCommandline
    >>> muscle_cline = MuscleCommandline(input="opuntia.fasta")
    >>> child = subprocess.Popen(str(muscle_cline),
    ...                          stdout=subprocess.PIPE,
    ...                          stderr=subprocess.PIPE,
    ...                          universal_newlines=True,
    ...                          shell=(sys.platform!="win32"))
    >>> from Bio import AlignIO
    >>> align = AlignIO.read(child.stdout, "fasta")
    >>> print(align)
    SingleLetterAlphabet() alignment with 7 rows and 906 columns
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273289|gb|AF191663.1|AF191663
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273291|gb|AF191665.1|AF191665
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273290|gb|AF191664.1|AF191664
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273287|gb|AF191661.1|AF191661
    TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273286|gb|AF191660.1|AF191660
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273285|gb|AF191659.1|AF191659
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273284|gb|AF191658.1|AF191658

MUSCLE using stdin and stdout
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We don’t actually *need* to have our FASTA input sequences prepared in a
file, because by default MUSCLE will read in the input sequence from
standard input! Note this is a bit more advanced and fiddly, so don’t
bother with this technique unless you need to.

First, we’ll need some unaligned sequences in memory as ``SeqRecord``
objects. For this demonstration I’m going to use a filtered version of
the original FASTA file (using a generator expression), taking just six
of the seven sequences:

.. doctest

.. code:: pycon

    >>> from Bio import SeqIO
    >>> records = (r for r in SeqIO.parse("opuntia.fasta", "fasta") if len(r) < 900)

Then we create the MUSCLE command line, leaving the input and output to
their defaults (stdin and stdout). I’m also going to ask for strict
ClustalW format as for the output.

.. doctest

.. code:: pycon

    >>> from Bio.Align.Applications import MuscleCommandline
    >>> muscle_cline = MuscleCommandline(clwstrict=True)
    >>> print(muscle_cline)
    muscle -clwstrict

Now for the fiddly bits using the ``subprocess`` module, stdin and
stdout:

.. code:: pycon

    >>> import subprocess
    >>> import sys
    >>> child = subprocess.Popen(str(cline),
    ...                          stdin=subprocess.PIPE,
    ...                          stdout=subprocess.PIPE,
    ...                          stderr=subprocess.PIPE,
    ...                          universal_newlines=True,
    ...                          shell=(sys.platform!="win32"))

That should start MUSCLE, but it will be sitting waiting for its FASTA
input sequences, which we must supply via its stdin handle:

.. code:: pycon

    >>> SeqIO.write(records, child.stdin, "fasta")
    6
    >>> child.stdin.close()

After writing the six sequences to the handle, MUSCLE will still be
waiting to see if that is all the FASTA sequences or not – so we must
signal that this is all the input data by closing the handle. At that
point MUSCLE should start to run, and we can ask for the output:

.. code:: pycon

    >>> from Bio import AlignIO
    >>> align = AlignIO.read(child.stdout, "clustal")
    >>> print(align)
    SingleLetterAlphabet() alignment with 6 rows and 900 columns
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273290|gb|AF191664.1|AF19166
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273289|gb|AF191663.1|AF19166
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273287|gb|AF191661.1|AF19166
    TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273286|gb|AF191660.1|AF19166
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273285|gb|AF191659.1|AF19165
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273284|gb|AF191658.1|AF19165

Wow! There we are with a new alignment of just the six records, without
having created a temporary FASTA input file, or a temporary alignment
output file. However, a word of caution: Dealing with errors with this
style of calling external programs is much more complicated. It also
becomes far harder to diagnose problems, because you can’t try running
MUSCLE manually outside of Biopython (because you don’t have the input
file to supply). There can also be subtle cross platform issues (e.g.
Windows versus Linux, Python 2 versus Python 3), and how you run your
script can have an impact (e.g. at the command line, from IDLE or an
IDE, or as a GUI script). These are all generic Python issues though,
and not specific to Biopython.

If you find working directly with ``subprocess`` like this scary, there
is an alternative. If you execute the tool with ``muscle_cline()`` you
can supply any standard input as a big string,
``muscle_cline(stdin=...)``. So, provided your data isn’t very big, you
can prepare the FASTA input in memory as a string using ``StringIO``
(see Section [sec:appendix-handles]):

.. doctest

.. code:: pycon

    >>> from Bio import SeqIO
    >>> records = (r for r in SeqIO.parse("opuntia.fasta", "fasta") if len(r) < 900)
    >>> from StringIO import StringIO
    >>> handle = StringIO()
    >>> SeqIO.write(records, handle, "fasta")
    6
    >>> data = handle.getvalue()

You can then run the tool and parse the alignment as follows:

.. code:: pycon

    >>> stdout, stderr = muscle_cline(stdin=data)
    >>> from Bio import AlignIO
    >>> align = AlignIO.read(StringIO(stdout), "clustal")
    >>> print(align)
    SingleLetterAlphabet() alignment with 6 rows and 900 columns
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273290|gb|AF191664.1|AF19166
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273289|gb|AF191663.1|AF19166
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273287|gb|AF191661.1|AF19166
    TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273286|gb|AF191660.1|AF19166
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273285|gb|AF191659.1|AF19165
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273284|gb|AF191658.1|AF19165

You might find this easier, but it does require more memory (RAM) for
the strings used for the input FASTA and output Clustal formatted data.

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

Let’s start by creating a complete ``needle`` command line object in one
go:

.. doctest

.. code:: pycon

    >>> from Bio.Emboss.Applications import NeedleCommandline
    >>> needle_cline = NeedleCommandline(asequence="alpha.faa", bsequence="beta.faa",
    ...                                  gapopen=10, gapextend=0.5, outfile="needle.txt")
    >>> print(needle_cline)
    needle -outfile=needle.txt -asequence=alpha.faa -bsequence=beta.faa -gapopen=10 -gapextend=0.5

Why not try running this by hand at the command prompt? You should see
it does a pairwise comparison and records the output in the file
``needle.txt`` (in the default EMBOSS alignment file format).

Even if you have EMBOSS installed, running this command may not work –
you might get a message about “command not found” (especially on
Windows). This probably means that the EMBOSS tools are not on your PATH
environment variable. You can either update your PATH setting, or simply
tell Biopython the full path to the tool, for example:

.. doctest

.. code:: pycon

    >>> from Bio.Emboss.Applications import NeedleCommandline
    >>> needle_cline = NeedleCommandline(r"C:\EMBOSS\needle.exe",
    ...                                  asequence="alpha.faa", bsequence="beta.faa",
    ...                                  gapopen=10, gapextend=0.5, outfile="needle.txt")

Remember in Python that for a default string ``\n`` or ``\t`` means a
new line or a tab – which is why we’re put a letter “r” at the start for
a raw string.

At this point it might help to try running the EMBOSS tools yourself by
hand at the command line, to familiarise yourself the other options and
compare them to the Biopython help text:

.. code:: pycon

    >>> from Bio.Emboss.Applications import NeedleCommandline
    >>> help(NeedleCommandline)
    ...

Note that you can also specify (or change or look at) the settings like
this:

.. doctest

.. code:: pycon

    >>> from Bio.Emboss.Applications import NeedleCommandline
    >>> needle_cline = NeedleCommandline()
    >>> needle_cline.asequence="alpha.faa"
    >>> needle_cline.bsequence="beta.faa"
    >>> needle_cline.gapopen=10
    >>> needle_cline.gapextend=0.5
    >>> needle_cline.outfile="needle.txt"
    >>> print(needle_cline)
    needle -outfile=needle.txt -asequence=alpha.faa -bsequence=beta.faa -gapopen=10 -gapextend=0.5
    >>> print(needle_cline.outfile)
    needle.txt

Next we want to use Python to run this command for us. As explained
above, for full control, we recommend you use the built in Python
``subprocess`` module, but for simple usage the wrapper object usually
suffices:

.. code:: pycon

    >>> stdout, stderr = needle_cline()
    >>> print(stdout + stderr)
    Needleman-Wunsch global alignment of two sequences

Next we can load the output file with ``Bio.AlignIO`` as discussed
earlier in this chapter, as the ``emboss`` format:

.. code:: pycon

    >>> from Bio import AlignIO
    >>> align = AlignIO.read("needle.txt", "emboss")
    >>> print(align)
    SingleLetterAlphabet() alignment with 2 rows and 149 columns
    MV-LSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTY...KYR HBA_HUMAN
    MVHLTPEEKSAVTALWGKV--NVDEVGGEALGRLLVVYPWTQRF...KYH HBB_HUMAN

In this example, we told EMBOSS to write the output to a file, but you
*can* tell it to write the output to stdout instead (useful if you don’t
want a temporary output file to get rid of – use ``stdout=True`` rather
than the ``outfile`` argument), and also to read *one* of the one of the
inputs from stdin (e.g. ``asequence=stdin``, much like in the MUSCLE
example in the section above).

This has only scratched the surface of what you can do with ``needle``
and ``water``. One useful trick is that the second file can contain
multiple sequences (say five), and then EMBOSS will do five pairwise
alignments.

Pairwise sequence alignment
---------------------------

Pairwise sequence alignment is the process of aligning two sequences to
each other by optimizing the similarity score between them. Biopython
includes two built-in pairwise aligners: the ’old’ ``Bio.pairwise2``
module and the new ``PairwiseAligner`` class within the ``Bio.Align``
module (since Biopython version 1.72). Both can perform global and local
alignments and offer numerous options to change the alignment
parameters. Although ``pairwise2`` has gained some speed and memory
enhancements recently, the new ``PairwiseAligner`` is much faster; so if
you need to make many alignments with larger sequences, the latter would
be the tool to choose. ``pairwise2``, on the contrary, is also able to
align lists, which can be useful if your sequences do not consist of
single characters only.

Given that the parameters and sequences are the same, both aligners will
return the same alignments and alignment score (if the number of
alignments is too high they may return different subsets of all valid
alignments).

pairwise2
~~~~~~~~~

``Bio.pairwise2`` contains essentially the same algorithms as ``water``
(local) and ``needle`` (global) from the
`EMBOSS <http://emboss.sourceforge.net/>`__ suite (see above) and should
return the same results. The ``pairwise2`` module has undergone some
optimization regarding speed and memory consumption recently (Biopython
versions >1.67) so that for short sequences (global alignments: ~2000
residues, local alignments ~600 residues) it’s faster (or equally fast)
to use ``pairwise2`` than calling EMBOSS’ ``water`` or ``needle`` via
the command line tools.

Suppose you want to do a global pairwise alignment between the same two
hemoglobin sequences from above (``HBA_HUMAN``, ``HBB_HUMAN``) stored in
``alpha.faa`` and ``beta.faa``:

.. doctest examples

.. code:: pycon

    >>> from Bio import pairwise2
    >>> from Bio import SeqIO
    >>> seq1 = SeqIO.read("alpha.faa", "fasta")
    >>> seq2 = SeqIO.read("beta.faa", "fasta")
    >>> alignments = pairwise2.align.globalxx(seq1.seq, seq2.seq)

As you see, we call the alignment function with ``align.globalxx``. The
tricky part are the last two letters of the function name (here:
``xx``), which are used for decoding the scores and penalties for
matches (and mismatches) and gaps. The first letter decodes the match
score, e.g. ``x`` means that a match counts 1 while mismatches have no
costs. With ``m`` general values for either matches or mismatches can be
defined (for more options see `Biopython’s
API <http://biopython.org/DIST/docs/api/Bio.pairwise2-module.html>`__).
The second letter decodes the cost for gaps; ``x`` means no gap costs at
all, with ``s`` different penalties for opening and extending a gap can
be assigned. So, ``globalxx`` means that only matches between both
sequences are counted.

Our variable ``alignments`` now contains a list of alignments (at least
one) which have the same optimal score for the given conditions. In our
example this are 80 different alignments with the score 72
(``Bio.pairwise2`` will return up to 1000 alignments). Have a look at
one of these alignments:

.. cont-doctest

.. code:: pycon

    >>> len(alignments)
    80

.. code:: pycon

    >>> print(alignments[0])
    ('MV-LSPADKTNV---K-A--A-WGKVGAHAG...YR-', 'MVHL-----T--PEEKSAVTALWGKV----...Y-H',
    72.0, 0, 217)

Each alignment is a tuple consisting of the two aligned sequences, the
score, the start and the end positions of the alignment (in global
alignments the start is always 0 and the end the length of the
alignment). ``Bio.pairwise2`` has a function ``format_alignment`` for a
nicer printout:

.. code:: pycon

    >>> print(pairwise2.format_alignment(*alignment[0]))
    MV-LSPADKTNV---K-A--A-WGKVGAHAG---EY-GA-EALE-RMFLSF----PTTK-TY--F...YR-
    || |     |     | |  | ||||        |  |  |||  |  |      |    |   |...|  
    MVHL-----T--PEEKSAVTALWGKV-----NVDE-VG-GEAL-GR--L--LVVYP---WT-QRF...Y-H
      Score=72

Better alignments are usually obtained by penalizing gaps: higher costs
for opening a gap and lower costs for extending an existing gap. For
amino acid sequences match scores are usually encoded in matrices like
``PAM`` or ``BLOSUM``. Thus, a more meaningful alignment for our example
can be obtained by using the BLOSUM62 matrix, together with a gap open
penalty of 10 and a gap extension penalty of 0.5 (using ``globalds``):

.. doctest examples

.. code:: pycon

    >>> from Bio import pairwise2
    >>> from Bio import SeqIO
    >>> from Bio.SubsMat.MatrixInfo import blosum62
    >>> seq1 = SeqIO.read("alpha.faa", "fasta")
    >>> seq2 = SeqIO.read("beta.faa", "fasta")
    >>> alignments = pairwise2.align.globalds(seq1.seq, seq2.seq, blosum62, -10, -0.5)
    >>> len(alignments)
    2

.. code:: pycon

    >>> print(pairwise2.format_alignment(*alignments[0]))
    MV-LSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTY...KYR
    || |.|..|..|.|.|||| ......|............|.......||.
    MVHLTPEEKSAVTALWGKV-NVDEVGGEALGRLLVVYPWTQRFF...KYH
      Score=292.5

This alignment has the same score that we obtained earlier with EMBOSS
needle using the same sequences and the same parameters.

Local alignments are called similarly with the function
``align.localXX``, where again XX stands for a two letter code for the
match and gap functions:

.. doctest

.. code:: pycon

    >>> from Bio import pairwise2
    >>> from Bio.SubsMat.MatrixInfo import blosum62
    >>> alignments = pairwise2.align.localds("LSPADKTNVKAA", "PEEKSAV", blosum62, -10, -1)
    >>> print(pairwise2.format_alignment(*alignments[0]))
    3 PADKTNV
      |..|..|
    1 PEEKSAV
      Score=16
    <BLANKLINE>

Note that local alignments must, as defined by Smith & Waterman, have a
positive score (>0). Thus, ``pairwise2`` may return no alignments if no
score >0 has been obtained. Also, ``pairwise2`` will not report
alignments which are the result of the addition of zero-scoring
extensions on either site. In the next example, the pairs
serin/aspartate (S/D) and lysin/asparagin (K/N) both have a match score
of 0. As you see, the aligned part has not been extended:

.. doctest

.. code:: pycon

    >>> from Bio import pairwise2
    >>> from Bio.SubsMat.MatrixInfo import blosum62
    >>> alignments = pairwise2.align.localds("LSSPADKTNVKKAA", "DDPEEKSAVNN", blosum62, -10, -1)
    >>> print(pairwise2.format_alignment(*alignments[0]))
    4 PADKTNV
      |..|..|
    3 PEEKSAV
      Score=16
    <BLANKLINE>

Instead of supplying a complete match/mismatch matrix, the match code
``m`` allows for easy defining general match/mismatch values. The next
example uses match/mismatch scores of 5/-4 and gap penalties
(open/extend) of 2/0.5 using ``localms``:

.. cont-doctest

.. code:: pycon

    >>> alignments = pairwise2.align.localms("AGAACT", "GAC", 5, -4, -2, -0.5)
    >>> print(pairwise2.format_alignment(*alignments[0]))
    2 GAAC
      | ||
    1 G-AC
      Score=13
    <BLANKLINE>

One useful keyword argument of the ``Bio.pairwise2.align`` functions is
``score_only``. When set to ``True`` it will only return the score of
the best alignment(s), but in a significantly shorter time. It will also
allow the alignment of longer sequences before a memory error is raised.
Another useful keyword argument is ``one_alignment_only=True`` which
will also result in some speed gain.

Unfortunately, ``Bio.pairwise2`` does not work with Biopython’s multiple
sequence alignment objects (yet). However, the module has some
interesting advanced features: you can define your own match and gap
functions (interested in testing affine logarithmic gap costs?), gap
penalties and end gaps penalties can be different for both sequences,
sequences can be supplied as lists (useful if you have residues that are
encoded by more than one character), etc. These features are hard (if at
all) to realize with other alignment tools. For more details see the
modules documentation in `Biopython’s
API <http://biopython.org/DIST/docs/api/Bio.pairwise2-module.html>`__.

PairwiseAligner
~~~~~~~~~~~~~~~

The new ``Bio.Align.PairwiseAligner`` implements the Needleman-Wunsch,
Smith-Waterman, Gotoh (three-state), and Waterman-Smith-Beyer global and
local pairwise alignment algorithms. We refer to Durbin *et al.*
:raw-latex:`\cite{durbin1998}` for in-depth information on sequence
alignment algorithms.

Basic usage
^^^^^^^^^^^

To generate pairwise alignments, first create a ``PairwiseAligner``
object:

.. doctest examples

.. code:: pycon

    >>> from Bio import Align
    >>> aligner = Align.PairwiseAligner()

The ``PairwiseAligner`` object ``aligner`` (see
Section [sec:pairwise-aligner]) stores the alignment parameters to be
used for the pairwise alignments.

Use the ``aligner.score`` method to calculate the alignment score
between two sequences:

.. cont-doctest

.. code:: pycon

    >>> seq1 = "GAACT"
    >>> seq2 = "GAT"
    >>> score = aligner.score(seq1, seq2)
    >>> score
    3.0

To see the actual alignments, use the ``aligner.align`` method and
iterate over the ``PairwiseAlignment`` objects returned:

.. cont-doctest

.. code:: pycon

    >>> alignments = aligner.align(seq1, seq2)
    >>> for alignment in alignments:
    ...     print(alignment)
    ... 
    GAACT
    ||--|
    GA--T
    <BLANKLINE>
    GAACT
    |-|-|
    G-A-T
    <BLANKLINE>

By default, a global pairwise alignment is performed, which finds the
optimal alignment over the whole length of ``seq1`` and ``seq2``.
Instead, a local alignment will find the subsequence of ``seq1`` and
``seq2`` with the highest alignment score. Local alignments can be
generated by setting ``aligner.mode`` to ``"local"``:

.. cont-doctest

.. code:: pycon

    >>> aligner.mode = 'local'
    >>> seq1 = "AGAACTC"
    >>> seq2 = "GAACT"
    >>> score = aligner.score(seq1, seq2)
    >>> score
    5.0
    >>> alignments = aligner.align(seq1, seq2)
    >>> for alignment in alignments:
    ...     print(alignment)
    ...
    AGAACTC
    .|||||.
    .GAACT.
    <BLANKLINE>

Note that there is some ambiguity in the definition of the best local
alignments if segments with a score 0 can be added to the alignment. We
follow the suggestion by Waterman & Eggert
:raw-latex:`\cite{waterman1987}` and disallow such extensions.

The pairwise aligner object
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``PairwiseAligner`` object stores all alignment parameters to be
used for the pairwise alignments. To see an overview of the values for
all parameters, use

.. cont-doctest

.. code:: pycon

    >>> print(aligner)
    Pairwise sequence aligner with parameters
      match_score: 1.000000
      mismatch_score: 0.000000
      target_open_gap_score: 0.000000
      target_extend_gap_score: 0.000000
      target_left_open_gap_score: 0.000000
      target_left_extend_gap_score: 0.000000
      target_right_open_gap_score: 0.000000
      target_right_extend_gap_score: 0.000000
      query_open_gap_score: 0.000000
      query_extend_gap_score: 0.000000
      query_left_open_gap_score: 0.000000
      query_left_extend_gap_score: 0.000000
      query_right_open_gap_score: 0.000000
      query_right_extend_gap_score: 0.000000
      mode: local
    <BLANKLINE>

See Sections [sec:pairwise-matchscores],
[sec:pairwise-affine-gapscores], and [sec:pairwise-general-gapscores]
below for the definition of these parameters. The attribute ``mode``
(described above in Section [sec:pairwise-basic]) can be set equal to
``"global"`` or ``"local"`` to specify global or local pairwise
alignment, respectively.

Depending on the gap scoring parameters (see
Sections [sec:pairwise-affine-gapscores] and
[sec:pairwise-general-gapscores]) and mode, a ``PairwiseAligner`` object
automatically chooses the appropriate algorithm to use for pairwise
sequence alignment. To verify the selected algorithm, use

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

Match and mismatch scores
^^^^^^^^^^^^^^^^^^^^^^^^^

The match and mismatch scores are stored as attributes of an
``PairwiseAligner`` object:

.. doctest examples

.. code:: pycon

    >>> from Bio import Align
    >>> aligner = Align.PairwiseAligner()
    >>> aligner.match_score
    1.0
    >>> aligner.mismatch_score
    0.0
    >>> score = aligner.score("AAA","AAC")
    >>> print(score)
    2.0
    >>> aligner.match_score = 2.0
    >>> score = aligner.score("AAA","AAC")
    >>> print(score)
    4.0

Alternatively, you can specify a substitution matrix as follows:

.. cont-doctest

.. code:: pycon

    >>> matrix = {('A','A'): 1.0, ('A','B'): -1.0, ('B','B'): 2.0}
    >>> aligner.substitution_matrix = matrix

The attributes ``aligner.match_score`` and ``aligner.mismatch_score``
are ignored if ``aligner.substitution_matrix`` is specified. Likewise,
after specifying ``aligner.match_score`` or ``aligner.mismatch_score``,
``aligner.substitution_matrix`` will be ignored.

Note that ``aligner.substitution_matrix`` will return a copy of the
substitution matrix stored in the ``PairwiseAligner`` object. Therefore,
modifying the substitution matrix directly has no effect:

.. cont-doctest

.. code:: pycon

    >>> aligner.substitution_matrix[('A','A')] = 5.0 # does nothing
    >>> aligner.substitution_matrix[('A','A')]
    1.0

In all cases, the character ``X`` is used to denote unknown characters,
which will always get a zero score in alignments, irrespective of the
match or mismatch score.

Affine gap scores
^^^^^^^^^^^^^^^^^

Affine gap scores are defined by a score to open a gap, and a score to
extend an existing gap:

:math:`\textrm{gap score} = \textrm{open gap score} + (n-1) \times \textrm{extend gap score}`,

where :math:`n` is the length of the gap. Biopython’s pairwise sequence
aligner allows fine-grained control over the gap scoring scheme by
specifying the following twelve attributes of a ``PairwiseAligner``
object:

+-----------------------------------------------------------------------------+----+
| **Opening scores & **Extending scores                                       |    |
| ``query_left_open_gap_score`` & ``query_left_extend_gap_score``             |    |
| ``query_internal_open_gap_score`` & ``query_internal_extend_gap_score``     |    |
| ``query_right_open_gap_score`` & ``query_right_extend_gap_score``           |    |
| ``target_left_open_gap_score`` & ``target_left_extend_gap_score``           |    |
| ``target_internal_open_gap_score`` & ``target_internal_extend_gap_score``   |    |
| ``target_right_open_gap_score`` & ``target_right_extend_gap_score``         |    |
| ****                                                                        |    |
+-----------------------------------------------------------------------------+----+

These attributes allow for different gap scores for internal gaps and on
either end of the sequence, as shown in this example:

+--------------------------------------------+----+----+
| **target & **query & **score               |    |    |
| A & - & query left open gap score          |    |    |
| C & - & query left extend gap score        |    |    |
| C & - & query left extend gap score        |    |    |
| G & G & match score                        |    |    |
| G & T & mismatch score                     |    |    |
| G & - & query internal open gap score      |    |    |
| A & - & query internal extend gap score    |    |    |
| A & - & query internal extend gap score    |    |    |
| T & T & match score                        |    |    |
| A & A & match score                        |    |    |
| G & - & query internal open gap score      |    |    |
| C & C & match score                        |    |    |
| - & C & target internal open gap score     |    |    |
| - & C & target internal extend gap score   |    |    |
| C & C & match score                        |    |    |
| T & G & mismatch score                     |    |    |
| C & C & match score                        |    |    |
| - & C & target internal open gap score     |    |    |
| A & A & match score                        |    |    |
| - & T & target right open gap score        |    |    |
| - & A & target right extend gap score      |    |    |
| - & A & target right extend gap score      |    |    |
| ******                                     |    |    |
+--------------------------------------------+----+----+

For convenience, ``PairwiseAligner`` objects have additional attributes
that refer to a number of these values collectively, as shown
(hierarchically) in Table [table:align-meta-attributes].

+-------------------------------------------------------------------------------------------------------------+----+
| **Meta-attribute & **Attributes it maps to                                                                  |    |
| ``gap_score`` & ``target_gap_score``, ``query_gap_score``                                                   |    |
| ``open_gap_score`` & ``target_open_gap_score``, ``query_open_gap_score``                                    |    |
| ``extend_gap_score`` & ``target_extend_gap_score``, ``query_extend_gap_score``                              |    |
| ``internal_gap_score`` & ``target_internal_gap_score``, ``query_internal_gap_score``                        |    |
| ``internal_open_gap_score`` & ``target_internal_open_gap_score``, ``query_internal_open_gap_score``         |    |
| ``internal_extend_gap_score`` & ``target_internal_extend_gap_score``, ``query_internal_extend_gap_score``   |    |
| ``end_gap_score`` & ``target_end_gap_score``, ``query_end_gap_score``                                       |    |
| ``end_open_gap_score`` & ``target_end_open_gap_score``, ``query_end_open_gap_score``                        |    |
| ``end_extend_gap_score`` & ``target_end_extend_gap_score``, ``query_end_extend_gap_score``                  |    |
| ``left_gap_score`` & ``target_left_gap_score``, ``query_left_gap_score``                                    |    |
| ``right_gap_score`` & ``target_right_gap_score``, ``query_right_gap_score``                                 |    |
| ``left_open_gap_score`` & ``target_left_open_gap_score``, ``query_left_open_gap_score``                     |    |
| ``left_extend_gap_score`` & ``target_left_extend_gap_score``, ``query_left_extend_gap_score``               |    |
| ``right_open_gap_score`` & ``target_right_open_gap_score``, ``query_right_open_gap_score``                  |    |
| ``right_extend_gap_score`` & ``target_right_extend_gap_score``, ``query_right_extend_gap_score``            |    |
| ``target_open_gap_score`` & ``target_internal_open_gap_score``, ``target_left_open_gap_score``,             |    |
| & ``target_right_open_gap_score``                                                                           |    |
| ``target_extend_gap_score`` & ``target_internal_extend_gap_score``, ``target_left_extend_gap_score``,       |    |
| & ``target_right_extend_gap_score``                                                                         |    |
| ``target_gap_score`` & ``target_open_gap_score``, ``target_extend_gap_score``                               |    |
| ``query_open_gap_score`` & ``query_internal_open_gap_score``, ``query_left_open_gap_score``,                |    |
| & ``query_right_open_gap_score``                                                                            |    |
| ``query_extend_gap_score`` & ``query_internal_extend_gap_score``, ``query_left_extend_gap_score``,          |    |
| & ``query_right_extend_gap_score``                                                                          |    |
| ``query_gap_score`` & ``query_open_gap_score``, ``query_extend_gap_score``                                  |    |
| ``target_internal_gap_score`` & ``target_internal_open_gap_score``, ``target_internal_extend_gap_score``    |    |
| ``target_end_gap_score`` & ``target_end_open_gap_score``, ``target_end_extend_gap_score``                   |    |
| ``target_end_open_gap_score`` & ``target_left_open_gap_score``, ``target_right_open_gap_score``             |    |
| ``target_end_extend_gap_score`` & ``target_left_extend_gap_score``, ``target_right_extend_gap_score``       |    |
| ``target_left_gap_score`` & ``target_left_open_gap_score``, ``target_left_extend_gap_score``                |    |
| ``target_right_gap_score`` & ``target_right_open_gap_score``, ``target_right_extend_gap_score``             |    |
| ``query_end_gap_score`` & ``query_end_open_gap_score``, ``query_end_extend_gap_score``                      |    |
| ``query_end_open_gap_score`` & ``query_left_open_gap_score``, ``query_right_open_gap_score``                |    |
| ``query_end_extend_gap_score`` & ``query_left_extend_gap_score``, ``query_right_extend_gap_score``          |    |
| ``query_internal_gap_score`` & ``query_internal_open_gap_score``, ``query_internal_extend_gap_score``       |    |
| ``query_left_gap_score`` & ``query_left_open_gap_score``, ``query_left_extend_gap_score``                   |    |
| ``query_right_gap_score`` & ``query_right_open_gap_score``, ``query_right_extend_gap_score``                |    |
| ****                                                                                                        |    |
+-------------------------------------------------------------------------------------------------------------+----+

Table: Meta-attributes of the pairwise aligner objects.

[table:align-meta-attributes]

General gap scores
^^^^^^^^^^^^^^^^^^

For even more fine-grained control over the gap scores, you can specify
a gap scoring function. For example, the gap scoring function below
disallows a gap after two nucleotides in the query sequence:

.. doctest

.. code:: pycon

    >>> from Bio import Align
    >>> aligner = Align.PairwiseAligner()
    >>> def my_gap_score_function(start, length):
    ...     if start==2:
    ...         return -1000
    ...     else:
    ...         return -1 * length
    ...
    >>> aligner.query_gap_score = my_gap_score_function
    >>> alignments = aligner.align("AACTT", "AATT")
    >>> for alignment in alignments:
    ...     print(alignment)
    ... 
    AACTT
    -|X||
    -AATT
    <BLANKLINE>
    AACTT
    |-X||
    A-ATT
    <BLANKLINE>
    AACTT
    ||X-|
    AAT-T
    <BLANKLINE>
    AACTT
    ||X|-
    AATT-
    <BLANKLINE>

Iterating over alignments
^^^^^^^^^^^^^^^^^^^^^^^^^

The ``alignments`` returned by ``aligner.align`` are a kind of immutable
iterable objects (similar to ``range``). While they appear similarto a
``tuple`` or ``list`` of ``PairwiseAlignment`` objects, they are
different in the sense that each ``PairwiseAlignment`` object is created
dynamically when it is needed. This approach was chosen because the
number of alignments can be extremely large, in particular for poor
alignments (see Section [sec:pairwise-examples] for an example).

You can perform the following operations on ``alignments``:

-  ``len(alignments)`` returns the number of alignments stored. This
   function returns quickly, even if the number of alignments is huge.
   If the number of alignments is extremely large (typically, larger
   than 9,223,372,036,854,775,807, which is the largest integer that can
   be stored as a ``long int`` on 64 bit machines), ``len(alignments)``
   will raise an ``OverflowError``. A large number of alignments
   suggests that the alignment quality is low.

   .. doctest examples

   .. code:: pycon

       >>> from Bio import Align
       >>> aligner = Align.PairwiseAligner()
       >>> alignments = aligner.align("AAA", "AA")
       >>> len(alignments)
       3

-  You can extract a specific alignment by index:

   .. doctest examples

   .. code:: pycon

       >>> from Bio import Align
       >>> aligner = Align.PairwiseAligner()
       >>> alignments = aligner.align("AAA", "AA")
       >>> print(alignments[2])
       AAA
       -||
       -AA
       <BLANKLINE>
       >>> print(alignments[0])
       AAA
       ||-
       AA-
       <BLANKLINE>

-  You can iterate over alignments, for example as in

   .. code:: pycon

       >>> for alignment in alignments:
       ...     print(alignment)
       ...

   Note that ``alignments`` can be reused, i.e. you can iterate over
   alignments multiple times:

   .. doctest examples

   .. code:: pycon

       >>> from Bio import Align
       >>> aligner = Align.PairwiseAligner()
       >>> alignments = aligner.align("AAA", "AA")
       >>> for alignment in alignments:
       ...     print(alignment)
       ...
       AAA
       ||-
       AA-
       <BLANKLINE>
       AAA
       |-|
       A-A
       <BLANKLINE>
       AAA
       -||
       -AA
       <BLANKLINE>
       >>> for alignment in alignments:
       ...     print(alignment)
       ...
       AAA
       ||-
       AA-
       <BLANKLINE>
       AAA
       |-|
       A-A
       <BLANKLINE>
       AAA
       -||
       -AA
       <BLANKLINE>

   You can also convert the ``alignments`` iterator into a ``list`` or
   ``tuple``:

   .. code:: pycon

       >>> alignments = list(alignments)

   It is wise to check the number of alignments by calling
   ``len(alignments)`` before attempting to call ``list(alignments)`` to
   save all alignments as a list.

-  The alignment score (which has the same value for each alignment in
   ``alignments``), is stored as an attribute. This allows you to check
   the alignment score before proceeding to extract individual
   alignments:

   .. cont-doctest

   .. code:: pycon

       >>> print(alignments.score)
       2.0

Alignment objects
^^^^^^^^^^^^^^^^^

The ``aligner.align`` method returns ``PairwiseAlignment`` objects, each
representing one alignment between the two sequences.

.. doctest

.. code:: pycon

    >>> from Bio import Align
    >>> aligner = Align.PairwiseAligner()
    >>> seq1 = "GAACT"
    >>> seq2 = "GAT"
    >>> alignments = aligner.align(seq1, seq2)
    >>> alignment = alignments[0]

.. code:: pycon

    >>> alignment
    <Bio.Align.PairwiseAlignment object at 0x10204d250>

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

Print the ``PairwiseAlignment`` object to show the alignment explicitly:

.. cont-doctest

.. code:: pycon

    >>> print(alignment)
    GAACT
    ||--|
    GA--T
    <BLANKLINE>

You can also represent the alignment as a string in PSL (Pattern Space
Layout, as generated by BLAT :raw-latex:`\cite{kent2002}`) format:

.. cont-doctest

.. code:: pycon

    >>> format(alignment, 'psl')
    '3\t0\t0\t0\t0\t0\t1\t2\t+\tquery\t3\t0\t3\ttarget\t5\t0\t5\t2\t2,1,\t0,2,\t0,4,\n'

Use the ``aligned`` property to find the start and end indices of
subsequences in the target and query sequence that were aligned to each
other. Generally, if the alignment between target (t) and query (q)
consists of N chunks, you get two tuples of length N:

.. code:: python

    (((t_start1, t_end1), (t_start2, t_end2), ..., (t_startN, t_endN)),
     ((q_start1, q_end1), (q_start2, q_end2), ..., (q_startN, q_endN)))

In the current example, ‘alignment.aligned‘ returns two tuples of length
2:

.. cont-doctest

.. code:: pycon

    >>> alignment.aligned
    (((0, 2), (4, 5)), ((0, 2), (2, 3)))

while for the alternative alignment, two tuples of length 3 are
returned:

.. cont-doctest

.. code:: pycon

    >>> alignment = alignments[1]
    >>> print(alignment)
    GAACT
    |-|-|
    G-A-T
    <BLANKLINE>
    >>> alignment.aligned
    (((0, 1), (2, 3), (4, 5)), ((0, 1), (1, 2), (2, 3)))

Note that different alignments may have the same subsequences aligned to
each other. In particular, this may occur if alignments differ from each
other in terms of their gap placement only:

.. cont-doctest

.. code:: pycon

    >>> aligner.mismatch_score = -10
    >>> alignments = aligner.align("AAACAAA", "AAAGAAA")
    >>> len(alignments)
    2
    >>> print(alignments[0])
    AAAC-AAA
    |||--|||
    AAA-GAAA
    <BLANKLINE>
    >>> alignments[0].aligned
    (((0, 3), (4, 7)), ((0, 3), (4, 7)))
    >>> print(alignments[1])
    AAA-CAAA
    |||--|||
    AAAG-AAA
    <BLANKLINE>
    >>> alignments[1].aligned
    (((0, 3), (4, 7)), ((0, 3), (4, 7)))

The ``aligned`` property can be used to identify alignments that are
identical to each other in terms of their aligned sequences.

Example
^^^^^^^

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

.. code:: pycon

    >>> print(alignment)
    MV-LS-PAD--KTN--VK-AA-WGKV-----GAHAGEYGAEALE-RMFLSF----P-TTKTY--FPHF--...
    ||-|--|----|----|--|--||||-----|---||--|--|--|--|------|-|------|--|--...
    MVHL-TP--EEK--SAV-TA-LWGKVNVDEVG---GE--A--L-GR--L--LVVYPWT----QRF--FES...

Better alignments are usually obtained by penalizing gaps: higher costs
for opening a gap and lower costs for extending an existing gap. For
amino acid sequences match scores are usually encoded in matrices like
``PAM`` or ``BLOSUM``. Thus, a more meaningful alignment for our example
can be obtained by using the BLOSUM62 matrix, together with a gap open
penalty of 10 and a gap extension penalty of 0.5:

.. doctest examples

.. code:: pycon

    >>> from Bio import Align
    >>> from Bio import SeqIO
    >>> from Bio.SubsMat.MatrixInfo import blosum62
    >>> seq1 = SeqIO.read("alpha.faa", "fasta")
    >>> seq2 = SeqIO.read("beta.faa", "fasta")
    >>> aligner = Align.PairwiseAligner()
    >>> aligner.open_gap_score = -10
    >>> aligner.extend_gap_score = -0.5
    >>> aligner.substitution_matrix = blosum62
    >>> score = aligner.score(seq1.seq, seq2.seq)
    >>> print(score)
    292.5
    >>> alignments = aligner.align(seq1.seq, seq2.seq)
    >>> len(alignments)
    2
    >>> print(alignments[0].score)
    292.5

.. code:: pycon

    >>> print(alignments[0])
    MV-LSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF-DLS-----HGSAQVKGHGKKV...
    ||-|X|XX|XX|X|X||||--XXX|X|X|||X|XXXXX|X|XXX|XX|-|||---- X|XXX||X|||||...
    MVHLTPEEKSAVTALWGKV--NVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKV...

This alignment has the same score that we obtained earlier with EMBOSS
needle using the same sequences and the same parameters.

To perform a local alignment, set ``aligner.mode`` to ``'local'``:

.. cont-doctest

.. code:: pycon

    >>> aligner.mode = 'local'
    >>> aligner.open_gap_score = -10
    >>> aligner.extend_gap_score = -1
    >>> alignments = aligner.align("LSPADKTNVKAA", "PEEKSAV")
    >>> print(len(alignments))
    1
    >>> alignment = alignments[0]
    >>> print(alignment)
    LSPADKTNVKAA
    ..|XX|XX|...
    ..PEEKSAV...
    <BLANKLINE>
    >>> print(alignment.score)
    16.0
