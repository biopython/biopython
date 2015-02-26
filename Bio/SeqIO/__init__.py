# Copyright 2006-2010 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
# Nice link:
# http://www.ebi.ac.uk/help/formats_frame.html

r"""Sequence input/output as SeqRecord objects.

Bio.SeqIO is also documented at SeqIO_ and by
a whole chapter in our tutorial:

  - `HTML Tutorial`_
  - `PDF Tutorial`_

.. _SeqIO: http://biopython.org/wiki/SeqIO
.. _`HTML Tutorial`: http://biopython.org/DIST/docs/tutorial/Tutorial.html
.. _`PDF Tutorial`: http://biopython.org/DIST/docs/tutorial/Tutorial.pdf

Input
-----
The main function is Bio.SeqIO.parse(...) which takes an input file handle
(or in recent versions of Biopython alternatively a filename as a string),
and format string.  This returns an iterator giving SeqRecord objects:

>>> from Bio import SeqIO
>>> for record in SeqIO.parse("Fasta/f002", "fasta"):
...     print("%s %i" % (record.id, len(record)))
gi|1348912|gb|G26680|G26680 633
gi|1348917|gb|G26685|G26685 413
gi|1592936|gb|G29385|G29385 471

Note that the parse() function will invoke the relevant parser for the
format with its default settings.  You may want more control, in which case
you need to create a format specific sequence iterator directly.

Input - Single Records
----------------------
If you expect your file to contain one-and-only-one record, then we provide
the following 'helper' function which will return a single SeqRecord, or
raise an exception if there are no records or more than one record:

>>> from Bio import SeqIO
>>> record = SeqIO.read("Fasta/f001", "fasta")
>>> print("%s %i" % (record.id, len(record)))
gi|3318709|pdb|1A91| 79

This style is useful when you expect a single record only (and would
consider multiple records an error).  For example, when dealing with GenBank
files for bacterial genomes or chromosomes, there is normally only a single
record.  Alternatively, use this with a handle when downloading a single
record from the internet.

However, if you just want the first record from a file containing multiple
record, use the next() function on the iterator (or under Python 2, the
iterator's next() method):

>>> from Bio import SeqIO
>>> record = next(SeqIO.parse("Fasta/f002", "fasta"))
>>> print("%s %i" % (record.id, len(record)))
gi|1348912|gb|G26680|G26680 633

The above code will work as long as the file contains at least one record.
Note that if there is more than one record, the remaining records will be
silently ignored.


Input - Multiple Records
------------------------
For non-interlaced files (e.g. Fasta, GenBank, EMBL) with multiple records
using a sequence iterator can save you a lot of memory (RAM).  There is
less benefit for interlaced file formats (e.g. most multiple alignment file
formats).  However, an iterator only lets you access the records one by one.

If you want random access to the records by number, turn this into a list:

>>> from Bio import SeqIO
>>> records = list(SeqIO.parse("Fasta/f002", "fasta"))
>>> len(records)
3
>>> print(records[1].id)
gi|1348917|gb|G26685|G26685

If you want random access to the records by a key such as the record id,
turn the iterator into a dictionary:

>>> from Bio import SeqIO
>>> record_dict = SeqIO.to_dict(SeqIO.parse("Fasta/f002", "fasta"))
>>> len(record_dict)
3
>>> print(len(record_dict["gi|1348917|gb|G26685|G26685"]))
413

However, using list() or the to_dict() function will load all the records
into memory at once, and therefore is not possible on very large files.
Instead, for *some* file formats Bio.SeqIO provides an indexing approach
providing dictionary like access to any record. For example,

>>> from Bio import SeqIO
>>> record_dict = SeqIO.index("Fasta/f002", "fasta")
>>> len(record_dict)
3
>>> print(len(record_dict["gi|1348917|gb|G26685|G26685"]))
413
>>> record_dict.close()

Many but not all of the supported input file formats can be indexed like
this. For example "fasta", "fastq", "qual" and even the binary format "sff"
work, but alignment formats like "phylip", "clustalw" and "nexus" will not.

In most cases you can also use SeqIO.index to get the record from the file
as a raw string (not a SeqRecord). This can be useful for example to extract
a sub-set of records from a file where SeqIO cannot output the file format
(e.g. the plain text SwissProt format, "swiss") or where it is important to
keep the output 100% identical to the input). For example,

>>> from Bio import SeqIO
>>> record_dict = SeqIO.index("Fasta/f002", "fasta")
>>> len(record_dict)
3
>>> print(record_dict.get_raw("gi|1348917|gb|G26685|G26685").decode())
>gi|1348917|gb|G26685|G26685 human STS STS_D11734.
CGGAGCCAGCGAGCATATGCTGCATGAGGACCTTTCTATCTTACATTATGGCTGGGAATCTTACTCTTTC
ATCTGATACCTTGTTCAGATTTCAAAATAGTTGTAGCCTTATCCTGGTTTTACAGATGTGAAACTTTCAA
GAGATTTACTGACTTTCCTAGAATAGTTTCTCTACTGGAAACCTGATGCTTTTATAAGCCATTGTGATTA
GGATGACTGTTACAGGCTTAGCTTTGTGTGAAANCCAGTCACCTTTCTCCTAGGTAATGAGTAGTGCTGT
TCATATTACTNTAAGTTCTATAGCATACTTGCNATCCTTTANCCATGCTTATCATANGTACCATTTGAGG
AATTGNTTTGCCCTTTTGGGTTTNTTNTTGGTAAANNNTTCCCGGGTGGGGGNGGTNNNGAAA
<BLANKLINE>
>>> print(record_dict["gi|1348917|gb|G26685|G26685"].format("fasta"))
>gi|1348917|gb|G26685|G26685 human STS STS_D11734.
CGGAGCCAGCGAGCATATGCTGCATGAGGACCTTTCTATCTTACATTATGGCTGGGAATC
TTACTCTTTCATCTGATACCTTGTTCAGATTTCAAAATAGTTGTAGCCTTATCCTGGTTT
TACAGATGTGAAACTTTCAAGAGATTTACTGACTTTCCTAGAATAGTTTCTCTACTGGAA
ACCTGATGCTTTTATAAGCCATTGTGATTAGGATGACTGTTACAGGCTTAGCTTTGTGTG
AAANCCAGTCACCTTTCTCCTAGGTAATGAGTAGTGCTGTTCATATTACTNTAAGTTCTA
TAGCATACTTGCNATCCTTTANCCATGCTTATCATANGTACCATTTGAGGAATTGNTTTG
CCCTTTTGGGTTTNTTNTTGGTAAANNNTTCCCGGGTGGGGGNGGTNNNGAAA
<BLANKLINE>
>>> record_dict.close()

Here the original file and what Biopython would output differ in the line
wrapping. Also note that under Python 3, the get_raw method will return a
bytes string, hence the use of decode to turn it into a (unicode) string.
This is uncessary on Python 2.

Also note that the get_raw method will preserve the newline endings. This
example FASTQ file uses Unix style endings (b"\n" only),

>>> from Bio import SeqIO
>>> fastq_dict = SeqIO.index("Quality/example.fastq", "fastq")
>>> len(fastq_dict)
3
>>> raw = fastq_dict.get_raw("EAS54_6_R1_2_1_540_792")
>>> raw.count(b"\n")
4
>>> raw.count(b"\r\n")
0
>>> b"\r" in raw
False
>>> len(raw)
78
>>> fastq_dict.close()

Here is the same file but using DOS/Windows new lines (b"\r\n" instead),

>>> from Bio import SeqIO
>>> fastq_dict = SeqIO.index("Quality/example_dos.fastq", "fastq")
>>> len(fastq_dict)
3
>>> raw = fastq_dict.get_raw("EAS54_6_R1_2_1_540_792")
>>> raw.count(b"\n")
4
>>> raw.count(b"\r\n")
4
>>> b"\r\n" in raw
True
>>> len(raw)
82
>>> fastq_dict.close()

Because this uses two bytes for each new line, the file is longer than
the Unix equivalent with only one byte.


Input - Alignments
------------------
You can read in alignment files as alignment objects using Bio.AlignIO.
Alternatively, reading in an alignment file format via Bio.SeqIO will give
you a SeqRecord for each row of each alignment:

>>> from Bio import SeqIO
>>> for record in SeqIO.parse("Clustalw/hedgehog.aln", "clustal"):
...     print("%s %i" % (record.id, len(record)))
gi|167877390|gb|EDS40773.1| 447
gi|167234445|ref|NP_001107837. 447
gi|74100009|gb|AAZ99217.1| 447
gi|13990994|dbj|BAA33523.2| 447
gi|56122354|gb|AAV74328.1| 447


Output
------
Use the function Bio.SeqIO.write(...), which takes a complete set of
SeqRecord objects (either as a list, or an iterator), an output file handle
(or in recent versions of Biopython an output filename as a string) and of
course the file format::

  from Bio import SeqIO
  records = ...
  SeqIO.write(records, "example.faa", "fasta")

Or, using a handle::

    from Bio import SeqIO
    records = ...
    with open("example.faa", "w") as handle:
      SeqIO.write(records, handle, "fasta")

You are expected to call this function once (with all your records) and if
using a handle, make sure you close it to flush the data to the hard disk.


Output - Advanced
-----------------
The effect of calling write() multiple times on a single file will vary
depending on the file format, and is best avoided unless you have a strong
reason to do so.

If you give a filename, then each time you call write() the existing file
will be overwriten. For sequential files formats (e.g. fasta, genbank) each
"record block" holds a single sequence.  For these files it would probably
be safe to call write() multiple times by re-using the same handle.


However, trying this for certain alignment formats (e.g. phylip, clustal,
stockholm) would have the effect of concatenating several multiple sequence
alignments together.  Such files are created by the PHYLIP suite of programs
for bootstrap analysis, but it is clearer to do this via Bio.AlignIO instead.


Conversion
----------
The Bio.SeqIO.convert(...) function allows an easy interface for simple
file format conversions. Additionally, it may use file format specific
optimisations so this should be the fastest way too.

In general however, you can combine the Bio.SeqIO.parse(...) function with
the Bio.SeqIO.write(...) function for sequence file conversion. Using
generator expressions or generator functions provides a memory efficient way
to perform filtering or other extra operations as part of the process.


File Formats
------------
When specifying the file format, use lowercase strings.  The same format
names are also used in Bio.AlignIO and include the following:

    - abif    - Applied Biosystem's sequencing trace format
    - ace     - Reads the contig sequences from an ACE assembly file.
    - embl    - The EMBL flat file format. Uses Bio.GenBank internally.
    - fasta   - The generic sequence file format where each record starts with
      an identifer line starting with a ">" character, followed by
      lines of sequence.
    - fastq   - A "FASTA like" format used by Sanger which also stores PHRED
      sequence quality values (with an ASCII offset of 33).
    - fastq-sanger - An alias for "fastq" for consistency with BioPerl and EMBOSS
    - fastq-solexa - Original Solexa/Illumnia variant of the FASTQ format which
      encodes Solexa quality scores (not PHRED quality scores) with an
      ASCII offset of 64.
    - fastq-illumina - Solexa/Illumina 1.3 to 1.7 variant of the FASTQ format
      which encodes PHRED quality scores with an ASCII offset of 64
      (not 33). Note as of version 1.8 of the CASAVA pipeline Illumina
      will produce FASTQ files using the standard Sanger encoding.
    - genbank - The GenBank or GenPept flat file format.
    - gb      - An alias for "genbank", for consistency with NCBI Entrez Utilities
    - ig      - The IntelliGenetics file format, apparently the same as the
      MASE alignment format.
    - imgt    - An EMBL like format from IMGT where the feature tables are more
      indented to allow for longer feature types.
    - phd     - Output from PHRED, used by PHRAP and CONSED for input.
    - pir     - A "FASTA like" format introduced by the National Biomedical
      Research Foundation (NBRF) for the Protein Information Resource
      (PIR) database, now part of UniProt.
    - seqxml  - SeqXML, simple XML format described in Schmitt et al (2011).
    - sff     - Standard Flowgram Format (SFF), typical output from Roche 454.
    - sff-trim - Standard Flowgram Format (SFF) with given trimming applied.
    - swiss   - Plain text Swiss-Prot aka UniProt format.
    - tab     - Simple two column tab separated sequence files, where each
      line holds a record's identifier and sequence. For example,
      this is used as by Aligent's eArray software when saving
      microarray probes in a minimal tab delimited text file.
    - qual    - A "FASTA like" format holding PHRED quality values from
      sequencing DNA, but no actual sequences (usually provided
      in separate FASTA files).
    - uniprot-xml - The UniProt XML format (replacement for the SwissProt plain
      text format which we call "swiss")

Note that while Bio.SeqIO can read all the above file formats, it cannot
write to all of them.

You can also use any file format supported by Bio.AlignIO, such as "nexus",
"phlip" and "stockholm", which gives you access to the individual sequences
making up each alignment as SeqRecords.
"""

from __future__ import print_function
from Bio._py3k import basestring

__docformat__ = "restructuredtext en"  # not just plaintext

# TODO
# - define policy on reading aligned sequences with gaps in
#   (e.g. - and . characters) including how the alphabet interacts
#
# - How best to handle unique/non unique record.id when writing.
#   For most file formats reading such files is fine; The stockholm
#   parser would fail.
#
# - MSF multiple alignment format, aka GCG, aka PileUp format (*.msf)
#   http://www.bioperl.org/wiki/MSF_multiple_alignment_format

"""
FAO BioPython Developers
------------------------
The way I envision this SeqIO system working as that for any sequence file
format we have an iterator that returns SeqRecord objects.

This also applies to interlaced fileformats (like clustal - although that
is now handled via Bio.AlignIO instead) where the file cannot be read record
by record.  You should still return an iterator, even if the implementation
could just as easily return a list.

These file format specific sequence iterators may be implemented as:
    - Classes which take a handle for __init__ and provide the __iter__ method
    - Functions that take a handle, and return an iterator object
    - Generator functions that take a handle, and yield SeqRecord objects

It is then trivial to turn this iterator into a list of SeqRecord objects,
an in memory dictionary, or a multiple sequence alignment object.

For building the dictionary by default the id propery of each SeqRecord is
used as the key.  You should always populate the id property, and it should
be unique in most cases. For some file formats the accession number is a good
choice.  If the file itself contains ambiguous identifiers, don't try and
dis-ambiguate them - return them as is.

When adding a new file format, please use the same lower case format name
as BioPerl, or if they have not defined one, try the names used by EMBOSS.

See also http://biopython.org/wiki/SeqIO_dev

--Peter
"""


from Bio.File import as_handle
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import Alphabet, AlphabetEncoder, _get_base_alphabet

from . import AbiIO
from . import AceIO
from . import FastaIO
from . import IgIO  # IntelliGenetics or MASE format
from . import InsdcIO  # EMBL and GenBank
from . import PdbIO
from . import PhdIO
from . import PirIO
from . import SeqXmlIO
from . import SffIO
from . import SwissIO
from . import TabIO
from . import QualityIO  # FastQ and qual files
from . import UniprotIO


# Convention for format names is "mainname-subtype" in lower case.
# Please use the same names as BioPerl or EMBOSS where possible.
#
# Note that this simple system copes with defining
# multiple possible iterators for a given format/extension
# with the -subtype suffix
#
# Most alignment file formats will be handled via Bio.AlignIO

_FormatToIterator = {"fasta": FastaIO.FastaIterator,
                     "gb": InsdcIO.GenBankIterator,
                     "genbank": InsdcIO.GenBankIterator,
                     "genbank-cds": InsdcIO.GenBankCdsFeatureIterator,
                     "embl": InsdcIO.EmblIterator,
                     "embl-cds": InsdcIO.EmblCdsFeatureIterator,
                     "imgt": InsdcIO.ImgtIterator,
                     "ig": IgIO.IgIterator,
                     "swiss": SwissIO.SwissIterator,
                     "pdb-atom": PdbIO.PdbAtomIterator,
                     "pdb-seqres": PdbIO.PdbSeqresIterator,
                     "phd": PhdIO.PhdIterator,
                     "ace": AceIO.AceIterator,
                     "tab": TabIO.TabIterator,
                     "pir": PirIO.PirIterator,
                     "fastq": QualityIO.FastqPhredIterator,
                     "fastq-sanger": QualityIO.FastqPhredIterator,
                     "fastq-solexa": QualityIO.FastqSolexaIterator,
                     "fastq-illumina": QualityIO.FastqIlluminaIterator,
                     "qual": QualityIO.QualPhredIterator,
                     "sff": SffIO.SffIterator,
                     # Not sure about this in the long run:
                     "sff-trim": SffIO._SffTrimIterator,
                     "uniprot-xml": UniprotIO.UniprotIterator,
                     "seqxml": SeqXmlIO.SeqXmlIterator,
                     "abi": AbiIO.AbiIterator,
                     "abi-trim": AbiIO._AbiTrimIterator,
                     }

_FormatToWriter = {"fasta": FastaIO.FastaWriter,
                   "gb": InsdcIO.GenBankWriter,
                   "genbank": InsdcIO.GenBankWriter,
                   "embl": InsdcIO.EmblWriter,
                   "imgt": InsdcIO.ImgtWriter,
                   "tab": TabIO.TabWriter,
                   "fastq": QualityIO.FastqPhredWriter,
                   "fastq-sanger": QualityIO.FastqPhredWriter,
                   "fastq-solexa": QualityIO.FastqSolexaWriter,
                   "fastq-illumina": QualityIO.FastqIlluminaWriter,
                   "phd": PhdIO.PhdWriter,
                   "qual": QualityIO.QualPhredWriter,
                   "sff": SffIO.SffWriter,
                   "seqxml": SeqXmlIO.SeqXmlWriter,
                   }

_BinaryFormats = ["sff", "sff-trim", "abi", "abi-trim"]


def write(sequences, handle, format):
    """Write complete set of sequences to a file.

        - sequences - A list (or iterator) of SeqRecord objects, or (if using
          Biopython 1.54 or later) a single SeqRecord.
        - handle    - File handle object to write to, or filename as string
          (note older versions of Biopython only took a handle).
        - format    - lower case string describing the file format to write.

    You should close the handle after calling this function.

    Returns the number of records written (as an integer).
    """
    from Bio import AlignIO

    # Try and give helpful error messages:
    if not isinstance(format, basestring):
        raise TypeError("Need a string for the file format (lower case)")
    if not format:
        raise ValueError("Format required (lower case string)")
    if format != format.lower():
        raise ValueError("Format string '%s' should be lower case" % format)

    if isinstance(sequences, SeqRecord):
        # This raised an exception in order version of Biopython
        sequences = [sequences]

    if format in _BinaryFormats:
        mode = 'wb'
    else:
        mode = 'w'

    with as_handle(handle, mode) as fp:
        # Map the file format to a writer class
        if format in _FormatToWriter:
            writer_class = _FormatToWriter[format]
            count = writer_class(fp).write_file(sequences)
        elif format in AlignIO._FormatToWriter:
            # Try and turn all the records into a single alignment,
            # and write that using Bio.AlignIO
            alignment = MultipleSeqAlignment(sequences)
            alignment_count = AlignIO.write([alignment], fp, format)
            assert alignment_count == 1, \
                "Internal error - the underlying writer " \
                " should have returned 1, not %s" % repr(alignment_count)
            count = len(alignment)
            del alignment_count, alignment
        elif format in _FormatToIterator or format in AlignIO._FormatToIterator:
            raise ValueError("Reading format '%s' is supported, but not writing"
                             % format)
        else:
            raise ValueError("Unknown format '%s'" % format)

        assert isinstance(count, int), "Internal error - the underlying %s " \
            "writer should have returned the record count, not %s" \
            % (format, repr(count))

    return count


def parse(handle, format, alphabet=None):
    r"""Turns a sequence file into an iterator returning SeqRecords.

        - handle   - handle to the file, or the filename as a string
          (note older versions of Biopython only took a handle).
        - format   - lower case string describing the file format.
        - alphabet - optional Alphabet object, useful when the sequence type
          cannot be automatically inferred from the file itself
          (e.g. format="fasta" or "tab")

    Typical usage, opening a file to read in, and looping over the record(s):

    >>> from Bio import SeqIO
    >>> filename = "Fasta/sweetpea.nu"
    >>> for record in SeqIO.parse(filename, "fasta"):
    ...    print("ID %s" % record.id)
    ...    print("Sequence length %i" % len(record))
    ...    print("Sequence alphabet %s" % record.seq.alphabet)
    ID gi|3176602|gb|U78617.1|LOU78617
    Sequence length 309
    Sequence alphabet SingleLetterAlphabet()

    For file formats like FASTA where the alphabet cannot be determined, it
    may be useful to specify the alphabet explicitly:

    >>> from Bio import SeqIO
    >>> from Bio.Alphabet import generic_dna
    >>> filename = "Fasta/sweetpea.nu"
    >>> for record in SeqIO.parse(filename, "fasta", generic_dna):
    ...    print("ID %s" % record.id)
    ...    print("Sequence length %i" % len(record))
    ...    print("Sequence alphabet %s" % record.seq.alphabet)
    ID gi|3176602|gb|U78617.1|LOU78617
    Sequence length 309
    Sequence alphabet DNAAlphabet()

    If you have a string 'data' containing the file contents, you must
    first turn this into a handle in order to parse it:

    >>> data = ">Alpha\nACCGGATGTA\n>Beta\nAGGCTCGGTTA\n"
    >>> from Bio import SeqIO
    >>> try:
    ...     from StringIO import StringIO # Python 2
    ... except ImportError:
    ...     from io import StringIO # Python 3
    ...
    >>> for record in SeqIO.parse(StringIO(data), "fasta"):
    ...     print("%s %s" % (record.id, record.seq))
    Alpha ACCGGATGTA
    Beta AGGCTCGGTTA

    Use the Bio.SeqIO.read(...) function when you expect a single record
    only.
    """
    # NOTE - The above docstring has some raw \n characters needed
    # for the StringIO example, hense the whole docstring is in raw
    # string mode (see the leading r before the opening quote).
    from Bio import AlignIO

    # Hack for SFF, will need to make this more general in future
    if format in _BinaryFormats:
        mode = 'rb'
    else:
        mode = 'rU'

    # Try and give helpful error messages:
    if not isinstance(format, basestring):
        raise TypeError("Need a string for the file format (lower case)")
    if not format:
        raise ValueError("Format required (lower case string)")
    if format != format.lower():
        raise ValueError("Format string '%s' should be lower case" % format)
    if alphabet is not None and not (isinstance(alphabet, Alphabet) or
                                     isinstance(alphabet, AlphabetEncoder)):
        raise ValueError("Invalid alphabet, %s" % repr(alphabet))

    with as_handle(handle, mode) as fp:
        # Map the file format to a sequence iterator:
        if format in _FormatToIterator:
            iterator_generator = _FormatToIterator[format]
            if alphabet is None:
                i = iterator_generator(fp)
            else:
                try:
                    i = iterator_generator(fp, alphabet=alphabet)
                except TypeError:
                    i = _force_alphabet(iterator_generator(fp), alphabet)
        elif format in AlignIO._FormatToIterator:
            # Use Bio.AlignIO to read in the alignments
            i = (r for alignment in AlignIO.parse(fp, format,
                                                  alphabet=alphabet)
                 for r in alignment)
        else:
            raise ValueError("Unknown format '%s'" % format)
        # This imposes some overhead... wait until we drop Python 2.4 to fix it
        for r in i:
            yield r


def _force_alphabet(record_iterator, alphabet):
    """Iterate over records, over-riding the alphabet (PRIVATE)."""
    # Assume the alphabet argument has been pre-validated
    given_base_class = _get_base_alphabet(alphabet).__class__
    for record in record_iterator:
        if isinstance(_get_base_alphabet(record.seq.alphabet),
                      given_base_class):
            record.seq.alphabet = alphabet
            yield record
        else:
            raise ValueError("Specified alphabet %s clashes with "
                             "that determined from the file, %s"
                             % (repr(alphabet), repr(record.seq.alphabet)))


def read(handle, format, alphabet=None):
    """Turns a sequence file into a single SeqRecord.

        - handle   - handle to the file, or the filename as a string
          (note older versions of Biopython only took a handle).
        - format   - string describing the file format.
        - alphabet - optional Alphabet object, useful when the sequence type
          cannot be automatically inferred from the file itself
          (e.g. format="fasta" or "tab")

    This function is for use parsing sequence files containing
    exactly one record.  For example, reading a GenBank file:

    >>> from Bio import SeqIO
    >>> record = SeqIO.read("GenBank/arab1.gb", "genbank")
    >>> print("ID %s" % record.id)
    ID AC007323.5
    >>> print("Sequence length %i" % len(record))
    Sequence length 86436
    >>> print("Sequence alphabet %s" % record.seq.alphabet)
    Sequence alphabet IUPACAmbiguousDNA()

    If the handle contains no records, or more than one record,
    an exception is raised.  For example:

    >>> from Bio import SeqIO
    >>> record = SeqIO.read("GenBank/cor6_6.gb", "genbank")
    Traceback (most recent call last):
        ...
    ValueError: More than one record found in handle

    If however you want the first record from a file containing
    multiple records this function would raise an exception (as
    shown in the example above).  Instead use:

    >>> from Bio import SeqIO
    >>> record = next(SeqIO.parse("GenBank/cor6_6.gb", "genbank"))
    >>> print("First record's ID %s" % record.id)
    First record's ID X55053.1

    Use the Bio.SeqIO.parse(handle, format) function if you want
    to read multiple records from the handle.
    """
    iterator = parse(handle, format, alphabet)
    try:
        first = next(iterator)
    except StopIteration:
        first = None
    if first is None:
        raise ValueError("No records found in handle")
    try:
        second = next(iterator)
    except StopIteration:
        second = None
    if second is not None:
        raise ValueError("More than one record found in handle")
    return first


def to_dict(sequences, key_function=None):
    """Turns a sequence iterator or list into a dictionary.

        - sequences  - An iterator that returns SeqRecord objects,
          or simply a list of SeqRecord objects.
        - key_function - Optional callback function which when given a
          SeqRecord should return a unique key for the dictionary.

    e.g. key_function = lambda rec : rec.name
    or,  key_function = lambda rec : rec.description.split()[0]

    If key_function is omitted then record.id is used, on the assumption
    that the records objects returned are SeqRecords with a unique id.

    If there are duplicate keys, an error is raised.

    Example usage, defaulting to using the record.id as key:

    >>> from Bio import SeqIO
    >>> filename = "GenBank/cor6_6.gb"
    >>> format = "genbank"
    >>> id_dict = SeqIO.to_dict(SeqIO.parse(filename, format))
    >>> print(sorted(id_dict))
    ['AF297471.1', 'AJ237582.1', 'L31939.1', 'M81224.1', 'X55053.1', 'X62281.1']
    >>> print(id_dict["L31939.1"].description)
    Brassica rapa (clone bif72) kin mRNA, complete cds.

    A more complex example, using the key_function argument in order to
    use a sequence checksum as the dictionary key:

    >>> from Bio import SeqIO
    >>> from Bio.SeqUtils.CheckSum import seguid
    >>> filename = "GenBank/cor6_6.gb"
    >>> format = "genbank"
    >>> seguid_dict = SeqIO.to_dict(SeqIO.parse(filename, format),
    ...               key_function = lambda rec : seguid(rec.seq))
    >>> for key, record in sorted(seguid_dict.items()):
    ...     print("%s %s" % (key, record.id))
    /wQvmrl87QWcm9llO4/efg23Vgg AJ237582.1
    BUg6YxXSKWEcFFH0L08JzaLGhQs L31939.1
    SabZaA4V2eLE9/2Fm5FnyYy07J4 X55053.1
    TtWsXo45S3ZclIBy4X/WJc39+CY M81224.1
    l7gjJFE6W/S1jJn5+1ASrUKW/FA X62281.1
    uVEYeAQSV5EDQOnFoeMmVea+Oow AF297471.1

    This approach is not suitable for very large sets of sequences, as all
    the SeqRecord objects are held in memory. Instead, consider using the
    Bio.SeqIO.index() function (if it supports your particular file format).
    """
    if key_function is None:
        key_function = lambda rec: rec.id

    d = dict()
    for record in sequences:
        key = key_function(record)
        if key in d:
            raise ValueError("Duplicate key '%s'" % key)
        d[key] = record
    return d


def index(filename, format, alphabet=None, key_function=None):
    """Indexes a sequence file and returns a dictionary like object.

        - filename - string giving name of file to be indexed
        - format   - lower case string describing the file format
        - alphabet - optional Alphabet object, useful when the sequence type
          cannot be automatically inferred from the file itself
          (e.g. format="fasta" or "tab")
        - key_function - Optional callback function which when given a
          SeqRecord identifier string should return a unique
          key for the dictionary.

    This indexing function will return a dictionary like object, giving the
    SeqRecord objects as values:

    >>> from Bio import SeqIO
    >>> records = SeqIO.index("Quality/example.fastq", "fastq")
    >>> len(records)
    3
    >>> sorted(records)
    ['EAS54_6_R1_2_1_413_324', 'EAS54_6_R1_2_1_443_348', 'EAS54_6_R1_2_1_540_792']
    >>> print(records["EAS54_6_R1_2_1_540_792"].format("fasta"))
    >EAS54_6_R1_2_1_540_792
    TTGGCAGGCCAAGGCCGATGGATCA
    <BLANKLINE>
    >>> "EAS54_6_R1_2_1_540_792" in records
    True
    >>> print(records.get("Missing", None))
    None
    >>> records.close()

    If the file is BGZF compressed, this is detected automatically. Ordinary
    GZIP files are not supported:

    >>> from Bio import SeqIO
    >>> records = SeqIO.index("Quality/example.fastq.bgz", "fastq")
    >>> len(records)
    3
    >>> print(records["EAS54_6_R1_2_1_540_792"].seq)
    TTGGCAGGCCAAGGCCGATGGATCA
    >>> records.close()

    Note that this pseudo dictionary will not support all the methods of a
    true Python dictionary, for example values() is not defined since this
    would require loading all of the records into memory at once.

    When you call the index function, it will scan through the file, noting
    the location of each record. When you access a particular record via the
    dictionary methods, the code will jump to the appropriate part of the
    file and then parse that section into a SeqRecord.

    Note that not all the input formats supported by Bio.SeqIO can be used
    with this index function. It is designed to work only with sequential
    file formats (e.g. "fasta", "gb", "fastq") and is not suitable for any
    interlaced file format (e.g. alignment formats such as "clustal").

    For small files, it may be more efficient to use an in memory Python
    dictionary, e.g.

    >>> from Bio import SeqIO
    >>> records = SeqIO.to_dict(SeqIO.parse("Quality/example.fastq", "fastq"))
    >>> len(records)
    3
    >>> sorted(records)
    ['EAS54_6_R1_2_1_413_324', 'EAS54_6_R1_2_1_443_348', 'EAS54_6_R1_2_1_540_792']
    >>> print(records["EAS54_6_R1_2_1_540_792"].format("fasta"))
    >EAS54_6_R1_2_1_540_792
    TTGGCAGGCCAAGGCCGATGGATCA
    <BLANKLINE>

    As with the to_dict() function, by default the id string of each record
    is used as the key. You can specify a callback function to transform
    this (the record identifier string) into your preferred key. For example:

    >>> from Bio import SeqIO
    >>> def make_tuple(identifier):
    ...     parts = identifier.split("_")
    ...     return int(parts[-2]), int(parts[-1])
    >>> records = SeqIO.index("Quality/example.fastq", "fastq",
    ...                       key_function=make_tuple)
    >>> len(records)
    3
    >>> sorted(records)
    [(413, 324), (443, 348), (540, 792)]
    >>> print(records[(540, 792)].format("fasta"))
    >EAS54_6_R1_2_1_540_792
    TTGGCAGGCCAAGGCCGATGGATCA
    <BLANKLINE>
    >>> (540, 792) in records
    True
    >>> "EAS54_6_R1_2_1_540_792" in records
    False
    >>> print(records.get("Missing", None))
    None
    >>> records.close()

    Another common use case would be indexing an NCBI style FASTA file,
    where you might want to extract the GI number from the FASTA identifer
    to use as the dictionary key.

    Notice that unlike the to_dict() function, here the key_function does
    not get given the full SeqRecord to use to generate the key. Doing so
    would impose a severe performance penalty as it would require the file
    to be completely parsed while building the index. Right now this is
    usually avoided.

    See also: Bio.SeqIO.index_db() and Bio.SeqIO.to_dict()
    """
    # Try and give helpful error messages:
    if not isinstance(filename, basestring):
        raise TypeError("Need a filename (not a handle)")
    if not isinstance(format, basestring):
        raise TypeError("Need a string for the file format (lower case)")
    if not format:
        raise ValueError("Format required (lower case string)")
    if format != format.lower():
        raise ValueError("Format string '%s' should be lower case" % format)
    if alphabet is not None and not (isinstance(alphabet, Alphabet) or
                                     isinstance(alphabet, AlphabetEncoder)):
        raise ValueError("Invalid alphabet, %s" % repr(alphabet))

    # Map the file format to a sequence iterator:
    from ._index import _FormatToRandomAccess  # Lazy import
    from Bio.File import _IndexedSeqFileDict
    try:
        proxy_class = _FormatToRandomAccess[format]
    except KeyError:
        raise ValueError("Unsupported format %r" % format)
    repr = "SeqIO.index(%r, %r, alphabet=%r, key_function=%r)" \
        % (filename, format, alphabet, key_function)
    return _IndexedSeqFileDict(proxy_class(filename, format, alphabet),
                               key_function, repr, "SeqRecord")


def index_db(index_filename, filenames=None, format=None, alphabet=None,
             key_function=None):
    """Index several sequence files and return a dictionary like object.

    The index is stored in an SQLite database rather than in memory (as in the
    Bio.SeqIO.index(...) function).

        - index_filename - Where to store the SQLite index
        - filenames - list of strings specifying file(s) to be indexed, or when
          indexing a single file this can be given as a string.
          (optional if reloading an existing index, but must match)
        - format   - lower case string describing the file format
          (optional if reloading an existing index, but must match)
        - alphabet - optional Alphabet object, useful when the sequence type
          cannot be automatically inferred from the file itself
          (e.g. format="fasta" or "tab")
        - key_function - Optional callback function which when given a
          SeqRecord identifier string should return a unique
          key for the dictionary.

    This indexing function will return a dictionary like object, giving the
    SeqRecord objects as values:

    >>> from Bio.Alphabet import generic_protein
    >>> from Bio import SeqIO
    >>> files = ["GenBank/NC_000932.faa", "GenBank/NC_005816.faa"]
    >>> def get_gi(name):
    ...     parts = name.split("|")
    ...     i = parts.index("gi")
    ...     assert i != -1
    ...     return parts[i+1]
    >>> idx_name = ":memory:" #use an in memory SQLite DB for this test
    >>> records = SeqIO.index_db(idx_name, files, "fasta", generic_protein, get_gi)
    >>> len(records)
    95
    >>> records["7525076"].description
    'gi|7525076|ref|NP_051101.1| Ycf2 [Arabidopsis thaliana]'
    >>> records["45478717"].description
    'gi|45478717|ref|NP_995572.1| pesticin [Yersinia pestis biovar Microtus str. 91001]'
    >>> records.close()

    In this example the two files contain 85 and 10 records respectively.

    BGZF compressed files are supported, and detected automatically. Ordinary
    GZIP compressed files are not supported.

    See also: Bio.SeqIO.index() and Bio.SeqIO.to_dict(), and the Python module
    glob which is useful for building lists of files.
    """
    # Try and give helpful error messages:
    if not isinstance(index_filename, basestring):
        raise TypeError("Need a string for the index filename")
    if isinstance(filenames, basestring):
        # Make the API a little more friendly, and more similar
        # to Bio.SeqIO.index(...) for indexing just one file.
        filenames = [filenames]
    if filenames is not None and not isinstance(filenames, list):
        raise TypeError(
            "Need a list of filenames (as strings), or one filename")
    if format is not None and not isinstance(format, basestring):
        raise TypeError("Need a string for the file format (lower case)")
    if format and format != format.lower():
        raise ValueError("Format string '%s' should be lower case" % format)
    if alphabet is not None and not (isinstance(alphabet, Alphabet) or
                                     isinstance(alphabet, AlphabetEncoder)):
        raise ValueError("Invalid alphabet, %s" % repr(alphabet))

    # Map the file format to a sequence iterator:
    from ._index import _FormatToRandomAccess  # Lazy import
    from Bio.File import _SQLiteManySeqFilesDict
    repr = "SeqIO.index_db(%r, filenames=%r, format=%r, alphabet=%r, key_function=%r)" \
               % (index_filename, filenames, format, alphabet, key_function)

    def proxy_factory(format, filename=None):
        """Given a filename returns proxy object, else boolean if format OK."""
        if filename:
            return _FormatToRandomAccess[format](filename, format, alphabet)
        else:
            return format in _FormatToRandomAccess

    return _SQLiteManySeqFilesDict(index_filename, filenames,
                                   proxy_factory, format,
                                   key_function, repr)


def convert(in_file, in_format, out_file, out_format, alphabet=None):
    """Convert between two sequence file formats, return number of records.

        - in_file - an input handle or filename
        - in_format - input file format, lower case string
        - out_file - an output handle or filename
        - out_format - output file format, lower case string
        - alphabet - optional alphabet to assume

    **NOTE** - If you provide an output filename, it will be opened which will
    overwrite any existing file without warning. This may happen if even
    the conversion is aborted (e.g. an invalid out_format name is given).

    For example, going from a filename to a handle:

    >>> from Bio import SeqIO
    >>> try:
    ...     from StringIO import StringIO # Python 2
    ... except ImportError:
    ...     from io import StringIO # Python 3
    ...
    >>> handle = StringIO("")
    >>> SeqIO.convert("Quality/example.fastq", "fastq", handle, "fasta")
    3
    >>> print(handle.getvalue())
    >EAS54_6_R1_2_1_413_324
    CCCTTCTTGTCTTCAGCGTTTCTCC
    >EAS54_6_R1_2_1_540_792
    TTGGCAGGCCAAGGCCGATGGATCA
    >EAS54_6_R1_2_1_443_348
    GTTGCTTCTGGCGTGGGTGGGGGGG
    <BLANKLINE>
    """
    # Hack for SFF, will need to make this more general in future
    if in_format in _BinaryFormats:
        in_mode = 'rb'
    else:
        in_mode = 'rU'

    # Don't open the output file until we've checked the input is OK?
    if out_format in ["sff", "sff_trim"]:
        out_mode = 'wb'
    else:
        out_mode = 'w'

    # This will check the arguments and issue error messages,
    # after we have opened the file which is a shame.
    from ._convert import _handle_convert  # Lazy import
    with as_handle(in_file, in_mode) as in_handle:
        with as_handle(out_file, out_mode) as out_handle:
            count = _handle_convert(in_handle, in_format,
                                    out_handle, out_format,
                                    alphabet)
    return count


# This helpful trick for testing no longer works with the
# local imports :(
#
# if __name__ == "__main__":
#    from Bio._utils import run_doctest
#    run_doctest()
