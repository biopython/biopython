# Copyright 2006-2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
#Nice link:
# http://www.ebi.ac.uk/help/formats_frame.html

"""Sequence input/output as SeqRecord objects.

Bio.SeqIO is also documented at U{http://biopython.org/wiki/SeqIO} and by
a whole chapter in our tutorial:
 - U{http://biopython.org/DIST/docs/tutorial/Tutorial.html}
 - U{http://biopython.org/DIST/docs/tutorial/Tutorial.pdf}  

Input
=====
The main function is Bio.SeqIO.parse(...) which takes an input file handle,
and format string.  This returns an iterator giving SeqRecord objects:

    >>> from Bio import SeqIO
    >>> handle = open("Fasta/f002", "rU")
    >>> for record in SeqIO.parse(handle, "fasta") :
    ...     print record.id, len(record)
    gi|1348912|gb|G26680|G26680 633
    gi|1348917|gb|G26685|G26685 413
    gi|1592936|gb|G29385|G29385 471
    >>> handle.close()

Note that the parse() function will all invoke the relevant parser for the
format with its default settings.  You may want more control, in which case
you need to create a format specific sequence iterator directly.

Input - Single Records
======================
If you expect your file to contain one-and-only-one record, then we provide
the following 'helper' function which will return a single SeqRecord, or
raise an exception if there are no records or more than one record:

    >>> from Bio import SeqIO
    >>> handle = open("Fasta/f001", "rU")
    >>> record = SeqIO.read(handle, "fasta")
    >>> handle.close()
    >>> print record.id, len(record)
    gi|3318709|pdb|1A91| 79

This style is useful when you expect a single record only (and would
consider multiple records an error).  For example, when dealing with GenBank
files for bacterial genomes or chromosomes, there is normally only a single
record.  Alternatively, use this with a handle when download a single record
from the internet.

However, if you just want the first record from a file containing multiple
record, use the iterator's next() method:

    >>> from Bio import SeqIO
    >>> handle = open("Fasta/f002", "rU")
    >>> record = SeqIO.parse(handle, "fasta").next()
    >>> handle.close()
    >>> print record.id, len(record)
    gi|1348912|gb|G26680|G26680 633

The above code will work as long as the file contains at least one record.
Note that if there is more than one record, the remaining records will be
silently ignored.


Input - Multiple Records
========================
For non-interlaced files (e.g. Fasta, GenBank, EMBL) with multiple records
using a sequence iterator can save you a lot of memory (RAM).  There is
less benefit for interlaced file formats (e.g. most multiple alignment file
formats).  However, an iterator only lets you access the records one by one.

If you want random access to the records by number, turn this into a list:

    >>> from Bio import SeqIO
    >>> handle = open("Fasta/f002", "rU")
    >>> records = list(SeqIO.parse(handle, "fasta"))
    >>> handle.close()
    >>> print records[1].id
    gi|1348917|gb|G26685|G26685

If you want random access to the records by a key such as the record id,
turn the iterator into a dictionary:

    >>> from Bio import SeqIO
    >>> handle = open("Fasta/f002", "rU")
    >>> record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    >>> handle.close()
    >>> len(record_dict)
    3
    >>> print len(record_dict["gi|1348917|gb|G26685|G26685"])
    413

However, using list() or the to_dict() function will load all the records
into memory at once, and therefore is not possible on very large files.
Instead, for *some* file formats Bio.SeqIO provides an indexing approach
providing dictionary like access to any record. For example,

    >>> from Bio import SeqIO
    >>> record_dict = SeqIO.index("Fasta/f002", "fasta")
    >>> len(record_dict)
    3
    >>> print len(record_dict["gi|1348917|gb|G26685|G26685"])
    413

Many but not all of the supported input file formats can be indexed like
this. For example "fasta", "fastq", "qual" and even the binary format "sff"
work, but alignment formats like "phylip", "clustalw" and "nexus" will not.

Input - Alignments
==================
You can read in alignment files as Alignment objects using Bio.AlignIO.
Alternatively, reading in an alignment file format via Bio.SeqIO will give
you a SeqRecord for each row of each alignment:

    >>> from Bio import SeqIO
    >>> handle = open("Clustalw/hedgehog.aln", "rU")
    >>> for record in SeqIO.parse(handle, "clustal") :
    ...     print record.id, len(record)
    gi|167877390|gb|EDS40773.1| 447
    gi|167234445|ref|NP_001107837. 447
    gi|74100009|gb|AAZ99217.1| 447
    gi|13990994|dbj|BAA33523.2| 447
    gi|56122354|gb|AAV74328.1| 447
    >>> handle.close()

Output
======
Use the function Bio.SeqIO.write(...), which takes a complete set of
SeqRecord objects (either as a list, or an iterator), an output file handle
and of course the file format::

    from Bio import SeqIO
    records = ...
    handle = open("example.faa", "w")
    SeqIO.write(records, handle, "fasta")
    handle.close()

In general, you are expected to call this function once (with all your
records) and then close the file handle.

Output - Advanced
=================
The effect of calling write() multiple times on a single file will vary
depending on the file format, and is best avoided unless you have a strong
reason to do so.

Trying this for certain alignment formats (e.g. phylip, clustal, stockholm)
would have the effect of concatenating several multiple sequence alignments
together.  Such files are created by the PHYLIP suite of programs for
bootstrap analysis.

For sequential files formats (e.g. fasta, genbank) each "record block" holds
a single sequence.  For these files it would probably be safe to call
write() multiple times.

Conversion
==========
The Bio.SeqIO.convert(...) function allows an easy interface for simple
file format conversions. Additionally, it may use file format specific
optimisations so this should be the fastest way too.

In general however, you can combine the Bio.SeqIO.parse(...) function with the
Bio.SeqIO.write(...) function for sequence file conversion. Using generator
expressions provides a memory efficient way to perform filtering or other
extra operations as part of the process.

File Formats
============
When specifying the file format, use lowercase strings.  The same format
names are also used in Bio.AlignIO and include the following:

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
 - fastq-illumina - Solexa/Illumnia 1.3+ variant of the FASTQ format which
             encodes PHRED quality scores with an ASCII offset of 64 (not 33).
 - genbank - The GenBank or GenPept flat file format.
 - gb      - An alias for "genbank", for consistency with NCBI Entrez Utilities
 - ig      - The IntelliGenetics file format, apparently the same as the
             MASE alignment format.
 - phd     - Output from PHRED, used by PHRAP and CONSED for input.
 - pir     - A "FASTA like" format introduced by the National Biomedical
             Research Foundation (NBRF) for the Protein Information Resource
             (PIR) database, now part of UniProt.
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

Note that while Bio.SeqIO can read all the above file formats, it cannot
write to all of them.

You can also use any file format supported by Bio.AlignIO, such as "nexus",
"phlip" and "stockholm", which gives you access to the individual sequences
making up each alignment as SeqRecords.
"""
__docformat__ = "epytext en" #not just plaintext

#TODO
# - define policy on reading aligned sequences with gaps in
#   (e.g. - and . characters) including how the alphabet interacts
#
# - Can we build the to_alignment(...) functionality
#   into the generic Alignment class instead?
#
# - How best to handle unique/non unique record.id when writing.
#   For most file formats reading such files is fine; The stockholm
#   parser would fail.
#
# - MSF multiple alignment format, aka GCG, aka PileUp format (*.msf)
#   http://www.bioperl.org/wiki/MSF_multiple_alignment_format 

"""
FAO BioPython Developers
========================
The way I envision this SeqIO system working as that for any sequence file
format we have an iterator that returns SeqRecord objects.

This also applies to interlaced fileformats (like clustal - although that
is now handled via Bio.AlignIO instead) where the file cannot be read record
by record.  You should still return an iterator, even if the implementation
could just as easily return a list.

These file format specific sequence iterators may be implemented as:
* Classes which take a handle for __init__ and provide the __iter__ method
* Functions that take a handle, and return an iterator object
* Generator functions that take a handle, and yield SeqRecord objects

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

import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Generic import Alignment
from Bio.Alphabet import Alphabet, AlphabetEncoder, _get_base_alphabet

import AceIO
import FastaIO
import IgIO #IntelliGenetics or MASE format
import InsdcIO #EMBL and GenBank
import PhdIO
import PirIO
import SffIO
import SwissIO
import TabIO
import QualityIO #FastQ and qual files


#Convention for format names is "mainname-subtype" in lower case.
#Please use the same names as BioPerl or EMBOSS where possible.
#
#Note that this simple system copes with defining
#multiple possible iterators for a given format/extension
#with the -subtype suffix
#
#Most alignment file formats will be handled via Bio.AlignIO

_FormatToIterator ={"fasta" : FastaIO.FastaIterator,
                    "gb" : InsdcIO.GenBankIterator,
                    "genbank" : InsdcIO.GenBankIterator,
                    "genbank-cds" : InsdcIO.GenBankCdsFeatureIterator,
                    "embl" : InsdcIO.EmblIterator,
                    "embl-cds" : InsdcIO.EmblCdsFeatureIterator,
                    "ig" : IgIO.IgIterator,
                    "swiss" : SwissIO.SwissIterator,
                    "phd" : PhdIO.PhdIterator,
                    "ace" : AceIO.AceIterator,
                    "tab" : TabIO.TabIterator,
                    "pir" : PirIO.PirIterator,
                    "fastq" : QualityIO.FastqPhredIterator,
                    "fastq-sanger" : QualityIO.FastqPhredIterator,
                    "fastq-solexa" : QualityIO.FastqSolexaIterator,
                    "fastq-illumina" : QualityIO.FastqIlluminaIterator,
                    "qual" : QualityIO.QualPhredIterator,
                    "sff": SffIO.SffIterator,
                    #Not sure about this in the long run:
                    "sff-trim": SffIO._SffTrimIterator,
                    }

_FormatToWriter ={"fasta" : FastaIO.FastaWriter,
                  "gb" : InsdcIO.GenBankWriter,
                  "genbank" : InsdcIO.GenBankWriter,
                  "tab" : TabIO.TabWriter,
                  "fastq" : QualityIO.FastqPhredWriter,
                  "fastq-sanger" : QualityIO.FastqPhredWriter,
                  "fastq-solexa" : QualityIO.FastqSolexaWriter,
                  "fastq-illumina" : QualityIO.FastqIlluminaWriter,
                  "phd" : PhdIO.PhdWriter,
                  "qual" : QualityIO.QualPhredWriter,
                  "sff" : SffIO.SffWriter,
                  }

def write(sequences, handle, format) :
    """Write complete set of sequences to a file.

     - sequences - A list (or iterator) of SeqRecord objects.
     - handle    - File handle object to write to.
     - format    - lower case string describing the file format to write.

    You should close the handle after calling this function.

    Returns the number of records written (as an integer).
    """
    from Bio import AlignIO

    #Try and give helpful error messages:
    if isinstance(handle, basestring) :
        raise TypeError("Need a file handle, not a string (i.e. not a filename)")
    if not isinstance(format, basestring) :
        raise TypeError("Need a string for the file format (lower case)")
    if not format :
        raise ValueError("Format required (lower case string)")
    if format != format.lower() :
        raise ValueError("Format string '%s' should be lower case" % format)
    if isinstance(sequences,SeqRecord):
        raise ValueError("Use a SeqRecord list/iterator, not just a single SeqRecord")

    #Map the file format to a writer class
    if format in _FormatToWriter :
        writer_class = _FormatToWriter[format]
        count = writer_class(handle).write_file(sequences)
    elif format in AlignIO._FormatToWriter :
        #Try and turn all the records into a single alignment,
        #and write that using Bio.AlignIO
        alignment = to_alignment(sequences)
        alignment_count = AlignIO.write([alignment], handle, format)
        assert alignment_count == 1, "Internal error - the underlying writer " \
           + " should have returned 1, not %s" % repr(alignment_count)
        count = len(alignment.get_all_seqs())
        del alignment_count, alignment
    elif format in _FormatToIterator or format in AlignIO._FormatToIterator :
        raise ValueError("Reading format '%s' is supported, but not writing" \
                         % format)
    else :
        raise ValueError("Unknown format '%s'" % format)

    assert isinstance(count, int), "Internal error - the underlying %s " \
           "writer should have returned the record count, not %s" \
           % (format, repr(count))
    return count
    
def parse(handle, format, alphabet=None) :
    r"""Turns a sequence file into an iterator returning SeqRecords.

     - handle   - handle to the file.
     - format   - lower case string describing the file format.
     - alphabet - optional Alphabet object, useful when the sequence type
                  cannot be automatically inferred from the file itself
                  (e.g. format="fasta" or "tab")

    Typical usage, opening a file to read in, and looping over the record(s):

    >>> from Bio import SeqIO
    >>> filename = "Nucleic/sweetpea.nu"
    >>> for record in SeqIO.parse(open(filename,"rU"), "fasta") :
    ...    print "ID", record.id
    ...    print "Sequence length", len(record)
    ...    print "Sequence alphabet", record.seq.alphabet
    ID gi|3176602|gb|U78617.1|LOU78617
    Sequence length 309
    Sequence alphabet SingleLetterAlphabet()

    For file formats like FASTA where the alphabet cannot be determined, it
    may be useful to specify the alphabet explicitly:

    >>> from Bio import SeqIO
    >>> from Bio.Alphabet import generic_dna
    >>> filename = "Nucleic/sweetpea.nu"
    >>> for record in SeqIO.parse(open(filename,"rU"), "fasta", generic_dna) :
    ...    print "ID", record.id
    ...    print "Sequence length", len(record)
    ...    print "Sequence alphabet", record.seq.alphabet
    ID gi|3176602|gb|U78617.1|LOU78617
    Sequence length 309
    Sequence alphabet DNAAlphabet()

    If you have a string 'data' containing the file contents, you must
    first turn this into a handle in order to parse it:

    >>> data = ">Alpha\nACCGGATGTA\n>Beta\nAGGCTCGGTTA\n"
    >>> from Bio import SeqIO
    >>> from StringIO import StringIO
    >>> for record in SeqIO.parse(StringIO(data), "fasta") :
    ...     print record.id, record.seq
    Alpha ACCGGATGTA
    Beta AGGCTCGGTTA

    Use the Bio.SeqIO.read(handle, format) function when you expect a single
    record only.
    """
    #NOTE - The above docstring has some raw \n characters needed
    #for the StringIO example, hense the whole docstring is in raw
    #string mode (see the leading r before the opening quote).
    from Bio import AlignIO

    #Try and give helpful error messages:
    if isinstance(handle, basestring) :
        raise TypeError("Need a file handle, not a string (i.e. not a filename)")
    if not isinstance(format, basestring) :
        raise TypeError("Need a string for the file format (lower case)")
    if not format :
        raise ValueError("Format required (lower case string)")
    if format != format.lower() :
        raise ValueError("Format string '%s' should be lower case" % format)
    if alphabet is not None and not (isinstance(alphabet, Alphabet) or \
                                     isinstance(alphabet, AlphabetEncoder)) :
        raise ValueError("Invalid alphabet, %s" % repr(alphabet))

    #Map the file format to a sequence iterator:    
    if format in _FormatToIterator :
        iterator_generator = _FormatToIterator[format]
        if alphabet is None :
            return iterator_generator(handle)
        try :
            return iterator_generator(handle, alphabet=alphabet)
        except :
            return _force_alphabet(iterator_generator(handle), alphabet)
    elif format in AlignIO._FormatToIterator :
        #Use Bio.AlignIO to read in the alignments
        #TODO - Can this helper function can be replaced with a generator expression,
        #or something from itertools?
        return _iterate_via_AlignIO(handle, format, alphabet)
    else :
        raise ValueError("Unknown format '%s'" % format)

#This is a generator function
def _iterate_via_AlignIO(handle, format, alphabet) :
    """Iterate over all records in several alignments (PRIVATE)."""
    from Bio import AlignIO
    for align in AlignIO.parse(handle, format, alphabet=alphabet) :
        for record in align :
            yield record

def _force_alphabet(record_iterator, alphabet) :
     """Iterate over records, over-riding the alphabet (PRIVATE)."""
     #Assume the alphabet argument has been pre-validated
     given_base_class = _get_base_alphabet(alphabet).__class__
     for record in record_iterator :
         if isinstance(_get_base_alphabet(record.seq.alphabet),
                       given_base_class) :
             record.seq.alphabet = alphabet
             yield record
         else :
             raise ValueError("Specified alphabet %s clashes with "\
                              "that determined from the file, %s" \
                              % (repr(alphabet), repr(record.seq.alphabet)))

def read(handle, format, alphabet=None) :
    """Turns a sequence file into a single SeqRecord.

     - handle   - handle to the file.
     - format   - string describing the file format.
     - alphabet - optional Alphabet object, useful when the sequence type
                  cannot be automatically inferred from the file itself
                  (e.g. format="fasta" or "tab")

    This function is for use parsing sequence files containing
    exactly one record.  For example, reading a GenBank file:

    >>> from Bio import SeqIO
    >>> record = SeqIO.read(open("GenBank/arab1.gb", "rU"), "genbank")
    >>> print "ID", record.id
    ID AC007323.5
    >>> print "Sequence length", len(record)
    Sequence length 86436
    >>> print "Sequence alphabet", record.seq.alphabet
    Sequence alphabet IUPACAmbiguousDNA()

    If the handle contains no records, or more than one record,
    an exception is raised.  For example:

    >>> from Bio import SeqIO
    >>> record = SeqIO.read(open("GenBank/cor6_6.gb", "rU"), "genbank")
    Traceback (most recent call last):
        ...
    ValueError: More than one record found in handle

    If however you want the first record from a file containing
    multiple records this function would raise an exception (as
    shown in the example above).  Instead use:

    >>> from Bio import SeqIO
    >>> record = SeqIO.parse(open("GenBank/cor6_6.gb", "rU"), "genbank").next()
    >>> print "First record's ID", record.id
    First record's ID X55053.1

    Use the Bio.SeqIO.parse(handle, format) function if you want
    to read multiple records from the handle.
    """
    iterator = parse(handle, format, alphabet)
    try :
        first = iterator.next()
    except StopIteration :
        first = None
    if first is None :
        raise ValueError("No records found in handle")
    try :
        second = iterator.next()
    except StopIteration :
        second = None
    if second is not None :
        raise ValueError("More than one record found in handle")
    return first

def to_dict(sequences, key_function=None) :
    """Turns a sequence iterator or list into a dictionary.

     - sequences  - An iterator that returns SeqRecord objects,
                    or simply a list of SeqRecord objects.
     - key_function - Optional callback function which when given a
                    SeqRecord should return a unique key for the dictionary.

    e.g. key_function = lambda rec : rec.name
    or,  key_function = lambda rec : rec.description.split()[0]

    If key_function is ommitted then record.id is used, on the assumption
    that the records objects returned are SeqRecords with a unique id.
    
    If there are duplicate keys, an error is raised.

    Example usage, defaulting to using the record.id as key:

    >>> from Bio import SeqIO
    >>> handle = open("GenBank/cor6_6.gb", "rU")
    >>> format = "genbank"
    >>> id_dict = SeqIO.to_dict(SeqIO.parse(handle, format))
    >>> print sorted(id_dict.keys())
    ['AF297471.1', 'AJ237582.1', 'L31939.1', 'M81224.1', 'X55053.1', 'X62281.1']
    >>> print id_dict["L31939.1"].description
    Brassica rapa (clone bif72) kin mRNA, complete cds.

    A more complex example, using the key_function argument in order to
    use a sequence checksum as the dictionary key:

    >>> from Bio import SeqIO
    >>> from Bio.SeqUtils.CheckSum import seguid
    >>> handle = open("GenBank/cor6_6.gb", "rU")
    >>> format = "genbank"
    >>> seguid_dict = SeqIO.to_dict(SeqIO.parse(handle, format),
    ...               key_function = lambda rec : seguid(rec.seq))
    >>> for key, record in sorted(seguid_dict.iteritems()) :
    ...     print key, record.id
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
    if key_function is None :
        key_function = lambda rec : rec.id

    d = dict()
    for record in sequences :
        key = key_function(record)
        if key in d :
            raise ValueError("Duplicate key '%s'" % key)
        d[key] = record
    return d

def index(filename, format, alphabet=None, key_function=None) :
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
    >>> sorted(records.keys())
    ['EAS54_6_R1_2_1_413_324', 'EAS54_6_R1_2_1_443_348', 'EAS54_6_R1_2_1_540_792']
    >>> print records["EAS54_6_R1_2_1_540_792"].format("fasta")
    >EAS54_6_R1_2_1_540_792
    TTGGCAGGCCAAGGCCGATGGATCA
    <BLANKLINE>
    >>> "EAS54_6_R1_2_1_540_792" in records
    True
    >>> print records.get("Missing", None)
    None

    Note that this psuedo dictionary will not support all the methods of a
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
    >>> records = SeqIO.to_dict(SeqIO.parse(open("Quality/example.fastq"), "fastq"))
    >>> len(records)
    3
    >>> sorted(records.keys())
    ['EAS54_6_R1_2_1_413_324', 'EAS54_6_R1_2_1_443_348', 'EAS54_6_R1_2_1_540_792']
    >>> print records["EAS54_6_R1_2_1_540_792"].format("fasta")
    >EAS54_6_R1_2_1_540_792
    TTGGCAGGCCAAGGCCGATGGATCA
    <BLANKLINE>

    As with the to_dict() function, by default the id string of each record
    is used as the key. You can specify a callback function to transform
    this (the record identifier string) into your prefered key. For example:

    >>> from Bio import SeqIO
    >>> def make_tuple(identifier) :
    ...     parts = identifier.split("_")
    ...     return int(parts[-2]), int(parts[-1])
    >>> records = SeqIO.index("Quality/example.fastq", "fastq",
    ...                       key_function=make_tuple)
    >>> len(records)
    3
    >>> sorted(records.keys())
    [(413, 324), (443, 348), (540, 792)]
    >>> print records[(540, 792)].format("fasta")
    >EAS54_6_R1_2_1_540_792
    TTGGCAGGCCAAGGCCGATGGATCA
    <BLANKLINE>
    >>> (540, 792) in records
    True
    >>> "EAS54_6_R1_2_1_540_792" in records
    False
    >>> print records.get("Missing", None)
    None

    Another common use case would be indexing an NCBI style FASTA file,
    where you might want to extract the GI number from the FASTA identifer
    to use as the dictionary key.

    Notice that unlike the to_dict() function, here the key_function does
    not get given the full SeqRecord to use to generate the key. Doing so
    would impose a severe performance penalty as it would require the file
    to be completely parsed while building the index. Right now this is
    usually avoided.
    """
    #Try and give helpful error messages:
    if not isinstance(filename, basestring) :
        raise TypeError("Need a filename (not a handle)")
    if not isinstance(format, basestring) :
        raise TypeError("Need a string for the file format (lower case)")
    if not format :
        raise ValueError("Format required (lower case string)")
    if format != format.lower() :
        raise ValueError("Format string '%s' should be lower case" % format)
    if alphabet is not None and not (isinstance(alphabet, Alphabet) or \
                                     isinstance(alphabet, AlphabetEncoder)) :
        raise ValueError("Invalid alphabet, %s" % repr(alphabet))

    #Map the file format to a sequence iterator:    
    import _index #Lazy import
    try :
        indexer = _index._FormatToIndexedDict[format]
    except KeyError :
        raise ValueError("Unsupported format '%s'" % format)
    return indexer(filename, alphabet, key_function)

def to_alignment(sequences, alphabet=None, strict=True) :
    """Returns a multiple sequence alignment (OBSOLETE).

     - sequences -An iterator that returns SeqRecord objects,
                  or simply a list of SeqRecord objects.  All
                  the record sequences must be the same length.
     - alphabet - Optional alphabet.  Stongly recommended.
     - strict   - Optional, defaults to True.  Should error checking
                  be done?

    Using this function is now discouraged.  Rather doing this:

    >>> from Bio import SeqIO
    >>> handle = open("Clustalw/protein.aln")
    >>> alignment = SeqIO.to_alignment(SeqIO.parse(handle, "clustal"))
    >>> handle.close()

    You are now encouraged to use Bio.AlignIO instead, e.g.

    >>> from Bio import AlignIO
    >>> handle = open("Clustalw/protein.aln")
    >>> alignment = AlignIO.read(handle, "clustal")
    >>> handle.close()
    """
    #TODO - Move this functionality into the Alignment class instead?
    from Bio.Alphabet import generic_alphabet
    from Bio.Alphabet import _consensus_alphabet
    if alphabet is None :
        sequences = list(sequences)
        alphabet = _consensus_alphabet([rec.seq.alphabet for rec in sequences \
                                        if rec.seq is not None])

    if not (isinstance(alphabet, Alphabet) or isinstance(alphabet, AlphabetEncoder)) :
        raise ValueError("Invalid alphabet")

    alignment_length = None
    alignment = Alignment(alphabet)
    for record in sequences :
        if strict :
            if alignment_length is None :
                alignment_length = len(record.seq)
            elif alignment_length != len(record.seq) :
                raise ValueError("Sequences must all be the same length")

            assert isinstance(record.seq.alphabet, Alphabet) \
            or isinstance(record.seq.alphabet, AlphabetEncoder), \
                "Sequence does not have a valid alphabet"

            #TODO - Move this alphabet comparison code into the Alphabet module/class?
            #TODO - Is a normal alphabet "ungapped" by default, or does it just mean
            #undecided?
            if isinstance(record.seq.alphabet, Alphabet) \
            and isinstance(alphabet, Alphabet) :
                #Comparing two non-gapped alphabets            
                if not isinstance(record.seq.alphabet, alphabet.__class__) :
                    raise ValueError("Incompatible sequence alphabet " \
                                     + "%s for %s alignment" \
                                     % (record.seq.alphabet, alphabet))
            elif isinstance(record.seq.alphabet, AlphabetEncoder) \
            and isinstance(alphabet, Alphabet) :
                raise ValueError("Sequence has a gapped alphabet, alignment does not")
            elif isinstance(record.seq.alphabet, Alphabet) \
            and isinstance(alphabet, Gapped) :
                #Sequence isn't gapped, alignment is.
                if not isinstance(record.seq.alphabet, alphabet.alphabet.__class__) :
                    raise ValueError("Incompatible sequence alphabet " \
                                     + "%s for %s alignment" \
                                     % (record.seq.alphabet, alphabet))
            else :
                #Comparing two gapped alphabets
                if not isinstance(record.seq.alphabet, alphabet.__class__) :
                    raise ValueError("Incompatible sequence alphabet " \
                                     + "%s for %s alignment" \
                                     % (record.seq.alphabet, alphabet))
                if record.seq.alphabet.gap_char != alphabet.gap_char :
                    raise ValueError("Sequence gap characters != alignment gap char")
            #ToDo, additional checks on the specified alignment...
            #Should we look at the alphabet.contains() method?
        if record.seq is None :
            raise TypeError("SeqRecord (id=%s) has None for its sequence." % record.id)
            
        #This is abusing the "private" records list,
        #we should really have a method like add_sequence
        #but which takes SeqRecord objects.  See also Bug 1944
        alignment._records.append(record)
    return alignment

def convert(in_file, in_format, out_file, out_format, alphabet=None) :
    """Convert between two sequence file formats, return number of records.

     - in_file - an input handle or filename
     - in_format - input file format, lower case string
     - out_file - an output handle or filename
     - out_format - output file format, lower case string
     - alphabet - optional alphabet to assume

    NOTE - If you provide an output filename, it will be opened which will
    overwrite any existing file without warning. This may happen if even
    the conversion is aborted (e.g. an invalid out_format name is given).

    For example, going from a filename to a handle:

    >>> from Bio import SeqIO
    >>> from StringIO import StringIO
    >>> handle = StringIO("")
    >>> SeqIO.convert("Quality/example.fastq", "fastq", handle, "fasta")
    3
    >>> print handle.getvalue()
    >EAS54_6_R1_2_1_413_324
    CCCTTCTTGTCTTCAGCGTTTCTCC
    >EAS54_6_R1_2_1_540_792
    TTGGCAGGCCAAGGCCGATGGATCA
    >EAS54_6_R1_2_1_443_348
    GTTGCTTCTGGCGTGGGTGGGGGGG
    <BLANKLINE>
    """
    #TODO - Add optimised versions of important conversions
    #For now just off load the work to SeqIO parse/write    
    if isinstance(in_file, basestring) :
        #Hack for SFF, will need to make this more general in future
        if in_format in ["sff", "sff-trim"] :
            in_handle = open(in_file, "rb")
        else :
            in_handle = open(in_file, "rU")
        in_close = True
    else :
        in_handle = in_file
        in_close = False
    #Don't open the output file until we've checked the input is OK?
    if isinstance(out_file, basestring) :
        if out_format in ["sff", "sff_trim"] :
            out_handle = open(out_file, "wb")
        else :
            out_handle = open(out_file, "w")
        out_close = True
    else :
        out_handle = out_file
        out_close = False
    #This will check the arguments and issue error messages,
    #after we have opened the file which is a shame.
    from _convert import _handle_convert #Lazy import
    count = _handle_convert(in_handle, in_format,
                            out_handle, out_format,
                            alphabet)
    #Must now close any handles we opened
    if in_close : in_handle.close()
    if out_close : out_handle.close()
    return count
           
def _test():
    """Run the Bio.SeqIO module's doctests.

    This will try and locate the unit tests directory, and run the doctests
    from there in order that the relative paths used in the examples work.
    """
    import doctest
    import os
    if os.path.isdir(os.path.join("..","..","Tests")) :
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join("..","..","Tests"))
        doctest.testmod()
        os.chdir(cur_dir)
        del cur_dir
        print "Done"
    elif os.path.isdir(os.path.join("Tests", "Fasta")) :
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join("Tests"))
        doctest.testmod()
        os.chdir(cur_dir)
        del cur_dir
        print "Done"

if __name__ == "__main__":
    #Run the doctests
    _test()
