# Copyright 2008-2010 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Multiple sequence alignment input/output as alignment objects.

The Bio.AlignIO interface is deliberately very similar to Bio.SeqIO, and in
fact the two are connected internally.  Both modules use the same set of file
format names (lower case strings).  From the user's perspective, you can read
in a PHYLIP file containing one or more alignments using Bio.AlignIO, or you
can read in the sequences within these alignmenta using Bio.SeqIO.

Bio.AlignIO is also documented at U{http://biopython.org/wiki/AlignIO} and by
a whole chapter in our tutorial:
 - U{http://biopython.org/DIST/docs/tutorial/Tutorial.html}
 - U{http://biopython.org/DIST/docs/tutorial/Tutorial.pdf}

Input
=====
For the typical special case when your file or handle contains one and only
one alignment, use the function Bio.AlignIO.read().  This takes an input file
handle (or in recent versions of Biopython a filename as a string), format
string and optional number of sequences per alignment.  It will return a single
MultipleSeqAlignment object (or raise an exception if there isn't just one
alignment):

    >>> from Bio import AlignIO
    >>> align = AlignIO.read("Phylip/interlaced.phy", "phylip")
    >>> print align
    SingleLetterAlphabet() alignment with 3 rows and 384 columns
    -----MKVILLFVLAVFTVFVSS---------------RGIPPE...I-- CYS1_DICDI
    MAHARVLLLALAVLATAAVAVASSSSFADSNPIRPVTDRAASTL...VAA ALEU_HORVU
    ------MWATLPLLCAGAWLLGV--------PVCGAAELSVNSL...PLV CATH_HUMAN

For the general case, when the handle could contain any number of alignments,
use the function Bio.AlignIO.parse(...) which takes the same arguments, but
returns an iterator giving MultipleSeqAlignment objects (typically used in a
for loop). If you want random access to the alignments by number, turn this
into a list:

    >>> from Bio import AlignIO
    >>> alignments = list(AlignIO.parse("Emboss/needle.txt", "emboss"))
    >>> print alignments[2]
    SingleLetterAlphabet() alignment with 2 rows and 120 columns
    -KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIVTKER...--- ref_rec
    LHIVVVDDDPGTCVYIESVFAELGHTCKSFVRPEAAEEYILTHP...HKE gi|94967506|receiver

Most alignment file formats can be concatenated so as to hold as many
different multiple sequence alignments as possible.  One common example
is the output of the tool seqboot in the PHLYIP suite.  Sometimes there
can be a file header and footer, as seen in the EMBOSS alignment output.

Output
======
Use the function Bio.AlignIO.write(...), which takes a complete set of
Alignment objects (either as a list, or an iterator), an output file handle
(or filename in recent versions of Biopython) and of course the file format::

    from Bio import AlignIO
    alignments = ...
    count = SeqIO.write(alignments, "example.faa", "fasta")

If using a handle make sure to close it to flush the data to the disk::

    from Bio import AlignIO
    alignments = ...
    handle = open("example.faa", "w")
    count = SeqIO.write(alignments, handle, "fasta")
    handle.close()

In general, you are expected to call this function once (with all your
alignments) and then close the file handle.  However, for file formats
like PHYLIP where multiple alignments are stored sequentially (with no file
header and footer), then multiple calls to the write function should work as
expected when using handles.

If you are using a filename, the repeated calls to the write functions will
overwrite the existing file each time.

Conversion
==========
The Bio.AlignIO.convert(...) function allows an easy interface for simple
alignnment file format conversions. Additionally, it may use file format
specific optimisations so this should be the fastest way too.

In general however, you can combine the Bio.AlignIO.parse(...) function with
the Bio.AlignIO.write(...) function for sequence file conversion. Using
generator expressions provides a memory efficient way to perform filtering or
other extra operations as part of the process.

File Formats
============
When specifying the file format, use lowercase strings.  The same format
names are also used in Bio.SeqIO and include the following:

 - clustal   - Ouput from Clustal W or X, see also the module Bio.Clustalw
               which can be used to run the command line tool from Biopython.
 - emboss    - EMBOSS tools' "pairs" and "simple" alignment formats.
 - fasta     - The generic sequence file format where each record starts with
               an identifer line starting with a ">" character, followed by
               lines of sequence.
 - fasta-m10 - For the pairswise alignments output by Bill Pearson's FASTA
               tools when used with the -m 10 command line option for machine
               readable output.
 - ig        - The IntelliGenetics file format, apparently the same as the
               MASE alignment format.
 - nexus     - Output from NEXUS, see also the module Bio.Nexus which can also
               read any phylogenetic trees in these files.
 - phylip    - Used by the PHLIP tools.
 - stockholm - A richly annotated alignment file format used by PFAM.

Note that while Bio.AlignIO can read all the above file formats, it cannot
write to all of them.

You can also use any file format supported by Bio.SeqIO, such as "fasta" or
"ig" (which are listed above), PROVIDED the sequences in your file are all the
same length.
"""

# For using with statement in Python 2.5 or Jython
from __future__ import with_statement

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

from Bio.Align import MultipleSeqAlignment
from Bio.Align.Generic import Alignment
from Bio.Alphabet import Alphabet, AlphabetEncoder, _get_base_alphabet
from Bio.File import as_handle

import StockholmIO
import ClustalIO
import NexusIO
import PhylipIO
import EmbossIO
import FastaIO

#Convention for format names is "mainname-subtype" in lower case.
#Please use the same names as BioPerl and EMBOSS where possible.

_FormatToIterator = {#"fasta" is done via Bio.SeqIO
                     "clustal" : ClustalIO.ClustalIterator,
                     "emboss" : EmbossIO.EmbossIterator,
                     "fasta-m10" : FastaIO.FastaM10Iterator,
                     "nexus" : NexusIO.NexusIterator,
                     "phylip" : PhylipIO.PhylipIterator,
                     "phylip-sequential" : PhylipIO.SequentialPhylipIterator,
                     "phylip-relaxed" : PhylipIO.RelaxedPhylipIterator,
                     "stockholm" : StockholmIO.StockholmIterator,
                     }

_FormatToWriter = {#"fasta" is done via Bio.SeqIO
                   #"emboss" : EmbossIO.EmbossWriter, (unfinished)
                   "nexus" : NexusIO.NexusWriter,
                   "phylip" : PhylipIO.PhylipWriter,
                   "phylip-sequential" : PhylipIO.SequentialPhylipWriter,
                   "phylip-relaxed" : PhylipIO.RelaxedPhylipWriter,
                   "stockholm" : StockholmIO.StockholmWriter,
                   "clustal" : ClustalIO.ClustalWriter,
                   }

def write(alignments, handle, format):
    """Write complete set of alignments to a file.

    Arguments:
     - alignments - A list (or iterator) of Alignment objects (ideally the
                   new MultipleSeqAlignment objects), or (if using Biopython
                   1.54 or later) a single alignment object.
     - handle    - File handle object to write to, or filename as string
                   (note older versions of Biopython only took a handle).
     - format    - lower case string describing the file format to write.

    You should close the handle after calling this function.

    Returns the number of alignments written (as an integer).
    """
    from Bio import SeqIO

    #Try and give helpful error messages:
    if not isinstance(format, basestring):
        raise TypeError("Need a string for the file format (lower case)")
    if not format:
        raise ValueError("Format required (lower case string)")
    if format != format.lower():
        raise ValueError("Format string '%s' should be lower case" % format)

    if isinstance(alignments, Alignment):
        #This raised an exception in older version of Biopython
        alignments = [alignments]

    with as_handle(handle, 'w') as fp:
        #Map the file format to a writer class
        if format in _FormatToWriter:
            writer_class = _FormatToWriter[format]
            count = writer_class(fp).write_file(alignments)
        elif format in SeqIO._FormatToWriter:
            #Exploit the existing SeqIO parser to the dirty work!
            #TODO - Can we make one call to SeqIO.write() and count the alignments?
            count = 0
            for alignment in alignments:
                if not isinstance(alignment, Alignment):
                    raise TypeError(\
                        "Expect a list or iterator of Alignment objects.")
                SeqIO.write(alignment, fp, format)
                count += 1
        elif format in _FormatToIterator or format in SeqIO._FormatToIterator:
            raise ValueError("Reading format '%s' is supported, but not writing" \
                             % format)
        else:
            raise ValueError("Unknown format '%s'" % format)

    assert isinstance(count, int), "Internal error - the underlying %s " \
           "writer should have returned the alignment count, not %s" \
           % (format, repr(count))

    return count

#This is a generator function!
def _SeqIO_to_alignment_iterator(handle, format, alphabet=None, seq_count=None):
    """Uses Bio.SeqIO to create an MultipleSeqAlignment iterator (PRIVATE).

    Arguments:
     - handle    - handle to the file.
     - format    - string describing the file format.
     - alphabet  - optional Alphabet object, useful when the sequence type
                   cannot be automatically inferred from the file itself
                   (e.g. fasta, phylip, clustal)
     - seq_count - Optional integer, number of sequences expected in each
                   alignment.  Recommended for fasta format files.

    If count is omitted (default) then all the sequences in the file are
    combined into a single MultipleSeqAlignment.
    """
    from Bio import SeqIO
    assert format in SeqIO._FormatToIterator

    if seq_count:
        #Use the count to split the records into batches.
        seq_record_iterator = SeqIO.parse(handle, format, alphabet)

        records = []
        for record in seq_record_iterator:
            records.append(record)
            if len(records) == seq_count:
                yield MultipleSeqAlignment(records, alphabet)
                records = []
        if len(records) > 0:
            raise ValueError("Check seq_count argument, not enough sequences?")
    else:
        #Must assume that there is a single alignment using all
        #the SeqRecord objects:
        records = list(SeqIO.parse(handle, format, alphabet))
        if records:
            yield MultipleSeqAlignment(records, alphabet)
    raise StopIteration

def _force_alphabet(alignment_iterator, alphabet):
    """Iterate over alignments, over-riding the alphabet (PRIVATE)."""
    #Assume the alphabet argument has been pre-validated
    given_base_class = _get_base_alphabet(alphabet).__class__
    for align in alignment_iterator:
        if not isinstance(_get_base_alphabet(align._alphabet),
                          given_base_class):
            raise ValueError("Specified alphabet %s clashes with "\
                             "that determined from the file, %s" \
                             % (repr(alphabet), repr(align._alphabet)))
        for record in align:
            if not isinstance(_get_base_alphabet(record.seq.alphabet),
                              given_base_class):
                raise ValueError("Specified alphabet %s clashes with "\
                                 "that determined from the file, %s" \
                           % (repr(alphabet), repr(record.seq.alphabet)))
            record.seq.alphabet = alphabet
        align._alphabet = alphabet
        yield align

def parse(handle, format, seq_count=None, alphabet=None):
    """Iterate over an alignment file as MultipleSeqAlignment objects.

    Arguments:
     - handle    - handle to the file, or the filename as a string
                   (note older verions of Biopython only took a handle).
     - format    - string describing the file format.
     - alphabet  - optional Alphabet object, useful when the sequence type
                   cannot be automatically inferred from the file itself
                   (e.g. fasta, phylip, clustal)
     - seq_count - Optional integer, number of sequences expected in each
                   alignment.  Recommended for fasta format files.

    If you have the file name in a string 'filename', use:

    >>> from Bio import AlignIO
    >>> filename = "Emboss/needle.txt"
    >>> format = "emboss"
    >>> for alignment in AlignIO.parse(filename, format):
    ...     print "Alignment of length", alignment.get_alignment_length()
    Alignment of length 124
    Alignment of length 119
    Alignment of length 120
    Alignment of length 118
    Alignment of length 125

    If you have a string 'data' containing the file contents, use:

    from Bio import AlignIO
    from StringIO import StringIO
    my_iterator = AlignIO.parse(StringIO(data), format)

    Use the Bio.AlignIO.read() function when you expect a single record only.
    """
    from Bio import SeqIO

    #Try and give helpful error messages:
    if not isinstance(format, basestring):
        raise TypeError("Need a string for the file format (lower case)")
    if not format:
        raise ValueError("Format required (lower case string)")
    if format != format.lower():
        raise ValueError("Format string '%s' should be lower case" % format)
    if alphabet is not None and not (isinstance(alphabet, Alphabet) or \
                                     isinstance(alphabet, AlphabetEncoder)):
        raise ValueError("Invalid alphabet, %s" % repr(alphabet))
    if seq_count is not None and not isinstance(seq_count, int):
        raise TypeError("Need integer for seq_count (sequences per alignment)")

    with as_handle(handle, 'rU') as fp:
        #Map the file format to a sequence iterator:
        if format in _FormatToIterator:
            iterator_generator = _FormatToIterator[format]
            if alphabet is None :
                i = iterator_generator(fp, seq_count)
            else:
                try:
                    #Initially assume the optional alphabet argument is supported
                    i = iterator_generator(fp, seq_count, alphabet=alphabet)
                except TypeError:
                    #It isn't supported.
                    i = _force_alphabet(iterator_generator(fp, seq_count),
                                        alphabet)

        elif format in SeqIO._FormatToIterator:
            #Exploit the existing SeqIO parser to the dirty work!
            i = _SeqIO_to_alignment_iterator(fp, format,
                                                alphabet=alphabet,
                                                seq_count=seq_count)
        else:
            raise ValueError("Unknown format '%s'" % format)

        #This imposes some overhead... wait until we drop Python 2.4 to fix it
        for a in i:
            yield a

def read(handle, format, seq_count=None, alphabet=None):
    """Turns an alignment file into a single MultipleSeqAlignment object.

    Arguments:
     - handle    - handle to the file, or the filename as a string
                   (note older verions of Biopython only took a handle).
     - format    - string describing the file format.
     - alphabet  - optional Alphabet object, useful when the sequence type
                   cannot be automatically inferred from the file itself
                   (e.g. fasta, phylip, clustal)
     - seq_count - Optional integer, number of sequences expected in each
                   alignment.  Recommended for fasta format files.

    If the handle contains no alignments, or more than one alignment,
    an exception is raised.  For example, using a PFAM/Stockholm file
    containing one alignment:

    >>> from Bio import AlignIO
    >>> filename = "Clustalw/protein.aln"
    >>> format = "clustal"
    >>> alignment = AlignIO.read(filename, format)
    >>> print "Alignment of length", alignment.get_alignment_length()
    Alignment of length 411

    If however you want the first alignment from a file containing
    multiple alignments this function would raise an exception.

    >>> from Bio import AlignIO
    >>> filename = "Emboss/needle.txt"
    >>> format = "emboss"
    >>> alignment = AlignIO.read(filename, format)
    Traceback (most recent call last):
        ...
    ValueError: More than one record found in handle

    Instead use:

    >>> from Bio import AlignIO
    >>> filename = "Emboss/needle.txt"
    >>> format = "emboss"
    >>> alignment = AlignIO.parse(filename, format).next()
    >>> print "First alignment has length", alignment.get_alignment_length()
    First alignment has length 124

    You must use the Bio.AlignIO.parse() function if you want to read multiple
    records from the handle.
    """
    iterator = parse(handle, format, seq_count, alphabet)
    try:
        first = iterator.next()
    except StopIteration:
        first = None
    if first is None:
        raise ValueError("No records found in handle")
    try:
        second = iterator.next()
    except StopIteration:
        second = None
    if second is not None:
        raise ValueError("More than one record found in handle")
    if seq_count:
        assert len(first)==seq_count
    return first

def convert(in_file, in_format, out_file, out_format, alphabet=None):
    """Convert between two alignment files, returns number of alignments.

     - in_file - an input handle or filename
     - in_format - input file format, lower case string
     - output - an output handle or filename
     - out_file - output file format, lower case string
     - alphabet - optional alphabet to assume

    NOTE - If you provide an output filename, it will be opened which will
    overwrite any existing file without warning. This may happen if even the
    conversion is aborted (e.g. an invalid out_format name is given).
    """
    #TODO - Add optimised versions of important conversions
    #For now just off load the work to SeqIO parse/write
    with as_handle(in_file, 'rU') as in_handle:
        #Don't open the output file until we've checked the input is OK:
        alignments = parse(in_handle, in_format, None, alphabet)

        #This will check the arguments and issue error messages,
        #after we have opened the file which is a shame.
        with as_handle(out_file, 'w') as out_handle:
            count = write(alignments, out_handle, out_format)

    return count

def _test():
    """Run the Bio.AlignIO module's doctests.

    This will try and locate the unit tests directory, and run the doctests
    from there in order that the relative paths used in the examples work.
    """
    import doctest
    import os
    if os.path.isdir(os.path.join("..", "..", "Tests")):
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join("..", "..", "Tests"))
        doctest.testmod()
        os.chdir(cur_dir)
        del cur_dir
        print "Done"
    elif os.path.isdir(os.path.join("Tests", "Fasta")):
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join("Tests"))
        doctest.testmod()
        os.chdir(cur_dir)
        del cur_dir
        print "Done"

if __name__ == "__main__":
    _test()
