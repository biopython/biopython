# Copyright 2008-2018 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Multiple sequence alignment input/output as alignment objects.

The Bio.AlignIO interface is deliberately very similar to Bio.SeqIO, and in
fact the two are connected internally.  Both modules use the same set of file
format names (lower case strings).  From the user's perspective, you can read
in a PHYLIP file containing one or more alignments using Bio.AlignIO, or you
can read in the sequences within these alignments using Bio.SeqIO.

Bio.AlignIO is also documented at http://biopython.org/wiki/AlignIO and by
a whole chapter in our tutorial:

* `HTML Tutorial`_
* `PDF Tutorial`_

.. _`HTML Tutorial`: http://biopython.org/DIST/docs/tutorial/Tutorial.html
.. _`PDF Tutorial`: http://biopython.org/DIST/docs/tutorial/Tutorial.pdf

Input
-----
For the typical special case when your file or handle contains one and only
one alignment, use the function Bio.AlignIO.read().  This takes an input file
handle (or in recent versions of Biopython a filename as a string), format
string and optional number of sequences per alignment.  It will return a single
MultipleSeqAlignment object (or raise an exception if there isn't just one
alignment):

>>> from Bio import AlignIO
>>> align = AlignIO.read("Phylip/interlaced.phy", "phylip")
>>> print(align)
Alignment with 3 rows and 384 columns
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
>>> print(alignments[2])
Alignment with 2 rows and 120 columns
-KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIVTKER...--- ref_rec
LHIVVVDDDPGTCVYIESVFAELGHTCKSFVRPEAAEEYILTHP...HKE gi|94967506|receiver

Most alignment file formats can be concatenated so as to hold as many
different multiple sequence alignments as possible.  One common example
is the output of the tool seqboot in the PHLYIP suite.  Sometimes there
can be a file header and footer, as seen in the EMBOSS alignment output.

Output
------
Use the function Bio.AlignIO.write(...), which takes a complete set of
Alignment objects (either as a list, or an iterator), an output file handle
(or filename in recent versions of Biopython) and of course the file format::

    from Bio import AlignIO
    alignments = ...
    count = SeqIO.write(alignments, "example.faa", "fasta")

If using a handle make sure to close it to flush the data to the disk::

    from Bio import AlignIO
    alignments = ...
    with open("example.faa", "w") as handle:
        count = SeqIO.write(alignments, handle, "fasta")

In general, you are expected to call this function once (with all your
alignments) and then close the file handle.  However, for file formats
like PHYLIP where multiple alignments are stored sequentially (with no file
header and footer), then multiple calls to the write function should work as
expected when using handles.

If you are using a filename, the repeated calls to the write functions will
overwrite the existing file each time.

Conversion
----------
The Bio.AlignIO.convert(...) function allows an easy interface for simple
alignment file format conversions. Additionally, it may use file format
specific optimisations so this should be the fastest way too.

In general however, you can combine the Bio.AlignIO.parse(...) function with
the Bio.AlignIO.write(...) function for sequence file conversion. Using
generator expressions provides a memory efficient way to perform filtering or
other extra operations as part of the process.

File Formats
------------
When specifying the file format, use lowercase strings.  The same format
names are also used in Bio.SeqIO and include the following:

  - clustal -   Output from Clustal W or X, see also the module Bio.Clustalw
    which can be used to run the command line tool from Biopython.
  - emboss    - EMBOSS tools' "pairs" and "simple" alignment formats.
  - fasta     - The generic sequence file format where each record starts with
    an identifier line starting with a ">" character, followed by
    lines of sequence.
  - fasta-m10 - For the pairwise alignments output by Bill Pearson's FASTA
    tools when used with the -m 10 command line option for machine
    readable output.
  - ig        - The IntelliGenetics file format, apparently the same as the
    MASE alignment format.
  - msf       - The GCG MSF alignment format, originally from PileUp tool.
  - nexus     - Output from NEXUS, see also the module Bio.Nexus which can also
    read any phylogenetic trees in these files.
  - phylip    - Interlaced PHYLIP, as used by the PHYLIP tools.
  - phylip-sequential - Sequential PHYLIP.
  - phylip-relaxed - PHYLIP like format allowing longer names.
  - stockholm - A richly annotated alignment file format used by PFAM.
  - mauve - Output from progressiveMauve/Mauve

Note that while Bio.AlignIO can read all the above file formats, it cannot
write to all of them.

You can also use any file format supported by Bio.SeqIO, such as "fasta" or
"ig" (which are listed above), PROVIDED the sequences in your file are all the
same length.
"""


# TODO
# - define policy on reading aligned sequences with gaps in
#   (e.g. - and . characters)
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
from Bio.File import as_handle

from . import StockholmIO
from . import ClustalIO
from . import NexusIO
from . import PhylipIO
from . import EmbossIO
from . import FastaIO
from . import MafIO
from . import MauveIO
from . import MsfIO

# Convention for format names is "mainname-subtype" in lower case.
# Please use the same names as BioPerl and EMBOSS where possible.

_FormatToIterator = {  # "fasta" is done via Bio.SeqIO
    "clustal": ClustalIO.ClustalIterator,
    "emboss": EmbossIO.EmbossIterator,
    "fasta-m10": FastaIO.FastaM10Iterator,
    "maf": MafIO.MafIterator,
    "mauve": MauveIO.MauveIterator,
    "msf": MsfIO.MsfIterator,
    "nexus": NexusIO.NexusIterator,
    "phylip": PhylipIO.PhylipIterator,
    "phylip-sequential": PhylipIO.SequentialPhylipIterator,
    "phylip-relaxed": PhylipIO.RelaxedPhylipIterator,
    "stockholm": StockholmIO.StockholmIterator,
}

_FormatToWriter = {  # "fasta" is done via Bio.SeqIO
    "clustal": ClustalIO.ClustalWriter,
    "maf": MafIO.MafWriter,
    "mauve": MauveIO.MauveWriter,
    "nexus": NexusIO.NexusWriter,
    "phylip": PhylipIO.PhylipWriter,
    "phylip-sequential": PhylipIO.SequentialPhylipWriter,
    "phylip-relaxed": PhylipIO.RelaxedPhylipWriter,
    "stockholm": StockholmIO.StockholmWriter,
}


def write(alignments, handle, format):
    """Write complete set of alignments to a file.

    Arguments:
     - alignments - A list (or iterator) of MultipleSeqAlignment objects,
       or a single alignment object.
     - handle    - File handle object to write to, or filename as string
       (note older versions of Biopython only took a handle).
     - format    - lower case string describing the file format to write.

    You should close the handle after calling this function.

    Returns the number of alignments written (as an integer).
    """
    from Bio import SeqIO

    # Try and give helpful error messages:
    if not isinstance(format, str):
        raise TypeError("Need a string for the file format (lower case)")
    if not format:
        raise ValueError("Format required (lower case string)")
    if format != format.lower():
        raise ValueError("Format string '%s' should be lower case" % format)

    if isinstance(alignments, MultipleSeqAlignment):
        # This raised an exception in older versions of Biopython
        alignments = [alignments]

    with as_handle(handle, "w") as fp:
        # Map the file format to a writer class
        if format in _FormatToWriter:
            writer_class = _FormatToWriter[format]
            count = writer_class(fp).write_file(alignments)
        elif format in SeqIO._FormatToWriter:
            # Exploit the existing SeqIO parser to do the dirty work!
            # TODO - Can we make one call to SeqIO.write() and count the alignments?
            count = 0
            for alignment in alignments:
                if not isinstance(alignment, MultipleSeqAlignment):
                    raise TypeError(
                        "Expect a list or iterator of MultipleSeqAlignment "
                        "objects, got: %r" % alignment
                    )
                SeqIO.write(alignment, fp, format)
                count += 1
        elif format in _FormatToIterator or format in SeqIO._FormatToIterator:
            raise ValueError(
                "Reading format '%s' is supported, but not writing" % format
            )
        else:
            raise ValueError("Unknown format '%s'" % format)

    if not isinstance(count, int):
        raise RuntimeError(
            "Internal error - the underlying %s "
            "writer should have returned the alignment count, not %r" % (format, count)
        )

    return count


# This is a generator function!
def _SeqIO_to_alignment_iterator(handle, format, seq_count=None):
    """Use Bio.SeqIO to create an MultipleSeqAlignment iterator (PRIVATE).

    Arguments:
     - handle    - handle to the file.
     - format    - string describing the file format.
     - seq_count - Optional integer, number of sequences expected in each
       alignment.  Recommended for fasta format files.

    If count is omitted (default) then all the sequences in the file are
    combined into a single MultipleSeqAlignment.
    """
    from Bio import SeqIO

    if format not in SeqIO._FormatToIterator:
        raise ValueError("Unknown format '%s'" % format)

    if seq_count:
        # Use the count to split the records into batches.
        seq_record_iterator = SeqIO.parse(handle, format)

        records = []
        for record in seq_record_iterator:
            records.append(record)
            if len(records) == seq_count:
                yield MultipleSeqAlignment(records)
                records = []
        if records:
            raise ValueError("Check seq_count argument, not enough sequences?")
    else:
        # Must assume that there is a single alignment using all
        # the SeqRecord objects:
        records = list(SeqIO.parse(handle, format))
        if records:
            yield MultipleSeqAlignment(records)


def parse(handle, format, seq_count=None):
    """Iterate over an alignment file as MultipleSeqAlignment objects.

    Arguments:
     - handle    - handle to the file, or the filename as a string
       (note older versions of Biopython only took a handle).
     - format    - string describing the file format.
     - seq_count - Optional integer, number of sequences expected in each
       alignment.  Recommended for fasta format files.

    If you have the file name in a string 'filename', use:

    >>> from Bio import AlignIO
    >>> filename = "Emboss/needle.txt"
    >>> format = "emboss"
    >>> for alignment in AlignIO.parse(filename, format):
    ...     print("Alignment of length %i" % alignment.get_alignment_length())
    Alignment of length 124
    Alignment of length 119
    Alignment of length 120
    Alignment of length 118
    Alignment of length 125

    If you have a string 'data' containing the file contents, use::

      from Bio import AlignIO
      from io import StringIO
      my_iterator = AlignIO.parse(StringIO(data), format)

    Use the Bio.AlignIO.read() function when you expect a single record only.
    """
    from Bio import SeqIO

    # Try and give helpful error messages:
    if not isinstance(format, str):
        raise TypeError("Need a string for the file format (lower case)")
    if not format:
        raise ValueError("Format required (lower case string)")
    if format != format.lower():
        raise ValueError("Format string '%s' should be lower case" % format)
    if seq_count is not None and not isinstance(seq_count, int):
        raise TypeError("Need integer for seq_count (sequences per alignment)")

    with as_handle(handle) as fp:
        # Map the file format to a sequence iterator:
        if format in _FormatToIterator:
            iterator_generator = _FormatToIterator[format]
            i = iterator_generator(fp, seq_count)

        elif format in SeqIO._FormatToIterator:
            # Exploit the existing SeqIO parser to the dirty work!
            i = _SeqIO_to_alignment_iterator(fp, format, seq_count=seq_count)
        else:
            raise ValueError("Unknown format '%s'" % format)

        yield from i


def read(handle, format, seq_count=None):
    """Turn an alignment file into a single MultipleSeqAlignment object.

    Arguments:
     - handle    - handle to the file, or the filename as a string
       (note older versions of Biopython only took a handle).
     - format    - string describing the file format.
     - seq_count - Optional integer, number of sequences expected in each
       alignment.  Recommended for fasta format files.

    If the handle contains no alignments, or more than one alignment,
    an exception is raised.  For example, using a PFAM/Stockholm file
    containing one alignment:

    >>> from Bio import AlignIO
    >>> filename = "Clustalw/protein.aln"
    >>> format = "clustal"
    >>> alignment = AlignIO.read(filename, format)
    >>> print("Alignment of length %i" % alignment.get_alignment_length())
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
    >>> alignment = next(AlignIO.parse(filename, format))
    >>> print("First alignment has length %i" % alignment.get_alignment_length())
    First alignment has length 124

    You must use the Bio.AlignIO.parse() function if you want to read multiple
    records from the handle.
    """
    iterator = parse(handle, format, seq_count)
    try:
        alignment = next(iterator)
    except StopIteration:
        raise ValueError("No records found in handle") from None
    try:
        next(iterator)
        raise ValueError("More than one record found in handle")
    except StopIteration:
        pass
    if seq_count:
        if len(alignment) != seq_count:
            raise RuntimeError(
                "More sequences found in alignment than specified in seq_count: %s."
                % seq_count
            )
    return alignment


def convert(in_file, in_format, out_file, out_format, molecule_type=None):
    """Convert between two alignment files, returns number of alignments.

    Arguments:
     - in_file - an input handle or filename
     - in_format - input file format, lower case string
     - output - an output handle or filename
     - out_file - output file format, lower case string
     - molecule_type - optional molecule type to apply, string containing
       "DNA", "RNA" or "protein".

    **NOTE** - If you provide an output filename, it will be opened which will
    overwrite any existing file without warning. This may happen if even the
    conversion is aborted (e.g. an invalid out_format name is given).

    Some output formats require the molecule type be specified where this
    cannot be determined by the parser. For example, converting to FASTA,
    Clustal, or PHYLIP format to NEXUS:

    >>> from io import StringIO
    >>> from Bio import AlignIO
    >>> handle = StringIO()
    >>> AlignIO.convert("Phylip/horses.phy", "phylip", handle, "nexus", "DNA")
    1
    >>> print(handle.getvalue())
    #NEXUS
    begin data;
    dimensions ntax=10 nchar=40;
    format datatype=dna missing=? gap=-;
    matrix
    Mesohippus   AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    Hypohippus   AAACCCCCCCAAAAAAAAACAAAAAAAAAAAAAAAAAAAA
    Archaeohip   CAAAAAAAAAAAAAAAACACAAAAAAAAAAAAAAAAAAAA
    Parahippus   CAAACAACAACAAAAAAAACAAAAAAAAAAAAAAAAAAAA
    Merychippu   CCAACCACCACCCCACACCCAAAAAAAAAAAAAAAAAAAA
    'M. secundu' CCAACCACCACCCACACCCCAAAAAAAAAAAAAAAAAAAA
    Nannipus     CCAACCACAACCCCACACCCAAAAAAAAAAAAAAAAAAAA
    Neohippari   CCAACCCCCCCCCCACACCCAAAAAAAAAAAAAAAAAAAA
    Calippus     CCAACCACAACCCACACCCCAAAAAAAAAAAAAAAAAAAA
    Pliohippus   CCCACCCCCCCCCACACCCCAAAAAAAAAAAAAAAAAAAA
    ;
    end;
    <BLANKLINE>
    """
    if molecule_type:
        if not isinstance(molecule_type, str):
            raise TypeError("Molecule type should be a string, not %r" % molecule_type)
        elif (
            "DNA" in molecule_type
            or "RNA" in molecule_type
            or "protein" in molecule_type
        ):
            pass
        else:
            raise ValueError("Unexpected molecule type, %r" % molecule_type)

    # TODO - Add optimised versions of important conversions
    # For now just off load the work to SeqIO parse/write
    # Don't open the output file until we've checked the input is OK:
    alignments = parse(in_file, in_format, None)

    if molecule_type:
        # Edit the records on the fly to set molecule type

        def over_ride(alignment):
            """Over-ride molecule in-place."""
            for record in alignment:
                record.annotations["molecule_type"] = molecule_type
            return alignment

        alignments = (over_ride(_) for _ in alignments)
    return write(alignments, out_file, out_format)


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
