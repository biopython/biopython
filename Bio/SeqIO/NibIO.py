# Copyright 2019 by Michiel de Hoon.  All rights reserved.
# Based on code contributed and copyright 2016 by Peter Cock.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Bio.SeqIO support for the binary Standard Flowgram Format (Nib) file format.

Nib was designed by 454 Life Sciences (Roche), the Whitehead Institute for
Biomedical Research and the Wellcome Trust Sanger Institute. Nib was also used
as the native output format from early versions of Ion Torrent's PGM platform
as well. You are expected to use this module via the Bio.SeqIO functions under
the format name "sff" (or "sff-trim" as described below).

For example, to iterate over the records in an Nib file,

    >>> from Bio import SeqIO
    >>> for record in SeqIO.parse("Roche/E3MFGYR02_random_10_reads.sff", "sff"):
    ...     print("%s %i %s..." % (record.id, len(record), record.seq[:20]))
    ...
    E3MFGYR02JWQ7T 265 tcagGGTCTACATGTTGGTT...
    E3MFGYR02JA6IL 271 tcagTTTTTTTTGGAAAGGA...
    E3MFGYR02JHD4H 310 tcagAAAGACAAGTGGTATC...
    E3MFGYR02GFKUC 299 tcagCGGCCGGGCCTCTCAT...
    E3MFGYR02FTGED 281 tcagTGGTAATGGGGGGAAA...
    E3MFGYR02FR9G7 261 tcagCTCCGTAAGAAGGTGC...
    E3MFGYR02GAZMS 278 tcagAAAGAAGTAAGGTAAA...
    E3MFGYR02HHZ8O 221 tcagACTTTCTTCTTTACCG...
    E3MFGYR02GPGB1 269 tcagAAGCAGTGGTATCAAC...
    E3MFGYR02F7Z7G 219 tcagAATCATCCACTTTTTA...

Each SeqRecord object will contain all the annotation from the Nib file,
including the PHRED quality scores.

    >>> print("%s %i" % (record.id, len(record)))
    E3MFGYR02F7Z7G 219
    >>> print("%s..." % record.seq[:10])
    tcagAATCAT...
    >>> print("%r..." % (record.letter_annotations["phred_quality"][:10]))
    [22, 21, 23, 28, 26, 15, 12, 21, 28, 21]...

Notice that the sequence is given in mixed case, the central upper case region
corresponds to the trimmed sequence. This matches the output of the Roche
tools (and the 3rd party tool sff_extract) for Nib to FASTA.

    >>> print(record.annotations["clip_qual_left"])
    4
    >>> print(record.annotations["clip_qual_right"])
    134
    >>> print(record.seq[:4])
    tcag
    >>> print("%s...%s" % (record.seq[4:20], record.seq[120:134]))
    AATCATCCACTTTTTA...CAAAACACAAACAG
    >>> print(record.seq[134:])
    atcttatcaacaaaactcaaagttcctaactgagacacgcaacaggggataagacaaggcacacaggggataggnnnnnnnnnnn

The annotations dictionary also contains any adapter clip positions
(usually zero), and information about the flows. e.g.

    >>> len(record.annotations)
    11
    >>> print(record.annotations["flow_key"])
    TCAG
    >>> print(record.annotations["flow_values"][:10])
    (83, 1, 128, 7, 4, 84, 6, 106, 3, 172)
    >>> print(len(record.annotations["flow_values"]))
    400
    >>> print(record.annotations["flow_index"][:10])
    (1, 2, 3, 2, 2, 0, 3, 2, 3, 3)
    >>> print(len(record.annotations["flow_index"]))
    219

Note that to convert from a raw reading in flow_values to the corresponding
homopolymer stretch estimate, the value should be rounded to the nearest 100:

    >>> print("%r..." % [int(round(value, -2)) // 100
    ...                  for value in record.annotations["flow_values"][:10]])
    ...
    [1, 0, 1, 0, 0, 1, 0, 1, 0, 2]...

If a read name is exactly 14 alphanumeric characters, the annotations
dictionary will also contain meta-data about the read extracted by
interpretting the name as a 454 Sequencing System "Universal" Accession
Number. Note that if a read name happens to be exactly 14 alphanumeric
characters but was not generated automatically, these annotation records
will contain nonsense information.

    >>> print(record.annotations["region"])
    2
    >>> print(record.annotations["time"])
    [2008, 1, 9, 16, 16, 0]
    >>> print(record.annotations["coords"])
    (2434, 1658)

As a convenience method, you can read the file with SeqIO format name "sff-trim"
instead of "sff" to get just the trimmed sequences (without any annotation
except for the PHRED quality scores and anything encoded in the read names):

    >>> from Bio import SeqIO
    >>> for record in SeqIO.parse("Roche/E3MFGYR02_random_10_reads.sff", "sff-trim"):
    ...     print("%s %i %s..." % (record.id, len(record), record.seq[:20]))
    ...
    E3MFGYR02JWQ7T 260 GGTCTACATGTTGGTTAACC...
    E3MFGYR02JA6IL 265 TTTTTTTTGGAAAGGAAAAC...
    E3MFGYR02JHD4H 292 AAAGACAAGTGGTATCAACG...
    E3MFGYR02GFKUC 295 CGGCCGGGCCTCTCATCGGT...
    E3MFGYR02FTGED 277 TGGTAATGGGGGGAAATTTA...
    E3MFGYR02FR9G7 256 CTCCGTAAGAAGGTGCTGCC...
    E3MFGYR02GAZMS 271 AAAGAAGTAAGGTAAATAAC...
    E3MFGYR02HHZ8O 150 ACTTTCTTCTTTACCGTAAC...
    E3MFGYR02GPGB1 221 AAGCAGTGGTATCAACGCAG...
    E3MFGYR02F7Z7G 130 AATCATCCACTTTTTAACGT...

Looking at the final record in more detail, note how this differs to the
example above:

    >>> print("%s %i" % (record.id, len(record)))
    E3MFGYR02F7Z7G 130
    >>> print("%s..." % record.seq[:10])
    AATCATCCAC...
    >>> print("%r..." % record.letter_annotations["phred_quality"][:10])
    [26, 15, 12, 21, 28, 21, 36, 28, 27, 27]...
    >>> len(record.annotations)
    3
    >>> print(record.annotations["region"])
    2
    >>> print(record.annotations["coords"])
    (2434, 1658)
    >>> print(record.annotations["time"])
    [2008, 1, 9, 16, 16, 0]

You might use the Bio.SeqIO.convert() function to convert the (trimmed) Nib
reads into a FASTQ file (or a FASTA file and a QUAL file), e.g.

    >>> from Bio import SeqIO
    >>> try:
    ...     from StringIO import StringIO # Python 2
    ... except ImportError:
    ...     from io import StringIO # Python 3
    ...
    >>> out_handle = StringIO()
    >>> count = SeqIO.convert("Roche/E3MFGYR02_random_10_reads.sff", "sff",
    ...                       out_handle, "fastq")
    ...
    >>> print("Converted %i records" % count)
    Converted 10 records

The output FASTQ file would start like this:

    >>> print("%s..." % out_handle.getvalue()[:50])
    @E3MFGYR02JWQ7T
    tcagGGTCTACATGTTGGTTAACCCGTACTGATT...

Bio.SeqIO.index() provides memory efficient random access to the reads in an
Nib file by name. Nib files can include an index within the file, which can
be read in making this very fast. If the index is missing (or in a format not
yet supported in Biopython) the file is indexed by scanning all the reads -
which is a little slower. For example,

    >>> from Bio import SeqIO
    >>> reads = SeqIO.index("Roche/E3MFGYR02_random_10_reads.sff", "sff")
    >>> record = reads["E3MFGYR02JHD4H"]
    >>> print("%s %i %s..." % (record.id, len(record), record.seq[:20]))
    E3MFGYR02JHD4H 310 tcagAAAGACAAGTGGTATC...
    >>> reads.close()

Or, using the trimmed reads:

    >>> from Bio import SeqIO
    >>> reads = SeqIO.index("Roche/E3MFGYR02_random_10_reads.sff", "sff-trim")
    >>> record = reads["E3MFGYR02JHD4H"]
    >>> print("%s %i %s..." % (record.id, len(record), record.seq[:20]))
    E3MFGYR02JHD4H 292 AAAGACAAGTGGTATCAACG...
    >>> reads.close()

You can also use the Bio.SeqIO.write() function with the "sff" format. Note
that this requires all the flow information etc, and thus is probably only
useful for SeqRecord objects originally from reading another Nib file (and
not the trimmed SeqRecord objects from parsing an Nib file as "sff-trim").

As an example, let's pretend this example Nib file represents some DNA which
was pre-amplified with a PCR primers AAAGANNNNN. The following script would
produce a sub-file containing all those reads whose post-quality clipping
region (i.e. the sequence after trimming) starts with AAAGA exactly (the non-
degenerate bit of this pretend primer):

    >>> from Bio import SeqIO
    >>> records = (record for record in
    ...            SeqIO.parse("Roche/E3MFGYR02_random_10_reads.sff", "sff")
    ...            if record.seq[record.annotations["clip_qual_left"]:].startswith("AAAGA"))
    ...
    >>> count = SeqIO.write(records, "temp_filtered.sff", "sff")
    >>> print("Selected %i records" % count)
    Selected 2 records

Of course, for an assembly you would probably want to remove these primers.
If you want FASTA or FASTQ output, you could just slice the SeqRecord. However,
if you want Nib output we have to preserve all the flow information - the trick
is just to adjust the left clip position!

    >>> from Bio import SeqIO
    >>> def filter_and_trim(records, primer):
    ...     for record in records:
    ...         if record.seq[record.annotations["clip_qual_left"]:].startswith(primer):
    ...             record.annotations["clip_qual_left"] += len(primer)
    ...             yield record
    ...
    >>> records = SeqIO.parse("Roche/E3MFGYR02_random_10_reads.sff", "sff")
    >>> count = SeqIO.write(filter_and_trim(records, "AAAGA"),
    ...                     "temp_filtered.sff", "sff")
    ...
    >>> print("Selected %i records" % count)
    Selected 2 records

We can check the results, note the lower case clipped region now includes the "AAAGA"
sequence:

    >>> for record in SeqIO.parse("temp_filtered.sff", "sff"):
    ...     print("%s %i %s..." % (record.id, len(record), record.seq[:20]))
    ...
    E3MFGYR02JHD4H 310 tcagaaagaCAAGTGGTATC...
    E3MFGYR02GAZMS 278 tcagaaagaAGTAAGGTAAA...
    >>> for record in SeqIO.parse("temp_filtered.sff", "sff-trim"):
    ...     print("%s %i %s..." % (record.id, len(record), record.seq[:20]))
    ...
    E3MFGYR02JHD4H 287 CAAGTGGTATCAACGCAGAG...
    E3MFGYR02GAZMS 266 AGTAAGGTAAATAACAAACG...
    >>> import os
    >>> os.remove("temp_filtered.sff")

For a description of the file format, please see the Roche manuals and:
http://www.ncbi.nlm.nih.gov/Traces/trace.cgi?cmd=show&f=formats&m=doc&s=formats

"""

from __future__ import print_function

from Bio.SeqIO.Interfaces import SequenceWriter
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import struct
import sys

from Bio._py3k import _bytes_to_string, _as_bytes


# This is a generator function!
def NibIterator(handle, alphabet=None):
    """Iterate over Standard Flowgram Format (SFF) reads (as SeqRecord objects).

        - handle - input file in the Nib file format as defibed by UCSC.
          This must be opened in binary mode!
        - alphabet - always ignored.

    The sequence of the resulting SeqRecord object should match the sequence
    generated by the nibFrag utility, except that it will be in upper case,
    while nibFrag uses lower case.

    This function is used internally via the Bio.SeqIO functions:

    >>> from Bio import SeqIO
    >>> record = SeqIO.read("Nib/test_bigendian.nib", "nib")
    >>> print("%s %i" % (record.seq, len(record)))
    ACGTAAACCGTACCCGTANANCANNNNACNANNANCN 37

    You can also call it directly:

    >>> with open("Nib/test_bigendian.nib", "rb") as handle:
    ...     for record in NibIterator(handle):
    ...         print("%s %i" % (record.seq, len(record)))
    ...
    ACGTAAACCGTACCCGTANANCANNNNACNANNANCN 37

    """
    if alphabet is not None:
        raise ValueError("Alphabets are ignored.")
    signature = handle.read(4).hex()
    if signature == '3a3de96b':
        byteorder = 'little' # little-endian
    elif signature == '6be93d3a':
        byteorder = 'big' # big-endian
    else:
        raise ValueError('unexpected signature in Nib header')
    number = handle.read(4)
    length = int.from_bytes(number, byteorder)
    indices = handle.read().hex()
    if length % 2 == 0:
        if len(indices) != length:
            raise ValueError('Unexpected file size')
    elif length % 2 == 1:
        if len(indices) != length + 1:
            raise ValueError('Unexpected file size')
        indices = indices[:length]
    if set(indices) != set('01234'):
        raise ValueError('Unexpected sequence data found in file')
    table = str.maketrans('01234','TCAGN')
    nucleotides = indices.translate(table)
    sequence = Seq(nucleotides)
    record = SeqRecord(sequence)
    yield record


class NibWriter(SequenceWriter):
    """Nib file writer."""

    def __init__(self, handle):
        """Initialize an Nib writer object.

        Arguments:
         - handle - Output handle, in binary write mode.
        """
        self.handle = handle
        byteorder = sys.byteorder
        if byteorder == 'little': # little-endian
            signature = '3a3de96b'
        elif byteorder == 'big': # big-endian
            signature = '6be93d3a'
        else:
            raise RuntimeError('unexpected system byte order %s' % byteorder)
        handle.write(bytes.fromhex(signature))

    def write_file(self, records):
        """Use this to write an entire file containing the given record."""
        count = 0
        for record in records:
            count += 1
        if count == 0:
            raise ValueError("Must have one sequence")
        if count > 1:
            raise ValueError('More than one sequence found')
        handle = self.handle
        sequence = record.seq
        nucleotides = str(sequence)
        length = len(sequence)
        handle.write(struct.pack('i', length))
        table = str.maketrans('TCAGNtcagn', '0123401234')
        padding = length % 2
        suffix = padding * 'T'
        nucleotides += suffix
        indices = nucleotides.translate(table)
        if set(indices) != set('01234'):
            raise ValueError('Sequence should contain A,C,G,T,N,a,c,g,t,n only')
        handle.write(bytes.fromhex(indices))
        return count


if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest(verbose=0)
