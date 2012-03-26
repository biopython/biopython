# Copyright 2009-2010 by Peter Cock.  All rights reserved.
# Based on code contributed and copyright 2009 by Jose Blanca (COMAV-UPV).
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Bio.SeqIO support for the binary Standard Flowgram Format (SFF) file format.

SFF was designed by 454 Life Sciences (Roche), the Whitehead Institute for
Biomedical Research and the Wellcome Trust Sanger Institute. You are expected
to use this module via the Bio.SeqIO functions under the format name "sff" (or
"sff-trim" as described below).

For example, to iterate over the records in an SFF file,

    >>> from Bio import SeqIO
    >>> for record in SeqIO.parse("Roche/E3MFGYR02_random_10_reads.sff", "sff"):
    ...     print record.id, len(record), record.seq[:20]+"..."
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

Each SeqRecord object will contain all the annotation from the SFF file,
including the PHRED quality scores.

    >>> print record.id, len(record)
    E3MFGYR02F7Z7G 219
    >>> print record.seq[:10], "..."
    tcagAATCAT ...
    >>> print record.letter_annotations["phred_quality"][:10], "..."
    [22, 21, 23, 28, 26, 15, 12, 21, 28, 21] ...

Notice that the sequence is given in mixed case, the central upper case region
corresponds to the trimmed sequence. This matches the output of the Roche
tools (and the 3rd party tool sff_extract) for SFF to FASTA.

    >>> print record.annotations["clip_qual_left"]
    4
    >>> print record.annotations["clip_qual_right"]
    134
    >>> print record.seq[:4]
    tcag
    >>> print record.seq[4:20], "...", record.seq[120:134]
    AATCATCCACTTTTTA ... CAAAACACAAACAG
    >>> print record.seq[134:]
    atcttatcaacaaaactcaaagttcctaactgagacacgcaacaggggataagacaaggcacacaggggataggnnnnnnnnnnn

The annotations dictionary also contains any adapter clip positions
(usually zero), and information about the flows. e.g.

    >>> print record.annotations["flow_key"]
    TCAG
    >>> print record.annotations["flow_values"][:10], "..."
    (83, 1, 128, 7, 4, 84, 6, 106, 3, 172) ...
    >>> print len(record.annotations["flow_values"])
    400
    >>> print record.annotations["flow_index"][:10], "..."
    (1, 2, 3, 2, 2, 0, 3, 2, 3, 3) ...
    >>> print len(record.annotations["flow_index"])
    219

Note that to convert from a raw reading in flow_values to the corresponding
homopolymer stretch estimate, the value should be rounded to the nearest 100:

    >>> print [int(round(value, -2)) // 100
    ...        for value in record.annotations["flow_values"][:10]], '...'
    [1, 0, 1, 0, 0, 1, 0, 1, 0, 2] ...

As a convenience method, you can read the file with SeqIO format name "sff-trim"
instead of "sff" to get just the trimmed sequences (without any annotation
except for the PHRED quality scores):

    >>> from Bio import SeqIO
    >>> for record in SeqIO.parse("Roche/E3MFGYR02_random_10_reads.sff", "sff-trim"):
    ...     print record.id, len(record), record.seq[:20]+"..."
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

    >>> print record.id, len(record)
    E3MFGYR02F7Z7G 130
    >>> print record.seq[:10], "..."
    AATCATCCAC ...
    >>> print record.letter_annotations["phred_quality"][:10], "..."
    [26, 15, 12, 21, 28, 21, 36, 28, 27, 27] ...
    >>> print record.annotations
    {}

You might use the Bio.SeqIO.convert() function to convert the (trimmed) SFF
reads into a FASTQ file (or a FASTA file and a QUAL file), e.g.

    >>> from Bio import SeqIO
    >>> from StringIO import StringIO
    >>> out_handle = StringIO()
    >>> count = SeqIO.convert("Roche/E3MFGYR02_random_10_reads.sff", "sff",
    ...                       out_handle, "fastq")
    >>> print "Converted %i records" % count
    Converted 10 records

The output FASTQ file would start like this:

    >>> print "%s..." % out_handle.getvalue()[:50]
    @E3MFGYR02JWQ7T
    tcagGGTCTACATGTTGGTTAACCCGTACTGATT...

Bio.SeqIO.index() provides memory efficient random access to the reads in an
SFF file by name. SFF files can include an index within the file, which can
be read in making this very fast. If the index is missing (or in a format not
yet supported in Biopython) the file is indexed by scanning all the reads -
which is a little slower. For example,

    >>> from Bio import SeqIO
    >>> reads = SeqIO.index("Roche/E3MFGYR02_random_10_reads.sff", "sff")
    >>> record = reads["E3MFGYR02JHD4H"]
    >>> print record.id, len(record), record.seq[:20]+"..."
    E3MFGYR02JHD4H 310 tcagAAAGACAAGTGGTATC...

Or, using the trimmed reads:

    >>> from Bio import SeqIO
    >>> reads = SeqIO.index("Roche/E3MFGYR02_random_10_reads.sff", "sff-trim")
    >>> record = reads["E3MFGYR02JHD4H"]
    >>> print record.id, len(record), record.seq[:20]+"..."
    E3MFGYR02JHD4H 292 AAAGACAAGTGGTATCAACG...

You can also use the Bio.SeqIO.write() function with the "sff" format. Note
that this requires all the flow information etc, and thus is probably only
useful for SeqRecord objects originally from reading another SFF file (and
not the trimmed SeqRecord objects from parsing an SFF file as "sff-trim").

As an example, let's pretend this example SFF file represents some DNA which
was pre-amplified with a PCR primers AAAGANNNNN. The following script would
produce a sub-file containing all those reads whose post-quality clipping
region (i.e. the sequence after trimming) starts with AAAGA exactly (the non-
degenerate bit of this pretend primer):

    >>> from Bio import SeqIO
    >>> records = (record for record in
    ...            SeqIO.parse("Roche/E3MFGYR02_random_10_reads.sff","sff")
    ...            if record.seq[record.annotations["clip_qual_left"]:].startswith("AAAGA"))
    >>> count = SeqIO.write(records, "temp_filtered.sff", "sff")
    >>> print "Selected %i records" % count
    Selected 2 records

Of course, for an assembly you would probably want to remove these primers.
If you want FASTA or FASTQ output, you could just slice the SeqRecord. However,
if you want SFF output we have to preserve all the flow information - the trick
is just to adjust the left clip position!

    >>> from Bio import SeqIO
    >>> def filter_and_trim(records, primer):
    ...     for record in records:
    ...         if record.seq[record.annotations["clip_qual_left"]:].startswith(primer):
    ...             record.annotations["clip_qual_left"] += len(primer)
    ...             yield record
    >>> records = SeqIO.parse("Roche/E3MFGYR02_random_10_reads.sff", "sff")
    >>> count = SeqIO.write(filter_and_trim(records,"AAAGA"),
    ...                     "temp_filtered.sff", "sff")
    >>> print "Selected %i records" % count
    Selected 2 records

We can check the results, note the lower case clipped region now includes the "AAAGA"
sequence:

    >>> for record in SeqIO.parse("temp_filtered.sff", "sff"):
    ...     print record.id, len(record), record.seq[:20]+"..."
    E3MFGYR02JHD4H 310 tcagaaagaCAAGTGGTATC...
    E3MFGYR02GAZMS 278 tcagaaagaAGTAAGGTAAA...
    >>> for record in SeqIO.parse("temp_filtered.sff", "sff-trim"):
    ...     print record.id, len(record), record.seq[:20]+"..."
    E3MFGYR02JHD4H 287 CAAGTGGTATCAACGCAGAG...
    E3MFGYR02GAZMS 266 AGTAAGGTAAATAACAAACG...
    >>> import os
    >>> os.remove("temp_filtered.sff")

For a description of the file format, please see the Roche manuals and:
http://www.ncbi.nlm.nih.gov/Traces/trace.cgi?cmd=show&f=formats&m=doc&s=formats

"""
from Bio.SeqIO.Interfaces import SequenceWriter
from Bio import Alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import struct
import sys

from Bio._py3k import _bytes_to_string, _as_bytes
_null = _as_bytes("\0")
_sff = _as_bytes(".sff")
_hsh = _as_bytes(".hsh")
_srt = _as_bytes(".srt")
_mft = _as_bytes(".mft")
_flag = _as_bytes("\xff")

def _sff_file_header(handle):
    """Read in an SFF file header (PRIVATE).

    Assumes the handle is at the start of the file, will read forwards
    though the header and leave the handle pointing at the first record.
    Returns a tuple of values from the header (header_length, index_offset,
    index_length, number_of_reads, flows_per_read, flow_chars, key_sequence)

    >>> handle = open("Roche/greek.sff", "rb")
    >>> values = _sff_file_header(handle)
    >>> print values[0]
    840
    >>> print values[1]
    65040
    >>> print values[2]
    256
    >>> print values[3]
    24
    >>> print values[4]
    800
    >>> values[-1]
    'TCAG'

    """
    if hasattr(handle,"mode") and "U" in handle.mode.upper():
        raise ValueError("SFF files must NOT be opened in universal new "
                         "lines mode. Binary mode is recommended (although "
                         "on Unix the default mode is also fine).")
    elif hasattr(handle,"mode") and "B" not in handle.mode.upper() \
    and sys.platform == "win32":
        raise ValueError("SFF files must be opened in binary mode on Windows")
    #file header (part one)
    #use big endiean encdoing   >
    #magic_number               I
    #version                    4B
    #index_offset               Q
    #index_length               I
    #number_of_reads            I
    #header_length              H
    #key_length                 H
    #number_of_flows_per_read   H
    #flowgram_format_code       B
    #[rest of file header depends on the number of flows and how many keys]
    fmt = '>4s4BQIIHHHB'
    assert 31 == struct.calcsize(fmt)
    data = handle.read(31)
    if not data:
        raise ValueError("Empty file.")
    elif len(data) < 13:
        raise ValueError("File too small to hold a valid SFF header.")
    magic_number, ver0, ver1, ver2, ver3, index_offset, index_length, \
    number_of_reads, header_length, key_length, number_of_flows_per_read, \
    flowgram_format = struct.unpack(fmt, data)
    if magic_number in [_hsh, _srt, _mft]:
        #Probably user error, calling Bio.SeqIO.parse() twice!
        raise ValueError("Handle seems to be at SFF index block, not start")
    if magic_number != _sff: # 779314790
        raise ValueError("SFF file did not start '.sff', but %s" \
                         % repr(magic_number))
    if (ver0, ver1, ver2, ver3) != (0, 0, 0, 1):
        raise ValueError("Unsupported SFF version in header, %i.%i.%i.%i" \
                         % (ver0, ver1, ver2, ver3))
    if flowgram_format != 1:
        raise ValueError("Flowgram format code %i not supported" \
                         % flowgram_format)
    if (index_offset!=0) ^ (index_length!=0):
        raise ValueError("Index offset %i but index length %i" \
                         % (index_offset, index_length))
    flow_chars = _bytes_to_string(handle.read(number_of_flows_per_read))
    key_sequence = _bytes_to_string(handle.read(key_length))
    #According to the spec, the header_length field should be the total number
    #of bytes required by this set of header fields, and should be equal to
    #"31 + number_of_flows_per_read + key_length" rounded up to the next value
    #divisible by 8.
    assert header_length % 8 == 0
    padding = header_length - number_of_flows_per_read - key_length - 31
    assert 0 <= padding < 8, padding
    if handle.read(padding).count(_null) != padding:
        raise ValueError("Post header %i byte padding region contained data" \
                         % padding)
    return header_length, index_offset, index_length, \
           number_of_reads, number_of_flows_per_read, \
           flow_chars, key_sequence

#This is a generator function!
def _sff_do_slow_index(handle):
    """Generates an index by scanning though all the reads in an SFF file (PRIVATE).

    This is a slow but generic approach if we can't parse the provided index
    (if present).

    Will use the handle seek/tell functions.
    """
    handle.seek(0)
    header_length, index_offset, index_length, number_of_reads, \
    number_of_flows_per_read, flow_chars, key_sequence \
        = _sff_file_header(handle)
    #Now on to the reads...
    read_header_fmt = '>2HI4H'
    read_header_size = struct.calcsize(read_header_fmt)
    #NOTE - assuming flowgram_format==1, which means struct type H
    read_flow_fmt = ">%iH" % number_of_flows_per_read
    read_flow_size = struct.calcsize(read_flow_fmt)
    assert 1 == struct.calcsize(">B")
    assert 1 == struct.calcsize(">s")
    assert 1 == struct.calcsize(">c")
    assert read_header_size % 8 == 0 #Important for padding calc later!
    for read in range(number_of_reads):
        record_offset = handle.tell()
        if record_offset == index_offset:
            #Found index block within reads, ignore it:
            offset = index_offset + index_length
            if offset % 8:
                offset += 8 - (offset % 8)
            assert offset % 8 == 0
            handle.seek(offset)
            record_offset = offset
        #assert record_offset%8 == 0 #Worth checking, but slow
        #First the fixed header
        data = handle.read(read_header_size)
        read_header_length, name_length, seq_len, clip_qual_left, \
        clip_qual_right, clip_adapter_left, clip_adapter_right \
            = struct.unpack(read_header_fmt, data)
        if read_header_length < 10 or read_header_length % 8 != 0:
            raise ValueError("Malformed read header, says length is %i:\n%s" \
                             % (read_header_length, repr(data)))
        #now the name and any padding (remainder of header)
        name = _bytes_to_string(handle.read(name_length))
        padding = read_header_length - read_header_size - name_length
        if handle.read(padding).count(_null) != padding:
            raise ValueError("Post name %i byte padding region contained data" \
                             % padding)
        assert record_offset + read_header_length == handle.tell()
        #now the flowgram values, flowgram index, bases and qualities
        size = read_flow_size + 3*seq_len
        handle.seek(size, 1)
        #now any padding...
        padding = size % 8
        if padding:
            padding = 8 - padding
            if handle.read(padding).count(_null) != padding:
                raise ValueError("Post quality %i byte padding region contained data" \
                                 % padding)
        #print read, name, record_offset
        yield name, record_offset
    if handle.tell() % 8 != 0:
        raise ValueError("After scanning reads, did not end on a multiple of 8")

def _sff_find_roche_index(handle):
    """Locate any existing Roche style XML meta data and read index (PRIVATE).

    Makes a number of hard coded assumptions based on reverse engineered SFF
    files from Roche 454 machines.

    Returns a tuple of read count, SFF "index" offset and size, XML offset
    and size, and the actual read index offset and size.

    Raises a ValueError for unsupported or non-Roche index blocks.
    """
    handle.seek(0)
    header_length, index_offset, index_length, number_of_reads, \
    number_of_flows_per_read, flow_chars, key_sequence \
        = _sff_file_header(handle)
    assert handle.tell() == header_length
    if not index_offset or not index_offset:
        raise ValueError("No index present in this SFF file")
    #Now jump to the header...
    handle.seek(index_offset)
    fmt = ">4s4B"
    fmt_size = struct.calcsize(fmt)
    data = handle.read(fmt_size)
    if not data:
        raise ValueError("Premature end of file? Expected index of size %i at offest %i, found nothing" \
                         % (index_length, index_offset))
    if len(data) < fmt_size:
        raise ValueError("Premature end of file? Expected index of size %i at offest %i, found %s" \
                         % (index_length, index_offset, repr(data)))
    magic_number, ver0, ver1, ver2, ver3 = struct.unpack(fmt, data)
    if magic_number == _mft: # 778921588
        #Roche 454 manifest index
        #This is typical from raw Roche 454 SFF files (2009), and includes
        #both an XML manifest and the sorted index.
        if (ver0, ver1, ver2, ver3) != (49, 46, 48, 48):
            #This is "1.00" as a string
            raise ValueError("Unsupported version in .mft index header, %i.%i.%i.%i" \
                             % (ver0, ver1, ver2, ver3))
        fmt2 = ">LL"
        fmt2_size = struct.calcsize(fmt2)
        xml_size, data_size = struct.unpack(fmt2, handle.read(fmt2_size))
        if index_length != fmt_size + fmt2_size + xml_size + data_size:
            raise ValueError("Problem understanding .mft index header, %i != %i + %i + %i + %i" \
                             % (index_length, fmt_size, fmt2_size, xml_size, data_size))
        return number_of_reads, header_length, \
               index_offset, index_length, \
               index_offset + fmt_size + fmt2_size, xml_size, \
               index_offset + fmt_size + fmt2_size + xml_size, data_size
    elif magic_number == _srt: #779317876
        #Roche 454 sorted index
        #I've had this from Roche tool sfffile when the read identifiers
        #had nonstandard lengths and there was no XML manifest.
        if (ver0, ver1, ver2, ver3) != (49, 46, 48, 48):
            #This is "1.00" as a string
            raise ValueError("Unsupported version in .srt index header, %i.%i.%i.%i" \
                             % (ver0, ver1, ver2, ver3))
        data = handle.read(4)
        if data != _null*4:
            raise ValueError("Did not find expected null four bytes in .srt index")
        return number_of_reads, header_length, \
               index_offset, index_length, \
               0, 0, \
               index_offset + fmt_size + 4, index_length - fmt_size - 4
    elif magic_number == _hsh:
        raise ValueError("Hash table style indexes (.hsh) in SFF files are "
                         "not (yet) supported")
    else:
        raise ValueError("Unknown magic number %s in SFF index header:\n%s" \
                         % (repr(magic_number), repr(data)))

def _sff_read_roche_index_xml(handle):
    """Reads any existing Roche style XML manifest data in the SFF "index" (PRIVATE, DEPRECATED).

    Will use the handle seek/tell functions. Returns a string.

    This has been replaced by ReadRocheXmlManifest. We would normally just
    delete an old private function without warning, but I believe some people
    are using this so we'll handle this with a deprecation warning.
    """
    import warnings
    warnings.warn("Private function _sff_read_roche_index_xml is deprecated. "
                  "Use new public function ReadRocheXmlManifest instead",
                  DeprecationWarning)
    return ReadRocheXmlManifest(handle)


def ReadRocheXmlManifest(handle):
    """Reads any Roche style XML manifest data in the SFF "index".

    The SFF file format allows for multiple different index blocks, and Roche
    took advantage of this to define their own index block wich also embeds
    an XML manifest string. This is not a publically documented extension to
    the SFF file format, this was reverse engineered.

    The handle should be to an SFF file opened in binary mode. This function
    will use the handle seek/tell functions and leave the handle in an
    arbitrary location.

    Any XML manifest found is returned as a Python string, which you can then
    parse as appropriate, or reuse when writing out SFF files with the
    SffWriter class.

    Returns a string, or raises a ValueError if an Roche manifest could not be
    found.
    """
    number_of_reads, header_length, index_offset, index_length, xml_offset, \
    xml_size, read_index_offset, read_index_size = _sff_find_roche_index(handle)
    if not xml_offset or not xml_size:
        raise ValueError("No XML manifest found")
    handle.seek(xml_offset)
    return _bytes_to_string(handle.read(xml_size))


#This is a generator function!
def _sff_read_roche_index(handle):
    """Reads any existing Roche style read index provided in the SFF file (PRIVATE).

    Will use the handle seek/tell functions.

    This works on ".srt1.00" and ".mft1.00" style Roche SFF index blocks.

    Roche SFF indices use base 255 not 256, meaning we see bytes in range the
    range 0 to 254 only. This appears to be so that byte 0xFF (character 255)
    can be used as a marker character to separate entries (required if the
    read name lengths vary).

    Note that since only four bytes are used for the read offset, this is
    limited to 255^4 bytes (nearly 4GB). If you try to use the Roche sfffile
    tool to combine SFF files beyound this limit, they issue a warning and
    omit the index (and manifest).
    """
    number_of_reads, header_length, index_offset, index_length, xml_offset, \
    xml_size, read_index_offset, read_index_size = _sff_find_roche_index(handle)
    #Now parse the read index...
    handle.seek(read_index_offset)
    fmt = ">5B"
    for read in range(number_of_reads):
        #TODO - Be more aware of when the index should end?
        data = handle.read(6)
        while True:
            more = handle.read(1)
            if not more:
                raise ValueError("Premature end of file!")
            data += more
            if more == _flag: break
        assert data[-1:] == _flag, data[-1:]
        name = _bytes_to_string(data[:-6])
        off4, off3, off2, off1, off0 = struct.unpack(fmt, data[-6:-1])
        offset = off0 + 255*off1 + 65025*off2 + 16581375*off3
        if off4:
            #Could in theory be used as a fifth piece of offset information,
            #i.e. offset =+ 4228250625L*off4, but testing the Roche tools this
            #is not the case. They simple don't support such large indexes.
            raise ValueError("Expected a null terminator to the read name.")
        yield name, offset
    if handle.tell() != read_index_offset + read_index_size:
        raise ValueError("Problem with index length? %i vs %i" \
                         % (handle.tell(), read_index_offset + read_index_size))


def _sff_read_seq_record(handle, number_of_flows_per_read, flow_chars,
                         key_sequence, alphabet, trim=False):
    """Parse the next read in the file, return data as a SeqRecord (PRIVATE)."""
    #Now on to the reads...
    #the read header format (fixed part):
    #read_header_length     H
    #name_length            H
    #seq_len                I
    #clip_qual_left         H
    #clip_qual_right        H
    #clip_adapter_left      H
    #clip_adapter_right     H
    #[rest of read header depends on the name length etc]
    read_header_fmt = '>2HI4H'
    read_header_size = struct.calcsize(read_header_fmt)
    read_flow_fmt = ">%iH" % number_of_flows_per_read
    read_flow_size = struct.calcsize(read_flow_fmt)

    read_header_length, name_length, seq_len, clip_qual_left, \
    clip_qual_right, clip_adapter_left, clip_adapter_right \
        = struct.unpack(read_header_fmt, handle.read(read_header_size))
    if clip_qual_left:
        clip_qual_left -= 1 #python counting
    if clip_adapter_left:
        clip_adapter_left -= 1 #python counting
    if read_header_length < 10 or read_header_length % 8 != 0:
        raise ValueError("Malformed read header, says length is %i" \
                         % read_header_length)
    #now the name and any padding (remainder of header)
    name = _bytes_to_string(handle.read(name_length))
    padding = read_header_length - read_header_size - name_length
    if handle.read(padding).count(_null) != padding:
        raise ValueError("Post name %i byte padding region contained data" \
                         % padding)
    #now the flowgram values, flowgram index, bases and qualities
    #NOTE - assuming flowgram_format==1, which means struct type H
    flow_values = handle.read(read_flow_size) #unpack later if needed
    temp_fmt = ">%iB" % seq_len # used for flow index and quals
    flow_index = handle.read(seq_len) #unpack later if needed
    seq = _bytes_to_string(handle.read(seq_len)) #TODO - Use bytes in Seq?
    quals = list(struct.unpack(temp_fmt, handle.read(seq_len)))
    #now any padding...
    padding = (read_flow_size + seq_len*3)%8
    if padding:
        padding = 8 - padding
        if handle.read(padding).count(_null) != padding:
            raise ValueError("Post quality %i byte padding region contained data" \
                             % padding)
    #Follow Roche and apply most aggressive of qual and adapter clipping.
    #Note Roche seems to ignore adapter clip fields when writing SFF,
    #and uses just the quality clipping values for any clipping.
    clip_left = max(clip_qual_left, clip_adapter_left)
    #Right clipping of zero means no clipping
    if clip_qual_right:
        if clip_adapter_right:
            clip_right = min(clip_qual_right, clip_adapter_right)
        else:
            #Typical case with Roche SFF files
            clip_right = clip_qual_right
    elif clip_adapter_right:
        clip_right = clip_adapter_right
    else:
        clip_right = seq_len
    #Now build a SeqRecord
    if trim:
        seq = seq[clip_left:clip_right].upper()
        quals = quals[clip_left:clip_right]
        #Don't record the clipping values, flow etc, they make no sense now:
        annotations = {}
    else:
        #This use of mixed case mimics the Roche SFF tool's FASTA output
        seq = seq[:clip_left].lower() + \
              seq[clip_left:clip_right].upper() + \
              seq[clip_right:].lower()
        annotations = {"flow_values":struct.unpack(read_flow_fmt, flow_values),
                       "flow_index":struct.unpack(temp_fmt, flow_index),
                       "flow_chars":flow_chars,
                       "flow_key":key_sequence,
                       "clip_qual_left":clip_qual_left,
                       "clip_qual_right":clip_qual_right,
                       "clip_adapter_left":clip_adapter_left,
                       "clip_adapter_right":clip_adapter_right}
    record = SeqRecord(Seq(seq, alphabet),
                       id=name,
                       name=name,
                       description="",
                       annotations=annotations)
    #Dirty trick to speed up this line:
    #record.letter_annotations["phred_quality"] = quals
    dict.__setitem__(record._per_letter_annotations,
                     "phred_quality", quals)
    #Return the record and then continue...
    return record


def _sff_read_raw_record(handle, number_of_flows_per_read):
    """Extract the next read in the file as a raw (bytes) string (PRIVATE)."""
    read_header_fmt = '>2HI'
    read_header_size = struct.calcsize(read_header_fmt)
    read_flow_fmt = ">%iH" % number_of_flows_per_read
    read_flow_size = struct.calcsize(read_flow_fmt)

    raw = handle.read(read_header_size)
    read_header_length, name_length, seq_len \
                        = struct.unpack(read_header_fmt, raw)
    if read_header_length < 10 or read_header_length % 8 != 0:
        raise ValueError("Malformed read header, says length is %i" \
                         % read_header_length)
    #now the four clip values (4H = 8 bytes), and read name
    raw += handle.read(8 + name_length)
    #and any padding (remainder of header)
    padding = read_header_length - read_header_size - 8 - name_length
    pad = handle.read(padding)
    if pad.count(_null) != padding:
        raise ValueError("Post name %i byte padding region contained data" \
                         % padding)
    raw += pad
    #now the flowgram values, flowgram index, bases and qualities
    raw += handle.read(read_flow_size + seq_len*3)
    padding = (read_flow_size + seq_len*3)%8
    #now any padding...
    if padding:
        padding = 8 - padding
        pad = handle.read(padding)
        if pad.count(_null) != padding:
            raise ValueError("Post quality %i byte padding region contained data" \
                             % padding)
        raw += pad
    #Return the raw bytes
    return raw

class _AddTellHandle(object):
    """Wrapper for handles which do not support the tell method (PRIVATE).

    Intended for use with things like network handles where tell (and reverse
    seek) are not supported. The SFF file needs to track the current offset in
    order to deal with the index block.
    """
    def __init__(self, handle):
        self._handle = handle
        self._offset = 0

    def read(self, length):
        data = self._handle.read(length)
        self._offset += len(data)
        return data

    def tell(self):
        return self._offset

    def seek(self, offset):
        if offset < self._offset:
            raise RunTimeError("Can't seek backwards")
        self._handle.read(offset - self._offset)

    def close(self):
        return self._handle.close()


#This is a generator function!
def SffIterator(handle, alphabet=Alphabet.generic_dna, trim=False):
    """Iterate over Standard Flowgram Format (SFF) reads (as SeqRecord objects).

    handle - input file, an SFF file, e.g. from Roche 454 sequencing.
             This must NOT be opened in universal read lines mode!
    alphabet - optional alphabet, defaults to generic DNA.
    trim - should the sequences be trimmed?

    The resulting SeqRecord objects should match those from a paired FASTA
    and QUAL file converted from the SFF file using the Roche 454 tool
    ssfinfo. i.e. The sequence will be mixed case, with the trim regions
    shown in lower case.

    This function is used internally via the Bio.SeqIO functions:

    >>> from Bio import SeqIO
    >>> handle = open("Roche/E3MFGYR02_random_10_reads.sff", "rb")
    >>> for record in SeqIO.parse(handle, "sff"):
    ...     print record.id, len(record)
    E3MFGYR02JWQ7T 265
    E3MFGYR02JA6IL 271
    E3MFGYR02JHD4H 310
    E3MFGYR02GFKUC 299
    E3MFGYR02FTGED 281
    E3MFGYR02FR9G7 261
    E3MFGYR02GAZMS 278
    E3MFGYR02HHZ8O 221
    E3MFGYR02GPGB1 269
    E3MFGYR02F7Z7G 219
    >>> handle.close()

    You can also call it directly:

    >>> handle = open("Roche/E3MFGYR02_random_10_reads.sff", "rb")
    >>> for record in SffIterator(handle):
    ...     print record.id, len(record)
    E3MFGYR02JWQ7T 265
    E3MFGYR02JA6IL 271
    E3MFGYR02JHD4H 310
    E3MFGYR02GFKUC 299
    E3MFGYR02FTGED 281
    E3MFGYR02FR9G7 261
    E3MFGYR02GAZMS 278
    E3MFGYR02HHZ8O 221
    E3MFGYR02GPGB1 269
    E3MFGYR02F7Z7G 219
    >>> handle.close()

    Or, with the trim option:

    >>> handle = open("Roche/E3MFGYR02_random_10_reads.sff", "rb")
    >>> for record in SffIterator(handle, trim=True):
    ...     print record.id, len(record)
    E3MFGYR02JWQ7T 260
    E3MFGYR02JA6IL 265
    E3MFGYR02JHD4H 292
    E3MFGYR02GFKUC 295
    E3MFGYR02FTGED 277
    E3MFGYR02FR9G7 256
    E3MFGYR02GAZMS 271
    E3MFGYR02HHZ8O 150
    E3MFGYR02GPGB1 221
    E3MFGYR02F7Z7G 130
    >>> handle.close()

    """
    if isinstance(Alphabet._get_base_alphabet(alphabet),
                  Alphabet.ProteinAlphabet):
        raise ValueError("Invalid alphabet, SFF files do not hold proteins.")
    if isinstance(Alphabet._get_base_alphabet(alphabet),
                  Alphabet.RNAAlphabet):
        raise ValueError("Invalid alphabet, SFF files do not hold RNA.")
    try:
        assert 0 == handle.tell()
    except AttributeError:
        #Probably a network handle or something like that
        handle = _AddTellHandle(handle)
    header_length, index_offset, index_length, number_of_reads, \
    number_of_flows_per_read, flow_chars, key_sequence \
        = _sff_file_header(handle)
    #Now on to the reads...
    #the read header format (fixed part):
    #read_header_length     H
    #name_length            H
    #seq_len                I
    #clip_qual_left         H
    #clip_qual_right        H
    #clip_adapter_left      H
    #clip_adapter_right     H
    #[rest of read header depends on the name length etc]
    read_header_fmt = '>2HI4H'
    read_header_size = struct.calcsize(read_header_fmt)
    read_flow_fmt = ">%iH" % number_of_flows_per_read
    read_flow_size = struct.calcsize(read_flow_fmt)
    assert 1 == struct.calcsize(">B")
    assert 1 == struct.calcsize(">s")
    assert 1 == struct.calcsize(">c")
    assert read_header_size % 8 == 0 #Important for padding calc later!
    #The spec allows for the index block to be before or even in the middle
    #of the reads. We can check that if we keep track of our position
    #in the file...
    for read in range(number_of_reads):
        if index_offset and handle.tell() == index_offset:
            offset = index_offset + index_length
            if offset % 8:
                offset += 8 - (offset % 8)
            assert offset % 8 == 0
            handle.seek(offset)
            #Now that we've done this, we don't need to do it again. Clear
            #the index_offset so we can skip extra handle.tell() calls:
            index_offset = 0
        yield _sff_read_seq_record(handle,
                                   number_of_flows_per_read,
                                   flow_chars,
                                   key_sequence,
                                   alphabet,
                                   trim)
    #The following is not essential, but avoids confusing error messages
    #for the user if they try and re-parse the same handle.
    if index_offset and handle.tell() == index_offset:
        offset = index_offset + index_length
        if offset % 8:
            offset += 8 - (offset % 8)
        assert offset % 8 == 0
        handle.seek(offset)
    #Should now be at the end of the file...
    if handle.read(1):
        raise ValueError("Additional data at end of SFF file")


#This is a generator function!
def _SffTrimIterator(handle, alphabet=Alphabet.generic_dna):
    """Iterate over SFF reads (as SeqRecord objects) with trimming (PRIVATE)."""
    return SffIterator(handle, alphabet, trim=True)


class SffWriter(SequenceWriter):
    """SFF file writer."""

    def __init__(self, handle, index=True, xml=None):
        """Creates the writer object.

        handle - Output handle, ideally in binary write mode.
        index - Boolean argument, should we try and write an index?
        xml - Optional string argument, xml manifest to be recorded in the index
              block (see function ReadRocheXmlManifest for reading this data).
        """
        if hasattr(handle,"mode") and "U" in handle.mode.upper():
            raise ValueError("SFF files must NOT be opened in universal new "
                             "lines mode. Binary mode is required")
        elif hasattr(handle,"mode") and "B" not in handle.mode.upper():
            raise ValueError("SFF files must be opened in binary mode")
        self.handle = handle
        self._xml = xml
        if index:
            self._index = []
        else:
            self._index = None

    def write_file(self, records):
        """Use this to write an entire file containing the given records."""
        try:
            self._number_of_reads = len(records)
        except TypeError:
            self._number_of_reads = 0 #dummy value
            if not hasattr(self.handle, "seek") \
            or not hasattr(self.handle, "tell"):
                raise ValueError("A handle with a seek/tell methods is "
                                 "required in order to record the total "
                                 "record count in the file header (once it "
                                 "is known at the end).")
        if self._index is not None and \
        not (hasattr(self.handle, "seek") and hasattr(self.handle, "tell")):
            import warnings
            warnings.warn("A handle with a seek/tell methods is required in "
                          "order to record an SFF index.")
            self._index = None
        self._index_start = 0
        self._index_length = 0
        if not hasattr(records, "next"):
            records = iter(records)
        #Get the first record in order to find the flow information
        #we will need for the header.
        try:
            record = records.next()
        except StopIteration:
            record = None
        if record is None:
            #No records -> empty SFF file (or an error)?
            #We can't write a header without the flow information.
            #return 0
            raise ValueError("Must have at least one sequence")
        try:
            self._key_sequence = _as_bytes(record.annotations["flow_key"])
            self._flow_chars = _as_bytes(record.annotations["flow_chars"])
            self._number_of_flows_per_read = len(self._flow_chars)
        except KeyError:
            raise ValueError("Missing SFF flow information")
        self.write_header()
        self.write_record(record)
        count = 1
        for record in records:
            self.write_record(record)
            count += 1
        if self._number_of_reads == 0:
            #Must go back and record the record count...
            offset = self.handle.tell()
            self.handle.seek(0)
            self._number_of_reads = count
            self.write_header()
            self.handle.seek(offset) #not essential?
        else:
            assert count == self._number_of_reads
        if self._index is not None:
            self._write_index()
        return count

    def _write_index(self):
        assert len(self._index)==self._number_of_reads
        handle = self.handle
        self._index.sort()
        self._index_start = handle.tell() #need for header
        #XML...
        if self._xml is not None:
            xml = _as_bytes(self._xml)
        else:
            from Bio import __version__
            xml = "<!-- This file was output with Biopython %s -->\n" % __version__
            xml += "<!-- This XML and index block attempts to mimic Roche SFF files -->\n"
            xml += "<!-- This file may be a combination of multiple SFF files etc -->\n"
            xml = _as_bytes(xml)
        xml_len = len(xml)
        #Write to the file...
        fmt = ">I4BLL"
        fmt_size = struct.calcsize(fmt)
        handle.write(_null*fmt_size + xml) #will come back later to fill this
        fmt2 = ">6B"
        assert 6 == struct.calcsize(fmt2)
        self._index.sort()
        index_len = 0 #don't know yet!
        for name, offset in self._index:
            #Roche files record the offsets using base 255 not 256.
            #See comments for parsing the index block. There may be a faster
            #way to code this, but we can't easily use shifts due to odd base
            off3 = offset
            off0 = off3 % 255
            off3 -= off0
            off1 = off3 % 65025
            off3 -= off1
            off2 = off3 % 16581375
            off3 -= off2
            assert offset == off0 + off1 + off2 + off3, \
                   "%i -> %i %i %i %i" % (offset, off0, off1, off2, off3)
            off3, off2, off1, off0 = off3//16581375, off2//65025, \
                                     off1//255, off0
            assert off0 < 255 and off1 < 255 and off2 < 255 and off3 < 255, \
                   "%i -> %i %i %i %i" % (offset, off0, off1, off2, off3)
            handle.write(name + struct.pack(fmt2, 0, \
                                            off3, off2, off1, off0, 255))
            index_len += len(name) + 6
        #Note any padding in not included:
        self._index_length = fmt_size + xml_len + index_len #need for header
        #Pad out to an 8 byte boundary (although I have noticed some
        #real Roche SFF files neglect to do this depsite their manual
        #suggesting this padding should be there):
        if self._index_length % 8:
            padding = 8 - (self._index_length%8)
            handle.write(_null*padding)
        else:
            padding = 0
        offset = handle.tell()
        assert offset == self._index_start + self._index_length + padding, \
               "%i vs %i + %i + %i"  % (offset, self._index_start, \
                                        self._index_length, padding)
        #Must now go back and update the index header with index size...
        handle.seek(self._index_start)
        handle.write(struct.pack(fmt, 778921588, #magic number
                                 49,46,48,48, #Roche index version, "1.00"
                                 xml_len, index_len) + xml)
        #Must now go back and update the header...
        handle.seek(0)
        self.write_header()
        handle.seek(offset) #not essential?

    def write_header(self):
        #Do header...
        key_length = len(self._key_sequence)
        #file header (part one)
        #use big endiean encdoing   >
        #magic_number               I
        #version                    4B
        #index_offset               Q
        #index_length               I
        #number_of_reads            I
        #header_length              H
        #key_length                 H
        #number_of_flows_per_read   H
        #flowgram_format_code       B
        #[rest of file header depends on the number of flows and how many keys]
        fmt = '>I4BQIIHHHB%is%is' % (self._number_of_flows_per_read, key_length)
        #According to the spec, the header_length field should be the total
        #number of bytes required by this set of header fields, and should be
        #equal to "31 + number_of_flows_per_read + key_length" rounded up to
        #the next value divisible by 8.
        if struct.calcsize(fmt) % 8 == 0:
            padding = 0
        else:
            padding = 8 - (struct.calcsize(fmt) % 8)
        header_length = struct.calcsize(fmt) + padding
        assert header_length % 8 == 0
        header = struct.pack(fmt, 779314790, #magic number 0x2E736666
                             0, 0, 0, 1, #version
                             self._index_start, self._index_length,
                             self._number_of_reads,
                             header_length, key_length,
                             self._number_of_flows_per_read,
                             1, #the only flowgram format code we support
                             self._flow_chars, self._key_sequence)
        self.handle.write(header + _null*padding)

    def write_record(self, record):
        """Write a single additional record to the output file.

        This assumes the header has been done.
        """
        #Basics
        name = _as_bytes(record.id)
        name_len = len(name)
        seq = _as_bytes(str(record.seq).upper())
        seq_len = len(seq)
        #Qualities
        try:
            quals = record.letter_annotations["phred_quality"]
        except KeyError:
            raise ValueError("Missing PHRED qualities information")
        #Flow
        try:
            flow_values = record.annotations["flow_values"]
            flow_index = record.annotations["flow_index"]
            if self._key_sequence != _as_bytes(record.annotations["flow_key"]) \
            or self._flow_chars != _as_bytes(record.annotations["flow_chars"]):
                raise ValueError("Records have inconsistent SFF flow data")
        except KeyError:
            raise ValueError("Missing SFF flow information")
        except AttributeError:
            raise ValueError("Header not written yet?")
        #Clipping
        try:
            clip_qual_left = record.annotations["clip_qual_left"]
            if clip_qual_left:
                clip_qual_left += 1
            clip_qual_right = record.annotations["clip_qual_right"]
            clip_adapter_left = record.annotations["clip_adapter_left"]
            if clip_adapter_left:
                clip_adapter_left += 1
            clip_adapter_right = record.annotations["clip_adapter_right"]
        except KeyError:
            raise ValueError("Missing SFF clipping information")

        #Capture information for index
        if self._index is not None:
            offset = self.handle.tell()
            #Check the position of the final record (before sort by name)
            #See comments earlier about how base 255 seems to be used.
            #This means the limit is 255**4 + 255**3 +255**2 + 255**1
            if offset > 4244897280:
                import warnings
                warnings.warn("Read %s has file offset %i, which is too large "
                              "to store in the Roche SFF index structure. No "
                              "index block will be recorded." % (name, offset))
                #No point recoring the offsets now
                self._index = None
            else:
                self._index.append((name, self.handle.tell()))

        #the read header format (fixed part):
        #read_header_length     H
        #name_length            H
        #seq_len                I
        #clip_qual_left         H
        #clip_qual_right        H
        #clip_adapter_left      H
        #clip_adapter_right     H
        #[rest of read header depends on the name length etc]
        #name
        #flow values
        #flow index
        #sequence
        #padding
        read_header_fmt = '>2HI4H%is' % name_len
        if struct.calcsize(read_header_fmt) % 8 == 0:
            padding = 0
        else:
            padding = 8 - (struct.calcsize(read_header_fmt) % 8)
        read_header_length = struct.calcsize(read_header_fmt) + padding
        assert read_header_length % 8 == 0
        data = struct.pack(read_header_fmt,
                           read_header_length,
                           name_len, seq_len,
                           clip_qual_left, clip_qual_right,
                           clip_adapter_left, clip_adapter_right,
                           name) + _null*padding
        assert len(data) == read_header_length
        #now the flowgram values, flowgram index, bases and qualities
        #NOTE - assuming flowgram_format==1, which means struct type H
        read_flow_fmt = ">%iH" % self._number_of_flows_per_read
        read_flow_size = struct.calcsize(read_flow_fmt)
        temp_fmt = ">%iB" % seq_len # used for flow index and quals
        data += struct.pack(read_flow_fmt, *flow_values) \
                + struct.pack(temp_fmt, *flow_index) \
                + seq \
                + struct.pack(temp_fmt, *quals)
        #now any final padding...
        padding = (read_flow_size + seq_len*3)%8
        if padding:
            padding = 8 - padding
        self.handle.write(data + _null*padding)


if __name__ == "__main__":
    print "Running quick self test"
    filename = "../../Tests/Roche/E3MFGYR02_random_10_reads.sff"
    metadata = ReadRocheXmlManifest(open(filename, "rb"))
    index1 = sorted(_sff_read_roche_index(open(filename, "rb")))
    index2 = sorted(_sff_do_slow_index(open(filename, "rb")))
    assert index1 == index2
    assert len(index1) == len(list(SffIterator(open(filename, "rb"))))
    from StringIO import StringIO
    try:
        #This is in Python 2.6+, and is essential on Python 3
        from io import BytesIO
    except ImportError:
        BytesIO = StringIO
    assert len(index1) == len(list(SffIterator(BytesIO(open(filename,"rb").read()))))

    if sys.platform != "win32":
        assert len(index1) == len(list(SffIterator(open(filename, "r"))))
        index2 = sorted(_sff_read_roche_index(open(filename)))
        assert index1 == index2
        index2 = sorted(_sff_do_slow_index(open(filename)))
        assert index1 == index2
        assert len(index1) == len(list(SffIterator(open(filename))))
        assert len(index1) == len(list(SffIterator(BytesIO(open(filename,"r").read()))))
        assert len(index1) == len(list(SffIterator(BytesIO(open(filename).read()))))

    sff = list(SffIterator(open(filename, "rb")))

    sff2 = list(SffIterator(open("../../Tests/Roche/E3MFGYR02_alt_index_at_end.sff", "rb")))
    assert len(sff) == len(sff2)
    for old, new in zip(sff, sff2):
        assert old.id == new.id
        assert str(old.seq) == str(new.seq)

    sff2 = list(SffIterator(open("../../Tests/Roche/E3MFGYR02_alt_index_at_start.sff", "rb")))
    assert len(sff) == len(sff2)
    for old, new in zip(sff, sff2):
        assert old.id == new.id
        assert str(old.seq) == str(new.seq)

    sff2 = list(SffIterator(open("../../Tests/Roche/E3MFGYR02_alt_index_in_middle.sff", "rb")))
    assert len(sff) == len(sff2)
    for old, new in zip(sff, sff2):
        assert old.id == new.id
        assert str(old.seq) == str(new.seq)

    sff2 = list(SffIterator(open("../../Tests/Roche/E3MFGYR02_index_at_start.sff", "rb")))
    assert len(sff) == len(sff2)
    for old, new in zip(sff, sff2):
        assert old.id == new.id
        assert str(old.seq) == str(new.seq)

    sff2 = list(SffIterator(open("../../Tests/Roche/E3MFGYR02_index_in_middle.sff", "rb")))
    assert len(sff) == len(sff2)
    for old, new in zip(sff, sff2):
        assert old.id == new.id
        assert str(old.seq) == str(new.seq)

    sff_trim = list(SffIterator(open(filename, "rb"), trim=True))

    print ReadRocheXmlManifest(open(filename, "rb"))

    from Bio import SeqIO
    filename = "../../Tests/Roche/E3MFGYR02_random_10_reads_no_trim.fasta"
    fasta_no_trim = list(SeqIO.parse(open(filename,"rU"), "fasta"))
    filename = "../../Tests/Roche/E3MFGYR02_random_10_reads_no_trim.qual"
    qual_no_trim = list(SeqIO.parse(open(filename,"rU"), "qual"))

    filename = "../../Tests/Roche/E3MFGYR02_random_10_reads.fasta"
    fasta_trim = list(SeqIO.parse(open(filename,"rU"), "fasta"))
    filename = "../../Tests/Roche/E3MFGYR02_random_10_reads.qual"
    qual_trim = list(SeqIO.parse(open(filename,"rU"), "qual"))

    for s, sT, f, q, fT, qT in zip(sff, sff_trim, fasta_no_trim,
                                   qual_no_trim, fasta_trim, qual_trim):
        #print
        print s.id
        #print s.seq
        #print s.letter_annotations["phred_quality"]

        assert s.id == f.id == q.id
        assert str(s.seq) == str(f.seq)
        assert s.letter_annotations["phred_quality"] == q.letter_annotations["phred_quality"]

        assert s.id == sT.id == fT.id == qT.id
        assert str(sT.seq) == str(fT.seq)
        assert sT.letter_annotations["phred_quality"] == qT.letter_annotations["phred_quality"]


    print "Writing with a list of SeqRecords..."
    handle = StringIO()
    w = SffWriter(handle, xml=metadata)
    w.write_file(sff) #list
    data = handle.getvalue()
    print "And again with an iterator..."
    handle = StringIO()
    w = SffWriter(handle, xml=metadata)
    w.write_file(iter(sff))
    assert data == handle.getvalue()
    #Check 100% identical to the original:
    filename = "../../Tests/Roche/E3MFGYR02_random_10_reads.sff"
    original = open(filename,"rb").read()
    assert len(data) == len(original)
    assert data == original
    del data
    handle.close()

    print "-"*50
    filename = "../../Tests/Roche/greek.sff"
    for record in SffIterator(open(filename,"rb")):
        print record.id
    index1 = sorted(_sff_read_roche_index(open(filename, "rb")))
    index2 = sorted(_sff_do_slow_index(open(filename, "rb")))
    assert index1 == index2
    try:
        print ReadRocheXmlManifest(open(filename, "rb"))
        assert False, "Should fail!"
    except ValueError:
        pass


    handle = open(filename, "rb")
    for record in SffIterator(handle):
        pass
    try:
        for record in SffIterator(handle):
            print record.id
        assert False, "Should have failed"
    except ValueError, err:
        print "Checking what happens on re-reading a handle:"
        print err


    """
    #Ugly code to make test files...
    index = ".diy1.00This is a fake index block (DIY = Do It Yourself), which is allowed under the SFF standard.\0"
    padding = len(index)%8
    if padding:
        padding = 8 - padding
    index += chr(0)*padding
    assert len(index)%8 == 0

    #Ugly bit of code to make a fake index at start
    records = list(SffIterator(open("../../Tests/Roche/E3MFGYR02_random_10_reads.sff", "rb")))
    out_handle = open("../../Tests/Roche/E3MFGYR02_alt_index_at_start.sff", "w")
    index = ".diy1.00This is a fake index block (DIY = Do It Yourself), which is allowed under the SFF standard.\0"
    padding = len(index)%8
    if padding:
        padding = 8 - padding
    index += chr(0)*padding
    w = SffWriter(out_handle, index=False, xml=None)
    #Fake the header...
    w._number_of_reads = len(records)
    w._index_start = 0
    w._index_length = 0
    w._key_sequence = records[0].annotations["flow_key"]
    w._flow_chars = records[0].annotations["flow_chars"]
    w._number_of_flows_per_read = len(w._flow_chars)
    w.write_header()
    w._index_start = out_handle.tell()
    w._index_length = len(index)
    out_handle.seek(0)
    w.write_header() #this time with index info
    w.handle.write(index)
    for record in records:
        w.write_record(record)
    out_handle.close()
    records2 = list(SffIterator(open("../../Tests/Roche/E3MFGYR02_alt_index_at_start.sff", "rb")))
    for old, new in zip(records, records2):
        assert str(old.seq)==str(new.seq)
    i = list(_sff_do_slow_index(open("../../Tests/Roche/E3MFGYR02_alt_index_at_start.sff", "rb")))

    #Ugly bit of code to make a fake index in middle
    records = list(SffIterator(open("../../Tests/Roche/E3MFGYR02_random_10_reads.sff", "rb")))
    out_handle = open("../../Tests/Roche/E3MFGYR02_alt_index_in_middle.sff", "w")
    index = ".diy1.00This is a fake index block (DIY = Do It Yourself), which is allowed under the SFF standard.\0"
    padding = len(index)%8
    if padding:
        padding = 8 - padding
    index += chr(0)*padding
    w = SffWriter(out_handle, index=False, xml=None)
    #Fake the header...
    w._number_of_reads = len(records)
    w._index_start = 0
    w._index_length = 0
    w._key_sequence = records[0].annotations["flow_key"]
    w._flow_chars = records[0].annotations["flow_chars"]
    w._number_of_flows_per_read = len(w._flow_chars)
    w.write_header()
    for record in records[:5]:
        w.write_record(record)
    w._index_start = out_handle.tell()
    w._index_length = len(index)
    w.handle.write(index)
    for record in records[5:]:
        w.write_record(record)
    out_handle.seek(0)
    w.write_header() #this time with index info
    out_handle.close()
    records2 = list(SffIterator(open("../../Tests/Roche/E3MFGYR02_alt_index_in_middle.sff", "rb")))
    for old, new in zip(records, records2):
        assert str(old.seq)==str(new.seq)
    j = list(_sff_do_slow_index(open("../../Tests/Roche/E3MFGYR02_alt_index_in_middle.sff", "rb")))

    #Ugly bit of code to make a fake index at end
    records = list(SffIterator(open("../../Tests/Roche/E3MFGYR02_random_10_reads.sff", "rb")))
    out_handle = open("../../Tests/Roche/E3MFGYR02_alt_index_at_end.sff", "w")
    w = SffWriter(out_handle, index=False, xml=None)
    #Fake the header...
    w._number_of_reads = len(records)
    w._index_start = 0
    w._index_length = 0
    w._key_sequence = records[0].annotations["flow_key"]
    w._flow_chars = records[0].annotations["flow_chars"]
    w._number_of_flows_per_read = len(w._flow_chars)
    w.write_header()
    for record in records:
        w.write_record(record)
    w._index_start = out_handle.tell()
    w._index_length = len(index)
    out_handle.write(index)
    out_handle.seek(0)
    w.write_header() #this time with index info
    out_handle.close()
    records2 = list(SffIterator(open("../../Tests/Roche/E3MFGYR02_alt_index_at_end.sff", "rb")))
    for old, new in zip(records, records2):
        assert str(old.seq)==str(new.seq)
    try:
        print ReadRocheXmlManifest(open("../../Tests/Roche/E3MFGYR02_alt_index_at_end.sff", "rb"))
        assert False, "Should fail!"
    except ValueError:
        pass
    k = list(_sff_do_slow_index(open("../../Tests/Roche/E3MFGYR02_alt_index_at_end.sff", "rb")))
    """

    print "Done"
