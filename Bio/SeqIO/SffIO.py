# Copyright 2009 by Peter Cock.  All rights reserved.
# Based on code contributed and copyright 2009 by Jose Blanca (COMAV-UPV).
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Bio.SeqIO support for the binary Standard Flowgram Format (SFF) file format.

SFF was designed by 454 Life Sciences (Roche), the Whitehead Institute for
Biomedical Research and the Wellcome Trust Sanger Institute. You are expected
to use this module via the Bio.SeqIO functions under the format name "sff".

For a description of the file format, please see:
http://www.ncbi.nlm.nih.gov/Traces/trace.cgi?cmd=show&f=formats&m=doc&s=formats

"""
from Interfaces import SequenceWriter
from Bio import Alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import struct
import sys

def _sff_file_header(handle) :
    """Read in an SFF file header (PRIVATE).

    Assumes the handle is at the start of the file, will read forwards
    though the header and leave the handle pointing at the first record.
    Returns a tuple of values from the header.
    """
    if hasattr(handle,"mode") and "U" in handle.mode.upper() :
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
    fmt = '>I4BQIIHHHB'
    assert 31 == struct.calcsize(fmt)
    magic_number, ver0, ver1, ver2, ver3, index_offset, index_length, \
    number_of_reads, header_length, key_length, number_of_flows_per_read, \
    flowgram_format = struct.unpack(fmt, handle.read(31))
    if magic_number != 779314790 :
        raise ValueError("Wrong SFF magic number in header")
    if (ver0, ver1, ver2, ver3) != (0,0,0,1) :
        raise ValueError("Unsupported SFF version in header, %i.%i.%i.%i" \
                         % (ver0, ver1, ver2, ver3))
    if flowgram_format != 1 :
        raise ValueError("Flowgram format code %i not supported" \
                         % flowgram_format)
    if (index_offset!=0) ^ (index_length!=0) :
        raise ValueError("Index offset %i but index length %i" \
                         % (index_offset, index_length))
    flow_chars = handle.read(number_of_flows_per_read)
    key_sequence = handle.read(key_length)
    #According to the spec, the header_length field should be the total number
    #of bytes required by this set of header fields, and should be equal to
    #"31 + number_of_flows_per_read + key_length" rounded up to the next value
    #divisible by 8.
    assert header_length % 8 == 0
    padding = header_length - number_of_flows_per_read - key_length - 31
    assert 0 <= padding < 8, padding
    if chr(0)*padding != handle.read(padding) :
        raise ValueError("Post header %i byte padding region contained data" \
                         % padding)
    return header_length, index_offset, index_length, \
           number_of_reads, number_of_flows_per_read, \
           flow_chars, key_sequence

#This is a generator function!
def _sff_do_slow_index(handle) :
    """Generates an index by scanning though all the reads in an SFF file (PRIVATE).

    This is a slow but generic approach if we can't parse the provided index (if
    present).

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
    for read in range(number_of_reads) :
        record_offset = handle.tell()
        #assert record_offset%8 == 0 #Worth checking, but slow
        #First the fixed header
        read_header_length, name_length, seq_len, clip_qual_left, \
        clip_qual_right, clip_adapter_left, clip_adapter_right \
            = struct.unpack(read_header_fmt, handle.read(read_header_size))
        if read_header_length < 10 or read_header_length%8 != 0 :
            raise ValueError("Malformed read header, says length is %i" \
                             % read_header_length)
        #now the name and any padding (remainder of header)
        name = handle.read(name_length)
        padding = read_header_length - read_header_size - name_length
        if chr(0)*padding != handle.read(padding) :
            raise ValueError("Post name %i byte padding region contained data" \
                             % padding)
        assert record_offset + read_header_length == handle.tell()
        #now the flowgram values, flowgram index, bases and qualities
        size = read_flow_size + 3*seq_len
        handle.seek(size,1)
        #now any padding...
        padding = size%8
        if padding :
            padding = 8 - padding
            if chr(0)*padding != handle.read(padding) :
                raise ValueError("Post quality %i byte padding region contained data" \
                                 % padding)
        #print read, name, record_offset
        yield name, record_offset
    if handle.tell() % 8 != 0 :
        raise ValueError("After scanning reads, did not end on a multiple of 8")

def _sff_find_roche_index(handle) :
    """Locate any existing Roche style XML meta data and read index (PRIVATE).

    Makes a number of hard coded assumptions based on reverse engineered SFF
    files from Roche 454 machines.

    Returns a tuple of read count, SFF "index" offset and size, XML offset and
    size, and the actual read index offset and size."""    
    handle.seek(0)
    header_length, index_offset, index_length, number_of_reads, \
    number_of_flows_per_read, flow_chars, key_sequence \
        = _sff_file_header(handle)
    assert handle.tell() == header_length
    if not index_offset or not index_offset :
        raise ValueError("No index present in this SFF file")
    #Now jump to the header...
    handle.seek(index_offset)
    fmt = ">4s4B"
    fmt_size = struct.calcsize(fmt)
    data = handle.read(fmt_size)
    magic_number, ver0, ver1, ver2, ver3 = struct.unpack(fmt, data)
    if magic_number == ".mft" : # 778921588
        #Roche 454 manifest index
        #This is typical from raw Roche 454 SFF files (2009), and includes
        #both an XML manifest and the sorted index.
        if (ver0, ver1, ver2, ver3) != (49,46,48,48) :
            #This is "1.00" as a string
            raise ValueError("Unsupported version in .mft index header, %i.%i.%i.%i" \
                             % (ver0, ver1, ver2, ver3))
        fmt2 = ">LL"
        fmt2_size = struct.calcsize(fmt2)
        xml_size, data_size = struct.unpack(fmt2, handle.read(fmt2_size))
        if index_length != fmt_size + fmt2_size + xml_size + data_size :
            raise ValueError("Problem understanding .mft index header, %i != %i + %i + %i + %i" \
                             % (index_length, fmt_size, fmt2_size, xml_size, data_size))
        return number_of_reads, header_length, \
               index_offset, index_length, \
               index_offset + fmt_size + fmt2_size, xml_size, \
               index_offset + fmt_size + fmt2_size + xml_size, data_size
    elif magic_number == ".srt" : #779317876
        #Roche 454 sorted index
        #I've had this from Roche tool sfffile when the read identifiers
        #had nonstandard lengths and there was no XML manifest.
        if (ver0, ver1, ver2, ver3) != (49,46,48,48) :
            #This is "1.00" as a string
            raise ValueError("Unsupported version in .srt index header, %i.%i.%i.%i" \
                             % (ver0, ver1, ver2, ver3))
        data = handle.read(4)
        if data != chr(0)*4 :
            raise ValueError("Did not find expected null four bytes in .srt index")
        return number_of_reads, header_length, \
               index_offset, index_length, \
               0, 0, \
               index_offset + fmt_size + 4, index_length - fmt_size - 4
    elif magic_number == ".hsh" :
        raise ValueError("Hash table style indexes (.hsh) in SFF files are "
                         "not (yet) supported")
    else :
        raise ValueError("Unknown magic number %s in SFF index header:\n%s" \
                         % (repr(magic_number), repr(data)))

def _sff_read_roche_index_xml(handle) :
    """Reads any existing Roche style XML manifest data in the SFF "index" (PRIVATE).

    Will use the handle seek/tell functions. Returns a string.
    """
    number_of_reads, header_length, index_offset, index_length, xml_offset, \
    xml_size, read_index_offset, read_index_size = _sff_find_roche_index(handle)
    if not xml_offset or not xml_size :
        raise ValueError("No XML manifest found")
    handle.seek(xml_offset)
    return handle.read(xml_size)


#This is a generator function!
def _sff_read_roche_index(handle) :
    """Reads any existing Roche style read index provided in the SFF file (PRIVATE).

    Will use the handle seek/tell functions.

    This works on ".srt1.00" and ".mft1.00" style Roche SFF index blocks.

    Roche SFF indices use base 255 not 256, meaning we see bytes in range the range
    0 to 254 only. This appears to be so that byte 0xFF (character 255) can be used
    as a marker character to separate entries (required if the read name lengths vary).

    Note that since only four bytes are used for the read offset, this is limited to
    255^4 bytes (nearly 4GB). If you try to use the Roche sfffile tool to combined SFF
    files beyound this limit, they issue a warning and ommit the index (and manifest).
    """
    number_of_reads, header_length, index_offset, index_length, xml_offset, \
    xml_size, read_index_offset, read_index_size = _sff_find_roche_index(handle)
    #Now parse the read index...    
    handle.seek(read_index_offset)
    start=chr(0)
    end=chr(255)
    fmt = ">5B"
    for read in range(number_of_reads) :
        #TODO - Be more aware of when the index should end?
        data = handle.read(6)
        while data[-1] != end :
            more = handle.read(1)
            if not more : raise ValueError("Premature end of file!")
            data += more
        name = data[:-6]
        off4, off3, off2, off1, off0 = struct.unpack(fmt, data[-6:-1])
        offset = off0 + 255*off1 + 65025*off2 + 16581375*off3
        if off4 :
            #Could in theory be used as a fifth piece of offset information,
            #i.e. offset =+ 4228250625L*off4, but testing the Roche tools this
            #is not the case. They simple don't support such large indexes.
            raise ValueError("Expected a null terminator to the read name.")
        yield name, offset
    if handle.tell() != read_index_offset + read_index_size :
        raise ValueError("Problem with index length? %i vs %i" \
                         % (handle.tell(), read_index_offset + read_index_size))

def _sff_read_seq_record(handle, number_of_flows_per_read, flow_chars,
                         key_sequence, alphabet, trim=False) :
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
    if clip_qual_left : clip_qual_left -= 1 #python counting
    if clip_adapter_left : clip_adapter_left -= 1 #python counting
    if read_header_length < 10 or read_header_length%8 != 0 :
        raise ValueError("Malformed read header, says length is %i" \
                         % read_header_length)
    #now the name and any padding (remainder of header)
    name = handle.read(name_length)
    padding = read_header_length - read_header_size - name_length
    if chr(0)*padding != handle.read(padding) :
        raise ValueError("Post name %i byte padding region contained data" \
                         % padding)
    #now the flowgram values, flowgram index, bases and qualities
    #NOTE - assuming flowgram_format==1, which means struct type H
    flow_values = struct.unpack(read_flow_fmt, handle.read(read_flow_size))
    temp_fmt = ">%iB" % seq_len # used for flow index and quals
    flow_index = struct.unpack(temp_fmt, handle.read(seq_len))
    seq = handle.read(seq_len)
    quals = list(struct.unpack(temp_fmt, handle.read(seq_len)))
    #now any padding...
    padding = (read_flow_size + seq_len*3)%8
    if padding :
        padding = 8 - padding
        if chr(0)*padding != handle.read(padding) :
            raise ValueError("Post quality %i byte padding region contained data" \
                             % padding)
    #Now build a SeqRecord
    if trim :
        seq = seq[clip_qual_left:clip_qual_right].upper()
        quals = quals[clip_qual_left:clip_qual_right]
        #Don't record the clipping values, flow etc, they make no sense now:
        annotations = {}
    else :
        #This use of mixed case mimics the Roche SFF tool's FASTA output
        seq = seq[:clip_qual_left].lower() + \
              seq[clip_qual_left:clip_qual_right].upper() + \
              seq[clip_qual_right:].lower()
        annotations = {"flow_values":flow_values,
                       "flow_index":flow_index,
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
    #TODO - adaptor clipping
    #Return the record and then continue...
    return record

#This is a generator function!
def SffIterator(handle, alphabet=Alphabet.generic_dna, trim=False) :
    """Iterate over Standard Flowgram Format (SFF) reads (as SeqRecord objects).

    handle - input file, an SFF file, e.g. from Roche 454 sequencing.
             This must NOT be opened in universal read lines mode!
    alphabet - optional alphabet, defaults to generic DNA.
    trim - should the sequences be trimmed?

    The resulting SeqRecord objects should match those from a paired
    FASTA and QUAL file converted from the SFF file using the Roche
    454 tool ssfinfo. i.e. The sequence will be mixed case, with the
    trim regions shown in lower case.
    """
    if isinstance(Alphabet._get_base_alphabet(alphabet),
                  Alphabet.ProteinAlphabet) :
        raise ValueError("Invalid alphabet, SFF files do not hold proteins.")
    if isinstance(Alphabet._get_base_alphabet(alphabet),
                  Alphabet.RNAAlphabet) :
        raise ValueError("Invalid alphabet, SFF files do not hold RNA.")
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
    for read in range(number_of_reads) :
        if index_offset and handle.tell() == index_offset :
            #import warnings
            #warnings.warn("Found SFF index before read %i (not at end)" % (read+1))
            offset = index_offset + index_length
            if offset % 8 :
                offset += 8 - (offset % 8)
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

#This is a generator function!
def _SffTrimIterator(handle, alphabet=Alphabet.generic_dna) :
    """Iterate over SFF reads (as SeqRecord objects) with trimming (PRIVATE)."""
    return SffIterator(handle, alphabet, trim=True)


class SffWriter(SequenceWriter) :
    def __init__(self, handle, index=True, xml=None):
        """Creates the writer object.

        handle - Output handle, ideally in binary write mode.
        index - Boolean argument, should we try and write an index?
        xml - Optional string argument, xml to be recorded in the index block.
        """
        if hasattr(handle,"mode") and "U" in handle.mode.upper() :
            raise ValueError("SFF files must NOT be opened in universal new "
                             "lines mode. Binary mode is recommended (although "
                             "on Unix the default mode is also fine).")
        elif hasattr(handle,"mode") and "B" not in handle.mode.upper() \
        and sys.platform == "win32":
            raise ValueError("SFF files must be opened in binary mode on Windows")
        self.handle = handle
        self._xml = xml
        if index :
            self._index = []
        else :
            self._index = None

    def write_file(self, records) :
        """Use this to write an entire file containing the given records."""
        try :
            self._number_of_reads = len(records)
        except TypeError :
            self._number_of_reads = 0 #dummy value
            if not hasattr(self.handle, "seek") \
            or not hasattr(self.handle, "tell") :
                raise ValueError("A handle with a seek/tell methods is required in "
                                 "order to record the total record count in the file "
                                 "header (once it is known at the end).")
        if self._index is not None and \
        not (hasattr(self.handle, "seek") and hasattr(self.handle, "tell")) :
            import warnings
            warnings.warn("A handle with a seek/tell methods is required in "
                          "order to record an SFF index.")
            self._index = None
        self._index_start = 0
        self._index_length = 0
        if not hasattr(records, "next") :
            records = iter(records)
        #Get the first record in order to find the flow information
        #we will need for the header.
        try :
            record = records.next()
        except StopIteration :
            record = None
        if record is None :
            #No records -> empty SFF file (or an error)?
            #We can't write a header without the flow information.
            #return 0
            raise ValueError("Need at least one record for SFF output")
        try :
            self._key_sequence = record.annotations["flow_key"]
            self._flow_chars = record.annotations["flow_chars"]
            self._number_of_flows_per_read = len(self._flow_chars)
        except KeyError :
            raise ValueError("Missing SFF flow information")
        self.write_header()
        self.write_record(record)
        count = 1
        for record in records :
            self.write_record(record)
            count += 1
        if self._number_of_reads == 0 :
            #Must go back and record the record count...
            offset = self.handle.tell()
            self.handle.seek(0)
            self._number_of_reads = count
            self.write_header()
            self.handle.seek(offset) #not essential?
        else :
            assert count == self._number_of_reads
        if self._index is not None :
            self._write_index()
        return count

    def _write_index(self) :
        assert len(self._index)==self._number_of_reads
        handle = self.handle
        self._index.sort()
        self._index_start = handle.tell() #need for header
        #XML...
        if isinstance(self._xml, str) :
            xml = self._xml
        else :
            from Bio import __version__
            xml = "<!-- This file was output with Biopython %s -->\n" % __version__
            xml += "<!-- This XML and index block attempts to mimic Roche SFF files -->\n"
            xml += "<!-- This file may be a combination of multiple SFF files etc -->\n"
        xml_len = len(xml)
        #Write to the file...
        fmt = ">I4BLL"
        fmt_size = struct.calcsize(fmt)
        handle.write(chr(0)*fmt_size + xml) #will come back later to fill this...
        fmt2 = ">6B"
        assert 6 == struct.calcsize(fmt2)
        self._index.sort()
        index_len = 0 #don't know yet!
        for name, offset in self._index :
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
            off3, off2, off1, off0 = off3//16581375, off2//65025, off1//255, off0
            assert off0 < 255 and off1 < 255 and off2 < 255 and off3 < 255, \
                   "%i -> %i %i %i %i" % (offset, off0, off1, off2, off3)
            handle.write(name + struct.pack(fmt2, 0, off3, off2, off1, off0, 255))
            index_len += len(name) + 6
        #Note any padding in not included:
        self._index_length = fmt_size + xml_len + index_len #need for header
        #Pad out to an 8 byte boundary (although I have noticed some
        #real Roche SFF files neglect to do this depsite their manual
        #suggesting this padding should be there):
        if self._index_length % 8 :
            padding = 8 - (self._index_length%8)
            handle.write(chr(0)*padding)
        else :
            padding = 0
        offset = handle.tell()
        assert offset == self._index_start + self._index_length + padding, \
               "%i vs %i + %i + %i"  %(offset, self._index_start, self._index_length, padding)
        #Must now go back and update the index header with index size...
        handle.seek(self._index_start)
        handle.write(struct.pack(fmt, 778921588, #magic number
                                 49,46,48,48, #Roche index version, "1.00"
                                 xml_len, index_len) + xml)
        #Must now go back and update the header...
        handle.seek(0)
        self.write_header()
        handle.seek(offset) #not essential?

    def write_header(self) :
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
        #According to the spec, the header_length field should be the total number
        #of bytes required by this set of header fields, and should be equal to
        #"31 + number_of_flows_per_read + key_length" rounded up to the next value
        #divisible by 8.
        if struct.calcsize(fmt) % 8 == 0 :
            padding = 0
        else :
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
        self.handle.write(header + chr(0)*padding)
        
    def write_record(self, record):
        """Write a single additional record to the output file.

        This assumes the header has been done.
        """
        #Basics
        name = record.id
        name_len = len(name)
        seq = str(record.seq).upper()
        seq_len = len(seq)
        #Qualities
        try :
            quals = record.letter_annotations["phred_quality"]
        except KeyError :
            raise ValueError("Missing PHRED qualities information")
        #Flow
        try :
            flow_values = record.annotations["flow_values"]
            flow_index = record.annotations["flow_index"]
            if self._key_sequence != record.annotations["flow_key"] \
            or self._flow_chars != record.annotations["flow_chars"] :
                raise ValueError("Records have inconsistent SFF flow data")
        except KeyError :
            raise ValueError("Missing SFF flow information")
        except AttributeError :
            raise ValueError("Header not written yet?")
        #Clipping
        try :
            clip_qual_left = record.annotations["clip_qual_left"]
            if clip_qual_left : clip_qual_left += 1
            clip_qual_right = record.annotations["clip_qual_right"]
            clip_adapter_left = record.annotations["clip_adapter_left"]
            if clip_adapter_left : clip_adapter_left += 1
            clip_adapter_right = record.annotations["clip_adapter_right"]
        except KeyError :
            raise ValueError("Missing SFF clipping information")

        #Capture information for index
        if self._index is not None :
            offset = self.handle.tell()
            #Check the position of the final record (before sort by name)
            #See comments earlier about how base 255 seems to be used.
            #This means the limit is 255**4 + 255**3 +255**2 + 255**1
            if offset > 4244897280L :
                import warnings
                warnings.warn("Read %s has file offset %i, which is too large to store "
                              "in the Roche SFF index structure. No index block will be "
                              "recorded." % (name, offset))
                #No point recoring the offsets now
                self._index = None
            else :
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
        if struct.calcsize(read_header_fmt) % 8 == 0 :
            padding = 0
        else :
            padding = 8 - (struct.calcsize(read_header_fmt) % 8)
        read_header_length = struct.calcsize(read_header_fmt) + padding
        assert read_header_length % 8 == 0
        data = struct.pack(read_header_fmt,
                           read_header_length,
                           name_len, seq_len,
                           clip_qual_left, clip_qual_right,
                           clip_adapter_left, clip_adapter_right,
                           name) + chr(0)*padding
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
        if padding :
            padding = 8 - padding
        self.handle.write(data + chr(0)*padding)


if __name__ == "__main__" :
    print "Running quick self test"
    filename = "../../Tests/Roche/E3MFGYR02_random_10_reads.sff"
    metadata = _sff_read_roche_index_xml(open(filename, "rb"))
    index1 = sorted(_sff_read_roche_index(open(filename, "rb")))
    index2 = sorted(_sff_do_slow_index(open(filename, "rb")))
    assert index1 == index2
    assert len(index1) == len(list(SffIterator(open(filename, "rb"))))
    from StringIO import StringIO
    assert len(index1) == len(list(SffIterator(StringIO(open(filename,"rb").read()))))

    if sys.platform != "win32" :
        assert len(index1) == len(list(SffIterator(open(filename, "r"))))
        index2 = sorted(_sff_read_roche_index(open(filename)))
        assert index1 == index2
        index2 = sorted(_sff_do_slow_index(open(filename)))
        assert index1 == index2
        assert len(index1) == len(list(SffIterator(open(filename))))
        assert len(index1) == len(list(SffIterator(StringIO(open(filename,"r").read()))))
        assert len(index1) == len(list(SffIterator(StringIO(open(filename).read()))))
                    
    sff = list(SffIterator(open(filename, "rb")))

    sff2 = list(SffIterator(open("../../Tests/Roche/E3MFGYR02_index_at_start.sff", "rb")))
    assert len(sff) == len(sff2)
    for old,new in zip(sff,sff2) :
        assert old.id == new.id

    sff2 = list(SffIterator(open("../../Tests/Roche/E3MFGYR02_index_in_middle.sff", "rb")))
    assert len(sff) == len(sff2)
    for old,new in zip(sff,sff2) :
        assert old.id == new.id
    
    sff_trim = list(SffIterator(open(filename, "rb"), trim=True))

    print _sff_read_roche_index_xml(open(filename, "rb"))

    from Bio import SeqIO
    filename = "../../Tests/Roche/E3MFGYR02_random_10_reads_no_trim.fasta"
    fasta_no_trim = list(SeqIO.parse(open(filename,"rU"), "fasta"))
    filename = "../../Tests/Roche/E3MFGYR02_random_10_reads_no_trim.qual"
    qual_no_trim = list(SeqIO.parse(open(filename,"rU"), "qual"))

    filename = "../../Tests/Roche/E3MFGYR02_random_10_reads.fasta"
    fasta_trim = list(SeqIO.parse(open(filename,"rU"), "fasta"))
    filename = "../../Tests/Roche/E3MFGYR02_random_10_reads.qual"
    qual_trim = list(SeqIO.parse(open(filename,"rU"), "qual"))

    for s, sT, f, q, fT, qT in zip(sff, sff_trim, fasta_no_trim, qual_no_trim, fasta_trim, qual_trim) :
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
    for record in SffIterator(open(filename,"rb")) :
        print record.id
    index1 = sorted(_sff_read_roche_index(open(filename, "rb")))
    index2 = sorted(_sff_do_slow_index(open(filename, "rb")))
    assert index1 == index2
    try :
        print _sff_read_roche_index_xml(open(filename, "rb"))
        assert False, "Should fail!"
    except ValueError :
        pass


    print "Done"
