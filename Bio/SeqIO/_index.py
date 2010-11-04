# Copyright 2009-2010 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Dictionary like indexing of sequence files (PRIVATE).

You are not expected to access this module, or any of its code, directly. This
is all handled internally by the Bio.SeqIO.index(...) function which is the
public interface for this functionality.

The basic idea is that we scan over a sequence file, looking for new record
markers. We then try and extract the string that Bio.SeqIO.parse/read would
use as the record id, ideally without actually parsing the full record. We
then use a subclassed Python dictionary to record the file offset for the
record start against the record id.

Note that this means full parsing is on demand, so any invalid or problem
record may not trigger an exception until it is accessed. This is by design.

This means our dictionary like objects have in memory ALL the keys (all the
record identifiers), which shouldn't be a problem even with second generation
sequencing. If this is an issue later on, storing the keys and offsets in a
temp lookup file might be one idea (e.g. using SQLite or an OBDA style index).
"""

import re
from Bio import SeqIO
from Bio import Alphabet

class _IndexedSeqFileDict(dict):
    """Read only dictionary interface to a sequential sequence file.

    Keeps the keys in memory, reads the file to access entries as
    SeqRecord objects using Bio.SeqIO for parsing them. This approach
    is memory limited, but will work even with millions of sequences.

    Note - as with the Bio.SeqIO.to_dict() function, duplicate keys
    (record identifiers by default) are not allowed. If this happens,
    a ValueError exception is raised.

    By default the SeqRecord's id string is used as the dictionary
    key. This can be changed by suppling an optional key_function,
    a callback function which will be given the record id and must
    return the desired key. For example, this allows you to parse
    NCBI style FASTA identifiers, and extract the GI number to use
    as the dictionary key.

    Note that this dictionary is essentially read only. You cannot
    add or change values, pop values, nor clear the dictionary.
    """
    def __init__(self, filename, format, alphabet, key_function):
        #Use key_function=None for default value
        dict.__init__(self) #init as empty dict!
        if format in SeqIO._BinaryFormats:
            mode = "rb"
        else:
            mode = "rU"
        self._handle = open(filename, mode)
        self._alphabet = alphabet
        self._format = format
        self._key_function = key_function
        if key_function:
            offset_iter = ((key_function(k),o) for (k,o) in self._build())
        else:
            offset_iter = self._build()
        for key, offset in offset_iter:
            if key in self:
                raise ValueError("Duplicate key '%s'" % key)
            else:
                dict.__setitem__(self, key, offset)
    
    def _build(self):
        """Actually scan the file identifying records and offsets (PRIVATE).
        
        Returns an iterator giving tuples of record names and their offsets.
        """
        pass

    def __repr__(self):
        return "SeqIO.index('%s', '%s', alphabet=%s, key_function=%s)" \
               % (self._handle.name, self._format,
                  repr(self._alphabet), self._key_function)

    def __str__(self):
        if self:
            return "{%s : SeqRecord(...), ...}" % repr(self.keys()[0])
        else:
            return "{}"

    if hasattr(dict, "iteritems"):
        #Python 2, use iteritems but not items etc
        def values(self):
            """Would be a list of the SeqRecord objects, but not implemented.

            In general you can be indexing very very large files, with millions
            of sequences. Loading all these into memory at once as SeqRecord
            objects would (probably) use up all the RAM. Therefore we simply
            don't support this dictionary method.
            """
            raise NotImplementedError("Due to memory concerns, when indexing a "
                                      "sequence file you cannot access all the "
                                      "records at once.")

        def items(self):
            """Would be a list of the (key, SeqRecord) tuples, but not implemented.

            In general you can be indexing very very large files, with millions
            of sequences. Loading all these into memory at once as SeqRecord
            objects would (probably) use up all the RAM. Therefore we simply
            don't support this dictionary method.
            """
            raise NotImplementedError("Due to memory concerns, when indexing a "
                                      "sequence file you cannot access all the "
                                      "records at once.")

        def iteritems(self):
            """Iterate over the (key, SeqRecord) items."""
            for key in self.__iter__():
                yield key, self.__getitem__(key)
    else:
        #Python 3 - define items and values as iterators
        def items(self):
            """Iterate over the (key, SeqRecord) items."""
            for key in self.__iter__():
                yield key, self.__getitem__(key)

        def values(self):
            """Iterate over the SeqRecord items."""
            for key in self.__iter__():
                yield self.__getitem__(key)


    def __getitem__(self, key):
        """x.__getitem__(y) <==> x[y]"""
        #Should be done by each sub-class
        raise NotImplementedError("Not implemented for this file format (yet).")

    def get(self, k, d=None):
        """D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None."""
        try:
            return self.__getitem__(k)
        except KeyError:
            return d

    def get_raw(self, key):
        """Similar to the get method, but returns the record as a raw string.

        If the key is not found, a KeyError exception is raised.

        NOTE - This functionality is not supported for every file format.
        """
        #Should be done by each sub-class (if possible)
        raise NotImplementedError("Not available for this file format.")

    def __setitem__(self, key, value):
        """Would allow setting or replacing records, but not implemented."""
        raise NotImplementedError("An indexed a sequence file is read only.")
    
    def update(self, **kwargs):
        """Would allow adding more values, but not implemented."""
        raise NotImplementedError("An indexed a sequence file is read only.")

    
    def pop(self, key, default=None):
        """Would remove specified record, but not implemented."""
        raise NotImplementedError("An indexed a sequence file is read only.")
    
    def popitem(self):
        """Would remove and return a SeqRecord, but not implemented."""
        raise NotImplementedError("An indexed a sequence file is read only.")

    
    def clear(self):
        """Would clear dictionary, but not implemented."""
        raise NotImplementedError("An indexed a sequence file is read only.")

    def fromkeys(self, keys, value=None):
        """A dictionary method which we don't implement."""
        raise NotImplementedError("An indexed a sequence file doesn't "
                                  "support this.")

    def copy(self):
        """A dictionary method which we don't implement."""
        raise NotImplementedError("An indexed a sequence file doesn't "
                                  "support this.")



####################
# Special indexers #
####################

# Anything where the records cannot be read simply by parsing from
# the record start. For example, anything requiring information from
# a file header - e.g. SFF files where we would need to know the
# number of flows.

class SffDict(_IndexedSeqFileDict) :
    """Indexed dictionary like access to a Standard Flowgram Format (SFF) file."""
    def _build(self):
        """Load any index block in the file, or build it the slow way (PRIVATE)."""
        if self._alphabet is None:
            self._alphabet = Alphabet.generic_dna
        handle = self._handle
        #Record the what we'll need for parsing a record given its offset
        header_length, index_offset, index_length, number_of_reads, \
        self._flows_per_read, self._flow_chars, self._key_sequence \
            = SeqIO.SffIO._sff_file_header(handle)
        if index_offset and index_length:
            #There is an index provided, try this the fast way:
            try :
                for name, offset in SeqIO.SffIO._sff_read_roche_index(handle) :
                    yield name, offset
                assert len(self) == number_of_reads, \
                       "Indexed %i records, expected %i" \
                       % (len(self), number_of_reads)
                return
            except ValueError, err :
                import warnings
                warnings.warn("Could not parse the SFF index: %s" % err)
                assert len(self)==0, "Partially populated index"
                handle.seek(0)
        else :
            #TODO - Remove this debug warning?
            import warnings
            warnings.warn("No SFF index, doing it the slow way")
        #Fall back on the slow way!
        for name, offset in SeqIO.SffIO._sff_do_slow_index(handle) :
            yield name, offset
        assert len(self) == number_of_reads, \
               "Indexed %i records, expected %i" % (len(self), number_of_reads)

    def __getitem__(self, key) :
        handle = self._handle
        handle.seek(dict.__getitem__(self, key))
        return SeqIO.SffIO._sff_read_seq_record(handle,
                                                self._flows_per_read,
                                                self._flow_chars,
                                                self._key_sequence,
                                                self._alphabet)


class SffTrimmedDict(SffDict) :
    def __getitem__(self, key) :
        handle = self._handle
        handle.seek(dict.__getitem__(self, key))
        return SeqIO.SffIO._sff_read_seq_record(handle,
                                                self._flows_per_read,
                                                self._flow_chars,
                                                self._key_sequence,
                                                self._alphabet,
                                                trim=True)


###################
# Simple indexers #
###################

class SequentialSeqFileDict(_IndexedSeqFileDict):
    """Indexed dictionary like access to most sequential sequence files."""
    def __init__(self, filename, format, alphabet, key_function):
        _IndexedSeqFileDict.__init__(self, filename, format, alphabet, key_function)
        #Load the parser class/function once an avoid the dict lookup in each
        #__getitem__ call:
        i = SeqIO._FormatToIterator[format]
        #The following alphabet code is a bit nasty... duplicates logic in
        #Bio.SeqIO.parse()
        if alphabet is None:
            def _parse():
                """Dynamically generated parser function (PRIVATE)."""
                return i(self._handle).next()
        else:
            #TODO - Detect alphabet support ONCE at __init__
            def _parse():
                """Dynamically generated parser function (PRIVATE)."""
                try:
                    return i(self._handle, alphabet=alphabet).next()
                except TypeError:
                    return SeqIO._force_alphabet(i(self._handle),
                                                 alphabet).next()
        self._parse = _parse

    def _build(self):
        handle = self._handle
        marker = {"ace" : "CO ",
                  "fasta": ">",
                  "phd" : "BEGIN_SEQUENCE",
                  "pir" : ">..;",
                  "qual": ">",
                   }[self._format]
        marker_re = re.compile("^%s" % marker)
        marker_offset = len(marker)
        self._marker_re = marker_re #saved for the get_raw method
        while True:
            offset = handle.tell()
            line = handle.readline()
            if marker_re.match(line):
                #Here we can assume the record.id is the first word after the
                #marker. This is generally fine... but not for GenBank, EMBL, Swiss
                yield line[marker_offset:].strip().split(None, 1)[0], offset
            elif not line:
                #End of file
                break

    def __getitem__(self, key):
        """x.__getitem__(y) <==> x[y]"""
        handle = self._handle
        handle.seek(dict.__getitem__(self, key))
        record = self._parse()
        if self._key_function:
            assert self._key_function(record.id) == key, \
                   "Requested key %s, found record.id %s which has key %s" \
                   % (repr(key), repr(record.id),
                      repr(self._key_function(record.id)))
        else:
            assert record.id == key, \
                   "Requested key %s, found record.id %s" \
                   % (repr(key), repr(record.id))
        return record

    def get_raw(self, key):
        """Similar to the get method, but returns the record as a raw string."""
        #For non-trivial file formats this must be over-ridden in the subclass
        handle = self._handle
        marker_re = self._marker_re
        handle.seek(dict.__getitem__(self, key))
        data = handle.readline()
        while True:
            line = handle.readline()
            if marker_re.match(line) or not line:
                #End of file, or start of next record => end of this record
                break
            data += line
        return data


#######################################
# Fiddly indexers: GenBank, EMBL, ... #
#######################################

class GenBankDict(SequentialSeqFileDict):
    """Indexed dictionary like access to a GenBank file."""
    def _build(self):
        handle = self._handle
        marker_re = re.compile("^LOCUS ")
        self._marker_re = marker_re #saved for the get_raw method
        while True:
            offset = handle.tell()
            line = handle.readline()
            if marker_re.match(line):
                #We cannot assume the record.id is the first word after LOCUS,
                #normally the first entry on the VERSION or ACCESSION line is used.
                key = None
                done = False
                while not done:
                    line = handle.readline()
                    if line.startswith("ACCESSION "):
                        key = line.rstrip().split()[1]
                    elif line.startswith("VERSION "):
                        version_id = line.rstrip().split()[1]
                        if version_id.count(".")==1 and version_id.split(".")[1].isdigit():
                            #This should mimic the GenBank parser...
                            key = version_id
                            done = True
                            break
                    elif line.startswith("FEATURES ") \
                    or line.startswith("ORIGIN ") \
                    or line.startswith("//") \
                    or marker_re.match(line) \
                    or not line:
                        done = True
                        break
                if not key:
                    raise ValueError("Did not find ACCESSION/VERSION lines")
                yield key, offset
            elif not line:
                #End of file
                break


class EmblDict(SequentialSeqFileDict):
    """Indexed dictionary like access to an EMBL file."""
    def _build(self):
        handle = self._handle
        marker_re = re.compile("^ID ")
        self._marker_re = marker_re #saved for the get_raw method
        while True:
            offset = handle.tell()
            line = handle.readline()
            if marker_re.match(line):
                #We cannot assume the record.id is the first word after ID,
                #normally the SV line is used.
                if line[2:].count(";") == 6:
                    #Looks like the semi colon separated style introduced in 2006
                    parts = line[3:].rstrip().split(";")
                    if parts[1].strip().startswith("SV "):
                        #The SV bit gives the version
                        key = "%s.%s" \
                              % (parts[0].strip(), parts[1].strip().split()[1])
                    else:
                        key = parts[0].strip()
                elif line[2:].count(";") == 3:
                    #Looks like the pre 2006 style, take first word only
                    key = line[3:].strip().split(None,1)[0]
                else:
                    raise ValueError('Did not recognise the ID line layout:\n' + line)
                while True:
                    line = handle.readline()
                    if line.startswith("SV "):
                        key = line.rstrip().split()[1]
                        break
                    elif line.startswith("FH ") \
                    or line.startswith("FT ") \
                    or line.startswith("SQ ") \
                    or line.startswith("//") \
                    or marker_re.match(line) \
                    or not line:
                        break
                yield key, offset
            elif not line:
                #End of file
                break


class SwissDict(SequentialSeqFileDict):
    """Indexed dictionary like access to a SwissProt file."""
    def _build(self):
        handle = self._handle
        marker_re = re.compile("^ID ")
        self._marker_re = marker_re #saved for the get_raw method
        while True:
            offset = handle.tell()
            line = handle.readline()
            if marker_re.match(line):
                #We cannot assume the record.id is the first word after ID,
                #normally the following AC line is used.
                line = handle.readline()
                assert line.startswith("AC ")
                key = line[3:].strip().split(";")[0].strip()
                yield key, offset
            elif not line:
                #End of file
                break


class UniprotDict(SequentialSeqFileDict):
    """Indexed dictionary like access to a UniProt XML file."""
    def _build(self):
        handle = self._handle
        marker_re = re.compile("^<entry ") #Assumes start of line!
        self._marker_re = marker_re #saved for the get_raw method
        while True:
            offset = handle.tell()
            line = handle.readline()
            if marker_re.match(line):
                #We expect the next line to be <accession>xxx</accession>
                #but allow it to be later on within the <entry>
                key = None
                done = False
                while True:
                    line = handle.readline()
                    if line.startswith("<accession>"):
                        assert "</accession>" in line, line
                        key = line[11:].split("<")[0]
                        break
                    elif "</entry>" in line:
                        break
                    elif marker_re.match(line) or not line:
                        #Start of next record or end of file
                        raise ValueError("Didn't find end of record")
                if not key:
                    raise ValueError("Did not find <accession> line")
                yield key, offset
            elif not line:
                #End of file
                break
    
    def get_raw(self, key):
        """Similar to the get method, but returns the record as a raw string."""
        handle = self._handle
        marker_re = self._marker_re
        handle.seek(dict.__getitem__(self, key))
        data = handle.readline()
        while True:
            line = handle.readline()
            i = line.find("</entry>")
            if i != -1:
                data += line[:i+8]
                break
            if marker_re.match(line) or not line:
                #End of file, or start of next record
                raise ValueError("Didn't find end of record")
            data += line
        return data

    def __getitem__(self, key) :
        #TODO - Can we handle this directly in the parser?
        #This is a hack - use get_raw for <entry>...</entry> and wrap it with
        #the apparently required XML header and footer.
        data = """<?xml version='1.0' encoding='UTF-8'?>
        <uniprot xmlns="http://uniprot.org/uniprot"
        xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
        xsi:schemaLocation="http://uniprot.org/uniprot
        http://www.uniprot.org/support/docs/uniprot.xsd">
        %s
        </uniprot>
        """ % self.get_raw(key)
        #TODO - For consistency, this function should not accept a string:
        return SeqIO.UniprotIO.UniprotIterator(data).next()


class IntelliGeneticsDict(SequentialSeqFileDict):
    """Indexed dictionary like access to a IntelliGenetics file."""
    def _build(self):
        handle = self._handle
        marker_re = re.compile("^;")
        self._marker_re = marker_re #saved for the get_raw method
        while True:
            offset = handle.tell()
            line = handle.readline()
            if marker_re.match(line):
                #Now look for the first line which doesn't start ";"
                while True:
                    line = handle.readline()
                    if line[0] != ";" and line.strip():
                        key = line.split()[0]
                        yield key, offset
                        break
                    if not line:
                        raise ValueError("Premature end of file?")
            elif not line:
                #End of file
                break


class TabDict(SequentialSeqFileDict):
    """Indexed dictionary like access to a simple tabbed file."""
    def _build(self):
        handle = self._handle
        while True:
            offset = handle.tell()
            line = handle.readline()
            if not line : break #End of file
            try:
                key = line.split("\t")[0]
            except ValueError, err:
                if not line.strip():
                    #Ignore blank lines
                    continue
                else:
                    raise err
            else:
                yield key, offset

    def get_raw(self, key):
        """Like the get method, but returns the record as a raw string."""
        handle = self._handle
        handle.seek(dict.__getitem__(self, key))
        return handle.readline()


##########################
# Now the FASTQ indexers #
##########################
         
class FastqDict(SequentialSeqFileDict):
    """Indexed dictionary like access to a FASTQ file (any supported variant).
    
    With FASTQ the records all start with a "@" line, but so can quality lines.
    Note this will cope with line-wrapped FASTQ files.
    """
    def _build(self):
        handle = self._handle
        pos = handle.tell()
        line = handle.readline()
        if not line:
            #Empty file!
            return
        if line[0] != "@":
            raise ValueError("Problem with FASTQ @ line:\n%s" % repr(line))
        while line:
            #assert line[0]=="@"
            #This record seems OK (so far)
            yield line[1:].rstrip().split(None, 1)[0], pos
            #Find the seq line(s)
            seq_len = 0
            while line:
                line = handle.readline()
                if line.startswith("+") : break
                seq_len += len(line.strip())
            if not line:
                raise ValueError("Premature end of file in seq section")
            #assert line[0]=="+"
            #Find the qual line(s)
            qual_len = 0
            while line:
                if seq_len == qual_len:
                    #Should be end of record...
                    pos = handle.tell()
                    line = handle.readline()
                    if line and line[0] != "@":
                        ValueError("Problem with line %s" % repr(line))
                    break
                else:
                    line = handle.readline()
                    qual_len += len(line.strip())
            if seq_len != qual_len:
                raise ValueError("Problem with quality section")
        #print "EOF"

    def get_raw(self, key):
        """Similar to the get method, but returns the record as a raw string."""
        #TODO - Refactor this and the __init__ method to reduce code duplication?
        handle = self._handle
        handle.seek(dict.__getitem__(self, key))
        line = handle.readline()
        data = line
        if line[0] != "@":
            raise ValueError("Problem with FASTQ @ line:\n%s" % repr(line))
        identifier = line[1:].rstrip().split(None, 1)[0]
        if self._key_function:
            identifier = self._key_function(identifier)
        if key != identifier:
            raise ValueError("Key did not match")
        #Find the seq line(s)
        seq_len = 0
        while line:
            line = handle.readline()
            data += line
            if line.startswith("+") : break
            seq_len += len(line.strip())
        if not line:
            raise ValueError("Premature end of file in seq section")
        assert line[0]=="+"
        #Find the qual line(s)
        qual_len = 0
        while line:
            if seq_len == qual_len:
                #Should be end of record...
                pos = handle.tell()
                line = handle.readline()
                if line and line[0] != "@":
                    ValueError("Problem with line %s" % repr(line))
                break
            else:
                line = handle.readline()
                data += line
                qual_len += len(line.strip())
        if seq_len != qual_len:
            raise ValueError("Problem with quality section")
        return data


###############################################################################

_FormatToIndexedDict = {"ace" : SequentialSeqFileDict,
                        "embl" : EmblDict,
                        "fasta" : SequentialSeqFileDict,
                        "fastq" : FastqDict, #Class handles all three variants
                        "fastq-sanger" : FastqDict, #alias of the above
                        "fastq-solexa" : FastqDict,
                        "fastq-illumina" : FastqDict,
                        "genbank" : GenBankDict,
                        "gb" : GenBankDict, #alias of the above
                        "ig" : IntelliGeneticsDict,
                        "imgt" : EmblDict,
                        "phd" : SequentialSeqFileDict,
                        "pir" : SequentialSeqFileDict,
                        "sff" : SffDict,
                        "sff-trim" : SffTrimmedDict,
                        "swiss" : SwissDict,
                        "tab" : TabDict,
                        "qual" : SequentialSeqFileDict,
                        "uniprot-xml" : UniprotDict,
                        }

