# Copyright 2009 by Peter Cock.  All rights reserved.
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
    def __init__(self, filename, alphabet, key_function, mode="rU"):
        #Use key_function=None for default value
        dict.__init__(self) #init as empty dict!
        self._handle = open(filename, mode)
        self._alphabet = alphabet
        self._format = ""
        self._key_function = key_function
        #Now scan it in a subclassed method, and set the format!

    def __repr__(self):
        return "SeqIO.index('%s', '%s', alphabet=%s, key_function=%s)" \
               % (self._handle.name, self._format,
                  repr(self._alphabet), self._key_function)

    def __str__(self):
        if self:
            return "{%s : SeqRecord(...), ...}" % repr(self.keys()[0])
        else:
            return "{}"

    def _record_key(self, identifier, seek_position):
        """Used by subclasses to record file offsets for identifiers (PRIVATE).

        This will apply the key_function (if given) to map the record id
        string to the desired key.

        This will raise a ValueError if a key (record id string) occurs
        more than once.
        """
        if self._key_function:
            key = self._key_function(identifier)
        else:
            key = identifier
        if key in self:
            raise ValueError("Duplicate key '%s'" % key)
        else:
            dict.__setitem__(self, key, seek_position)

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

    def __getitem__(self, key):
        """x.__getitem__(y) <==> x[y]"""
        #For non-trivial file formats this must be over-ridden in the subclass
        handle = self._handle
        handle.seek(dict.__getitem__(self, key))
        record = SeqIO.parse(handle, self._format, self._alphabet).next()
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

    def get(self, k, d=None):
        """D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None."""
        try:
            return self.__getitem__(k)
        except KeyError:
            return d

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

###################
# Simple indexers #
###################

class _SequentialSeqFileDict(_IndexedSeqFileDict):
    """Subclass for easy cases (PRIVATE)."""
    def __init__(self, filename, alphabet, key_function, format, marker):
        _IndexedSeqFileDict.__init__(self, filename, alphabet, key_function)
        self._format = format
        handle = self._handle
        marker_re = re.compile("^%s" % marker)
        marker_offset = len(marker)
        while True:
            offset = handle.tell()
            line = handle.readline()
            if not line : break #End of file
            if marker_re.match(line):
                #Here we can assume the record.id is the first word after the
                #marker. This is generally fine... but not for GenBank, EMBL, Swiss
                self._record_key(line[marker_offset:].strip().split(None,1)[0], offset)

class FastaDict(_SequentialSeqFileDict):
    """Indexed dictionary like access to a FASTA file."""
    def __init__(self, filename, alphabet, key_function):
        _SequentialSeqFileDict.__init__(self, filename, alphabet, key_function,
                                        "fasta", ">")

class QualDict(_SequentialSeqFileDict):
    """Indexed dictionary like access to a QUAL file."""
    def __init__(self, filename, alphabet, key_function):
        _SequentialSeqFileDict.__init__(self, filename, alphabet, key_function,
                                        "qual", ">")

class PirDict(_SequentialSeqFileDict):
    """Indexed dictionary like access to a PIR/NBRF file."""
    def __init__(self, filename, alphabet, key_function):
        _SequentialSeqFileDict.__init__(self, filename, alphabet, key_function,
                                        "pir", ">..;")

class PhdDict(_SequentialSeqFileDict):
    """Indexed dictionary like access to a PHD (PHRED) file."""
    def __init__(self, filename, alphabet, key_function):
        _SequentialSeqFileDict.__init__(self, filename, alphabet, key_function,
                                        "phd", "BEGIN_SEQUENCE")

class AceDict(_SequentialSeqFileDict):
    """Indexed dictionary like access to an ACE file."""
    def __init__(self, filename, alphabet, key_function):
        _SequentialSeqFileDict.__init__(self, filename, alphabet, key_function,
                                        "ace", "CO ")


#######################################
# Fiddly indexers: GenBank, EMBL, ... #
#######################################

class GenBankDict(_IndexedSeqFileDict):
    """Indexed dictionary like access to a GenBank file."""
    def __init__(self, filename, alphabet, key_function):
        _IndexedSeqFileDict.__init__(self, filename, alphabet, key_function)
        self._format = "genbank"
        handle = self._handle
        marker_re = re.compile("^LOCUS ")
        while True:
            offset = handle.tell()
            line = handle.readline()
            if not line : break #End of file
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
                            #This should mimics the GenBank parser...
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
                self._record_key(key, offset)

class EmblDict(_IndexedSeqFileDict):
    """Indexed dictionary like access to an EMBL file."""
    def __init__(self, filename, alphabet, key_function):
        _IndexedSeqFileDict.__init__(self, filename, alphabet, key_function)
        self._format = "embl"
        handle = self._handle
        marker_re = re.compile("^ID ")
        while True:
            offset = handle.tell()
            line = handle.readline()
            if not line : break #End of file
            if marker_re.match(line):
                #We cannot assume the record.id is the first word after ID,
                #normally the SV line is used.
                parts = line[3:].rstrip().split(";")
                if parts[1].strip().startswith("SV "):
                    #The SV bit gives the version
                    key = "%s.%s" % (parts[0].strip(),parts[1].strip().split()[1])
                else:
                    key = parts[0].strip()
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
                self._record_key(key, offset)

class SwissDict(_IndexedSeqFileDict):
    """Indexed dictionary like access to a SwissProt file."""
    def __init__(self, filename, alphabet, key_function):
        _IndexedSeqFileDict.__init__(self, filename, alphabet, key_function)
        self._format = "swiss"
        handle = self._handle
        marker_re = re.compile("^ID ")
        while True:
            offset = handle.tell()
            line = handle.readline()
            if not line : break #End of file
            if marker_re.match(line):
                #We cannot assume the record.id is the first word after ID,
                #normally the following AC line is used.
                line = handle.readline()
                assert line.startswith("AC ")
                key = line[3:].strip().split(";")[0].strip()
                self._record_key(key, offset)

class IntelliGeneticsDict(_IndexedSeqFileDict):
    """Indexed dictionary like access to a IntelliGenetics file."""
    def __init__(self, filename, alphabet, key_function):
        _IndexedSeqFileDict.__init__(self, filename, alphabet, key_function)
        self._format = "ig"
        handle = self._handle
        marker_re = re.compile("^;")
        while True:
            offset = handle.tell()
            line = handle.readline()
            if not line : break #End of file
            if marker_re.match(line):
                #Now look for the first line which doesn't start ";"
                while True:
                    line = handle.readline()
                    if not line:
                        raise ValueError("Premature end of file?")
                    if line[0] != ";" and line.strip():
                        key = line.split()[0]
                        self._record_key(key, offset)
                        break

class TabDict(_IndexedSeqFileDict):
    """Indexed dictionary like access to a simple tabbed file."""
    def __init__(self, filename, alphabet, key_function):
        _IndexedSeqFileDict.__init__(self, filename, alphabet, key_function)
        self._format = "tab"
        handle = self._handle
        while True:
            offset = handle.tell()
            line = handle.readline()
            if not line : break #End of file
            try:
                key, rest = line.split("\t")
            except ValueError, err:
                if not line.strip():
                    #Ignore blank lines
                    continue
                else:
                    raise err
            else:
                self._record_key(key, offset)

##########################
# Now the FASTQ indexers #
##########################
         
class _FastqSeqFileDict(_IndexedSeqFileDict):
    """Subclass for easy cases (PRIVATE).

    With FASTQ the records all start with a "@" line, but so too can some
    quality lines. Note this will cope with line-wrapped FASTQ files.
    """
    def __init__(self, filename, alphabet, key_function, fastq_format):
        _IndexedSeqFileDict.__init__(self, filename, alphabet, key_function)
        self._format = fastq_format
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
            self._record_key(line[1:].rstrip().split(None,1)[0],pos)
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
                    if line and line[0]!="@":
                        ValueError("Problem with line %s" % repr(line))
                    break
                else:
                    line = handle.readline()
                    qual_len += len(line.strip())
            if seq_len != qual_len:
                raise ValueError("Problem with quality section")
        #print "EOF"

class FastqSangerDict(_FastqSeqFileDict):
    """Indexed dictionary like access to a standard Sanger FASTQ file."""
    def __init__(self, filename, alphabet, key_function):
        _FastqSeqFileDict.__init__(self, filename, alphabet, key_function,
                                   "fastq-sanger")

class FastqSolexaDict(_FastqSeqFileDict):
    """Indexed dictionary like access to a Solexa (or early Illumina) FASTQ file."""
    def __init__(self, filename, alphabet, key_function):
        _FastqSeqFileDict.__init__(self, filename, alphabet, key_function,
                                   "fastq-solexa")

class FastqIlluminaDict(_FastqSeqFileDict):
    """Indexed dictionary like access to a Illumina 1.3+ FASTQ file."""
    def __init__(self, filename, alphabet, key_function):
        _FastqSeqFileDict.__init__(self, filename, alphabet, key_function,
                                   "fastq-illumina")

###############################################################################

_FormatToIndexedDict = {"ace" : AceDict,
                        "embl" : EmblDict,
                        "fasta" : FastaDict,
                        "fastq" : FastqSangerDict,
                        "fastq-sanger" : FastqSangerDict, #alias of the above
                        "fastq-solexa" : FastqSolexaDict,
                        "fastq-illumina" : FastqIlluminaDict,
                        "genbank" : GenBankDict,
                        "gb" : GenBankDict, #alias of the above
                        "ig" : IntelliGeneticsDict,
                        "phd" : PhdDict,
                        "pir" : PirDict,
                        "swiss" : SwissDict,
                        "tab" : TabDict,
                        "qual" : QualDict
                        }

