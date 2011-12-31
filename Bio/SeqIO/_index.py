# Copyright 2009-2011 by Peter Cock.  All rights reserved.
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

import os
try:
    from collections import UserDict as _dict_base
except ImportError:
    from UserDict import DictMixin as _dict_base
import re
import itertools
from StringIO import StringIO

try:
    from sqlite3 import dbapi2 as _sqlite
    from sqlite3 import IntegrityError as _IntegrityError
    from sqlite3 import OperationalError as _OperationalError
except ImportError:
    #Not expected to be present on Python 2.4, ignore it
    #and at least offer Bio.SeqIO.index() functionality
    _sqlite = None
    pass

from Bio._py3k import _bytes_to_string, _as_bytes, _as_string

from Bio import SeqIO
from Bio import Alphabet

class _IndexedSeqFileDict(_dict_base):
    """Read only dictionary interface to a sequential sequence file.

    Keeps the keys and associated file offsets in memory, reads the file to
    access entries as SeqRecord objects using Bio.SeqIO for parsing them.
    This approach is memory limited, but will work even with millions of
    sequences.

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
        try:
            proxy_class = _FormatToRandomAccess[format]
        except KeyError:
            raise ValueError("Unsupported format '%s'" % format)
        random_access_proxy = proxy_class(filename, format, alphabet)
        self._proxy = random_access_proxy
        self._key_function = key_function
        if key_function:
            offset_iter = ((key_function(k),o,l) for (k,o,l) in random_access_proxy)
        else:
            offset_iter = random_access_proxy
        offsets = {}
        for key, offset, length in offset_iter:
            #Note - we don't store the length because I want to minimise the
            #memory requirements. With the SQLite backend the length is kept
            #and is used to speed up the get_raw method (by about 3 times).
            if key in offsets:
                self._proxy._handle.close()
                raise ValueError("Duplicate key '%s'" % key)
            else:
                offsets[key] = offset
        self._offsets = offsets
    
    def __repr__(self):
        return "SeqIO.index(%r, %r, alphabet=%r, key_function=%r)" \
               % (self._proxy._handle.name, self._proxy._format,
                  self._proxy._alphabet, self._key_function)

    def __str__(self):
        if self:
            return "{%s : SeqRecord(...), ...}" % repr(self.keys()[0])
        else:
            return "{}"

    def __contains__(self, key) :
        return key in self._offsets
        
    def __len__(self):
        """How many records are there?"""
        return len(self._offsets)

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

        def keys(self) :
            """Return a list of all the keys (SeqRecord identifiers)."""
            #TODO - Stick a warning in here for large lists? Or just refuse?
            return self._offsets.keys()

        def itervalues(self):
            """Iterate over the SeqRecord) items."""
            for key in self.__iter__():
                yield self.__getitem__(key)

        def iteritems(self):
            """Iterate over the (key, SeqRecord) items."""
            for key in self.__iter__():
                yield key, self.__getitem__(key)
        
        def iterkeys(self):
            """Iterate over the keys."""
            return self.__iter__()

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

        def keys(self):
            """Iterate over the keys."""
            return self.__iter__()

    def __iter__(self):
        """Iterate over the keys."""
        return iter(self._offsets)
        
    def __getitem__(self, key):
        """x.__getitem__(y) <==> x[y]"""
        #Pass the offset to the proxy
        record = self._proxy.get(self._offsets[key])
        if self._key_function:
            key2 = self._key_function(record.id)
        else:
            key2 = record.id
        if key != key2:
            raise ValueError("Key did not match (%s vs %s)" % (key, key2))
        return record

    def get(self, k, d=None):
        """D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None."""
        try:
            return self.__getitem__(k)
        except KeyError:
            return d

    def get_raw(self, key):
        """Similar to the get method, but returns the record as a raw string.

        If the key is not found, a KeyError exception is raised.

        Note that on Python 3 a bytes string is returned, not a typical
        unicode string.

        NOTE - This functionality is not supported for every file format.
        """
        #Pass the offset to the proxy
        return self._proxy.get_raw(self._offsets[key])

    def __setitem__(self, key, value):
        """Would allow setting or replacing records, but not implemented."""
        raise NotImplementedError("An indexed a sequence file is read only.")
    
    def update(self, *args, **kwargs):
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


class _SQLiteManySeqFilesDict(_IndexedSeqFileDict):
    """Read only dictionary interface to many sequential sequence files.

    Keeps the keys, file-numbers and offsets in an SQLite database. To access
    a record by key, reads from the offset in the approapriate file using
    Bio.SeqIO for parsing.
    
    There are OS limits on the number of files that can be open at once,
    so a pool are kept. If a record is required from a closed file, then
    one of the open handles is closed first.
    """
    def __init__(self, index_filename, filenames, format, alphabet,
                 key_function, max_open=10):
        random_access_proxies = {}
        #TODO? - Don't keep filename list in memory (just in DB)?
        #Should save a chunk of memory if dealing with 1000s of files.
        #Furthermore could compare a generator to the DB on reloading
        #(no need to turn it into a list)
        if not _sqlite:
            #Hack for Python 2.4 (of if Python is compiled without it)
            from Bio import MissingPythonDependencyError
            raise MissingPythonDependencyError("Requires sqlite3, which is "
                                               "included Python 2.5+")
        if filenames is not None:
            filenames = list(filenames) #In case it was a generator
        if os.path.isfile(index_filename):
            #Reuse the index.
            con = _sqlite.connect(index_filename)
            self._con = con
            #Check the count...
            try:
                count, = con.execute("SELECT value FROM meta_data WHERE key=?;",
                                     ("count",)).fetchone()
                self._length = int(count)
                if self._length == -1:
                    con.close()
                    raise ValueError("Unfinished/partial database")
                count, = con.execute("SELECT COUNT(key) FROM offset_data;").fetchone()
                if self._length <> int(count):
                    con.close()
                    raise ValueError("Corrupt database? %i entries not %i" \
                                     % (int(count), self._length))
                self._format, = con.execute("SELECT value FROM meta_data WHERE key=?;",
                                           ("format",)).fetchone()
                if format and format != self._format:
                    con.close()
                    raise ValueError("Index file says format %s, not %s" \
                                     % (self._format, format))
                self._filenames = [row[0] for row in \
                                  con.execute("SELECT name FROM file_data "
                                              "ORDER BY file_number;").fetchall()]
                if filenames and len(filenames) != len(self._filenames):
                    con.close()
                    raise ValueError("Index file says %i files, not %i" \
                                     % (len(self.filenames) != len(filenames)))
                if filenames and filenames != self._filenames:
                    con.close()
                    raise ValueError("Index file has different filenames")
            except _OperationalError, err:
                con.close()
                raise ValueError("Not a Biopython index database? %s" % err)
            #Now we have the format (from the DB if not given to us),
            try:
                proxy_class = _FormatToRandomAccess[self._format]
            except KeyError:
                con.close()
                raise ValueError("Unsupported format '%s'" % self._format)
        else:
            self._filenames = filenames
            self._format = format
            if not format or not filenames:
                raise ValueError("Filenames to index and format required")
            try:
                proxy_class = _FormatToRandomAccess[format]
            except KeyError:
                raise ValueError("Unsupported format '%s'" % format)
            #Create the index
            con = _sqlite.connect(index_filename)
            self._con = con
            #print "Creating index"
            # Sqlite PRAGMA settings for speed
            con.execute("PRAGMA synchronous='OFF'")
            con.execute("PRAGMA locking_mode=EXCLUSIVE")
            #Don't index the key column until the end (faster)
            #con.execute("CREATE TABLE offset_data (key TEXT PRIMARY KEY, "
            # "offset INTEGER);")
            con.execute("CREATE TABLE meta_data (key TEXT, value TEXT);")
            con.execute("INSERT INTO meta_data (key, value) VALUES (?,?);",
                        ("count", -1))
            con.execute("INSERT INTO meta_data (key, value) VALUES (?,?);",
                        ("format", format))
            #TODO - Record the alphabet?
            #TODO - Record the file size and modified date?
            con.execute("CREATE TABLE file_data (file_number INTEGER, name TEXT);")
            con.execute("CREATE TABLE offset_data (key TEXT, file_number INTEGER, offset INTEGER, length INTEGER);")
            count = 0
            for i, filename in enumerate(filenames):
                con.execute("INSERT INTO file_data (file_number, name) VALUES (?,?);",
                            (i, filename))
                random_access_proxy = proxy_class(filename, format, alphabet)
                if key_function:
                    offset_iter = ((key_function(k),i,o,l) for (k,o,l) in random_access_proxy)
                else:
                    offset_iter = ((k,i,o,l) for (k,o,l) in random_access_proxy)
                while True:
                    batch = list(itertools.islice(offset_iter, 100))
                    if not batch: break
                    #print "Inserting batch of %i offsets, %s ... %s" \
                    # % (len(batch), batch[0][0], batch[-1][0])
                    con.executemany("INSERT INTO offset_data (key,file_number,offset,length) VALUES (?,?,?,?);",
                                    batch)
                    con.commit()
                    count += len(batch)
                if len(random_access_proxies) < max_open:
                    random_access_proxies[i] = random_access_proxy
                else:
                    random_access_proxy._handle.close()
            self._length = count
            #print "About to index %i entries" % count
            try:
                con.execute("CREATE UNIQUE INDEX IF NOT EXISTS "
                            "key_index ON offset_data(key);")
            except _IntegrityError, err:
                self._proxies = random_access_proxies
                self.close()
                con.close()
                raise ValueError("Duplicate key? %s" % err)
            con.execute("PRAGMA locking_mode=NORMAL")
            con.execute("UPDATE meta_data SET value = ? WHERE key = ?;",
                        (count, "count"))
            con.commit()
            #print "Index created"
        self._proxies = random_access_proxies
        self._max_open = max_open
        self._index_filename = index_filename
        self._alphabet = alphabet
        self._key_function = key_function
    
    def __repr__(self):
        return "SeqIO.index_db(%r, filenames=%r, format=%r, alphabet=%r, key_function=%r)" \
               % (self._index_filename, self._filenames, self._format,
                  self._alphabet, self._key_function)

    def __contains__(self, key):
        return bool(self._con.execute("SELECT key FROM offset_data WHERE key=?;",
                                      (key,)).fetchone())

    def __len__(self):
        """How many records are there?"""
        return self._length
        #return self._con.execute("SELECT COUNT(key) FROM offset_data;").fetchone()[0]

    def __iter__(self):
        """Iterate over the keys."""
        for row in self._con.execute("SELECT key FROM offset_data;"):
            yield str(row[0])

    if hasattr(dict, "iteritems"):
        #Python 2, use iteritems but not items etc
        #Just need to override this...
        def keys(self) :
            """Return a list of all the keys (SeqRecord identifiers)."""
            return [str(row[0]) for row in \
                    self._con.execute("SELECT key FROM offset_data;").fetchall()]

    def __getitem__(self, key):
        """x.__getitem__(y) <==> x[y]"""
        #Pass the offset to the proxy
        row = self._con.execute("SELECT file_number, offset FROM offset_data WHERE key=?;",
                                (key,)).fetchone()
        if not row: raise KeyError
        file_number, offset = row
        proxies = self._proxies
        if file_number in proxies:
            record = proxies[file_number].get(offset)
        else:
            if len(proxies) >= self._max_open:
                #Close an old handle...
                proxies.popitem()[1]._handle.close()
            #Open a new handle...
            proxy = _FormatToRandomAccess[self._format]( \
                        self._filenames[file_number],
                        self._format, self._alphabet)
            record = proxy.get(offset)
            proxies[file_number] = proxy
        if self._key_function:
            key2 = self._key_function(record.id)
        else:
            key2 = record.id
        if key != key2:
            raise ValueError("Key did not match (%s vs %s)" % (key, key2))
        return record

    def get(self, k, d=None):
        """D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None."""
        try:
            return self.__getitem__(k)
        except KeyError:
            return d

    def get_raw(self, key):
        """Similar to the get method, but returns the record as a raw string.

        If the key is not found, a KeyError exception is raised.

        Note that on Python 3 a bytes string is returned, not a typical
        unicode string.

        NOTE - This functionality is not supported for every file format.
        """
        #Pass the offset to the proxy
        row = self._con.execute("SELECT file_number, offset, length FROM offset_data WHERE key=?;",
                                (key,)).fetchone()
        if not row: raise KeyError
        file_number, offset, length = row
        proxies = self._proxies
        if file_number in proxies:
            if length:
                #Shortcut if we have the length
                h = proxies[file_number]._handle
                h.seek(offset)
                return h.read(length)
            else:
                return proxies[file_number].get_raw(offset)
        else:
            #This code is duplicated from __getitem__ to avoid a function call
            if len(proxies) >= self._max_open:
                #Close an old handle...
                proxies.popitem()[1]._handle.close()
            #Open a new handle...
            proxy = _FormatToRandomAccess[self._format]( \
                        self._filenames[file_number],
                        self._format, self._alphabet)
            proxies[file_number] = proxy
            if length:
                #Shortcut if we have the length
                h = proxy._handle
                h.seek(offset)
                return h.read(length)
            else:
                return proxy.get_raw(offset)

    def close(self):
        """Close any open file handles."""
        proxies = self._proxies
        while proxies:
            proxies.popitem()[1]._handle.close()
        

##############################################################################

class SeqFileRandomAccess(object):
    def __init__(self, filename, format, alphabet):
        self._handle = open(filename, "rb")
        self._alphabet = alphabet
        self._format = format
        #Load the parser class/function once an avoid the dict lookup in each
        #__getitem__ call:
        i = SeqIO._FormatToIterator[format]
        #The following alphabet code is a bit nasty... duplicates logic in
        #Bio.SeqIO.parse()
        if alphabet is None:
            def _parse(handle):
                """Dynamically generated parser function (PRIVATE)."""
                return i(handle).next()
        else:
            #TODO - Detect alphabet support ONCE at __init__
            def _parse(handle):
                """Dynamically generated parser function (PRIVATE)."""
                try:
                    return i(handle, alphabet=alphabet).next()
                except TypeError:
                    return SeqIO._force_alphabet(i(handle),
                                                 alphabet).next()
        self._parse = _parse

    def __iter__(self):
        """Returns (id,offset) tuples."""
        raise NotImplementedError("Subclass should implement this")

    def get(self, offset):
        """Returns SeqRecord."""
        #Should be overriden for binary file formats etc:
        return self._parse(StringIO(_bytes_to_string(self.get_raw(offset))))

    def get_raw(self, offset):
        """Returns bytes string (if implemented for this file format)."""
        #Should be done by each sub-class (if possible)
        raise NotImplementedError("Not available for this file format.")




####################
# Special indexers #
####################

# Anything where the records cannot be read simply by parsing from
# the record start. For example, anything requiring information from
# a file header - e.g. SFF files where we would need to know the
# number of flows.

class SffRandomAccess(SeqFileRandomAccess):
    """Random access to a Standard Flowgram Format (SFF) file."""
    def __init__(self, filename, format, alphabet):
        SeqFileRandomAccess.__init__(self, filename, format, alphabet)
        header_length, index_offset, index_length, number_of_reads, \
        self._flows_per_read, self._flow_chars, self._key_sequence \
            = SeqIO.SffIO._sff_file_header(self._handle)

    def __iter__(self):
        """Load any index block in the file, or build it the slow way (PRIVATE)."""
        if self._alphabet is None:
            self._alphabet = Alphabet.generic_dna
        handle = self._handle
        handle.seek(0)
        #Alread did this in __init__ but need handle in right place
        header_length, index_offset, index_length, number_of_reads, \
        self._flows_per_read, self._flow_chars, self._key_sequence \
            = SeqIO.SffIO._sff_file_header(handle)
        if index_offset and index_length:
            #There is an index provided, try this the fast way:
            count = 0
            try :
                for name, offset in SeqIO.SffIO._sff_read_roche_index(handle) :
                    yield name, offset, 0
                    count += 1
                assert count == number_of_reads, \
                       "Indexed %i records, expected %i" \
                       % (count, number_of_reads)
                return
            except ValueError, err :
                import warnings
                warnings.warn("Could not parse the SFF index: %s" % err)
                assert count==0, "Partially populated index"
                handle.seek(0)
        #We used to give a warning in this case, but Ion Torrent's
        #SFF files don't have an index so that would be annoying.
        #Fall back on the slow way!
        count = 0
        for name, offset in SeqIO.SffIO._sff_do_slow_index(handle) :
            yield name, offset, 0
            count += 1
        assert count == number_of_reads, \
               "Indexed %i records, expected %i" % (count, number_of_reads)

    def get(self, offset) :
        handle = self._handle
        handle.seek(offset)
        return SeqIO.SffIO._sff_read_seq_record(handle,
                                                self._flows_per_read,
                                                self._flow_chars,
                                                self._key_sequence,
                                                self._alphabet)

    def get_raw(self, offset):
        handle = self._handle
        handle.seek(offset)
        return SeqIO.SffIO._sff_read_raw_record(handle, self._flows_per_read)


class SffTrimedRandomAccess(SffRandomAccess) :
    def get(self, offset) :
        handle = self._handle
        handle.seek(offset)
        return SeqIO.SffIO._sff_read_seq_record(handle,
                                                self._flows_per_read,
                                                self._flow_chars,
                                                self._key_sequence,
                                                self._alphabet,
                                                trim=True)


###################
# Simple indexers #
###################

class SequentialSeqFileRandomAccess(SeqFileRandomAccess):
    def __init__(self, filename, format, alphabet):
        SeqFileRandomAccess.__init__(self, filename, format, alphabet)
        marker = {"ace" : "CO ",
                  "embl" : "ID ",
                  "fasta" : ">",
                  "genbank" : "LOCUS ",
                  "gb": "LOCUS ",
                  "imgt" : "ID ",
                  "phd" : "BEGIN_SEQUENCE",
                  "pir" : ">..;",
                  "qual": ">",
                  "qual": ">",
                  "swiss" : "ID ",
                  "uniprot-xml" : "<entry ",
                   }[format]
        self._marker = marker
        self._marker_re = re.compile(_as_bytes("^%s" % marker))
        
    def __iter__(self):
        """Returns (id,offset) tuples."""
        marker_offset = len(self._marker)
        marker_re = self._marker_re
        handle = self._handle
        handle.seek(0)
        #Skip and header before first record
        while True:
            start_offset = handle.tell()
            line = handle.readline()
            if marker_re.match(line) or not line:
                break
        #Should now be at the start of a record, or end of the file
        while marker_re.match(line):
            #Here we can assume the record.id is the first word after the
            #marker. This is generally fine... but not for GenBank, EMBL, Swiss
            id = line[marker_offset:].strip().split(None, 1)[0]
            while True:
                line = handle.readline()
                if marker_re.match(line) or not line:
                    end_offset = handle.tell() - len(line)
                    yield _bytes_to_string(id), start_offset, end_offset - start_offset
                    start_offset = end_offset
                    break
        assert not line, repr(line)

    def get_raw(self, offset):
        """Similar to the get method, but returns the record as a raw string."""
        #For non-trivial file formats this must be over-ridden in the subclass
        handle = self._handle
        marker_re = self._marker_re
        handle.seek(offset)
        lines = [handle.readline()]
        while True:
            line = handle.readline()
            if marker_re.match(line) or not line:
                #End of file, or start of next record => end of this record
                break
            lines.append(line)
        return _as_bytes("").join(lines)


#######################################
# Fiddly indexers: GenBank, EMBL, ... #
#######################################

class GenBankRandomAccess(SequentialSeqFileRandomAccess):
    """Indexed dictionary like access to a GenBank file."""
    def __iter__(self):
        handle = self._handle
        handle.seek(0)
        marker_re = self._marker_re
        dot_char = _as_bytes(".")
        accession_marker = _as_bytes("ACCESSION ")
        version_marker = _as_bytes("VERSION ")
        #Skip and header before first record
        while True:
            start_offset = handle.tell()
            line = handle.readline()
            if marker_re.match(line) or not line:
                break
        #Should now be at the start of a record, or end of the file
        while marker_re.match(line):
            #We cannot assume the record.id is the first word after LOCUS,
            #normally the first entry on the VERSION or ACCESSION line is used.
            key = None
            while True:
                line = handle.readline()
                if marker_re.match(line) or not line:
                    if not key:
                        raise ValueError("Did not find ACCESSION/VERSION lines")
                    end_offset = handle.tell() - len(line)
                    yield _bytes_to_string(key), start_offset, end_offset - start_offset
                    start_offset = end_offset
                    break
                elif line.startswith(accession_marker):
                    key = line.rstrip().split()[1]
                elif line.startswith(version_marker):
                    version_id = line.rstrip().split()[1]
                    if version_id.count(dot_char)==1 and version_id.split(dot_char)[1].isdigit():
                        #This should mimic the GenBank parser...
                        key = version_id
        assert not line, repr(line)


class EmblRandomAccess(SequentialSeqFileRandomAccess):
    """Indexed dictionary like access to an EMBL file."""
    def __iter__(self):
        handle = self._handle
        handle.seek(0)
        marker_re = self._marker_re
        semi_char = _as_bytes(";")
        dot_char = _as_bytes(".")
        sv_marker = _as_bytes("SV ")
        #Skip any header before first record
        while True:
            start_offset = handle.tell()
            line = handle.readline()
            if marker_re.match(line) or not line:
                break
        #Should now be at the start of a record, or end of the file
        while marker_re.match(line):
            #We cannot assume the record.id is the first word after ID,
            #normally the SV line is used.
            if line[2:].count(semi_char) == 6:
                #Looks like the semi colon separated style introduced in 2006
                parts = line[3:].rstrip().split(semi_char)
                if parts[1].strip().startswith(sv_marker):
                    #The SV bit gives the version
                    key = parts[0].strip() + dot_char + parts[1].strip().split()[1]
                else:
                    key = parts[0].strip()
            elif line[2:].count(semi_char) == 3:
                #Looks like the pre 2006 style, take first word only
                key = line[3:].strip().split(None,1)[0]
            else:
                raise ValueError('Did not recognise the ID line layout:\n' + line)
            while True:
                line = handle.readline()
                if marker_re.match(line) or not line:
                    end_offset = handle.tell() - len(line)
                    yield _bytes_to_string(key), start_offset, end_offset - start_offset
                    start_offset = end_offset
                    break
                elif line.startswith(sv_marker):
                    key = line.rstrip().split()[1]
        assert not line, repr(line)


class SwissRandomAccess(SequentialSeqFileRandomAccess):
    """Random access to a SwissProt file."""
    def __iter__(self):
        handle = self._handle
        handle.seek(0)
        marker_re = self._marker_re
        semi_char = _as_bytes(";")
        #Skip any header before first record
        while True:
            start_offset = handle.tell()
            line = handle.readline()
            if marker_re.match(line) or not line:
                break
        #Should now be at the start of a record, or end of the file
        while marker_re.match(line):
            #We cannot assume the record.id is the first word after ID,
            #normally the following AC line is used.
            line = handle.readline()
            assert line.startswith(_as_bytes("AC "))
            key = line[3:].strip().split(semi_char)[0].strip()
            while True:
                line = handle.readline()
                if marker_re.match(line) or not line:
                    end_offset = handle.tell() - len(line)
                    yield _bytes_to_string(key), start_offset, end_offset - start_offset
                    start_offset = end_offset
                    break
        assert not line, repr(line)


class UniprotRandomAccess(SequentialSeqFileRandomAccess):
    """Random access to a UniProt XML file."""
    def __iter__(self):
        handle = self._handle
        handle.seek(0)
        marker_re = self._marker_re
        start_acc_marker = _as_bytes("<accession>")
        end_acc_marker = _as_bytes("</accession>")
        end_entry_marker = _as_bytes("</entry>")
        #Skip any header before first record
        while True:
            start_offset = handle.tell()
            line = handle.readline()
            if marker_re.match(line) or not line:
                break
        #Should now be at the start of a record, or end of the file
        while marker_re.match(line):
            #We expect the next line to be <accession>xxx</accession>
            #(possibly with leading spaces)
            #but allow it to be later on within the <entry>
            key = None
            done = False
            while True:
                line = handle.readline()
                if key is None and start_acc_marker in line:
                    assert end_acc_marker in line, line
                    key = line[line.find(start_acc_marker)+11:].split(_as_bytes("<"))[0]
                elif end_entry_marker in line:
                    end_offset = handle.tell() - len(line) \
                               + line.find(end_entry_marker) + 8
                    break
                elif marker_re.match(line) or not line:
                    #Start of next record or end of file
                    raise ValueError("Didn't find end of record")
            if not key:
                raise ValueError("Did not find <accession> line in bytes %i to %i" \
                                 % (start_offset, end_offset))
            yield _bytes_to_string(key), start_offset, end_offset - start_offset
            #Find start of next record
            while not marker_re.match(line) and line:
                start_offset = handle.tell()
                line = handle.readline()
        assert not line, repr(line)
    
    def get_raw(self, offset):
        """Similar to the get method, but returns the record as a raw string."""
        handle = self._handle
        marker_re = self._marker_re
        end_entry_marker = _as_bytes("</entry>")
        handle.seek(offset)
        data = handle.readline()
        while True:
            line = handle.readline()
            i = line.find(end_entry_marker)
            if i != -1:
                data += line[:i+8]
                break
            if marker_re.match(line) or not line:
                #End of file, or start of next record
                raise ValueError("Didn't find end of record")
            data += line
        return data

    def get(self, offset) :
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
        """ % _bytes_to_string(self.get_raw(offset))
        #TODO - For consistency, this function should not accept a string:
        return SeqIO.UniprotIO.UniprotIterator(data).next()


class IntelliGeneticsRandomAccess(SeqFileRandomAccess):
    """Random access to a IntelliGenetics file."""
    def __init__(self, filename, format, alphabet):
        SeqFileRandomAccess.__init__(self, filename, format, alphabet)
        self._marker_re = re.compile(_as_bytes("^;"))

    def __iter__(self):
        handle = self._handle
        handle.seek(0)
        marker_re = self._marker_re
        semi_char = _as_bytes(";")
        while True:
            offset = handle.tell()
            line = handle.readline()
            if marker_re.match(line):
                #Now look for the first line which doesn't start ";"
                while True:
                    line = handle.readline()
                    if line[0:1] != semi_char and line.strip():
                        key = line.split()[0]
                        yield _bytes_to_string(key), offset, 0
                        break
                    if not line:
                        raise ValueError("Premature end of file?")
            elif not line:
                #End of file
                break


    def get_raw(self, offset):
        handle = self._handle
        handle.seek(offset)
        marker_re = self._marker_re
        lines = []
        line = handle.readline()
        semi_char = _as_bytes(";")
        while line.startswith(semi_char):
            lines.append(line)
            line = handle.readline()
        while line and not line.startswith(semi_char):
            lines.append(line)
            line = handle.readline()
        return _as_bytes("").join(lines)

class TabRandomAccess(SeqFileRandomAccess):
    """Random access to a simple tabbed file."""
    def __iter__(self):
        handle = self._handle
        handle.seek(0)
        start_offset = handle.tell()
        tab_char = _as_bytes("\t")
        while True:
            line = handle.readline()
            if not line : break #End of file
            try:
                key = line.split(tab_char)[0]
            except ValueError, err:
                if not line.strip():
                    #Ignore blank lines
                    start_offset = handle.tell()
                    continue
                else:
                    raise err
            else:
                end_offset = handle.tell()
                yield _bytes_to_string(key), start_offset, end_offset - start_offset
                start_offset = end_offset

    def get_raw(self, offset):
        """Like the get method, but returns the record as a raw string."""
        handle = self._handle
        handle.seek(offset)
        return handle.readline()


##########################
# Now the FASTQ indexers #
##########################
         
class FastqRandomAccess(SeqFileRandomAccess):
    """Random access to a FASTQ file (any supported variant).
    
    With FASTQ the records all start with a "@" line, but so can quality lines.
    Note this will cope with line-wrapped FASTQ files.
    """
    def __iter__(self):
        handle = self._handle
        handle.seek(0)
        id = None
        start_offset = handle.tell()
        line = handle.readline()
        if not line:
            #Empty file!
            return
        at_char = _as_bytes("@")
        plus_char = _as_bytes("+")
        if line[0:1] != at_char:
            raise ValueError("Problem with FASTQ @ line:\n%s" % repr(line))
        while line:
            #assert line[0]=="@"
            #This record seems OK (so far)
            id = line[1:].rstrip().split(None, 1)[0]
            #Find the seq line(s)
            seq_len = 0
            while line:
                line = handle.readline()
                if line.startswith(plus_char) : break
                seq_len += len(line.strip())
            if not line:
                raise ValueError("Premature end of file in seq section")
            #assert line[0]=="+"
            #Find the qual line(s)
            qual_len = 0
            while line:
                if seq_len == qual_len:
                    #Should be end of record...
                    line = handle.readline()
                    if line and line[0:1] != at_char:
                        ValueError("Problem with line %s" % repr(line))
                    break
                else:
                    line = handle.readline()
                    qual_len += len(line.strip())
            if seq_len != qual_len:
                raise ValueError("Problem with quality section")
            end_offset = handle.tell() - len(line)
            yield _bytes_to_string(id), start_offset, end_offset - start_offset
            start_offset = end_offset
        #print "EOF"

    def get_raw(self, offset):
        """Similar to the get method, but returns the record as a raw string."""
        #TODO - Refactor this and the __init__ method to reduce code duplication?
        handle = self._handle
        handle.seek(offset)
        line = handle.readline()
        data = line
        at_char = _as_bytes("@")
        plus_char = _as_bytes("+")
        if line[0:1] != at_char:
            raise ValueError("Problem with FASTQ @ line:\n%s" % repr(line))
        identifier = line[1:].rstrip().split(None, 1)[0]
        #Find the seq line(s)
        seq_len = 0
        while line:
            line = handle.readline()
            data += line
            if line.startswith(plus_char) : break
            seq_len += len(line.strip())
        if not line:
            raise ValueError("Premature end of file in seq section")
        assert line[0:1] == plus_char
        #Find the qual line(s)
        qual_len = 0
        while line:
            if seq_len == qual_len:
                #Should be end of record...
                pos = handle.tell()
                line = handle.readline()
                if line and line[0:1] != at_char:
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

_FormatToRandomAccess = {"ace" : SequentialSeqFileRandomAccess,
                        "embl" : EmblRandomAccess,
                        "fasta" : SequentialSeqFileRandomAccess,
                        "fastq" : FastqRandomAccess, #Class handles all three variants
                        "fastq-sanger" : FastqRandomAccess, #alias of the above
                        "fastq-solexa" : FastqRandomAccess,
                        "fastq-illumina" : FastqRandomAccess,
                        "genbank" : GenBankRandomAccess,
                        "gb" : GenBankRandomAccess, #alias of the above
                        "ig" : IntelliGeneticsRandomAccess,
                        "imgt" : EmblRandomAccess,
                        "phd" : SequentialSeqFileRandomAccess,
                        "pir" : SequentialSeqFileRandomAccess,
                        "sff" : SffRandomAccess,
                        "sff-trim" : SffTrimedRandomAccess,
                        "swiss" : SwissRandomAccess,
                        "tab" : TabRandomAccess,
                        "qual" : SequentialSeqFileRandomAccess,
                        "uniprot-xml" : UniprotRandomAccess, 
                        }

