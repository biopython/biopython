# Copyright 2009-2011 by Peter Cock.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

#TODO: factor out this module with SeqIO's _index to stay DRY

"""Custom indexing for Bio.SearchIO objects (PRIVATE)."""

import itertools
import os
from StringIO import StringIO
try:
    from sqlite3 import dbapi2 as sqlite
    from sqlite3 import IntegrityError, OperationalError
except ImportError:
    # apparently jython2.5 may not support sqlite
    sqlite = None
try:
    from collections import UserDict as _dict_base
except ImportError:
    from UserDict import DictMixin as _dict_base

from Bio._py3k import _bytes_to_string
from Bio.SearchIO import _INDEXER_MAP
from Bio.SearchIO._utils import get_processor


__all__ = ['SearchIndexer']


class SearchIndexer(object):

    """Iterator returning file positions of results in a search output file."""

    def __init__(self, filename, **kwargs):
        self._handle = open(filename, 'rb')
        self._handle.seek(0) 
        self._kwargs = kwargs

    def _parse(self, handle):
        return iter(self._parser(handle, **self._kwargs)).next()

    def get(self, offset):
        return self._parse(StringIO(_bytes_to_string(self.get_raw(offset))))


class _IndexedSearch(_dict_base):

    """Dictionary-like object for implementing Search indexing."""

    def __init__(self, filename, format, key_function=None, **kwargs):
        """Initializes _IndexedSearch instance.

        filename -- The source filename as string.
        format -- Lower case string denoting one of the supported formats.
        key_function -- Optional callbak function which when given a Result
                        should return a unique key for the dictionary.

        """
        self._filename = filename
        self._format = format
        self._key_function = key_function

        indexer_class = get_processor(format, _INDEXER_MAP)
        indexed_obj = indexer_class(filename, **kwargs)
        self._indexer = indexed_obj

        # default key function is lambda rec: rec.id
        if key_function:
            offset_iter = ((key_function(key), offset, length) for \
                    (key, offset, length) in indexed_obj)
        else:
            offset_iter = indexed_obj

        index = {}
        for key, offset, length in offset_iter:
            if key in index:
                self._indexer._handle.close()
                raise ValueError("Duplicate key %r" % key)
            else:
                index[key] = offset

        self._index = index

    def __iter__(self):
        return iter(self._index)

    def __contains__(self, key):
        return key in self._index

    def __len__(self):
        return len(self._index)

    def __repr__(self):
        return "SearchIO.index(%r, %r, key_function=%r)" % \
                (self._filename, self._format, self._key_function)

    def __str__(self):
        if self:
            return "{%r: QueryResult(...), ...}" % repr(self.keys()[0])
        else:
            return "{}"

    # dynamic method assignments to deal with different dict behavior in
    # python 2 and python 3
    # case 1: python2
    if hasattr(dict, 'iteritems'):
        def values(self):
            raise NotImplementedError("Due to memory concerns, when indexing "
                    "a search output file you cannot access all the records "
                    "at once.")

        def items(self):
            raise NotImplementedError("Due to memory concerns, when indexing "
                    "a search output file you cannot access all the records "
                    "at once.")

        def keys(self):
            return self._index.keys()

        def itervalues(self):
            for key in self.__iter__():
                yield self.__getitem__(key)

        def iteritems(self):
            for key in self.__iter__():
                yield key, self.__getitem__(key)

        def iterkeys(self):
            return self.__iter__()

    # case 2: python3, where dict.{items,keys,values} returns a generator-like
    # object
    else:
        def values(self):
            for key in self.__iter__():
                yield self.__getitem__(key)

        def items(self):
            for key in self.__iter__():
                yield key, self.__getitem__(key)

        def keys(self):
            return self.__iter__()

    def __getitem__(self, key):
        result = self._indexer.get(self._index[key])
        if self._key_function:
            key2 = self._key_function(result.id)
        else:
            key2 = result.id

        if key != key2:
            raise ValueError("Key did not match (%r vs %r)" % (key, key2))

        return result

    def get(self, key, default=None):
        try:
            return self.__getitem__(key)
        except KeyError:
            return default

    def get_raw(self, key):
        """Similar to get, but returns the raw string of the Result object.

        Note that any text or information in the search output file located
        prior to the first result is ignored.

        """
        return self._indexer.get_raw(self._index[key])

    def __delitem__(self, key, value):
        raise NotImplementedError("An indexed search output file is read-only.")

    def __setitem__(self, key, value):
        raise NotImplementedError("An indexed search output file is read-only.")

    def clear(self, *args, **kwargs):
        raise NotImplementedError("An indexed search output file is read-only.")

    def pop(self, key, default=None):
        raise NotImplementedError("An indexed search output file is read-only.")

    def popitem(self):
        raise NotImplementedError("An indexed search output file is read-only.")

    def update(self, *args, **kwargs):
        raise NotImplementedError("An indexed search output file is read-only.")

    def copy(self, *args, **kwargs):
        raise NotImplementedError("An indexed search output file does not "
                "support this.")

    def fromkeys(self, *args, **kwargs):
        raise NotImplementedError("An indexed search output file does not "
                "support this.")


class _DbIndexedSearch(_IndexedSearch):

    """Dictionary-like object for implementing storable Search indexing."""

    def __init__(self, index_filename, filenames, format, key_function,
            max_open=10, overwrite=False, **kwargs):
        """Initializes a _DbIndexedSearch instance.

        index_filename -- The SQLite filename.
        filenames -- List of strings specifying file(s) to be indexed, or when
                     indexing a single file this can be given as a string.
                     (optional if reloading an existing index, but must match)
        format -- Lower case string denoting one of the supported formats.
                  (optional if reloading an existing index, but must match)
        indexer -- Format-specific Indexer class.
        key_function - Optional callback function which when given a
                       Result identifier string should return a unique
                       key for the dictionary.
        max_open -- Integer of maximum open file objects allowed.
        overwrite -- Boolean, whether to overwrite existing index database
                     (if it exists) or not.

        """
        # COMPAT: for Jython, which may not have sqlite3 baked in
        if not sqlite:
            from Bio import MissingPythonDependencyError
            raise MissingPythonDependencyError("Requires sqlite3, which is "
                                               "included Python 2.5+")
        indexer_proxies = {}
        indexer_class = get_processor(format, _INDEXER_MAP)
        self._indexer_class = indexer_class
        self._kwargs = kwargs

        # remove index_filename if overwrite is True and the file exists
        if overwrite and os.path.isfile(index_filename):
            os.remove(index_filename)

        if os.path.isfile(index_filename):
            con = sqlite.connect(index_filename)
            self._con = con
            try:
                # get the # of result offsets stored in the database
                count, = con.execute("SELECT value FROM meta_data WHERE "
                        "key=?;", ('count',)).fetchone()
                self._length = int(count)
                if self._length == -1:
                    con.close()
                    raise ValueError("Unfinished/partial database")

                # count the # of result offsets stored in the database
                # for cross-checking
                count, = con.execute("SELECT COUNT(key) FROM "
                        "offset_data;").fetchone()
                if self._length != int(count):
                    con.close()
                    raise ValueError("Corrupt database? %i entries not %i" %
                            (int(count), self._length))

                # check if the database format is the same as the given format
                self._format, = con.execute("SELECT value FROM meta_data "
                        "WHERE key=?;", ('format',)).fetchone()
                if format and format != self._format:
                    con.close()
                    raise ValueError("Incorrect format specified: '%s', "
                            "expected: '%s'" % (format, self._format))

                # check filenames # and names
                self._filenames = [row[0] for row in \
                        con.execute("SELECT name FROM file_data ORDER BY "
                            "file_number;").fetchall()]
                if filenames and len(filenames) != len(self._filenames):
                    con.close()
                    raise ValueError("Index file says %i files, not %i" %
                            (len(self._filenames), len(filenames)))
                if filenames and filenames != self._filenames:
                    con.close()
                    raise ValueError("Index file has different filenames")

            except OperationalError, err:
                con.close()
                raise ValueError("Not a Biopython index database? %s" % err)
        else:
            self._filenames = filenames
            self._format = format

            # create the index db file
            con = sqlite.connect(index_filename)
            self._con = con

            # speed optimization
            con.execute("PRAGMA synchronous=OFF")
            con.execute("PRAGMA locking_mode=EXCLUSIVE")

            # create the tables
            # meta_data: for storing # of results and file format
            # file_data: for storing source filenames
            # offset_data: for storing index, offset pair
            con.execute("CREATE TABLE meta_data (key TEXT, value TEXT);")
            con.execute("INSERT INTO meta_data (key, value) VALUES (?,?);",
                    ('count', -1))
            con.execute("INSERT INTO meta_data (key, value) VALUES (?,?);",
                    ('format', format))
            con.execute("CREATE TABLE file_data (file_number INTEGER, "
                    "name TEXT);")
            con.execute("CREATE TABLE offset_data (key TEXT, file_number "
                    "INTEGER, offset INTEGER, length INTEGER);")
            
            # and fill them up
            count = 0
            for idx, filename in enumerate(filenames):
                # fill the file_data
                con.execute("INSERT INTO file_data(file_number, name) VALUES "
                        "(?,?);", (idx, filename))
                indexed_obj = indexer_class(filename, **self._kwargs)

                if key_function:
                    offset_iter = ((key_function(key), idx, offset, length) for
                            (key, offset, length) in indexed_obj)
                else:
                    offset_iter = ((key, idx, offset, length) for (key, offset,
                            length) in indexed_obj)

                # and the results in a file into offset_data
                while True:
                    batch = list(itertools.islice(offset_iter, 100))
                    if not batch:
                        break
                    con.executemany("INSERT INTO offset_data "
                            "(key,file_number,offset,length) VALUES "
                            "(?,?,?,?);", batch)
                    con.commit()
                    count += len(batch)

                # check if we can still open more file objects
                if len(indexer_proxies) < max_open:
                    indexer_proxies[idx] = indexed_obj
                else:
                    indexed_obj._handle.close()

            self._length = count

            try:
                con.execute("CREATE UNIQUE INDEX IF NOT EXISTS "
                        "key_index ON offset_data(key);")
            except IntegrityError, err:
                self._proxies = indexer_proxies
                self.close()
                con.close()
                raise ValueError("Duplicate key? %s" % err)

            con.execute("PRAGMA locking_mode=NORMAL")
            con.execute("UPDATE meta_data SET value = ? WHERE key = ?;",
                    (count, "count"))

            con.commit()

        self._proxies = indexer_proxies
        self._max_open = max_open
        self._index_filename = index_filename
        self._key_function = key_function

    def __repr__(self):
        return "SearchIO.index_db(%r, %r, sources=%r, key_function=%r)" % \
                (self._index_filename, self._format, self._filenames,
                self._key_function)

    def __contains__(self, key):
        return bool(self._con.execute("SELECT key FROM offset_data WHERE "
            "key=?;", (key,)).fetchone())

    def __len__(self):
        return self._length

    def __iter__(self):
        for row in self._con.execute("SELECT key FROM offset_data;"):
            yield str(row[0])

    # handle python2
    if hasattr(dict, 'iteritems'):
        def keys(self):
            return [str(row[0]) for row in
                self._con.execute("SELECT key FROM offset_data;").fetchall()]

    def __getitem__(self, key):
        row = self._con.execute("SELECT file_number, offset FROM offset_data "
                "WHERE key=?;", (key,)).fetchone()
        if not row:
            raise KeyError
        file_number, offset = row
        proxies = self._proxies

        if file_number in proxies:
            result = proxies[file_number].get(offset)
        else:
            if len(proxies) >= self._max_open:
                proxies.popitem()[1]._handle.close()
            # open a new handle
            proxy = self._indexer_class(self._filenames[file_number],
                    **self._kwargs)
            result = proxy.get(offset)
            proxies[file_number] = proxy

        if self._key_function:
            key2 = self._key_function(result.id)
        else:
            key2 = result.id
        if key != key2:
            raise ValueError("Key does not match (%r vs %r)" % (key, key2))

        return result

    def get_raw(self, key):
        """Similar to get, but returns the raw string of the Result object.

        Note that any text or information in the search output file located
        prior to the first result is ignored.

        """
        row = self._con.execute("SELECT file_number, offset, length FROM "
                "offset_data WHERE key=?;", (key,)).fetchone()

        if not row:
            raise KeyError("Key %r does not point to any result in "
                    "index." % key)

        file_number, offset, length = row
        proxies = self._proxies
        if file_number in proxies:
            if length:
                handle = proxies[file_number]._handle
                handle.seek(offset)
                return handle.read(length)
            else:
                return proxies[file_number].get_raw(offset)
        else:
            if len(proxies) >= self._max_open:
                proxies.popitem()[1]._handle.close()
            # open a new handle
            proxy = self._indexer_class(self._filenames[file_number],
                    **self._kwargs)
            proxies[file_number] = proxy
            if length:
                handle = proxy._handle
                handle.seek(offset)
                return handle.read(length)
            else:
                return proxy.get_raw(offset)

    def close(self):
        """Close any open file handles."""
        proxies = self._proxies
        while proxies:
            proxies.popitem()[1]._handle.close()
