# Copyright 2009-2011 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Lazy parsing and feature indexing of sequence files (PRIVATE).

You are not expected to access this module, or any of its code,
directly. This is all handled internally by the following functions
when accessed with the lazy=True kwarg.

    Bio.SeqIO.read(..., lazy=True)
    Bio.SeqIO.parse(.... lazy=True)
    Bio.SeqIO.index(..., lazy=True)

The lazy loading parsers will read an entire record efficiently storing
and parsing only the mimimum amount of information to provide fast
lookup. Lazy records returned from parse(), index_db(), and read() will
act as a proxy to a regular SeqRecord class. Internally they will be
quite different but externally accessing any property shared by the
typical SeqRecord will return the correct value. Because the lazy record
inherits from SeqRecord, all SeqRecord methods are available.

The lazy loading strategy is, for very large sequences, faster to
initialize than SeqIO.parse() and SeqIO.index_db() but also more
memory efficient. The when fully accessing SeqRecord attributes
the lazy loading parser reaches performance parity with both standard
parsers. Partial record access, such as only using the seq attribute
or slicing the record beore use, gives significant performance
advantages to using the lazy parser.

The lazy loader will partially parse a sequence file noting the
sequence span and file position of sequence annoations. These
annotations will be stored for efficient retrieval from the file when
requested. Sequence data will also be efficiently pulled from
structured regions and file so that long sequences aren't necessarily
loaded into memory fully.

To take full advantage of the lazy loading strategy, operate slices on
the lazy record whenever possible instead of the record attributes. See
the example below using a lazy loading record 'bigrecord'

    #getting sequence information
    bigrecord.seq[2000000:3000000]  # bad: this parses the whole seq
    bigrecord[2000000:3000000].seq  # good: only 10M bases are loaded

Because all parsing is on demand, improperly formatted records may not
trigger an exception until the problem region is requested.
"""

from copy import copy
import re
from math import floor, ceil, log
from os.path import isfile, basename
try:
    from sqlite3 import dbapi2 as _sqlite
except ImportError:
    #Even without sqlite we still want to offer in-memory indexing.
    _sqlite = None

from Bio._py3k import _is_int_or_long, _bytes_to_string, _as_string
from ..SeqRecord import SeqRecord

def _get_db_connection(dbfile):
    """do some basic checking and return a connection to the db"""
    if not _sqlite:
        # Hack for Jython (or if Python is compiled without it)
        from Bio import MissingPythonDependencyError
        raise MissingPythonDependencyError("Requires sqlite3, which is "
                                           "included Python 2.5+")
    if not isfile(dbfile):
        raise ValueError("A valid database or None must be provided")
    #case: connection does not exist yet
    con = _sqlite.connect(dbfile)
    return con

class SeqProxyIndexManager(object):
    """Class to manage the index and any sqlite db-IO required

    From the perspective of SeqRecordProxyBase derived classes,
    this index manager will provide the following public attriubes
    and method.

    Attributes:
      record (dict): stores the main record index
    """

    def __init__(self, filefmt, indexdb=None, recordkey=None, handlename=None):
        """Initialize SeqProxyIndexManager

        The SeqProxyIndexManager may be initialized with or without an
        index database file, this init must adjust behavior of the class
        to account for both cases.

        Args:
          filefmt (str): The format of files to be indexed.
          indexdb (str): A string with the full path to the SQLite
            index bein used
          recordkey (str or int): When fetching from the database, the
            `recordkey` is used to grab a specific record. int type keys
            will fetch a record by the order in which it was added.
          handlename (str): The basename of the file a record is
            associated with. This is used for both storage and lookup.
        """
        #recordkey should never be passed in the absence of a valid indexdb
        if recordkey is not None and indexdb is None:
            raise ValueError("SeqProxyIndexManager requires an indexdb"
                             " to fetch by key.")
        elif indexdb is not None and handlename is None:
            raise ValueError("SeqProxyIndexManager requires handlename"
                             " to store records.")
        #public attributes
        self.record = {}
        #private attributes
        self._format = filefmt
        self._recordkey = recordkey
        self._features = None
        self._handlename = handlename
        self._indexdb = indexdb
        self._handle_id = None
        self._index_id = None
        self.index_exists = False

        if self._indexdb is not None and recordkey is not None:
            self.record.update(self._load_record_index())
            self.index_exists = True

    def _load_record_index(self):
        """Returns the index dictionary from memory or from database

        When the database (_indexdb) is None:
        This function will preferentially supply an index from memory.
        Since the _indexdb is None, it will provide None for the _index
        indicating to the lazy loading proxy that an index must be
        created.

        When the database (_indexdb) is set:
        This function will still preferentially supply an index from
        memory. If the database is valid and empty it will provide None
        for the _index indicating to the lazy loading proxy that an
        index must be created. If an indexdb is present and valid,
        it will be queried to produce the _index dictionary which
        will be loaded into memory.
        """

        #this will connect to the db
        con = _get_db_connection(self._indexdb)
        cursor = con.cursor()

        #case: indexdb is provided and has records
        if _is_int_or_long(self._recordkey):
            recordindex = cursor.execute("SELECT main_index.* "
               "FROM main_index "
               "INNER JOIN indexed_files "
               "ON main_index.fileid = indexed_files.fileid "
               "WHERE indexed_files.filename=? "
               "ORDER BY main_index.recordoffsetstart "
               "LIMIT 2 OFFSET ?;",
               (self._handlename, self._recordkey))
        elif isinstance(self._recordkey, str):
            recordindex = cursor.execute("SELECT main_index.* "
               "FROM main_index "
               "INNER JOIN indexed_files "
               "ON main_index.fileid = indexed_files.fileid "
               "WHERE indexed_files.filename=? AND "
               "main_index.id=?;",
               (self._handlename, self._recordkey))
        #description is always set, so this is safe to do first
        record_keys = [key[0] for key in recordindex.description]
        record_index = recordindex.fetchone()
        con.commit()
        con.close()
        if record_index is None:
            raise KeyError("key {0} not found".format(self._recordkey))
        record_indexdict = dict((record_keys[i], record_index[i]) for \
                                    i in range(len(record_keys)))
        # save some usefulreturn the index
        self._index_id = record_indexdict.pop("indexid")
        self._handle_id = record_indexdict.pop("fileid")
        #the return
        return record_indexdict


    def set_record_index(self, indexdict):
        """Save a new index dict to the database

        Given a new indexdict, the database will be loaded with the
        record unless the database already contains the given record.
        Should a conflict be found, a ValueError will be raised.

        When no database exists, using this just sets the dict in
        memory for later retrieval.

        Args:
          indexdict (dict): A dictionary corresponding to an indexed
            sequence record.
        """
        if not isinstance(indexdict, dict):
            raise ValueError("New index must be a dictionary")
        if self._indexdb is None:
            self.record.update(indexdict)
            index_exists = True
            return

        con = _get_db_connection(self._indexdb)
        #case: indexdb provided but it is empty
        if not self._is_valid_db_with_tables(con):
            self._create_tables(con, indexdict)
            con.commit()
        cursor = con.cursor()
        fileid = self._get_handle_id(con, write=True)

        #case: indexdb is provided but already contains
        # this specific record. This operation should raise
        cursor.execute("SELECT main_index.* "
                       "FROM main_index "
                       "INNER JOIN indexed_files "
                       "ON main_index.fileid = indexed_files.fileid "
                       "WHERE indexed_files.filename=? AND "
                       "main_index.recordoffsetstart=?;",
                       (basename(self._handlename),
                       indexdict["recordoffsetstart"]))
        samefileposition = cursor.fetchone()
        cursor.execute("SELECT main_index.* "
                       "FROM main_index "
                       "INNER JOIN indexed_files "
                       "ON main_index.fileid = indexed_files.fileid "
                       "WHERE indexed_files.filename=? AND "
                       "main_index.id=?;",
                       (basename(self._handlename),
                       indexdict["id"]))
        con.commit()
        sameid = cursor.fetchone()
        if samefileposition is not None or sameid is not None:
            raise ValueError("indexdb already contains a similar record")

        #create main index input table
        # using 'insecure' string format query generation
        # due to the lack of flexible substitution syntax
        keys, placehold, values = "fileid", "?", [fileid]
        for key, value in indexdict.items():
            keys += ", " + key
            placehold += ",?"
            values.append(value)
        values = tuple(values)
        cur = con.cursor()
        cur.execute("INSERT INTO main_index ({0})".format(keys) +\
                    "VALUES ({0})".format(placehold), values)
        indexid = cur.lastrowid
        #increment the file record counter
        cur.execute("UPDATE indexed_files SET count=count + 1 "
                    "WHERE filename=?", (basename(self._handlename),))
        con.commit()
        con.close()
        self._index_id = indexid
        self.record.update(indexdict)
        index_exists = True

    def _create_tables(self, con, indexdict):
        """create metadata table and populate format

        For a blank database, a valid indexdict is used as the template
        for creating the correct rows in the record index.

        Args:
          con (sqlite db connection): a connection to the sqlite db
          indexdict (dict): A dictionary corresponding to an indexed
            sequence record. This is a template for the main_index
            table
        """
        #set locking mode and synchronous
        con.execute("PRAGMA synchronous=0")
        con.execute("PRAGMA locking_mode=EXCLUSIVE")
        con.commit()

        con.execute( \
            "CREATE TABLE meta_data(key TEXT UNIQUE, value TEXT);")
        con.execute("INSERT INTO meta_data (key, value) VALUES (?,?);",
                    ("format", self._format))
        #create file table
        con.execute("CREATE TABLE indexed_files(fileid INTEGER PRIMARY "
                    "KEY, filename TEXT UNIQUE, count INTEGER);")

        #create basic index table:
        # The (less secure) string manipulation is used because sqlite3
        # API does not provide rich logic in the parameter substitution
        # query interface. The variant nature of each file's index
        # requires that index db's are format specific and may contain
        # a variable number of fields
        rows = sorted(list(indexdict.keys()))
        indextab = ", ".join([str(key) + " INTEGER" for key in rows \
                              if key != "id"])
        indextab = "CREATE TABLE main_index(" +\
                   "indexid INTEGER PRIMARY KEY, " +\
                   "id TEXT, fileid INTEGER, " + indextab +\
                   ", FOREIGN KEY(fileid) REFERENCES indexed_files(fileid));"
        con.execute(indextab)

        #create features table
        feattab = "CREATE TABLE features(" +\
                  "fileid INTEGER, " +\
                  "indexid INTEGER, " +\
                  "offsetbegin INTEGER, " +\
                  "offsetend INTEGER, " +\
                  "seqbegin INTEGER, " +\
                  "seqend INTEGER, " +\
                  "qualifier TEXT, " +\
                  "FOREIGN KEY(fileid) REFERENCES indexed_files(fileid) " +\
                  "FOREIGN KEY(indexid) REFERENCES main_index(indexid));"
        con.execute(feattab)
        con.commit()

    def _get_handle_id(self, con, write=False):
        """Return the file id from the DB, write an entry if necessary

        returns fileid for a filename. If write is False a fileid not
        found in the database will return None while setting write to
        True will require that the file is written prior to returning
        the fileid

        Args:
          con (sqlite db connection): a connection to the sqlite db
          write (bool): if write is True, the handle id will be added
            to the indexed_files table.
        """
        if self._handle_id is not None:
            return self._handle_id
        if self._handle_id is None:
            cursor = con.cursor()
            name = cursor.execute("SELECT fileid "
                           "FROM indexed_files WHERE "
                           "filename=?;",
                           (self._handlename,))
            name = name.fetchone()
            if name is not None:
                self._handle_id = name[0]
                return name[0]
            elif write:
                cursor.execute("INSERT INTO indexed_files "
                               "(filename, count) VALUES (?, ?);",
                               (self._handlename, 0))
                name = cursor.execute("SELECT fileid "
                           "FROM indexed_files WHERE "
                           "filename=?;",
                           (self._handlename,)).fetchone()
                self._handle_id = name[0]
                con.commit()
                return name[0]
            else:
                return None

    def _is_valid_db_with_tables(self, con):
        """if db is empty False, if correct tables True, else raise

        Args:
          con (sqlite db connection): a connection to the sqlite db
        """
        cursor = con.cursor()
        #case: database is empty
        tablenames = cursor.execute(
            "SELECT name FROM sqlite_master WHERE "
            "type='table' ORDER BY name;").fetchall()
        if len(tablenames) == 0:
            return False
        #case: database is valid and has correct format
        tablenames = [name[0] for name in tablenames]
        expectedtables = ['features', 'indexed_files',
                          'main_index', 'meta_data']
        if tablenames == expectedtables:
            #quick check: format matches self._format
            format, = cursor.execute(
                "SELECT value FROM meta_data WHERE key=?;",
                ("format",)).fetchone()
            if format != self._format:
                raise ValueError("provided database is for % files"%(format,))
            return True
        #something was wrong with the database
        else:
            raise ValueError("provided database is incomplete or corrupt")


    def set_feature_index(self, feature_index_list):
        """Save feature index to the database and set to memory

         Args:
          feature_index_list (FeatureBinCollection): All features made
            by the feature indexer
        """

        if not isinstance(feature_index_list, FeatureBinCollection):
            raise ValueError("_feature_index must be a FeatureBinCollection")
        self._features = feature_index_list
        if self._indexdb is None:
            return

        #from this point, do database interaction
        feature_index_list = feature_index_list[:]
        con = _get_db_connection(self._indexdb)
        if not self._is_valid_db_with_tables(con):
            raise ValueError("_index must be set before _feature_index")

        cursor = con.cursor()
        fileid = self._get_handle_id(con, write=False)

        #this will raise a NameError if the _index_id has not been set
        assert self._index_id is not None
        cursor.execute("SELECT features.* "
                       "FROM features "
                       "INNER JOIN main_index mi "
                       "ON mi.indexid = features.indexid "
                       "WHERE features.indexid=? "
                       "LIMIT 1;",
                       (self._index_id,))
        samefileposition = cursor.fetchone()

        if samefileposition is not None:
            raise ValueError("indexdb already contains a similar record")

        #create featuresdict table
        # using 'insecure' string format query generation
        # due to the lack of flexible substitution syntax
        for feature in feature_index_list:
            begin, end, beginoffset, endoffset, qualifier = feature
            #these assertions are done because SQLite is not static typed
            assert _is_int_or_long(begin)
            assert _is_int_or_long(end)
            assert _is_int_or_long(beginoffset)
            assert _is_int_or_long(endoffset)
            cursor.execute("INSERT INTO features "
                           "(indexid, fileid, offsetbegin, offsetend, "
                           "seqbegin, seqend, qualifier) "
                           "VALUES (?, ?, ?, ?, ?, ?, ?);",
                           (self._index_id, fileid, beginoffset,
                           endoffset, begin, end, qualifier))
        con.commit()
        con.close()

    def get_features(self, begin, end):
        """Return all features contained the defined range

        Args:
          begin (int): inclusive beginning of features to fetch
          end (int): exclusive endpoint of features to fetch
        """
        truebegin, trueend = 0, self.record["seqlen"]

        if isinstance(self._features, FeatureBinCollection):
            return self._features

        #case: index is not loaded and no db has been used
        if self._indexdb is None:
            #case: no index key provided thus no record will be searched
            if self._recordkey is None:
                return None

        #case: retrieve from db because a key was given and a db exists
        #      this will connect to the db (this is implicit since the
        #      _index getter must have already run successfully)
        con = _get_db_connection(self._indexdb)
        cursor = con.cursor()

        #case: empty db, this is already dealt with in _index getter

        #case: indexdb is provided and has records
        featureindexes = cursor.execute("SELECT f.seqbegin, "
            "f.seqend, f.offsetbegin, f.offsetend, f.qualifier "
            "FROM features f "
            "WHERE f.indexid = ? "
            "AND f.seqbegin >= ? "
            "AND f.seqend <= ? "
            "ORDER BY f.seqbegin;",
            (self._index_id, begin, end))

        container = FeatureBinCollection()
        for featureindex in featureindexes:
            container.insert(featureindex)

        con.commit()
        con.close()

        if begin == truebegin and end == trueend:
            self._features = container
        return container

class SeqRecordProxyBase(SeqRecord):
    """A SeqRecord object holds a sequence and information about it.

    Main attributes:
     - id          - Identifier such as a locus tag (string)
     - seq         - The sequence itself (Seq object or similar)

    Additional attributes:
     - name        - Sequence name, e.g. gene name (string)
     - description - Additional text (string)
     - dbxrefs     - List of database cross references.
     - features    - Any (sub)features, as SeqFeature objects
     - annotations - Further information about the whole sequence.
                     Most entries are strings, or lists of strings.
     - letter_annotations - Per letter/symbol annotation (restricted
                     dictionary). This holds Python sequences (lists,
                     strings or tuples) whose length matches that of
                     the sequence. A typical use would be to hold a
                     list of integers representing sequencing quality
                     scores, or representing the secondary structure.

    You will typically use Bio.SeqIO to read in sequences from files
    as SeqRecord objects.  However, you may want to create your own
    SeqRecord objects directly (see the __init__ method for details):

    #subclassing to-do list
    attributes that must be managed or utilized
     -   _seq         ::  sequence type, set by _read_seq
     -   name
     -   dbxrefs      ::  [] db cross references
     -   id           ::  id string
     -   _features    ::  list of SeqFeatures
     -   _annotations ::  dict of annotations
     -   _letter_an...::  per letter information

     -   _index_begin ::  defines begin index w/r/t global seq position
     -   _index_end   ::  defines end index w/r/t global seq position

    methods need implementation by derived class:
     - _format               :: a text string of the format
     - _read_seq             :: gets the sequence from file and sets it
     - _make_record_index    :: index the nce portion
              the record index is primarily format specific
              but it must contain the keys "id" and "seqlen"
     - _read_features        :: read features from ft. index
     - _make_feature_index   :: make feature index
     - _load_non_lazy_values :: load all values from lazy

    """

    _seq = None
    _format = None
    _features = None
    _feature_index = None
    _per_letter_annotations = {}
    annotations = {}

    def __init__(self, handle, startoffset=None, indexdb=None,
                 indexkey=None, alphabet=None):

        self._handle = handle
        self._alphabet = alphabet
        self._indexdb = indexdb
        self._indexkey = indexkey
        if hasattr(handle, "name"):
            name = basename(handle.name)
        else:
            name = None
        self._index = SeqProxyIndexManager(self._format,
            indexdb=indexdb, recordkey=indexkey,
            handlename=name)

        # create the index or load an existing one from file
        if self._index.index_exists is False:
            # make main index and sets _index
            new_index = {"recordoffsetstart":startoffset}
            new_index = self._make_record_index(new_index)
            self._index.set_record_index(new_index)
            # make feature index and sets _feature_index
            feature_index = FeatureBinCollection()
            feature_index = self._make_feature_index(feature_index)
            self._index.set_feature_index(feature_index)

        #set base values
        self._load_non_lazy_values()
        if self.id != self._index.record["id"]:
            raise ValueError("The index does not contain the correct id.")
        self._index_begin = 0
        self._index_end = self._index.record["seqlen"]

    def _load_non_lazy_values(self):
        """this is implemented to save the non-lazy values from a seq file"""
        raise NotImplementedError( \
            "_load_non_lazy_values must be implemented in the derived class")

    def _make_record_index(self, new_index):
        """this is implemented to set the index from a seq file"""
        raise NotImplementedError( \
            "_make_record_index must be implemented in the derived class")

    def _make_feature_index(self, new_index):
        """this is implemented to set the feature index from a seq file"""
        raise NotImplementedError( \
            "_make_feature_index must be implemented in the derived class")

    def _return_seq(self):
        """(private) removes getter logic from _read_seq"""
        if self._seq is None:
            self._read_seq()
        return self._seq

    def _set_seq_by_user(self, newseq):
        """allow the user to set a sequence, avoid some errors"""
        if newseq is None:
            raise ValueError("Setting seq to 'None' is reserved "
                             "for lazy proxy.")
        if len(newseq) != len(self):
            raise ValueError("New sequence must match record length.")
        self._seq = newseq

    def _read_seq(self):
        """this is implemented to handle file access for setting _seq"""
        raise NotImplementedError( \
            "_read_seq must be implemented in the derived class")

    seq = property(fget=_return_seq,
                   fset=_set_seq_by_user,
                   doc="The sequence itself, as a Seq or MutableSeq object.")

    def _return_features(self):
        """(private) removes getter logic from _read_seq"""
        if self._features is None:
            self._read_features()
        return self._features

    def _set_features_by_user(self, newfeatures):
        """allow users to change _features"""
        if newfeatures is None:
            raise ValueError("Setting features to 'None' is reserved "
                             "for lazy proxy")
        self._features = newfeatures

    def _read_features(self):
        """this is implemented in the derived class to set _seq"""
        raise NotImplementedError( \
            "_read_properties must be implemented in the derived class")

    features = property(fget=_return_features,
                        doc="list of SeqFeature objects")

    def __getitem__(self, index):
        """Returns a sub-sequence or an individual letter.

        This will need to cleverly either return a sequence letter
        or a copy of itself with new marker indices.
        """
        if isinstance(index, int):
            #The recursive call is required to prevent full parsing when
            # only calling a single resiude
            return self[index:index+1].seq[0]

        elif isinstance(index, slice):
            #Raise an error if there is no seq to access.
            parent_length = len(self)
            start, stop, step = index.indices(parent_length)
            if parent_length <= 0:
                raise ValueError("Cannot slice sequence length <= 0.")
            #For step values != 1, return a SeqRecord instance mimicing
            # the return behavior of SeqRecord. Some lazy loading is still
            # performed by pre-narrowing the sequence window
            if step != 1:
                indexsmall = min(start, stop)
                indexlarge = max(start, stop)
                if self.seq is None:
                    raise ValueError("If the sequence is None, "
                                     "we cannot slice it.")
                return SeqRecord(self[indexsmall:indexlarge].seq[::step],
                                 id=self.id,
                                 name=self.name,
                                 description=self.description)
            elif step == 1:
                #this is what will be returned
                seq_proxy_copy = copy(self)
                #do some index math
                seq_proxy_copy._index_begin = self._index_begin + start
                seq_proxy_copy._index_end = self._index_begin + stop
                #fix _seq property when set
                if self._seq:
                    seq_proxy_copy._seq = self._seq[start:stop]

                #fix _features property when set
                if self._features:
                    seq_proxy_copy._features = []
                    for f in self._features:
                        if f.ref or f.ref_db:
                            #TODO - Implement this (with lots of tests)?
                            import warnings
                            warnings.warn("When slicing SeqRecord objects, "
                                  "any SeqFeature referencing other "
                                  "sequences (e.g. from segmented GenBank"
                                  " records) are ignored.")
                            continue
                        if start <= f.location.nofuzzy_start \
                            and f.location.nofuzzy_end <= stop:
                            seq_proxy_copy._features.append(f._shift(-start))
                return seq_proxy_copy

        raise ValueError("Invalid index")

    def __len__(self):
        """Returns the length of the sequence.

        The derived class must set attribute _seq_len
        """
        return self._index_end - self._index_begin

    def upper(self):
        """Returns copy of the record with an upper case sequence. """
        if not self._seq:
            self._read_seq()
        newself = copy(self)
        newself._seq = self._seq.upper()
        return newself


    def lower(self):
        """Returns copy of the record with a lower case sequence. """
        if not self._seq:
            self._read_seq()
        newself = copy(self)
        newself._seq = self._seq.lower()
        return newself

    def __repr__(self):
        """this is a shortened repr for the lazy loading parsers

        modification of the repr() magic method is required to prevent
        simple command line options from unnecessarily invoking full
        parsing of a file.

        This strategy of redefining __repr__ contrasts with the __str__
        method that is left to the base class since invoking the str()
        is less used and more likely in output focused programs.
        """
        idrepr = "id=%s" % str(self.id)
        namerepr = "name=%s" % str(self.name)
        descriptionrepr = "description=%s" % str(self.description)
        dbxrefsrepr = "dbxrefs=%s" % repr(self.dbxrefs)
        if self._seq is None:
            seqrepr = "seq=NOT_READ"
        else:
            seqrepr = "seq=%s" % repr(self.seq)
        return self.__class__.__name__ \
         + "(%s, %s, %s, %s, %s)" \
         % (seqrepr, idrepr, namerepr, descriptionrepr, dbxrefsrepr)

    def next_record_offset(self):
        """return the offset of the next record"""
        return self._index.record["nextrecordoffset"]

    def get_raw(self):
        """Get the raw text of the lazy record"""
        begin_offset = self._index.record["recordoffsetstart"]
        recordoffsetlength = self._index.record["recordoffsetlength"]
        self._handle.seek(begin_offset)
        return self._handle.read(recordoffsetlength)




    # All methods tagged below are implemented in the base class
    #
    #def __bool__(self):
    #__nonzero__= __bool__
    #def reverse_complement(self, id=False, name=False, description=False,
    #def __radd__(self, other):
    #def format(self, format):
    #def __format__(self, format_spec):
    #def __add__(self, other):
    #    returns non-proxy class where possible
    #def __radd__(self, other):
    #    returns non-proxy class where possible
    #def __str__(self):
    #def __repr__(self):
    #def __contains__(self, char):
    #def __iter__(self):

def inherit_lazy_method_doc(meth):
    """Use as a decorator to fix docstrings of lazy derived methods"""
    methodname = meth.__name__
    if hasattr(SeqRecordProxyBase, methodname):
        meth.__doc__ = getattr(SeqRecordProxyBase, methodname).__doc__
    else:
        raise ValueError("'{0}' is not defined in SeqRecordProxyBase"\
                         .format(methodname))
    return meth

class HandleQueueLRU(object):
    """A simplified LRU cache that serves handles to HandleWrapper

    Because this LRU cache uses a simple list to store the cached
    handles, this cache operates at O(n) with respect to the number
    of cached handles. For a deep cache a different implementation
    would be needed, but for this specific application where the
    number of concurrently open files is very small, this is fine
    """
    def __init__(self, filenames, cachelen=5):
        if cachelen > 10:
            raise ValueError("cachelen must not exceed 10")
        if isinstance(filenames, str):
            raise ValueError("filenames must be a list, not a str")
        self.cachelen = min(len(filenames), cachelen)

        #initialize namedict and fill the cache
        self.namedict = {}
        self.stack = []
        stackinitializerlen = 0
        for fname in filenames:
            value = {"fullpath":fname, "handle":None}
            self.namedict[basename(fname)] = value
            if stackinitializerlen < self.cachelen:
                value["handle"] = open(value["fullpath"], 'rb')
                self.stack.append(value)
                stackinitializerlen += 1

    def __getitem__(self, key):
        """Get a file handle from the LRU cache

        This if the handle is at the top of the cache, it is returned
        quickly without making any changes. If the handle is deeper in
        the cache it is returned after being pushed to the top.

        Handles not in the cache will be created (via opening in binary
        writing mode) and they will be added to the top of the cache.
        The handle at the bottom of the cache will be closed and the
        reference will be terminated.
        """
        if key not in self.namedict:
            raise KeyError("HandleQueueLRU does not contain '{}'".format(key))
        record = self.namedict[key]
        #case: record is the most recent
        if record == self.stack[-1]:
            pass
        #case: record isn't in the stack
        elif record["handle"] is None:
            oldrec = self.stack.pop(0)
            oldrec["handle"].close()
            oldrec["handle"] = None
            record["handle"] = open(record["fullpath"], 'rb')
            self.stack.append(record)
        #case: record is in the stack but not at the top
        else:
            recordpos = self.stack.index(record) #this is O(n) for stack len
            self.stack.pop(recordpos) #this is O(n) for stack len
            self.stack.append(record)
        return record["handle"]

    def __contains__(self, key):
        if key in self.namedict:
            return True
        return False

class HandleWrapper(object):
    """Simple tool to wrap handles being provided by the queue

    This will only work when duck typing is used since __class__
    will return HandleWrapper no.
    """
    def __init__(self, filekey, handle_queue):
        self.handle_queue = handle_queue
        self.filekey = filekey

    def __getattr__(self, requested_attribute):
        if self.filekey in self.handle_queue:
            realhandle = self.handle_queue[self.filekey]
            return getattr(realhandle, requested_attribute)

class LazyIterator(object):
    """wrapper for most lazy-loading/indexing access routes for SeqIO

    This class is used to combine most workflows for the lazy loading
    module. Use of this will allow SeqIO.parse, and SeqIO.index to
    use the contents of the lazy loading module through sparse code
    additions
    """

    def __init__(self, files, return_class, index=True,
                 alphabet=None, asdict=False, key_function=None):
        #set up files and file keys
        if isinstance(files, list) and len(files) > 0:
            self.files = files
            self.filekeys = [basename(f) for f in files]
        else:
            raise ValueError("files must be a non-empty list")

        #set key_function
        if not hasattr(key_function, "__call__"):
            key_function = lambda keyin: keyin

        #other default values
        self.return_class = return_class
        self.alphabet = alphabet
        self.index = index
        self.format = return_class._format
        self.asdict = asdict

        #check and fix index existence
        if isinstance(index, bool):
            self.use_an_index = False
        elif isfile(index):
            self.use_an_index = True
        else:
            #make a new file
            temp = open(index, 'wb')
            temp.close()
            self.use_an_index = True

        if self.use_an_index:
            #since an index is being used, the HandleQueue is needed
            for f in files:
                assert isfile(f), "File '{}' doesn't exist".format(f)
            handle_queue = HandleQueueLRU(files)
            self.handles = {}
            for f in files:
                fkey = basename(f)
                self.handles[fkey] = HandleWrapper(fkey, handle_queue)
            #Check if the database is empty and which files are indexed
            con = _get_db_connection(self.index)
            table_exists = con.execute("SELECT name FROM sqlite_master "
                "WHERE type='table' AND name='indexed_files'")
            table_exists = True and [val for val in table_exists]
            if table_exists:
                filesquery = con.execute("SELECT filename FROM indexed_files;")
                return_keys = [_as_string(key[0]) for key in filesquery]
            else:
                return_keys = []
            con.close()
            for fname in self.files:
                key = basename(fname)
                if key not in return_keys:
                    temphandle = open(fname, 'rb')
                    try:
                        _make_index_db(handle=temphandle,
                                       return_class=self.return_class,
                                       indexdb=self.index,
                                       format=self.format,
                                       alphabet=self.alphabet)
                    except ValueError as e:
                        if "corrupt" in str(e):
                            raise ValueError("database is corrupt, please"
                                " delete it and make a new one", e)
                        else:
                            raise e
                    temphandle.close()
            #set keys in light of key_function
            realkeylist = self._get_keys_from_db()
            self._keymap = {}
            self._keys = []
            for key in realkeylist:
                modifiedkey = key_function(key)
                if modifiedkey in self._keymap:
                    raise KeyError("key_function induces key conflict; " +
                        "{0}, and {1} both produce {2}"\
                        .format(self._keymap[modifiedkey], key, modifiedkey))
                self._keymap[key_function(key)] = key
                self._keys.append(modifiedkey)

    def __iter__(self):
        """This can do both dict-type key iteration or list iteration

        when this class is invoked with asdict=True, this iterator
        will iterate over keys, if this is invoked with asdict=False
        this iterator will iterate over records either pulling from
        the captive db, or by indexing on invocation.
        """
        return_class = self.return_class
        if self.use_an_index:
            for k in self._keys:
                if self.asdict:
                    yield k
                else:
                    yield self[k]
        else:
            for f in self.files:
                handle = open(f, 'rb')
                record_offset = _get_first_record_start_offset(handle,
                    file_format=self.format)
                while record_offset is not None:
                    result = return_class(handle, startoffset=record_offset,
                        indexdb=None, indexkey=None, alphabet=self.alphabet)
                    record_offset = result.next_record_offset()
                    yield result

    def _get_keys_from_db(self):
        """For a list of handle name, return all record keys"""
        if not self.use_an_index:
            raise RuntimeError("No dict-type access without a db-index")
        else:
            self.record_to_file = {}
            keys = []
            con = _get_db_connection(self.index)
            for f in self.files:
                fname = basename(f)
                tempkeys = con.execute("SELECT idx.id "
                    "FROM main_index idx "
                    "INNER JOIN indexed_files idxf "
                    "ON idxf.fileid = idx.fileid "
                    "WHERE idxf.filename = ?",
                    (fname,))
                tempkeys = [_as_string(k[0]) for k in tempkeys]
                for key in tempkeys:
                    self.record_to_file[key] = fname
                keys.extend(tempkeys)
            con.close()
            if len(keys) != len(set(keys)):
                raise KeyError("Files contain overlapping records and "
                               "can only be loaded separately.")
            return keys

    def get(self, key, d=None):
        """D.get(k[,d]) -> D[k] if k in D, else d. d defaults to None."""
        if not self.use_an_index:
            raise RuntimeError("An index file is required to use 'get'")
        if key in self._keys:
            return self[key]
        else:
            return d

    def get_raw(self, key):
        """D.get_raw(k), returns raw file data for record [k]."""
        if not self.use_an_index:
            raise RuntimeError("An index file is required to use 'get_raw'")
        rec = self[key]
        return rec.get_raw()

    def items(self):
        """D.items() -> a set-like object providing a view on D's items"""
        if not self.use_an_index:
            raise RuntimeError("An index file is required to use 'items'")
        for k in self._keys:
            yield (k, self[k])

    iteritems = items

    def keys(self):
        """"D.keys() -> a set-like object providing a view on D's keys"""
        if not self.use_an_index:
            raise RuntimeError("An index file is required to use 'keys'")
        return self._keys

    def __contains__(self, key):
        if not self.use_an_index:
            raise RuntimeError("An index file is required for '__contains__'")
        if key in self.record_to_file:
            return True
        return False

    def __getitem__(self, key):
        if not self.use_an_index:
            raise RuntimeError("An index file is required for '__getitem__'")
        recordid = self._keymap[key]
        filekey = self.record_to_file[recordid]
        handle = self.handles[filekey]
        return self.return_class(handle, indexkey=recordid,
                                 indexdb=self.index,
                                 alphabet=self.alphabet)

def _make_index_db(handle, return_class, indexdb, format, alphabet=None):
    """(private) populate index db with file's index information"""
    record_offset = _get_first_record_start_offset(handle, format)
    while record_offset is not None:
        temp = return_class(handle, startoffset=record_offset,
                            indexdb=indexdb, indexkey=None, alphabet=alphabet)
        record_offset = temp.next_record_offset()

def _get_first_record_start_offset(handle, file_format=None):
    """(private) return the offset of the first record, or None"""
    marker = {"fasta": ">",
              "genbank": "LOCUS",
              "gb": "LOCUS",
              "uniprot-xml": "<",
              "embl": "ID   "}[file_format]
    marker_re = re.compile(marker)

    # Set handle to beginning and set offsets to valid initial values
    handle.seek(0)
    current_offset, working_offset = -1, 0
    while True:
        # Advance the current_offset markers, using handle.readline().
        # unlike iterators, reaching the end of the handle and calling
        # readline will not raise StopIteration, instead it will return
        # a an empty string and the file offset will not advance
        current_offset = working_offset
        currentline = _bytes_to_string(handle.readline())
        working_offset = handle.tell()

        # Check if we encounter the first record then save the offset
        if marker_re.match(currentline):
            return current_offset
        if not currentline:
            return None

class FeatureBinCollection(object):
    """Manage efficient creation and retrieval of feature indices

    This class is used to organize feature data in a quickly
    retrievable data structure. The feature data must be added as a
    tuple containing at least two indices: first annotated residue
    and the last as a half open half closed interval [first, last).
    The indices are assumed to be the first two elements of the stored
    tuple, but they may be re-assigned on instantiation via the
    beginindex and endindex kwarks.

    EXAMPLE
    -------
    defined below is this 3-tuple: (beginindex, endindex, fileidx),
    three features are added to a newly initialized featurebin.

    >>> ft0 = (5574, 5613, 2300)
    >>> ft1 = (0, 18141, 1300 )
    >>> ft2 = (5298, 6416, 3540)
    >>> featurebin = FeatureBinCollection(bounded_only_returns=False)
    >>> featurebin.insert( ft0 )
    >>> featurebin.insert( ft1 )
    >>> featurebin.insert( ft2 )
    >>> len(featurebin)
    3

    Now that the 'featurebin' instance has some features, they can be
    retrieved with a standard getter using single integer indices or
    slice notation.

    >>> featurebin[1]
    [(0, 18141, 1300)]
    >>> sliceresult = featurebin[5200:5300]
    >>> sliceresult.sort()
    >>> sliceresult
    [(0, 18141, 1300), (5298, 6416, 3540)]


    BACKGROUND:
    -----------
    The basic idea of using feature bins is to group features into
    bins organized by their span and sequence location. These bins then
    allow only likely candidate features to be queried rather than all
    features. The example below illustrated with Figure 1 shows a
    similar scheme where feature1 is stored in bin-0, feature2 in
    bin-4, and feature3 in bin-2. Each sequence is stored in the
    smallest bin that will fully contain the sequence. A query of all
    features in the region denoted by query1 could be quickly performed
    by only searching through bins 0, 2, 5, and 6. Were this data
    structure many levels deep, the performance savings would be large.
    It is also useful to see that any bin fully bound by a query can
    have the full contents returned without individual comparison
    operations.

    ___Figure 1_________________________________________________
    |                                                           |
    |    feature1  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~              |
    |    feature2  |       ~~~~                  |              |
    |    feature3  |       |  |               ~~~~~~~~~~~~~     |
    |              |       |  |               |  |        |     |
    | bins:        |       |  |               |  |        |     |
    |    0_________|_______|__|_______._______|__|____.___|_    |
    |    1_________________|__|__   2_:_______|_______:___|_    |
    |    3__________  4____|__|__   5_:________  6____:_____    |
    |                                 :               :         |
    |                                 :               :         |
    |    query1                       [ ? ? ? ? ? ? ? ]         |
    |...........................................................|

    Further reading on the math behind the idea can be found in:
        Journal:  Bioinformatics Vol 27 no. 5 2011, pages 718-719
        Article:  "Tabix: fast retrieval of sequence features from
                  generic tab delimited files"
        Author:   Heng Li

    The implementation by Li has, as its largest bin ~500 million
    (2^29) and its smallest bin ~16000 (2^14). Each level of binning
    is separated by a factor of 8 (2^3). The implementation herein
    abandons a static binning scheme and instead starts with the
    smallest and largest bins as 256 and 8 million respectively. These
    bins can then be dynamically expanded increasing by a factor of 8
    every time new data is found to be larger than the largest bin.
    As a practical matter of sanity checking, bin sizes are capped at
    2.2 trillion residues (2^41).

    Under some circumstances the exact size of a sequence and all
    related annotations is known beforehand. If this is the case the
    length kwarg allows the binning object to be solidified on
    instantiation at the correct length. Attempts to add data outside
    of this range will raise an error.

    This structure knows nothing about the global sequence index and
    is indexed at zero. Any index transformation must be done at a
    higher level. It is important that all sequences and features
    stored here are indexed to zero.
    """

    def __init__(self, length=None, beginindex=0, endindex=1,
                 bounded_only_returns=True):
        """ initialize the class and set standard attributes

        kwargs:

          length:
            when length == None, the bins are dynamically sized.
            when length is a positive integer, the appropriate bin
            size is selected and locked. Exceeding this value will
            cause exceptions when the max bin size is locked

          beginindex:
            the index of the first residue within the tuple that will
            be stored with the FeatureBinCollection.

          endindex:
            the index of the last residue (as a open interval) inside
            the tuple that will be stored with FeatureBinCollection
        """
        #these should not be changed
        self._bin_level_count = 6
        self._bins = [[] for i in range(37449)]

        # this defines the indices of the begin and end sequence info
        # in the tuple structures stored in the bins
        self.bounded_only_returns = bounded_only_returns
        self._beginindex = beginindex
        self._endindex = endindex

        #default action: start small (8M) and allow expansion
        self._sorted = False
        self._dynamic_size = True
        if length is None:
            self._set_max_bin_power(23)

        #alternate action if a sequence length is provided
        # set to smallest power able to fully contain
        elif _is_int_or_long(length) and length > 0:
            default_powers = [23, 26, 29, 32, 35, 38, 41]
            for power in default_powers:
                if length <= 2**power:
                    self._set_max_bin_power(power)
                    self._dynamic_size = False
                    break
            if self._dynamic_size: #this should have been set to False
                error_string = "Sequence length is {0}:".format(length) +\
                               " must be less than 2^41"
                raise ValueError(error_string)

    def _increase_bin_sizes(self):
        """increase max bin size 8x and reorganize existing binned data

        In order to increase the total maximum bin size, the set with
        smallest bins is merged with the next level of larger bins.
        Next, the entire set moved down a single level in the hierarchy
        and the top level bin (the largest) is opened as an empty list.

        An assertion in this routine blocks sequences larger than 2**41
        from being created. This is a soft limit based on the expected
        largest possible sequence. Larger sequences can be supported
        by modifiying this assertion.
        """
        oldsizepower = self._max_bin_power
        newsizepower = oldsizepower + 3
        assert newsizepower <= 41
        self._set_max_bin_power(newsizepower)

        # first, remove the lowest level
        # by merging it up to the previous level
        level = 5
        oL = int((2**(3*level) - 1)/7.0)
        new_level = 4
        new_oL = int((2**(3*new_level) - 1)/7.0)
        old_size = 2**(oldsizepower - 3*level)
        new_size = 2**(oldsizepower - 3*new_level)
        for k in range(4681, 37449):
            bin_begin_old = (k-oL)*old_size
            k_new = int(floor(new_oL + (bin_begin_old/new_size)))
            #extend required to save existing data
            self._bins[k_new].extend(self._bins[k])
            self._bins[k] = []

        #then, move everything down.
        for k_inverse in range(4681):
            k = 4680 - k_inverse
            level = int(floor(log((7*k + 1), 2)/3.0))
            new_level = level + 1
            oL = int((2**(3*level) - 1)/7.0)
            new_oL = int((2**(3*new_level) - 1)/7.0)

            new_index = k - oL + new_oL
            self._bins[new_index] = self._bins[k]
            self._bins[k] = []

    def _set_max_bin_power(self, power):
        """Set the maximum bin power and fix other attributes"""

        self._max_bin_power = power
        self._min_bin_power = power - 3*(self._bin_level_count-1)
        self._size_list = [2**(self._min_bin_power+3*n) \
                          for n in range(self._bin_level_count)]
        self._max_sequence_length = self._size_list[-1]

    def insert(self, feature_tuple):
        """inserts a tuple with a sequence range into the feature bins

        data is assumed to be somewhat scrubbed, coming from a parser
        or a parser consumer."""

        beginindex = self._beginindex
        endindex = self._endindex

        #reset sorted quality
        self._sorted = False

        begin = feature_tuple[beginindex]
        end = feature_tuple[endindex]
        #use _py3k module later
        assert _is_int_or_long(begin)
        assert _is_int_or_long(end)
        assert begin <= end
        span = end-begin

        bin_index = self._calculate_bin_index(begin, span)
        self._bins[bin_index].append(feature_tuple)

    def __len__(self):
        """return the count of elements contained in all bins"""
        return sum(len(binn) for binn in self._bins)

    def sort(self):
        """perform bin-centric sorting for faster retrieval"""
        #bins must be sorted by the begin index, this is fastest
        if self._beginindex == 0:
            for i in range(len(self._bins)):
                self._bins[i].sort()
        #this is a bit slower but accomodates diverse data structures
        else:
            for i in range(len(self._bins)):
                self._bins[i].sort(key=lambda tup: tup[self._beginindex])
        #reset sorted quality
        self._sorted = True

    def __getitem__(self, key):
        """Efficiently retrieves the required entries, return as a list

        This getter primarily works as expected with the exception
        of is the treatment of slices indices where the start is
        greater than the stop. Rather than just throwing calculated
        output (generally no matches would be returned for sequence
        based data sets) an IndexError is raised.
        """
        #set some locals
        beginindex = self._beginindex
        endindex = self._endindex

        #check that it is a slice and it has no step property (or step==1)
        if not isinstance(key, slice):
            if _is_int_or_long(key):
                key = slice(key, key+1)
            else:
                raise TypeError("lookups in the feature bin"
                                " must use slice or int keys")

        #fix begin or end index for slicing of forms: bins[50:] or bins[:]
        keystart, keystop, keystep = key.indices(self._max_sequence_length)
        if keystart > keystop:
            raise IndexError("key not valid, slice.start > slice.stop")

        #any integers are just converted to a 'len() == 1' slice
        if keystep is not None and keystep != 1:
            raise KeyError("lookups in the feature bin may not use"
                           " slice stepping ex. bins[0:50:2]")

        #pre-sort if necessary
        if not self._sorted:
            self.sort()

        #code taken from self._calculate_bin_index(), comments removed
        return_entries = []
        bin_level_count = self._bin_level_count
        max_bin_power = self._max_bin_power
        for l_inverse in range(bin_level_count):
            L = bin_level_count - 1 - l_inverse
            oL = (2**(3*L) - 1)/7
            sL = float(2**(max_bin_power-3*L))
            k1 = int(floor(oL + (keystart/sL)))
            #k2 is incremented since range is used
            k2 = int(ceil(oL - 1 + (keystop)/sL)) + 1
            if k2-k1 > 2:
                for entry in range(k1+1, k2-1):
                    return_entries.extend(self._bins[entry])
            for binn in set([k1, k2-1]):
                #for binn in range(k1,k2):
                for feature in self._bins[binn]:
                    if self.bounded_only_returns:
                        #this covers fully bound sequence and left overlap
                        if keystart <= feature[beginindex] and \
                                       feature[endindex] <= keystop:
                            return_entries.append(feature)
                    else:
                        #this covers fully bound sequence and left overlap
                        if keystart <= feature[beginindex] < keystop:
                            return_entries.append(feature)
                        #this covers left sequence right sequence overlap
                        elif keystart < feature[endindex] <= keystop:
                            return_entries.append(feature)
                        #this covers seqyebces fully bound by a feature
                        elif keystart > feature[beginindex] and \
                             keystop < feature[endindex]:
                            return_entries.append(feature)
                        if keystop < feature[beginindex]:
                            break
        return return_entries

    def _calculate_bin_index(self, begin, span):
        """ Returns a bin index given a (begin, span) interval

        The equations in the Tabix paper (see class docstring)
        are used to determine the bin index

        This function should only be used privately in the context
        of having no easier relationship to assign bin index. Placing
        this in a loop for any task other than arbitrary assignments
        is a bad idea since many other the common tasks can provide
        bin index through more direct relationships.
        _increase_bin_sizes is an example of a routine that would
        suffer performance penalties were it to use this method
        yet can be run efficiently using alternate bin index
        relationships.
        """

        #take care base cases with the span parameter
        assert span >= 0
        #determination of bin location for zero length seq's
        span = max(1, span)

        # fix bin sizes if needed. Also, if the bin size is
        # larger than expected, do some self-consistency checks
        while begin+span > self._max_sequence_length:
            if self._dynamic_size:
                # len(seq) > 2.19 trillion is not reasonable
                assert begin+span <= 2**41
            elif not self._dynamic_size \
                and begin+span > 2**self._max_bin_power:
                error_string = "feature index at {0}:".format(begin+span) +\
                    " must be less than 2^{0}".format(self._max_bin_power)
                raise ValueError(error_string)
            self._increase_bin_sizes()

        #run the assignment loop
        bin_level_count = self._bin_level_count
        max_bin_power = self._max_bin_power
        for l_inverse in range(bin_level_count):
            level = bin_level_count - 1 - l_inverse
            # calculate offset (oL) of the list at level L
            offset_at_L = (2**(3*level) - 1)/7
            group_length = (2**(3*(level+1)) - 1)/7
            #calculate size (sL) of the list: the number of residues in width
            size_at_L = float(2**(max_bin_power-3*level))
            # interval[0] >= (k - oL)*sL
            # rearrange to form
            # k =< oL + (interval[0])/sL
            k1 = int(floor(offset_at_L + (begin/size_at_L)))
            # interval[1] < (k - oL + 1)*sL
            # rearrange to form
            #k > 1 + oL + (interval[1])/sL
            #k2 = oL - 1 + (begin+span)/sL
            k2 = int(ceil(offset_at_L - 1 + (begin+span)/size_at_L))
            if k1 == k2 and k1 < group_length:
                return k1

        assert False # the assignment loop failed
from xml.parsers import expat

def xml_index_iter(filename, targetfield, tagstoparse=None, returndict=False):
    """A xml file iter that returns indexes for sequential tags

    targetfield is a text field that defines the tag over which this
    function iterates. All targetfield tags must be sequential and any
    tag that interupts a set of targetfield tags will result in the end
    of iteration.

    tagstoparse
    """
    if not hasattr(filename, "read"):
        handle = open(filename, 'rb')
    else:
        handle = filename

    if tagstoparse is None:
        tagstoparse = []
    if targetfield not in tagstoparse:
        tagstoparse.append(targetfield)

    position = 0
    handler = ExpatHandler(handle, targetfield, tagstoparse)
    while True:
        root = handler.parse_from_position(position)
        if returndict:
            yield root.first_child().flatten_to_dict()
        else:
            yield root.first_child()
        if root.lastrecord is True:
            break
        position = root.nextelementoffset


class LinkedElement(object):
    """A simple to use LinkedElement class for indexing

    The LinkedElement is an element class that references it's
    parent and makes explicit the use of index positions.

    LinkedElement is similar to xml.etree.Element but with a simplified
    set of interactions. The main reason this was created is to
    add a reference to the parent element. Addding this reference
    to the default Element requires significant reworking of multiple
    operations so this custom implementation will actually reduce
    """
    def __init__(self, tag, attributes=None, begin=None, end=None):
        self.parent = None
        self.tag = tag
        self.text = ""
        if attributes is None:
            self.attributes = {}
        self.children = []
        self.indexbegin = begin
        self.indexend = end

        #used for looking at next element
        self.lastrecord = True
        self.nextelementoffset = None

    def append(self, child):
        child.parent = self
        self.children.append(child)

    def find_all_children_by_tag(self, tag):
        validchildren = []
        for child in self.children:
            if child.tag == tag:
                validchildren.append(child)
            validchildren.extend(child.find_children_by_tag(tag))
        return validchildren

    def find_children_by_tag(self, tag):
        validchildren = []
        for child in self.children:
            if child.tag == tag:
                validchildren.append(child)
        return validchildren

    def first_child(self):
        if not self.children:
            return None
        else:
            return self.children[0]

    def __repr__(self):
        return "< LinkedElement tag={0}, position=({1},{2}) >"\
                .format(self.tag, self.indexbegin, self.indexend)

    def depth(self):
        """Recursive call; mostly useful for debugging"""
        if self.parent is None:
            return 0
        else:
            return self.parent.depth() + 1

    def extract_from_handle(self, handle):
        handle.seek(self.indexbegin)
        segment = handle.read(self.indexend - self.indexbegin)
        if hasattr(segment, "decode"):
            segment = _bytes_to_string(segment)
        return segment

    def flatten_to_dict(self):
        asdict = {"tag":self.tag,
           "text":self.text,
           "attributes":self.attributes,
           "children":[child.flatten_to_dict() for child in self.children],
           "file_offset_begin":self.indexbegin,
           "file_offset_length":(self.indexend - self.indexbegin)}
        return asdict

class ExpatHandler(object):
    """ExpatHandler class will return an indexed LinkedElement tree

    The  targetfield attribute is the fundamental unit of indexing,
    these xml tags must be a sequential list with no tags in between
    having a different identity. The tagstoparse list contains tags
    that will be extracted into the LinkedElement tree.

    ExpatHandler assumes compilant well formmated XML, several types
    of formatting errors will result in difficult to decipher
    errors while other sorts of errors will not be detected. Best
    practice is to use a secondary parser for actual parsing and
    validation of data.
    """

    def __init__(self, handle, targetfield, tagstoparse,
                 parser_class=expat.ParserCreate):
        self._handle = handle
        self._parser_class = parser_class

        self.targetfield = targetfield
        self.tagstoparse = tagstoparse

        #values used in self.parse_from_position
        self.currentelem = None
        self.rootelem = None
        self.savetext = False
        self.baseposition = None
        self._parser = None

    def parse_from_position(self, position=0):
        """Parse XML from a position and return an indexed root node"""
        handle = self._handle
        #initialize the parser
        parser = self._parser = self._parser_class()
        parser.StartElementHandler = self.start_element
        parser.EndElementHandler = self.end_element
        parser.CharacterDataHandler = self.char_data
        handle.seek(position)
        self.baseposition = position

        rootelem = LinkedElement(tag="ROOT", begin=position)
        self.rootelem = rootelem
        self.currentelem = rootelem

        #make the index
        try:
            parser.ParseFile(handle)
        except StopIteration:
            pass
        # fix index for end-tags and next file begin
        handle.seek(self.rootelem.indexend - 1)
        readlen = 100
        padding = 0
        endfound = False
        given_end = self.rootelem.indexend - 1
        rootschild = rootelem.first_child()
        char = None
        while True:
            if not char:
                endregion = _bytes_to_string(handle.read(readlen))
                if not endregion or len(endregion) == 0:
                    raise ValueError( \
                        "file does not contain end tag on/after line {0}"\
                        .format(parser.CurrentLineNumber))
                char = endregion[padding%readlen]
            if char == "<":
                nexttag = given_end + padding
                self.rootelem.nextelementoffset = nexttag
                rootschild.nextelementoffset = nexttag
                if not endfound:
                    raise ValueError( \
                        "file does not contain end tag on/after line {0}"\
                        .format(parser.CurrentLineNumber))
                break
            padding += 1
            if char == ">":
                self.rootelem.indexend = given_end + padding
                rootschild.indexend = given_end + padding
                endfound = True

            char = endregion[padding%readlen]

        #check the next tag
        handle.seek(self.rootelem.nextelementoffset+1)
        beginregion = _bytes_to_string(handle.read(len(self.targetfield)))
        if self.targetfield == beginregion:
            self.rootelem.lastrecord = False
            rootschild.lastrecord = False


        if len(rootelem.children) == 0:
            raise ValueError("The XML @ offset={0} didn't contain a '{1}' tag"\
                              .format(position, self.targetfield))

        return rootelem

    def start_element(self, tag, attrs):
        """handle expat sart element"""
        if self.currentelem.indexend is True:
            self._finish_linked_element()

        if tag in self.tagstoparse:
            self.savetext = True
            byteindex = self._parser.CurrentByteIndex + self.baseposition
            new_element = LinkedElement(tag, begin=byteindex)
            new_element.attributes = attrs
            self.currentelem.append(new_element)
            self.currentelem = new_element
            if tag == self.targetfield:
                self.rootelem.indexbegin = byteindex
        else:
            self.savetext = False

    def end_element(self, tag):
        """handle expat end element"""
        if tag == self.targetfield:
            #for a compact xml file, this will produce the index of the end
            # tag without the trailing '>'. The parser fixes this.
            end = self._parser.CurrentByteIndex + len(self.targetfield) + 3 \
                  + self.baseposition
            self.currentelem.indexend = end
            self.rootelem.indexend = end
            raise StopIteration()
        if self.currentelem.indexend is True:
            self._finish_linked_element()
        if tag == self.currentelem.tag:
            self.currentelem.indexend = True

    def char_data(self, data):
        """handle expat character data"""
        if data.strip() and self.savetext:
            self.currentelem.text += data.strip()

    def _finish_linked_element(self):
        """ any LinkedElement eligible for finishing is saved here

        An LinkedElement has ended; fix the end byte index and fetch
        the parent node fixing it to currentelem.
        """
        assert self.currentelem.indexend is True
        position = self._parser.CurrentByteIndex + self.baseposition
        self.currentelem.indexend = position
        self.currentelem = self.currentelem.parent
