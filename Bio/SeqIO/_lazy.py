# Copyright 2009-2011 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Lazy parsing and feature indexing of sequence files (PRIVATE).

You are not expected to access this module, or any of its code, directly. This
is all handled internally by the following functions when accessed with the
lazy=True kwarg.
Bio.SeqIO.parse(.... lazy=True)
Bio.SeqIO.index(..., lazy=True)
Bio.SeqIO.read(..., lazy=True) 


The lazy loading parsers will read an entire record efficiently storing and 
parsing only the mimimum amount of information to provide fast lookup. Records
returned from parse(), index(), and read() will act as a proxy to a regular 
SeqRecord class.

The lazy loading strategy is, for large sequences, faster to initialize than
SeqIO.parse() and more memory efficient. It is slower to initialize than
SeqIO.index() but also more memory efficient.

The lazy loader will partially parse a sequence file noting the sequence
span and file position of sequence annoations. These annotations will
be stored for efficient retrieval from the file when requested. Sequence 
data will also be efficiently pulled from structured regions
and file

Note that this means all parsing is on demand, so improperly formatted
records may not trigger an exception unless the problem region is requested.
This is by design.
"""

from Bio._py3k import _is_int_or_long, _bytes_to_string, _as_bytes
from ..SeqRecord import SeqRecord, _RestrictedDict
from copy import copy
import re
from math import floor, ceil, log
from os.path import isfile, basename
try:
    from sqlite3 import dbapi2 as _sqlite
    from sqlite3 import IntegrityError as _IntegrityError
    from sqlite3 import OperationalError as _OperationalError
except ImportError:
    #Not present on Jython, but should be included in Python 2.5
    #or later (unless compiled from source without its dependencies)
    #Still want to offer in-memory indexing.
    _sqlite = None
    pass


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
     -   _seq         ::   sequence type, set by _read_seq
     -   name
     -   dbxrefs      ::   [] db cross references
     -   id           ::   id string
     -   _features    ::   list of SeqFeatures
     -   _annotations ::   dict of annotations
     -   _letter_an...::   per letter information

     -   _index_begin ::   defines begin index w/r/t global seq position
     -   _index_end   ::   defines end index w/r/t global seq position

     

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
    _features = None
    _feature_index = None

    def __init__(self, handle, startoffset=None, length=None,\
                 indexdb = None, indexkey = None, alphabet=None):

        self._handle = handle
        self._alphabet = alphabet
        self._indexdb = indexdb
        self._indexkey = indexkey
        
        # create the index or load an existing one from file 
        if self._index is None:
            # make main index and sets _index
            new_index = {"recordoffsetstart":startoffset, 
                           "recordoffsetlength":length}
            self._make_record_index(new_index)
            # make feature index and sets _feature_index
            feature_index = FeatureBinCollection()
            self._make_feature_index(feature_index)

        #set base values
        self._load_non_lazy_values()
        self._index_begin = 0
        self._index_end = self._index["seqlen"]
    
    def _return_seq(self):
        """(private) removes getter logic from _read_seq"""
        if not self._seq:
            self._read_seq()
        return self._seq
    
    def _read_seq(self):
        """this is implemented to handle file access for setting _seq"""
        raise NotImplementedError( \
            "_read_seq must be implemented in the derived class")

    seq = property(fget=_return_seq,
                   doc="The sequence itself, as a Seq or MutableSeq object.")

    def _return_features(self):
        """(private) removes getter logic from _read_seq"""
        if not self._features:
            self._read_features()
        return self._features
    
    def _read_features(self):
        """this is implemented to handle file access for setting _seq"""
        raise NotImplementedError( \
            "_read_properties must be implemented in the derived class")

    features = property(fget=_return_features,
                        doc="list of SeqFeature objects")

    #TODO - letter annotations?
    letter_annotations = None

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
                indexsmall = min(start,stop)
                indexlarge = max(start,stop)
                if self.seq is None:
                    raise ValueError("If the sequence is None, we cannot slice it.")
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

                return seq_proxy_copy

        raise ValueError("Invalid index")
    
    
    def __len__(self):
        """Returns the length of the sequence.
        
        The derived class must set attribute _seq_len
        """
        return self._index_end - self._index_begin


    def upper(self):
        """Returns a copy of the record with an upper case sequence.

        
        """
        if not self._seq:
            self._read_seq()
        newSelf = copy(self)
        newSelf._seq = self._seq.upper()
        return newSelf


    def lower(self):
        """Returns a copy of the record with a lower case sequence.
        """
        if not self._seq:
            self._read_seq()
        newSelf = copy(self)
        newSelf._seq = self._seq.lower()
        return newSelf

    def __repr__(self):
        """this is a shortened repr for the lazy loading parsers

        modification of the repr() magic method is required to prevent
        simple command line options from unnecessarily invoking full
        parsing of a file. 

        This strategy of redefining __repr__ contrasts with the __str__
        method that is left to the base class since invoking the str() is
        less used and more likely in output focused programs. 
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

    #
    # the following section defines behavior specific to index database io
    #

    __index = None
    __file_id = None

    def __load_and_return_index(self):
        """Returns the index dictionary from memory or from database
        
        When the database is not provided:
        This function will preferentially supply an index from memory.
        Since the _indexdb is None, it will provide None for the _index
        indicating to the lazy loading proxy that an index must be
        created.
        
        When the database (_indexdb) is provided:
        This function will still preferentially supply an index from memory.
        if the database is valid and empty it will provide None for the _index
        indicating to the lazy loading proxy that an index must be created. If
        an indexdb is present and valid, it will be queried to produce the
        _index dictionary which will be loaded into memory"""
        
        #case: index already exists in memory
        if isinstance(self.__index, dict):
            return self.__index
        
        #case: index is not loaded and no db has been used
        if self.__index is None:
            #case: no index key provided thus no record will be searched
            if self._indexkey is None:
                return None
            #case: a key request was made but no db was given
            if self._indexdb is None and self._indexkey is not None:
                raise KeyError("cannot retrieve by key from no db")
        #case: retrieve from db because a key was given and a db exists

        #this will connect to the db
        con = self.__connect_db_return_connection()
        cursor = con.cursor()
        
        #case: indexdb provided but it is empty
        if not self.__is_valid_db_with_tables(con):
            raise KeyError("_indexkey must match a stored record")
        #case: indexdb is provided and has records
        if _is_int_or_long(self._indexkey):
            recordindex = cursor.execute("SELECT main_index.* " +\
                           "FROM main_index " +\
                           "INNER JOIN indexed_files " +\
                           "ON main_index.fileid = indexed_files.fileid " +\
                           "WHERE indexed_files.filename=? " +\
                           "ORDER BY main_index.recordoffsetstart " + \
                           "LIMIT 2 OFFSET ?;",
                           (basename(self._handle.name),self._indexkey))
        elif isinstance(self._indexkey, str):
            recordindex = cursor.execute("SELECT main_index.* " +\
                           "FROM main_index " +\
                           "INNER JOIN indexed_files " +\
                           "ON main_index.fileid = indexed_files.fileid " +\
                           "WHERE indexed_files.filename=? AND " +\
                           "main_index.id=?;", \
                           (basename(self._handle.name), self._indexkey))
        record_keys = [key[0] for key in recordindex.description]
        record_index = recordindex.fetchone()
        con.commit()
        con.close()
        if not record_index:
            con.close()
            raise KeyError("_indexkey must match a stored record")
        record_indexdict = dict((record_keys[i], record_index[i]) for \
                                    i in range(len(record_keys)))
        del(record_indexdict["fileid"])
        self.__index = record_indexdict
        return record_indexdict                     


    def __set_and_save_index(self, indexdict):
        """Save a new index dict to the database
        
        (stub)
        """
        if isinstance(indexdict, dict):
            self.__index = indexdict
            #case, no indexdb
            if self._indexdb is None \
                or self._indexkey is not None:
                return
        else:
            raise ValueError("_index must be a dictionary")
        
        con = self.__connect_db_return_connection()
        #case: indexdb provided but it is empty
        if not self.__is_valid_db_with_tables(con):
            self.__create_tables(con, indexdict)
            con.commit()
        cursor = con.cursor()
        fileid = self.__get_fileid(con, write=True)
        
        #case: indexdb is provided but already contains
        # this specific record. This operation should raise
        cursor.execute("SELECT main_index.* " +\
                       "FROM main_index " +\
                       "INNER JOIN indexed_files " +\
                       "ON main_index.fileid = indexed_files.fileid " +\
                       "WHERE indexed_files.filename=? AND " +\
                       "main_index.recordoffsetstart=?;", \
                       (basename(self._handle.name), \
                       indexdict["recordoffsetstart"]))
        samefileposition = cursor.fetchone()
        cursor.execute("SELECT main_index.* " +\
                       "FROM main_index " +\
                       "INNER JOIN indexed_files " +\
                       "ON main_index.fileid = indexed_files.fileid " +\
                       "WHERE indexed_files.filename=? AND " +\
                       "main_index.id=?;", \
                       (basename(self._handle.name), \
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
        con.execute("INSERT INTO main_index ({0})".format(keys) +\
                       "VALUES ({0})".format(placehold), values)
        #increment the file record counter
        con.execute("UPDATE indexed_files SET count=count + 1 "
                    "WHERE filename=?", (basename(self._handle.name),))
        con.commit()
        con.close()
        #case: indexdb provided and has correct format/table
        #self.__insert_index("implementation pending")
    
    def __create_tables(self, con, indexdict):
        #create metadata table and populate format
        #set locking mode and synchronous
        con.execute("PRAGMA synchronous=0")
        con.execute("PRAGMA locking_mode=EXCLUSIVE")
        con.commit()
        
        con.execute( \
                "CREATE TABLE meta_data(key TEXT UNIQUE, value TEXT);")
        con.execute("INSERT INTO meta_data (key, value) VALUES (?,?);",
                      ("format", self._format))
        #create file table
        con.execute("CREATE TABLE indexed_files(fileid INTEGER PRIMARY " +\
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
        indextab = "CREATE TABLE main_index(" + \
                   "id TEXT, fileid INTEGER, "+ indextab + \
                   ", FOREIGN KEY(fileid) REFERENCES indexed_files(fileid));"
        con.execute(indextab)
        
        #create features table
        feattab = "CREATE TABLE features(" + \
                  "fileid INTEGER, " + \
                  "featurenumber INTEGER, " + \
                  "offsetbegin INTEGER, " + \
                  "offsetend INTEGER, " +\
                  "seqbegin INTEGER, " + \
                  "seqend INTEGER, " + \
                  "meta TEXT, " + \
                  "FOREIGN KEY(fileid) REFERENCES indexed_files(fileid));"
        con.execute(feattab)
        con.commit()
        
    def __connect_db_return_connection(self):
        """do some basic checking and return a connection to the db"""
        if not _sqlite:
            # Hack for Jython (or if Python is compiled without it)
            from Bio import MissingPythonDependencyError
            raise MissingPythonDependencyError("Requires sqlite3, which is "
                                               "included Python 2.5+")
        if not isfile(self._indexdb):
            raise ValueError("A valid database or None must be provided")
        #case: connection does not exist yet
        con = _sqlite.connect(self._indexdb)
        return con
    
    def __get_fileid(self, con, write=False):
        """Return the file id from the DB, write an entry if necessary
        
        returns fileid for a filename. If write is False a fileid not
        found in the database will return None while setting write to 
        True will require that the file is written prior to returning
        the fileid
        """
        if self.__file_id is not None:
            return self.__file_id 
        if self.__file_id is None:
            cursor = con.cursor()
            name = cursor.execute("SELECT fileid " +\
                           "FROM indexed_files WHERE " +\
                           "filename=?;", \
                           (basename(self._handle.name),))
            name = name.fetchone()
            if name is not None:
                self.__file_id = name[0]
                return name[0]
            elif write:
                cursor.execute("INSERT INTO indexed_files " +\
                               "(filename, count) VALUES (?, ?);",
                               (basename(self._handle.name), 0))
                name = cursor.execute("SELECT fileid " +\
                           "FROM indexed_files WHERE " +\
                           "filename=?;", \
                           (basename(self._handle.name),)).fetchone()
                self.__file_id = name[0]
                con.commit()
                return name[0]
            else:
                return None
            
                
    def __is_valid_db_with_tables(self, con):
        """if db is empty False, if correct tables True, else raise"""
        cursor = con.cursor()
        tempcursor = con.cursor()
        #case: database is empty
        tablenames = cursor.execute(
                    "SELECT name FROM sqlite_master WHERE " + \
                    "type='table' ORDER BY name;").fetchall()
        if len(tablenames) == 0:
            return False
        #case: database is valid and has correct format
        tablenames = [name[0] for name in tablenames]
        expectedtables = ['features', 'indexed_files', 'main_index', 'meta_data']
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

    _index = property(fget = __load_and_return_index,
                      fset = __set_and_save_index,
                      doc="the index may be loaded from a sqlite db")


    def __set_and_save_feature_index(self, feature_index_list):
        """Save a new index dict to the database

        (stub)
        """
        if isinstance(feature_index_list, list):
            self.__feature_index = feature_index_list
            #case, no indexdb
            if self._indexdb is None \
                or self._indexkey is not None:
                return
        else:
            raise ValueError("_feature_index must be a list")

        #case: indexdb provided but it is empty
        con = self.__connect_db_return_connection()
        if not self.__is_valid_db_with_tables(con):
            raise ValueError("_index must be set before _feature_index")

        cursor = con.cursor()
        fileid = self.__get_fileid(con, write=True)

        #--------------------- unedited below -------------------------------------------------------
        #case: indexdb is provided but already contains
        # this specific record. This operation should raise
        features_begin = self.__index.get("featuresoffsetstart")
        if features_begin:
            cursor.execute("SELECT features.* " +\
                           "FROM features " +\
                           "INNER JOIN indexed_files " +\
                           "ON features.fileid = indexed_files.fileid " +\
                           "INNER JOIN main_index " +\
                           "ON main_index.fileid = features.fileid " +\
                           "WHERE indexed_files.filename=? AND " +\
                           "main_index.featuresoffsetstart=?" +\
                           "LIMIT 2;", \
                           (basename(self._handle.name), \
                           features_begin))
            samefileposition = cursor.fetchone()
        else:
            raise ValueError("cannot save feature indexes without idx. begin")
        if samefileposition is not None:
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
        con.execute("INSERT INTO main_index ({0})".format(keys) +\
                       "VALUES ({0})".format(placehold), values)
        #increment the file record counter
        con.execute("UPDATE indexed_files SET count=count + 1 "
                    "WHERE filename=?", (basename(self._handle.name),))
        con.commit()
        con.close()
        #case: indexdb provided and has correct format/table
        #self.__insert_index("implementation pending")
    _feature_index = property(fget = __load_and_return_feature_index,
                              fset = __set_and_save_feature_index)

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

def lazy_iterator(handle, return_class=None, format = None,
                  alphabet=None, index=None):
    """A general iterator for sequential sequence files

    This function steps through a sequentially oriented sequence
    record file and returns SeqRecordProxy objects correspoding
    to the suported file formats.
    """
    if index is None:
        offsetiterator = sequential_record_offset_iterator(handle, format=format)
        for offset_begin, length in offsetiterator:
            yield return_class(handle, startoffset=offset_begin, \
                               length=length, indexdb = None, \
                               indexkey = None, alphabet=alphabet)
    else:
        #test that the db is populated by pulling the first record
        try:
            return_class(handle, indexdb = index, \
                         indexkey = 0, alphabet=alphabet)
        except KeyError:
            _make_index_db(handle, return_class =return_class,
                            indexdb = index, format=format,
                            alphabet=alphabet)

        #find the number of records
        if not _sqlite:
            # Hack for Jython (or if Python is compiled without it)
            from Bio import MissingPythonDependencyError
            raise MissingPythonDependencyError("Requires sqlite3, which is "
                                               "included Python 2.5+")
        con = _sqlite.connect(index)
        count = con.execute("SELECT count " +\
                            "FROM indexed_files WHERE " +\
                            "filename=?;", \
                            (basename(handle.name),)).fetchone()[0]
        con.close()
        for idx in range(count):
            yield return_class(handle, indexdb = index, \
                                indexkey = idx, alphabet=alphabet)
        

def _make_index_db(handle, return_class, indexdb, format, alphabet=None):
    """(private) populate index db with file's index information"""
    #this is used privately, but ensure it was used correctly
    offsetiterator = sequential_record_offset_iterator(handle, format=format)
    for offset_begin, length in offsetiterator:
        temp = return_class(handle, startoffset=offset_begin, \
                           length=length, indexdb = indexdb, \
                           indexkey = None, alphabet=alphabet)


def sequential_record_offset_iterator(handle, format = None):
    """Steps through a sequence file and returns offset pairs

    This function will step through a sequentially oriented sequence
    record file and return offset tuples for each format compliant
    record.

    The yield format is a 2-tuple representing the following:
    (record_begin_offset, record_offset_length)

    record_begin_offset = the index of the first character
    record_offset_length = the length (bytes) of the record,
                           accounting for padding and/or EOF
    """
    
    marker = {"fasta": ">",
              "genbank": "LOCUS",
              "gb": "LOCUS"}[format]
    marker_re = re.compile(marker)

    # Set handle to beginning and set offsets to valid initial values
    handle.seek(0)
    current_offset, next_offset = -1, 0
    end_offset = None
    padding = 0
    while True:
        # Advance the current_offset markers, using handle.readline().
        # unlike iterators, reaching the end of the handle and calling
        # readline will not raise StopIteration, instead it will return
        # a an empty string and the file offset will not advance
        current_offset = next_offset
        currentline = _bytes_to_string(handle.readline())
        next_offset = handle.tell()

        # Check if we encounter the first record then save the offset
        if marker_re.match(currentline) and not end_offset:
            current_record_offset = current_offset
        # Encountering any record but the first will return the most recent
        # fully-read record
        elif marker_re.match(currentline) or not currentline:
            length = current_offset - padding - current_record_offset
            yield current_record_offset, length
            current_record_offset = current_offset
            if not currentline:
                break
            # One cannot assume the handle will retain position.
            handle.seek(next_offset)
        elif currentline.strip() == "":
            padding += len(currentline)
            end_offset = current_offset
        else:
            #save position if no marker is found in order to denote
            #the end of the feature when a new marker is found
            padding = 0
            end_offset = current_offset

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

    def __init__(self, length = None,
                 beginindex=0, endindex=1,
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
            default_powers = [23,26,29,32,35,38,41]
            for power in default_powers:
                if length <= 2**power:
                    self._set_max_bin_power(power)
                    self._dynamic_size = False
                    break
            if self._dynamic_size: #this should have been set to False
                error_string = "Sequence length is {0}: must be less than 2^41".format(length)
                raise ValueError(error_string)
        
    def _increase_bin_sizes(self):
        """increase max bin size 8x (2**3) and re-organize existing binned data
        
        In order to increase the total maximum bin size, the lowest set
        of bins must be merged up one level (first step) then the entire set
        must be moved down one level without disturbing the organization scheme.
        
        An assertion in this routine blocks sequences larger than 2**41 from
        being created.
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
            level = int( floor( log((7*k + 1),2)/3.0 ) )
            new_level = level + 1
            oL = int((2**(3*level) - 1)/7.0)
            new_oL = int((2**(3*new_level) - 1)/7.0)
            
            new_index = k - oL + new_oL 
            self._bins[new_index] = self._bins[k]
            self._bins[k] = []
               
    def _set_max_bin_power(self, power):
        """sets the maximum bin power and fixes other necessary attributes"""
        
        self._max_bin_power = power
        self._min_bin_power = self._max_bin_power - 3*(self._bin_level_count-1)
        self._size_list = [2**(self._min_bin_power+3*n) for n in range(self._bin_level_count)]
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
        return sum(len(bin) for bin in self._bins) 

    def sort(self):
        """this performs bin-centric sorting, necessary for faster retrieval"""
        #bins must be sorted by the begin index, this is fastest
        if self._beginindex == 0:
            for i in range(len(self._bins)):
                self._bins[i].sort()
        #this is a bit slower but accomodates diverse data structures
        else:
            for i in range(len(self._bins)):
                self.bins[i].sort(key = lambda tup: tup[self._beginindex])
        #reset sorted quality
        self._sorted = True
            
    def __getitem__(self, key):
        """This getter efficiently retrieves the required entries
        
        This getter primarily works as expected and in a pythonic
        fashion, one exception is it's treatment of slices indices
        where the start is greater than the stop. Rather than just 
        throwing calculated output, an IndexError is raised.
        """    
        #set some locals
        beginindex = self._beginindex
        endindex = self._endindex
        
        #check that it is a slice and it has no step property (or step==1)
        if not isinstance(key, slice):
            if _is_int_or_long(key):
                key = slice(key, key+1)
            else:
                raise TypeError("lookups in the feature bin must use slice or int keys")
        
        #any integers are just converted to a 'len() == 1' slice
        if key.step is not None and key.step != 1:
            raise KeyError("lookups in the feature bin may not use slice stepping ex. bins[0:50:2]")
        
        #fix begin or end index for slicing of forms: bins[50:] or bins[:50] or even bins[:]
        keystart, keystop, keystep = key.indices(self._max_sequence_length)
        if keystart > keystop:
            raise IndexError("key not valid, slice.start > slice.stop")
        
        #pre-sort if necessary
        if not self._sorted:
            self.sort()
            
        #code taken from self._calculate_bin_index(), comments removed
        return_entries = []
        possible_entries = []
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
                for bin in range(k1+1, k2-1):
                    return_entries.extend( self._bins[bin]) 
            for binn in set([k1,k2-1]):
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

    def _calculate_bin_index(self, begin,span):
        """ This function returns a bin index given a (begin, span) interval
        
        The equations for determination of bin index derived from this publication: 
           Journal:  Bioinformatics Vol 27 no. 5 2011, pages 718-719
           Article:  "Tabix: fast retrieval of sequence features from generic
                      tab delimited files"
           Author:   Heng Li
        
        This function should only be used privately in the context of having no
        easier relationship to assign bin index. Placing this in a loop for any
        task other than arbitrary assignments is a bad idea since many of the
        common tasks can provide bin index through other relationships.
        _increase_bin_sizes is an example of a routine that would suffer
        performance penalties were it to use this routine yet can be run
        efficiently using alternate bin index relationships.
        """
        
        #take care base cases with the span parameter
        assert span >= 0
        #this is required for determination of bin location for zero length seq's
        span = max(1, span)
        
        # fix bin sizes if needed. Also, if the bin size is
        # larger than expected, do some self-consistency checks
        while begin+span > self._max_sequence_length:
            if self._dynamic_size:
                assert begin+span <= 2**41   # len(seq) > 2.19 trillion is not reasonable
            elif not self._dynamic_size and begin+span > 2**self._max_bin_power:
                error_string = "feature index at {0}: must be less than 2^{1}".format \
                                                (begin+span, self._max_bin_power)
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