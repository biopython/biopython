# Copyright 2002 by Jeffrey Chang, Brad Chapman.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# The SQL and Corba was modified from an original implementation by
# Brad Chapman.

"""Implements Registry to access databases.  These objects access
databases using a dictionary-like interface, where the key is the ID
of the thing to look up, and the value returned is the data associated
with the key.

Classes:
DBRegistry     Accesses databases with a dictionary-like interface.
DBObject       Base class for Registry objects for databases.
DBGroup        Groups DBObjects.

CGIDB          Accesses CGI databases.
EUtilsDB       Accesses NCBI using EUtils.
BioSQLDB       Accesses a BioSQL database.
BioCorbaDB     Accesses a BioCorba database.
IndexedFileDB  Accesses a Mindy Indexed file.
"""
from Bio.config.Registry import *

class DBRegistry(Registry):
    """This implements a dictionary-like interface to databases.

    """
    def __init__(self, name, load_path=None):
        Registry.__init__(self, name, load_path=load_path)

# Create a registry for access to databases.
db = DBRegistry("db", "Bio.dbdefs")

def _clean_abbrev(abbrev):
    return abbrev.replace("-", "_")

class DBObject(RegisterableObject):
    """This is a base class for dictionary-like interfaces to
    databases.

    Methods:
    get                  Lookup a key in a database, with a default value.
    get_as               Lookup a key and convert to an object.
    __getitem__          Lookup a key in a database.
    
        THE FOLLOWING SHOULD BE IMPLEMENTED IN A DERIVED CLASS.
    _get                 Return the data indicated by key.
    _convert_to          Convert the data to another object.
        IMPLEMENT THESE ONLY IF TIMEOUT OR CONCURRENT ACCESS IS NEEDED.
    _make_pickleable     Make the object returned by _get to a pickleable.
    _unmake_pickleable   Turn the pickleable object back into the original

    """
    def __init__(self, name, abbrev=None, doc=None, delay=None):
        """DBObject(name[, abbrev][, doc][, delay])"""
        import _support
        abbrev = _clean_abbrev(abbrev or name)
        RegisterableObject.__init__(self, name, abbrev, doc)
        if delay is not None:
            x = _support.make_rate_limited_function(self._get, delay)
            setattr(self, "_get", x)

    def set(self, key, data):
        self._set(key, data)

    def get(self, key, default=None):
        """S.get(key[, default]) -> data"""
        try:
            results = self[key]
        except KeyError:
            results = default
        return results

    def get_as(self, key, to_io=None, default=None):
        """S.get_as(key[, to_io][, default]) -> object"""
        data = self.get(key, default=default)
        return self._convert_to(data, to_io)

    def __getitem__(self, key):
        try:
            return self._get(key)
        except IOError, x:
            if str(x) == "timed out":
                raise KeyError, x
            raise

    # THESE FUNCTIONS CAN BE OVERLOADED IN A DERIVED CLASS.
        
    def _get(self, key):
        """S._get(key) -> data"""
        # Look up a key in the DB and return the data.
        raise NotImplementedError, "Please implement in a derived class."
    def _convert_to(self, data, to_io):
        """S._convert_to(data, to_io) -> another data type"""
        # Convert the data returned by _get to the type specified by
        # to_io, which is a FormatIO object.
    def _set(self, key, data):
        """S._set(key, data)"""
        # Not used.  May be used in the future to support caching.
        raise NotImplementedError, "Caching not supported here."
    def _make_pickleable(self, data):
        """S._make_pickleable(key, data) -> pickleable_obj"""
        # Make the handle a pickle-able python object.
        # Only need to implement if supporting timeout or concurrent
        # access.
        raise NotImplementedError, "pickling not supported."
    def _unmake_pickleable(self, pickleable_obj):
        """S._unmake_pickleable(key, pickleable_obj) -> data"""
        # Turn the pickle-able python object back into a handle.
        # Only need to implement if supporting timeout or concurrent
        # access.
        raise NotImplementedError, "pickling not supported."

class DBGroup(RegisterableGroup):
    """Groups DBObjects that return the same kind of data.

    """
    def __init__(self, name, abbrev=None, doc=None, cache=None):
        """DBGroup(name[, abbrev][, doc])

        name is the name of the object, and abbrev is an abbreviation
        for the name.
        """
        abbrev = _clean_abbrev(abbrev or name)
        RegisterableGroup.__init__(self, name, abbrev, doc)
        self._last_object_used = None

    def __getitem__(self, key):
        for obj in self.objs:
            try:
                handle = obj[key]
            except SystemError, KeyboardInterrupt:
                raise
            except Exception, x:
                continue
            else:
                self._last_object_used = obj
                return handle
        raise KeyError, "I could not get any results."

    def get(self, key, default=None):
        try:
            data = self[key]
        except KeyError:
            data = default
        return data

    def get_as(self, key, to_io=None, default=None):
        """S.get_as(key[, to_io][, default]) -> object"""
        data = self.get(key, default=default)
        return self._last_object_used._convert_to(data, to_io)

class TextLikeMixin:
    """Mixin class with useful functionality for retrival of text files.

    This implements some useful helper functions and overrides of DBObject
    for those implementations which need to retrieve text, check for errors in
    the retrieve text, and then convert that text to other formats.
    """
    def _check_for_errors(self, handle, failure_cases):
        from Martel import Parser
        from Bio import StdHandler
        from Bio.EUtils.ReseekFile import ReseekFile
        
        if not failure_cases:
            return handle
        handle = ReseekFile(handle)
        pos = handle.tell()
        for expression, errormsg in failure_cases:
            handle.seek(pos)
            parser = expression.make_parser()
            handler = StdHandler.RecognizeHandler()
            parser.setContentHandler(handler)
            parser.setErrorHandler(handler)
            try:
                parser.parseFile(handle)
            except Parser.ParserException:
                pass
            if handler.recognized:
                raise KeyError, errormsg
        handle.seek(pos)
        return handle

    def _convert_to(self, handle, to_io):
        from Bio import FormatIO
        x = to_io.read(handle)
        if isinstance(x, FormatIO.FormatIOIterator):
            i = 0
            for rec in x:
                if i > 0:
                    raise AssertionError, "Multiple records returned"
                i += 1
        else:
            rec = x
        return rec

class CGIDB(DBObject, TextLikeMixin):
    """This class implements DBObject for accessing CGI databases.

    """
    def __init__(self, name, cgi, url=None, key=None, params=None, 
                 abbrev=None, doc=None, delay=None, timeout=None,
                 getmethod=1, failure_cases=None):
        """CGIDB(name, cgi[, url][, key][, params][, abbrev][, doc]
        [, delay][, timeout][, getmethod][, failure_cases])

        name is the name of the object, abbrev is an abbreviation for
        the name, and doc is some documentation describing the object.
        
        cgi is the URL for the cgi script.  url points to the
        human-readable URL of the form.

        params is a list of (key, value) tuples indicating the
        parameters that should be passed to the CGI script.  key is
        the name of the parameter for the CGI script whose value is
        the ID of the object to retrieve.

        getmethod is a boolean describing whether a GET or POST should
        be used.  By default, GET is used.

        failure_cases is a list of (Martel Expression, error message)
        describing patterns of errors in the text returned by the
        script.

        """
        import _support
        DBObject.__init__(self, name=name, abbrev=abbrev,
                          doc=doc, delay=delay, timeout=timeout)
        self.cgi = cgi
        self.key = key or ''
        self.params = params or []
        self.url = url
        self.getmethod = getmethod
        self.failure_cases = []
        for exp, message in failure_cases or []:
            exp = _support.make_cached_expression(exp)
            self.failure_cases.append((exp, message))

    def _normalize_params(self, key):
        return self.params + [(self.key, key)]
    
    def _get(self, key):
        handle = self._cgiopen(key)
        handle = self._check_for_errors(handle, self.failure_cases)
        return handle

    def _cgiopen(self, key):
        import urllib
        params = self._normalize_params(key)
        options = _my_urlencode(params)
        if self.getmethod:
            fullcgi = self.cgi
            if options:
                fullcgi = "%s?%s" % (self.cgi, options)
            handle = urllib.urlopen(fullcgi)
        else:    # do a POST
            handle = urllib.urlopen(self.cgi, options)
        return handle

    def _make_pickleable(self, handle):
        return handle.read()

    def _unmake_pickleable(self, obj):
        import StringIO
        return StringIO.StringIO(obj)

class EUtilsDB(DBObject, TextLikeMixin):
    """Implement DBObject for accessing EUtils databases at NCBI.
    """
    def __init__(self, name, db, rettype, abbrev = None, doc = None,
                 failure_cases = None, delay = None, timeout = None):
        """Initialize an EUtilsDB connection for retrieval.

        name is the name of the object, abbrev is an abbreviation for
        the name, and doc is some documentation describing the object.

        db is the name of the database at NCBI you want to retrieve from
        (ie. protein, nucleotide, pubmed)

        rettype is the type of information to return
        (ie. gp, gb, fasta, medline)

        failure_cases is a list of (Martel Expression, error message)
        describing patterns of errors in the text returned by the
        script.
        """
        import _support
        DBObject.__init__(self, name=name, abbrev=abbrev,
                          doc=doc, delay=delay, timeout=timeout)
        self.db = db
        self.rettype = rettype
        self.failure_cases = []
        for exp, message in failure_cases or []:
            exp = _support.make_cached_expression(exp)
            self.failure_cases.append((exp, message))

    def _get(self, key):
        """Implementation of retrieval -- used DBIds client from EUtils.
        """
        from Bio.EUtils import DBIds
        from Bio.EUtils import DBIdsClient
        db_id = DBIds(self.db, [key])
        eutils_client = DBIdsClient.from_dbids(db_id)
        handle = eutils_client.efetch(retmode = "text", rettype =
                self.rettype)
        handle = self._check_for_errors(handle, self.failure_cases)
        return handle

class BioSQLDB(DBObject):
    """Represent a BioSQL-style database to retrieve SeqRecord objects.

    This returns a SeqRecord-like object from _get() instead of a
    handle (since BioSQL is not going to give you a handle).
    
    """
    def __init__(self, name, doc = "", db_host = 'localhost', db_port = '',
                 db_user = 'root', db_passwd = '', sql_db = '',
                 namespace_db = '', db_type = 'mysql'):
        """Intialize with information for connecting to the BioSQL db.
        """
        DBObject.__init__(self, name=name, doc=doc)
        self.db_host = db_host
        self.db_port = db_port
        self.db_user = db_user
        self.db_passwd = db_passwd
        self.sql_db = sql_db
        self.namespace_db = namespace_db
        self.db_type = db_type

    def _get_db_module(self, db_type):
        """Retrieve the appropriate module to use for connecting to a database

        This parses a description of the database and tries to determine
        which module is appropriate for that database type.
        """
        if db_type in ['mysql']:
            return 'MySQLdb'
        elif db_type in ['pg', 'postgres', 'postgresql']:
            raise ValueError("Postgres not supported yet. Sorry.")
        else:
            raise ValueError("Unknown database type: %s" % db_type)
   
    def _get(self, key):
        # do the import here to prevent circular import problems
        from BioSQL import BioSeqDatabase

        # for params, we expect to get something like
        # [('accession', 'AB030760')]. We don't worry about what the id
        # is called right now, and just try to find it in the database
        # any way we can
        find_id = key

        db_driver = self._get_db_module(self.db_type)
        open_args = {"user" : self.db_user,
                     "passwd" : self.db_passwd,
                     "host" : self.db_host,
                     "db" : self.sql_db,
                     "driver" : db_driver}
        if self.db_port:
            open_args["port"] = self.db_port
        server = BioSeqDatabase.open_database( *(), **open_args)
        db = server[self.namespace_db]
        # try our different id choices to test the query
        item = None
        for possible_id_type in ["accession", "display_id"]:
            try:
                item = db.lookup( *(), **{possible_id_type : find_id})
            except IndexError:
                pass
        if item is None:
            raise KeyError("Could not get item with id: %s" % find_id)
        return item
    
    def _convert_to(self, data, to_io):
        from Bio import SeqRecord
        if to_io != SeqRecord.io:
            raise ValueError, "format %s not supported" % to_io.name
        return data
    
    def _make_pickleable(self, item):
        return item
    def _unmake_pickleable(self, item):
        return item
    
class BioCorbaDB(DBObject):
    """Represent a BioCorba BioSequenceCollection for SeqRecord objects.

    Returns SeqRecord-like objects.
    
    """
    def __init__(self, name, ior_ref, server_type=None, doc=""):
        """Intialize with IOR reference for a BioCorba Collection.
        
        ior_ref is a URL or file reference to an IOR string. The IOR
        should reference a BioSequenceCollection. This is the top level
        BioCorba object we should use for making objects available.

        server_type is a hack parameter which might be necessary if there
        are server/client issues (ie. as with Perl ORBit) that we need
        to muck around with. If not set, we just use a standard retriever.
        """
        DBObject.__init__(self, name=name, doc=doc)
        self.retriever = self._get_retriever(server_type)
        self.ior_ref = ior_ref
        self.corba_dict = None

    def _get_retriever(self, server_type):
        """Return a BioCorba retriever object based on the specified server.

        This returns a ready-to-go client retriever which can be used to
        connect to a BioCorba server.
        """
        # do the BioCorba imports here, so we don't have to have it
        # installed to use this module
        from BioCorba.Client.BiocorbaConnect import PerlCorbaClient, \
          PythonCorbaClient, JavaCorbaClient, GenericCorbaClient
        from BioCorba.Client.Seqcore.CorbaCollection import \
          BioSequenceCollection

        if server_type is None:
            client_type = GenericCorbaClient
        else:
            server_type = server_type.lower()
            if server_type.find("python") >= 0:
                client_type = PythonCorbaClient
            elif server_type.find("java") >= 0:
                client_type = JavaCorbaClient
            elif server_type.find("perl") >= 0:
                client_type = PerlCorbaClient
            else:
                raise ValueError("Unexpected server type specified: %s" % 
                                 server_type)

        retriever = client_type(BioSequenceCollection)
        return retriever

    def _get_corba_client(self, ior_ref, retriever):
        """Get a connection to the CORBA server based on the ior_ref
        """
        # do the imports here so we don't need BioCorba for whole module
        from BioCorba.Bio import GenBank
        
        if ior_ref.find("http") >= 0: # assume it is a url
            client = retriever.from_url_ior(ior_ref)
        else: # assume it is a file
            client = retriever.from_file_ior(ior_ref)

        return GenBank.Dictionary(client, GenBank.FeatureParser())

    def _get(self, key):
        # get the corba dictionary only once when fetched
        if self.corba_dict is None:
            self.corba_dict = self._get_corba_client(self.ior_ref, 
                                                     self.retriever)
        return self.corba_dict[key]
    
    def _convert_to(self, data, to_io):
        from Bio import SeqRecord
        if to_io != SeqRecord.io:
            raise ValueError, "format %s not supported" % to_io.name
        return data

class IndexedFileDB(DBObject):
    """Return SeqRecord objects from an indexed file.

    This module deals with both flat file and BerkeleyDB indexes.
    These indexed files can be created by any of the compliant indexing
    implementations from Biopython, BioPerl, BioJava, etc...
    
    """
    def __init__(self, name, dbname, doc = ""):
        """Intialize with information about loading the database.

        dbname is the name of the database to open. This will likely
        be a filesystem path to a database directory.
        """
        DBObject.__init__(self, name=name, doc=doc)
        self.db = self._load_database(dbname)

    def _load_database(self, name):
        """Get a connection with the given database.
        """
        from Bio import Mindy
        db = Mindy.open(dbname = name)
        return db

    def _get_check_names(self, given_name, db):
        """Get a list of all namespaces to search for the file under.

        If given_name is a valid key, then it is returned as the only
        thing to check. Otherwise, we go forward and check all possible
        namespaces.
        """
        if given_name is not None and given_name in db.keys():
            return [given_name]
        else:
            return db.keys()

    def _get(self, key):
        """Do the database retrieval of the sequence, returning a handle.
        """
        # XXX jchang: how does this namespace/key stuff work?  can we
        # get rid of namespace?
        import operator
        import StringIO
        if not operator.isSequenceType(key) or len(key) != 2:
            raise ValueError, "Key should be tuple of (namespace, key)"
        namespace, key = key
        names_to_check = self._get_check_names(namespace, self.db)
        for check_name in names_to_check:
            location = self.db.lookup( *(), **{check_name : key})
            if len(location) >= 1:
                break
        assert len(location) == 1, "Got multiple hits: %s" % location
        return StringIO(location[0].text)

    def _convert_to(self, handle, to_io):
        from Bio import FormatIO
        x = to_io.read(handle)
        if isinstance(x, FormatIO.FormatIOIterator):
            i = 0
            for rec in x:
                if i > 0:
                    raise AssertionError, "Multiple records returned"
                i += 1
        else:
            rec = x
        return rec

def _my_urlencode(params):
    # urllib only handles key=value pairs.  However, some CGI
    # scripts also contain parameters that are passed without the
    # key= part.  Thus, search through the params for empty
    # strings (or None), and handle these myself.

    # params could be a dictionary of key->value or a list of
    # (key,value) pairs.  If it's a dictionary, convert it to a list.
    import operator
    import urllib

    if operator.isMappingType(params) and hasattr(params, "items"):
        params = params.items()

    paramlist = []
    for key, value in params:
        if key:
            paramlist.append(urllib.urlencode([(key, value)]))
        else:
            paramlist.append(urllib.quote_plus(value))
    return '&'.join(paramlist)
