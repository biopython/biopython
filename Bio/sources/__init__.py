"""Base code for representing a database source of sequences.

This module contains code which can be used to describe a specialized
location where sequence information is available from. The code here
covers the basic cases (ie. a CGI script to retrieve information from)
as classes (CGI, in our example) which need to be supplied the database
specific information. So, if you want to programmatically describe
a database, this is the place to start looking for useful classes.
"""
import time
import urllib
import operator
import StringIO

from Bio import StdHandler
from Bio.Tools.MultiProc.copen import copen_fn
from Bio.ReseekFile import ReseekFile
from Bio.WWW import RequestLimiter
from Bio import SeqRecord

from Martel import Parser

class Source:
    """Base class for representing a database source.

    This is the class to inherit from when you want to create
    a general class to represent a specific type of database source
    (ie. CGI scripts). When deriving, you must implement the
    _rawget() function which takes a list of parameters and returns
    results. These results might be a handle or any other kind of
    result. You can optionally derive from _post_process() to do
    any sort of clean-up on the data.
    
    (XXX This should return a SeqFeature
    object, me thinks, and only a handle if you don't know how to
    parse the results that come back. Or something like that.)
    """
    def __init__(self, name, delay, timeout, doc, failure_cases):
        self.name = name
        self.delay = delay
        self.timeout = timeout
        self.doc = doc
        if self.delay is not None:
            self.limiter = RequestLimiter(delay)
        self.failure_cases = failure_cases or []
    def __call_rawget(self, params):
        return self._rawget(params)
    def __do_post_processing(self, results):
        return self._post_process(results)
    def __serialize_rawget(self, params):
        handle = self._rawget(params)
        return handle.read()
    def __call_rawget_forked(self, params):
        end_time = time.time() + self.timeout
        handle = copen_fn(self.__serialize_rawget, params)
        while time.time() < end_time:
            if handle.poll():
                break
            time.sleep(0.01)
        else:
            handle.close()
            raise IOError, "timed out"
        return handle
    def __check_for_errors(self, handle, more_failure_cases):
        # don't put it inside a RessekFile handle unless needed
        if len(self.failure_cases + more_failure_cases) == 0:
            return handle
        handle = ReseekFile(handle)
        pos = handle.tell()
        # XXX I can optimize this by precompiling the expressions.
        for expression, errormsg in self.failure_cases + more_failure_cases:
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
                raise IOError, errormsg
        handle.seek(pos)
        return handle
    def get(self, params, failure_cases=[]):
        if hasattr(self, "limiter"):
            self.limiter.wait()
        if self.timeout is None:
            handle = self.__call_rawget(params)
        else:
            handle = self.__call_rawget_forked(params)
        handle = self.__check_for_errors(handle, failure_cases)
        results = self.__do_post_processing(handle)
        return results
    def set(self, key, handle):
        self._rawset(key, handle)
    def _rawget(self, params):
        raise NotImplementedError, "Please implement in a derived class."
    def _post_process(self, results):
        """Can be hooked into by a derived class. Does nothing by default.
        """
        return results
    def _rawset(self, key, handle):
        raise NotImplementedError, "Caching not supported here."

class CGI(Source):
    # cgi
    # url
    # getmethod
    def __init__(self, name=None, delay=None, timeout=None, doc=None,
                 failure_cases=None, cgi=None, url=None, getmethod=1):
        Source.__init__(self, name, delay, timeout, doc, failure_cases)
        self.cgi = cgi
        self.url = url
        self.getmethod = getmethod
    def _rawget(self, params):
        options = _my_urlencode(params)
        if self.getmethod:
            fullcgi = self.cgi
            if options:
                fullcgi = "%s?%s" % (self.cgi, options)
            handle = urllib.urlopen(fullcgi)
        else:    # do a POST
            handle = urllib.urlopen(self.cgi, options)
        return handle

    def _post_process(self, handle):
        """Try to parse the results of a query into a SeqRecord.

        This uses Martel auto-detection to try to determine the
        format of a file. If the format cannot be determined, the
        handle itself is returned.
        """
        try:
            records = SeqRecord.io.read(handle)
            # don't do this check now; not really a list, an Iterator
            # assert len(records) == 1, "Expected one record, got %s" \
            #                           % (str(records))
            return records[0]
        except TypeError: # could not determine file type
            return handle

class BioSQL(Source):
    """Represent a BioSQL-style database to retrieve SeqRecord objects.

    XXX This returns a SeqRecord-like object from get() instead of
    a handle (since BioSQL is not going to give you a handle). This
    should be standardized.
    """
    def __init__(self, name, doc = "", db_host = 'localhost',
                 db_user = 'root', db_passwd = '', sql_db = '',
                 namespace_db = ''):
        """Intialize with information for connecting to the BioSQL db.
        """
        Source.__init__(self, name, None, None, doc, None)
        self.db_host = db_host
        self.db_user = db_user
        self.db_passwd = db_passwd
        self.sql_db = sql_db
        self.namespace_db = namespace_db
    
    def _rawget(self, params):
        # do the import here to prevent circular import problems
        from BioSQL import BioSeqDatabase

        # for params, we expect to get something like
        # [('accession', 'AB030760')]. We don't worry about what the id
        # is called right now, and just try to find it in the database
        # any way we can
        assert len(params) == 1, "Expected only one parameter, got %s" % \
                                 (str(params))
        assert len(params[0]) == 2, "Expected two item parameter got %s" % \
                                    (str(params[0]))
        find_id = params[0][1]
        server = BioSeqDatabase.open_database(user = self.db_user,
                   passwd = self.db_passwd, host = self.db_host,
                   db = self.sql_db)
        db = server[self.namespace_db]
        # try our different id choices to test the query
        item = None
        for possible_id_type in ["accession", "display_id"]:
            try:
                item = apply(db.lookup, (), {possible_id_type : find_id})
            except IndexError:
                pass
        if item is None:
            raise KeyError("Could not get item with id: %s" % find_id)
        
        return item

class BioCorba(Source):
    """Represent a BioCorba BioSequenceCollection for SeqRecord objects.

    XXX This has the same SeqRecord-returning style as BioSQL.
    """
    def __init__(self, name, doc = "", ior_ref = None, server_type = None):
        """Intialize with IOR reference for a BioCorba Collection.
        
        ior_ref is a URL or file reference to an IOR string. The IOR
        should reference a BioSequenceCollection. This is the top level
        BioCorba object we should use for making objects available.

        server_type is a hack parameter which might be necessary if there
        are server/client issues (ie. as with Perl ORBit) that we need
        to muck around with. If not set, we just use a standard retriever.
        """
        Source.__init__(self, name, None, None, doc, None)
        retriever = self._get_retriever(server_type)
        self.corba_dict = self._get_corba_client(ior_ref, retriever)

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

    def _rawget(self, params):
        # for params, we expect to get something like
        # [('accession', 'AB030760')]. We don't worry about what the id
        # is called right now, and just try to find it in the database
        # any way we can
        assert len(params) == 1, "Expected only one parameter, got %s" % \
                                 (str(params))
        assert len(params[0]) == 2, "Expected two item parameter got %s" % \
                                    (str(params[0]))
        find_id = params[0][1]

        return self.corba_dict[find_id]

class IndexedFile(Source):
    """Return SeqRecord objects from an indexed file.

    This module deals with both flat file and BerkeleyDB indexes.
    These indexed files can be created by any of the compliant indexing
    implementations from Biopython, BioPerl, BioJava, etc...
    """
    def __init__(self, name, doc = "", dbname = ""):
        """Intialize with information about loading the database.

        dbname is the name of the database to open. This will likely
        be a filesystem path to a database directory.
        """
        Source.__init__(self, name, None, None, doc, None)
        self.db = self._load_database(dbname)

    def _load_database(self, name):
        """Get a connection with the given database.
        """
        from Bio import Mindy
        db = Mindy.open(dbname = name)
        return db

    def _rawget(self, params):
        """Do the database retrieval of the sequence, returning a handle.
        """
        assert len(params) == 1, "Expected one parameter, got %s" % \
                                 str(params)
        namespace, key = params[0]
        location = apply(self.db.lookup, (), {namespace : key})
        assert len(location) == 1, "Got multiple hits: %s" % location
        return StringIO.StringIO(location[0].text)

    def _post_process(self, handle):
        """Try to parse the results of a query into a SeqRecord.

        This uses Martel auto-detection to try to determine the
        format of a file. If the format cannot be determined, the
        handle itself is returned.
        """
        try:
            records = SeqRecord.io.read(handle)
            # don't do this check now; not really a list, an Iterator
            # assert len(records) == 1, "Expected one record, got %s" \
            #                           % (str(records))
            return records[0]
        except TypeError: # could not determine file type
            return handle

##class PythonFunction:
##    function name

##class Application:
##    location

def _my_urlencode(params):
    # urllib only handles key=value pairs.  However, some CGI
    # scripts also contain parameters that are passed without the
    # key= part.  Thus, search through the params for empty
    # strings (or None), and handle these myself.

    # params could be a dictionary of key->value or a list of
    # (key,value) pairs.  If it's a dictionary, convert it to a list.
    if operator.isMappingType(params):
        params = params.items()

    paramlist = []
    for x in params:
        if x[0]:
            paramlist.append(urllib.urlencode([x]))
        else:
            paramlist.append(urllib.quote_plus(x[1]))
    return '&'.join(paramlist)
        
