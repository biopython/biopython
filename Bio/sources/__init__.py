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
            records = SeqRecord.io.read(f)
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
        print params # what will I get here?
        server = BioSeqDatabase.open_database(user = self.db_user,
                   passwd = self.db_passwd, host = self.db_host,
                   db = self.sql_db)
        db = server[self.namespace_db]
        # try our different id choices to test the query
        item = None
        for possible_id_type in ["accession", "display_id"]:
            try:
                item = apply(self.db.lookup, (), {possible_id_type : find_id})
            except IndexError:
                pass
        if item is None:
            raise KeyError("Could not get item with id: %s" % find_id)
        
        return item

class BioCorba(Source):
    """Represent a BioCorba BioSequenceCollection for SeqRecord objects.

    XXX This has the same SeqRecord-returning style as BioSQL.
    """
    def __init__(self, name, doc = "", ior_ref = None):
        """Intialize with IOR reference for a BioCorba Collection.
        
        ior_ref is a URL or file reference to an IOR string. The IOR
        should reference a BioSequenceCollection. This is the top level
        BioCorba object we should use for making objects available.
        """
        Source.__init__(self, name, None, None, doc, None)

##class PythonFunction:
##    function name

##class Application:
##    location

##class BerkeleyDBIndex:
##    index_filename

##class FlatIndex:
##    filename

##class SQL:
##    pass





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
        
