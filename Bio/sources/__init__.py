import time
import urllib
import operator

from Bio import StdHandler
from Bio.Tools.MultiProc.copen import copen_fn
from Bio.ReseekFile import ReseekFile
from Bio.WWW import RequestLimiter

from Martel import Parser

class Source:
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
        return handle
    def set(self, key, handle):
        self._rawset(key, handle)
    def _rawget(self, params):
        raise NotImplementedError, "Please implement in a derived class."
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

##class PythonFunction:
##    function name

##class Application:
##    location

##class BerkeleyDBIndex:
##    index_filename

##class FlatIndex:
##    filename

##class CORBA:
##    pass

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
        
