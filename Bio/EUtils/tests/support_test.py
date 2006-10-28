import urlparse, cgi, urllib2, sys
import cPickle as pickle
import cStringIO as StringIO
from Bio.EUtils import ThinClient

def _assert_same_urls(url1, url2):
    # The URL may have key/value parameters after the '?'
    # The order of these parameters is not defined
    split1 = urlparse.urlsplit(url1)
    split2 = urlparse.urlsplit(url2)
    assert split1[:3] == split2[:3], (split1[:3], split2[:3])
    d1 = dict(cgi.parse_qsl(split1[3]))
    d2 = dict(cgi.parse_qsl(split2[3]))
    for d in (d1, d2):
        if "tool" in d: del d["tool"]
        if "email" in d: del d["email"]
    assert d1 == d2, (d1, d2)
    

def _getname(name, prefix):
    n = 2
    if name is not None:
        if isinstance(name, int):
            assert name > 0
            n = n + name
        else:
            return name

    f = sys._getframe()
    for i in range(n):
        f = getattr(f, "f_back")
    name = f.f_code.co_name
    if prefix is not None:
        name = prefix + "." + name
    return name

class OpenFromPickleStore:
    def __init__(self, picklestring):
        infile = StringIO.StringIO(picklestring)
        self.unpickler = pickle.Unpickler(infile)
    def open(self, url, query = None):
        saved_url, saved_query = self.unpickler.load()
        _assert_same_urls(saved_url, url)
        assert saved_query == query

        return StringIO.StringIO(self.unpickler.load())

class UsePickleStore:
    def __init__(self, filename):
        self.data = pickle.load(open(filename))
    def client(self, name = None, prefix = None):
        name = _getname(name, prefix)
        
        s = self.data[name]
        return ThinClient.ThinClient(
            opener = OpenFromPickleStore(s))
    def done(self):
        pass

class NoPickleStore:
    def __init__(self, filename = None):
        pass
    def client(self, name = None, prefix = None):
        return ThinClient.ThinClient()
    def done(self):
        pass

class AddToPickleStore:
    def __init__(self):
        self.outfile = StringIO.StringIO()
        self.pickler = pickle.Pickler(self.outfile, 1)
        self.opener = urllib2.build_opener()
        
    def open(self, url, query = None):
        self.pickler.dump( (url, query) )
        s = self.opener.open(url, query).read()
        self.pickler.dump(s)
        return StringIO.StringIO(s)
    def getvalue(self):
        return self.outfile.getvalue()
    
class CreatePickleStore:
    def __init__(self, filename):
        self.filename = filename
        self.data = {}
    def client(self, name = None, prefix = None):
        name = _getname(name, prefix)
        if name in self.data:
            raise TypeError("multiple %r" % (name,))
        self.data[name] = AddToPickleStore()
        return ThinClient.ThinClient(opener = self.data[name])
    
    def done(self):
        outfile = open(self.filename, "wb")
        d = {}
        for k, v in self.data.items():
            d[k] = v.getvalue()
        pickle.dump(d, outfile, 1)
