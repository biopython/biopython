"""Search and retreive information using the EUtils history.

EUtils has two major modes.  One uses history while the other uses
database identifiers.  This is a high-level interface for working with
the history.  You should use this module if you expect to work with
large or an unknown number of identifiers.

See DBIdsClient if you want to get information about a set of known
database identifiers.

>>> from Bio.EUtils import HistoryClient
>>> client = HistoryClient.HistoryClient()
>>> cancer = client.search("cancer")
>>> print len(cancer)
1458353
>>> 

That's quite a few hits.  Most people would like to see the first few
records then try to refine the search.

>>> print cancer[:5].efetch(retmode = "text", rettype = "docsum").read()

1:  Seow-Choen F.
Author's reply: Adjuvant therapy for rectal cancer cannot be based on the
results of other surgeons (Br J Surg 2002; 89: 946-947).
Br J Surg. 2003 Jan;90(1):121-122.
PMID: 12520589 [PubMed - as supplied by publisher]

2:  Mortensen N, Lindsey I.
Adjuvant therapy for rectal cancer cannot be based on the results of other
surgeons (Br J Surg 2002; 89: 946-947).
Br J Surg. 2003 Jan;90(1):121.
PMID: 12520588 [PubMed - as supplied by publisher]

3:  Osugi H, Takemura M, Higashino M, Takada N, Lee S, Kinoshita H.
A comparison of video-assisted thoracoscopic oesophagectomy and radical lymph
node dissection for squamous cell cancer of the oesophagus with open operation.
Br J Surg. 2003 Jan;90(1):108-13.
PMID: 12520585 [PubMed - in process]

4:  Tanaka M, Kitajima Y, Sato S, Miyazaki K.
Combined evaluation of mucin antigen and E-cadherin expression may help select
patients with gastric cancer suitable for minimally invasive therapy.
Br J Surg. 2003 Jan;90(1):95-101.
PMID: 12520583 [PubMed - in process]

5:  Diaz De Liano A, Oteiza Martinez F, Ciga MA, Aizcorbe M, Cobo F, Trujillo R.
Impact of surgical procedure for gastric cancer on quality of life.
Br J Surg. 2003 Jan;90(1):91-4.
PMID: 12520582 [PubMed - in process]

>>>

Now refine the query to publications in the last day

>>> from Bio import EUtils
>>> recent_cancer = client.search("#%s" % (cancer.query_key,),
...                               daterange = EUtils.WithinNDays(1))
>>> len(recent_cancer)
106
>>>

Still quite a few.  What's the last one about?
>>> for k, v in recent_cancer[-1].summary().dataitems.allitems():
...     print k, "=", v
...

PubDate = 2002/12/01
Source = Nippon Shokakibyo Gakkai Zasshi
Authors = Kuroki T
Title = [Strategy against cancer in 21 century, with emphasis of cancer prevention and refractory cancer]
Volume = 99
Pages = 1423-7
EntrezDate = 2003/01/10
PubMedId = 12518389
MedlineId = 22406828
Lang = Japanese
PubType =
RecordStatus = PubMed - in process
Issue = 12
SO = 2002 Dec;99(12):1423-7
DOI =
JTA = KJY
ISSN = 0446-6586
PubId =
PubStatus = 4
Status = 6
HasAbstract = 0
ArticleIds = {'MedlineUID': u'22406828', 'PubMedId': u'12518389'}
>>> 

Here's an interesting one.  Which articles are related to this one but
are not about cancer?  First, get the related articles.


>>> neighbors = recent_cancer[-1].neighbor_links()
>>> dbids = neighbors.linksetdbs["pubmed_pubmed"].dbids
>>> len(dbids)
10296
>>> 

Upload that back to the server

>>> related_result = client.post(dbids)
>>>
>>> non_cancer = client.search("#%s NOT #%s" % (related_result.query_key,
...                                             cancer.query_key))
>>> len(non_cancer)
4000
>>>

The HistoryClient instance has an attribute named 'query_history'
which stores the searches done so far, keyed by the query_key value
assigned by the server.  The history on the server can expire.  If
that is detected during a search then previous results are invalidated
and removed from the query_history.  Future requests from invalidated
results will raise an error.

If a request is made from a search which has not been invalidated but
whose history has expired then queries like 'summary' will raise an
error.  Some other request (like 'dbids') may return success but
contain undefined information.

""" #"
    
import types
import ThinClient, parse, Datatypes, Mixins, Config

class HistoryCookie:
    """Data needed to get back to the history"""
    def __init__(self, db, webenv_ref, query_key):
        self.db = db
        self.webenv_ref = webenv_ref
        self.query_key = query_key

class HistoryLookup(object):
    """Look up information about a search in history

    To get the list of dbids by fetching the server's "uilist",
    use the "dbids" attribute.
    """
    def __init__(self, eutils, cookie, retstart, retmax):
        self.eutils = eutils
        self.cookie = cookie
        self.retstart = retstart
        self.retmax = retmax
        self.db = cookie.db
        self.query_key = cookie.query_key

    def _check_invalid(self):
        # Check if we can get data from this history
        if self.cookie.query_key is None:
            raise NotImplementedError("empty data set")
        if self.query_key is None:
            raise Datatypes.EUtilsError(
                "query history no longer available on server")

    def esummary(self, retmode = 'xml', rettype = None):
        """Request the eSummary for this history; returns the socket handle"""
        self._check_invalid()
        infile = self.eutils.esummary_using_history(
            webenv = self.cookie.webenv_ref[0],
            db = self.cookie.db,
            query_key = self.cookie.query_key,
            retstart = self.retstart,
            retmax = self.retmax)
        return infile

    def summary(self):
        """the Datatypes.Summary for this history"""
        return parse.parse_summary_xml(self.esummary("xml"))

    def elink(self,
              db = "pubmed",
              cmd = "neighbor",
              term = None,
              field = None,
              daterange = None):
        """Request an eLink for this history; returns the socket handle"""
        self._check_invalid()
        return self.eutils.elink_using_history(
            webenv = self.cookie.webenv_ref[0],
            query_key = self.cookie.query_key,
            db = db,
            dbfrom = self.cookie.db,
            cmd = cmd,
            retstart = self.retstart,
            retmax = self.retmax,
            daterange = daterange,
            term = term,
            field = field,
            )

    def _get_dbids(self):
        infile = self.efetch(retmode = "text", rettype = "uilist")
        ids = parse.parse_fetch_identifiers(infile)
        return Datatypes.DBIds(self.cookie.db, ids)
    dbids = property(_get_dbids, None, None,
        "The DBIds for this results set, fetched from the server's 'uilist'")

class HistoryRecord(HistoryLookup):
    """Get information about a single record in a history"""
    def __init__(self, eutils, cookie, offset):
        HistoryLookup.__init__(self, eutils, cookie, offset, 1)
    def summary(self):
        """the Datatypes.Summary for this history record"""
        return HistoryLookup.summary(self)[0]

class SequenceHistoryFetchMixin:
    def efetch(self, retmode = 'xml', rettype = None,
               seq_start = None, seq_stop = None, strand = None,
               complexity = None):
        self._check_invalid()
        if strand not in (None, 1, 2):
            raise TypeError("Strand can only be 1 (plus, default) or 2 (minus)")
        return self.eutils.efetch_using_history(
            webenv = self.cookie.webenv_ref[0],
            db = self.cookie.db,
            query_key = self.cookie.query_key,
            retstart = self.retstart,
            retmax = self.retmax,
            retmode = retmode,
            rettype = rettype,
            seq_start = seq_start,
            seq_stop = seq_stop,
            strand = strand,
            complexity = complexity)

class SequenceHistoryRecord(Mixins.SequenceFetchMixin,
                            SequenceHistoryFetchMixin,
                            HistoryRecord):
    pass

class PublicationHistoryFetchMixin:
    def efetch(self, retmode = "xml", rettype = None):
        self._check_invalid()
        return self.eutils.efetch_using_history(
            webenv = self.cookie.webenv_ref[0],
            db = self.cookie.db,
            query_key = self.cookie.query_key,
            retstart = self.retstart,
            retmax = self.retmax,
            retmode = retmode,
            rettype = rettype)

class PublicationHistoryRecord(Mixins.PublicationFetchMixin,
                               PublicationHistoryFetchMixin,
                               HistoryRecord):
    pass

class BaseHistoryRecordSet(HistoryLookup):
    def __init__(self, eutils, cookie, retstart, retmax, metadata = None):
        HistoryLookup.__init__(self, eutils, cookie, retstart, retmax)
        self.metadata = metadata

    def __len__(self):
        return self.retmax
    
    def __getitem__(self, i):
        if isinstance(i, types.SliceType):
            if i.step is not None:
                raise TypeError("cannot set step size in slice")
            # Don't pass metadata downwards
            start = i.start
            if start is None: start = 0
            stop = i.stop
            if stop is None: stop = self.retmax
            # Potentially expensive, but this is the easy way
            # to get the semantics correct.
            x = range(self.retstart, self.retstart + self.retmax)[start:stop]
            if x:
                retstart = x[0]
                retmax = x[-1] - x[0] + 1
            else:
                retstart = 0
                retmax = 0
            return self.__class__(self.eutils, self.cookie, retstart,
                                  retmax)
        if 0 <= i < self.retmax:
            pos = self.retstart + i
        elif 1 <= -i <= self.retmax:
            pos = self.retstart + i + self.retmax
        else:
            raise IndexError(i)
        return self._record_class(self.eutils, self.cookie, pos)

class SequenceHistoryRecordSet(Mixins.SequenceFetchMixin,
                               SequenceHistoryFetchMixin,
                               BaseHistoryRecordSet):
    _record_class = SequenceHistoryRecord

class PublicationHistoryRecordSet(Mixins.PublicationFetchMixin,
                                  PublicationHistoryFetchMixin,
                                  BaseHistoryRecordSet):
    _record_class = PublicationHistoryRecord

def _get_recordset_constructor(db, dbtype):
    dbtype = Config.databases.gettype(db, dbtype)
    if dbtype == Config.SEQUENCE_TYPE:
        return SequenceHistoryRecordSet
    elif dbtype == Config.PUBLICATION_TYPE:
        return PublicationHistoryRecordSet
    else:
        raise TypeError("unknown database type: %r" % (dbtype,))
    return Hrec_set, Hrec

class HistoryClient:
    def __init__(self, eutils = None):
        if eutils is None:
            eutils = ThinClient.ThinClient()
        self.eutils = eutils
        self.webenv_ref = [None]
        self.query_history = {}

    def _check_for_cache_reset(self, query_key):
        # If this is a repeat then it's because the history
        # has expired.  Set the existing results to have a
        # query_key of None and reset the cache.
        if query_key not in self.query_history:
            # New value, so there was no reset
            return

        for v in self.query_history.values():
            v.query_key = None
        
        self.query_history.clear()

    def search(self,
               term,
               db = "pubmed",
               field = None,
               daterange = None,
               dbtype = None
               ):
        
        set_klass = _get_recordset_constructor(db, dbtype)
        
        infile = self.eutils.esearch(
            term = term,
            db = db,
            field = field,

            retstart = 0,
            retmax = 0,
            
            daterange = daterange,
            
            usehistory = 1,
            webenv = self.webenv_ref[0],
            )
        searchinfo = parse.parse_search(infile, self.webenv_ref)

        if searchinfo.query_key is not None:
            cookie = HistoryCookie(db, self.webenv_ref, searchinfo.query_key)
        else:
            assert searchinfo.count == 0
            cookie = HistoryCookie(db, None, None)

        recordset = set_klass(self.eutils, cookie, 0, searchinfo.count,
                              searchinfo)
        # won't have a query_key if the search turned up empty
        if searchinfo.query_key is not None:
            self._check_for_cache_reset(searchinfo.query_key)
            self.query_history[searchinfo.query_key] = recordset

        return recordset

    def post(self, dbids, dbtype = None):
        set_klass = _get_recordset_constructor(dbids.db, dbtype)

        infile = self.eutils.epost(dbids,
                                   webenv = self.webenv_ref[0])
        # Extract the webenv_ref since it may change.
        postinfo = parse.parse_post(infile, self.webenv_ref)

        # Were there any invalid identifiers?
        n = len(dbids) - len(postinfo.invalid_ids)

        cookie = HistoryCookie(dbids.db, self.webenv_ref,
                               postinfo.query_key)
        recordset = set_klass(self.eutils, cookie, 0, n, postinfo)
        self._check_for_cache_reset(postinfo.query_key)
        self.query_history[postinfo.query_key] = recordset
        return recordset

    from_dbids = post  # alias for similarity to DBIdsClient
