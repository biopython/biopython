"""Search and retrieve information given a set of database identifiers.

EUtils has two major modes.  One uses history while the other uses
database identifiers.  This is a high-level interface for working with
identifiers.  You should use this module to get information about a
set of known database identifiers.

See HistoryClient if you want to work with a large number of
identifiers or potentially large search results.

>>> from Bio import EUtils
>>> from Bio.EUtils import DBIdsClient
>>> client = DBIdsClient.DBIdsClient()
>>> result = client.search("dalke", retmax = 100)
>>> len(result)
30
>>> print result[0].efetch(retmode = "text", rettype = "abstract").read()

1: Pac Symp Biocomput  1997;:85-96

Using Tcl for molecular visualization and analysis.

Dalke A, Schulten K.

Beckman Institute, Urbana, IL 61801, USA.

Reading and manipulating molecular structure data is a standard task in every
molecular visualization and analysis program, but is rarely available in a form
readily accessible to the user. Instead, the development of new methods for
analysis, display, and interaction is often achieved by writing a new program,
rather than building on pre-existing software. We present the Tcl-based script
language used in our molecular modeling program, VMD, and show how it can access
information about the molecular structure, perform analysis, and graphically
display and animate the results. The commands are available to the user and make
VMD a useful environment for studying biomolecules.


PMID: 9390282 [PubMed - indexed for MEDLINE]

>>>


Find sequences similar to GI:4579714 which were published in 2002.

>>> protein = DBIdsClient.from_dbids(EUtils.DBIds("protein", "4579714"))
>>> neighbors = protein.neighbor_links("protein",
...        daterange = EUtils.DateRange("2002/01/01", "2002/12/31", "pdat"))
>>> dbids = neighbors.linksetdbs["protein_protein"].dbids
>>> len(dbids)
28
>>> print dbids
DBIds(u'protein', [u'4579714', u'25298947', u'24158913', u'24158914', u'24158915', u'17942993', u'17942994', u'17942995', u'20150921', u'20150922', u'20151159', u'25298949', u'19716034', u'20663737', u'20663738', u'20663741', u'24987328', u'25533128', u'25298946', u'25298948', u'23008597', u'20219020', u'21218340', u'21218344', u'19075395', u'21218338', u'21218342', u'21311795'])
>>> 
>>> print client.from_dbids(dbids[:5]).efetch(retmode="text",
...                                           rettype="summary").read()

1: BAA75200
Bacteriorhodopsin [Halobacterium sp.]
gi|4579714|dbj|BAA75200.1|[4579714]


2: H84300
bacteriorhodopsin [imported] - Halobacterium sp. NRC-1
gi|25298947|pir||H84300[25298947]


3: 1M0KA
Chain A, Bacteriorhodopsin K Intermediate At 1.43 A Resolution
gi|24158913|pdb|1M0K|A[24158913]


4: 1M0LA
Chain A, BacteriorhodopsinLIPID COMPLEX AT 1.47 A RESOLUTION
gi|24158914|pdb|1M0L|A[24158914]


5: 1M0MA
Chain A, Bacteriorhodopsin M1 Intermediate At 1.43 A Resolution
gi|24158915|pdb|1M0M|A[24158915]

>>>

"""

import types
import parse, Mixins, Config, ThinClient, Datatypes

class DBIdsLookup(object):
    """Look up information about a DBIds

    To get the list of dbids, as interpreted by fetching the
    server's "uilist", use the "dbids" attribute.
    """
    def __init__(self, eutils, records_dbids):
        self.eutils = eutils
        self.records_dbids = records_dbids

    def esummary(self, retmode = 'xml', rettype = None):
        """call esummary on this DBIds; returns the socket handle"""
        return self.eutils.esummary_using_dbids(
            dbids = self.records_dbids)

    def summary(self):
        """get the summary for these DBIds, parsed into a Datatypes.Summary"""
        return parse.parse_summary_xml(self.esummary("xml"))

    def elink(self,
              db = "pubmed",
              cmd = "neighbor",
              term = None,
              field = None,
              daterange = None):
        """call elink on this DBIds; returns the socket handle"""
        return self.eutils.elink_using_dbids(
            dbids = self.dbids,
            db = db,
            cmd = cmd,
            daterange = daterange,
            term = term,
            field = field,
            )

    def _get_dbids(self):
        infile = self.efetch(retmode = "text", rettype = "uilist")
        ids = parse.parse_fetch_identifiers(infile)
        return Datatypes.DBIds(self.records_dbids.db, ids)
    dbids = property(_get_dbids, None, None,
        "The DBIds for this results set, validated from the server's 'uilist'")
    
    
class DBIdsRecord(DBIdsLookup):
    """A single record on the server"""
    def summary(self):
        return DBIdsLookup.summary(self)[0]

class SequenceDBIdsFetchMixin:
    """Support 'efetch' for sequence records"""
    def efetch(self, retmode = 'xml', rettype = None,
               seq_start = None, seq_stop = None, strand = None,
               complexity = None):
        if strand not in (None, 1, 2):
            raise TypeError("Strand can only be 1 (plus, default) or 2 (minus)")
        return self.eutils.efetch_using_dbids(
            dbids = self.records_dbids,
            retmode = retmode,
            rettype = rettype,
            seq_start = seq_start,
            seq_stop = seq_stop,
            strand = strand,
            complexity = complexity)

class SequenceDBIdsRecord(Mixins.SequenceFetchMixin,
                          SequenceDBIdsFetchMixin,
                          DBIdsRecord):
    """a single sequence record, referenced by database identifier"""
    pass

class PublicationDBIdsFetchMixin:
    """Support 'efetch' for publication records"""
    def efetch(self, retmode = "xml", rettype = None):
        return self.eutils.efetch_using_dbids(
            dbids = self.records_dbids,
            retmode = retmode,
            rettype = rettype)

class PublicationDBIdsRecord(Mixins.PublicationFetchMixin,
                             PublicationDBIdsFetchMixin,
                             DBIdsRecord):
    """a single publication record, referenced by database identifier"""
    pass

class BaseDBIdsRecordSet(DBIdsLookup):
    """Base class for dealing with a set of records, reference by identifier"""
    def __init__(self, eutils, records_dbids, metadata = None):
        DBIdsLookup.__init__(self, eutils, records_dbids)
        self.metadata = metadata

    def __len__(self):
        """Number of records referenced by this RecordSet"""
        return len(self.records_dbids)

    def __getitem__(self, i):
        """Return subset of the records"""
        if isinstance(i, types.SliceType):
            # Metadata is not passed downwards
            if i.step is None:
                return self.__class__(
                    self.eutils,
                    self.records_dbids[i.start:i.stop])
            return self.__class__(
                self.eutils,
                self.records_dbids[i.start:i.stop:i.step])

        return self._record_class(self.eutils, self.records_dbids.item(i))
        
class SequenceDBIdsRecordSet(Mixins.SequenceFetchMixin,
                             SequenceDBIdsFetchMixin,
                             BaseDBIdsRecordSet):
    """a set of sequence records, referenced by database identifier"""
    _record_class = SequenceDBIdsRecord

class PublicationDBIdsRecordSet(Mixins.PublicationFetchMixin,
                                PublicationDBIdsFetchMixin,
                                BaseDBIdsRecordSet):
    """a set of publication records, referenced by database identifier"""
    _record_class = PublicationDBIdsRecord


def _get_recordset_constructor(db, dbtype):
    """get the right DataSet constructor for a database"""
    dbtype = Config.databases.gettype(db, dbtype)
    if dbtype == Config.SEQUENCE_TYPE:
        return SequenceDBIdsRecordSet
    elif dbtype == Config.PUBLICATION_TYPE:
        return PublicationDBIdsRecordSet
    else:
        raise TypeError("Unknown database type: %r" % (dbtype,))

def from_dbids(dbids, dbtype = None, eutils = None):
    """create a RecordSet interface for the set of database identifiers

    Parameters are:
      dbids -- a DBIds
      dbtype -- the dbtype to use (EUtils.Config.{SEQUENCE,PUBLIATION}_TYPE)
           in case dbids.db isn't in the list of know NCBI databases.
           Defaults to None.
      eutils -- the ThinClient to use, defaults to creating a new
           ThinClient.ThinClient()
    """
    return DBIdsClient(eutils).from_dbids(dbids, dbtype)

class DBIdsClient:
    """Create a RecordSet either from a search or a set of dbids

    The constructor takes an optional ThinClient to use for
    connecting to NCBI.
    """
    def __init__(self, eutils = None):
        if eutils is None:
            eutils = ThinClient.ThinClient()
        self.eutils = eutils

    def from_dbids(self, dbids, dbtype = None):
        """Return a RecordSet given the DBIds

        This RecordSet can be used to fetch data from NCBI
        related to the given DBIds.
        """
        set_klass = _get_recordset_constructor(dbids.db, dbtype)
        return set_klass(self.eutils, dbids, None)

    def search(self,
               term,
               db = "pubmed",
               field = None,

               retstart = 0,
               retmax = 20,

               daterange = None,
               dbtype = None,
               ):
        """do an Entrez search

        The parameters are:
          'term' -- the query string in the Entrez query language; see
             http://www.ncbi.nlm.nih.gov/entrez/query/static/help/pmhelp.html
          'db' -- the database to search

          'field' -- the field to use for unqualified words
                  Eg, "dalke[au] AND gene" with field==None becomes
                    dalke[au] AND (genes[MeSH Terms] OR gene[Text Word]
                  and "dalke[au] AND gene" with field=="au" becomes
                    dalke[au] AND genes[Author]
                 (Yes, I think the first "au" should be "Author" too)

          'retstart' -- include identifiers in the output, starting with
                   position 'retstart' (normally starts with 0)
          'retmax' -- return at most 'retmax' identifiers in the output
                   (if not specified, NCBI returns 20 identifiers)
          'daterange' -- a date restriction; either WithinNDays or DateRange
          
          'dbtype' -- (optional) the database type (Config.PUBLICATION_TYPE
                  or SEQUENCE_TYPE).  Overrides the type based on the 'db'
        """
        set_klass = _get_recordset_constructor(db, dbtype)
        infile = self.eutils.esearch(
            term = term,
            db = db,
            field = field,
            retstart = retstart,
            retmax = retmax,
            daterange = daterange)
        searchinfo = parse.parse_search(infile, [None])

        dbids = Datatypes.DBIds(db, searchinfo.ids)
        return set_klass(self.eutils, dbids, searchinfo)
