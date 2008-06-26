"""Low-level interface to NCBI's EUtils for Entrez search and retrieval.

For higher-level interfaces, see DBIdsClient (which works with a set
of database identifiers) and HistoryClient (which does a much better
job of handling history).

There are five classes of services:
  ESearch - search a database
  EPost - upload a list of indicies for further use
  ESummary - get document summaries for a given set of records
  EFetch - get the records translated to a given format
  ELink - find related records in other databases

You can find more information about them at
  http://www.ncbi.nlm.nih.gov/entrez/query/static/eutils_help.html
but that document isn't very useful.  Perhaps the following is better.

EUtils offers a structured way to query Entrez, get the results in
various formats, and get information about related documents.  The way
to start off is create an EUtils object.

>>> from Bio import EUtils
>>> from Bio.EUtils.ThinClient import ThinClient
>>> eutils = ThinClient.ThinClient()
>>> 

You can search Entrez with the "esearch" method.  This does a query on
the server, which generates a list of identifiers for records that
matched the query.  However, not all the identifiers are returned.
You can request only a subset of the matches (using the 'retstart' and
'retmax') terms.  This is useful because searches like 'cancer' can
have over 1.4 million matches.  Most people would rather change the
query or look at more details about the first few hits than wait to
download all the identifiers before doing anything else.

The esearch method, and indeed all these methods, returns a
'urllib.addinfourl' which is an HTTP socket connection that has
already parsed the HTTP header and is ready to read the data from the
server.

For example, here's a query and how to use it

  Search in PubMed for the term cancer for the entrez date from the
  last 60 days and retrieve the first 10 IDs and translations using
  the history parameter.

>>> infile = eutils.esearch("cancer",
...                         daterange = EUtils.WithinNDays(60, "edat"),
...                         retmax = 10)
>>>
>>> print infile.read()
<?xml version="1.0"?>
<!DOCTYPE eSearchResult PUBLIC "-//NLM//DTD eSearchResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eSearch_020511.dtd">
<eSearchResult>
        <Count>7228</Count>
        <RetMax>10</RetMax>
        <RetStart>0</RetStart>
        <IdList>
                <Id>12503096</Id>
                <Id>12503075</Id>
                <Id>12503073</Id>
                <Id>12503033</Id>
                <Id>12503030</Id>
                <Id>12503028</Id>
                <Id>12502932</Id>
                <Id>12502925</Id>
                <Id>12502881</Id>
                <Id>12502872</Id>
        </IdList>
        <TranslationSet>
                <Translation>
                        <From>cancer%5BAll+Fields%5D</From>
                        <To>(%22neoplasms%22%5BMeSH+Terms%5D+OR+cancer%5BText+Word%5D)</To>
                </Translation>
        </TranslationSet>
        <TranslationStack>
                <TermSet>
                        <Term>"neoplasms"[MeSH Terms]</Term>
                        <Field>MeSH Terms</Field>
                        <Count>1407151</Count>
                        <Explode>Y</Explode>
                </TermSet>
                <TermSet>
                        <Term>cancer[Text Word]</Term>
                        <Field>Text Word</Field>
                        <Count>382919</Count>
                        <Explode>Y</Explode>
                </TermSet>
                <OP>OR</OP>
                <TermSet>
                        <Term>2002/10/30[edat]</Term>
                        <Field>edat</Field>
                        <Count>-1</Count>
                        <Explode>Y</Explode>
                </TermSet>
                <TermSet>
                        <Term>2002/12/29[edat]</Term>
                        <Field>edat</Field>
                        <Count>-1</Count>
                        <Explode>Y</Explode>
                </TermSet>
                <OP>RANGE</OP>
                <OP>AND</OP>
        </TranslationStack>
</eSearchResult>

>>>

You get a raw XML input stream which you can process in many ways.
(The appropriate DTDs are included in the subdirectory "DTDs" and see
also the included POM reading code.)

    WARNING! As of this writing (2002/12/3) NCBI returns their
    XML encoded as Latin-1 but their processing instruction says
    it is UTF-8 because they leave out the "encoding" attribute.
    Until they fix it you will need to recode the input stream
    before processing it with XML tools, like this

        import codecs
        infile = codecs.EncodedFile(infile, "utf-8", "iso-8859-1")


The XML fields are mostly understandable:
  Count -- the total number of matches from this search
  RetMax -- the number of <ID> values returned in this subset
  RetStart -- the start position of this subset in the list of
      all matches

  IDList and ID -- the identifiers in this subset

  TranslationSet / Translation -- if the search field is not
      explicitly specified ("qualified"), then the server will
      apply a set of hueristics to improve the query.  Eg, in
      this case "cancer" is first parsed as
        cancer[All Fields]
      then turned into the query
        "neoplasms"[MeSH Terms] OR cancer[Text Word]

      Note that these terms are URL escaped.
      For details on how the translation is done, see
http://www.ncbi.nlm.nih.gov/entrez/query/static/help/pmhelp.html#AutomaticTermMapping

  TranslationStack -- The (possibly 'improved' query) fully
      parsed out and converted into a postfix (RPN) notation.
      The above example is written in the Entrez query language as

        ("neoplasms"[MeSH Terms] OR cancer[Text Word]) AND
                     2002/10/30:2002/12/29[edat]
      Note that these terms are *not* URL escaped.  Nothing like
      a bit of inconsistency for the soul.

      The "Count" field shows how many matches were found for each
      term of the expression.  I don't know what "Explode" does.


Let's get more information about the first record, which has an id of
12503096.  There are two ways to query for information, one uses a set
of identifiers and the other uses the history.  I'll talk about the
history one in a bit.  To use a set of identifiers you need to make a
DBIds object containing the that list.

>>> dbids = EUtils.DBIds("pubmed", ["12503096"])
>>>

Now get the summary using dbids

>>> infile = eutils.esummary_using_dbids(dbids)
>>> print infile.read()
<?xml version="1.0"?>
<!DOCTYPE eSummaryResult PUBLIC "-//NLM//DTD eSummaryResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eSummary_020511.dtd">
<eSummaryResult>
<DocSum>
        <Id>12503096</Id>
        <Item Name="PubDate" Type="Date">2003 Jan 30</Item>
        <Item Name="Source" Type="String">Am J Med Genet</Item>
        <Item Name="Authors" Type="String">Coyne JC, Kruus L, Racioppo M, Calzone KA, Armstrong K</Item>
        <Item Name="Title" Type="String">What do ratings of cancer-specific distress mean among women at high risk of breast and ovarian cancer?</Item>
        <Item Name="Volume" Type="String">116</Item>
        <Item Name="Pages" Type="String">222-8</Item>
        <Item Name="EntrezDate" Type="Date">2002/12/28 04:00</Item>
        <Item Name="PubMedId" Type="Integer">12503096</Item>
        <Item Name="MedlineId" Type="Integer">22390532</Item>
        <Item Name="Lang" Type="String">English</Item>
        <Item Name="PubType" Type="String"></Item>
        <Item Name="RecordStatus" Type="String">PubMed - in process</Item>
        <Item Name="Issue" Type="String">3</Item>
        <Item Name="SO" Type="String">2003 Jan 30;116(3):222-8</Item>
        <Item Name="DOI" Type="String">10.1002/ajmg.a.10844</Item>
        <Item Name="JTA" Type="String">3L4</Item>
        <Item Name="ISSN" Type="String">0148-7299</Item>
        <Item Name="PubId" Type="String"></Item>
        <Item Name="PubStatus" Type="Integer">4</Item>
        <Item Name="Status" Type="Integer">5</Item>
        <Item Name="HasAbstract" Type="Integer">1</Item>
        <Item Name="ArticleIds" Type="List">
                <Item Name="PubMedId" Type="String">12503096</Item>
                <Item Name="DOI" Type="String">10.1002/ajmg.a.10844</Item>
                <Item Name="MedlineUID" Type="String">22390532</Item>
        </Item>
</DocSum>
</eSummaryResult>
>>>

This is just a summary.  To get the full details, including an
abstract (if available) use the 'efetch' method.  I'll only print a
bit to convince you it has an abstract.

>>> s = eutils.efetch_using_dbids(dbids).read()
>>> print s[587:860]
<ArticleTitle>What do ratings of cancer-specific distress mean among women at high risk of breast and ovarian cancer?</ArticleTitle>
<Pagination>
<MedlinePgn>222-8</MedlinePgn>
</Pagination>
<Abstract>
<AbstractText>Women recruited from a hereditary cancer registry provided
>>>

Suppose instead you want the data in a text format.  Different
databases have different text formats.  For example, PubMed has a
"docsum" format which gives just the summary of a document and
"medline" format as needed for a citation database.  To get these, use
a "text" "retmode" ("return mode") and select the appropriate
"rettype" ("return type").

Here are examples of those two return types

>>> print eutils.efetch_using_dbids(dbids, "text", "docsum").read()[:497]
1:  Coyne JC, Kruus L, Racioppo M, Calzone KA, Armstrong K.
What do ratings of cancer-specific distress mean among women at high risk of breast and ovarian cancer?
Am J Med Genet. 2003 Jan 30;116(3):222-8.
PMID: 12503096 [PubMed - in process]
>>> print eutils.efetch_using_dbids(dbids, "text", "medline").read()[:369]
UI  - 22390532
PMID- 12503096
DA  - 20021227
IS  - 0148-7299
VI  - 116
IP  - 3
DP  - 2003 Jan 30
TI  - What do ratings of cancer-specific distress mean among women at high risk
      of breast and ovarian cancer?
PG  - 222-8
AB  - Women recruited from a hereditary cancer registry provided ratings of
      distress associated with different aspects of high-risk status
>>> 

It's also possible to get a list of records related to a given
article.  This is done through the "elink" method.  For example,
here's how to get the list of PubMed articles related to the above
PubMed record.  (Again, truncated because otherwise there is a lot of
data.)

>>> print eutils.elink_using_dbids(dbids).read()[:590]
<?xml version="1.0"?>
<!DOCTYPE eLinkResult PUBLIC "-//NLM//DTD eLinkResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eLink_020511.dtd">
<eLinkResult>
<LinkSet>
        <DbFrom>pubmed</DbFrom>
        <IdList>
                <Id>12503096</Id>
        </IdList>
        <LinkSetDb>
                <DbTo>pubmed</DbTo>
                <LinkName>pubmed_pubmed</LinkName>
                <Link>
                        <Id>12503096</Id>
                        <Score>2147483647</Score>
                </Link>
                <Link>
                        <Id>11536413</Id>
                        <Score>30817790</Score>
                </Link>
                <Link>
                        <Id>11340606</Id>
                        <Score>29939219</Score>
                </Link>
                <Link>
                        <Id>10805955</Id>
                        <Score>29584451</Score>
                </Link>
>>>

For a change of pace, let's work with the protein database to learn
how to work with history.  Suppose I want to do a multiple sequene
alignment of bacteriorhodopsin with all of its neighbors, where
"neighbors" is defined by NCBI.  There are good programs for this -- I
just need to get the records in the right format, like FASTA.

The bacteriorhodopsin I'm interested in is BAA75200, which is
GI:4579714, so I'll start by asking for its neighbors.

>>> results = eutils.elink_using_dbids(
...             EUtils.DBIds("protein", ["4579714"]),
...             db = "protein").read()
>>> print results[:454]
<?xml version="1.0"?>
<!DOCTYPE eLinkResult PUBLIC "-//NLM//DTD eLinkResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eLink_020511.dtd">
<eLinkResult>
<LinkSet>
        <DbFrom>protein</DbFrom>
        <IdList>
                <Id>4579714</Id>
        </IdList>
        <LinkSetDb>
                <DbTo>protein</DbTo>
                <LinkName>protein_protein</LinkName>
                <Link>
                        <Id>4579714</Id>
                        <Score>2147483647</Score>
                </Link>
                <Link>
                        <Id>11277596</Id>
                        <Score>1279</Score>
                </Link>
>>>

Let's get all the <Id> fields.  (While the following isn't a good way
to parse XML, it is easy to understand and works well enough for this
example.)  Note that I remove the first <Id> because that's from the
query and not from the results.

>>> import re
>>> ids = re.findall(r"<Id>(\d+)</Id>", results)
>>> ids = ids[1:]
>>> len(ids)
222
>>> dbids = EUtils.DBIds("protein", ids)
>>> 

That's a lot of records.  I could use 'efetch_using_dbids' but there's
a problem with that.  Efetch uses the HTTP GET protocol to pass
information to the EUtils server.  ("GET" is what's used when you type
a URL in the browser.)  Each id takes about 9 characters, so the URL
would be over 2,000 characters long.  This may not work on some
systems, for example, some proxies do not support long URLs.  (Search
for "very long URLs" for examples.)

Instead, we'll upload the list to the server then fetch the FASTA
version using the history.

The first step is to upload the data.  We want to put that into the
history so we set 'usehistory' to true.  There's no existing history
so the webenv string is None.


>>> print eutils.epost(dbids, usehistory = 1, webenv = None).read()
<?xml version="1.0"?>
<!DOCTYPE ePostResult PUBLIC "-//NLM//DTD ePostResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/ePost_020511.dtd">
<ePostResult>
        <QueryKey>1</QueryKey>
        <WebEnv>%7BPgTHRHFBsJfC%3C%5C%5C%5B%3EAfJCKQ%5Ey%60%3CGkH%5DH%5E%3DJHGBKAJ%3F%40CbCiG%3FE%3C</WebEnv>
</ePostResult>

>>>

This says that the identifiers were saved as query #1, which will be
used later on as the "query_key" field.  The WebEnv is a cookie (or
token) used to tell the server where to find that query.  The WebEnv
changes after every history-enabled ESearch or EPost so you'll need to
parse the output from those to get the new WebEnv field.  You'll also
need to unquote it since it is URL-escaped.

Also, you will need to pass in the name of the database used for the
query in order to access the history.  Why?  I don't know -- I figure
the WebEnv and query_key should be enough to get the database name.

>>> import urllib
>>> webenv = urllib.unquote("%7BPgTHRHFBsJfC%3C%5C%5C%5B%3EAfJCKQ%5Ey%60%3CGkH%5DH%5E%3DJHGBKAJ%3F%40CbCiG%3FE%3C")
>>> print webenv
{PgTHRHFBsJfC<\\[>AfJCKQ^y`<GkH]H^=JHGBKAJ?@CbCiG?E<
>>>

Okay, now to get the data in FASTA format.  Notice that I need the
'retmax' in order to include all the records in the result.  (The
default is 20 records.)

>>> fasta = eutils.efetch_using_history("protein", webenv, query_key = "1",
...                                     retmode = "text", rettype = "fasta",
...                                     retmax = len(dbids)).read()
>>> fasta.count(">")
222
>>> print fasta[:694]
>gi|14194475|sp|O93742|BACH_HALSD Halorhodopsin (HR)
MMETAADALASGTVPLEMTQTQIFEAIQGDTLLASSLWINIALAGLSILLFVYMGRNLEDPRAQLIFVAT
LMVPLVSISSYTGLVSGLTVSFLEMPAGHALAGQEVLTPWGRYLTWALSTPMILVALGLLAGSNATKLFT
AVTADIGMCVTGLAAALTTSSYLLRWVWYVISCAFFVVVLYVLLAEWAEDAEVAGTAEIFNTLKLLTVVL
WLGYPIFWALGAEGLAVLDVAVTSWAYSGMDIVAKYLFAFLLLRWVVDNERTVAGMAAGLGAPLARCAPA
DD
>gi|14194474|sp|O93741|BACH_HALS4 Halorhodopsin (HR)
MRSRTYHDQSVCGPYGSQRTDCDRDTDAGSDTDVHGAQVATQIRTDTLLHSSLWVNIALAGLSILVFLYM
ARTVRANRARLIVGATLMIPLVSLSSYLGLVTGLTAGPIEMPAAHALAGEDVLSQWGRYLTWTLSTPMIL
LALGWLAEVDTADLFVVIAADIGMCLTGLAAALTTSSYAFRWAFYLVSTAFFVVVLYALLAKWPTNAEAA
GTGDIFGTLRWLTVILWLGYPILWALGVEGFALVDSVGLTSWGYSLLDIGAKYLFAALLLRWVANNERTI
AVGQRSGRGAIGDPVED
>>> 

To round things out, here's a query which refines the previous query.
I want to get all records from the first search which also have the
word "Structure" in them.  (My background was originally structural
biophysics, whaddya expect?  :)

>>> print eutils.search("#1 AND structure", db = "protein", usehistory = 1,
...                     webenv = webenv).read()
<?xml version="1.0"?>
<!DOCTYPE eSearchResult PUBLIC "-//NLM//DTD eSearchResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eSearch_020511.dtd">
<eSearchResult>
        <Count>67</Count>
        <RetMax>20</RetMax>
        <RetStart>0</RetStart>
        <QueryKey>2</QueryKey>
        <WebEnv>UdvMf%3F%60G%3DIE%60bG%3DGec%3E%3D%3Cbc_%5DgBAf%3EAi_e%5EAJcHgDi%3CIqGdE%7BmC%3C</WebEnv>
        <IdList>
                <Id>461608</Id>
                <Id>114808</Id>
                <Id>1364150</Id>
                <Id>1363466</Id>
                <Id>1083906</Id>
                <Id>99232</Id>
                <Id>99212</Id>
                <Id>81076</Id>
                <Id>114811</Id>
                <Id>24158915</Id>
                <Id>24158914</Id>
                <Id>24158913</Id>
                <Id>1168615</Id>
                <Id>114812</Id>
                <Id>114809</Id>
                <Id>17942995</Id>
                <Id>17942994</Id>
                <Id>17942993</Id>
                <Id>20151159</Id>
                <Id>20150922</Id>
        </IdList>
        <TranslationSet>
        </TranslationSet>
        <TranslationStack>
                <TermSet>
                        <Term>#1</Term>
                        <Field>All Fields</Field>
                        <Count>222</Count>
                        <Explode>Y</Explode>
                </TermSet>
                <TermSet>
                        <Term>structure[All Fields]</Term>
                        <Field>All Fields</Field>
                        <Count>142002</Count>
                        <Explode>Y</Explode>
                </TermSet>
                <OP>AND</OP>
        </TranslationStack>
</eSearchResult>

>>> 

One last thing about history.  It doesn't last very long -- perhaps an
hour or so.  (Untested.)  You may be able to toss it some keep-alive
signal every once in a while.  Or you may want to keep 

The known 'db' fields and primary IDs (if known) are
  genome -- GI number
  nucleotide -- GI number
  omim  -- MIM number
  popset -- GI number
  protein -- GI number
  pubmed  -- PMID
  sequences (not available; this will combine all sequence databases)
  structure -- MMDB ID
  taxonomy -- TAXID

The 'field' parameter is different for different databases.  The
fields for PubMed are listed at

http://www.ncbi.nlm.nih.gov/entrez/query/static/help/pmhelp.html#SearchFieldDescriptionsandTags

  Affiliation -- AD
  All Fields -- All
  Author -- AU
  EC/RN Number -- RN
  Entrez Date -- EDAT  (also valid for 'datetype')
  Filter -- FILTER
  Issue -- IP
  Journal Title -- TA
  Language -- LA
  MeSH Date -- MHDA  (also valid for 'datetype')
  MeSH Major Topic -- MAJR
  MeSH Subheadings -- SH
  MeSH Terms -- MH
  Pagination -- PG
  Personal Name as Subject -- PS
  Publication Date -- DP  (also valid for 'datetype')
  Publication Type -- PT
  Secondary Source ID -- SI
  Subset -- SB
  Substance Name -- NM
  Text Words -- TW
  Title -- TI
  Title/Abstract -- TIAB
  Unique Identifiers -- UID
  Volume -- VI

The fields marked as 'datetype' can also be used for date searches.
Date searches can be done in the query (for example, as

   1990/01/01:1999/12/31[edat]

or by passing a WithinNDays or DateRange field to the 'date' parameter
of the search.


Please pay attention to the usage limits!  The are listed at
  http://www.ncbi.nlm.nih.gov/entrez/query/static/eutils_help.html

At the time of this writing they are:
    * Run retrieval scripts on weekends or between 9 PM and 5 AM ET
            weekdays for any series of more than 100 requests.
    * Make no more than one request every 3 seconds.
    * Only 5000 PubMed records may be retrieved in a single day.

    * NCBI's Disclaimer and Copyright notice must be evident to users
      of your service.  NLM does not hold the copyright on the PubMed
      abstracts the journal publishers do.  NLM provides no legal
      advice concerning distribution of copyrighted materials, consult
      your legal counsel.

(Their disclaimer is at
       http://www.ncbi.nlm.nih.gov/About/disclaimer.html )


""" # "  # Emacs cruft

import urllib, urllib2, cStringIO
import time

DUMP_URL = 0
DUMP_RESULT = 0

# These tell NCBI who is using the tool.  They are meant to provide
# hints to NCBI about how their service is being used and provide a
# means of getting ahold of the author.
#
# To use your own values, pass them in to the EUtils constructor.
#
TOOL = "EUtils_Python_client"
EMAIL = "biopython-dev@biopython.org"

assert " " not in TOOL
assert " " not in EMAIL

def _dbids_to_id_string(dbids):
    """Internal function: convert a list of ids to a comma-seperated string"""
    # NOTE: the server strips out non-numeric characters
    # Eg, "-1" is treated as "1".  So do some sanity checking.
    # XXX Should I check for non-digits?
    # Are any of the IDs non-integers?
    if not dbids:
        raise TypeError("dbids list must have at least one term")
    for x in dbids.ids:
        if "," in x:
            raise TypeError("identifiers cannot contain a comma: %r " %
                            (x,))
    id_string = ",".join(dbids.ids)
    assert id_string.count(",") == len(dbids.ids)-1, "double checking"
    return id_string

#Record the time at module level, in case the user has multiple copies
#of the ThinClient class in operation at once.
_open_previous = time.time()

class ThinClient:
    """Client-side interface to the EUtils services

    See the module docstring for much more complete information.
    """
    def __init__(self,
                 opener = None,
                 tool = TOOL,
                 email = EMAIL,
                 baseurl = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/"):
        """opener = None, tool = TOOL, email = EMAIL, baseurl = ".../eutils/"

        'opener' -- an object which implements the 'open' method like a
             urllib2.OpenDirector.  Defaults to urllib2.build_opener()

        'tool' -- the term to use for the 'tool' field, used by NCBI to
             track which programs use their services.  If you write your
             own tool based on this package, use your own tool name.

        'email' -- a way for NCBI to contact you (the developer, not
             the user!) if there are problems and to tell you about
             updates or changes to their system.

        'baseurl' -- location of NCBI's EUtils directory.  Shouldn't need
             to change this at all.
        """
        
        if tool is not None and " " in tool:
            raise TypeError("No spaces allowed in 'tool'")
        if email is not None and " " in email:
            raise TypeError("No spaces allowed in 'email'")

        if opener is None:
            opener = urllib2.build_opener()

        self.opener = opener
        self.tool = tool
        self.email = email
        self.baseurl = baseurl

    def _fixup_query(self, query):
        """Internal function to add and remove fields from a query"""
        q = query.copy()

        # Set the 'tool' and 'email' fields
        q["tool"] = self.tool
        q["email"] = self.email

        # Kinda cheesy -- shouldn't really do this here.
        # If 'usehistory' is true, use the value of 'Y' instead.
        # Otherwise, don't use history
        if "usehistory" in q:
            if q["usehistory"]:
                q["usehistory"] = "y"
            else:
                q["usehistory"] = None

        # This will also remove the history, email, etc. fields
        # if they are set to None.
        for k, v in q.items():
            if v is None:
                del q[k]

        # Convert the query into the form needed for a GET.
        return urllib.urlencode(q)

    def _wait(self, delay = 3.0) :
        """Enforce the NCBI requirement of one request every three seconds.

        Ideally the calling code would have respected the 3 second rule,
        but as this often hasn't happenend we check this here.

        wait - number of seconds between queries."""
        global _open_previous
        wait = _open_previous + delay - time.time()
        if wait > 0:
            time.sleep(wait)
        _open_previous = time.time()
    
    def _get(self, program, query):
        """Internal function: send the query string to the program as GET"""
        # NOTE: epost uses a different interface

        self._wait()

        q = self._fixup_query(query)
        url = self.baseurl + program + "?" + q
        if DUMP_URL:
            print "Opening with GET:", url
        if DUMP_RESULT:
            print " ================== Results ============= "
            s = self.opener.open(url).read()
            print s
            print " ================== Finished ============ "
            return cStringIO.StringIO(s)
        return self.opener.open(url)
            
    def esearch(self,
                term,           # In Entrez query language
                db = "pubmed",  # Required field, default to PubMed
                field = None,   # Field to use for unqualified words
                daterange = None,    # Date restriction

                retstart = 0,
                retmax = 20,     # Default from NCBI is 20, so I'll use that

                usehistory = 0,  # Enable history tracking
                webenv = None,   # If given, add to an existing history
                ):
        
        """term, db="pubmed", field=None, daterange=None, retstart=0, retmax=20, usehistory=0, webenv=none

        Search the given database for records matching the query given
        in the 'term'.  See the module docstring for examples.

        'term' -- the query string in the Entrez query language; see
           http://www.ncbi.nlm.nih.gov/entrez/query/static/help/pmhelp.html
        'db' -- the database to search

        'field' -- the field to use for unqualified words
                  Eg, "dalke[au] AND gene" with field==None becomes
                    dalke[au] AND (genes[MeSH Terms] OR gene[Text Word]
                  and "dalke[au] AND gene" with field=="au" becomes
                    dalke[au] AND genes[Author]
                 (Yes, I think the first "au" should be "Author" too)

        'daterange' -- a date restriction; either WithinNDays or DateRange
        'retstart' -- include identifiers in the output, starting with
                  position 'retstart' (normally starts with 0)
        'retmax' -- return at most 'retmax' identifiers in the output
                  (if not specified, NCBI returns 20 identifiers)

        'usehistory' -- flag to enable history tracking
        'webenv' -- if this string is given, add the search results
                  to an existing history. (WARNING: the history disappers
                  after about an hour of non-use.)

        You will need to parse the output XML to get the new QueryKey
        and WebEnv fields.

        Returns an input stream from an HTTP request.  The stream
        contents are in XML.
        """
        query = {"term": term,
                 "db": db,
                 "field": field,
                 "retstart": retstart,
                 "retmax": retmax,
                 "usehistory": usehistory,
                 "WebEnv": webenv,
                 }
        if daterange is not None:
            query.update(daterange.get_query_params())
        
        return self._get(program = "esearch.fcgi", query = query)
    
    def epost(self,
              dbids,

              webenv = None,    # If given, add to an existing history
              ):
        """dbids, webenv = None

        Create a new collection in the history containing the given
        list of identifiers for a database.

        'dbids' -- a DBIds, which contains the database name and
             a list of identifiers in that database
        'webenv' -- if this string is given, add the collection
                  to an existing history. (WARNING: the history disappers
                  after about an hour of non-use.)

        You will need to parse the output XML to get the new QueryKey
        and WebEnv fields.  NOTE: The order of the IDs on the server
        is NOT NECESSARILY the same as the upload order.

        Returns an input stream from an HTTP request.  The stream
        contents are in XML.
        """
        id_string = _dbids_to_id_string(dbids)

        # Looks like it will accept *any* ids.  Wonder what that means.
        program = "epost.fcgi"
        query = {"id": id_string,
                 "db": dbids.db,
                 "WebEnv": webenv,
                 }
        q = self._fixup_query(query)

        self._wait()

        # Need to use a POST since the data set can be *very* long;
        # even too long for GET.
        if DUMP_URL:
            print "Opening with POST:", self.baseurl + program + "?" + q
        if DUMP_RESULT:
            print " ================== Results ============= "
            s = self.opener.open(self.baseurl + program, q).read()
            print s
            print " ================== Finished ============ "
            return cStringIO.StringIO(s)
        return self.opener.open(self.baseurl + program, q)

    def esummary_using_history(self,
                               db,  # This is required.  Don't use a
                                    # default here because it must match
                                    # that of the webenv
                               webenv,
                               query_key,
                               retstart = 0,
                               retmax = 20,
                               retmode = "xml",  # any other modes?
                               ):
        """db, webenv, query_key, retstart = 0, retmax = 20, retmode = "xml"

        Get the summary for a collection of records in the history

        'db' -- the database containing the history/collection
        'webenv' -- the WebEnv cookie for the history
        'query_key' -- the collection in the history
        'retstart' -- get the summaries starting with this position
        'retmax' -- get at most this many summaries
        'retmode' -- can only be 'xml'.  (Are there others?)

        Returns an input stream from an HTTP request.  The stream
        contents are in 'retmode' format.
        """
        return self._get(program = "esummary.fcgi",
                         query = {"db": db,
                                  "WebEnv": webenv,
                                  "query_key": query_key,
                                  "retstart": retstart,
                                  "retmax": retmax,
                                  "retmode": retmode,
                                  })
                        
    def esummary_using_dbids(self,
                             dbids,
                             retmode = "xml",  # any other modes?
                             ):
        """dbids, retmode = "xml"

        Get the summary for records specified by identifier

        'dbids' -- a DBIds containing the database name and list
                       of record identifiers
        'retmode' -- can only be 'xml'

        Returns an input stream from an HTTP request.  The stream
        contents are in 'retmode' format.
        """
        
        id_string = _dbids_to_id_string(dbids)
        return self._get(program = "esummary.fcgi",
                         query = {"id": id_string,
                                  "db": dbids.db,
                                  # "retmax": len(dbids.ids), # needed?
                                  "retmode": retmode,
                                  })

    def efetch_using_history(self,
                             db,
                             webenv,
                             query_key,

                             retstart = 0,
                             retmax = 20,

                             retmode = None,
                             rettype = None,

                             # sequence only
                             seq_start = None,
                             seq_stop = None,
                             strand = None,
                             complexity = None,
                             ):
        """db, webenv, query_key, retstart=0, retmax=20, retmode=None, rettype=None, seq_start=None, seq_stop=None, strand=None, complexity=None

        Fetch information for a collection of records in the history,
        in a variety of formats.

        'db' -- the database containing the history/collection
        'webenv' -- the WebEnv cookie for the history
        'query_key' -- the collection in the history
        'retstart' -- get the formatted data starting with this position
        'retmax' -- get data for at most this many records

        These options work for sequence databases

        'seq_start' -- return the sequence starting at this position.
               The first position is numbered 1
        'seq_stop' -- return the sequence ending at this position
               Includes the stop position, so seq_start = 1 and
               seq_stop = 5 returns the first 5 bases/residues.
        'strand' -- strand.  Use EUtils.PLUS_STRAND (== 1) for plus
               strand and EUtils.MINUS_STRAND (== 2) for negative
        'complexity' -- regulates the level of display.  Options are
            0 - get the whole blob
            1 - get the bioseq for gi of interest (default in Entrez)
            2 - get the minimal bioseq-set containing the gi of interest
            3 - get the minimal nuc-prot containing the gi of interest
            4 - get the minimal pub-set containing the gi of interest
        
        http://www.ncbi.nlm.nih.gov/entrez/query/static/efetchseq_help.html

        The valid retmode and rettype values are

        For publication databases (omim, pubmed, journals) the
        retmodes are 'xml', 'asn.1', 'text', and 'html'.
        
          If retmode == xml     ---> XML (default)
          if retmode == asn.1   ---> ASN.1

          The following rettype values work for retmode == 'text'.
           
             docsum    ----> author / title / cite / PMID
             brief     ----> a one-liner up to about 66 chars
             abstract  ----> cite / title / author / dept /
                                    full abstract / PMID
             citation  ----> cite / title / author / dept /
                                    full abstract / MeSH terms /
                                    substances / PMID
             medline   ----> full record in medline format
             asn.1     ----> full record in one ASN.1 format
             mlasn1    ----> full record in another ASN.1 format
             uilist    ----> list of uids, one per line
             sgml      ----> same as retmode="xml"

        Sequence databases (genome, protein, nucleotide, popset)
        also have retmode values of 'xml', 'asn.1', 'text', and
        'html'.

          If retmode == 'xml'   ---> XML (default; only supports
                                        rettype == 'native')
          If retmode == 'asn.1' ---> ASN.1 text (only works for rettype
                                        of 'native' and 'sequin')

          The following work with a retmode of 'text' or 'html' 

             native    ----> Default format for viewing sequences
             fasta     ----> FASTA view of a sequence
             gb        ---->  GenBank view for sequences, constructed sequences
                          will be shown as contigs (by pointing to its parts).
                          Valid for nucleotides.
             gbwithparts --> GenBank view for sequences, the sequence will
                          always be shown.  Valid for nucleotides.
             est       ----> EST Report.  Valid for sequences from
                          dbEST database.
             gss       ----> GSS Report.  Valid for sequences from dbGSS
                          database.
             gp        ----> GenPept view.  Valid for proteins.
             seqid     ----> To convert list of gis into list of seqids
             acc       ----> To convert list of gis into list of accessions

             # XXX TRY THESE
             fasta_xml
             gb_xml
             gi  (same as uilist?)

             
        
        A retmode of 'file' is the same as 'text' except the data is
        sent with a Content-Type of application/octet-stream, which tells
        the browser to save the data to a file.

        A retmode of 'html' is the same as 'text' except a HTML header
        and footer are added and special character are properly escaped.

        Returns an input stream from an HTTP request.  The stream
        contents are in the requested format.
        """

        # NOTE: found the list of possible values by sending illegal
        # parameters, to see which comes up as an error message.  Used
        # that to supplement the information from the documentation.
        # Looks like efetch is based on pmfetch code and uses the same
        # types.
        
        # if retmax is specified and larger than 500, NCBI only returns 
        # 500 sequences. Removing it from the URL relieves this constraint.
        # To get around this, if retstart is 0 and retmax is greater than 500,
        # we set retmax to be None.
        if retstart == 0 and retmax > 500:
            retmax = None
        return self._get(program = "efetch.fcgi",
                         query = {"db": db,
                                  "WebEnv": webenv,
                                  "query_key": query_key,
                                  "retstart": retstart,
                                  "retmax": retmax,
                                  "retmode": retmode,
                                  "rettype": rettype,
                                  "seq_start": seq_start,
                                  "seq_stop": seq_stop,
                                  "strand": strand,
                                  "complexity": complexity,
                                  })

    def efetch_using_dbids(self,
                       dbids,
                       retmode = None,
                       rettype = None,

                       # sequence only
                       seq_start = None,
                       seq_stop = None,
                       strand = None,
                       complexity = None,
                       ):
        """dbids, retmode = None, rettype = None, seq_start = None, seq_stop = None, strand = None, complexity = None

        Fetch information for records specified by identifier

        'dbids' -- a DBIds containing the database name and list
                       of record identifiers
        'retmode' -- See the docstring for 'efetch_using_history'
        'rettype' -- See the docstring for 'efetch_using_history'

        These options work for sequence databases

        'seq_start' -- return the sequence starting at this position.
               The first position is numbered 1
        'seq_stop' -- return the sequence ending at this position
               Includes the stop position, so seq_start = 1 and
               seq_stop = 5 returns the first 5 bases/residues.
        'strand' -- strand.  Use EUtils.PLUS_STRAND (== 1) for plus
               strand and EUtils.MINUS_STRAND (== 2) for negative
        'complexity' -- regulates the level of display.  Options are
            0 - get the whole blob
            1 - get the bioseq for gi of interest (default in Entrez)
            2 - get the minimal bioseq-set containing the gi of interest
            3 - get the minimal nuc-prot containing the gi of interest
            4 - get the minimal pub-set containing the gi of interest
        
        Returns an input stream from an HTTP request.  The stream
        contents are in the requested format.
        """
        id_string = _dbids_to_id_string(dbids)
        return self._get(program = "efetch.fcgi",
                         query = {"id": id_string,
                                  "db": dbids.db,
                                  # "retmax": len(dbids.ids), # needed?
                                  "retmode": retmode,
                                  "rettype": rettype,
                                  "seq_start": seq_start,
                                  "seq_stop": seq_stop,
                                  "strand": strand,
                                  "complexity": complexity,
                                  })

    def elink_using_history(self,
                            dbfrom,
                            webenv,
                            query_key,

                            db = "pubmed",

                            retstart = 0,
                            retmax = 20,

                            cmd = "neighbor",
                            retmode = None,

                            term = None,
                            field = None,
                            # "Date limits not valid for link commands"
                            daterange = None,
                            ):
        """dbfrom, webenv, query_key, db="pubmed", retstart=0, retmax=20, cmd="neighbor", retmode=None, term=None, field=None, daterange=None, 

        Find records related (in various ways) to a collection of
        records in the history.

        'dbfrom' -- this is the name of the database containing the
                      collection of record.  NOTE!  For the other methods
                      this is named 'db'.  But I'm keeping NCBI's notation.
                      This is where the records come FROM.
        'webenv' -- the WebEnv cookie for the history
        'query_key' -- the collection in the history

        'db' -- Where the records link TO.  This is where you want to
                  find the new records.  For example, if you want to
                  find PubMed records related to a protein then 'dbfrom'
                  is 'protein' and 'db' is 'pubmed'

        'cmd'-- one of the following (unless specified, retmode is the
                      default value, which returns data in XML)
           neighbor:  Display neighbors and their scores by database and ID.
                         (This is the default 'cmd'.)
           prlinks: List the hyperlink to the primary LinkOut provider
                         for multiple IDs and database.
                    When retmode == 'ref' this URL redirects the browser
                         to the primary LinkOut provider for a single ID
                         and database.
           llinks:  List LinkOut URLs and Attributes for multiple IDs
                         and database.
           lcheck:  Check for the existence (Y or N) of an external
                         link in for multiple IDs and database.
           ncheck:  Check for the existence of a neighbor link for
                         each ID, e.g., Related Articles in PubMed.

        'retstart' -- get the formatted data starting with this position
        'retmax' -- get data for at most this many records

        'retmode' -- only used with 'prlinks'

        'term' -- restrict results to records which also match this
                         Entrez search
        'field' -- the field to use for unqualified words

        'daterange' -- restrict results to records which also match this
                         date criteria; either WithinNDays or DateRange
                         NOTE: DateRange must have both mindate and maxdate

        Some examples:
          In PubMed, to get a list of "Related Articles"
            dbfrom = pubmed
            cmd = neighbor

          To get MEDLINE index only related article
            dbfrom = pubmed
            db = pubmed
            term = medline[sb]
            cmd = neighbor

          Given a PubMed record, find the related nucleotide records
            dbfrom = pubmed
            db = nucleotide  (or "protein" for related protein records)
            cmd = neighbor

          To get "LinkOuts" (external links) for a PubMed record set
            dbfrom = pubmed
            cmd = llinks

          Get the primary link information for a PubMed document; includes
              various hyperlinks, image URL for the provider, etc.
            dbfrom = pubmed
            cmd = prlinks
            (optional) retmode = "ref" (causes a redirect to the privder)
        
        Returns an input stream from an HTTP request.  The stream
        contents are in XML unless 'retmode' is 'ref'.
        """
        query = {"WebEnv": webenv,
                 "query_key": query_key,
                 "db": db,
                 "dbfrom": dbfrom,
                 "cmd": cmd,
                 "retstart": retstart,
                 "retmax": retmax,
                 "retmode": retmode,
                 "term": term,
                 "field": field,
                 }
        if daterange is not None:
            if daterange.mindate is None or daterange.maxdate is None:
                raise TypeError("Both mindate and maxdate must be set for eLink")
            query.update(daterange.get_query_params())
        return self._get(program = "elink.fcgi", query = query)

    def elink_using_dbids(self,
                          dbids,
                          db = "pubmed",

                          cmd = "neighbor",
                          
                          retmode = None,
                          term = None,
                          field = None,

                          daterange = None,
                          
                          ):
        """dbids, db="pubmed", cmd="neighbor", retmode=None, term=None, daterange=None

        Find records related (in various ways) to a set of records
        specified by identifier.

        'dbids' -- a DBIds containing the database name and list
                       of record identifiers
        'db' -- Where the records link TO.  This is where you want to
                  find the new records.  For example, if you want to
                  find PubMed records related to a protein then 'db'
                  is 'pubmed'.  (The database they are from is part
                  of the DBIds object.)

        'cmd' -- see the docstring for 'elink_using_history'
        'retmode' -- see 'elink_using_history'
        'term' -- see 'elink_using_history'
        'daterange' -- see 'elink_using_history'

        Returns an input stream from an HTTP request.  The stream
        contents are in XML unless 'retmode' is 'ref'.
        """
        id_string = _dbids_to_id_string(dbids)
        query = {"id": id_string,
                 "db": db,
                 "dbfrom": dbids.db,
                 "cmd": cmd,
                 "retmode": retmode,
                 "field" : field,
                 "term": term,
                 }
        if daterange is not None:
            import Datatypes
            if isinstance(daterange, Datatypes.DateRange) and \
               (daterange.mindate is None or daterange.maxdate is None):
                raise TypeError("Both mindate and maxdate must be set for eLink")
            query.update(daterange.get_query_params())
        
        return self._get(program = "elink.fcgi", query = query)
