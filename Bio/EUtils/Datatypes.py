"""Various EUtils datatypes."""

import re, types

class EUtilsError(Exception):
    """Base class for all EUtils-specific errors

    Contains a single error string -- use str(err) to get it.
    """
    pass

class EUtilsSearchError(EUtilsError):
    """Used when the ESearch XML says there is an ERROR

    The main error is in err.errmsg but more information
    may be available in err.errors or err.warnings.  Eg,
    the error message is often "Can't run executor" but
    you can get more information from the list of errors.
    
    """
    def __init__(self, errmsg, errors = None, warnings = None):
        EUtilsError.__init__(self, errmsg)

        if errors is None: errors = []
        if warnings is None: warnings = []

        self.errmsg = errmsg
        self.errors = errors
        self.warnings = warnings
    def __repr__(self):
        return "%s(%r, %r, %r)" % (self.__class__.__name__,
                                   self.errmsg, self.errors, self.warnings)
    def __str__(self):
        s = self.errmsg
        if self.errors:
            s = s + "; ERRORS: " + ", ".join(map(str, self.errors))
        if self.warnings:
            s = s + "; WARNINGS: " + ", ".join(map(str, self.warnings))
        return s.encode("latin1")
    


####################################
class DBIds:
    """Store a list of identifiers for a database

    This is used as input for the '*_using_dbids' functions.

    Constructed with the database name and list of identifier strings.
    
    """
    def __init__(self, db, ids):
        """db, ids

        'db' -- the database for those identifiers
        'ids' -- a list of identifiers for the given database
        """
        self.db = db
        self.ids = ids
    def __len__(self):
        """number of identifers"""
        return len(self.ids)
    def __getitem__(self, i):
        """get an identifier or a subset of the DBIds"""
        if isinstance(i, types.SliceType):
            # XXX Python 2.3 fixes this, I think
            # Either that, or I'm doing something wrong?
            step = i.step
            start = i.start
            if start is None: start = 0
            stop = i.stop
            if stop is None: stop = len(self.ids)
            if step is None:
                return self.__class__(self.db, self.ids[start:stop])
            else:
                return self.__class__(self.db, self.ids[start:stop:step])
        # XXX Should this return a DBIds as well?  Because of this, I nee
        # the 'item' method
        return self.ids[i]
    def item(self, i):
        """Get a DBIds containing the item at position i

        Can't use dbids[i] since that returns only the identifier.
        This returns a DBIds, which can be used for another request.
        """
        return self.__class__(self.db, [self.ids[i]])
    
    def __iter__(self):
        """Iterate over the list of identifiers"""
        return iter(self.ids)
    def __repr__(self):
        return "DBIds(%r, %r)" % (self.db, self.ids)
    def __eq__(self, other):
        """does this DBIds equal the other?

        The database names must match, but the identifiers
        themselves can be in any order.
        """
        if self.ids == other.ids:
            return self.db == other.db
        if self.db != other.db:
            return 0
        # Could be in a different order, and there may be non-unique
        # keys.  XXX use a sets.Set from Python 2.3?  But then
        # there won't be a simple mapping from id to efetch results.
        d1 = {}
        for x in self.ids:
            d1[x] = 0
        d2 = {}
        for x in other.ids:
            d2[x] = 0
        return d1 == d2
    def __ne__(self, other):
        """check if this isn't equal to the other DBIds"""
        return not self == other

    def __sub__(self, other):
        """DBIds of the identifiers in this set which aren't in the other"""
        if self.db != other.db:
            raise TypeError("Different databases: %r and %r" % (
                self.db, other.db))
        other_d = {}
        for x in other.ids:
            other_d[x] = 0
        new_ids = [x for x in self.ids if x not in other_d]
        return DBIds(self.db, new_ids)

class WithinNDays:
    """Restrict a search to matches in the last N days

    Eg, to see what's been published in PubMed about rabies
    in the last 20 days.

    client.search("rabies", daterange = WithinNDays(20, "pdat")
    """
    def __init__(self, ndays, datetype = None):
        """ndays, datetype = None

        'ndays' -- within this many days of now (the 'reldate' field
               of a search)
        'datetype' -- the date field to use (defaults to Entrez date,
               which is "edat")
        """
        self.ndays = ndays
        self.datetype = datetype
    def get_query_params(self):
        """returns the fields to add to the EUtils query

        This is an internal implementation feature you can ignore.
        """
        return {"reldate": self.ndays,
                "datetype": self.datetype}

# Could actually check the month and day fields...
_date_re_match = re.compile(r"\d{4}(/\d\d(/\d\d)?)?$").match

class DateRange:
    """Restrict a search to matches within a date range

    Some examples:
        matches between 1995 and 2000 -- DateRange("1995", "1999/12/31")
        matches before 1990 -- DateRange(maxdate = "1990/01/01")
        matches in 2002 or later -- DateRange(mindate = "2002/01/01")
        matches in June or July of 2001 -- DateRange("2001/06", "2001/07")
                 
    """
    def __init__(self, mindate = None, maxdate = None, datetype = None):
        """mindate = None, maxdate = None, datetype = None

        'mindate' -- matches must be on or after this date
        'maxdate' -- matches must be on or before this date
        'datetype' -- the date field to use for the search (defaults
             to Entrez date, which is "edat")

        At least one of mindate or maxdate must be specified.
        If mindate is omitted, all results on or before maxdate are returned.
        If maxdate is omitted, all results on or after mindate are returned.

        Dates must be formatted as 'YYYY/MM/DD', 'YYYY/MM', or 'YYYY'.
        """
        if mindate is None and maxdate is None:
            raise TypeError("Must specify at least one of mindate or maxdate")

        errinfo = None
        if mindate is not None and _date_re_match(mindate) is None:
            errinfo = ("mindate", mindate)
        elif maxdate is not None and _date_re_match(maxdate) is None:
            errinfo = ("maxdate", maxdate)
        if errinfo:
            raise TypeError(
                "%s is not in YYYY/MM/DD format (month and "
                "day are optional): %r" % errinfo)
        self.mindate = mindate
        self.maxdate = maxdate
        self.datetype = datetype

    def get_query_params(self):
        """returns the fields to add to the EUtils query

        This is an internal implementation feature you can ignore.
        """
        return {"mindate": str(self.mindate),
                "maxdate": str(self.maxdate),
                "datetype": self.datetype}

####################################

class Expression:
    """Base class for the Expression given in the eSearch output

    NCBI does some processing on the request.  They return the
    translated expression as part of the search results.  To get the
    expression as an Entrez string, use str(expression).

    iter(expression) traverses the expression tree in postfix order.
    """
    def __and__(self, other):
        """intersection of two expressions"""
        return And(self, other)
    def __or__(self, other):
        """union of two expressions"""
        return Or(self, other)
    def __iter__(self):
        """Traverse the tree in postfix order"""
        raise NotImplementedError

class Term(Expression):
    """Information about an Expression Term, which is the leaf node

    The fields are:
       term -- a word from the search term
       field -- the field searched by this term
       count -- the number of records matching this word
       explode -- no idea
    """
    def __init__(self, term, field, count, explode):
        self.term = term
        self.field = field
        self.count = count
        self.explode = explode
    def __str__(self):
        return self.term
    def __iter__(self):
        """Traverse the tree in postfix order"""
        yield self

class BinaryOp(Expression):
    """Base class for binary expressions.  Has a left and a right child"""
    def __init__(self, left, right):
        self.left = left
        self.right = right
    def __iter__(self):
        """Traverse the tree in postfix order"""
        for x in self.left:
            yield x
        for x in self.right:
            yield x
        yield self
        
# NCBI processes booleans left to right (no precedence)
# I'm not going to worry about using minimal parens,
# I'll just always put them around them
class And(BinaryOp):
    """intersection of two subexpressions"""
    def __str__(self):
        return "(%s AND %s)" % (self.left, self.right)

class Or(BinaryOp):
    """union two subexpressions"""
    def __str__(self):
        return "(%s OR %s)" % (self.left, self.right)

# NOT and BUTNOT
class Not(BinaryOp):
    """the set of the left child without elements from the right child

    This is used for something like "poliovirus NOT polio"
    """
    def __str__(self):
        return "(%s NOT %s)" % (self.left, self.right)

class Range(BinaryOp):
    """Used to store a date range"""
    def __init__(self, left, right):
        if left.field != right.field:
            raise TypeError("dates must have the same field: %r and %r" %
                            (left.field, right.field))
        BinaryOp.__init__(self, left, right)

    def __str__(self):
        i = self.left.term.rfind("[")
        if i == -1:
            i = len(self.left.term)
        x = self.left.term[:i]

        i = self.right.term.rfind("[")
        if i == -1:
            i = len(self.right.term)
        y = self.right.term[:i]
        
        return "%s:%s[%s]" % (x, y, self.left.field)

##################

class SearchResult:
    """Store results from a database search

    Attributes are:
       count -- total number of matches to the query
       retmax -- total number of identifiers requested
       retstart -- a search can return a portion of the total
           number of results.  retstart is the offset into this list
       ids -- matching identifiers (may be a subset of the full list)
       translation_set -- dict mapping an input name to the canonical
           form prefered by NCBI
       expression -- the full equery as understood by NCBI
       webenv -- the WebEnv string (if use_history is set)
       query_key -- the query_key (if use_history is set)
       errors -- list of Problems in the ErrorList
       warnings -- list of Problems in the WarningList
       timestamp -- timestamp (from time.time()) when this record
           was received from the server.

    Returns a list of identifers instead of a DBIds because the output
    from NCBI's eSearch doesn't include the database name.
    """
    def __init__(self,
                 count, retmax, retstart, ids,
                 translation_set, expression,
                 webenv, query_key, errors,
                 warnings, timestamp):
        self.count = count
        self.retmax = retmax
        self.retstart = retstart
        self.ids = ids
        self.translation_set = translation_set
        self.expression = expression
        self.webenv = webenv
        self.query_key = query_key
        self.errors = errors
        self.warnings = warnings
        self.timestamp = timestamp

class PostResult:
    """Store the results of a Post

    Attributes are:
      webenv -- the WebEnv string
      query_key -- the query_ket
      timestamp -- timestamp (from time.time()) when this record
           was received from the server.
    """
    def __init__(self, webenv, query_key, invalid_ids, timestamp):
        self.webenv = webenv
        self.query_key = query_key
        self.invalid_ids = invalid_ids
        self.timestamp = timestamp

class Summary:
    """Store information from calling eSummary

    Attributes are:
      id -- the identifier string for this record
      dataitems -- an OrderedDictList containing the parsed Item
         elements for this Summary.
    """
    def __init__(self, id, dataitems):
        self.id = id
        self.dataitems = dataitems
    def __repr__(self):
        return "Summary(%r, %r)" % (self.id, self.dataitems)
    def __str__(self):
        return "<Summary id=%s, %s>" % (self.id, self.dataitems)

# XXX Use the new 'datetime' module when 2.3 is out!
class Date:
    """Allow simple Date storage

    Parameters and attributes are 'year', 'month', and 'day'
    """
    def __init__(self, year, month, day):
        self.year = year
        self.month = month
        self.day = day
    def __repr__(self):
        return "%s(%r, %r, %r)" % (self.__class__.__name__,
                                   self.year, self.month, self.day)
    def __str__(self):
        return "%4d/%02d/%02d" % (self.year, self.month, self.day)
    def timetuple(self):
        """Return the 9-tuple needed by various time functions"""
        # NOTE: I don't yet deal with the last three fields
        # (day of week, day of year, isDST)
        return (self.year, self.month, self.day, 0, 0, 0, 0, 0, -1)
    def __eq__(self, other):
        """Are these two times equal?"""
        return (self.year == other.year and
                self.month == other.month and
                self.day == other.day)
    def __ne__(self, other):
        """Are these two times dissimilar?"""
        return not self == other


#     possible errors from eSearch
# <!ELEMENT        ErrorList      (PhraseNotFound*,FieldNotFound*)>
# <!ELEMENT        WarningList    (PhraseIgnored*,
#                                 QuotedPhraseNotFound*,
#                                 OutputMessage*)>

class Problem:
    """Base class for Search Errors or Warnings

    A problem has:
      text -- the text of the problem
      severity -- either Problem.ERROR or Problem.WARNING
      category -- how NCBI categorizes this problem
    """
    ERROR = "ERROR"
    WARNING = "WARNING"
    def __init__(self, text):
        self.text = text
    def __eq__(self, other):
        return (self.text == other.text and
                self.severity == other.severity and
                self.category == other.category)
    def __ne__(self, other):
        return not self == other
    def __repr__(self):
        return "%s(%r)" % (self.__class__.__name__, self.text)
    def __str__(self):
        return str(self.text)
    
class ErrorProblem(Problem):
    severity = Problem.ERROR

class WarningProblem(Problem):
    severity = Problem.WARNING
    
class PhraseNotFound(ErrorProblem):
    category = "PhraseNotFound"

class FieldNotFound(ErrorProblem):
    severity = Problem.ERROR
    category = "FieldNotFound"
                     
class PhraseIgnored(WarningProblem):
    category = "PhraseIgnored"
    
class QuotedPhraseNotFound(WarningProblem):
    category = "QuotedPhraseNotFound"
    
class OutputMessage(WarningProblem):
    category = "OutputMessage"

def _build_problem_mapping():
    """Internal: make a map from category name (in XML) to the right class"""
    mapping = {}
    for v in globals().values():
        try:
            if issubclass(v, Problem) and hasattr(v, "category"):
                mapping[v.category] = v
        except TypeError:
            pass
    return mapping

problem_category_mapping = _build_problem_mapping()


# elinks with cmd=="neighbor"
class Link:
    """Store neighbor Link information for a given record

    Attributes are;
      id -- the identifier used as the input for the neighbor request
      score -- the amount of similarity, high numbers are better
    """
    def __init__(self, id, score = None):
        self.id = id                      # in everything
        self.score = score                # in neighbor
    def __eq__(self, other):
        return self.id == other.id and self.score == other.score
    def __ne__(self, other):
        return not self == other
    def __repr__(self):
        return "Link(%r, %r)" % (self.id, self.score)

class IdCheck:
    """Store results from an lcheck link

    Attributes are:
       id -- the id of the requested record
       has_linkout -- boolean, either it does or doesn't
       has_neighbor -- boolean, either it does or doesn't
    """
    def __init__(self, id, has_linkout = 0, has_neighbor = 0):
        self.id = id
        self.has_linkout = has_linkout
        self.has_neighbor = has_neighbor
    def __eq__(self, other):
        return (self.id == other.id and
                self.has_linkout == other.has_linkout and
                self.has_neighbor == other.has_neighbor)
    def __ne__(self, other):
        return not self == other
    def __repr__(self):
        return "IdCheck(%r, %r, %r)" % (self.id, self.has_linkout, self.has_neighbor)

class LinkSetDb(object):
    """Used in eLink with cmd == neighbor

    Attributes are:
      dbto -- the links are TO this database name
      linkname -- the name for this set (eg, "pubmed_protein")
      links -- list of Links, one per matching record (includes score)
         List order is the sames as the XML, which is ordered from
         most likely to least.  The identifer is from 'dbto'
      info -- ignored; this is only used as a warning when there is
         an empty list

    You can also use
      dbids -- get a DBIds of dbto and the identifiers in each Link
    """
    def __init__(self, dbto, linkname, links = None, info = None):
        if links is None:
            if info is None:
                raise TypeError("At least one of 'links' and 'info' must be set")
            links = []
        self.dbto = dbto
        self.linkname = linkname
        self.links = links

    def _get_dbids(self):
        return DBIds(self.dbto, [link.id for link in self.links])
    dbids = property(_get_dbids)
    
    def __eq__(self, other):
        return (self.dbto == other.dbto and
                self.linkname == other.linkname and
                self.links == other.links)
    def __ne__(self, other):
        return not self == other
    def __repr__(self):
        return "LinkSetDb(%r, %r, %r)" % (self.dbto, self.linkname, self.links)
        
class NeighborLinkSet:
    """Results from an eLink neighbor search

    Attributes are:
      dbids -- the DBIds of the *REQUESTED* identifiers
      linksetdbs -- an OrderedMultiDict of LinkSetDb objects

    """
    def __init__(self, dbids, linksetdbs):
        self.dbids = dbids
        self.linksetdbs = linksetdbs
    def __eq__(self, other):
        return (self.dbids == other.dbids and
                self.linksetdbs == other.linksetdbs)
    def __ne__(self, other):
        return not self == other

    def __repr__(self):
        return "NeighborLinkSet(%r, %r)" % (self.dbids, self.linksetdbs)

# elinks with cmd in ("ncheck", "lcheck")
class CheckLinkSet(object):
    """Results from 'ncheck' and 'lcheck' searches

    This is used to check if a set of records has neighbors
    or links.

    Attributes are:
      dbfrom -- the database containing those records
      idchecks -- list of IdCheck objects, one per id

      dbids -- the DBIds make from dbfrom and the idchecks
    """
    def __init__(self, dbfrom, idchecks):
        self.dbfrom = dbfrom
        self.idchecks = idchecks

    def _get_dbids(self):
        return DBIds(self.dbfrom, [idcheck.id for idcheck in self.idchecks])
    dbids = property(_get_dbids)
    
    def __eq__(self, other):
        return (self.dbfrom == other.dbfrom and
                self.idchecks == other.idchecks)
    def __ne__(self, other):
        return not self == other
    def __repr__(self):
        return "CheckLinkSet(%r, %r)" % (self.dbfrom, self.idchecks)
                

# elinks with cmd == "llinks"
class Provider:
    """The Provider, as listed in 'llinks' (LinkOut)

    Attributes are:
        name -- name of the provider
        name_abbr -- an abbreviated name for the provider
        id -- a unique id for the provider
        url -- where to go for more information about the provider
        icon_url -- a small image to use for the provider
    
    """
    def __init__(self, name, name_abbr, id,
                 url = None, icon_url = None):
        self.name = name
        self.name_abbr = name_abbr
        self.id = id
        self.url = url
        self.icon_url = icon_url
    def __eq__(self, other):
        return (self.name == other.name and
                self.name_abbr == other.name_abbr and
                self.id == other.id and
                self.url == other.url and
                self.icon_url == other.icon_url)
    def __ne__(self, other):
        return not self == other
    def __repr__(self):
        return "Provider(%r, %r, %r, %r, %r)" % (
            self.name, self.name_abbr, self.id, self.url, self.icon_url)
    

class ObjUrl:
    """The ObjUrl containing LinkOut information for a record

    Attributes are:
      subject_types -- list of strings describing this link (0 or more)
      provider -- a Provider instance
      linkname -- a name used to categorize this link (optional)
      attributes -- list of attributes (text strings), (0 or more)
      url -- URL of the link (optional)
      iconurl -- URL containing image for this link (optional)
    """
    def __init__(self, subject_types, provider,
                 linkname = None, url = None, attributes = None):
        assert isinstance(subject_types, list)
        self.subject_types = subject_types
        self.provider = provider
        self.linkname = linkname
        if attributes is None:
            attributes = []
        self.url = url
        self.attributes = attributes
    def __eq__(self, other):
        return (self.linkname == other.linkname and
                self.subject_types == other.subject_types and
                self.url == other.url and
                self.attributes == other.attributes and
                self.provider == other.provider)
    def __ne__(self, other):
        return not self == other
    def __repr__(self):
        return "ObjUrl(%r, %r, %r, %r, %r)" % (
            self.subject_types, self.provider, self.linkname,
            self.url, self.attributes)
        
class IdUrlSet:
    """Set of ObjUrls for the record with the given 'id'"""
    def __init__(self, id, objurls):
        self.id = id
        self.objurls = objurls
    def __eq__(self, other):
        return (self.id == other.id and
                self.objurls == other.objurls)
    def __ne__(self, other):
        return not self == other
    def __repr__(self):
        return "IdUrlSet(%r, %r)" % (self.id, self.objurls)

class LinksLinkSet:
    """Results of an 'llink' (LinkOut) search

    Finds links from records in a given database to external
    resources.

    Fields are:
      dbfrom -- the database in which search started
      idurlset -- a list of IdUrlSet, one for each identifier
    """
      
    def __init__(self, dbfrom, idurlset):
        self.dbfrom = dbfrom
        self.idurlset = idurlset
    def __eq__(self, other):
        return (self.dbfrom == other.dbfrom and
                self.idurlset == other.idurlset)
    def __ne__(self, other):
        return not self == other
    def __repr__(self):
        return "LinksLinkSet(%r, %r)" % (self.dbfrom, self.idurlset)
