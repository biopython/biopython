import codecs, time, urllib, re, htmlentitydefs
from xml.sax import xmlreader, SAXException
import Datatypes, ReseekFile, MultiDict
from xml.sax.handler import feature_external_ges
import POM

def _construct_pattern():
    n = max([len(x) for x in htmlentitydefs.entitydefs.keys()])
    entity_pattern = re.compile(r"&([a-zA-Z]{1,%d});" % n)

    defs = {}
    for k, v in htmlentitydefs.entitydefs.items():
        if len(v) == 1:
            defs[k] = unicode(v, "latin-1")
        elif v[:2] == "&#" and v[-1] == ";":
            defs[k] = unichr(int(v[2:-1]))
        else:
            raise AssertionError("Unexpected entitydef value: %r" % v)

    return entity_pattern, defs

_entity_pattern, entitydefs = _construct_pattern()

def _load_module(name):
    mod = __import__(name)
    for term in name.split(".")[1:]:
        mod = getattr(mod, term)
    return mod

class GetObject:
    def __call__(self, obj):
        self.obj = obj

class UsePOMParser:
    def __init__(self, module_name):
        self.module_name = "Bio.EUtils.DTDs." + module_name

    def parse_using_dtd(self, file):
        module = _load_module(self.module_name)
        cb = GetObject()
        parser = POM.get_parser(callback = cb, module = module)
        # This tells the parser to not resolve the NCBI DTDs
        try:
            parser.setFeature(feature_external_ges, 0)
            SAXException
        except SAXException:
            pass
        parser.parse(file)
        return cb.obj

# Pull out the "ERROR", "ErrorList", and "WarningList" terms
def _check_for_errors(pom):
    errmsg = None
    errors = []
    warnings = []

    err = pom.get("ERROR", None)
    if err is not None:
        errmsg = err.tostring()

    for x in pom.get("ErrorList", []):
        errors.append(
            Datatypes.problem_category_mapping[x.__class__.__name__](
                      x.tostring()))
    for x in pom.get("WarningList", []):
        warnings.append(
            Datatypes.problem_category_mapping[x.__class__.__name__](
                        x.tostring()))
    return errmsg, errors, warnings


def _check_for_bad_input_stream(infile, force_encoding = 1):
    reseekfile = ReseekFile.ReseekFile(infile)
    s = reseekfile.read(500)
    reseekfile.seek(0)
    reseekfile.nobuffer()
    
    lines = s.split("\n")
    if len(lines) > 3:
        if lines[0] == "<Html>":
            if lines[2].find("<h2>Error occured:") != 1:
                s = re.findall(r"Error occured:([^<]+)", lines[2])[0]
                s = urllib.unquote(s)
                raise Datatypes.EUtilsError(s)
            raise Datatypes.EUtilsError("Unknown error:\n" +
                                        reseekfile.read(1000))

        # On error, fetch can return a valid XML document, but not one
        # which matches the DTD.  Rather than change the DTD (which is
        # pubmed_020114.dtd) I'll check it here to raise the error.
        # XXX HACK!
        if lines[2] == "<pmFetchResult>":
            # <pmFetchResult>
            # \t<ERROR>Empty id list - nothing todo</ERROR>
            # </pmFetchResult>
            s = "Unable to parse pmFetchResult error message"
            if len(lines) > 4:
                s = re.findall(r"<ERROR>([^>]+)</ERROR>", lines[3])[0]
            raise Datatypes.EUtilsError(s)

        # This happens when you choose a database which doesn't exist
        # Are there other reasons?  Probably yes, if you choose
        # other illegal parameters.
        if lines[0].startswith("<!doctype"):
            raise Datatypes.EUtilsError("Parameter not allowed")

    if force_encoding and lines[0].startswith('<?xml version="1.0"?>'):
        # Doesn't use an encoding, which means the XML is supposed
        # to be in UTF-8 encoding.  However, it seems NCBI uses
        # Latin-1 so we need to translate the Latin-1 input to
        # UTF-8 output else the XML parsers will fail for non-ASCII
        # characters.
        reseekfile = codecs.EncodedFile(reseekfile, "utf-8", "iso-8859-1")

    return reseekfile

##############################

def parse_search(infile, webenv_ref = [None]):
    # Need to pull out the webenv from the input stream
    infile = _check_for_bad_input_stream(infile)
    
    xml_parser = UsePOMParser("eSearch_020511")
    pom = xml_parser.parse_using_dtd(infile)

    errmsg, errors, warnings = _check_for_errors(pom)
    #        ErrorList      (PhraseNotFound*,FieldNotFound*)>
    #        WarningList    (PhraseIgnored*,
    #                            QuotedPhraseNotFound*,
    #                            OutputMessage*)>

    # If it's only "PhraseNotFound" erros, with an
    # OutputMessage of "No items found." then personally
    # think that should be considered the same as a search
    # which returned no results.

    # Set things up for an empty match
    webenv = None
    query_key = None
    count = 0
    retmax = 0
    retstart = 0
    ids = []
    translation_set = {}
    expression = None

    nothing_matched = 0
    if errmsg == "Can't run executor":
        # Check that the error list only contains PhraseNotFound terms
        flg = 1
        for x in errors:
            if x.category != "PhraseNotFound":
                flg = 0
                break
        if flg:
            # Okay, only PhraseNotFound.  Make sure there is
            # only one OutputMessage, with the text "No items found."
            # (Eg, an OutputMessage of 'Query syntax error.' means
            # there was a real problem.)
            msgs = [x for x in warnings if x.category == "OutputMessage"]
            if len(msgs) == 1 and msgs[0].text == "No items found.":
                nothing_matched = 1

        if not nothing_matched:
            # This is an error
            raise Datatypes.EUtilsSearchError(errmsg,
                                              errors,
                                              warnings)
            
    # In other words, check if something matched
    if not nothing_matched:
        ## Get WebEnv, if it exists
        if pom.get_element("WebEnv") is not None:
            s = pom["WebEnv"].tostring()
            webenv = urllib.unquote(s)
            # ONLY change webenv_ref if there's a new one
            webenv_ref[0] = webenv

        # Other simple fields
        if pom.get_element("QueryKey") is not None:
            query_key = pom["QueryKey"].tostring()

        count = int(pom["Count"].tostring())
        retmax = int(pom["RetMax"].tostring())
        retstart = int(pom["RetStart"].tostring())

        # The identifiers (if any)
        # NOTE: not a DBIds because the search result doesn't list the
        # database searched!
        ids = [x.tostring() for x in pom["IdList"].find_elements("Id")]

        # TranslationSet
        translation_set = {}
        for ele in pom["TranslationSet"]:
            translation_set[urllib.unquote_plus(ele["From"].tostring())] = \
                            urllib.unquote_plus(ele["To"].tostring())

        # Convert the RPN TranslationStack into an Expression
        stack = []
        try:
            translation_stack = pom["TranslationStack"]
        except IndexError:
            translation_stack = []
        for ele in translation_stack:
            if ele.__class__.__name__ == "TermSet":
                stack.append(Datatypes.Term(
                    term = urllib.unquote_plus(ele["Term"].tostring()),
                    field = urllib.unquote_plus(ele["Field"].tostring()),
                    count = int(ele["Count"].tostring()),
                    explode = ele["Explode"].tostring()))
            elif ele.__class__.__name__ == "OP":
                s = ele.tostring().strip()
                if s == "AND":
                    stack[-2:] = [stack[-2] & stack[-1]]
                elif s == "OR":
                    stack[-2:] = [stack[-2] | stack[-1]]
                elif s == "RANGE":
                    stack[-2:] = [Datatypes.Range(stack[-2], stack[-1])]
                elif s == "NOT":
                    stack[-2:] = [Datatypes.Not(stack[-2], stack[-1])]
                elif s == "GROUP":
                    # GROUP doesn't appear to do any more than put an extra
                    # parenthesis around ANDs and ORs -- can't find any
                    # specific documentation on its role
                    # So right now it is redundant and just ignore it
                    pass
                else:
                    raise TypeError("Unknown OP code: %r" % (s,))
            else:
                raise TypeError("Unknown TranslationStack element: %r" %
                                (ele.__class__.__name__,))
        
        # hack -- it appears as if the translation stack is sometimes missing
        # an AND at the end, which I guess is supposed to be implicit. For
        # instance, doing a text word search plus date range leaves off a
        # trailing and to link the final elements.
        if len(stack) == 2:
            stack[-2:] = [stack[-2] & stack[-1]]

        if len(stack) > 1:
            raise TypeError("Incomplete TranslationStack: %r" % stack)
        elif not stack:
            stack = [None]
        
        expression = stack[0]

    # Return either our synthesized query or 
    search_result = Datatypes.SearchResult(count, retmax, retstart, ids,
                                           translation_set, expression,
                                           webenv, query_key, errors,
                                           warnings, time.time())

    return search_result


###########################

def parse_post(infile, webenv_ref):
    # It doesn't look like I need check for a bad input stream
    # since I can only generate two types of error messages
    
    # ePost_020511.dtd
    xml_parser = UsePOMParser("ePost_020511")
    pom = xml_parser.parse_using_dtd(infile)

    # If there was an ERROR, raise it now
    errmsg, errors, warnings = _check_for_errors(pom)
    if errmsg is not None:
        raise Datatypes.EUtilsError(errmsg)

    # Get any invalid identifies
    invalid_ids = [x.tostring() for x in pom.get("InvalidIdList", [])]

    # Otherwise, get the WebEnv string
    s = pom["WebEnv"].tostring()
    webenv = urllib.unquote(s)
    webenv_ref[0] = webenv

    query_key = pom["QueryKey"].tostring()

    return Datatypes.PostResult(webenv, query_key, invalid_ids, time.time())

###############################

    
#  PubDate:  '2000 Feb 1'  or  '1975 Jun' or '1995'
#  BLAH! PubMed 8318652 also has "1993 May-Jun" for 
_pubdate_format1 = re.compile(
    r"(?P<year>\d{4})( (?P<month>[A-Za-z]{3})( (?P<day>\d+))?)?$")
_pubdate_format2 = re.compile(
    r"(?P<year>\d{4}) (?P<month1>[A-Za-z]{3})-(?P<month2>[A-Za-z]{3})")

_month_names_to_number = {
    None: 1,
    "Jan":  1,
    "Feb":  2,
    "Mar":  3,
    "Apr":  4,
    "May":  5,
    "Jun":  6,
    "Jul":  7,
    "Aug":  8,
    "Sep":  9,
    "Oct": 10,
    "Nov": 11,
    "Dec": 12,
    }
    
    

# Ignoring the hour and minute parts -- they seem to be either
# midnight or 09:00 and since I don't know the timezone it seems
# rather pointless
#  EntrezDate: 2000/02/17 09:00
_entrezdate_format = re.compile(r"(?P<year>\d+)/(?P<month>\d+)/(?P<day>\d+)")

# This may not be the right way to do this.
# Perhaps should keep the string and only translate upon request
# to a given time format?
def convert_summary_Date(x):
    return convert_summary_Date_string(x.tostring())

def convert_summary_Date_string(s):
    # Can be in one of several different formats
    m = _pubdate_format1.match(s)
    if m is not None:
        # 2000 Feb 15
        d = {}
        d["year"] = int(m.group("year"))
        d["month"] = _month_names_to_number[m.group("month")]
        try:
            d["day"] = int(m.group("day"))
        except TypeError:  # if this is None
            d["day"] = 1
        return Datatypes.Date(**d)
    m = _pubdate_format2.match(s)
    if m is not None:
        # 1993 May-Jun
        d = {}
        d["year"] = int(m.group("year"))
        d["month"] = _month_names_to_number[m.group("month1")]
        d["day"] = 1
        return Datatypes.Date(**d)

    m = _entrezdate_format.match(s)
    if m is not None:
        return Datatypes.Date(year = int(m.group("year"),),
                              month = int(m.group("month")),
                              day = int(m.group("day")))
                                
    raise TypeError("Unknown date format: %s" % (s,))

def unescape_entities(s):
    if "&" not in s:
        return unicode(s)
    
    terms = []
    i = 0
    defs = entitydefs
    for m in _entity_pattern.finditer(s):
        terms.append(s[i:m.start()])
        try:
            terms.append(defs[m.group(1)])
        except KeyError:
            terms.append(m.group(0))
        i = m.end()
    terms.append(s[i:])
    return "".join(terms)

def convert_summary_String(x):
    # The text may have HTML entity definitions .. convert as needed
    #
    # XXX Is this correct?  Most other characters are properly
    # encoded.  This may mean that that data provider messed up and
    # sent data in the wrong format.
    return unescape_entities(x.tostring())

def convert_summary_Integer(x):
    return int(x.tostring())

def convert_summary_Unknown(x):
    return x.tostring()

def convert_summary_List(x):
    # XXX I'm not doing this as a list.. Should I?
    return convert_summary_Items(x.find_elements("Item"))

def convert_summary_Items(x):
    d = MultiDict.OrderedMultiDict()
    for item in x:
        name = item.Name
        if name in d:
            print "Found multiple Items named %r!" % (name,)
        d[name] = summary_type_parser_table[item.Type](item)
    return d

summary_type_parser_table = {
    "String": convert_summary_String,
    "Integer": convert_summary_Integer,
    "Unknown": convert_summary_Unknown,
    "Date": convert_summary_Date,
    "List": convert_summary_List,
    }

def parse_summary_xml(infile):
    infile = _check_for_bad_input_stream(infile)
    xml_parser = UsePOMParser("eSummary_020511")
    pom = xml_parser.parse_using_dtd(infile)
    errmsg, errors, warnings = _check_for_errors(pom)
    if errmsg is not None:
        raise Datatypes.EUtilsError(errmsg)

    results = []
    for docsum in pom:
        id = docsum["Id"].tostring()
        d = convert_summary_Items(docsum.find_elements("Item"))
        results.append(Datatypes.Summary(id, d))
    
    return results

###############################

# XML
def parse_fetch_publication_xml(infile):
    infile = _check_for_bad_input_stream(infile, force_encoding = 0)
    xml_parser = UsePOMParser("pubmed_020114")
    return xml_parser.parse_using_dtd(infile)

def parse_fetch_sequence_xml(infile):
    raise NotImplementedError

# Identifer list ("\n" separated)
# Useful for "uilist", "acc", and a few others
def parse_fetch_identifiers(infile):
    infile = _check_for_bad_input_stream(infile)
    return [x.strip() for x in infile.readlines() if x != "\n"]

###############################
def _check_for_link_errors(pom):
    if not pom.has_key("LinkSet"):
        if pom.has_key("ERROR"):
            raise Datatypes.EUtilsError(pom["ERROR"].tostring())
        raise Datatypes.EUtilsError("Server failed to process request")
    if len(pom.find_elements("LinkSet")) != 1:
        raise AssertionError(
            "Did not expect to find more than one LinkSet in the XML")
    linkset = pom["LinkSet"]
    if linkset.has_key("ERROR"):
        raise Datatypes.EUtilsError(linkset["ERROR"].tostring())

def _parse_link(infile):
    #infile = _check_for_bad_input_stream(infile)
    # Need this, as seen in
    #  http://www.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&cmd=llinks&db=pubmed&id=10611131%2C12085853
    # which has an a-with-umlaut in Latin-1 encoding
    
    infile = codecs.EncodedFile(infile, "utf-8", "iso-8859-1")
    xml_parser = UsePOMParser("eLink_020511")
    pom = xml_parser.parse_using_dtd(infile)
    _check_for_link_errors(pom)
    return pom

def parse_neighbor_links(infile):
    pom = _parse_link(infile)
    pom_linkset = pom["LinkSet"]
    dbfrom = pom_linkset["DbFrom"].tostring().lower()
    idlist = [x.tostring() for x in pom_linkset["IdList"].find_elements("Id")]
    
    linksetdbs = MultiDict.OrderedMultiDict()
    for pom_linksetdb in pom_linkset.find_elements("LinkSetDb"):
        if pom_linksetdb.has_key("ERROR"):
            raise Datatypes.EUtilsError(pom_linksetdb["ERROR"].tostring())
        dbto = pom_linksetdb["DbTo"].tostring().lower()
        linkname = pom_linksetdb["LinkName"].tostring()
        links = []
        for pom_link in pom_linksetdb.find_elements("Link"):
            score = pom_link.get("Score")
            if score is not None:
                score = int(score.tostring())
            links.append(Datatypes.Link(pom_link["Id"].tostring(), score))
        linksetdbs[linkname] = Datatypes.LinkSetDb(dbto, linkname, links)
    return Datatypes.NeighborLinkSet(Datatypes.DBIds(dbfrom.lower(), idlist),
                                     linksetdbs)


def parse_lcheck(infile):
    pom = _parse_link(infile)
    pom_linkset = pom["LinkSet"]
    dbfrom = pom_linkset["DbFrom"].tostring().lower()
    idchecks = []
    for ele in pom_linkset["IdCheckList"].find_elements("Id"):
        has_linkout = getattr(ele, "HasLinkOut", "N")
        has_linkout = {"Y": 1}.get(has_linkout, 0)
        has_neighbor = getattr(ele, "HasNeighbor", "N")
        has_neighbor = {"Y": 1}.get(has_neighbor, 0)
        idchecks.append(Datatypes.IdCheck(ele.tostring(),
                                          has_linkout,
                                          has_neighbor))
    return Datatypes.CheckLinkSet(dbfrom, idchecks)

parse_ncheck = parse_lcheck

def _get_opt_string(ele, name):
    x = ele.get(name)
    if x is None:
        return None
    s = x.tostring()
    if not s:
        return None
    return s

def parse_llinks(infile):
    pom = _parse_link(infile)
    pom_linkset = pom["LinkSet"]
    dbfrom = pom_linkset["DbFrom"].tostring().lower()
    idurlsets = []
    for ele in pom_linkset["IdUrlList"].find_elements("IdUrlSet"):
        id = ele["Id"].tostring()
        objurls = []
        for pom_objurl in ele.find_elements("ObjUrl"):
            url = _get_opt_string(pom_objurl, "Url")
            linkname = _get_opt_string(pom_objurl, "LinkName")
            subject_types = [x.tostring() for x in
                                   pom_objurl.find_elements("SubjectType")]
            attributes = [s.tostring() for s in pom_objurl.find_elements("Attribute")]

            pom_provider = pom_objurl["Provider"]
            provider_name = pom_provider["Name"].tostring()
            provider_name_abbr = pom_provider["NameAbbr"].tostring()
            provider_id = pom_provider["Id"].tostring()
            provider_url = _get_opt_string(pom_provider, "Url")
            provider_icon_url = _get_opt_string(pom_provider, "IconUrl")
                
            provider = Datatypes.Provider(provider_name,
                                          provider_name_abbr,
                                          provider_id,
                                          provider_url,
                                          provider_icon_url)
            objurl = Datatypes.ObjUrl(subject_types, provider,
                                      linkname, url, attributes)
            objurls.append(objurl)
            
        idurlsets.append(Datatypes.IdUrlSet(id, objurls))
                         
    return Datatypes.LinksLinkSet(dbfrom, idurlsets)

parse_prlinks = parse_llinks

def parse_link_xml(infile):
    infile = _check_for_bad_input_stream(infile)
    xml_parser = UsePOMParser("eLink_020511")
    return xml_parser.parse_using_dtd(infile)
