import unittest, cStringIO, re
from Bio.EUtils import parse, Datatypes, MultiDict

#  ======================  Test parsing Search results

search1 = """\
<?xml version="1.0"?>
<!DOCTYPE eSearchResult PUBLIC "-//NLM//DTD eSearchResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eSearch_020511.dtd">
<eSearchResult>
        <Count>30</Count>
        <RetMax>20</RetMax>
        <RetStart>0</RetStart>
        <IdList>
                <Id>12269969</Id>
                <Id>10912950</Id>
                <Id>9454215</Id>
                <Id>9454196</Id>
                <Id>9454186</Id>
                <Id>9390282</Id>
                <Id>9303476</Id>
                <Id>9300720</Id>
                <Id>8763495</Id>
                <Id>8744570</Id>
                <Id>8566008</Id>
                <Id>7648552</Id>
                <Id>7858482</Id>
                <Id>7909813</Id>
                <Id>8318652</Id>
                <Id>1644686</Id>
                <Id>1641862</Id>
                <Id>1599914</Id>
                <Id>1310101</Id>
                <Id>2007738</Id>
        </IdList>
        <TranslationSet>
        </TranslationSet>
        <TranslationStack>
                <TermSet>
                        <Term>dalke[All Fields]</Term>
                        <Field>All Fields</Field>
                        <Count>30</Count>
                        <Explode>Y</Explode>
                </TermSet>
        </TranslationStack>
</eSearchResult>
"""
search1_expected = Datatypes.SearchResult(
    30, 20, 0,
    ["12269969", "10912950", "9454215", "9454196", "9454186",
     "9390282", "9303476", "9300720", "8763495", "8744570",
     "8566008", "7648552", "7858482", "7909813", "8318652",
     "1644686", "1641862", "1599914", "1310101", "2007738"],
    {}, "dalke[All Fields]",
    None, None,
    [], [], None)

search1_terms = [("All Fields", 30, "Y")]


search2 = """\
<?xml version="1.0"?>
<!DOCTYPE eSearchResult PUBLIC "-//NLM//DTD eSearchResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eSearch_020511.dtd">
<eSearchResult>
        <Count>134325</Count>
        <RetMax>2</RetMax>
        <RetStart>7</RetStart>
        <IdList>
                <Id>9338476</Id>
                <Id>9338475</Id>
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
                        <Term>1995/02[EDAT]</Term>
                        <Field>EDAT</Field>
                        <Count>-1</Count>
                        <Explode>Y</Explode>
                </TermSet>
                <TermSet>
                        <Term>1997/10/24[EDAT]</Term>
                        <Field>EDAT</Field>
                        <Count>-1</Count>
                        <Explode>Y</Explode>
                </TermSet>
                <OP>RANGE</OP>
                <OP>AND</OP>
        </TranslationStack>
</eSearchResult>
"""
search2_expected = Datatypes.SearchResult(
    134325, 2, 7, ["9338476", "9338475"],
    {"cancer[All Fields]": '("neoplasms"[MeSH Terms] OR cancer[Text Word])'},
    '(("neoplasms"[MeSH Terms] OR cancer[Text Word]) AND '
                                    '1995/02:1997/10/24[EDAT])',
    None, None,
    [], [], None)

search2_terms = [("MeSH Terms", 1407151, "Y"),
                 ("Text Word", 382919, "Y"),
                 ("EDAT", -1, "Y"),
                 ("EDAT", -1, "Y")]
            
    
# Search: dalke BUTNOT VMD
search3 = """\
<?xml version="1.0"?>
<!DOCTYPE eSearchResult PUBLIC "-//NLM//DTD eSearchResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eSearch_020511.dtd">
<eSearchResult>
        <Count>28</Count>
        <RetMax>3</RetMax>
        <RetStart>0</RetStart>
        <QueryKey>1</QueryKey>
        <WebEnv>OPqbh%40H%3E%5EfG_A%3DEftI_JCJ%5CQKL%5D%3CIEueha%3CFgWB%3D_iJ%3Fg%60%3F%3EfAK%5C_</WebEnv>
        <IdList>
                <Id>12269969</Id>
                <Id>10912950</Id>
                <Id>9454215</Id>
        </IdList>
        <TranslationSet>
        </TranslationSet>
        <TranslationStack>
                <TermSet>
                        <Term>dalke[All Fields]</Term>
                        <Field>All Fields</Field>
                        <Count>30</Count>
                        <Explode>Y</Explode>
                </TermSet>
                <TermSet>
                        <Term>VMD[All Fields]</Term>
                        <Field>All Fields</Field>
                        <Count>68</Count>
                        <Explode>Y</Explode>
                </TermSet>
                <OP>NOT</OP>
        </TranslationStack>
</eSearchResult>
"""

search3_expected = Datatypes.SearchResult(
    28, 3, 0, ["12269969", "10912950", "9454215"],
    {}, "(dalke[All Fields] NOT VMD[All Fields])",
    "OPqbh@H>^fG_A=EftI_JCJ\\QKL]<IEueha<FgWB=_iJ?g`?>fAK\\_", "1",
    [], [], None)

search3_terms = [("All Fields", 30, "Y"),
                 ("All Fields", 68, "Y")]

class ParseSearchWebEnv(unittest.TestCase):
    def testSearchHistory(self):
        expected_webenv = ('OPqbh@H>^fG_A=EftI_JCJ\\QKL]<IEueha'
                           '<FgWB=_iJ?g`?>fAK\\_')

        webenv_ref = [None]
        search = parse.parse_search(cStringIO.StringIO(search3), webenv_ref)
        self.assertEquals(webenv_ref, [expected_webenv])
        self.assertEquals(search.webenv, expected_webenv)

    def testSearchNoHistory(self):
        webenv_ref = ["qwe"]
        search = parse.parse_search(cStringIO.StringIO(search2), webenv_ref)
        self.assertEquals(webenv_ref, ["qwe"])
        self.assertEquals(search.webenv, None)

## Various errors or problematical queries

# This was from an ill-formed query of: 'qwerty asdfg'
# instead of: "qwerty asdfg"
# Because apostrophes are allowed (eg, for "O'Henry")
# I judge this should be a valid result with no matches
# and *NOT* an error.
searcherr1 = """\
<eSearchResult>
        <ERROR>Can't run executor</ERROR>
        <ErrorList>
                <PhraseNotFound>'qwerty[All Fields]</PhraseNotFound>
                <PhraseNotFound>asdfg'[All Fields]</PhraseNotFound>
        </ErrorList>
        <WarningList>
                <OutputMessage>No items found.</OutputMessage>
        </WarningList>
</eSearchResult>
"""
searcherr1_expected = Datatypes.SearchResult(
    0, 0, 0, [], {}, None, None, None,
    [Datatypes.PhraseNotFound("'qwerty[All Fields]"),
         Datatypes.PhraseNotFound("asdfg'[All Fields]")],
    [Datatypes.OutputMessage("No items found.")],
    None)

searcherr1_terms = []
    

# This is a legitimate query!
# Search: qweqweqweqwe and qwerewtre
searcherr2 = """\
<eSearchResult>
        <ERROR>Can't run executor</ERROR>
        <ErrorList>
                <PhraseNotFound>qweqweqweqwe[All Fields]</PhraseNotFound>
                <PhraseNotFound>qwerewtre[All Fields]</PhraseNotFound>
        </ErrorList>
        <WarningList>
                <PhraseIgnored>and</PhraseIgnored>
                <OutputMessage>No items found.</OutputMessage>
        </WarningList>
</eSearchResult>
"""
searcherr2_expected = Datatypes.SearchResult(
    0, 0, 0, [], {}, None, None, None,
    [Datatypes.PhraseNotFound("qweqweqweqwe[All Fields]"),
              Datatypes.PhraseNotFound("qwerewtre[All Fields]")],
    [Datatypes.PhraseIgnored("and"),
              Datatypes.OutputMessage("No items found.")],
    None)

searcherr2_terms = []

# Search: qwertya[XYZ]
searcherr3 = """\
<eSearchResult>
        <ERROR>Can't run executor</ERROR>
        <ErrorList>
                <PhraseNotFound>qwertya[All Fields]</PhraseNotFound>
                <FieldNotFound>XYZ</FieldNotFound>
        </ErrorList>
        <WarningList>
                <OutputMessage>No items found.</OutputMessage>
        </WarningList>
</eSearchResult>
"""
searcherr3_expected = Datatypes.EUtilsSearchError(
    "Can't run executor",
    [Datatypes.PhraseNotFound("qwertya[All Fields]"),
     Datatypes.FieldNotFound("XYZ")],
    [Datatypes.OutputMessage("No items found.")])

searcherr3_terms = []

# Search: "
# (Yes, a single character, which is a double quote.)
searcherr4 = """\
<eSearchResult>
        <ERROR>Can't run executor</ERROR>
        <WarningList>
                <OutputMessage>Your request produced no results</OutputMessage>
                <OutputMessage>Query syntax error.</OutputMessage>
        </WarningList>
</eSearchResult>
"""
searcherr4_expected = Datatypes.EUtilsSearchError(
    "Can't run executor",
    [],
    [Datatypes.OutputMessage("Your request produced no results"),
     Datatypes.OutputMessage("Query syntax error.")])

searcherr4_terms = []

# Search: NOT VMD
searcherr5 = """\
<eSearchResult>
        <ERROR>Can't run executor</ERROR>
        <WarningList>
                <OutputMessage>Syntax error in query: improper beginning of expression</OutputMessage>
        </WarningList>
</eSearchResult>
"""
searcherr5_expected = Datatypes.EUtilsSearchError(
    "Can't run executor", [],
    [Datatypes.OutputMessage(
           "Syntax error in query: improper beginning of expression")])

searcherr5_terms = []

searcherr6_error = """\
<eSearchResult>
        <ERROR>Can't run executor</ERROR>
        <WarningList>
                <OutputMessage>Wrong subquery 2.</OutputMessage>
                <OutputMessage>No items found.</OutputMessage>
        </WarningList>
</eSearchResult>
"""
searcherr6_expected = Datatypes.EUtilsSearchError(
    "Can't run executor", [],
    [Datatypes.OutputMessage("Wrong subquery2."),
     Datatypes.OutputMessage("No items found."),])
    
searcherr6_terms = []

searchwarn1 = """\
<?xml version="1.0"?>
<!DOCTYPE eSearchResult PUBLIC "-//NLM//DTD eSearchResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eSearch_020511.dtd">
<eSearchResult>
        <Count>0</Count>
        <RetMax>0</RetMax>
        <RetStart>0</RetStart>
        <IdList>
        </IdList>
        <TranslationSet>
        </TranslationSet>
        <TranslationStack>
                <TermSet>
                        <Term>qwerty[All Fields]</Term>
                        <Field>All Fields</Field>
                        <Count>9</Count>
                        <Explode>Y</Explode>
                </TermSet>
                <TermSet>
                        <Term>asdfg[All Fields]</Term>
                        <Field>All Fields</Field>
                        <Count>1</Count>
                        <Explode>Y</Explode>
                </TermSet>
                <OP>AND</OP>
        </TranslationStack>
        <WarningList>
                <QuotedPhraseNotFound>"qwerty asdfg"[All Fields]</QuotedPhraseNotFound>
        </WarningList>
</eSearchResult>
"""
searchwarn1_expected = Datatypes.SearchResult(
    0, 0, 0, [], {},
    "(qwerty[All Fields] AND asdfg[All Fields])",  # XXX should not do this
    None, None,
    [],
    [Datatypes.QuotedPhraseNotFound('"qwerty asdfg"[All Fields]')],
    None)

searchwarn1_terms = [("All Fields", 9, "Y"),
                     ("All Fields", 1, "Y")]

# Search: cancer and aswfgh
# with retmax = 2
searchwarn2 = """\
<?xml version="1.0"?>
<!DOCTYPE eSearchResult PUBLIC "-//NLM//DTD eSearchResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eSearch_020511.dtd">
<eSearchResult>
        <Count>1456738</Count>
        <RetMax>20</RetMax>
        <RetStart>0</RetStart>
        <IdList>
                <Id>12508380</Id>
                <Id>12508360</Id>
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
                        <Count>1407473</Count>
                        <Explode>Y</Explode>
                </TermSet>
                <TermSet>
                        <Term>cancer[Text Word]</Term>
                        <Field>Text Word</Field>
                        <Count>383244</Count>
                        <Explode>Y</Explode>
                </TermSet>
                <OP>OR</OP>
        </TranslationStack>
        <ErrorList>
                <PhraseNotFound>aswfgh[All Fields]</PhraseNotFound>
        </ErrorList>
        <WarningList>
                <PhraseIgnored>and</PhraseIgnored>
        </WarningList>
</eSearchResult>
"""
searchwarn2_expected = Datatypes.SearchResult(
    1456738, 20, 0, ["12508380", "12508360"],
    {'cancer[All Fields]': '("neoplasms"[MeSH Terms] OR cancer[Text Word])'},
    '("neoplasms"[MeSH Terms] OR cancer[Text Word])',
    None, None,
    [Datatypes.PhraseNotFound("aswfgh[All Fields]")],
    [Datatypes.PhraseIgnored("and")],
    None)

searchwarn2_terms = [("MeSH Terms", 1407473, "Y"),
                     ("Text Word", 383244, "Y")]

# This was a search for: cancer and aswfgh OR ("qwertyu pijngq")
# with retmax set to 2
searchwarn3 = """\
<?xml version="1.0"?>
<!DOCTYPE eSearchResult PUBLIC "-//NLM//DTD eSearchResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eSearch_020511.dtd">
<eSearchResult>
        <Count>1456738</Count>
        <RetMax>2</RetMax>
        <RetStart>0</RetStart>
        <IdList>
                <Id>12508380</Id>
                <Id>12508360</Id>
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
                        <Count>1407473</Count>
                        <Explode>Y</Explode>
                </TermSet>
                <TermSet>
                        <Term>cancer[Text Word]</Term>
                        <Field>Text Word</Field>
                        <Count>383244</Count>
                        <Explode>Y</Explode>
                </TermSet>
                <OP>OR</OP>
                <TermSet>
                        <Term>"qwertyu pijngq"[All Fields]</Term>
                        <Field>All Fields</Field>
                        <Count>0</Count>
                        <Explode>Y</Explode>
                </TermSet>
                <OP>OR</OP>
        </TranslationStack>
        <ErrorList>
                <PhraseNotFound>aswfgh[All Fields]</PhraseNotFound>
                <PhraseNotFound>qwertyu[All Fields]</PhraseNotFound>
                <PhraseNotFound>pijngq[All Fields]</PhraseNotFound>
        </ErrorList>
        <WarningList>
                <PhraseIgnored>and</PhraseIgnored>
                <QuotedPhraseNotFound>"qwertyu pijngq"[All Fields]</QuotedPhraseNotFound>
        </WarningList>
</eSearchResult>
"""
searchwarn3_expected = Datatypes.SearchResult(
    1456738, 2, 0, ["12508380", "12508360"],
    {'cancer[All Fields]': '("neoplasms"[MeSH Terms] OR cancer[Text Word])'},
    '(("neoplasms"[MeSH Terms] OR cancer[Text Word])'
                           ' OR "qwertyu pijngq"[All Fields])',
     None, None,
    [Datatypes.PhraseNotFound("aswfgh[All Fields]"),
     Datatypes.PhraseNotFound("qwertyu[All Fields]"),
     Datatypes.PhraseNotFound("pijngq[All Fields]")],
    [Datatypes.PhraseIgnored("and"),
     Datatypes.QuotedPhraseNotFound('"qwertyu pijngq"[All Fields]')],
    None)

searchwarn3_terms = [("MeSH Terms", 1407473, "Y"),
                     ("Text Word", 383244, "Y"),
                     ("All Fields", 0, "Y")]

# Search for: dalke OR ("qwertyu pijngq")
# retmax = 0
searchwarn4 = """\
<?xml version="1.0"?>
<!DOCTYPE eSearchResult PUBLIC "-//NLM//DTD eSearchResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eSearch_020511.dtd">
<eSearchResult>
        <Count>30</Count>
        <RetMax>0</RetMax>
        <RetStart>0</RetStart>
        <IdList>
        </IdList>
        <TranslationSet>
        </TranslationSet>
        <TranslationStack>
                <TermSet>
                        <Term>dalke[All Fields]</Term>
                        <Field>All Fields</Field>
                        <Count>30</Count>
                        <Explode>Y</Explode>
                </TermSet>
                <TermSet>
                        <Term>"qwertyu pijngq"[All Fields]</Term>
                        <Field>All Fields</Field>
                        <Count>0</Count>
                        <Explode>Y</Explode>
                </TermSet>
                <OP>OR</OP>
        </TranslationStack>
        <ErrorList>
                <PhraseNotFound>qwertyu[All Fields]</PhraseNotFound>
                <PhraseNotFound>pijngq[All Fields]</PhraseNotFound>
        </ErrorList>
        <WarningList>
                <QuotedPhraseNotFound>"qwertyu pijngq"[All Fields]</QuotedPhraseNotFound>
        </WarningList>
</eSearchResult>
"""
searchwarn4_expected = Datatypes.SearchResult(
    30, 0, 0, [], {},
    '(dalke[All Fields] OR "qwertyu pijngq"[All Fields])',
    None, None,
    [Datatypes.PhraseNotFound("qwertyu[All Fields]"),
     Datatypes.PhraseNotFound("pijngq[All Fields]")],
    [Datatypes.QuotedPhraseNotFound('"qwertyu pijngq"[All Fields]')],
    None)

searchwarn4_terms = [("All Fields", 30, "Y"),
                     ("All Fields", 0, "Y")]


# Search: dalke OR ("qwerty asdfg")
# retmax = 0
searchwarn5 = """\
<?xml version="1.0"?>
<!DOCTYPE eSearchResult PUBLIC "-//NLM//DTD eSearchResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eSearch_020511.dtd">
<eSearchResult>
        <Count>30</Count>
        <RetMax>0</RetMax>
        <RetStart>0</RetStart>
        <IdList>
        </IdList>
        <TranslationSet>
        </TranslationSet>
        <TranslationStack>
                <TermSet>
                        <Term>dalke[All Fields]</Term>
                        <Field>All Fields</Field>
                        <Count>30</Count>
                        <Explode>Y</Explode>
                </TermSet>
                <TermSet>
                        <Term>qwerty[All Fields]</Term>
                        <Field>All Fields</Field>
                        <Count>9</Count>
                        <Explode>Y</Explode>
                </TermSet>
                <TermSet>
                        <Term>asdfg[All Fields]</Term>
                        <Field>All Fields</Field>
                        <Count>1</Count>
                        <Explode>Y</Explode>
                </TermSet>
                <OP>AND</OP>
                <OP>OR</OP>
        </TranslationStack>
        <WarningList>
                <QuotedPhraseNotFound>"qwerty asdfg"[All Fields]</QuotedPhraseNotFound>
        </WarningList>
</eSearchResult>
"""
searchwarn5_expected = Datatypes.SearchResult(
    30, 0, 0, [], {},
    "(dalke[All Fields] OR (qwerty[All Fields] AND asdfg[All Fields]))",
    None, None,
    [],
    [Datatypes.QuotedPhraseNotFound('"qwerty asdfg"[All Fields]')],
    None)

searchwarn5_terms = [("All Fields", 30, "Y"),
                     ("All Fields", 9, "Y"),
                     ("All Fields", 1, "Y")]

def _get_expression_terms2(expr):
    if isinstance(expr, Datatypes.Term):
        yield expr
    else:
        for x in _get_expression_terms2(expr.left):
            yield x
        for x in _get_expression_terms2(expr.right):
            yield x

def _get_expression_terms(expr):
    if expr is None:
        return []
    return _get_expression_terms2(expr)


class CompareSearchAgainstExpected(unittest.TestCase):
    def _compare_search_error(self, got, expected):
        self.assertEquals(got.errmsg, expected.errmsg)
        self._compare_problems(got.errors, expected.errors)
        self._compare_problems(got.warnings, expected.warnings)
    
    def _compare_problems(self, got_list, expected_list):
        self.assertEquals(len(got_list), len(expected_list))
        for got, expected in zip(got_list, expected_list):
            # They must be of the same class and have the same text
            self.assertEquals(got.category, expected.category)
            self.assertEquals(got.text, expected.text)

    def _compare_result(self, got, expected):
        for name in ("count", "retmax", "retstart", "ids",
                     "translation_set", "webenv", "query_key"):
            if getattr(got, name) != getattr(expected, name):
                raise AssertionError("Different %r: %r and %r" %
                                     (name, getattr(got, name),
                                      getattr(expected, name)))
        # the 'expected' expression is a simple string (HACK!) or None
        if expected.expression is None:
            self.assertEquals(got.expression, expected.expression)
        else:
            self.assertEquals(str(got.expression), expected.expression)

        self._compare_problems(got.errors, expected.errors)
        self._compare_problems(got.warnings, expected.warnings)

    def _compare_terms(self, got, expected):
        got_terms = [x for x in _get_expression_terms(got.expression)]
        self.assertEquals(len(got_terms), len(expected))
        for got_term, expected_term in zip(got_terms, expected):
            self.assertEquals(got_term.field, expected_term[0])
            self.assertEquals(got_term.count, expected_term[1])
            self.assertEquals(got_term.explode, expected_term[2])
    
    def _test(self, s, expected, terms):
        if isinstance(expected, Datatypes.SearchResult):
            got = parse.parse_search(cStringIO.StringIO(s))
            self._compare_result(got, expected)
            self._compare_terms(got, terms)
        else:
            self.assertEquals(terms, [])
            try:
                parse.parse_search(cStringIO.StringIO(s))
            except Datatypes.EUtilsSearchError, got:
                self._compare_search_error(got, expected)
            else:
                raise AssertionError("Expecting an exception")

def _add_search_test_methods():
    g = globals()
    for k, v in g.items():
        if (isinstance(v, type("")) and
            k.startswith("search")):
            i = 0
            if k + "_expected" in g:
                i = i + 1
            if k + "_terms" in g:
                i = i + 1
            if i == 0:
                break
            if i != 2:
                raise AssertionError(
                    "have %s without both %s_expected and %s_terms" %
                    (k, k, k))

            code = "def test_%s(self): self._test(%s, %s_expected, %s_terms)"\
                     % (k, k, k, k)
            d = {}
            exec code in globals(), d
            name = "test_" + k
            setattr(CompareSearchAgainstExpected, name, d[name])

_add_search_test_methods()


##  ====================== Test parsing Post results

posttest1 = """\
<?xml version="1.0"?>
<!DOCTYPE ePostResult PUBLIC "-//NLM//DTD ePostResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/ePost_020511.dtd">
<ePostResult>
        <QueryKey>1</QueryKey>
        <WebEnv>PB%60qj%3FGbNY_%3FbFJ%3FX%3CBE%3D%5Cb%5Cg%60%3DIfb%5Cab%3EI%40HI%3D%5E%3EdA%3DtCG%3EH%3DF%5Cb</WebEnv>
</ePostResult>
"""

posttest2 = """\
<?xml version="1.0"?>
<!DOCTYPE ePostResult PUBLIC "-//NLM//DTD ePostResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/ePost_020511.dtd">
<ePostResult>
        <QueryKey>2</QueryKey>
        <WebEnv>%3FPUV%3C%3F%3C%3FBNe%5C%3C%3C%5Ce%5DB%3F%3F%3CeGAAyFFDfs%5CgeD%40gfECDieGnhDCfAG%5C%60</WebEnv>
</ePostResult>
"""

posttest3 = """\
<?xml version="1.0"?>
<!DOCTYPE ePostResult PUBLIC "-//NLM//DTD ePostResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/ePost_020511.dtd">
<ePostResult>
  <InvalidIdList>
  <Id>14000000</Id> 
  <Id>15000000</Id> 
  </InvalidIdList>
  <QueryKey>1</QueryKey> 
  <WebEnv>pxTK%60eD%3CDkeIB%5EcFX%3CS%40KCJrg%3DgGgIJ_IEB%60%3DLD%3DDIFKtH%3E%40</WebEnv>
</ePostResult>
"""

class ParsePostBase:
    def _parse(self, s, webenv_ref = [None]):
        return parse.parse_post(cStringIO.StringIO(s), webenv_ref)
    def _parse_fail(self, s, expect, webenv_ref = [None]):
        try:
            parse.parse_post(cStringIO.StringIO(s), webenv_ref)
            raise AssertionError("Should have raised a parse error")
        except Datatypes.EUtilsError, err:
            self.assertEquals(str(err), expect)

class ParsePost(ParsePostBase, unittest.TestCase):
    def testParsePost(self):
        webenv_ref = [None]
        post1 = self._parse(posttest1, webenv_ref)
        assert post1.query_key == "1"
        assert webenv_ref[0] == post1.webenv == \
               u'PB`qj?GbNY_?bFJ?X<BE=\\b\\g`=Ifb\\ab>I@HI=^>dA=tCG>H=F\\b'

        post2 = self._parse(posttest2, webenv_ref)
        assert post2.query_key == "2"
        assert webenv_ref[0] == post2.webenv == \
               u'?PUV<?<?BNe\\<<\\e]B??<eGAAyFFDfs\\geD@gfECDieGnhDCfAG\\`'

    def testParsePostWithInvalidIds(self):
        webenv_ref = [None]
        post = self._parse(posttest3, webenv_ref)
        self.assertEquals(post.invalid_ids, ["14000000", "15000000"])

errpost1 = """\
<ePostResult>
        <ERROR>Wrong DB name</ERROR>
</ePostResult>
"""

class ParsePostError(ParsePostBase, unittest.TestCase):
    def testParsePostError(self):
        self._parse_fail(errpost1, "Wrong DB name")

##  ====================== Test parsing Summary results

summary1 = """\
<?xml version="1.0"?>
<!DOCTYPE eSummaryResult PUBLIC "-//NLM//DTD eSummaryResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eSummary_020511.dtd">
<eSummaryResult>
<DocSum>
        <Id>1</Id>
        <Item Name="PubDate" Type="Date">1975 Jun</Item>
        <Item Name="Source" Type="String">Biochem Med</Item>
        <Item Name="Authors" Type="String">Makar AB, McMartin KE, Palese M, Tephly TR</Item>
        <Item Name="Title" Type="String">Formate assay in body fluids: application in methanol poisoning.</Item>
        <Item Name="Volume" Type="String">13</Item>
        <Item Name="Pages" Type="String">117-26</Item>
        <Item Name="EntrezDate" Type="Date">1975/06/01 00:00</Item>
        <Item Name="PubMedId" Type="Integer">1</Item>
        <Item Name="MedlineId" Type="Integer">76061540</Item>
        <Item Name="Lang" Type="String">English</Item>
        <Item Name="PubType" Type="String"></Item>
        <Item Name="RecordStatus" Type="String">PubMed - indexed for MEDLINE</Item>
        <Item Name="Issue" Type="String">2</Item>
        <Item Name="SO" Type="String">1975 Jun;13(2):117-26</Item>
        <Item Name="DOI" Type="String"></Item>
        <Item Name="JTA" Type="String">9YW</Item>
        <Item Name="ISSN" Type="String">0006-2944</Item>
        <Item Name="PubId" Type="String"></Item>
        <Item Name="PubStatus" Type="Integer">4</Item>
        <Item Name="Status" Type="Integer">4</Item>
        <Item Name="HasAbstract" Type="Integer">0</Item>
        <Item Name="ArticleIds" Type="List">
                <Item Name="PubMedId" Type="String">1</Item>
                <Item Name="MedlineUID" Type="String">76061540</Item>
        </Item>
</DocSum>
<DocSum>
        <Id>2</Id>
        <Item Name="PubDate" Type="Date">1975 Oct 27</Item>
        <Item Name="Source" Type="String">Biochem Biophys Res Commun</Item>
        <Item Name="Authors" Type="String">Bose KS, Sarma RH</Item>
        <Item Name="Title" Type="String">Delineation of the intimate details of the backbone conformation of pyridine nucleotide coenzymes in aqueous solution.</Item>
        <Item Name="Volume" Type="String">66</Item>
        <Item Name="Pages" Type="String">1173-9</Item>
        <Item Name="EntrezDate" Type="Date">1975/10/27 00:00</Item>
        <Item Name="PubMedId" Type="Integer">2</Item>
        <Item Name="MedlineId" Type="Integer">76061562</Item>
        <Item Name="Lang" Type="String">English</Item>
        <Item Name="PubType" Type="String"></Item>
        <Item Name="RecordStatus" Type="String">PubMed - indexed for MEDLINE</Item>
        <Item Name="Issue" Type="String">4</Item>
        <Item Name="SO" Type="String">1975 Oct 27;66(4):1173-9</Item>
        <Item Name="DOI" Type="String"></Item>
        <Item Name="JTA" Type="String">9Y8</Item>
        <Item Name="ISSN" Type="String">0006-291X</Item>
        <Item Name="PubId" Type="String"></Item>
        <Item Name="PubStatus" Type="Integer">4</Item>
        <Item Name="Status" Type="Integer">4</Item>
        <Item Name="HasAbstract" Type="Integer">0</Item>
        <Item Name="ArticleIds" Type="List">
                <Item Name="PubMedId" Type="String">2</Item>
                <Item Name="MedlineUID" Type="String">76061562</Item>
        </Item>
</DocSum>
</eSummaryResult>
"""

summary1_expected = [
    Datatypes.Summary("1", MultiDict.OrderedMultiDict( (
                       ("PubDate", Datatypes.Date(1975, 6, 1)),
                       ("Source", "Biochem Med"),
                       ("Authors", "Makar AB, McMartin KE, Palese M, Tephly TR"),
                       ("Title", "Formate assay in body fluids: application in methanol poisoning."),
                       ("Volume", "13"),
                       ("Pages", "117-26"),
                       ("EntrezDate", Datatypes.Date(1975, 6, 1)),
                       ("PubMedId", 1),
                       ("MedlineId", 76061540),
                       ("Lang", "English"),
                       ("PubType", ""),
                       ("RecordStatus", "PubMed - indexed for MEDLINE"),
                       ("Issue", "2"),
                       ("SO", "1975 Jun;13(2):117-26"),
                       ("DOI", ""),
                       ("JTA", "9YW"),
                       ("ISSN", "0006-2944"),
                       ("PubId", ""),
                       ("PubStatus", 4), # ??
                       ("Status", 4),  # ??
                       ("HasAbstract", 0),
                       ("ArticleIds", MultiDict.OrderedMultiDict( (
                                         ("PubMedId", "1"),
                                         ("MedlineUID", "76061540") )))))),
    Datatypes.Summary("2",
                      MultiDict.OrderedMultiDict( (
                       ("PubDate", Datatypes.Date(1975, 10, 27)),
                       ("Source", "Biochem Biophys Res Commun"),
                       ("Authors", "Bose KS, Sarma RH"),
                       ("Title", "Delineation of the intimate details of the backbone conformation of pyridine nucleotide coenzymes in aqueous solution."),
                       ("Volume", "66"),
                       ("Pages", "1173-9"),
                       ("EntrezDate", Datatypes.Date(1975, 10, 27)),
                       ("PubMedId", 2),
                       ("MedlineId", 76061562),
                       ("Lang", "English"),
                       ("PubType", ""),
                       ("RecordStatus", "PubMed - indexed for MEDLINE"),
                       ("Issue", "4"),
                       ("SO", "1975 Oct 27;66(4):1173-9"),
                       ("DOI", ""),
                       ("JTA", "9Y8"),
                       ("ISSN", "0006-291X"),
                       ("PubId", ""),
                       ("PubStatus", 4),
                       ("Status", 4),
                       ("HasAbstract", 0),
                       ("ArticleIds", MultiDict.OrderedMultiDict((
                                         ("PubMedId", "2"),
                                         ("MedlineUID", "76061562"))))))),
    ]


# This one has non-ASCII characters
summary2 = """\
<?xml version="1.0"?>
<!DOCTYPE eSummaryResult PUBLIC "-//NLM//DTD eSummaryResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eSummary_020511.dtd">
<eSummaryResult>
<DocSum>
        <Id>9178602</Id>
        <Item Name="PubDate" Type="Date">1997 Jul</Item>
        <Item Name="Source" Type="String">Microb Ecol</Item>
        <Item Name="Authors" Type="String">Jeppesen E, Erlandsen M, S&amp;oslash;ndergaard M</Item>
        <Item Name="Title" Type="String">Can Simple Empirical Equations Describe the Seasonal Dynamics of Bacterioplankton in Lakes: An Eight-Year Study in Shallow Hypertrophic and Biologically Highly Dynamic Lake S&amp;oslash;byg&amp;aring;rd, Denmark</Item>
        <Item Name="Volume" Type="String">34</Item>
        <Item Name="Pages" Type="String">11-26</Item>
        <Item Name="EntrezDate" Type="Date">1997/07/01 00:00</Item>
        <Item Name="PubMedId" Type="Integer">9178602</Item>
        <Item Name="MedlineId" Type="Integer">0</Item>
        <Item Name="Lang" Type="String">English</Item>
        <Item Name="PubType" Type="String"></Item>
        <Item Name="RecordStatus" Type="String">PubMed - as supplied by publisher</Item>
        <Item Name="Issue" Type="String">1</Item>
        <Item Name="SO" Type="String">1997 Jul;34(1):11-26</Item>
        <Item Name="DOI" Type="String"></Item>
        <Item Name="JTA" Type="String"></Item>
        <Item Name="ISSN" Type="String">0095-3628</Item>
        <Item Name="PubId" Type="String">MECOJKF96-29</Item>
        <Item Name="PubStatus" Type="Integer">4</Item>
        <Item Name="Status" Type="Integer">4</Item>
        <Item Name="HasAbstract" Type="Integer">0</Item>
        <Item Name="ArticleIds" Type="List">
                <Item Name="PII" Type="String">MECOJKF96-29</Item>
                <Item Name="PubMedId" Type="String">9178602</Item>
        </Item>
</DocSum>
</eSummaryResult>
"""
OSLASH = u"\N{LATIN SMALL LETTER O WITH STROKE}"
ARING = u"\N{LATIN SMALL LETTER A WITH RING ABOVE}"
summary2_expected = [
    Datatypes.Summary("9178602",
                      MultiDict.OrderedMultiDict( (
                             ("PubDate", Datatypes.Date(1997, 7, 1)),
                             ("Source", "Microb Ecol"),
                             ("Authors", ("Jeppesen E, Erlandsen M, S" +
                                   OSLASH + "ndergaard M")),
                             ("Title", 
    "Can Simple Empirical Equations Describe the Seasonal Dynamics " +
    "of Bacterioplankton in Lakes: An Eight-Year Study in Shallow " +
    "Hypertrophic and Biologically Highly Dynamic Lake " +
    "S" + OSLASH + "byg" + ARING + "rd, Denmark"),
                             ("Volume", "34"),
                             ("Pages", "11-26"),
                             ("EntrezDate", Datatypes.Date(1997, 7, 1)),
                             ("PubMedId", 9178602),
                             ("MedlineId", 0),
                             ("Lang", "English"),
                             ("PubType", ""),
                             ("RecordStatus", "PubMed - as supplied by publisher"),
                             ("Issue", "1"),
                             ("SO", "1997 Jul;34(1):11-26"),
                             ("DOI", ""),
                             ("JTA", ""),
                             ("ISSN", "0095-3628"),
                             ("PubId", "MECOJKF96-29"),
                             ("PubStatus", 4),
                             ("Status", 4),
                             ("HasAbstract", 0),
                             ("ArticleIds", MultiDict.OrderedMultiDict( (
                                      ("PII", "MECOJKF96-29"),
                                      ("PubMedId",  "9178602"))))) ))
    ]
    
    

# In this case there were no matches
summary3 = """\
<?xml version="1.0"?>
<!DOCTYPE eSummaryResult PUBLIC "-//NLM//DTD eSummaryResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eSummary_020511.dtd">
<eSummaryResult>
</eSummaryResult>
"""
summary3_expected = []

# In this case the database doesn't exist
summaryerr1 = """\
<?xml version="1.0"?>
<!DOCTYPE eSummaryResult PUBLIC "-//NLM//DTD eSummaryResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eSummary_020511.dtd">
<eSummaryResult>
        <ERROR>Wrong DB name</ERROR>
</eSummaryResult>
"""
summaryerr1_expected = Datatypes.EUtilsError("Wrong DB name")


class ParseSummary(unittest.TestCase):
    def _test(self, summary, summary_expected):
        got = parse.parse_summary_xml(cStringIO.StringIO(summary))
        expected = summary_expected
        self.assertEquals(len(got), len(expected))
        for got_summary, expected_summary in zip(got, expected):
            self.assertEquals(got_summary.id, expected_summary.id)

            got_list = got_summary.dataitems.items()
            got_list.sort()
            expected_list = expected_summary.dataitems.items()
            expected_list.sort()

            for x, y in zip(got_list, expected_list):
                if isinstance(x[1], Datatypes.Date):
                    self.assertEquals(x[0], y[0])
                    self.assertEquals(x[1].year, y[1].year)
                    self.assertEquals(x[1].month, y[1].month)
                    self.assertEquals(x[1].day, y[1].day)
                else:
                    self.assertEquals(x, y)
                    
    def testParseSummary1(self):
        self._test(summary1, summary1_expected)

    def testParseSummary2(self):
        self._test(summary2, summary2_expected)

    def testParseSummary3(self):
        self._test(summary3, summary3_expected)

    def _errtest(self, summary, err_expected):
        try:
            parse.parse_summary_xml(cStringIO.StringIO(summary))
        except Datatypes.EUtilsError, err:
            self.assertEquals(str(err), str(err_expected))
        else:
            raise AssertionError("Should not be able to parse that")

    def testParseSummaryErr1(self):
        self._errtest(summaryerr1, summaryerr1_expected)

##  ====================== Test parsing Fetch results

fetch1 = u"""\
<?xml version="1.0"?>
<!DOCTYPE PubmedArticleSet PUBLIC "-//NLM//DTD PubMedArticle, 14th January 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/pubmed_020114.dtd">
<PubmedArticleSet>

<PubmedArticle>
<MedlineCitation Owner="NLM" Status="Completed">
<MedlineID>96365665</MedlineID>
<PMID>8769842</PMID>
<DateCreated>
<Year>1996</Year>
<Month>12</Month>
<Day>20</Day>
</DateCreated>
<DateCompleted>
<Year>1996</Year>
<Month>12</Month>
<Day>20</Day>
</DateCompleted>
<DateRevised>
<Year>2000</Year>
<Month>12</Month>
<Day>18</Day>
</DateRevised>
<Article>
<Journal>
<ISSN>0002-9513</ISSN>
<JournalIssue>
<Volume>270</Volume>
<Issue>1 Pt 2</Issue>
<PubDate>
<Year>1996</Year>
<Month>Jan</Month>
</PubDate>
</JournalIssue>
</Journal>
<ArticleTitle>Adrenal steroids stimulate thiazide-sensitive NaCl transport by rat renal distal tubules.</ArticleTitle>
<Pagination>
<MedlinePgn>F211-9</MedlinePgn>
</Pagination>
<Abstract>
<AbstractText>The current experiments were designed to test the hypothesis that adrenal steroids increase thiazide-sensitive Na and Cl transport by the mammalian renal distal convoluted tubule (DCT). Male Sprague-Dawley rats were adrenalectomized and received steroid hormones by osmotic pumps. Six groups of animals were studied as follows: group I, no hormones; group II, replacement levels of dexamethasone only; group III, replacement levels of aldosterone only; group IV, replacement levels of both hormones; group V; replacement levels of aldosterone and high levels of dexamethasone; and group VI, replacement levels of dexamethasone and high levels of aldosterone. Circulating levels of both hormones were found to be in the high physiological range when infused at the high rate. In vivo microperfusion of distal tubules was performed to determine rates of Na and Cl transport. Chlorothiazide was used to assess the magnitude of electroneutral Na-Cl cotransport. Both aldosterone and dexamethasone stimulated thiazide-sensitive Na and Cl transport by the distal tubule by more than fivefold. [3H]metolazone binding was measured to assess the number of thiazide-sensitive Na-Cl cotransporters in renal cortex. Each steroid also increased the number of [3H]metolazone binding sites in kidney cortex more than threefold. The results are consistent with the presence of both mineralocorticoid and glucocorticoid receptors in the mammalian DCT. Physiological changes in circulating levels of adrenal steroids may affect renal NaCl excretion in part by regulating the rate of electroneutral Na-Cl absorption by the DCT.</AbstractText>
</Abstract>
<Affiliation>Yale University School of Medicine, New Haven 06510, USA.</Affiliation>
<AuthorList CompleteYN="Y">
<Author>
<LastName>Vel\N{LATIN SMALL LETTER A WITH ACUTE}zquez</LastName>
<ForeName>H</ForeName>
<Initials>H</Initials>
</Author>
<Author>
<LastName>Bartiss</LastName>
<ForeName>A</ForeName>
<Initials>A</Initials>
</Author>
<Author>
<LastName>Bernstein</LastName>
<ForeName>P</ForeName>
<Initials>P</Initials>
</Author>
<Author>
<LastName>Ellison</LastName>
<ForeName>D H</ForeName>
<Initials>DH</Initials>
</Author>
</AuthorList>
<Language>eng</Language>
<PublicationTypeList>
<PublicationType>Journal Article</PublicationType>
</PublicationTypeList>
</Article>
<MedlineJournalInfo>
<Country>UNITED STATES</Country>
<MedlineTA>Am J Physiol</MedlineTA>
<NlmUniqueID>0370511</NlmUniqueID>
</MedlineJournalInfo>
<ChemicalList>
<Chemical>
<RegistryNumber>0</RegistryNumber>
<NameOfSubstance>Carrier Proteins</NameOfSubstance>
</Chemical>
<Chemical>
<RegistryNumber>0</RegistryNumber>
<NameOfSubstance>Drug Combinations</NameOfSubstance>
</Chemical>
<Chemical>
<RegistryNumber>0</RegistryNumber>
<NameOfSubstance>sodium-chloride cotransporter</NameOfSubstance>
</Chemical>
<Chemical>
<RegistryNumber>17560-51-9</RegistryNumber>
<NameOfSubstance>Metolazone</NameOfSubstance>
</Chemical>
<Chemical>
<RegistryNumber>50-02-2</RegistryNumber>
<NameOfSubstance>Dexamethasone</NameOfSubstance>
</Chemical>
<Chemical>
<RegistryNumber>52-39-1</RegistryNumber>
<NameOfSubstance>Aldosterone</NameOfSubstance>
</Chemical>
<Chemical>
<RegistryNumber>58-94-6</RegistryNumber>
<NameOfSubstance>Chlorothiazide</NameOfSubstance>
</Chemical>
</ChemicalList>
<CitationSubset>IM</CitationSubset>
<MeshHeadingList>
<MeshHeading>
<DescriptorName MajorTopicYN="N">Adrenalectomy</DescriptorName>
</MeshHeading>
<MeshHeading>
<DescriptorName MajorTopicYN="N">Aldosterone</DescriptorName>
<QualifierName MajorTopicYN="Y">pharmacology</QualifierName>
</MeshHeading>
<MeshHeading>
<DescriptorName MajorTopicYN="N">Animal</DescriptorName>
</MeshHeading>
<MeshHeading>
<DescriptorName MajorTopicYN="N">Binding Sites</DescriptorName>
</MeshHeading>
<MeshHeading>
<DescriptorName MajorTopicYN="N">Carrier Proteins</DescriptorName>
<QualifierName MajorTopicYN="Y">metabolism</QualifierName>
</MeshHeading>
<MeshHeading>
<DescriptorName MajorTopicYN="N">Chlorothiazide</DescriptorName>
<QualifierName MajorTopicYN="Y">pharmacology</QualifierName>
</MeshHeading>
<MeshHeading>
<DescriptorName MajorTopicYN="N">Dexamethasone</DescriptorName>
<QualifierName MajorTopicYN="Y">pharmacology</QualifierName>
</MeshHeading>
<MeshHeading>
<DescriptorName MajorTopicYN="N">Drug Combinations</DescriptorName>
</MeshHeading>
<MeshHeading>
<DescriptorName MajorTopicYN="N">Kidney Cortex</DescriptorName>
<QualifierName MajorTopicYN="N">metabolism</QualifierName>
</MeshHeading>
<MeshHeading>
<DescriptorName MajorTopicYN="N">Kidney Tubules, Distal</DescriptorName>
<QualifierName MajorTopicYN="Y">metabolism</QualifierName>
</MeshHeading>
<MeshHeading>
<DescriptorName MajorTopicYN="N">Male</DescriptorName>
</MeshHeading>
<MeshHeading>
<DescriptorName MajorTopicYN="N">Metolazone</DescriptorName>
<QualifierName MajorTopicYN="N">metabolism</QualifierName>
</MeshHeading>
<MeshHeading>
<DescriptorName MajorTopicYN="N">Rats</DescriptorName>
</MeshHeading>
<MeshHeading>
<DescriptorName MajorTopicYN="N">Rats, Sprague-Dawley</DescriptorName>
</MeshHeading>
<MeshHeading>
<DescriptorName MajorTopicYN="N">Support, Non-U.S. Gov't</DescriptorName>
</MeshHeading>
<MeshHeading>
<DescriptorName MajorTopicYN="N">Support, U.S. Gov't, Non-P.H.S.</DescriptorName>
</MeshHeading>
</MeshHeadingList>
</MedlineCitation>
<PubmedData>
        <History>
                <PubMedPubDate PubStatus="pubmed">
                        <Year>1996</Year>
                        <Month>1</Month>
                        <Day>1</Day>
                </PubMedPubDate>
                <PubMedPubDate PubStatus="medline">
                        <Year>1996</Year>
                        <Month>1</Month>
                        <Day>1</Day>
                        <Hour>0</Hour>
                        <Minute>1</Minute>
                </PubMedPubDate>
        </History>
        <PublicationStatus>ppublish</PublicationStatus>
        <ArticleIdList>
                <ArticleId IdType="pubmed">0008769842</ArticleId>
                <ArticleId IdType="medline">96365665</ArticleId>
        </ArticleIdList>
</PubmedData>
</PubmedArticle>


</PubmedArticleSet>
"""


# empty list or the first identifier is not in the database
# Eg, 
#  http://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?\
#                                 retmode=xml&db=pubmed&id=X

fetcherror1 = """\
<?xml version="1.0"?>
<!DOCTYPE PubmedArticleSet PUBLIC "-//NLM//DTD PubMedArticle, 14th January 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/pubmed_020114.dtd">
<pmFetchResult>
        <ERROR>Empty id list - nothing todo</ERROR>
</pmFetchResult>"""
fetcherror1_expected = "Empty id list - nothing todo"


class ParseFetch(unittest.TestCase):
    def testParseFetchError(self):
        try:
            parse.parse_fetch_publication_xml(cStringIO.StringIO(fetcherror1))
        except Datatypes.EUtilsError, err:
            self.assertEquals(str(err), str(fetcherror1_expected))
        else:
            raise AssertionError("Should not be able to parse that")


##  ====================== Test parsing ELink results

def _build_od(data):
    od = MultiDict.OrderedMultiDict()
    for x in data:
        od[x.linkname] = x
    return od

# PubMed neighbors to 3453362
link1 = """\
<?xml version="1.0"?>
<!DOCTYPE eLinkResult PUBLIC "-//NLM//DTD eLinkResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eLink_020511.dtd">
<eLinkResult>
<LinkSet>
        <DbFrom>pubmed</DbFrom>
        <IdList>
                <Id>3453362</Id>
        </IdList>
        <LinkSetDb>
                <DbTo>pubmed</DbTo>
                <LinkName>pubmed_pubmed</LinkName>
                <Link>
                        <Id>3453362</Id>
                        <Score>2147483647</Score>
                </Link>
                <Link>
                        <Id>3595635</Id>
                        <Score>64783724</Score>
                </Link>
                <Link>
                        <Id>3568116</Id>
                        <Score>60580017</Score>
                </Link>
                <Link>
                        <Id>749683</Id>
                        <Score>15740164</Score>
                </Link>
        </LinkSetDb>
</LinkSet>
</eLinkResult>
"""
link1_expected = Datatypes.NeighborLinkSet(
    Datatypes.DBIds("pubmed", ["3453362"]),
    _build_od([Datatypes.LinkSetDb(dbto = "pubmed",
                                   linkname = "pubmed_pubmed",
                                   links = [Datatypes.Link("3453362", 2147483647),
                                            Datatypes.Link("3595635", 64783724),
                                            Datatypes.Link("3568116", 60580017),
                                            Datatypes.Link("749683", 15740164)]),
     ]))

# PubMed neighbors to 3453362 and 9457642
link2 = """\
<?xml version="1.0"?>
<!DOCTYPE eLinkResult PUBLIC "-//NLM//DTD eLinkResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eLink_020511.dtd">
<eLinkResult>
<LinkSet>
        <DbFrom>pubmed</DbFrom>
        <IdList>
                <Id>3453362</Id>
                <Id>9457642</Id>
        </IdList>
        <LinkSetDb>
                <DbTo>pubmed</DbTo>
                <LinkName>pubmed_pubmed</LinkName>
                <Link>
                        <Id>9457642</Id>
                        <Score>2147483647</Score>
                </Link>
                <Link>
                        <Id>3453362</Id>
                        <Score>2147483647</Score>
                </Link>
                <Link>
                        <Id>3595635</Id>
                        <Score>64783724</Score>
                </Link>
                <Link>
                        <Id>3568116</Id>
                        <Score>60580017</Score>
                </Link>
        </LinkSetDb>
</LinkSet>
</eLinkResult>
"""
link2_expected = Datatypes.NeighborLinkSet(
    Datatypes.DBIds("pubmed", ["3453362", "9457642"]),
    _build_od([Datatypes.LinkSetDb("pubmed", "pubmed_pubmed",
                                   [Datatypes.Link("9457642", 2147483647),
                                    Datatypes.Link("3453362", 2147483647),
                                    Datatypes.Link("3595635", 64783724),
                                    Datatypes.Link("3568116", 60580017)]),
     ]))

# neighbor
# In this case, the 'to' database (db=) does not exist
link3 = """\
<?xml version="1.0"?>
<!DOCTYPE eLinkResult PUBLIC "-//NLM//DTD eLinkResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eLink_020511.dtd">
<eLinkResult>
<LinkSet>
        <DbFrom>pubmed</DbFrom>
        <IdList>
                <Id>3453362</Id>
                <Id>9457642</Id>
        </IdList>
</LinkSet>
</eLinkResult>
"""
link3_expected = Datatypes.NeighborLinkSet(
    Datatypes.DBIds("pubmed", ["3453362", "9457642"]),
    MultiDict.OrderedMultiDict())

# Or should I return an error instead?
##link3_expected = Datatypes.EUtilsError("illegal parameter?")

# neighbor
# In this case, the dbfrom database does not exist
link4 = """\
<?xml version="1.0"?>
<!DOCTYPE eLinkResult PUBLIC "-//NLM//DTD eLinkResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eLink_020511.dtd">
<eLinkResult>
<LinkSet>
        <DbFrom>qweqwe</DbFrom>
        <ERROR>Wrong DB name</ERROR>
</LinkSet>
</eLinkResult>
"""
link4_expected = Datatypes.EUtilsError("Wrong DB name")

# Here the identifier doesn't exist
# Though I actually wanted '34533622342342349457642' -- notice that it
#    got converted to the value of maxint.
link5 = """\
<?xml version="1.0"?>
<!DOCTYPE eLinkResult PUBLIC "-//NLM//DTD eLinkResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eLink_020511.dtd">
<eLinkResult>
<LinkSet>
        <DbFrom>PubMed</DbFrom>
        <IdList>
                <Id>2147483647</Id>
        </IdList>
        <LinkSetDb>
                <DbTo>protein</DbTo>
                <LinkName>pubmed_protein</LinkName>
                <Info>Empty result</Info>
        </LinkSetDb>
</LinkSet>
</eLinkResult>
"""
link5_expected = Datatypes.NeighborLinkSet(
    Datatypes.DBIds("pubmed", ["2147483647"]),
    _build_od([Datatypes.LinkSetDb(dbto = "protein",
                                   linkname = "pubmed_protein",
                                   info = "Empty result")]))


link6 = """\
<?xml version="1.0"?>
<!DOCTYPE eLinkResult PUBLIC "-//NLM//DTD eLinkResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eLink_020511.dtd">
<eLinkResult>
<LinkSet>
        <DbFrom>protein</DbFrom>
        <IdList>
                <Id>123123</Id>
        </IdList>
        <LinkSetDb>
                <DbTo>protein</DbTo>
                <LinkName>protein_protein</LinkName>
                <Link>
                        <Id>123123</Id>
                        <Score>2147483647</Score>
                </Link>
                <Link>
                        <Id>68751</Id>
                        <Score>1786</Score>
                </Link>

                <Link>
                        <Id>431353</Id>
                        <Score>96</Score>
                </Link>
        </LinkSetDb>
        <LinkSetDb>
                <DbTo>protein</DbTo>
                <LinkName>protein_protein_cdart_summary</LinkName>
                <Link>
                        <Id>123123</Id>
                        <Score>2147483647</Score>
                </Link>
                <Link>
                        <Id>123123</Id>
                        <Score>8719</Score>
                </Link>
        </LinkSetDb>
</LinkSet>
</eLinkResult>
"""
link6_expected = Datatypes.NeighborLinkSet(
    Datatypes.DBIds("protein", ["123123"]),
    _build_od([Datatypes.LinkSetDb("protein", "protein_protein",
                                   [Datatypes.Link("123123", 2147483647),
                                    Datatypes.Link("68751", 1786),
                                    Datatypes.Link("431353", 96)]),
               Datatypes.LinkSetDb("protein", "protein_protein_cdart_summary",
                                   [Datatypes.Link("123123", 2147483647),
                                    Datatypes.Link("123123", 8719)])]))

# llinks
link7 = """\
<?xml version="1.0"?>
<!DOCTYPE eLinkResult PUBLIC "-//NLM//DTD eLinkResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eLink_020511.dtd">
<eLinkResult>
<LinkSet>
        <DbFrom>protein</DbFrom>
        <IdUrlList>
                <IdUrlSet>
                        <Id>123123</Id>
                        <ObjUrl>
                                <Url>http://www.ncbi.nlm.nih.gov/Structure/lexington/lexington.cgi?cmd=seq&amp;FILTER=on&amp;EXPECT=0.01&amp;DATALIB=oasis_sap&amp;fasta=123123</Url>
                                <LinkName>Domain Neighbors</LinkName>
                                <SubjectType>protein identification/characterization</SubjectType>
                                <Provider>
                                        <Name>Domain Architecture Retrieval Tool</Name>
                                        <NameAbbr>DART</NameAbbr>
                                        <Id>3240</Id>
                                        <Url></Url>
                                </Provider>
                        </ObjUrl>
                </IdUrlSet>
        </IdUrlList>
</LinkSet>
</eLinkResult>
"""

link7_expected = Datatypes.LinksLinkSet(
    dbfrom = "protein",
    idurlset = [Datatypes.IdUrlSet("123123",
                                 objurls = [
    Datatypes.ObjUrl(
        url = "http://www.ncbi.nlm.nih.gov/Structure/lexington/lexington.cgi?cmd=seq&FILTER=on&EXPECT=0.01&DATALIB=oasis_sap&fasta=123123",
        linkname = "Domain Neighbors",
        subject_types = ["protein identification/characterization"],
        provider = Datatypes.Provider(
                    name = "Domain Architecture Retrieval Tool",
                    name_abbr = "DART",
                    id = "3240")
        ),
    ])])

# elink_using_dbids(DBIds("protein", ["123123"]), cmd = "neighbor",
#                         db = "taxonomy")
link8 = """\
<?xml version="1.0"?>
<!DOCTYPE eLinkResult PUBLIC "-//NLM//DTD eLinkResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eLink_020511.dtd">
<eLinkResult>
<LinkSet>
        <DbFrom>protein</DbFrom>
        <IdList>
                <Id>123123</Id>
        </IdList>
        <LinkSetDb>
                <DbTo>taxonomy</DbTo>
                <LinkName>protein_taxonomy</LinkName>
                <Link>
                        <Id>10243</Id>
                </Link>
        </LinkSetDb>
</LinkSet>
</eLinkResult>
"""

link8_expected = Datatypes.NeighborLinkSet(
    Datatypes.DBIds("protein", ["123123"]),
    _build_od([Datatypes.LinkSetDb("taxonomy", "protein_taxonomy",
                                   [Datatypes.Link("10243")])]))
    

# elink_using_dbids(DBIds("protein", ["123123"]), cmd = "lcheck", db = "taxonomy")
# lcheck doesn't seem to be affected by dbto)
# Looks like lcheck and ncheck have the same format
link9 = """\
<?xml version="1.0"?>
<!DOCTYPE eLinkResult PUBLIC "-//NLM//DTD eLinkResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eLink_020511.dtd">
<eLinkResult>
<LinkSet>
        <DbFrom>protein</DbFrom>
        <IdCheckList>
                <Id HasLinkOut="Y">123123</Id>
        </IdCheckList>
</LinkSet>
</eLinkResult>
"""
link9_expected = Datatypes.CheckLinkSet(
    "protein", [Datatypes.IdCheck("123123", has_linkout = 1)])

# How did I generate this?
link10 = """\
<?xml version="1.0"?>
<!DOCTYPE eLinkResult PUBLIC "-//NLM//DTD eLinkResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eLink_020511.dtd">
<eLinkResult>
</eLinkResult>
"""
link10_expected = Datatypes.EUtilsError("Server failed to process request")

# mindate is after maxdate
link11 = """\
<?xml version="1.0"?>
<!DOCTYPE eLinkResult PUBLIC "-//NLM//DTD eLinkResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eLink_020511.dtd">
<eLinkResult>
<LinkSet>
        <DbFrom>pubmed</DbFrom>
        <ERROR>Can't apply term/daterange</ERROR>
</LinkSet>
</eLinkResult>
"""
link11_expected = Datatypes.EUtilsError("Can't apply term/daterange")

#  To retrieve IDs from nucleotide for GI 18250303, 18250301, 18250299 to protein:
#   dbfrom=nucleotide&db=protein&id=18250303,18250307
# Note: there is no Score
link12 = """\
<?xml version="1.0"?>
<!DOCTYPE eLinkResult PUBLIC "-//NLM//DTD eLinkResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eLink_020511.dtd">
<eLinkResult>
<LinkSet>
	<DbFrom>nucleotide</DbFrom>
	<IdList>
		<Id>18250303</Id>
		<Id>18250307</Id>
	</IdList>
	<LinkSetDb>
		<DbTo>protein</DbTo>
		<LinkName>nucleotide_protein</LinkName>
		<Link>
			<Id>18250308</Id>
		</Link>
		<Link>
			<Id>18250304</Id>
		</Link>
	</LinkSetDb>
</LinkSet>
</eLinkResult>
"""
link12_expected = Datatypes.NeighborLinkSet(
    Datatypes.DBIds("nucleotide", ["18250303", "18250307"]),
    _build_od([Datatypes.LinkSetDb("protein", "nucleotide_protein",
                                   [Datatypes.Link("18250308"),
                                    Datatypes.Link("18250304")])]))


# To list all available links in pubmed for PMIDs 12085856 and 12085853:
# dbfrom=pubmed&id=12085856,12085853&cmd=llinks
link13 = """\
<?xml version="1.0"?>
<!DOCTYPE eLinkResult PUBLIC "-//NLM//DTD eLinkResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eLink_020511.dtd">
<eLinkResult>
<LinkSet>
  <DbFrom>pubmed</DbFrom>
  <IdUrlList>
    <IdUrlSet>
      <Id>12085856</Id>
      <ObjUrl>
        <Url>http://symptomresearch.nih.gov/chapter_1/index.htm</Url>
        <SubjectType>online tutorials/courses</SubjectType>
        <Provider>
          <Name>New England Research Institutes Inc.</Name>
          <NameAbbr>NERI</NameAbbr>
          <Id>3291</Id>
          <Url>http://www.symptomresearch.com</Url>
        </Provider>
      </ObjUrl>
      <ObjUrl>
        <Url>http://www.nlm.nih.gov/medlineplus/coronarydisease.html</Url>
        <LinkName>Coronary Disease</LinkName>
        <SubjectType>consumer health</SubjectType>
        <Provider>
          <Name>MEDLINEplus Health Information</Name>
          <NameAbbr>MEDPLUS</NameAbbr>
          <Id>3162</Id>
          <Url>http://medlineplus.gov/</Url>
          <IconUrl>http://www.nlm.nih.gov/medlineplus/images/linkout_sm.gif</IconUrl>
        </Provider>
      </ObjUrl>
      <ObjUrl>
        <Url>http://www.nlm.nih.gov/medlineplus/heartbypasssurgeryangioplasty.html</Url>
        <LinkName>Heart Bypass Surgery/Angioplasty</LinkName>
        <SubjectType>consumer health</SubjectType>
        <Provider>
          <Name>MEDLINEplus Health Information</Name>
          <NameAbbr>MEDPLUS</NameAbbr>
          <Id>3162</Id>
          <Url>http://medlineplus.gov/</Url>
          <IconUrl>http://www.nlm.nih.gov/medlineplus/images/linkout_sm.gif</IconUrl>
        </Provider>
      </ObjUrl>
    </IdUrlSet>
    <IdUrlSet>
      <Id>12085853</Id>
      <ObjUrl>
        <Url>http://www.nlm.nih.gov/medlineplus/arrhythmia.html</Url>
        <LinkName>Arrhythmia</LinkName>
        <SubjectType>consumer health</SubjectType>
        <Provider>
          <Name>MEDLINEplus Health Information</Name>
          <NameAbbr>MEDPLUS</NameAbbr>
          <Id>3162</Id>
          <Url>http://medlineplus.gov/</Url>
          <IconUrl>http://www.nlm.nih.gov/medlineplus/images/linkout_sm.gif</IconUrl>
        </Provider>
      </ObjUrl>
      <ObjUrl>
        <Url>http://www.nlm.nih.gov/medlineplus/exercisephysicalfitness.html</Url>
        <LinkName>Exercise/Physical Fitness</LinkName>
        <SubjectType>consumer health</SubjectType>
        <Provider>
          <Name>MEDLINEplus Health Information</Name>
          <NameAbbr>MEDPLUS</NameAbbr>
          <Id>3162</Id>
          <Url>http://medlineplus.gov/</Url>
          <IconUrl>http://www.nlm.nih.gov/medlineplus/images/linkout_sm.gif</IconUrl>
        </Provider>
      </ObjUrl>
    </IdUrlSet>
  </IdUrlList>
</LinkSet>
</eLinkResult>
"""
link13_expected = Datatypes.LinksLinkSet(
    "pubmed", [Datatypes.IdUrlSet("12085856",
                                  [Datatypes.ObjUrl(
    subject_types = ["online tutorials/courses"],
    provider = Datatypes.Provider("New England Research Institutes Inc.", "NERI",
                                  "3291", "http://www.symptomresearch.com"),
    url = "http://symptomresearch.nih.gov/chapter_1/index.htm"),

                                   Datatypes.ObjUrl(
    linkname = "Coronary Disease",
    subject_types = ["consumer health"],
    url = "http://www.nlm.nih.gov/medlineplus/coronarydisease.html",
    provider = Datatypes.Provider(name = "MEDLINEplus Health Information",
                                  name_abbr = "MEDPLUS",
                                  id = "3162",
                                  url = "http://medlineplus.gov/",
                                  icon_url = "http://www.nlm.nih.gov/medlineplus/images/linkout_sm.gif")),

                                   Datatypes.ObjUrl(
    ["consumer health"],
    Datatypes.Provider("MEDLINEplus Health Information", "MEDPLUS", "3162",
                       "http://medlineplus.gov/",
                       "http://www.nlm.nih.gov/medlineplus/images/linkout_sm.gif"),
    "Heart Bypass Surgery/Angioplasty",
    "http://www.nlm.nih.gov/medlineplus/heartbypasssurgeryangioplasty.html"),
                                   ]),

               Datatypes.IdUrlSet(id = "12085853",
                                  objurls = [
    Datatypes.ObjUrl(
         url = "http://www.nlm.nih.gov/medlineplus/arrhythmia.html",
         linkname = "Arrhythmia",
         subject_types = ["consumer health"],
         provider = Datatypes.Provider(
             name = "MEDLINEplus Health Information",
             name_abbr = "MEDPLUS",
             id = "3162",
             url = "http://medlineplus.gov/",
             icon_url = "http://www.nlm.nih.gov/medlineplus/images/linkout_sm.gif")),
    
    Datatypes.ObjUrl(
         url = "http://www.nlm.nih.gov/medlineplus/exercisephysicalfitness.html",
         linkname = "Exercise/Physical Fitness",
         subject_types = ["consumer health"],
         provider = Datatypes.Provider(
             name = "MEDLINEplus Health Information",
             name_abbr = "MEDPLUS",
             id = "3162",
             url = "http://medlineplus.gov/",
             icon_url = "http://www.nlm.nih.gov/medlineplus/images/linkout_sm.gif"))
    ])])
    
# elink.fcgi?dbfrom=pubmed&id=10611131&cmd=prlinks
link14 = """\
<?xml version="1.0"?>
<!DOCTYPE eLinkResult PUBLIC "-//NLM//DTD eLinkResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eLink_020511.dtd">
<eLinkResult>
<LinkSet>
	<DbFrom>pubmed</DbFrom>
	<IdUrlList>
		<IdUrlSet>
			<Id>10611131</Id>
			<ObjUrl>
				<Url>http://brain.oupjournals.org/cgi/pmidlookup?view=full&amp;pmid=10611131</Url>

				<SubjectType>publishers/providers</SubjectType>
				<Attribute>publisher of information in URL</Attribute>
				<Attribute>full-text online</Attribute>
				<Provider>
					<Name>HighWire Press</Name>
					<NameAbbr>HighWire</NameAbbr>

					<Id>3051</Id>
					<Url>http://highwire.stanford.edu</Url>
					<IconUrl>http://highwire.stanford.edu/icons/externalservices/pubmed/highwirepress.jpg</IconUrl>
				</Provider>
			</ObjUrl>
		</IdUrlSet>
	</IdUrlList>

</LinkSet>
</eLinkResult>
"""
link14_expected = Datatypes.LinksLinkSet(
    "pubmed", [Datatypes.IdUrlSet("10611131",
                                  [Datatypes.ObjUrl(
    
    url = "http://brain.oupjournals.org/cgi/pmidlookup?view=full&pmid=10611131",
    subject_types = ["publishers/providers"],
    attributes = ["publisher of information in URL", "full-text online"],
    provider = Datatypes.Provider(name = "HighWire Press",
                                  name_abbr = "HighWire",
                                  id = "3051",
                                  url = "http://highwire.stanford.edu",
                                  icon_url = "http://highwire.stanford.edu/icons/externalservices/pubmed/highwirepress.jpg"))])])

# prlinks for protein 4579714, has no primary links
link15 = """\
<?xml version="1.0"?>
<!DOCTYPE eLinkResult PUBLIC "-//NLM//DTD eLinkResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eLink_020511.dtd">
<eLinkResult>
<LinkSet>
        <DbFrom>protein</DbFrom>
        <IdUrlList>
                <IdUrlSet>
                        <Id>4579714</Id>
                        <Info>No primary links</Info>
                </IdUrlSet>
        </IdUrlList>
</LinkSet>
</eLinkResult>
"""
link15_expected = Datatypes.LinksLinkSet(
    "protein", [Datatypes.IdUrlSet("4579714", [])])


# This seems to be an intermittant error
# Seen when doing a neighbor link
link16 = """\
<?xml version="1.0"?>
<!DOCTYPE eLinkResult PUBLIC "-//NLM//DTD eLinkResult, 11 May 2002//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/eLink_020511.dtd">
<eLinkResult>
<LinkSet>
        <DbFrom>nucleotide</DbFrom>
        <IdList>
                <Id>18250303</Id>
        </IdList>
        <LinkSetDb>
                <DbTo>protein</DbTo>
                <LinkName>nucleotide_protein</LinkName>
                <ERROR>CPmCSResult::Fetch() : Can not fetch result</ERROR>
        </LinkSetDb>
</LinkSet>
</eLinkResult>
"""
link16_expected = Datatypes.EUtilsError(
    "CPmCSResult::Fetch() : Can not fetch result")

class ParseLink(unittest.TestCase):
    def _test(self, f, s, expected):
        if isinstance(expected, Datatypes.EUtilsError):
            try:
                f(cStringIO.StringIO(s))
            except Datatypes.EUtilsError, err:
                self.assertEquals(str(err), str(expected))
            else:
                raise AssertionError("Unexpected success")
        else:
            x = f(cStringIO.StringIO(s))
            if (isinstance(x, Datatypes.NeighborLinkSet) or
                isinstance(x, Datatypes.LinksLinkSet) or
                isinstance(x, Datatypes.CheckLinkSet)):
                self.assertEquals(x, expected)
            else:
                raise AssertionError("Unexpected return class")
    
    def testParseLink1(self): self._test(parse.parse_neighbor_links,
                                         link1, link1_expected)
    def testParseLink2(self): self._test(parse.parse_neighbor_links,
                                         link2, link2_expected)
    def testParseLink3(self): self._test(parse.parse_neighbor_links,
                                         link3, link3_expected)
    def testParseLink4(self): self._test(parse.parse_neighbor_links,
                                         link4, link4_expected)
    def testParseLink5(self): self._test(parse.parse_neighbor_links,
                                         link5, link5_expected)
    def testParseLink6(self): self._test(parse.parse_neighbor_links,
                                         link6, link6_expected)
    def testParseLink7(self): self._test(parse.parse_llinks,
                                         link7, link7_expected)
    def testParseLink8(self): self._test(parse.parse_neighbor_links,
                                         link8, link8_expected)
    def testParseLink9(self):
        self._test(parse.parse_lcheck, link9, link9_expected)
        self._test(parse.parse_ncheck, link9, link9_expected)
    
    def testParseLink10(self): self._test(parse.parse_neighbor_links,
                                          link10, link10_expected)
    def testParseLink11(self): self._test(parse.parse_neighbor_links,
                                          link11, link11_expected)
    def testParseLink12(self): self._test(parse.parse_neighbor_links,
                                          link12, link12_expected)
    def testParseLink13(self): self._test(parse.parse_llinks,
                                          link13, link13_expected)
    def testParseLink14(self): self._test(parse.parse_prlinks,
                                          link14, link14_expected)
    def testParseLink15(self): self._test(parse.parse_prlinks,
                                          link15, link15_expected)
    def testParseLink16(self): self._test(parse.parse_neighbor_links,
                                          link16, link16_expected)
        

##  ====================== Test generic errors caused by bad parameters


# gave the wrong database (this text is truncated)
error2 = '''\
<!doctype html public "-//w3c//dtd html 4.0 transitional//en">
<html>
<head>
   <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">
   <meta name="GENERATOR" content="Microsoft FrontPage 5.0">
   <title>EFetch Entrez Utility</title>
<script LANGUAGE="JavaScript">

   </script>
<link rel="stylesheet" href="http://www.ncbi.nlm.nih.gov/corehtml/ncbi2.css">
</head>
<body text="#000000" bgcolor="#FFFFFF" link="#CC6600" vlink="#CC6600">
'''

error2_expected = "Parameter not allowed"

class ParseErrors(unittest.TestCase):
    error_cases = [(error2, error2_expected)]
    def _test(self, f):
        for s, expect in self.error_cases:
            try:
                f(cStringIO.StringIO(s))
            except Datatypes.EUtilsError, err:
                self.assertEquals(str(err), expect)
                
    def testSearchError(self):
        self._test(parse.parse_search)

    ## Don't need to test this -- never raises those errors
##    def testPostError(self):
##        self._test(parse.parse_post)
        
    def testSummary(self):
        self._test(parse.parse_summary_xml)
        
    def testFetch(self):
        self._test(parse.parse_fetch_publication_xml)
        
    def testELink(self):
        self._test(parse.parse_link_xml)


class SummaryConversion(unittest.TestCase):
    def testDates(self):
        for s, d in (("2000 Feb 3", Datatypes.Date(2000, 2, 3)),
                     ("1975 Jun", Datatypes.Date(1975, 6, 1)),
                     ("1999", Datatypes.Date(1999, 1, 1)),
                     ("2000/02/17 09:00", Datatypes.Date(2000, 2, 17)),
                     ("1993 May-Jun", Datatypes.Date(1993, 5, 1))):
            self.assertEquals(parse.convert_summary_Date_string(s), d)
                     
if __name__ == "__main__":
    unittest.main()
