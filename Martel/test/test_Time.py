import re, time, string
from xml.sax import handler

import Martel
from Martel import Time

def test_terms(s, matches, non_matches):
    pat = Time.make_pattern(s, tag_format = None)
    match = re.compile(str(pat) + "$").match

    exp = Time.make_expression(s)
    exp_match = exp.make_parser().parseString
    
    for m in matches:
        if match(m) is None:
            raise AssertionError("%re: s %s does not match %s" % \
                                 (repr(s), repr(pat), repr(m)))
        try:
            exp_match(m)
        except Martel.Parser.ParserException:
            raise AssertionError("exp: %s %s does not match %s" % \
                                 (repr(s), repr(pat), repr(m)))
            
        
    for m in non_matches:
        if match(m) is not None:
            raise AssertionError("re: %s %s should not match %s" % \
                                 (repr(s), repr(pat), repr(m)))
        try:
            exp_match(m)
        except Martel.Parser.ParserException:
            pass
        else:
            raise AssertionError("exp: %s %s should not match %s" % \
                                 (repr(s), repr(pat), repr(m)))

def test_syntax():
    for pat in ("%A", "%%", "%A %A %B %C%%", "%", "Test%%", "Test%",
                "%Test", "%Test%", "%%Test%%%", "nothing", "",
                "%(A)", "%(day)"):
        Time.make_pattern(pat)
        Time.make_expression(pat)

    for bad_pat in ("%9A", "%9", "%(", "%()", "%(missing)", "%!", "%Q"):
        try:
            Time.make_pattern(bad_pat)
        except (TypeError, KeyError):
            pass
        else:
            raise AssertionError("Should not have allowed: %s" %
                                 (repr(bad_pat),))

def test_times():
    now = time.time()
    table = []
    for c in Time._time_table.keys():
        if len(c) != 1:
            continue
        pat = Time.make_pattern("%" + c, tag_format = None)
        exp = Time.make_expression("%" + c)
        table.append( (c, pat, re.compile(pat + "$").match,
                       exp.make_parser().parseString) )
    
    for j in range(-5, 5):
        for i in range(25):
            t = time.localtime(now + (i + j*100)*23*3603)
            for c, pat, match, parse in table:
                s = time.strftime("%" + c, t)
                m = match(s)
                if m is None or m.endpos != len(s):
                    raise AssertionError("Bad re %" + \
                                         "%s %s %s" % (c, pat, s))
                try:
                    parse(s)
                except Martel.Parser.ParserException:
                    raise AssertionError("Bad exp %" + \
                                         "%s %s %s" % (c, pat, s))

def test_expand():
    # make sure tag name expansion works correctly
    class Capture:
        def __init__(self):
            self.capture = []
        def __mod__(self, s):
            self.capture.append(s)
            return s
    cap = Capture()
    exp = Time.make_expression("%m-%Y", cap)
    parser = exp.make_parser()
    parser.parseString("05-1921")
    x = cap.capture
    x.sort()
    assert x == ["month", "year"]
    

def _find_quoted_words(line):
    words = []
    i = 0
    while 1:
        i = string.find(line, '"', i)
        if i == -1:
            break
        j = string.find(line, '"', i+1)
        assert i < j, line
        words.append(line[i+1:j])
        i = j + 1
    return words

def test_docstring():
    check_pattern = 0
    for line in Time.__doc__.split("\n"):
        if line[:1] == "%":
            pattern = line.split()[0]
            check_pattern = 1
            examples = []
            element_name = None
            element_attrs = None
            has_continue = "Definition"
            continue
        if not check_pattern:
            continue
        if line.find("Pattern:") != -1:
            has_continue = "Pattern"
            continue
        if line.find("Example:") != -1:
            has_continue = "Example"
            examples.extend(_find_quoted_words(line))
            continue
        if line.find("Element name:") != -1:
            element_name, = _find_quoted_words(line)
            has_continue = 0
            continue

        if line.find("Element attributes:") != -1:
            element_attrs = _find_quoted_words(line)
            has_continue = 0
            if not element_attrs:
                assert line.find("no attributes") != -1, line
            else:
                assert len(element_attrs) == 2, line
                assert element_attrs[0] == "type", (line, element_attrs[0])
        elif line.find("Element:") != -1:
            # Neither element name nor attrs may be specified here
            has_continue = 0
            assert element_name is None, element_name
            assert element_attrs is None, element_attrs
        else:
            if has_continue in ("Pattern", "Definition"):
                continue
            if has_continue == "Example":
                examples.extend(_find_quoted_words(line))
                continue
            raise AssertionError(line)

        assert (element_name is None) == (element_attrs is None), \
               (line, element_name, element_attrs)

        exp = Time.make_expression(pattern)
        parser = exp.make_parser()

        class StoreAttrs(handler.ContentHandler):
            def startElement(self, name, attrs):
                self.name = name
                self.attrs = attrs
        store = StoreAttrs()
        parser.setContentHandler(store)
        parseString = parser.parseString

        for example in examples:
            parseString(example)
            if element_name is not None:
                assert element_name == store.name, (element_name, store.name)
                if element_attrs:
                    if store.attrs[element_attrs[0]] != element_attrs[1]:
                        raise AssertionError(
                            "pattern = %s ; text = %s: %s != %s" % \
                            (repr(pattern), repr(example),
                             repr(element_attrs[0]), repr(element_attrs[1])))
        check_pattern = 0
        pattern = None
            

def test():
    test_docstring()
    test_expand()
    test_syntax()

    test_times()

    test_terms("%(Jan)",
               ("Jan", "JAN", "jan", "Feb", "MAR", "apR", "mAy",
                "jUN", "JUL", "aug", "sEP", "Oct", "NOv", "DeC"),
               ("Jan.", "October", ""))

    test_terms("%(January)",
               ("May", "March", "January", "October"),
               ("Jan.", "JAN", "jan", "", "Mayo"))

    test_terms("%(Mon)day",
               ("Monday", "TUEday", "Wedday", "fRiday", "THuday"),
               ("MON", "Sun", "wednesday", ""))

    test_terms("%(Monday)",
               ("Monday", "TUESDAY", "WedNeSdAY"),
               ("MON", "Mon", "mon", "wedday", ""))

    test_terms("%(second)",
               ("00", "09", "10", "01", "60", "61", "00"),
               ("-1", "0", "9", " 0", "62", "A", "1A", "A1"))

    test_terms("%(minute)",
               ("00", "01", "09", "10", "59"),
               ("60", "1", "A", "1A", "A1"))  # Allow "1" and " 1"?

    test_terms("%(hour)",
               ("00", "01", "02", "03", "04", "05", "06", "07", "08", "09",
                "10", "11", "12", "13", "14", "15", "16", "17", "18", "19",
                "20", "21", "22", "23", " 1", " 2", " 3", " 4", " 5", " 6",
                " 7", " 8", " 9", "1", "2", "3", "4", "5", "6", "7", "8", "9"),
               ("24", "-1", "A", "1A", "A1", "123", "  1"))

    test_terms("%(12-hour)",
               ("01", "02", "03", "04", "05", "06", "07", "08", "09",
                "10", "11", "12", " 1", " 2", " 3", " 4", " 5", " 6", " 7",
                " 8", " 9", "1", "2", "3", "4", "5", "6", "7", "8", "9"),
               ("00", "13", "14", "15", "24", "20", "-1", "A",
                "1A", "A1", "123", "  1"))    

    test_terms("%(24-hour)",
               ("00", "01", "02", "03", "04", "05", "06", "07", "08", "09",
                "10", "11", "12", "13", "14", "15", "16", "17", "18", "19",
                "20", "21", "22", "23", " 1", " 2", " 3", " 4", " 5", " 6",
                " 7", " 8", " 9", "1", "2", "3", "4", "5", "6", "7", "8", "9"),
               ("24", "-1", "A", "1A", "A1", "123", "  1"))

    test_terms("%(day)",
               (      "01", "02", "03", "04", "05", "06", "07", "08", "09",
                "10", "11", "12", "13", "14", "15", "16", "17", "18", "19",
                "20", "21", "22", "23", "24", "25", "26", "27", "28", "29",
                "30", "31", " 1", " 2", " 3", " 4", " 5", " 6", " 7", " 8",
                " 9", "1", "2", "3", "4", "5", "6", "7", "8", "9"),
               ("00", "32", "3 ", "A", "3A"))
               
    test_terms("%(month)",
               (      "01", "02", "03", "04", "05", "06", "07", "08", "09",
                "10", "11", "12", " 1", " 2", " 3", " 4", " 5", " 6", " 7",
                " 8", " 9", "1", "2", "3", "4", "5", "6", "7", "8", "9"),
               ("00", "13", "21", "20", "F", "1F"))

    test_terms("%(YY)",
               ("00", "01", "10", "99"),
               ("100", "00A", "2000", "AA"))

    test_terms("%(YYYY)",
               ("2000", "2001", "1970", "1492"),
               ("01", "ABCD", "1997A"))

    test_terms("%(year)",
               ("00", "01", "10", "99", "2000", "2001", "1970", "1492"),
               ("100", "00A", "AA", "ABCD", "1997A", ""))

    test_terms("%(Jan) %(day), %(YYYY)",
               ("OCT 31, 2000", "Jan 01, 1900", "Jan 1, 1900"),
               ("Aug 99, 2000", "Aug 22, 12345", "June 1 1900"))

    test_terms("%j",
               ("001", "002", "009", "099", "197", "226", "300",
                "301", "309", "310", "320", "330", "340", "350",
                "360", "365", "366"),
               ("000", "0", "00", "1", "12", "367", "1A", "200B"))

    test_terms("%Z",
               ("MST", "GMT", "Pacific Standard Time", "GRNLNDST",
                "MET DST", "New Zealand Standard Time", "NZST",
                "SAST", "GMT+0200", "IDT"),
               ("(MST)", "GMT0200", "+1000"))

    test_terms("%A%%%a",
               ("Monday%Mon", "Tuesday%FRI"),
               ("Thu%THU", "Thursday%%Sunday", "Tuesday%FRI!"))

    test_terms("%",
               ("%",),
               ("%%", ""))

    # make sure terms like "." and "[]" are escaped before turning into
    # a regular expression
    test_terms("A.B][",
               ("A.B][",),
               ("A-B][",))
    
if __name__ == "__main__":
    test()
