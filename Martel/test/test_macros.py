import Martel
from xml.sax import handler, saxutils
import string, StringIO

# Various "macro" definitions so you don't need to come up with your
# own regular expressions.

def must_parse(test_name, parseString, term):
    try:
        parseString(term)
    except Martel.Parser.ParserException:
        raise AssertionError("%s: Cannot parse %s" % \
                             (test_name, repr(term)))

def must_not_parse(test_name, parseString, term):
    try:
        parseString(term)
    except Martel.Parser.ParserException:
        pass
    else:
        raise AssertionError("%s: Could parse %s" % \
                             (test_name, repr(term)))

class Capture(handler.ContentHandler):
    def startDocument(self):
        self.capture = []
    def startElement(self, name, attrs):
        if not self.capture:
            self.capture.append( (name, attrs) )  # dup first item
        self.capture.append( (name, attrs) )
    def endElement(self, name):
        self.capture.pop()

def has_group(format, term, name, attr):
    parser = format.make_parser()
    cap = Capture()
    parser.setContentHandler(cap)
    try:
        parser.parseString(term)
    except Martel.Parser.ParserException:
        raise AssertionError("Cannot parse %s" % (repr(term),))
    assert len(cap.capture) == 1, "Should only have dup of first item"
    assert cap.capture[0][0] == name, (cap.capture[0], name)
    assert cap.capture[0][1]["x"] == attr, (cap.capture[0][1], attr)

def has_no_group(format, term):
    parser = format.make_parser()
    cap = Capture()
    parser.setContentHandler(cap)
    try:
        parser.parseString(term)
    except Martel.Parser.ParserException:
        raise AssertionError("Cannot parse %s" % (repr(term),))
    assert not cap.capture, "Should have nothing, not %s" % (cap.capture,)



def test_Float():
    parseString = Martel.Float().make_parser().parseString
    
    for head in ("", "-", "+", "-1", "+2", "3"):
        for tail in ("", "E0", "E+0", "E-0", "E4", "e+5", "e-6",
                     "E10", "E-19", "e+28"):
            for middle in (".1", "5.", "7.6", "989", ".0001"):
                must_parse("Float", parseString, head + middle + tail)

    for term in ("1E", ".E", "1.E", "1/", "E0", "1.2E0K",
                 "=1", "+-1", ".", "e", "-e", "-e0"):
        must_not_parse("not Float", parseString, term)

    has_group(Martel.Float("spam", {"x": "spot"}), "1.0", "spam", "spot")
    has_group(Martel.Float("eggs", {"x": "SPOT"}), "0.8", "eggs", "SPOT")
    has_no_group(Martel.Float(), "+1")

def test_Digits():
    parseString = Martel.Digits().make_parser().parseString
    for term in ("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                 "20", "99", "453", "34653", "34359739467623"):
        must_parse("Digits", parseString, term)

    for term in ("A", "1A", "123123123T", "-1"):
        must_not_parse("not Digits", parseString, term)

    has_group(Martel.Digits("spam", {"x": "this"}), "5", "spam", "this")
    has_group(Martel.Digits("eggs", {"x": "that"}), "9", "eggs", "that")
    has_no_group(Martel.Digits(), "00")

def test_Word():
    parseString = Martel.Word().make_parser().parseString
    for term in ("Andrew", "Dalke", "was_here", "test12", "12df"):
        must_parse("Word", parseString, term)
    for term in ("*", "", "this-that"):
        must_not_parse("not Word", parseString, term)

    has_group(Martel.Word("spam", {"x": "fly"}), "_", "spam", "fly")
    has_group(Martel.Word("eggs", {"x": "boy"}), "9", "eggs", "boy")
    has_no_group(Martel.Word(), "__")

def test_Spaces():
    parseString = Martel.Spaces().make_parser().parseString
    for term in (" ", "\t", "  ", " \t  \t\t  "):
        must_parse("Spaces", parseString, term)
    for term in ("\n", " \n", " X ", ""):
        must_not_parse("not Spaces", parseString, term)
    has_group(Martel.Spaces("spam", {"x": "pick"}), " "*100,
              "spam", "pick")
    has_group(Martel.Spaces("eggs", {"x": "name"}), "\t"*200,
              "eggs", "name")
    has_no_group(Martel.Spaces(), " ")
    
def test_Unprintable():
    parseString = Martel.Unprintable().make_parser().parseString
    unprintables = []
    for i in range(0, 256):
        c = chr(i)
        if c in string.printable:
            must_not_parse("not Unprintable", parseString, c)
        else:
            must_parse("Unprintable", parseString, c)
            unprintables.append(c)

    has_group(Martel.Unprintable("spam", {"x": "import"}),
              unprintables[0], "spam", "import")
    has_group(Martel.Unprintable("eggs", {"x": "export"}),
              unprintables[-1], "eggs", "export")
    has_no_group(Martel.Unprintable(), unprintables[1])

def test_Punctuation():
    parseString = Martel.Punctuation().make_parser().parseString
    for i in range(0, 256):
        c = chr(i)
        if c in string.punctuation:
            must_parse("Punctuation", parseString, c)
        else:
            must_not_parse("not Punctuation", parseString, c)

    has_group(Martel.Punctuation("spam", {"x": "Iran"}),
              string.punctuation[0], "spam", "Iran")
    has_group(Martel.Punctuation("eggs", {"x": "Iraq"}),
              string.punctuation[-1], "eggs", "Iraq")
    has_no_group(Martel.Punctuation(), string.punctuation[1])

def test_ToEol():
    parser = Martel.ToEol("SantaFe").make_parser()
    parseString = parser.parseString
    must_parse("ToEol", parseString, "Testing, 1, 2, 3!\n")
    must_parse("ToEol", parseString, "Andrew\n")
    must_not_parse("ToEol", parseString, "Dalke")
    must_not_parse("ToEol", parseString, "This\nis")
    must_not_parse("ToEol", parseString, "This\nis a test\n")

    file = StringIO.StringIO()
    parser.setContentHandler(saxutils.XMLGenerator(file))
    parser.parseString("This is a test.\n")
    s = file.getvalue()
    expect = "<SantaFe>This is a test.</SantaFe>\n"
    assert string.find(s, expect) != -1, ("Got: %s" % (repr(s),))

def test_ToSep():
    exp = Martel.Group("test",
                       Martel.ToSep("colon", ":") + \
                       Martel.ToSep("space", " ") + \
                       Martel.ToSep("empty", "!"))
    parser = exp.make_parser()

    file = StringIO.StringIO()
    parser.setContentHandler(saxutils.XMLGenerator(file))
    parser.parseString("q:wxy !")
    s = file.getvalue()
    expect = "<test><colon>q</colon>:<space>wxy</space> <empty></empty>!</test>"
    assert string.find(s, expect) != -1, ("Got: %s" % (repr(s),))
    

def test_DelimitedFields():
    exp = Martel.Group("test", Martel.DelimitedFields("Field", "/"))
    parser = exp.make_parser()

    file = StringIO.StringIO()
    parser.setContentHandler(saxutils.XMLGenerator(file))
    parser.parseString("a/b/cde/f//\n")
    s = file.getvalue()
    expect = "<test><Field>a</Field>/<Field>b</Field>/<Field>cde</Field>/" \
             "<Field>f</Field>/<Field></Field>/<Field></Field>\n</test>"
    assert string.find(s, expect) != -1, ("Got: %s" % (repr(s),))


def test():
    test_Digits()
    test_Float()
    test_Word()
    test_Spaces()
    test_Unprintable()
    test_Punctuation()
    test_ToEol()
    test_ToSep()
    test_DelimitedFields()

if __name__ == "__main__":
    test()
