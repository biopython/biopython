import Martel
from Martel import RecordReader, Parser

from xml.sax import handler, saxutils
from StringIO import StringIO

# NOTE: Do not write formats like this, eg, with "(.|\n)*".  Those
# depend on the implementation using a RecordReader-like interface.
# Instead, you need to write them so they could be used even as one
# huge regexp.

def test_reader_parser():
    record = Martel.Group("start", Martel.Rep(Martel.Str("abc"))) + \
             Martel.Group("end", Martel.Rep(Martel.Str("xyz")))
    parser = record.make_parser()

    parser = Parser.Parser(parser.tagtable)
    parser.setErrorHandler(handler.ErrorHandler())

    parser.parseString("abc" * 10 + "xyz")

    try:
        parser.parseString("abc" * 10 + "xyzQ")
    except Parser.ParserPositionException:
        pass
    else:
        raise AssertionError, "didn't get a position exception"

    try:
        parser.parseString("abc" * 10 + "x")
    except Parser.ParserPositionException:
        pass
    else:
        raise AssertionError, "didn't get a position exception"

class CountErrors(handler.ErrorHandler):
    def __init__(self):
        self.error_count = 0
        self.fatal_error_count = 0
    def error(self, exception):
        self.error_count = self.error_count + 1
    def fatalError(self, exception):
        self.fatal_error_count = self.fatal_error_count + 1

class CountRecords(handler.ContentHandler):
    def __init__(self, tag):
        self.tag = tag
        self.count = 0
    def startElement(self, tag, attrs):
        if tag == self.tag:
            self.count = self.count + 1
        
def test_record_parser():
    record = Martel.Group("A", Martel.Str("X\n") + Martel.Re("a*\n"))
    p = record.make_parser()

    parser = Parser.RecordParser("blah", {}, p.tagtable, (0, 1, {}),
                                 RecordReader.StartsWith, ("X",))

    err = CountErrors()
    parser.setErrorHandler(err)
    count = CountRecords("A")
    parser.setContentHandler(count)
    
    parser.parseString("X\na\nX\nb\nX\naaa\nX\naaaa\nX\nq\nX\na\n")
    
    assert err.fatal_error_count == 0, err.fatal_error_count
    assert err.error_count == 2, err.error_count
    assert count.count == 4, count.count

def test_header_footer1():
    s = """\
header
XX
record 1
//
record 2
//
record 3
//
footer
"""
    gold = """\
<?xml version="1.0" encoding="iso-8859-1"?>
<hf><header>header
XX
</header><record>record 1
//
</record><record>record 2
//
</record><record>record 3
//
</record><footer>footer
</footer></hf>"""
    
    debug_level = 1
    
    # Don't use regexps like these in your code - for testing only!
    header = Martel.Group("header", Martel.Re(r"header(.|\n)*"))
    record = Martel.Group("record", Martel.Re(r"rec(.|\n)*"))
    footer = Martel.Group("footer", Martel.Re(r"footer(.|\n)*"))

    header = header.make_parser(debug_level = debug_level)
    record = record.make_parser(debug_level = debug_level)
    footer = footer.make_parser(debug_level = debug_level)

    hf = Parser.HeaderFooterParser(
        "hf", {},
        RecordReader.EndsWith, ("XX\n", ), header.tagtable,
        RecordReader.EndsWith, ("//\n", ), record.tagtable,
        RecordReader.StartsWith, ("f", ), footer.tagtable,
        (0, debug_level, {}))

    outfile = StringIO()
    hf.setContentHandler(saxutils.XMLGenerator(outfile))
    hf.setErrorHandler(handler.ErrorHandler())
    hf.parseFile(StringIO(s))

    result = outfile.getvalue()
    assert result == gold, (result, gold)


def test_header_footer2():
    # Have a header but no footer
    s = """
This is some misc. header text
that goes on until the end.
ID 1
This is some data
ID 2
This is some more data
"""
    gold = """\
<?xml version="1.0" encoding="iso-8859-1"?>
<hf><header>
This is some misc. header text
that goes on until the end.
</header><record>ID 1
This is some data
</record><record>ID 2
This is some more data
</record></hf>"""
    
    # Don't use a regexp like this in your code - for testing only!
    header = Martel.Group("header", Martel.Re(r"(.|\n)*"))
    record = Martel.Group("record", Martel.Re(r"ID \d+(.|\n)*"))

    header = header.make_parser()
    record = record.make_parser()

    hf = Parser.HeaderFooterParser(
        "hf", {},
        RecordReader.Until, ("ID", ), header.tagtable,
        RecordReader.StartsWith, ("ID", ), record.tagtable,
        RecordReader.Nothing, (), (),
        (0, 1, {}))

    outfile = StringIO()
    hf.setContentHandler(saxutils.XMLGenerator(outfile))
    hf.setErrorHandler(handler.ErrorHandler())
    hf.parseFile(StringIO(s))

    text = outfile.getvalue()
    assert text == gold, (text, gold)

def test_header_footer3():
    # Have a footer but no header
    s = """\
ID 1
This is some data
//
ID 2
This is some more data
//
Okay, that was all of the data.
"""
    gold = """\
<?xml version="1.0" encoding="iso-8859-1"?>
<hf><record>ID 1
This is some data
//
</record><record>ID 2
This is some more data
//
</record><footer>Okay, that was all of the data.
</footer></hf>"""
    
    # Don't use a regexp like this in your code - for testing only!
    record = Martel.Group("record", Martel.Re(r"ID \d+(.|\n)*"))
    # Require at least 5 characters (just to be safe)
    footer = Martel.Group("footer", Martel.Re(r".....(.|\n)*"))

    record = record.make_parser()
    footer = footer.make_parser()

    hf = Parser.HeaderFooterParser(
        "hf", {},
        RecordReader.Nothing, (), (),
        RecordReader.EndsWith, ("//\n", ), record.tagtable,
        RecordReader.Everything, (), footer.tagtable,
        (0, 1, {}))

    outfile = StringIO()
    hf.setContentHandler(saxutils.XMLGenerator(outfile))
    hf.setErrorHandler(handler.ErrorHandler())
    hf.parseFile(StringIO(s))

    text = outfile.getvalue()
    assert text == gold, (text, gold)
    
def test_header_footer4():
    # Have a header but no footer - and not footer reader
    s = """
This is some misc. header text
that goes on until the end.
ID 1
This is some data
ID 2
This is some more data
"""
    gold = """\
<?xml version="1.0" encoding="iso-8859-1"?>
<hf><header>
This is some misc. header text
that goes on until the end.
</header><record>ID 1
This is some data
</record><record>ID 2
This is some more data
</record></hf>"""

    # Don't use a regexp like this in your code - for testing only!
    header = Martel.Group("header", Martel.Re(r"(.|\n)*"))
    record = Martel.Group("record", Martel.Re(r"ID \d+(.|\n)*"))

    header = header.make_parser()
    record = record.make_parser()

    hf = Parser.HeaderFooterParser(
        "hf", {},
        RecordReader.Until, ("ID", ), header.tagtable,
        RecordReader.StartsWith, ("ID", ), record.tagtable,
        None, (), (),
        (0, 1, {}))
    outfile = StringIO()
    hf.setContentHandler(saxutils.XMLGenerator(outfile))
    hf.setErrorHandler(handler.ErrorHandler())
    hf.parseFile(StringIO(s))

    text = outfile.getvalue()
    assert text == gold, (text, gold)

def test_header_footer5():
    # Make sure I can skip records when there are not footer records
    s = """
This is some misc. header text
that goes on until the end.
ID 1
This is some data
ID A
This is some more data
ID 3
This is again some more data
ID Q
This blah
ID W
QWE
ID 987
To be
ID 897
Or not to be
"""
    header = Martel.Group("header", Martel.Re(r"(.|\n)*"))
    record = Martel.Group("record", Martel.Re(r"ID \d+(.|\n)*"))

    header = header.make_parser()
    record = record.make_parser()

    hf = Parser.HeaderFooterParser(
        "hf", {},
        RecordReader.Until, ("ID", ), header.tagtable,
        RecordReader.StartsWith, ("ID", ), record.tagtable,
        None, (), (),
        (0, 1, {}))
    count = CountRecords("record")
    hf.setContentHandler(count)
    err = CountErrors()
    hf.setErrorHandler(err)
    hf.parseFile(StringIO(s))

    assert err.error_count == 3, err.error_count
    assert err.fatal_error_count == 0, err.fatal_error_count
    assert count.count == 4, count.count

def test_header_footer6():
    # Make sure I can skip records when there are footer records
    s = """
This is some misc. header text
that goes on until the end.
ID 1
This is some data
//
ID A
This is some more data
//
ID 3
This is again some more data
//
ID Q
This blah
//
ID W
QWE
//
ID 987
To be
//
ID 897
Or not to be
//
FOOTER
"""
    header = Martel.Group("header", Martel.Re(r"(.|\n)*"))
    record = Martel.Group("record", Martel.Re(r"ID \d+(.|\n)*"))
    footer = Martel.Group("footer", Martel.Re("FOOTER(.|\n)*"))

    header = header.make_parser()
    record = record.make_parser()
    footer = footer.make_parser()

    hf = Parser.HeaderFooterParser(
        "hf", {},
        RecordReader.Until, ("ID", ), header.tagtable,
        RecordReader.EndsWith, ("//", ), record.tagtable,
        RecordReader.StartsWith, ("FOOTER", ), footer.tagtable,
        (0, 1, {}))
    count = CountRecords("record")
    hf.setContentHandler(count)
    err = CountErrors()
    hf.setErrorHandler(err)
    hf.parseFile(StringIO(s))

    assert err.error_count == 3, err.error_count
    assert err.fatal_error_count == 0, err.fatal_error_count
    assert count.count == 4, count.count

def test_header_footer7():
    # header and footer but with no record data
    s = """\
This is some misc. header text
that goes on until the end.
FOOTER
"""
    header = Martel.Group("header", Martel.Re(r"(.|\n)*"))
    record = Martel.Group("record", Martel.Re(r"ID \d+(.|\n)*"))
    footer = Martel.Group("footer", Martel.Re("FOOTER(.|\n)*"))

    header = header.make_parser()
    record = record.make_parser()
    footer = footer.make_parser()

    hf = Parser.HeaderFooterParser(
        "hf", {},
        RecordReader.CountLines, (2, ), header.tagtable,
        RecordReader.EndsWith, ("//", ), record.tagtable,
        RecordReader.StartsWith, ("FOOTER", ), footer.tagtable,
        (0, 1, {}))
    count = CountRecords("record")
    hf.setContentHandler(count)
    err = CountErrors()
    hf.setErrorHandler(err)
    hf.parseFile(StringIO(s))

    assert err.error_count == 0, err.error_count
    assert err.fatal_error_count == 0, err.fatal_error_count
    assert count.count == 0, count.count

def test_header_footer8():
    # header, record and footer, but with extra data
    s1 = """Two lines in
the header.
Data 1
Data 2
Data Q
Data 4
FOOTER Abc
FOOTER B
"""
    s2 = """Two lines in
the header.
Data 1
Data 2
Data Q
Data 4
FOOTER Abc
"""
    s3 = """Two lines in
the header.
Data 1
Data 4
FOOTER Abc
"""
    s4 = """Two lines in
the header.
Data Q
FOOTER Abc
"""
    s5 = """Two lines in
the header.
FOOTER Abc
"""
    dataset = ( (s1, 3, 1, 1),
                (s2, 3, 1, 0),
                (s3, 2, 0, 0),
                (s4, 0, 1, 0),
                (s5, 0, 0, 0),
                )
                
    header = Martel.Group("header", Martel.Re(r"(.|\n)*"))
    record = Martel.Group("record", Martel.Re(r"Data \d+\n"))
    footer = Martel.Group("footer", Martel.Re("FOOTER \w+\n"))

    header = header.make_parser()
    record = record.make_parser()
    footer = footer.make_parser()

    hf = Parser.HeaderFooterParser(
        "hf", {},
        RecordReader.CountLines, (2, ), header.tagtable,
        RecordReader.CountLines, (1, ), record.tagtable,
        RecordReader.CountLines, (1, ), footer.tagtable,
        (0, 1, {}))
    for s, rec_count, err_count, fatal_count in dataset:
        count = CountRecords("record")
        hf.setContentHandler(count)
        err = CountErrors()
        hf.setErrorHandler(err)
        hf.parseFile(StringIO(s))

        assert err.error_count == err_count, (s, err.error_count, err_count)
        assert err.fatal_error_count == fatal_count, \
               (s, err.fatal_error_count, fatal_count)
        assert count.count == rec_count, (s, count.count, rec_count)

def test():
    test_reader_parser()
    test_record_parser()
    test_header_footer1()
    test_header_footer2()
    test_header_footer3()
    test_header_footer4()
    test_header_footer5()
    test_header_footer7()
    test_header_footer8()

if __name__ == "__main__":
    test()
