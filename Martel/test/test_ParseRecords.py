# Test that the ParseRecords object produces the same code as the
# underlying expression.

from xml.sax import handler, saxutils
from cStringIO import StringIO
import Martel
from testformats import swissprot38

import test_swissprot38

text = test_swissprot38.record1 + test_swissprot38.record2

def test():
    s1 = StringIO()
    parser = swissprot38.format_expression.make_parser()
    parser.setErrorHandler(handler.ErrorHandler())
    parser.setContentHandler(saxutils.XMLGenerator(s1))
    parser.parseString(text)

    s2 = StringIO()
    parser = swissprot38.format.make_parser()
    parser.setErrorHandler(handler.ErrorHandler())
    parser.setContentHandler(saxutils.XMLGenerator(s2))
    parser.parseString(text)

    s3 = StringIO()
    parser = swissprot38.format.expression.make_parser()
    parser.setErrorHandler(handler.ErrorHandler())
    parser.setContentHandler(saxutils.XMLGenerator(s3))
    parser.parseString(text)

    assert s1.getvalue() == s2.getvalue() == s3.getvalue()

if __name__ == "__main__":
    test()
    
