import os
from xml.sax import handler

import Martel
from testformats import delimiter

sample_dir = "samples"

expected_records = [
    ["Name", "Language"],
    ["Andrew", "Python"],
    ["Jeff", "Python"],
    ["Ewan", "Perl"],
    ["Jason", "Java"],
]

class CatchFields(handler.ContentHandler):
    def __init__(self):
        handler.ContentHandler.__init__(self)
        self._records = None
        self._fields = None
        self._text = None
        
    def startElement(self, name, attrs):
        if name == "delimited":
            self._records = []
        elif name == "record":
            self._fields = []
        elif name == "field":
            self._text = ""
            
    def characters(self, s):
        if self._text is not None:
            self._text = self._text + s
        
    def endElement(self, name):
        if name == "field":
            self._fields.append(self._text)
            self._text = None
        elif name == "record":
            self._records.append(self._fields)
            self._fields = None

def do_test(format, infile):
    parser = format.make_parser()
    catch = CatchFields()
    parser.setContentHandler(catch)
    parser.parseFile(infile)
    assert expected_records == catch._records
        
def test():
    infile = open(os.path.join(sample_dir, "sample.tab"))
    do_test(delimiter.tabformat, infile)

    infile = open(os.path.join(sample_dir, "sample.comma"))
    do_test(delimiter.commaformat, infile)

if __name__ == "__main__":
    test()
