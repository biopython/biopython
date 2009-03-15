import urllib
from xml.sax import handler, make_parser, expatreader
from xml.sax.expatreader import ExpatParser
from xml.sax._exceptions import SAXParseException

class Record(list):
    """\
This record is a list containing the search results returned by
ScanProsite. The record also contains the data members n_match, n_seq,
capped, and warning."""

    def __init__(self):
        self.n_match = None
        self.n_seq = None
        self.capped = None
        self.warning = None


def scan(seq="", mirror='http://www.expasy.org', output='xml', **keywords):
    """Execute a ScanProsite search.

    mirror:      The ScanProsite mirror to be used
                 (default: http://www.expasy.org).
    seq:         The query sequence, or UniProtKB (Swiss-Prot,
                 TrEMBL) accession
    output:      Format of the search results
                 (default: xml)
    
    Further search parameters can be passed as keywords; see the
    documentation for programmatic access to ScanProsite at
    http://www.expasy.org/tools/scanprosite/ScanPrositeREST.html
    for a description of such parameters.

    This function returns a handle to the search results returned by
    ScanProsite. Search results in the XML format can be parsed into a
    Python object, by using the Bio.ExPASy.ScanProsite.read function.
    """
    parameters = {'seq': seq,
                  'output': output}
    for key, value in keywords.iteritems():
        if value is not None:
            parameters[key] = value
    command = urllib.urlencode(parameters)
    url = "%s/cgi-bin/prosite/PSScan.cgi?%s" % (mirror, command)
    handle = urllib.urlopen(url)
    return handle

def read(handle):
    "Parse search results returned by ScanProsite into a Python object"
    content_handler = ContentHandler()
    saxparser = Parser()
    saxparser.setContentHandler(content_handler)
    saxparser.parse(handle)
    record = content_handler.record
    return record

# The functions below are considered private

class Parser(ExpatParser):

    def __init__(self):
        ExpatParser.__init__(self)
        self.firsttime = True

    def feed(self, data, isFinal = 0):
        # Error messages returned by the ScanProsite server are formatted as
        # as plain text instead of an XML document. To catch such error
        # messages, we override the feed method of the Expat parser.
        # The error message is (hopefully) contained in the data that was just
        # fed to the parser.
        if self.firsttime:
            if data[:5]!="<?xml":
                raise ValueError, data
        self.firsttime = False 
        return ExpatParser.feed(self, data, isFinal)


class ContentHandler(handler.ContentHandler):
    integers = ("start", "stop")
    strings = ("sequence_ac", 
               "sequence_id",
               "sequence_db",
               "signature_ac",
               "level",
               "level_tag")
    def __init__(self):
        self.element = []
    def startElement(self, name, attrs):
        self.element.append(name)
        self.content = ""
        if self.element==["matchset"]:
            self.record = Record()
            self.record.n_match = int(attrs["n_match"])
            self.record.n_seq = int(attrs["n_seq"])
        elif self.element==["matchset", "match"]:
            match = {}
            self.record.append(match)
    def endElement(self, name):
        assert name==self.element.pop()
        name = str(name)
        if self.element==["matchset", "match"]:
            match = self.record[-1]
            if name in ContentHandler.integers:
                match[name] = int(self.content)
            elif name in ContentHandler.strings:
                match[name] = self.content
            else:
                # Unknown type, treat it as a string
                match[name] = self.content
    def characters(self, content):
        self.content += content
