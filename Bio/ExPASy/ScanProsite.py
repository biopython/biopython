# Copyright 2009 by Michiel de Hoon. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Code for calling and parsing ScanProsite from ExPASy."""

# Importing these functions with leading underscore as not intended for reuse
from urllib.request import urlopen
from urllib.parse import urlencode

from xml.sax import handler
from xml.sax.expatreader import ExpatParser


class Record(list):
    """Represents search results returned by ScanProsite.

    This record is a list containing the search results returned by
    ScanProsite. The record also contains the data members n_match,
    n_seq, capped, and warning.
    """

    def __init__(self):
        """Initialize the class."""
        self.n_match = None
        self.n_seq = None
        self.capped = None
        self.warning = None


def scan(seq="", mirror="https://www.expasy.org", output="xml", **keywords):
    """Execute a ScanProsite search.

    Arguments:
     - mirror:   The ScanProsite mirror to be used
                 (default: https://www.expasy.org).
     - seq:      The query sequence, or UniProtKB (Swiss-Prot,
                 TrEMBL) accession
     - output:   Format of the search results
                 (default: xml)

    Further search parameters can be passed as keywords; see the
    documentation for programmatic access to ScanProsite at
    https://www.expasy.org/tools/scanprosite/ScanPrositeREST.html
    for a description of such parameters.

    This function returns a handle to the search results returned by
    ScanProsite. Search results in the XML format can be parsed into a
    Python object, by using the Bio.ExPASy.ScanProsite.read function.

    """
    parameters = {"seq": seq, "output": output}
    for key, value in keywords.items():
        if value is not None:
            parameters[key] = value
    command = urlencode(parameters)
    url = "%s/cgi-bin/prosite/PSScan.cgi?%s" % (mirror, command)
    handle = urlopen(url)
    return handle


def read(handle):
    """Parse search results returned by ScanProsite into a Python object."""
    content_handler = ContentHandler()
    saxparser = Parser()
    saxparser.setContentHandler(content_handler)
    saxparser.parse(handle)
    record = content_handler.record
    return record


# The classes below are considered private


class Parser(ExpatParser):
    """Process the result from a ScanProsite search (PRIVATE)."""

    def __init__(self):
        """Initialize the class."""
        ExpatParser.__init__(self)
        self.firsttime = True

    def feed(self, data, isFinal=0):
        """Raise an Error if plain text is received in the data.

        This is to show the Error messages returned by ScanProsite.
        """
        # Error messages returned by the ScanProsite server are formatted as
        # as plain text instead of an XML document. To catch such error
        # messages, we override the feed method of the Expat parser.
        # The error message is (hopefully) contained in the data that was just
        # fed to the parser.
        if self.firsttime:
            if data[:5].decode("utf-8") != "<?xml":
                raise ValueError(data)
        self.firsttime = False
        return ExpatParser.feed(self, data, isFinal)


class ContentHandler(handler.ContentHandler):
    """Process and fill in the records, results of the search (PRIVATE)."""

    integers = ("start", "stop")
    strings = (
        "sequence_ac",
        "sequence_id",
        "sequence_db",
        "signature_ac",
        "level",
        "level_tag",
    )

    def __init__(self):
        """Initialize the class."""
        self.element = []

    def startElement(self, name, attrs):
        """Define the beginning of a record and stores the search record."""
        self.element.append(name)
        self.content = ""
        if self.element == ["matchset"]:
            self.record = Record()
            self.record.n_match = int(attrs["n_match"])
            self.record.n_seq = int(attrs["n_seq"])
        elif self.element == ["matchset", "match"]:
            match = {}
            self.record.append(match)

    def endElement(self, name):
        """Define the end of the search record."""
        assert name == self.element.pop()
        name = str(name)
        if self.element == ["matchset", "match"]:
            match = self.record[-1]
            if name in ContentHandler.integers:
                match[name] = int(self.content)
            elif name in ContentHandler.strings:
                match[name] = self.content
            else:
                # Unknown type, treat it as a string
                match[name] = self.content

    def characters(self, content):
        """Store the record content."""
        self.content += content
