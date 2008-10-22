# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
This module provides code to work with the keywlist.txt file from
SwissProt.
http://www.expasy.ch/sprot/sprot-top.html


Classes:
Record            Stores the information about one keyword or one category
                  in the keywlist.txt file.

Functions:
parse             Parses the keywlist.txt file and returns an iterator to
                  the records it contains.

DEPRECATED:

Classes:
ListParser        Parses a keywlist.txt file into a list of keywords.

_Scanner          Scans the keywlist.txt file.
_ListConsumer     Consumes keywlist data to a list.


Functions:
extract_keywords  Return the keywords from a keywlist.txt file.

"""


class Record(dict):
    """
    This record stores the information of one keyword or category in the
    keywlist.txt as a Python dictionary. The keys in this dictionary are
    the line codes that can appear in the keywlist.txt file:

    ---------  ---------------------------     ----------------------
    Line code  Content                         Occurrence in an entry
    ---------  ---------------------------     ----------------------
    ID         Identifier (keyword)            Once; starts a keyword entry
    IC         Identifier (category)           Once; starts a category entry
    AC         Accession (KW-xxxx)             Once
    DE         Definition                      Once or more
    SY         Synonyms                        Optional; once or more
    GO         Gene ontology (GO) mapping      Optional; once or more
    HI         Hierarchy                       Optional; once or more
    WW         Relevant WWW site               Optional; once or more
    CA         Category                        Once per keyword entry; absent
                                               in category entries
    """
    def __init__(self):
        dict.__init__(self)
        for keyword in ("DE", "SY", "GO", "HI", "WW"):
            self[keyword] = []
    
def parse(handle):
    # First, skip the header
    for line in handle:
        if line.startswith("______________________________________"):
            break
    # Now parse the records
    record = Record()
    for line in handle:
        if line.startswith("-------------------------------------"):
            # We have reached the footer
            break
        key = line[:2]
        if key=="//":
            record["DE"] = " ".join(record["DE"])
            record["SY"] = " ".join(record["SY"])
            yield record
            record = Record()
        else:
            value = line[5:].strip()
            if key in ("ID", "IC", "AC", "CA"):
                record[key] = value
            elif key in ("DE", "SY", "GO", "HI", "WW"):
                record[key].append(value)
    # Read the footer and throw it away
    for line in handle:
        pass


# Everything below is deprecated.

from types import *

from Bio import File
from Bio.ParserSupport import *

class ListParser(AbstractParser):
    """Parses keywlist.txt data into a list of keywords.

    """
    def __init__(self):
        import warnings
        warnings.warn("Bio.SwissProt.KeyWList.ListParser is deprecated. Please use the function Bio.SwissProt.KeyWList.parse instead to parse the keywlist.txt file.  In case of any problems, please contact the Biopython developers (biopython-dev@biopython.org).",
              DeprecationWarning)
        self._scanner = _Scanner()
        self._consumer = _ListConsumer()

    def parse(self, handle):
        self._scanner.feed(handle, self._consumer)
        return self._consumer.keywords


class _Scanner:
    """Scan the keywlist.txt file included with the SwissProt distribution.

    Tested with:
    Release 37
    Release 38
    """

    def __init__(self):
        import warnings
        warnings.warn("Bio.SwissProt.KeyWList._Scanner is deprecated. Please use the function Bio.SwissProt.KeyWList.parse instead to parse the keywlist.txt file.  In case of any problems, please contact the Biopython developers (biopython-dev@biopython.org).",
              DeprecationWarning)

    def feed(self, handle, consumer):
        """feed(self, handle, consumer)

        Feed in the keywlist.txt file for scanning.  handle is a file-like
        object that contains keyword information.  consumer is a
        Consumer object that will receive events as the report is scanned.

        """
        if isinstance(handle, File.UndoHandle):
            uhandle = handle
        else:
            uhandle = File.UndoHandle(handle)
        
        self._scan_header(uhandle, consumer)
        self._scan_keywords(uhandle, consumer)
        self._scan_footer(uhandle, consumer)

    def _scan_header(self, uhandle, consumer):
        consumer.start_header()
        
        read_and_call(uhandle, consumer.noevent, start='----')
        read_and_call(uhandle, consumer.noevent, blank=1)
        read_and_call(uhandle, consumer.noevent, contains="SWISS-PROT")
        read_and_call(uhandle, consumer.noevent, contains="Release")
        read_and_call(uhandle, consumer.noevent, blank=1)
        read_and_call(uhandle, consumer.noevent, start='----')

        read_and_call(uhandle, consumer.noevent, blank=1)
        read_and_call(uhandle, consumer.noevent, start='List of keywords')
        read_and_call(uhandle, consumer.noevent, blank=1)
        read_and_call(uhandle, consumer.noevent, start='----')

        while 1:
            if attempt_read_and_call(uhandle, consumer.noevent, start='----'):
                break
            read_and_call(uhandle, consumer.noevent, blank=0)

        read_and_call(uhandle, consumer.noevent, start='Document name')
        read_and_call(uhandle, consumer.noevent, start='----')
        read_and_call(uhandle, consumer.noevent, blank=1)
        
        consumer.end_header()

    def _scan_keywords(self, uhandle, consumer):
        consumer.start_keywords()

        # SwissProt38 starts with lines:
        # Keyword
        # ______________________________________
        #
        # Check and see if it's release 38, and parse it.
        if attempt_read_and_call(uhandle, consumer.noevent, start='Keyword'):
            read_and_call(uhandle, consumer.noevent, start='____')

        while 1:
            if not attempt_read_and_call(uhandle, consumer.keyword, blank=0):
                break
        read_and_call(uhandle, consumer.noevent, blank=1)
        
        consumer.end_keywords()

    def _scan_footer(self, uhandle, consumer):
        consumer.start_footer()

        read_and_call(uhandle, consumer.noevent, start='----')
        while 1:
            if attempt_read_and_call(uhandle, consumer.noevent, start='----'):
                break
            read_and_call(uhandle, consumer.copyright, blank=0)

        consumer.end_footer()

class _ListConsumer(AbstractConsumer):
    """Consumer that converts a keywlist.txt file into a list of keywords.

    Members:
    keywords    List of keywords.

    """
    def __init__(self):
        import warnings
        warnings.warn("Bio.SwissProt.KeyWList._ListConsumer is deprecated. Please use the function Bio.SwissProt.KeyWList.parse instead to parse the keywlist.txt file.  In case of any problems, please contact the Biopython developers (biopython-dev@biopython.org).",
              DeprecationWarning)
        self.keywords = None

    def start_keywords(self):
        self.keywords = []

    def keyword(self, line):
        self.keywords.append(string.rstrip(line))

def extract_keywords(keywlist_handle):
    """extract_keywords(keywlist_handle) -> list of keywords

    Return the keywords from a keywlist.txt file.

    """
    import warnings
    warnings.warn("Bio.SwissProt.KeyWList.extract_keywords is deprecated. Please use the function Bio.SwissProt.KeyWList.parse instead to parse the keywlist.txt file.  In case of any problems, please contact the Biopython developers (biopython-dev@biopython.org).",
              DeprecationWarning)
    if type(keywlist_handle) is not FileType and \
       type(keywlist_handle) is not InstanceType:
        raise ValueError("I expected a file handle or file-like object")
    return ListParser().parse(keywlist_handle)
