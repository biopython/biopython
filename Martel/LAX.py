"""A simple way to read lists of fields from flat XML records.

Many XML formats are very simple: all the fields are needed, there is
no tree hierarchy, all the text inside of the tags is used, and the
text is short (it can easily fit inside of memory).  SAX is pretty
good for this but it's still somewhat complicated to use.  DOM is
designed to handle tree structures so is a bit too much for a simple
flat data structure.

This module implements a new, simpler API, which I'll call LAX.  It
only works well when the elements are small and non-hierarchical.  LAX
has three callbacks.

  start() -- the first method called

  element(tag, attrs, text) -- called once for each element, after the

    element has been fully read.  (Ie, called when the endElement
    would be called.)  The 'tag' is the element name, the attrs is the
    attribute object that would be used in a startElement, and the
    text is all the text between the two tags.  The text is the
    concatenation of all the characters() calls.
    
  end() -- the last method called (unless there was an error)

LAX.LAX is an content handler which converts the SAX events to
LAX events.  Here is an example use:

  >>> from Martel import Word, Whitespace, Group, Integer, Rep1, AnyEol
  >>> format = Rep1(Group("line", Word("name") + Whitespace() +
  ...                             Integer("age")) + AnyEol())
  >>> parser = format.make_parser()
  >>>
  >>> from Martel import LAX
  >>> class PrintFields(LAX.LAX):
  ...     def element(self, tag, attrs, text):
  ...         print tag, "has", repr(text)
  ...
  >>> parser.setContentHandler(PrintFields())
  >>> text = "Maggie 3\nPorter 1\n"
  >>> parser.parseString(text)
  name has 'Maggie'
  age has '3'
  line has 'Maggie 3'
  name has 'Porter'
  age has '1'
  line has 'Porter 1'
  >>>

Callbacks take some getting used to.  Many people prefer an iterative
solution which returns all of the fields of a given record at one
time.  The default implementation of LAX.LAX helps this case.
The 'start' method initializes a local variable named 'groups', which
is dictionary.  When the 'element' method is called, the information
is added to groups; the key is the element name and the value is the
list of text strings.  It's a list because the same field name may
occur multiple times.

If you need the element attributes as well as the name, use the
LAX.LAXAttrs class, which stores a list of 2-ples (text, attrs)
instead of just the text.

For examples:

  >>> iterator = format.make_iterator("line")
  >>> for record in iterator.iterateString(text, LAX.LAX()):
  ...     print record.groups["name"][0], "is", record.groups["age"][0]
  ...
  Maggie is 3
  Porter is 1
  >>>

If you only want a few fields, you can pass the list to constructor,
as in:

  >>> lax = LAX.LAX(["name", "sequence"])
  >>>

"""

import string
from xml.sax import handler

# Used to simplify the check if 

class _IsIn:
    def __contains__(self, obj):
        return 1

class LAX(handler.ContentHandler, dict):
    def __init__(self, fields = None):
        handler.ContentHandler.__init__(self)
        dict.__init__(self)
        if fields is None:
            fields = _IsIn()
        self.__fields = fields

    def __getattr__(self, name):
        if name == "document":
            return self
        raise AttributeError(name)

    def uses_tags(self):
        if isinstance(self.__fields, _IsIn):
            return None
        return self.__fields

    
    def startDocument(self):
        self.__capture = []
        self.__expect = None
        self.__pos = 0
        self.start()
        
    def start(self):
        self.clear()

    def startElement(self, tag, attrs):
        if tag in self.__fields:
            self.__capture.append( (tag, attrs, [], self.__pos) )
            self.__expect = tag
            
    def characters(self, s):
        self.__pos += len(s)
        for term in self.__capture:
            term[2].append(s)
            
    def endElement(self, tag):
        if tag == self.__expect:
            cap, attrs, text_items, start = self.__capture.pop()
            self.element(tag, attrs, string.join(text_items, ""),
                         start, self.__pos)
            if self.__capture:
                self.__expect = self.__capture[-1][0]
            else:
                self.__expect = None

    def element(self, tag, attrs, text, startpos, endpos):
        self.setdefault(tag, []).append(text)
            
    def endDocument(self):
        if self.__capture:
            missing = []
            for term in self.__capture:
                missing.append(term[0])
            raise TypeError("Looking for endElements for %s" % \
                            string.join(missing, ","))
        self.end()

    def end(self):
        pass


# Also stores the attributes
class LAXAttrs(LAX):
    def element(self, tag, attrs, text, startpos, endpos):
        self.setdefault(tag, []).append( (text, attrs) )

# Stores attributes and positions
class ElementInfo:
    def __init__(self, text, attrs, startpos, endpos):
        self.text = text
        self.attrs = attrs
        self.startpos = startpos
        self.endpos = endpos

class LAXPositions(LAXAttrs):
    def element(self, tag, attrs, text, startpos, endpos):
        self.setdefault(tag, []).append(
            ElementInfo(text, attrs, startpos, endpos) )
