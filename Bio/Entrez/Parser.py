# Copyright 2008-2014 by Michiel de Hoon.  All rights reserved.
# Revisions copyright 2008-2015 by Peter Cock. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Parser for XML results returned by NCBI's Entrez Utilities.

This parser is used by the read() function in Bio.Entrez, and is not
intended be used directly.

The question is how to represent an XML file as Python objects. Some
XML files returned by NCBI look like lists, others look like dictionaries,
and others look like a mix of lists and dictionaries.

My approach is to classify each possible element in the XML as a plain
string, an integer, a list, a dictionary, or a structure. The latter is a
dictionary where the same key can occur multiple times; in Python, it is
represented as a dictionary where that key occurs once, pointing to a list
of values found in the XML file.

The parser then goes through the XML and creates the appropriate Python
object for each element. The different levels encountered in the XML are
preserved on the Python side. So a subelement of a subelement of an element
is a value in a dictionary that is stored in a list which is a value in
some other dictionary (or a value in a list which itself belongs to a list
which is a value in a dictionary, and so on). Attributes encountered in
the XML are stored as a dictionary in a member .attributes of each element,
and the tag name is saved in a member .tag.

To decide which kind of Python object corresponds to each element in the
XML, the parser analyzes the DTD referred at the top of (almost) every
XML file returned by the Entrez Utilities. This is preferred over a hand-
written solution, since the number of DTDs is rather large and their
contents may change over time. About half the code in this parser deals
with parsing the DTD, and the other half with the XML itself.
"""
import sys
import re
import os
import warnings
from xml.parsers import expat
from io import BytesIO
import xml.etree.ElementTree as ET

# Importing these functions with leading underscore as not intended for reuse
from Bio._py3k import urlopen as _urlopen
from Bio._py3k import urlparse as _urlparse
from Bio._py3k import unicode


# The following four classes are used to add a member .attributes to integers,
# strings, lists, and dictionaries, respectively.


class IntegerElement(int):
    def __repr__(self):
        text = int.__repr__(self)
        try:
            attributes = self.attributes
        except AttributeError:
            return text
        return "IntegerElement(%s, attributes=%s)" % (text, repr(attributes))


class StringElement(str):
    def __repr__(self):
        text = str.__repr__(self)
        try:
            attributes = self.attributes
        except AttributeError:
            return text
        return "StringElement(%s, attributes=%s)" % (text, repr(attributes))


class UnicodeElement(unicode):
    def __repr__(self):
        text = unicode.__repr__(self)
        try:
            attributes = self.attributes
        except AttributeError:
            return text
        return "UnicodeElement(%s, attributes=%s)" % (text, repr(attributes))


class ListElement(list):
    def __repr__(self):
        text = list.__repr__(self)
        try:
            attributes = self.attributes
        except AttributeError:
            return text
        return "ListElement(%s, attributes=%s)" % (text, repr(attributes))


class DictionaryElement(dict):
    def __repr__(self):
        text = dict.__repr__(self)
        try:
            attributes = self.attributes
        except AttributeError:
            return text
        return "DictElement(%s, attributes=%s)" % (text, repr(attributes))


# A StructureElement is like a dictionary, but some of its keys can have
# multiple values associated with it. These values are stored in a list
# under each key.
class StructureElement(dict):
    def __init__(self, keys):
        """Initialize the class."""
        dict.__init__(self)
        for key in keys:
            dict.__setitem__(self, key, [])
        self.listkeys = keys

    def __setitem__(self, key, value):
        if key in self.listkeys:
            self[key].append(value)
        else:
            dict.__setitem__(self, key, value)

    def __repr__(self):
        text = dict.__repr__(self)
        try:
            attributes = self.attributes
        except AttributeError:
            return text
        return "DictElement(%s, attributes=%s)" % (text, repr(attributes))


class NotXMLError(ValueError):
    def __init__(self, message):
        """Initialize the class."""
        self.msg = message

    def __str__(self):
        return "Failed to parse the XML data (%s). Please make sure that the input data are in XML format." % self.msg


class CorruptedXMLError(ValueError):
    def __init__(self, message):
        """Initialize the class."""
        self.msg = message

    def __str__(self):
        return "Failed to parse the XML data (%s). Please make sure that the input data are not corrupted." % self.msg


class ValidationError(ValueError):
    """XML tag found which was not defined in the DTD.

    Validating parsers raise this error if the parser finds a tag in the XML
    that is not defined in the DTD. Non-validating parsers do not raise this
    error. The Bio.Entrez.read and Bio.Entrez.parse functions use validating
    parsers by default (see those functions for more information).
    """

    def __init__(self, name):
        """Initialize the class."""
        self.name = name

    def __str__(self):
        return ("Failed to find tag '%s' in the DTD. To skip all tags that "
                "are not represented in the DTD, please call Bio.Entrez.read "
                "or Bio.Entrez.parse with validate=False." % self.name)


class DataHandler(object):

    import platform
    if platform.system() == 'Windows':
        directory = os.path.join(os.getenv("APPDATA"), "biopython")
    else:  # Unix/Linux/Mac
        home = os.path.expanduser('~')
        directory = os.path.join(home, '.config', 'biopython')
        del home
    local_dtd_dir = os.path.join(directory, 'Bio', 'Entrez', 'DTDs')
    local_xsd_dir = os.path.join(directory, 'Bio', 'Entrez', 'XSDs')
    del directory
    del platform
    try:
        os.makedirs(local_dtd_dir)  # use exist_ok=True on Python >= 3.2
    except OSError as exception:
        # Check if local_dtd_dir already exists, and that it is a directory.
        # Trying os.makedirs first and then checking for os.path.isdir avoids
        # a race condition.
        if not os.path.isdir(local_dtd_dir):
            raise exception
    try:
        os.makedirs(local_xsd_dir)  # use exist_ok=True on Python >= 3.2
    except OSError as exception:
        if not os.path.isdir(local_xsd_dir):
            raise exception

    from Bio import Entrez
    global_dtd_dir = os.path.join(str(Entrez.__path__[0]), "DTDs")
    global_xsd_dir = os.path.join(str(Entrez.__path__[0]), "XSDs")
    del Entrez

    def __init__(self, validate):
        """Initialize the class."""
        self.stack = []
        self.errors = []
        self.integers = []
        self.strings = []
        self.lists = []
        self.dictionaries = []
        self.structures = {}
        self.items = []
        self.dtd_urls = []
        self.validating = validate
        self.parser = expat.ParserCreate(namespace_separator=" ")
        self.parser.SetParamEntityParsing(expat.XML_PARAM_ENTITY_PARSING_ALWAYS)
        self.parser.XmlDeclHandler = self.xmlDeclHandler
        self.is_schema = False

    def read(self, handle):
        """Set up the parser and let it parse the XML results."""
        # HACK: remove Bio._py3k handle conversion, since the Entrez XML parser
        # expects binary data
        if handle.__class__.__name__ == 'EvilHandleHack':
            handle = handle._handle
        if handle.__class__.__name__ == 'TextIOWrapper':
            handle = handle.buffer
        if hasattr(handle, "closed") and handle.closed:
            # Should avoid a possible Segmentation Fault, see:
            # http://bugs.python.org/issue4877
            raise IOError("Can't parse a closed handle")
        if sys.version_info[0] >= 3:
            # Another nasty hack to cope with a unicode StringIO handle
            # since the Entrez XML parser expects binary data (bytes)
            from io import StringIO
            if isinstance(handle, StringIO):
                from Bio._py3k import _as_bytes
                handle = BytesIO(_as_bytes(handle.read()))
        try:
            self.parser.ParseFile(handle)
        except expat.ExpatError as e:
            if self.parser.StartElementHandler:
                # We saw the initial <!xml declaration, so we can be sure that
                # we are parsing XML data. Most likely, the XML file is
                # corrupted.
                raise CorruptedXMLError(e)
            else:
                # We have not seen the initial <!xml declaration, so probably
                # the input data is not in XML format.
                raise NotXMLError(e)
        try:
            return self.object
        except AttributeError:
            if self.parser.StartElementHandler:
                # We saw the initial <!xml declaration, and expat didn't notice
                # any errors, so self.object should be defined. If not, this is
                # a bug.
                raise RuntimeError("Failed to parse the XML file correctly, possibly due to a bug in Bio.Entrez. Please contact the Biopython developers at biopython-dev@biopython.org for assistance.")
            else:
                # We did not see the initial <!xml declaration, so probably
                # the input data is not in XML format.
                raise NotXMLError("XML declaration not found")

    def parse(self, handle):
        BLOCK = 1024
        while True:
            # Read in another block of the file...
            text = handle.read(BLOCK)
            if not text:
                # We have reached the end of the XML file
                if self.stack:
                    # No more XML data, but there is still some unfinished
                    # business
                    raise CorruptedXMLError("Premature end of XML stream")
                try:
                    for record in self.object:
                        yield record
                except AttributeError:
                    if self.parser.StartElementHandler:
                        # We saw the initial <!xml declaration, and expat
                        # didn't notice any errors, so self.object should be
                        # defined. If not, this is a bug.
                        raise RuntimeError("Failed to parse the XML file correctly, possibly due to a bug in Bio.Entrez. Please contact the Biopython developers at biopython-dev@biopython.org for assistance.")
                    else:
                        # We did not see the initial <!xml declaration, so
                        # probably the input data is not in XML format.
                        raise NotXMLError("XML declaration not found")
                self.parser.Parse("", True)
                self.parser = None
                return

            try:
                self.parser.Parse(text, False)
            except expat.ExpatError as e:
                if self.parser.StartElementHandler:
                    # We saw the initial <!xml declaration, so we can be sure
                    # that we are parsing XML data. Most likely, the XML file
                    # is corrupted.
                    raise CorruptedXMLError(e)
                else:
                    # We have not seen the initial <!xml declaration, so
                    # probably the input data is not in XML format.
                    raise NotXMLError(e)

            if not self.stack:
                # Haven't read enough from the XML file yet
                continue

            records = self.stack[0]
            if not isinstance(records, list):
                raise ValueError("The XML file does not represent a list. Please use Entrez.read instead of Entrez.parse")
            while len(records) > 1:  # Then the top record is finished
                record = records.pop(0)
                yield record

    def xmlDeclHandler(self, version, encoding, standalone):
        # XML declaration found; set the handlers
        self.parser.StartElementHandler = self.startElementHandler
        self.parser.EndElementHandler = self.endElementHandler
        self.parser.CharacterDataHandler = self.characterDataHandler
        self.parser.ExternalEntityRefHandler = self.externalEntityRefHandler
        self.parser.StartNamespaceDeclHandler = self.startNamespaceDeclHandler

    def startNamespaceDeclHandler(self, prefix, un):
        # This is an xml schema
        if "Schema" in un:
            self.is_schema = True
        else:
            raise NotImplementedError("The Bio.Entrez parser cannot handle XML data that make use of XML namespaces")

    def startElementHandler(self, name, attrs):
        # preprocessing the xml schema
        if self.is_schema:
            if len(attrs) == 1:
                schema = list(attrs.values())[0]
                handle = self.open_xsd_file(os.path.basename(schema))
                # if there is no local xsd file grab the url and parse the file
                if not handle:
                    handle = _urlopen(schema)
                    text = handle.read()
                    self.save_xsd_file(os.path.basename(schema), text)
                    handle.close()
                    self.parse_xsd(ET.fromstring(text))
                else:
                    self.parse_xsd(ET.fromstring(handle.read()))
                    handle.close()
        self.content = ""
        if name in self.lists:
            object = ListElement()
        elif name in self.dictionaries:
            object = DictionaryElement()
        elif name in self.structures:
            object = StructureElement(self.structures[name])
        elif name in self.items:  # Only appears in ESummary
            name = str(attrs["Name"])  # convert from Unicode
            del attrs["Name"]
            itemtype = str(attrs["Type"])  # convert from Unicode
            del attrs["Type"]
            if itemtype == "Structure":
                object = DictionaryElement()
            elif name in ("ArticleIds", "History"):
                object = StructureElement(["pubmed", "medline"])
            elif itemtype == "List":
                object = ListElement()
            else:
                object = StringElement()
            object.itemname = name
            object.itemtype = itemtype
        elif name in self.strings + self.errors + self.integers:
            self.attributes = attrs
            return
        else:
            # Element not found in DTD
            if self.validating:
                raise ValidationError(name)
            else:
                # this will not be stored in the record
                object = ""
        if object != "":
            object.tag = name
            if attrs:
                object.attributes = dict(attrs)
            if len(self.stack) != 0:
                current = self.stack[-1]
                try:
                    current.append(object)
                except AttributeError:
                    current[name] = object
        self.stack.append(object)

    def endElementHandler(self, name):
        value = self.content
        if name in self.errors:
            if value == "":
                return
            else:
                raise RuntimeError(value)
        elif name in self.integers:
            value = IntegerElement(value)
        elif name in self.strings:
            # Convert Unicode strings to plain strings if possible
            try:
                value = StringElement(value)
            except UnicodeEncodeError:
                value = UnicodeElement(value)
        elif name in self.items:
            self.object = self.stack.pop()
            if self.object.itemtype in ("List", "Structure"):
                return
            elif self.object.itemtype == "Integer" and value:
                value = IntegerElement(value)
            else:
                # Convert Unicode strings to plain strings if possible
                try:
                    value = StringElement(value)
                except UnicodeEncodeError:
                    value = UnicodeElement(value)
            name = self.object.itemname
        else:
            self.object = self.stack.pop()
            value = re.sub(r"[\s]+", "", value)
            if self.is_schema and value:
                self.object.update({'data': value})
            return
        value.tag = name
        if self.attributes:
            value.attributes = dict(self.attributes)
            del self.attributes
        current = self.stack[-1]
        if current != "":
            try:
                current.append(value)
            except AttributeError:
                current[name] = value

    def characterDataHandler(self, content):
        self.content += content

    def parse_xsd(self, root):
        is_dictionary = False
        name = ""
        for child in root:
            for element in child.getiterator():
                if "element" in element.tag:
                    if "name" in element.attrib:
                        name = element.attrib['name']
                if "attribute" in element.tag:
                    is_dictionary = True
            if is_dictionary:
                self.dictionaries.append(name)
                is_dictionary = False
            else:
                self.lists.append(name)

    def elementDecl(self, name, model):
        """Call a call-back function for each element declaration in a DTD.

        This is used for each element declaration in a DTD like::

            <!ELEMENT       name          (...)>

        The purpose of this function is to determine whether this element
        should be regarded as a string, integer, list, dictionary, structure,
        or error.
        """
        if name.upper() == "ERROR":
            self.errors.append(name)
            return
        if name == 'Item' and model == (expat.model.XML_CTYPE_MIXED,
                                        expat.model.XML_CQUANT_REP,
                                        None, ((expat.model.XML_CTYPE_NAME,
                                                expat.model.XML_CQUANT_NONE,
                                                'Item',
                                                ()
                                                ),
                                               )
                                        ):
            # Special case. As far as I can tell, this only occurs in the
            # eSummary DTD.
            self.items.append(name)
            return
        # First, remove ignorable parentheses around declarations
        while (model[0] in (expat.model.XML_CTYPE_SEQ,
                            expat.model.XML_CTYPE_CHOICE) and
               model[1] in (expat.model.XML_CQUANT_NONE,
                           expat.model.XML_CQUANT_OPT) and
               len(model[3]) == 1):
            model = model[3][0]
        # PCDATA declarations correspond to strings
        if model[0] in (expat.model.XML_CTYPE_MIXED,
                        expat.model.XML_CTYPE_EMPTY):
            self.strings.append(name)
            return
        # List-type elements
        if (model[0] in (expat.model.XML_CTYPE_CHOICE,
                         expat.model.XML_CTYPE_SEQ) and
            model[1] in (expat.model.XML_CQUANT_PLUS,
                         expat.model.XML_CQUANT_REP)):
            self.lists.append(name)
            return
        # This is the tricky case. Check which keys can occur multiple
        # times. If only one key is possible, and it can occur multiple
        # times, then this is a list. If more than one key is possible,
        # but none of them can occur multiple times, then this is a
        # dictionary. Otherwise, this is a structure.
        # In 'single' and 'multiple', we keep track which keys can occur
        # only once, and which can occur multiple times.
        single = []
        multiple = []
        # The 'count' function is called recursively to make sure all the
        # children in this model are counted. Error keys are ignored;
        # they raise an exception in Python.

        def count(model):
            quantifier, name, children = model[1:]
            if name is None:
                if quantifier in (expat.model.XML_CQUANT_PLUS,
                                  expat.model.XML_CQUANT_REP):
                    for child in children:
                        multiple.append(child[2])
                else:
                    for child in children:
                        count(child)
            elif name.upper() != "ERROR":
                if quantifier in (expat.model.XML_CQUANT_NONE,
                                  expat.model.XML_CQUANT_OPT):
                    single.append(name)
                elif quantifier in (expat.model.XML_CQUANT_PLUS,
                                    expat.model.XML_CQUANT_REP):
                    multiple.append(name)
        count(model)
        if len(single) == 0 and len(multiple) == 1:
            self.lists.append(name)
        elif len(multiple) == 0:
            self.dictionaries.append(name)
        else:
            self.structures.update({name: multiple})

    def open_dtd_file(self, filename):
        path = os.path.join(DataHandler.local_dtd_dir, filename)
        try:
            handle = open(path, "rb")
        except IOError:
            pass
        else:
            return handle
        path = os.path.join(DataHandler.global_dtd_dir, filename)
        try:
            handle = open(path, "rb")
        except IOError:
            pass
        else:
            return handle
        return None

    def open_xsd_file(self, filename):
        path = os.path.join(DataHandler.local_xsd_dir, filename)
        try:
            handle = open(path, "rb")
        except IOError:
            pass
        else:
            return handle
        path = os.path.join(DataHandler.global_xsd_dir, filename)
        try:
            handle = open(path, "rb")
        except IOError:
            pass
        else:
            return handle
        return None

    def save_dtd_file(self, filename, text):
        path = os.path.join(DataHandler.local_dtd_dir, filename)
        try:
            handle = open(path, "wb")
        except IOError:
            warnings.warn("Failed to save %s at %s" % (filename, path))
        else:
            handle.write(text)
            handle.close()

    def save_xsd_file(self, filename, text):
        path = os.path.join(DataHandler.local_xsd_dir, filename)
        try:
            handle = open(path, "wb")
        except IOError:
            warnings.warn("Failed to save %s at %s" % (filename, path))
        else:
            handle.write(text)
            handle.close()

    def externalEntityRefHandler(self, context, base, systemId, publicId):
        """Handle external entiry reference in order to cache DTD locally.

        The purpose of this function is to load the DTD locally, instead
        of downloading it from the URL specified in the XML. Using the local
        DTD results in much faster parsing. If the DTD is not found locally,
        we try to download it. If new DTDs become available from NCBI,
        putting them in Bio/Entrez/DTDs will allow the parser to see them.
        """
        urlinfo = _urlparse(systemId)
        # Following attribute requires Python 2.5+
        # if urlinfo.scheme=='http':
        if urlinfo[0] in ['http', 'https', 'ftp']:
            # Then this is an absolute path to the DTD.
            url = systemId
        elif urlinfo[0] == '':
            # Then this is a relative path to the DTD.
            # Look at the parent URL to find the full path.
            try:
                source = self.dtd_urls[-1]
            except IndexError:
                # Assume the default URL for DTDs if the top parent
                # does not contain an absolute path
                source = "http://www.ncbi.nlm.nih.gov/dtd/"
            else:
                source = os.path.dirname(source)
            # urls always have a forward slash, don't use os.path.join
            url = source.rstrip("/") + "/" + systemId
        else:
            raise ValueError("Unexpected URL scheme %r" % (urlinfo[0]))
        self.dtd_urls.append(url)
        # First, try to load the local version of the DTD file
        location, filename = os.path.split(systemId)
        handle = self.open_dtd_file(filename)
        if not handle:
            # DTD is not available as a local file. Try accessing it through
            # the internet instead.
            try:
                handle = _urlopen(url)
            except IOError:
                raise RuntimeError("Failed to access %s at %s" % (filename, url))
            text = handle.read()
            handle.close()
            self.save_dtd_file(filename, text)
            handle = BytesIO(text)

        parser = self.parser.ExternalEntityParserCreate(context)
        parser.ElementDeclHandler = self.elementDecl
        parser.ParseFile(handle)
        handle.close()
        self.dtd_urls.pop()
        return 1
