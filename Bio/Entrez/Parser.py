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
from collections import Counter
from xml.parsers import expat
from io import BytesIO
import xml.etree.ElementTree as ET
from xml.sax.saxutils import escape

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


class Consumer(object):

    def __init__(self, name, attrs):
        """Create a do-nothing Consumer object."""
        return

    def startElementHandler(self, name, attrs):
        return False

    def endElementHandler(self, name):
        return False

    def consume(self, content):
        return

    def store(self, key, value):
        return

    @property
    def value(self):
        return


class ErrorConsumer(Consumer):

    def __init__(self, name, attrs):
        """Create a Consumer for ERROR messages in the XML data."""
        self.data = []

    def consume(self, content):
        self.data.append(content)

    @property
    def value(self):
        value = "".join(self.data)
        if value == "":
            return None
        else:
            raise RuntimeError(value)


class StringConsumer(Consumer):

    consumable = set()

    def __init__(self, name, attrs):
        """Create a Consumer for plain text elements in the XML data."""
        self.tag = name
        self.attributes = dict(attrs)
        self.data = []

    def startElementHandler(self, name, attrs, prefix=None):
        if prefix:
            key = "%s:%s" % (prefix, name)
        else:
            key = name
        if key in self.consumable:
            tag = "<%s" % name
            for key, value in attrs.items():
                tag += ' %s="%s"' % (key, value)
            tag += ">"
            self.data.append(tag)
            return True
        return False

    def endElementHandler(self, name, prefix=None):
        if prefix:
            key = "%s:%s" % (prefix, name)
        else:
            key = name
        if key in self.consumable:
            tag = "</%s>" % name
            self.data.append(tag)
            return True
        return False

    def consume(self, content):
        self.data.append(content)

    @property
    def value(self):
        value = "".join(self.data)
        # Convert Unicode strings to plain strings if possible
        try:
            value = StringElement(value)
        except UnicodeEncodeError:
            value = UnicodeElement(value)
        value.tag = self.tag
        if self.attributes:
            value.attributes = self.attributes
        return value


class IntegerConsumer(Consumer):

    def __init__(self, name, attrs):
        """Create a Consumer for integer elements in the XML data."""
        self.tag = name
        self.attributes = dict(attrs)
        self.data = []

    def consume(self, content):
        self.data.append(content)

    @property
    def value(self):
        value = int("".join(self.data))
        value = IntegerElement(value)
        value.tag = self.tag
        value.attributes = self.attributes
        return value


class ListConsumer(Consumer):

    keys = None

    def __init__(self, name, attrs):
        """Create a Consumer for list elements in the XML data."""
        data = ListElement()
        data.tag = name
        if attrs:
            data.attributes = dict(attrs)
        self.data = data

    def store(self, key, value):
        if self.keys is not None and key not in self.keys:
            raise ValueError("Unexpected item '%s' in list" % key)
        self.data.append(value)

    @property
    def value(self):
        return self.data


class DictionaryConsumer(Consumer):

    multiple = None

    def __init__(self, name, attrs):
        """Create a Consumer for dictionary elements in the XML data."""
        data = DictionaryElement()
        data.tag = name
        data.attributes = dict(attrs)
        for key in self.multiple:
            data[key] = []
        self.data = data

    def store(self, key, value):
        if key in self.multiple:
            self.data[key].append(value)
        else:
            self.data[key] = value

    @property
    def value(self):
        return self.data


def select_item_consumer(name, attrs):
    assert name == 'Item'
    name = str(attrs["Name"])  # convert from Unicode
    del attrs["Name"]
    itemtype = str(attrs["Type"])  # convert from Unicode
    del attrs["Type"]
    if itemtype == "Structure":
        cls = type(name,
                   (DictionaryConsumer,),
                   {"multiple": set()})
        consumer = cls(name, attrs)
    elif name in ("ArticleIds", "History"):
        cls = type(name,
                   (DictionaryConsumer,),
                   {"multiple": set(["pubmed", "medline"])})
        consumer = cls(name, attrs)
    elif itemtype == "List":
        # Keys are unknown in this case
        consumer = ListConsumer(name, attrs)
    elif itemtype == "Integer":
        consumer = IntegerConsumer(name, attrs)
    elif itemtype in ("String", "Unknown", "Date"):
        consumer = StringConsumer(name, attrs)
    else:
        raise ValueError("Unknown item type %s" % name)
    return consumer


class DataHandler(object):

    from Bio import Entrez
    global_dtd_dir = os.path.join(str(Entrez.__path__[0]), "DTDs")
    global_xsd_dir = os.path.join(str(Entrez.__path__[0]), "XSDs")
    local_dtd_dir = ''
    local_xsd_dir = ''

    del Entrez

    def __init__(self, validate, escape):
        """Create a DataHandler object."""
        self.dtd_urls = []
        self.classes = {}
        self.consumer = None
        self.validating = validate
        self.escaping = escape
        self.parser = expat.ParserCreate(namespace_separator=" ")
        self.parser.SetParamEntityParsing(expat.XML_PARAM_ENTITY_PARSING_ALWAYS)
        self.parser.XmlDeclHandler = self.xmlDeclHandler
        self.schema_namespace = None
        self.namespace_level = Counter()
        self.namespace_prefix = {}
        self._directory = None

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
            return self.record
        except AttributeError:
            if self.parser.StartElementHandler:
                # We saw the initial <!xml declaration, and expat didn't notice
                # any errors, so self.record should be defined. If not, this is
                # a bug.
                raise RuntimeError("Failed to parse the XML file correctly, possibly due to a bug in Bio.Entrez. Please contact the Biopython developers via the mailing list or GitHub for assistance.")
            else:
                # We did not see the initial <!xml declaration, so probably
                # the input data is not in XML format.
                raise NotXMLError("XML declaration not found")

    def parse(self, handle):
        BLOCK = 1024
        while True:
            # Read in another block of the file...
            text = handle.read(BLOCK)
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
            try:
                records = self.record
            except AttributeError:
                if self.parser.StartElementHandler:
                    # We saw the initial <!xml declaration, and expat
                    # didn't notice any errors, so self.record should be
                    # defined. If not, this is a bug.
                    raise RuntimeError("Failed to parse the XML file correctly, possibly due to a bug in Bio.Entrez. Please contact the Biopython developers via the mailing list or GitHub for assistance.")
                else:
                    # We did not see the initial <!xml declaration, so
                    # probably the input data is not in XML format.
                    raise NotXMLError("XML declaration not found")

            if not isinstance(records, list):
                raise ValueError("The XML file does not represent a list. Please use Entrez.read instead of Entrez.parse")

            while len(records) >= 1:  # Then the top record is finished
                record = records.pop(0)
                yield record

            if not text:
                sys.stdout.flush()
                self.parser = None
                if self.consumer:
                    # We have reached the end of the XML file
                    # No more XML data, but there is still some unfinished
                    # business
                    raise CorruptedXMLError("Premature end of XML stream")
                return

    def xmlDeclHandler(self, version, encoding, standalone):
        # XML declaration found; set the handlers
        self.parser.StartElementHandler = self.startElementHandler
        self.parser.EndElementHandler = self.endElementHandler
        if self.escaping:
            self.parser.CharacterDataHandler = self.characterDataHandlerEscape
        else:
            self.parser.CharacterDataHandler = self.characterDataHandlerRaw
        self.parser.ExternalEntityRefHandler = self.externalEntityRefHandler
        self.parser.StartNamespaceDeclHandler = self.startNamespaceDeclHandler
        self.parser.EndNamespaceDeclHandler = self.endNamespaceDeclHandler

    def startNamespaceDeclHandler(self, prefix, uri):
        if prefix == 'xsi':
            # This is an xml schema
            self.schema_namespace = uri
            self.parser.StartElementHandler = self.schemaHandler
        else:
            # Note that the DTD for MathML specifies a default attribute
            # that declares the namespace for each MathML element. This means
            # that MathML element in the XML has an invisible MathML namespace
            # declaration that triggers a call to startNamespaceDeclHandler
            # and endNamespaceDeclHandler. Therefore we need to count how often
            # startNamespaceDeclHandler and endNamespaceDeclHandler were called
            # to find out their first and last invocation for each namespace.
            self.namespace_level[prefix] += 1
            self.namespace_prefix[uri] = prefix

    def endNamespaceDeclHandler(self, prefix):
        if prefix != 'xsi':
            self.namespace_level[prefix] -= 1
            if self.namespace_level[prefix] == 0:
                for key, value in self.namespace_prefix.items():
                    if value == prefix:
                        break
                else:
                    raise RuntimeError("Failed to find namespace prefix")
                del self.namespace_prefix[key]

    def schemaHandler(self, name, attrs):
        # process the XML schema before processing the element
        key = "%s noNamespaceSchemaLocation" % self.schema_namespace
        schema = attrs[key]
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
        # continue handling the element
        self.startElementHandler(name, attrs)
        # reset the element handler
        self.parser.StartElementHandler = self.startElementHandler

    def startElementHandler(self, name, attrs):
        # check if the name is in a namespace
        prefix = None
        if self.namespace_prefix:
            try:
                uri, name = name.split()
            except ValueError:
                pass
            else:
                prefix = self.namespace_prefix[uri]
                if self.namespace_level[prefix] == 1:
                    attrs = {'xmlns': uri}
        # First, check if the current consumer can use the tag
        if self.consumer is not None:
            if prefix:
                consumed = self.consumer.startElementHandler(name, attrs, prefix)
            else:
                consumed = self.consumer.startElementHandler(name, attrs)
            if consumed:
                return
        cls = self.classes.get(name)
        if cls is None:
            # Element not found in DTD
            if self.validating:
                raise ValidationError(name)
            else:
                # this will not be stored in the record
                consumer = Consumer(name, attrs)
        else:
            consumer = cls(name, attrs)
        consumer.parent = self.consumer
        if self.consumer is None:
            # This is relevant only for Entrez.parse, not for Entrez.read.
            # If self.consumer is None, then this is the first start tag we
            # encounter, and it should refer to a list. Store this list in
            # the record attribute, so that Entrez.parse can iterate over it.
            # The record attribute will be set again at the last end tag;
            # However, it doesn't hurt to set it twice.
            value = consumer.value
            if value is not None:
                self.record = value
        self.consumer = consumer

    def endElementHandler(self, name):
        prefix = None
        if self.namespace_prefix:
            try:
                uri, name = name.split()
            except ValueError:
                pass
            else:
                prefix = self.namespace_prefix[uri]
        consumer = self.consumer
        # First, check if the current consumer can use the tag
        if consumer is not None:
            if prefix:
                consumed = consumer.endElementHandler(name, prefix)
            else:
                consumed = consumer.endElementHandler(name)
            if consumed:
                return
        self.consumer = consumer.parent
        value = consumer.value
        if self.consumer is None:
            self.record = value
        elif value is not None:
            name = value.tag
            self.consumer.store(name, value)

    def characterDataHandlerRaw(self, content):
        self.consumer.consume(content)

    def characterDataHandlerEscape(self, content):
        content = escape(content)
        self.consumer.consume(content)

    def parse_xsd(self, root):
        name = ""
        for child in root:
            is_dictionary = False
            multiple = []
            for element in child.getiterator():
                if "element" in element.tag:
                    if "name" in element.attrib:
                        name = element.attrib['name']
                if "attribute" in element.tag:
                    is_dictionary = True
                if "sequence" in element.tag:
                    for grandchild in element:
                        key = grandchild.attrib['ref']
                        multiple.append(key)
            if is_dictionary:
                bases = (DictionaryConsumer,)
                multiple = set(multiple)
                self.classes[name] = type(str(name),
                                          bases,
                                          {"multiple": multiple})
                is_dictionary = False
            else:
                self.classes[name] = ListConsumer

    def elementDecl(self, name, model):
        """Call a call-back function for each element declaration in a DTD.

        This is used for each element declaration in a DTD like::

            <!ELEMENT       name          (...)>

        The purpose of this function is to determine whether this element
        should be regarded as a string, integer, list, dictionary, structure,
        or error.
        """
        if name.upper() == "ERROR":
            self.classes[name] = ErrorConsumer
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
            self.classes[name] = select_item_consumer
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
            if model[1] == expat.model.XML_CQUANT_REP:
                tags = []
                children = model[3]
                for child in children:
                    tag = child[2]
                    tags.append(tag)
                    if tag in self.classes:
                        try:
                            keys = self.classes[tag].keys
                        except AttributeError:
                            continue
                        tags.extend(keys)
                bases = (StringConsumer, )
                self.classes[name] = type(str(name), bases, {'consumable': tags})
            else:
                self.classes[name] = StringConsumer
            return
        # List-type elements
        if (model[0] in (expat.model.XML_CTYPE_CHOICE,
                         expat.model.XML_CTYPE_SEQ) and
            model[1] in (expat.model.XML_CQUANT_PLUS,
                         expat.model.XML_CQUANT_REP)):
            children = model[3]
            if model[0] == expat.model.XML_CTYPE_SEQ:
                assert len(children) == 1
            keys = set([child[2] for child in children])
            bases = (ListConsumer,)
            self.classes[name] = type(str(name), bases, {'keys': keys})
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
            quantifier, key, children = model[1:]
            if key is None:
                if quantifier in (expat.model.XML_CQUANT_PLUS,
                                  expat.model.XML_CQUANT_REP):
                    for child in children:
                        multiple.append(child[2])
                else:
                    for child in children:
                        count(child)
            elif key.upper() != "ERROR":
                if quantifier in (expat.model.XML_CQUANT_NONE,
                                  expat.model.XML_CQUANT_OPT):
                    single.append(key)
                elif quantifier in (expat.model.XML_CQUANT_PLUS,
                                    expat.model.XML_CQUANT_REP):
                    multiple.append(key)
        count(model)
        if len(single) == 0 and len(multiple) == 1:
            keys = set(multiple)
            bases = (ListConsumer, )
            self.classes[name] = type(str(name), bases, {'keys': keys})
        else:
            multiple = set(multiple)
            bases = (DictionaryConsumer,)
            self.classes[name] = type(str(name),
                                      bases,
                                      {"multiple": multiple})

    def open_dtd_file(self, filename):
        self._initialize_directory()
        path = os.path.join(self.local_dtd_dir, filename)
        try:
            handle = open(path, "rb")
        except IOError:
            pass
        else:
            return handle
        path = os.path.join(self.global_dtd_dir, filename)
        try:
            handle = open(path, "rb")
        except IOError:
            pass
        else:
            return handle
        return None

    def open_xsd_file(self, filename):
        self._initialize_directory()
        path = os.path.join(self.local_xsd_dir, filename)
        try:
            handle = open(path, "rb")
        except IOError:
            pass
        else:
            return handle
        path = os.path.join(self.global_xsd_dir, filename)
        try:
            handle = open(path, "rb")
        except IOError:
            pass
        else:
            return handle
        return None

    def save_dtd_file(self, filename, text):
        self._initialize_directory()
        path = os.path.join(self.local_dtd_dir, filename)
        try:
            handle = open(path, "wb")
        except IOError:
            warnings.warn("Failed to save %s at %s" % (filename, path))
        else:
            handle.write(text)
            handle.close()

    def save_xsd_file(self, filename, text):
        self._initialize_directory()
        path = os.path.join(self.local_xsd_dir, filename)
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

    def _initialize_directory(self):
        """Initialize the local DTD/XSD directories (PRIVATE).

        Added to allow for custom directory (cache) locations,
        for example when code is deployed on AWS Lambda.
        """
        # If user hasn't set a custom cache location, initialize it.
        if self.directory is None:
            import platform
            if platform.system() == 'Windows':
                self.directory = os.path.join(os.getenv("APPDATA"), "biopython")
            else:  # Unix/Linux/Mac
                home = os.path.expanduser('~')
                self.directory = os.path.join(home, '.config', 'biopython')
                del home
            del platform
        # Create DTD local directory
        self.local_dtd_dir = os.path.join(self.directory, 'Bio', 'Entrez', 'DTDs')
        try:
            os.makedirs(self.local_dtd_dir)  # use exist_ok=True on Python >= 3.2
        except OSError as exception:
            # Check if local_dtd_dir already exists, and that it is a directory.
            # Trying os.makedirs first and then checking for os.path.isdir avoids
            # a race condition.
            if not os.path.isdir(self.local_dtd_dir):
                raise exception
        # Create XSD local directory
        self.local_xsd_dir = os.path.join(self.directory, 'Bio', 'Entrez', 'XSDs')
        try:
            os.makedirs(self.local_xsd_dir)  # use exist_ok=True on Python >= 3.2
        except OSError as exception:
            if not os.path.isdir(self.local_xsd_dir):
                raise exception

    @property
    def directory(self):
        return self._directory

    @directory.setter
    def directory(self, directory):
        """Allow user to set a custom directory, also triggering subdirectory initialization."""
        self._directory = directory
        self._initialize_directory()
