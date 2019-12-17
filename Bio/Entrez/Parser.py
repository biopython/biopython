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
from Bio._py3k import raise_from as _raise_from


# The following four classes are used to add a member .attributes to integers,
# strings, lists, and dictionaries, respectively.


class NoneElement:
    """NCBI Entrez XML element mapped to None."""

    def __init__(self, tag, attributes, key=None):
        """Create a NoneElement."""
        self.tag = tag
        if key is None:
            self.key = tag
        else:
            self.key = key
        self.attributes = attributes

    def __eq__(self, other):
        """Define equality with other None objects."""
        if other is None:
            return True
        elif other.__eq__(None):
            return True
        else:
            return False

    def __ne__(self, other):
        """Define non-equality."""
        if other is None:
            return False
        elif other.__eq__(None):
            return False
        else:
            return True

    def __repr__(self):
        """Return a string representation of the object."""
        try:
            attributes = self.attributes
        except AttributeError:
            return "NoneElement"
        return "NoneElement(attributes=%s)" % repr(attributes)


class IntegerElement(int):
    """NCBI Entrez XML element mapped to an integer."""

    def __new__(cls, value, tag, attributes, key=None):
        """Create an IntegerElement."""
        self = int.__new__(cls, value)
        self.tag = tag
        if key is None:
            self.key = tag
        else:
            self.key = key
        self.attributes = attributes
        return self

    def __repr__(self):
        """Return a string representation of the object."""
        text = int.__repr__(self)
        try:
            attributes = self.attributes
        except AttributeError:
            return text
        return "IntegerElement(%s, attributes=%s)" % (text, repr(attributes))


class StringElement(str):
    """NCBI Entrez XML element mapped to a string."""

    def __new__(cls, value, tag, attributes, key=None):
        """Create a StringElement."""
        self = str.__new__(cls, value)
        self.tag = tag
        if key is None:
            self.key = tag
        else:
            self.key = key
        self.attributes = attributes
        return self

    def __repr__(self):
        """Return a string representation of the object."""
        text = str.__repr__(self)
        attributes = self.attributes
        if not attributes:
            return text
        return "StringElement(%s, attributes=%s)" % (text, repr(attributes))


class UnicodeElement(unicode):
    """NCBI Entrez XML element mapped to a unicode string."""

    def __new__(cls, value, tag, attributes, key=None):
        """Create a UnicodeElement."""
        self = unicode.__new__(cls, value)
        self.tag = tag
        if key is None:
            self.key = tag
        else:
            self.key = key
        self.attributes = attributes
        return self

    def __repr__(self):
        """Return a string representation of the object."""
        text = unicode.__repr__(self)
        attributes = self.attributes
        if not attributes:
            return text
        return "UnicodeElement(%s, attributes=%s)" % (text, repr(attributes))


class ListElement(list):
    """NCBI Entrez XML element mapped to a list."""

    def __init__(self, tag, attributes, allowed_tags, key=None):
        """Create a ListElement."""
        self.tag = tag
        if key is None:
            self.key = tag
        else:
            self.key = key
        self.attributes = attributes
        self.allowed_tags = allowed_tags

    def __repr__(self):
        """Return a string representation of the object."""
        text = list.__repr__(self)
        attributes = self.attributes
        if not attributes:
            return text
        return "ListElement(%s, attributes=%s)" % (text, repr(attributes))

    def store(self, value):
        """Append an element to the list, checking tags."""
        key = value.key
        if self.allowed_tags is not None and key not in self.allowed_tags:
            raise ValueError("Unexpected item '%s' in list" % key)
        self.append(value)


class DictionaryElement(dict):
    """NCBI Entrez XML element mapped to a dictionaray."""

    def __init__(self, tag, attrs, allowed_tags, repeated_tags=None, key=None):
        """Create a DictionaryElement."""
        self.tag = tag
        if key is None:
            self.key = tag
        else:
            self.key = key
        self.attributes = attrs
        self.allowed_tags = allowed_tags
        self.repeated_tags = repeated_tags
        if repeated_tags:
            for key in repeated_tags:
                self[key] = []

    def __repr__(self):
        """Return a string representation of the object."""
        text = dict.__repr__(self)
        attributes = self.attributes
        if not attributes:
            return text
        return "DictElement(%s, attributes=%s)" % (text, repr(attributes))

    def store(self, value):
        """Add an entry to the dictionary, checking tags."""
        key = value.key
        tag = value.tag
        if self.allowed_tags is not None and tag not in self.allowed_tags:
            raise ValueError("Unexpected item '%s' in dictionary" % key)
        if self.repeated_tags and key in self.repeated_tags:
            self[key].append(value)
        else:
            self[key] = value


class NotXMLError(ValueError):
    """Failed to parse file as XML."""

    def __init__(self, message):
        """Initialize the class."""
        self.msg = message

    def __str__(self):
        """Return a string summary of the exception."""
        return (
            "Failed to parse the XML data (%s). Please make sure that the input data "
            "are in XML format." % self.msg
        )


class CorruptedXMLError(ValueError):
    """Corrupted XML."""

    def __init__(self, message):
        """Initialize the class."""
        self.msg = message

    def __str__(self):
        """Return a string summary of the exception."""
        return (
            "Failed to parse the XML data (%s). Please make sure that the input data "
            "are not corrupted." % self.msg
        )


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
        """Return a string summary of the exception."""
        return (
            "Failed to find tag '%s' in the DTD. To skip all tags that "
            "are not represented in the DTD, please call Bio.Entrez.read "
            "or Bio.Entrez.parse with validate=False." % self.name
        )


class DataHandler(object):
    """Data handler for parsing NCBI XML from Entrez."""

    from Bio import Entrez

    global_dtd_dir = os.path.join(str(Entrez.__path__[0]), "DTDs")
    global_xsd_dir = os.path.join(str(Entrez.__path__[0]), "XSDs")
    local_dtd_dir = ""
    local_xsd_dir = ""

    del Entrez

    def __init__(self, validate, escape):
        """Create a DataHandler object."""
        self.dtd_urls = []
        self.element = None
        self.level = 0
        self.data = []
        self.attributes = None
        self.allowed_tags = None
        self.strings = {}
        self.lists = {}
        self.dictionaries = {}
        self.items = set()
        self.errors = set()
        self.validating = validate
        self.parser = expat.ParserCreate(namespace_separator=" ")
        self.parser.SetParamEntityParsing(expat.XML_PARAM_ENTITY_PARSING_ALWAYS)
        self.parser.XmlDeclHandler = self.xmlDeclHandler
        self.schema_namespace = None
        self.namespace_level = Counter()
        self.namespace_prefix = {}
        self._directory = None
        if escape:
            self.characterDataHandler = self.characterDataHandlerEscape
        else:
            self.characterDataHandler = self.characterDataHandlerRaw

    def read(self, handle):
        """Set up the parser and let it parse the XML results."""
        # HACK: remove Bio._py3k handle conversion, since the Entrez XML parser
        # expects binary data
        if handle.__class__.__name__ == "EvilHandleHack":
            handle = handle._handle
        if handle.__class__.__name__ == "TextIOWrapper":
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
                _raise_from(CorruptedXMLError(e), None)
            else:
                # We have not seen the initial <!xml declaration, so probably
                # the input data is not in XML format.
                _raise_from(NotXMLError(e), None)
        try:
            return self.record
        except AttributeError:
            if self.parser.StartElementHandler:
                # We saw the initial <!xml declaration, and expat didn't notice
                # any errors, so self.record should be defined. If not, this is
                # a bug.
                _raise_from(
                    RuntimeError(
                        "Failed to parse the XML file correctly, possibly due to a bug "
                        "in Bio.Entrez. Please contact the Biopython developers via "
                        "the mailing list or GitHub for assistance."
                    ),
                    None,
                )
            else:
                # We did not see the initial <!xml declaration, so probably
                # the input data is not in XML format.
                _raise_from(NotXMLError("XML declaration not found"), None)

    def parse(self, handle):
        """Parse the XML in the given file handle."""
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
                    _raise_from(CorruptedXMLError(e), None)
                else:
                    # We have not seen the initial <!xml declaration, so
                    # probably the input data is not in XML format.
                    _raise_from(NotXMLError(e), None)
            try:
                records = self.record
            except AttributeError:
                if self.parser.StartElementHandler:
                    # We saw the initial <!xml declaration, and expat
                    # didn't notice any errors, so self.record should be
                    # defined. If not, this is a bug.
                    _raise_from(
                        RuntimeError(
                            "Failed to parse the XML file correctly, possibly due to a "
                            "bug in Bio.Entrez. Please contact the Biopython "
                            "developers via the mailing list or GitHub for assistance."
                        ),
                        None,
                    )
                else:
                    # We did not see the initial <!xml declaration, so
                    # probably the input data is not in XML format.
                    _raise_from(NotXMLError("XML declaration not found"), None)

            if not isinstance(records, list):
                raise ValueError(
                    "The XML file does not represent a list. Please use Entrez.read "
                    "instead of Entrez.parse"
                )

            if not text:
                break

            while len(records) >= 2:
                # Then the first record is finished, while the second record
                # is still a work in progress.
                record = records.pop(0)
                yield record

        # We have reached the end of the XML file
        self.parser = None
        if self.element is not None:
            # No more XML data, but there is still some unfinished business
            raise CorruptedXMLError("Premature end of XML stream")

        # Send out the remaining records
        for record in records:
            yield record

    def xmlDeclHandler(self, version, encoding, standalone):
        """Set XML handlers when an XML declaration is found."""
        self.parser.StartElementHandler = self.startElementHandler
        self.parser.CharacterDataHandler = self.characterDataHandler
        self.parser.ExternalEntityRefHandler = self.externalEntityRefHandler
        self.parser.StartNamespaceDeclHandler = self.startNamespaceDeclHandler
        self.parser.EndNamespaceDeclHandler = self.endNamespaceDeclHandler

    def startNamespaceDeclHandler(self, prefix, uri):
        """Handle start of an XML namespace declaration."""
        if prefix == "xsi":
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
            assert uri == "http://www.w3.org/1998/Math/MathML"
            assert prefix == "mml"

    def endNamespaceDeclHandler(self, prefix):
        """Handle end of an XML namespace declaration."""
        if prefix != "xsi":
            self.namespace_level[prefix] -= 1
            if self.namespace_level[prefix] == 0:
                for key, value in self.namespace_prefix.items():
                    if value == prefix:
                        break
                else:
                    raise RuntimeError("Failed to find namespace prefix")
                del self.namespace_prefix[key]

    def schemaHandler(self, name, attrs):
        """Process the XML schema (before processing the element)."""
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

    def startElementHandler(self, tag, attrs):
        """Handle start of an XML element."""
        if tag in self.items:
            assert tag == "Item"
            name = str(attrs["Name"])  # convert from Unicode
            itemtype = str(attrs["Type"])  # convert from Unicode
            del attrs["Type"]
            if itemtype == "Structure":
                del attrs["Name"]
                element = DictionaryElement(
                    name, attrs, allowed_tags=None, repeated_tags=None
                )
                parent = self.element
                element.parent = parent
                # For consistency with lists below, store the element here
                if parent is None:
                    self.record = element
                else:
                    parent.store(element)
                self.element = element
                self.parser.EndElementHandler = self.endElementHandler
                self.parser.CharacterDataHandler = self.skipCharacterDataHandler
            elif name in ("ArticleIds", "History"):
                del attrs["Name"]
                allowed_tags = None  # allowed tags are unknown
                repeated_tags = frozenset(["pubmed", "medline"])
                element = DictionaryElement(
                    tag,
                    attrs,
                    allowed_tags=allowed_tags,
                    repeated_tags=repeated_tags,
                    key=name,
                )
                parent = self.element
                element.parent = parent
                # For consistency with lists below, store the element here
                if parent is None:
                    self.record = element
                else:
                    parent.store(element)
                self.element = element
                self.parser.EndElementHandler = self.endElementHandler
                self.parser.CharacterDataHandler = self.skipCharacterDataHandler
            elif itemtype == "List":
                del attrs["Name"]
                allowed_tags = None  # allowed tags are unknown
                element = ListElement(tag, attrs, allowed_tags, name)
                parent = self.element
                element.parent = parent
                if self.element is None:
                    # Set self.record here to let Entrez.parse iterate over it
                    self.record = element
                else:
                    parent.store(element)
                self.element = element
                self.parser.EndElementHandler = self.endElementHandler
                self.parser.CharacterDataHandler = self.skipCharacterDataHandler
            elif itemtype == "Integer":
                self.parser.EndElementHandler = self.endIntegerElementHandler
                self.parser.CharacterDataHandler = self.characterDataHandler
                self.attributes = attrs
            elif itemtype in ("String", "Unknown", "Date", "Enumerator"):
                assert self.attributes is None
                self.attributes = attrs
                self.parser.StartElementHandler = self.startRawElementHandler
                self.parser.EndElementHandler = self.endStringElementHandler
                self.parser.CharacterDataHandler = self.characterDataHandler
            else:
                raise ValueError("Unknown item type %s" % name)
        elif tag in self.errors:
            self.parser.EndElementHandler = self.endErrorElementHandler
            self.parser.CharacterDataHandler = self.characterDataHandler
        elif tag in self.strings:
            self.parser.StartElementHandler = self.startRawElementHandler
            self.parser.EndElementHandler = self.endStringElementHandler
            self.parser.CharacterDataHandler = self.characterDataHandler
            assert self.allowed_tags is None
            self.allowed_tags = self.strings[tag]
            assert self.attributes is None
            self.attributes = attrs
        elif tag in self.dictionaries:
            allowed_tags, repeated_tags = self.dictionaries[tag]
            element = DictionaryElement(tag, attrs, allowed_tags, repeated_tags)
            parent = self.element
            element.parent = parent
            # For consistency with lists below, store the element here
            if parent is None:
                self.record = element
            else:
                parent.store(element)
            self.element = element
            self.parser.EndElementHandler = self.endElementHandler
            self.parser.CharacterDataHandler = self.skipCharacterDataHandler
        elif tag in self.lists:
            allowed_tags = self.lists[tag]
            element = ListElement(tag, attrs, allowed_tags)
            parent = self.element
            element.parent = parent
            if parent is None:
                # Set self.record here to let Entrez.parse iterate over it
                self.record = element
            else:
                parent.store(element)
            self.element = element
            self.parser.EndElementHandler = self.endElementHandler
            self.parser.CharacterDataHandler = self.skipCharacterDataHandler
        else:
            # Element not found in DTD
            if self.validating:
                raise ValidationError(tag)
            else:
                # this will not be stored in the record
                self.parser.StartElementHandler = self.startSkipElementHandler
                self.parser.EndElementHandler = self.endSkipElementHandler
                self.parser.CharacterDataHandler = self.skipCharacterDataHandler
                self.level = 1

    def startRawElementHandler(self, name, attrs):
        """Handle start of an XML raw element."""
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
                    attrs = {"xmlns": uri}
        if prefix:
            key = "%s:%s" % (prefix, name)
        else:
            key = name
        # self.allowed_tags is ignored for now. Anyway we know what to do
        # with this tag.
        tag = "<%s" % name
        for key, value in attrs.items():
            tag += ' %s="%s"' % (key, value)
        tag += ">"
        self.data.append(tag)
        self.parser.EndElementHandler = self.endRawElementHandler
        self.level += 1

    def startSkipElementHandler(self, name, attrs):
        """Handle start of an XML skip element."""
        self.level += 1

    def endStringElementHandler(self, tag):
        """Handle end of an XML string element."""
        element = self.element
        if element is not None:
            self.parser.StartElementHandler = self.startElementHandler
            self.parser.EndElementHandler = self.endElementHandler
            self.parser.CharacterDataHandler = self.skipCharacterDataHandler
        value = "".join(self.data)
        self.data = []
        attributes = self.attributes
        self.attributes = None
        if tag in self.items:
            assert tag == "Item"
            key = str(attributes["Name"])  # convert from Unicode
            del attributes["Name"]
        else:
            key = tag
        # Convert Unicode strings to plain strings if possible
        try:
            value = StringElement(value, tag, attributes, key)
        except UnicodeEncodeError:
            value = UnicodeElement(value, tag, attributes, key)
        if element is None:
            self.record = element
        else:
            element.store(value)
        self.allowed_tags = None

    def endRawElementHandler(self, name):
        """Handle start of an XML raw element."""
        self.level -= 1
        if self.level == 0:
            self.parser.EndElementHandler = self.endStringElementHandler
        if self.namespace_prefix:
            uri, name = name.split()
        tag = "</%s>" % name
        self.data.append(tag)

    def endSkipElementHandler(self, name):
        """Handle start of an XML skip element."""
        self.level -= 1
        if self.level == 0:
            self.parser.StartElementHandler = self.startElementHandler
            self.parser.EndElementHandler = self.endElementHandler

    def endErrorElementHandler(self, name):
        """Handle start of an XML error element."""
        if self.data:
            # error found:
            value = "".join(self.data)
            raise RuntimeError(value)
        # no error found:
        if self.element is not None:
            self.parser.EndElementHandler = self.endElementHandler
            self.parser.CharacterDataHandler = self.skipCharacterDataHandler

    def endElementHandler(self, name):
        """Handle end of an XML element."""
        element = self.element
        self.element = element.parent
        del element.parent

    def endIntegerElementHandler(self, tag):
        """Handle end of an XML integer element."""
        attributes = self.attributes
        self.attributes = None
        assert tag == "Item"
        key = str(attributes["Name"])  # convert from Unicode
        del attributes["Name"]
        if self.data:
            value = int("".join(self.data))
            self.data = []
            value = IntegerElement(value, tag, attributes, key)
        else:
            value = NoneElement(tag, attributes, key)
        element = self.element
        if element is None:
            self.record = value
        else:
            self.parser.EndElementHandler = self.endElementHandler
            self.parser.CharacterDataHandler = self.skipCharacterDataHandler
            if value is None:
                return
            element.store(value)

    def characterDataHandlerRaw(self, content):
        """Handle character data as-is (raw)."""
        self.data.append(content)

    def characterDataHandlerEscape(self, content):
        """Handle character data by encoding it."""
        content = escape(content)
        self.data.append(content)

    def skipCharacterDataHandler(self, content):
        """Handle character data by skipping it."""
        return

    def parse_xsd(self, root):
        """Parse an XSD file."""
        prefix = "{http://www.w3.org/2001/XMLSchema}"
        for element in root:
            isSimpleContent = False
            attribute_keys = []
            keys = []
            multiple = []
            assert element.tag == prefix + "element"
            name = element.attrib["name"]
            assert len(element) == 1
            complexType = element[0]
            assert complexType.tag == prefix + "complexType"
            for component in complexType:
                tag = component.tag
                if tag == prefix + "attribute":
                    # we could distinguish by type; keeping string for now
                    attribute_keys.append(component.attrib["name"])
                elif tag == prefix + "sequence":
                    maxOccurs = component.attrib.get("maxOccurs", "1")
                    for key in component:
                        assert key.tag == prefix + "element"
                        ref = key.attrib["ref"]
                        keys.append(ref)
                        if maxOccurs != "1" or key.attrib.get("maxOccurs", "1") != "1":
                            multiple.append(ref)
                elif tag == prefix + "simpleContent":
                    assert len(component) == 1
                    extension = component[0]
                    assert extension.tag == prefix + "extension"
                    assert extension.attrib["base"] == "xs:string"
                    for attribute in extension:
                        assert attribute.tag == prefix + "attribute"
                        # we could distinguish by type; keeping string for now
                        attribute_keys.append(attribute.attrib["name"])
                    isSimpleContent = True
            allowed_tags = frozenset(keys)
            if len(keys) == 1 and keys == multiple:
                assert not isSimpleContent
                self.lists[name] = allowed_tags
            elif len(keys) >= 1:
                assert not isSimpleContent
                repeated_tags = frozenset(multiple)
                self.dictionaries[name] = (allowed_tags, repeated_tags)
            else:
                self.strings[name] = allowed_tags

    def elementDecl(self, name, model):
        """Call a call-back function for each element declaration in a DTD.

        This is used for each element declaration in a DTD like::

            <!ELEMENT       name          (...)>

        The purpose of this function is to determine whether this element
        should be regarded as a string, integer, list, dictionary, structure,
        or error.
        """
        if name.upper() == "ERROR":
            self.errors.add(name)
            return
        if name == "Item" and model == (
            expat.model.XML_CTYPE_MIXED,
            expat.model.XML_CQUANT_REP,
            None,
            ((expat.model.XML_CTYPE_NAME, expat.model.XML_CQUANT_NONE, "Item", ()),),
        ):
            # Special case. As far as I can tell, this only occurs in the
            # eSummary DTD.
            self.items.add(name)
            return
        # First, remove ignorable parentheses around declarations
        while (
            model[0] in (expat.model.XML_CTYPE_SEQ, expat.model.XML_CTYPE_CHOICE)
            and model[1] in (expat.model.XML_CQUANT_NONE, expat.model.XML_CQUANT_OPT)
            and len(model[3]) == 1
        ):
            model = model[3][0]
        # PCDATA declarations correspond to strings
        if model[0] in (expat.model.XML_CTYPE_MIXED, expat.model.XML_CTYPE_EMPTY):
            if model[1] == expat.model.XML_CQUANT_REP:
                children = model[3]
                allowed_tags = frozenset(child[2] for child in children)
            else:
                allowed_tags = frozenset()
            self.strings[name] = allowed_tags
            return
        # List-type elements
        if model[0] in (
            expat.model.XML_CTYPE_CHOICE,
            expat.model.XML_CTYPE_SEQ,
        ) and model[1] in (expat.model.XML_CQUANT_PLUS, expat.model.XML_CQUANT_REP):
            children = model[3]
            if model[0] == expat.model.XML_CTYPE_SEQ:
                assert len(children) == 1
            allowed_tags = frozenset(child[2] for child in children)
            self.lists[name] = allowed_tags
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
                if quantifier in (
                    expat.model.XML_CQUANT_PLUS,
                    expat.model.XML_CQUANT_REP,
                ):
                    for child in children:
                        multiple.append(child[2])
                else:
                    for child in children:
                        count(child)
            elif key.upper() != "ERROR":
                if quantifier in (
                    expat.model.XML_CQUANT_NONE,
                    expat.model.XML_CQUANT_OPT,
                ):
                    single.append(key)
                elif quantifier in (
                    expat.model.XML_CQUANT_PLUS,
                    expat.model.XML_CQUANT_REP,
                ):
                    multiple.append(key)

        count(model)
        if len(single) == 0 and len(multiple) == 1:
            allowed_tags = frozenset(multiple)
            self.lists[name] = allowed_tags
        else:
            allowed_tags = frozenset(single + multiple)
            repeated_tags = frozenset(multiple)
            self.dictionaries[name] = (allowed_tags, repeated_tags)

    def open_dtd_file(self, filename):
        """Open specified DTD file."""
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
        """Open specified XSD file."""
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
        """Save DTD file to cache."""
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
        """Save XSD file to cache."""
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
        """Handle external entity reference in order to cache DTD locally.

        The purpose of this function is to load the DTD locally, instead
        of downloading it from the URL specified in the XML. Using the local
        DTD results in much faster parsing. If the DTD is not found locally,
        we try to download it. If new DTDs become available from NCBI,
        putting them in Bio/Entrez/DTDs will allow the parser to see them.
        """
        urlinfo = _urlparse(systemId)
        # Following attribute requires Python 2.5+
        # if urlinfo.scheme=='http':
        if urlinfo[0] in ["http", "https", "ftp"]:
            # Then this is an absolute path to the DTD.
            url = systemId
        elif urlinfo[0] == "":
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
                _raise_from(
                    RuntimeError("Failed to access %s at %s" % (filename, url)), None
                )
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

            if platform.system() == "Windows":
                self.directory = os.path.join(os.getenv("APPDATA"), "biopython")
            else:  # Unix/Linux/Mac
                home = os.path.expanduser("~")
                self.directory = os.path.join(home, ".config", "biopython")
                del home
            del platform
        # Create DTD local directory
        self.local_dtd_dir = os.path.join(self.directory, "Bio", "Entrez", "DTDs")
        try:
            os.makedirs(self.local_dtd_dir)  # use exist_ok=True on Python >= 3.2
        except OSError as exception:
            # Check if local_dtd_dir already exists, and that it is a directory.
            # Trying os.makedirs first and then checking for os.path.isdir avoids
            # a race condition.
            if not os.path.isdir(self.local_dtd_dir):
                _raise_from(exception, None)
        # Create XSD local directory
        self.local_xsd_dir = os.path.join(self.directory, "Bio", "Entrez", "XSDs")
        try:
            os.makedirs(self.local_xsd_dir)  # use exist_ok=True on Python >= 3.2
        except OSError as exception:
            if not os.path.isdir(self.local_xsd_dir):
                _raise_from(exception, None)

    @property
    def directory(self):
        """Directory for caching XSD and DTD files."""
        return self._directory

    @directory.setter
    def directory(self, directory):
        """Allow user to set a custom directory, also triggering subdirectory initialization."""
        self._directory = directory
        self._initialize_directory()
