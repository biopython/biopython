# Copyright 2008 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Parser for XML results returned by NCBI's Entrez Utilities. This
parser is used by the read() function in Bio.Entrez, and is not intended
be used directly.
"""

# The question is how to represent an XML file as Python objects. Some
# XML files returned by NCBI look like lists, others look like dictionaries,
# and others look like a mix of lists and dictionaries.
#
# My approach is to classify each possible element in the XML as a plain
# string, an integer, a list, a dictionary, or a structure. The latter is a
# dictionary where the same key can occur multiple times; in Python, it is
# represented as a dictionary where that key occurs once, pointing to a list
# of values found in the XML file.
#
# The parser then goes through the XML and creates the appropriate Python
# object for each element. The different levels encountered in the XML are
# preserved on the Python side. So a subelement of a subelement of an element
# is a value in a dictionary that is stored in a list which is a value in
# some other dictionary (or a value in a list which itself belongs to a list
# which is a value in a dictionary, and so on). Attributes encountered in 
# the XML are stored as a dictionary in a member .attributes of each element,
# and the tag name is saved in a member .tag.
#
# To decide which kind of Python object corresponds to each element in the
# XML, the parser analyzes the DTD referred at the top of (almost) every
# XML file returned by the Entrez Utilities. This is preferred over a hand-
# written solution, since the number of DTDs is rather large and their
# contents may change over time. About half the code in this parser deals
# wih parsing the DTD, and the other half with the XML itself.


import os.path
from xml.parsers import expat

# The following four classes are used to add a member .attributes to integers,
# strings, lists, and dictionaries, respectively.

class IntegerElement(int): pass

class StringElement(str): pass

class UnicodeElement(unicode): pass

class ListElement(list): pass

class DictionaryElement(dict): pass

# A StructureElement is like a dictionary, but some of its keys can have
# multiple values associated with it. These values are stored in a list
# under each key.
class StructureElement(dict):
    def __init__(self, keys):
        dict.__init__(self)
        for key in keys:
            dict.__setitem__(self, key, [])
        self.listkeys = keys
    def __setitem__(self, key, value):
        if key in self.listkeys:
            self[key].append(value)
        else:
            dict.__setitem__(self, key, value)

class DataHandler:

    def __init__(self, dtd_dir):
        self.stack = []
	self.errors = []
	self.integers = []
	self.strings = []
	self.lists = []
        self.dictionaries = []
        self.structures = {}
        self.items = []
        self.dtd_dir = dtd_dir

    def run(self, handle):
        """Set up the parser and let it parse the XML results"""
        self.parser = expat.ParserCreate()
        self.parser.SetParamEntityParsing(expat.XML_PARAM_ENTITY_PARSING_ALWAYS)
        self.parser.StartElementHandler = self.startElement
        self.parser.EndElementHandler = self.endElement
        self.parser.CharacterDataHandler = self.characters
        self.parser.ExternalEntityRefHandler = self.external_entity_ref_handler
        self.parser.ParseFile(handle)
        self.parser = None
        return self.object

    def parse(self, handle):
        BLOCK = 1024
        self.parser = expat.ParserCreate()
        self.parser.SetParamEntityParsing(expat.XML_PARAM_ENTITY_PARSING_ALWAYS)
        self.parser.StartElementHandler = self.startElement
        self.parser.EndElementHandler = self.endElement
        self.parser.CharacterDataHandler = self.characters
        self.parser.ExternalEntityRefHandler = self.external_entity_ref_handler

        while True :

            #Read in another block of the file...
            text = handle.read(BLOCK)
            if not text:
                # We have reached the end of the XML file
                for record in self.object:
                    yield record
                self.parser.Parse("", True)
                self.parser = None
                return

            self.parser.Parse(text, False)        

            if not self.stack:
                # Haven't read enough from the XML file yet
                continue

            records = self.stack[0]
            while len(records) > 1: # Then the top record is finished
                try:
                    record = records[0]
                except TypeError:
                    raise ValueError, "The XML file does not represent a list. Please use Entrez.read instead of Entrez.parse"
                records[:] = records[1:]
                yield record

    def startElement(self, name, attrs):
        self.content = ""
        if name in self.lists:
            object = ListElement()
        elif name in self.dictionaries:
            object = DictionaryElement()
        elif name in self.structures:
            object = StructureElement(self.structures[name])
        elif name in self.items:	# Only appears in ESummary
            name = str(attrs["Name"]) # convert from Unicode
            del attrs["Name"]
            itemtype = str(attrs["Type"]) # convert from Unicode
            del attrs["Type"]
            if itemtype=="Structure":
                object = DictionaryElement()
            elif name in ("ArticleIds", "History"):
                object = StructureElement(["pubmed", "medline"])
            elif itemtype=="List":
                object = ListElement()
            else:
                object = StringElement()
            object.itemname = name
            object.itemtype = itemtype
        elif name in self.strings + self.errors + self.integers:
            self.attributes = attrs
            return
        else:
            # Element not found in DTD; this will not be stored in the record
            object = ""
        if object!="":
            object.tag = name
            if attrs:
                object.attributes = dict(attrs)
            if len(self.stack)!=0:
                current = self.stack[-1]
                try:
                    current.append(object)
                except AttributeError:
                    current[name] = object
        self.stack.append(object)

    def endElement(self, name):
        value = self.content
        if name in self.errors:
            if value=="":
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
            elif self.object.itemtype=="Integer":
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
            return
        value.tag = name
        if self.attributes:
            value.attributes = dict(self.attributes)
            del self.attributes
        current = self.stack[-1]
        try:
            current.append(value)
        except AttributeError:
            current[name] = value

    def characters(self, content):
        self.content += content

    def elementDecl(self, name, model):
        """This callback function is called for each element declaration:
        <!ELEMENT       name          (...)>
        encountered in a DTD. The purpose of this function is to determine
        whether this element should be regarded as a string, integer, list
        dictionary, structure, or error."""
        if name.upper()=="ERROR":
            self.errors.append(name)
            return
        if name=='Item' and model==(expat.model.XML_CTYPE_MIXED,
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
                            expat.model.XML_CTYPE_CHOICE)
          and model[1] in (expat.model.XML_CQUANT_NONE,
                           expat.model.XML_CQUANT_OPT)
          and len(model[3])==1):
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
            if name==None:
                if quantifier in (expat.model.XML_CQUANT_PLUS,
                                  expat.model.XML_CQUANT_REP):
                    for child in children:
                        multiple.append(child[2])
                else:
                    for child in children:
                        count(child)
            elif name.upper()!="ERROR":
                if quantifier in (expat.model.XML_CQUANT_NONE,
                                  expat.model.XML_CQUANT_OPT):
                    single.append(name)
                elif quantifier in (expat.model.XML_CQUANT_PLUS,
                                    expat.model.XML_CQUANT_REP):
                    multiple.append(name)
        count(model)
        if len(single)==0 and len(multiple)==1:
            self.lists.append(name)
        elif len(multiple)==0:
            self.dictionaries.append(name)
        else:
            self.structures.update({name: multiple})

    def external_entity_ref_handler(self, context, base, systemId, publicId):
        """The purpose of this function is to load the DTD locally, instead
        of downloading it from the URL specified in the XML. Using the local
        DTD results in much faster parsing. If the DTD is not found locally,
        we try to download it. In practice, this may fail though, if the XML
        relies on many interrelated DTDs. If new DTDs appear, putting them in
        Bio/Entrez/DTDs will allow the parser to see them."""
        location, filename = os.path.split(systemId)
        path = os.path.join(self.dtd_dir, filename)
        try:
            handle = open(path)
        except IOError:
            message = """\
Unable to load DTD file %s.

Bio.Entrez uses NCBI's DTD files to parse XML files returned by NCBI Entrez.
Though most of NCBI's DTD files are included in the Biopython distribution,
sometimes you may find that a particular DTD file is missing. In such a
case, you can download the DTD file from NCBI and install it manually.

Usually, you can find missing DTD files at either
    http://www.ncbi.nlm.nih.gov/dtd/
or
    http://eutils.ncbi.nlm.nih.gov/entrez/query/DTD/
If you cannot find %s there, you may also try to search
for it with a search engine such as Google.

Please save %s in the directory
%s
in order for Bio.Entrez to find it.
Alternatively, you can save %s in the directory
Bio/Entrez/DTDs in the Biopython distribution, and reinstall Biopython.

Please also inform the Biopython developers by sending an email to
biopython-dev@biopython.org to inform us about this missing DTD, so that we
can include it with the next release of Biopython.
""" % (filename, filename, filename, self.dtd_dir, filename)
            raise RuntimeError(message)
            
        parser = self.parser.ExternalEntityParserCreate(context)
        parser.ElementDeclHandler = self.elementDecl
        parser.ParseFile(handle)
        return 1






