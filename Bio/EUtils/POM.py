#!/usr/bin/env python

"""
This module implements the XML POM -- the Python Object Model for XML. It is
something like DOM, but more Python-ic, and easier to use. These base classes
are used to build POM source files which are self-validating python-based XML
constructor objects. The major parts of the dtd2py command line tool are also
here.

"""

import sys, os, re, string

try:
	True
except NameError:
	True = 1
	False = 0

class ValidationError(ValueError):
	"""ValidationError
	This exception is raised when an attempt is made to construct an XML POM
	tree that would be invalid.

	"""
	pass

#########################################################
# XML generating classes
# These classes are used to generate XML documents, similar to DOM. But, this
# interface is simpler and more Python-ic.
#########################################################

# plain text data to be added to a GenericNode.
# this class needs to emulate much of the ElementNode interface.
class IndentedText(str):
	def __init__(self, data=""):
		self.data = unescape(unicode(data))
		self._level = 0
		self._parent = None
	def set_text(self, data):
		self.data = unescape(unicode(data))
	def get_text(self):
		return self.data
	def insert(self, data):
		self.data = unescape(unicode(data)) + self.data
	def add_text(self,data):
		self.data = self.data + unescape(unicode(data))
	append = add_text
	def __str__(self):
		return "%s%s" % ("\t"*self._level, escape(self.data))
	def __unicode__(self):
		return u"%s%s" % ("\t"*self._level, escape(self.data))
	def __repr__(self):
		return "%s(%r)" % (self.__class__.__name__, escape(self.data))
	def set_level(self, level):
		self._level = level
#	def __len__(self):
#		return len(self.data)
#	def __getslice__(self, start, end):
#		return self.data[start:end]
#	def __setslice__(self, start, end, v):
#		self.data[start:end] = v
#	def __delslice__(self, start, end):
#		del self.data[start:end]
	def get_escape_length(self):
		return len(escape(self.data))
	def destroy(self):
		self.data = None
		self._parent = None
	def fullpath(self):
		if self._parent:
			return "%s = %r" % (self._parent.fullpath(), self.data)
		else:
			return `self.data`
	def matchpath(self, pe):
		return 0
	def has_children(self):
		return 0
	def has_attributes(self):
		return 0
	def full_repr(self):
		return repr(self)


class Text(IndentedText):
	def __str__(self):
		return escape(self.data)
	def __unicode__(self):
		return escape(self.data)

class Comment(IndentedText):
	def __init__(self, data=""):
		self.data = unicode(data)
		self._level = 0
		self._parent = None
	def __str__(self):
		return "%s<!-- %s -->" % ("\t"*self._level, self._fix(self.data))
	def __unicode__(self):
		return u"%s<!-- %s -->" % ("\t"*self._level, self._fix(self.data))
	def set_text(self, data):
		self.data = unicode(data)
	def get_text(self):
		return self.data
	def insert(self, data):
		self.data = unicode(data) + self.data
	def add_text(self,data):
		self.data = self.data + unicode(data)
	append = add_text
	
	def _fix(self, data):
		if data.find(u"--") != -1:
			data = data.replace(u"--", u"- ")
		return data



# abstract base class for generic XML node generation. 
# Create an XML node by subclassing this and defining allowed attribute names
# in ATTLIST. CONTENTMODEL holds the content specification from the DTD.
# Then name of the subclass should exactly match the name of the XML element.

class ElementNode:
	ATTLIST = None
	CONTENTMODEL = None
	_acquired = { "_indented":1, "_namespace":None } # default acquired values
	def __init__(self, **attribs):
		self._attribs = {}
		for key, value in attribs.items():
			if self._validate_attribute(key):
				self._attribs[key] = value
			else:
				raise ValidationError, "invalid attribute name for this element"
		self._children = []
		self._parent = None
		self._level = 0
		self._inline = 0
		#self._indented = 1 # may be acquired.
		#self._namespace = None # may be acquired.
		# you can force element names to a particular case, regardless of
		# subclass name. This is sometimes needed overcome clashes with Python
		# keyword names.
		self._name = self.__class__.__name__

 	# check if attribute name is defined for this element
	def _validate_attribute(self, name):
		if self.ATTLIST:
			for xmlattr in self.ATTLIST:
				if name == xmlattr.name:
					return True
		return False
	
	def _verify_attributes(self):
		if not self.ATTLIST:
			return None
		for attr in self.ATTLIST:
			aval = self._attribs.get(attr.name, None)
			if aval is None:
				if attr.a_decl == REQUIRED:
					raise ValidationError, "required attribute not present: " + attr.name
			else:
				attr.verify(aval)


	def get_parent(self):
		return self._parent
	
	def reparent(self, newparent):
		if self._parent:
			i = self._parent.index(self)
			del self._parent[i]
		newparent.append(self)
	
	def detach(self):
		self._parent = None
		self._level = 0

	def destroy(self):
		"""destroy() Remove this node and all child node references."""
		# remove parent _children list reference
		if self._parent:
			i = self._parent.index(self)
			del self._parent[i]
		self._parent = None
		for n in self._children:
			n.detach()
		self._children = None
	
	def set_level(self, level):
		self._level = int(level)

	def set_inline(self, tf=1):
		self._inline = not not tf # force to boolean

	def set_indented(self, tf=1):
		self._indented = not not tf # force to boolean

	def inherit_indent(self):
		"clears indentation flag so that it may be acquired from parent."
		try:
			del self.__dict__["_indented"]
		except KeyError:
			pass

	def set_namespace(self, ns):
		self._namespace = ns

	# some ugly stuff for case-insensitive XHTML
	def use_lowercase(self):
		self._name = self.__class__.__name__.lower()

	def use_uppercase(self):
		self._name = self.__class__.__name__.upper()
	
	def use_truecase(self):
		self._name = self.__class__.__name__

	def index(self, obj):
		objid = id(obj)
		i = 0
		for o in self._children:
			if id(o) == objid:
				return i
			i += 1
		raise ValueError, "ElementNode: Object not contained here."

	def append(self, obj):
		obj.set_level(self._level+1)
		obj._parent = self
		self._children.append(obj)
	
	def extend(self, objlist):
		for obj in objlist:
			self.append(obj)

	def insert(self, index, obj):
		obj.set_level(self._level+1)
		obj._parent = self
		self._children.insert(index, obj)

	def add(self, klass, **kwargs):
		obj = klass( *(), **kwargs)
		self.append(obj)
		return obj

	def get_children(self):
		return self._children[:]

	def __iter__(self):
		return iter(self._children)

	def add_text(self, text):
		"Adding text to elements is so common, there is a special method for it."
		if self.has_children() and isinstance(self._children[-1], IndentedText):
			self._children[-1].add_text(text)
		else:
			t = Text(text)
			t.set_level(0)
			self.append(t)
	
	def replace_text(self, text):
		if self._children:
			del self._children[-1]
		self.append(Text(text))

	def __len__(self):
		return len(self._children)
	
	# The truth is, we exist.
	def __nonzero__(self):
		return 1

	def hasAttributes(self):
		return len(self._attribs)
	has_attributes = hasAttributes
	
	def has_attribute(self, name):
		if name in self._attribs.keys():
			return 1
		else:
			return 0

	def attributes(self):
		return map(lambda o: o.name, self.ATTLIST)

	def has_children(self):
		return len(self._children)

	def set_attribute(self, name, val):
		"""set_attribute(name, value)
		This exists to set attributes that have names with illegal Python
		identifier characters.

		"""
		if self._validate_attribute(name):
			self._attribs[name] = val

	def get_attribute(self, name):
		"""get_attribute(name)
		This exists to set attributes that have names with illegal Python
		identifier characters.

		"""
		return self._attribs[name]

	def __setattr__(self, name, value):
		if self._validate_attribute(name):
			self._attribs[name] = value
		else:
			self.__dict__[name] = value

	# this plus the _parent and _acquired attributes implement "acquisiton", 
	# or run-time inheritance.
	def __getattr__(self, name):
		try:
			return self._attribs[name]
		except KeyError:
			pass
		try:
			return self._acquire(name)
		except:
			pass
		raise AttributeError, "AttributeError: %s has no attribute '%s'" % (self._name, name)

	def _acquire(self, name):
		if self._parent:
			try:
				return self._parent.__dict__[name]
			except KeyError:
				pass
			return self._parent._acquire(name)
		else:
			try:
				return self._acquired[name]
			except KeyError:
				pass
		raise AttributeError

	def __delattr__(self, name):
		del self._attribs[name]

	def _find_index(self, index):
		if type(index) is str:
			for i in xrange(len(self._children)):
				if self._children[i].matchpath(index):
					return i
			raise IndexError, "no elements match"
		else:
			return index

	def __getitem__(self, index):
		if type(index) is str:
			el =  self.get_element(index)
			if el is None:
				raise IndexError, "no item matches"
			else:
				return el
		else:
			return self._children[index]

	def get(self, index, default = None):
		if isinstance(index, str):
			el = self.get_element(index)
			if el is None:
				return default
			return el
		return self._children[index]

	def has_key(self, index):
		if isinstance(index, str):
			return self.get_element(index) is not None
		raise TypeError("Can only use has_key on a string")
			
	
	def __setitem__(self, index, obj):
		index = self._find_index(index)
		obj.set_level(self._level+1)
		obj._parent = self
		self._children[index] = obj
	
	def __delitem__(self, index):
		index = self._find_index(index)
#		self._children[index].destroy()
		del self._children[index]

	def __repr__(self):
		attrs = map(lambda t: '%s=%r' % t, self._attribs.items())
		return "%s(%s)" % (self.__class__, ", ".join(attrs))

	def __str__(self):
		self._verify_attributes()
		if not self.CONTENTMODEL or self.CONTENTMODEL.is_empty():
			return self._empty_str()
		else:
			return self._non_empty_str()

	def __unicode__(self):
		self._verify_attributes()
		if not self.CONTENTMODEL or self.CONTENTMODEL.is_empty():
			return self._empty_unistr()
		else:
			return self._non_empty_unistr()
	
	def full_repr(self):
		s = ["n%d = %r" % ( self._level, self)]
		s.append("n%d.set_level(%d)" % (self._level, self._level+1))
		for c in self._children:
			if not c.has_children():
				s.append("n%d.append(%r)" % (self._level, c))
			else:
				s.append(c.full_repr())
				s.append("n%d.append(n%d)" % (self._level, self._level+1))
				s.append("del n%d" % (self._level+1))
		return "\n".join(s)

	def _tabs(self):
		return "\t"*(self._level*self._indented)

	def _get_ns(self):
		return IF(self._namespace, "%s:" % self._namespace, "")

	def _non_empty_str(self):
		s = ["%s<%s%s%s>" % (self._tabs(), self._get_ns(), self._name, self._attr_str())]
		map(s.append, map(str, self._children))
		s.append("%s</%s%s>" % (IF(self._inline, "", self._tabs()), self._get_ns(), self._name))
		if self._inline:
			return "".join(s)
		else:
			return "\n".join(s)

	def _empty_str(self):
		return "%s<%s%s%s />" % (self._tabs(), self._get_ns(), self._name, self._attr_str())
	
	def _attr_str(self):
		attrs = map(lambda t: ' %s="%s"' % t, map(lambda t: (t[0], escape(str(t[1]))), filter(lambda t: t[1] is not None, self._attribs.items())))
		return "".join(attrs)

	def _non_empty_unistr(self):
		s = [u"%s<%s%s%s>" % (self._tabs(), self._get_ns(), self._name, self._attr_unistr())]
		map(s.append, map(unicode, self._children))
		s.append(u"%s</%s%s>" % (IF(self._inline, "", self._tabs()), self._get_ns(), self._name))
		if self._inline:
			return u"".join(s)
		else:
			return u"\n".join(s)

	def _empty_unistr(self):
		return u"%s<%s%s%s />" % (self._tabs(), self._get_ns(), self._name, self._attr_unistr())
	
	def _attr_unistr(self):
		attrs = map(lambda t: u' %s="%s"' % t, map(lambda t: (t[0], escape(unicode(t[1]))), filter(lambda t: t[1] is not None, self._attribs.items())))
		return u"".join(attrs)

	# methods for node path manipulation
	def pathname(self):
		"""pathname() returns the ElementNode as a string in xpath format."""
		if self._attribs:
			s = map(lambda i: "@%s='%s'" % (i[0],i[1]), self._attribs.items())
			return "%s[%s]" % (self.__class__.__name__, " and ".join(s))
		else:
			return self.__class__.__name__

	def fullpath(self):
		"""fullpath() returns the ElementNode's full path as a string in xpath format."""
		if self._parent:
			base = self._parent.fullpath()
		else:
			base = ""
		return "%s/%s" % (base, self.pathname() )

	def matchpath(self, pathelement):
		if "[" not in pathelement:
			return pathelement == self._name
		else:
			xpath_re = re.compile(r'(\w*)(\[.*])')
			mo = xpath_re.match(pathelement)
			if mo:
				name, match = mo.groups()
				match = match.replace("@", "self.")
				match = match.replace("=", "==")
				return (name == self._name and eval(match[1:-1]))
			else:
				raise ValueError, "ivalid path element"

	def find_elements(self, pathelement):
		rv = []
		for child in self._children:
			if child.matchpath(pathelement):
				rv.append(child)
		return rv
	
	def get_element(self, pathelement):
		for child in self._children:
			if child.matchpath(pathelement):
				return child
		return None

	def _find_node(self, eltype, collect=None):
		if collect is None:
			collection = []
		else:
			collection = collect # should be a list
		for el in self._children:
			if el.has_children():
				el._find_node(eltype, collection)
			if isinstance(el, eltype):
				collection.append(el)
		return collection

	def find(self, elclass, **attribs):
		for obj in self._children:
			if isinstance(obj, elclass):
				if self._attribs_match(obj, attribs):
					return obj
		return None

	def getall(self, elclass, depth=0, collect=None):
		if collect is None:
			rv = []
		else:
			rv = collect # should be a list
		for el in self._children:
			if isinstance(el, elclass):
				rv.append(el)
			if depth > 0:
				el.getall(elclass, depth-1, rv)
		return rv

	def _attribs_match(self, obj, attribdict):
		for tname, tval in attribdict.items():
			try:
				if getattr(obj, tname) != tval:
					return 0
			except AttributeError:
				return 0
		return 1

	def tostring(self):
		return "".join([x.get_text() for x in self.text()])

	# XPath-like functions
	def comment(self):
		return self._find_node(Comment)

	def text(self):
		return self._find_node(IndentedText)

	def processing_instruction(self):
		return self._find_node(ProcessingInstruction)
	
	def node(self):
		return self._find_node(ElementNode)



class Fragments(ElementNode):
	"""Fragments is a special holder class to hold 'loose' markup fragments.
	That is, bits of markup that don't have a common container.  It is invisible."""
	
	def __str__(self):
		s = []
		map(s.append, map(str, self._children))
		if self._inline:
			return "".join(s)
		else:
			return "\n".join(s)

	def __unicode__(self):
		s = []
		map(s.append, map(str, self._children))
		if self._inline:
			return u"".join(s)
		else:
			return u"\n".join(s)



# base class for whole POM documents, including Header.
class POMDocument:
	HEADER = '<?xml version="1.0" encoding="iso-8859-1"?>\n'
	def __init__(self, dtd=None):
		self.dtd = dtd
		self.root = None
		self.parser = None
		self.dirty = 0

	def __str__(self):
		return self.HEADER + str(self.root) + "\n"

	def __unicode__(self):
		return self.HEADER + unicode(self.root) + "\n"

	def set_dirty(self, val=1):
		self.dirty = val

	def get_parser(self, handlerclass=None, module=None):
		mod = module or self.dtd
		self.parser = get_parser(handlerclass, self._callback, mod)
		return self.parser
	
	def del_parser(self):
		self.parser = None

	def _callback(self, doc):
		self.root = doc
		self.dirty = 0
	
	def parse(self, url, handlerclass=None, module=None):
		mod = module or self.dtd
		if not self.parser:
			self.get_parser(handlerclass, mod)
		self.parser.parse(url)
		self.del_parser()
	
	def parseFile(self, fo, handlerclass=None, module=None):
		mod = module or self.dtd
		if not self.parser:
			self.get_parser(handlerclass, mod)
		self.parser.parseFile(fo)
		self.del_parser()

	def write_xmlfile(self, filename=None):
		filename = filename or self.filename
		if filename:
			fo = open(os.path.expanduser(filename), "w")
			try:
				fo.write(str(self))
			finally:
				fo.close()
		self.dirty = 0
	writefile = write_xmlfile

	def writefileobject(self, fo):
		fo.write(str(self))

	def get_document(self, filename, dtdmodule):
		self.get_parser(module=dtdmodule)
		self.parse(filename)
		self.filename = filename
	
	def getnode(self, path):
		"""getnode(path) Returns an ElementNode addressed by the path."""
		elements = path.split("/")
		while not elements[0]: # eat empty first element
			elements.pop(0)
		node = self.root
		pathelement = elements.pop(0)
		if node.matchpath(pathelement):
			while elements:
				pathelement = elements.pop(0)
				node = node.get_element(pathelement)
				if node is None:
					raise IndexError, "path element not found"
			return node
		else:
			raise IndexError, "first path element not found"

	def setnode(self, path, text):
		node = self.getnode(path)
		node.replace_text(text)
	
	def delnode(self, path):
		els = path.split("/")
		path, endnode = "/".join(els[:-1]), els[-1]
		node = self.getnode(path)
		del node[endnode]
	
	def addnode(self, basepath, newnode):
		node = self.getnode(basepath)
		node.append(newnode)

	def add_text(self, basepath, text):
		node = self.getnode(basepath)
		node.add_text(text)

	def _write_text(self, fo, node):
		for n in node:
			if isinstance(n, IndentedText):
				fo.write(n.fullpath())
				fo.write("\n")
			else:
				self._write_text(fo, n)
		
	def write_repr(self, fo):
		realfile = 0
		if type(fo) is str:
			fo = open(fo, "w")
			realfile = 1
		fo.write(self.root.full_repr())
		if realfile:
			fo.close()

	def read_repr(self, filename, localdict=None):
		localdict = localdict or {}
		execfile(filename, globals(), localdict)
		self.root = localdict["n0"]

	def write_paths(self, fileobject):
		realfile = 0
		if type(fileobject) is str:
			fileobject = open(fileobject, "w")
			realfile = 1
		self._write_text(fileobject, self.root)
		if realfile:
			fileobject.close()
	


# parses XML files into a POM object model. A callback function is then called 
# with this object model as a paramter.
class ObjectParserHandler:
	def __init__(self, callback, module=None):
		self.stack = []
		self.msg = None
		self.callback = callback # gets called when message fully parsed. The
		                         # argument is the toplevel message object.
		self.modules = []
		if module is not None:
			if type(module) is list:
				self.modules.extend(module)
			else:
				self.modules.append(module)

	def add_module(self, module):
		self.modules.append(module)

	def _get_class(self, name):
		klass = None
		for mod in self.modules:
			try:
				klass = getattr(mod, name)
			except AttributeError:
				continue
			if klass:
				return klass
		raise AttributeError


	def startDocument(self):
		self.stack = []

	def endDocument(self):
		if self.stack: # stack should be empty now
			raise ValidationError, "unbalanced document!"
		self.callback(self.msg)
		self.msg = None

	def startElement(self, name, atts):
		"Handle an event for the beginning of an element."
		try:
			klass = self._get_class(name)
		except AttributeError:
			raise ValidationError, "Undefined element tag: "+name
		attr = {} # atts is a instance with unicode keys.. must convert to str..
		def fixatts(t):
			attr[str(t[0])] = unescape(str(t[1]))
		map(fixatts, atts.items())
		obj = klass( *(), **attr)
		obj.set_level(len(self.stack))
		self.stack.append(obj)

	def endElement(self, name):
		"Handle an event for the end of an element."
		obj = self.stack.pop()
		if self.stack:
			self.stack[-1].append(obj)
		else:
			self.msg = obj

	def characters(self, text):
		if self.stack:
			text = text.strip()
			if text:
				self.stack[-1].append(Text(text))
		
	def ignorableWhitespace(self, ch, start, length):
		pass
	def processingInstruction(self, target, data):
		"Handle a processing instruction event."
		print "unhandled processing instruction:", target, data
	def setDocumentLocator(self, locator):
		"Receive an object for locating the origin of SAX document events."
		pass


def _default_parser_callback(obj):
	print obj

def get_parser(handlerclass=None, callback=None, module=None):
	import xml.sax
	hc = handlerclass or ObjectParserHandler
	cb = callback or _default_parser_callback
	mod = module or sys.modules[__name__]
	handler = hc(cb, mod)
	parser  = xml.sax.make_parser()
	parser.setContentHandler(handler)
	return parser

#from xml.parsers.xmlproc.xmlapp import DTDConsumer
def get_dtd_compiler(fo, mixinmodule=None, toupper=0):
	global sourcegen
	import sourcegen
	from xml.parsers.xmlproc.dtdparser import DTDParser
	generator = sourcegen.get_sourcefile(fo)
	dh = DTDConsumerForSourceGeneration(generator, mixinmodule, toupper)
	parser = DTDParser()
	parser.set_dtd_consumer(dh)
	return parser



# xml helper classes, used in both generation and operation
# The are instantiated during compilation to generate themselves. 
# Then, when imported by the user from the dtds package, are used normally.
class ContentModel:
	"""Represents and validates a content model.

	"""
	def __init__(self, rawmodel=None):
		self.model = rawmodel # XXX

	def __repr__(self):
		return "%s(%r)" % (self.__class__, self.model)

	def is_empty(self):
		return not self.model


class _ContentModelGenerator:
	"""_ContentModelGenerator(rawmodel)
	The DTD parser generated and final content model are so different that a
	different content model generator is used for this object.

	"""
	def __init__(self, rawmodel=None):
		tm_type = type(rawmodel)
		if tm_type is str:
			if rawmodel == "EMPTY":
				self.model = EMPTY
			elif rawmodel == "#PCDATA":
				self.model = PCDATA
			elif rawmodel == "ANY":
				self.model = ANY
			else:
				raise ValidationError, "ContentModelGenerator: unknown special type"
		elif tm_type is tuple:
			self.model = rawmodel # XXX
		elif tm_type is type(None):
			self.model = None
		else:
			raise RuntimeError, "unknown content model format"

	def __repr__(self):
		return "%s(%r)" % (ContentModel, self.model)


class Enumeration(list):
	pass
# XXX

class AttributeList(list):
	def __repr__(self):
		return "%s(%r)" % (self.__class__, self.data)
	def __str__(self):
		return " ".join(map(str, self.data))
	def __unicode__(self):
		return u" ".join(map(str, self.data))

class _AttributeType(str):
	def __repr__(self):
		return "%s('%s')" % (self.__class__.__name__, self)

class IDREFS(AttributeList):
	def add_ref(self, value):
		self.data.append(IDREF(value))

class ENTITIES(AttributeList):
	pass
class NMTOKENS(AttributeList):
	pass

class CDATA(_AttributeType):
	pass
class ID(_AttributeType):
	pass
class IDREF(_AttributeType):
	pass
class NMTOKEN(_AttributeType):
	pass
class ENTITY(_AttributeType):
	pass


PCDATA = Text
ANY = True
EMPTY = None

# enumerations
AT_CDATA = 1
AT_ID = 2
AT_IDREF = 3
AT_IDREFS = 4
AT_ENTITY = 5
AT_ENTITIES = 6
AT_NMTOKEN = 7
AT_NMTOKENS = 8

REQUIRED = 11
IMPLIED = 12
DEFAULT = 13
FIXED = 14

_ATTRTYPEMAP = {
	"CDATA": AT_CDATA,
	"ID": AT_ID,
	"IDREF": AT_IDREF,
	"IDREFS": AT_IDREFS,
	"ENTITY": AT_ENTITY,
	"ENTITIES": AT_ENTITIES,
	"NMTOKEN": AT_NMTOKEN,
	"NMTOKENS": AT_NMTOKENS
}

_ATTRCLASSMAP = {
	AT_CDATA: CDATA,
	AT_ID: ID,
	AT_IDREF: IDREF,
	AT_IDREFS: IDREFS,
	AT_ENTITY: ENTITY,
	AT_ENTITIES: ENTITIES,
	AT_NMTOKEN: NMTOKEN,
	AT_NMTOKENS: NMTOKENS
}

_DEFAULTMAP = {
	u'#REQUIRED': REQUIRED,
	u'#IMPLIED': IMPLIED,
	u'#DEFAULT': DEFAULT,
	u'#FIXED': FIXED,
}

class XMLAttribute:
	def __init__(self, name, a_type, a_decl, default=None):
		self.name = str(name)
		a_type_type = type(a_type)
		#a_decl_type = type(a_decl)
		if a_type_type is unicode: # from the parser
			self.a_type = _ATTRTYPEMAP.get(str(a_type), a_type)
#		elif a_type_type is tuple or a_type_type is list:
#			self.a_type = a_type # XXX
		elif a_type_type is int: # from the generated file
			self.a_type = _ATTRCLASSMAP.get(a_type, a_type)
		elif a_type_type is list:
			self.a_type = Enumeration(map(str, a_type))
		else:
			self.a_type = a_type
		# declaration
		# convert string to int value when generating, just use the int when using.
		self.a_decl = _DEFAULTMAP.get(a_decl, a_decl)
		self.default = default
		# save the type to speed verify
		self.a_type_type = type(self.a_type)

	def __repr__(self):
		return "%s(%r, %r, %r, %r)" % (self.__class__, self.name, self.a_type, self.a_decl, self.default)

	def verify(self, value):
		if type(self.a_type) is list:
			if value not in self.a_type:
				raise ValidationError, "Enumeration has wrong value. %s is not one of %r." % (value, self.a_type)
		
		
		


# this DTD parser consumer generates the Python source code from the DTD. 
class DTDConsumerForSourceGeneration:
	def __init__(self, generator, mixins=None, toupper=0):
		self.generator = generator
		self.elements = {}
		self.parameter_entities = {}
		self.general_entities = {}
		self.toupper = toupper # should element names be converted to all caps?
		self.mixins = mixins # should be a module object

	def dtd_start(self):
		print "Starting to parse DTD...",
		self.generator.add_comment("This file generated by a program. do not edit.")
		self.generator.add_import(sys.modules[__name__])
		if self.mixins:
			self.generator.add_import(self.mixins)

	def dtd_end(self):
		print "done parsing. Writing file."
		self.generator.write()
	
	def new_element_type(self, elem_name, elem_cont):
		"Receives the declaration of an element type."
		try:
			element = self.elements[elem_name]
		except KeyError:
			parents = [ElementNode]
			mixinname = "%sMixin" % ( elem_name )
			if self.mixins and hasattr(self.mixins, mixinname):
				parents.insert(0, getattr(self.mixins, mixinname))
			ch = self.generator.add_class(IF(self.toupper, elem_name.upper(), elem_name), tuple(parents))
			ch.add_attribute("CONTENTMODEL", _ContentModelGenerator(elem_cont))
			self.elements[elem_name] = ch
			
	def new_attribute(self, elem, attr, a_type, a_decl, a_def):
		"Receives the declaration of a new attribute."
		try:
			element = self.elements[elem]
		except KeyError:
			raise ValidationError, "attribute defined before element!"
		try:
			attlist = element.get_attribute("ATTLIST")
		except KeyError:
			element.add_attribute("ATTLIST", AttributeList())
			attlist = element.get_attribute("ATTLIST")
		attlist.append(XMLAttribute(attr, a_type, a_decl, a_def))

	def handle_comment(self, contents):
		"Receives the contents of a comment."
		self.generator.add_comment(contents)

	def new_parameter_entity(self,name,val):
		"Receives internal parameter entity declarations."
		# these are handled internally by the DTD parser. but.. save it anyway.
		self.parameter_entities[name] = val
	
	def new_external_pe(self, name, pubid, sysid):
		"Receives external parameter entity declarations."
		# these are handled internally by the DTD parser.
	
	def new_general_entity(self,name,val):
		"Receives internal general entity declarations."
		self.general_entities[name] = val
		# XXX do we need to handle this?
		#print "XXX general entity:"
		#print name, val

	def new_external_entity(self, ent_name, pub_id, sys_id, ndata):
		"""Receives external general entity declarations. 'ndata' is the
		empty string if the entity is parsed."""
		# XXX do we need to handle this?
		print "XXX external entity:"
		print ent_name, pub_id, sys_id, ndata

	def new_notation(self,name,pubid,sysid):
		"Receives notation declarations."
		# XXX do we need to handle this?
		print "XXX unhandled notation:",
		print name, pubid, sysid

	def handle_pi(self, target, data):
		"Receives the target and data of processing instructions."
		# XXX do we need to handle this?
		print "XXX unhandled PI:",
		print target, data

#########################################################
# Utility functions
#########################################################

def IF(test, tv, fv=None):
	if test:
		return tv
	else:
		return fv

def get_mod_file(sourcefilename):
	"""get_mod_file(sourcefilename)
	Converts a file name into a file name inside the dtds package. This file
	name is the destination for generated python files.
	"""
	import DTDs as dtds
	modname = os.path.splitext(os.path.split(sourcefilename)[1])[0]
	return os.path.join(dtds.__path__[0], modname.translate(string.maketrans("-.", "__"))+".py")


def _find_element(elname, modules):
	for mod in modules:
		try:
			return getattr(mod, elname)
		except AttributeError:
			continue
	return None

def _construct_node(name, modules):
	if "[" not in name:
		nc = _find_element(name, modules)
		if nc is None:
			raise ValidationError, "no such element name in modules"
		return nc() # node
	else:
		xpath_re = re.compile(r'(\w*)(\[.*])')
		mo = xpath_re.match(name)
		if mo:
			attdict = {}
			ename, attribs = mo.groups()
			nc = _find_element(ename, modules)
			if nc is None:
				raise ValidationError, "no such element name in modules"
			attribs = attribs[1:-1].split("and") # chop brackets and split on 'and'
			attribs = map(string.strip, attribs) # strip whitespace
			for att in attribs:                  # dict elememnts are name and vaue
				name, val = att.split("=")
				attdict[name[1:]] = val[1:-1]
		return nc( *(), **attdict)



def make_node(path, modules, value=None):
	"""make_Node(path, modules, [value])
	Makes a node or an XML fragment given a path, element module list, and an
	optional value.
	"""
	if type(modules) is not list:
		modules = [modules]
	pathelements = path.split("/")
	if not pathelements[0]: # delete possible empty root node
		del pathelements[0]
	rootnode = current = _construct_node(pathelements[0], modules)
	for element in pathelements[1:]:
		new = _construct_node(element, modules)
		current.append(new)
		current = new
	current.set_inline()
	if value is not None:
		current.add_text(value)
	return rootnode
	

def unescape(s):
	if '&' not in s:
		return s
	s = s.replace("&lt;", "<")
	s = s.replace("&gt;", ">")
#	s = s.replace("&apos;", "'")
	s = s.replace("&quot;", '"')
	s = s.replace("&amp;", "&") # Must be last
	return s

def escape(s):
	s = s.replace("&", "&amp;") # Must be first
	s = s.replace("<", "&lt;")
	s = s.replace(">", "&gt;")
#	s = s.replace("'", "&apos;")
	s = s.replace('"', "&quot;")
	return s

# self test
if __name__ == "__main__":
	import os
	FILE = os.path.join(os.environ["PAF_HOME"], "etc", "dtd", "WCSinvalidation.dtd")
	outfilename = get_mod_file(FILE)
	argc = len(sys.argv)
#	outfile = open(outfilename, "w")

	# note: running this script as __main__ will not generate valid source code. 
	# Use the dtd2py script for that.
	dtdp = get_dtd_compiler(sys.stdout)
	dtdp.parse_resource(FILE)
#	outfile.close()
	print Comment("some ------- comment-")
	print "+++++++"
	import dtds.pvsystem as pvs
	n = make_node("/pvsystem[@major='2' and @dot='0' and @minor='0']/pvac/httpOwsPort", pvs, 8080)
	print n
	print "+++++++"
	print make_node('/pvsystem[@major="2" and @minor="0" and @dot="0"]/globals/enableMonitor', pvs, "true")
	print "+++++++"
	print make_node('globals/enableMonitor', pvs, "true")
	print "+++++++"
	print make_node('enableMonitor', pvs, "true")
	print "+++++++"
	print make_node('enableMonitor', pvs)

