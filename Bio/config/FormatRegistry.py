# Copyright 2002 by Jeffrey Chang, Andrew Dalke.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This is based on some older code by Andrew Dalke.

"""Implements a Registry to store Martel-type format expressions.

Classes:
FormatRegistry   Holds Biopython formats in a dictionary-like interface.
FormatObject     Describes a Biopython file format.
FormatGroup      Describes a group of Biopython file formats.

"""
# Private Functions:
# _parses_file     Return whether an expression can parse a file.
# _parses_string   Return whether an expression can parse a string.
# _normalize_expression   Turn an expression or path into an expression.
# _load_first_existing    Return the first format that loads successfully.
# _load_expression        Load a Martel expression.
# _load_object            Load a Python object.

import string
import weakref
try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO
import operator

from Bio import StdHandler
from Bio.config.Registry import *
import Martel
from Martel import Parser

import _support


class FormatRegistry(Registry):
    """This implements a dictionary-like interface to Biopython file
    formats.

    Methods:
    find_builder  Find a builder that converts from a format to an object.
    find_writer   Find a writer that can write an object to a format.

    """
    def __init__(self, name, load_path=None,
                 builder_path="Bio.builders", writer_path="Bio.writers"):
        Registry.__init__(self, name, load_path=load_path)
        self._builder_path = builder_path
        self._writer_path = writer_path
        
    def normalize(self, name_or_format): # XXX appropriate?
        if isinstance(name_or_format, type("")):
            # It's a name
            return self[name_or_format]
        return name_or_format

    def _build_parent_path(self, format, visited=None):
        if visited is None:
            visited = {}
        if visited.has_key(format.name):
            return []
        format_list = [format]
        for parent in format._parents:
            format_list.extend(self._build_parent_path(parent, visited))
        return format_list

    def _build_child_path(self, format, visited=None):
        if visited is None:
            visited = {}
        if visited.has_key(format.name):
            return []
        format_list = [format]
        for child in getattr(format, 'objs', []):
            format_list.extend(self._build_child_path(parent, visited))
        return format_list
        
    def find_builder(self, from_format, to_io):
        # The directory of the builders is organized according to:
        # builders/io/format
        basemodulename = "%s.%s" % (self._builder_path, to_io.abbrev)

        # Search through the formats in the order of most specific to
        # most general.
        all_formats = self._build_parent_path(from_format)
        for format in all_formats:
            name = basemodulename + "." + format.abbrev
            module = _support.safe_load_module(name)
            if module is not None:
                break
        else:
            raise TypeError("Cannot find builder for %r" % to_io.abbrev)
        return module.make_builder()

    def find_writer(self, from_io, to_format, outfile):
        # The directory of the writers is organized according to:
        # writers/io/format
        basemodulename = "%s.%s" % (self._writer_path, from_io.abbrev)
        
        # Search through the formats in the order of most general to
        # most specific.
        all_formats = self._build_child_path(to_format)
        for format in all_formats:
            name = basemodulename + "." + format.abbrev
            module = _support.safe_load_module(name)
            if module is not None:
                break
        else:
            raise TypeError("Cannot find writer for %r" % from_io.abbrev)
        return module.make_writer(outfile)

formats = FormatRegistry("formats", "Bio.formatdefs")


class FormatObject(RegisterableObject):
    """This object stores Biopython file formats and provides methods
    to work on them.

    Methods:
    identify        Identify the format at a URL.
    identifyFile    Identify the format of a file.
    identifyString  Identify the format of a string.

    make_parser     Make a parser that can parse the format.
    make_iterator   Make an iterator over files of this format.

    """
    def __init__(self, name, expression, abbrev=None, doc=None,
                 filter=None, multirecord=1):
        """FormatObject(name, expression[, abbrev][, doc]
        [, filter][, multirecord])

        name is the name of the object, abbrev is an abbreviation for
        the name, and doc is some documentation describing the object.

        expression is a Martel.Expression that can parse this format.
        filter is an optional Martel.Expression that can be used to
        quickly determine whether some input is parseable by this
        format.

        multirecord is either 0/1 indicating whether this format can
        be used to parse multiple records.  By default, it is 1.

        """
        RegisterableObject.__init__(self, name, abbrev, doc)
        self.expression = _normalize_expression(expression)
        self.filter = _normalize_expression(filter) or self.expression
        self.filter = _support.make_cached_expression(self.filter)
        self.multirecord = operator.truth(multirecord)
        self._parser_cache = {}
        self._iterator_cache = {}
        self._parents = []
        
    def identifyFile(self, infile, debug_level=0):
        """S.identifyFile(infile[, debug_level]) -> FormatObject or None"""
        if _parses_file(self.filter, infile, debug_level):
            return self
        return None
    
    def identifyString(self, s, debug_level=0):
        """S.identifyString(s[, debug_level]) -> FormatObject or None"""
        if _parses_string(self.filter, s, debug_level):
            return self
        return None
    
    def identify(self, source, debug_level=0):
        """S.identify(source[, debug_level]) -> FormatObject or None"""
        source = ReseekFile.prepare_input_source(source)
        f = source.getCharacterStream() or source.getByteStream()
        return self.identifyFile(f, debug_level)

    def make_parser(self, select_names=None, debug_level=0):
        """S.make_parser([select_names][, debug_level]) -> parser"""
        if select_names is not None:
            select_names = list(select_names)
            select_names.sort()
            key = tuple(select_names), debug_level
        else:
            key = None, debug_level

        if not self._parser_cache.has_key(key):
            exp = self.expression
            if select_names is not None:
                exp = Martel.select_names(exp, select_names)
            p = exp.make_parser(debug_level = debug_level)
            self._parser_cache[key] = p
        return self._parser_cache[key].copy()
    
    def make_iterator(self, tag="record", select_names=None, debug_level=0):
        """S.make_iterator([tag][, select_names][, debug_level]) -> iterator"""
        if select_names is not None:
            select_names = list(select_names)
            select_names.sort()
            key = tuple(select_names), debug_level
        else:
            key = None, debug_level

        if not self._iterator_cache.has_key(key):
            exp = self.expression
            if select_names is not None:
                exp = Martel.select_names(exp, select_names)
            p = exp.make_iterator(tag, debug_level = debug_level)
            self._iterator_cache[key] = p
        return self._iterator_cache[key].copy()
    
class FormatGroup(RegisterableGroup):
    """This object holds a group of FormatObjects.

    Methods:
    identify        Identify the format at a URL.
    identifyFile    Identify the format of a file.
    identifyString  Identify the format of a string.

    """
    def __init__(self, name, abbrev=None, filter=None, multirecord=1):
        """FormatGroup(name[, abbrev][, filter][, multirecord])

        name is the name of the object, abbrev is an abbreviation for
        the name.

        filter is an optional Martel.Expression that can be used to
        quickly determine whether some input is parseable by this
        group.

        multirecord is either 0/1 indicating whether this format can
        be used to parse multiple records.  By default, it is 1.

        """
        RegisterableGroup.__init__(self, name, abbrev, None)
        self.filter = _normalize_expression(filter)
        if filter is not None:
            self.filter = _support.make_cached_expression(self.filter)
        self.multirecord = multirecord
        self._parents = []
        
    def identifyFile(self, infile, debug_level=0):
        """S.identifyFile(infile[, debug_level]) -> FormatObject or None"""
        # See if the filter test weeds things out
        if self.filter:
            if not _parses_file(self.filter, infile, debug_level):
                return None
        for obj in self.objs:
            format = obj.identifyFile(infile, debug_level=debug_level)
            if format is not None:
                return format
        return None

    def identifyString(self, s, debug_level=0):
        """S.identifyString(s[, debug_level]) -> FormatObject or None"""
        return self.identifyFile(StringIO(s), debug_level)

    def identify(self, source, debug_level=0):
        """S.identify(source[, debug_level]) -> FormatObject or None"""
        source = ReseekFile.prepare_input_source(source)
        f = source.getCharacterStream() or source.getByteStream()
        return self.identifyFile(f, debug_level)

    def add(self, obj, *args, **keywds):
        RegisterableGroup.add(self, obj, *args, **keywds)
        obj._parents.append(weakref.proxy(self))
        

def _parses_file(expression, infile, debug_level):
    # Return a boolean indicating whether expression can parse infile.
    parser = expression.make_parser(debug_level)
    handler = StdHandler.RecognizeHandler()
    parser.setErrorHandler(handler)
    parser.setContentHandler(handler)
    pos = infile.tell()
    try:
        try:
            parser.parseFile(infile)
        except Parser.ParserException:
            pass
    finally:
        infile.seek(pos)
    return handler.recognized
 
def _parses_string(expression, s, debug_level):
    return _parses_string(expression, StringIO(s), debug_level)

def _normalize_expression(expression_or_path):
    if expression_or_path is None:
        return None
    if type(expression_or_path) != type(""):
        return expression_or_path
    return _load_expression(expression_or_path)

def _load_expression(path):
    from Martel import Expression
    x = _load_object(path)
    if x is not None:
        if not isinstance(x, Expression.Expression):
            try:
                klass = x.__class__.__name__
            except AttributeError:
                klass = type(x)
                raise TypeError("%r should be a Martel Expression but " \
                                "is a %r" % (path, klass))
        return x
 
    # Expression not found; make a useful error message
    msg = "Could not find %r\n" % (path,)
    msg = msg + "(You may need to add the top-level module to the PYTHONPATH)"
    raise TypeError(msg)

def _load_object(path):
    terms = string.split(path, ".")
    s = terms[0]
    # Import all the needed modules
    # (Don't know which are modules and which are classes, so simply
    # stop when imports fail.)
    # The order of appends is correct, since the last element cannot
    # be a module.
    x = __import__(s)
    prev_term = s
    for term in terms[1:]:
        try:
            __import__(s)
        except SyntaxError, exc:
##            raise SyntaxError("%s during import of %r" % (exc, s)), \
##                  None, sys.exc_info()[2]
            raise
        except ImportError, exc:
            # This is the only way I know to tell if the module
            # could not be loaded because it doesn't exist.
            error_text = str(exc)
            if error_text.find("No module named %s" % prev_term) == -1:
                raise
            break
        if not term:
            raise TypeError("There's a '.' in the wrong place: %r" % \
                            (path,))
        s = s + "." + term
        prev_term = term

    # Get the requested object
    s = terms[0]
    for term in terms[1:]:
        try:
            x = getattr(x, term)
        except AttributeError:
            raise AttributeError("%s object (%r) has no attribute %r" % \
                                 (type(x).__name__, s, term))
        s = s + "." + term
    return x

def _load_first_existing(basemodulename, possible_formats):
    for format in possible_formats:
        try:
            module = _support.load_module(basemodulename + "." + format.abbrev)
        except ImportError, exc:
            # This is the only way I know to tell if the module
            # could not be loaded because it doesn't exist.
            error_text = str(exc)
            if error_text.find("No module named %s" % format.abbrev) == -1:
                raise
            continue
        return module
    return None
