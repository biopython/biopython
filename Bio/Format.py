"""Define a Format for the Bioformats projects.

These are the allowed fields.

 o name -- the primary, canonical name for the format (required)

This is the externally visible string that people should use.  It must
be unique.  The grammer for a name is not yet written formalized, but
here are some examples:

  swissprot
  swissprot/38
  pdb-generic
  pdb-generic/2.1
  pdb/XPLOR-3.1

 o abbrev -- an abbreviated name used internally (optional)

This value should not be visible to most people.  It is used as part
of the resolution system to find builders, converters, and writers
without having to define everything at the beginning.  The abbrev name
must be usable as a Python module name, which means it must be of the
form: [a-zA-Z_][a-zA-Z0-9_]*

If not given, the 'name' is used.

 o expression -- Martel expression to parse the format (optional)

An expression need to be defined at the bottom-most level of the
recognition tree.

 o filter -- Martel expression to test if the format matches (optional)

The object needs to support the 'make_parser' method.  The created
parser is tested against the input text.  If the parse is successful,
or generates a ParserIncompleteException, then the expression is said
to accept the text.

The filter must not accept empty input.

If not given, 'expression' is used.

 o multirecord -- specifies if the format returns more than one record (optional)

Some formats contain a single record, others returns a list of records.
Different parts of the API change to handle these as expected.  This value
should be 0 or 1.  The default is 1, which means there may be
more than one record in the file.

 o description -- free form text which describes this format (optional)


 == These fields concern the original provider of the format
 This may have nothing to do with the Bioformats format definition.

 o provider_url -- URL for the site providing the format (optional)
     Example: http://www.expasy.ch/

 o provider_doc -- URL for the primary documentation on the format (optional)
     Example: http://www.expasy.ch/sprot/userman.html

There may not be a permanent URL for all versions of a given
documentation.  Use the site you think is the most appropriate.


 == These fields concern the people who wrote this definition

 o author_name -- name of the person who wrote this definition (optional)
     Example: Andrew Dalke

(The name may be of a person, organization, non-corporeal being, etc.)

 o author_email -- email address of the person who wrote this definition (optional)
     Example: dalke@dalkescientific.com

 o author_url -- URL for the person who wrote this definition (optional)

In case you want people to know more about you.


 == These fields concern the people who wrote this definition

The intent is to allow the author fields to contain information about
the person who wrote the definition, and the maintainer fields to
contain information about who to bug when things go wrong.  They may
be different people.  It may even be that the author is a person in a
company, and the company is listed as the maintainer.

 o maintainer_name -- name of the person who maintains the definition (optional)

If not given, the default is the same as the author_name

 o maintainer_email -- email address of the person who maintains the definition (optional)

If not given, the default is the same as the author_email

 o maintainer_url -- URL for whoever wrote this format definition (optional)
     Example: http://www.dalkescientific.com/

If not given, the default is the same as author_url.

"""
import re, urllib, sys
from xml.sax import handler, saxutils
import Martel
from Martel import Parser

import ReseekFile, _FmtUtils, StdHandler


# This helps ensure the documentation and the list of allowed
# parameters stay in sync.

_allowed_terms = {
}
_required_terms = []
def _parse_docstring():
    open_paren = r"\("  # dealing with emacs cruft
    close_paren = r"\)"
    
    info_pattern = re.compile(\
                   r"^ o (?P<name>[a-zA-Z][a-zA-Z0-9_]*)\s*--\s*" + \
                   r"(?P<text>[^%s]+)" % close_paren + open_paren + \
                   r"(?P<required>optional|required)" + close_paren + \
                   r"\s*$")

    for line in re.split(r"[\r\n]+", __doc__):
        if line[:3] == " o ":
            m = info_pattern.match(line)
            if m is None:
                raise TypeError("Incorrect format: %r" % (line,))
            name = m.group("name")
            text = m.group("text")
            required = m.group("required")
            assert name is not None
            assert text is not None
            assert required in ("optional", "required")
            assert not _allowed_terms.has_key(name)
            if required == "required":
                _required_terms.append(name)
            _allowed_terms[name] = text

_parse_docstring()

_legal_abbrev = re.compile(r"[a-zA-Z][a-zA-Z0-9_]*$")
check_abbrev = _legal_abbrev.match

def _build_parent_path(format, visited, format_list):
    if not visited.has_key(format.name):
        visited[format.name] = 1
        format_list.append(format)
        for parent in format._parents:
            # depth first traversal is correct order
            _build_parent_path(parent, visited, format_list)

def _build_child_path(format, visited, format_list):
    if not visited.has_key(format.name):
        visited[format.name] = 1
        format_list.append(format)
        for child in getattr(format, "_children", ()):
            # depth first traversal is correct order
            _build_child_path(child, visited, format_list)

class Trace:
    def __init__(self, outfile = sys.stdout):
        self.outfile = outfile
        self._indent = 0
    def indent(self):
        self._indent += 1
    def dedent(self):
        self._indent -= 1
        assert self._indent >= 0, "too many dedents"
    def writeln(self, s):
        if self.outfile is None:
            return
        pre = "  " * self._indent
        for line in s.split("\n"):
            self.outfile.write(pre + line + "\n")
def _Trace(trace):
    if isinstance(trace, Trace):
        return trace
    if trace == 0:
        return Trace(None)
    return Trace()


def check_parser_file(expression, infile, trace, debug_level):
    trace.indent()
    if expression is None:
        trace.writeln("No filter defined")
        trace.dedent()
        return 1
    parser = expression.make_parser(debug_level = debug_level)
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
        
    if trace:
        if handler.recognized:
            trace.writeln("Check passed")
        else:
            trace.writeln("Check failed because:")
            trace.indent()
            trace.writeln(str(handler.exc))
            trace.dedent()

    trace.dedent()
    return handler.recognized
 
def check_parser_string(expression, s, debug_level):
    if expression is None:
        trace.writeln("No filter defined")
        return 1
    parser = expression.make_parser(debug_level = debug_level)
    handler = StdHandler.RecognizeHandler()
    parser.setContentHandler(handler)
    parser.setErrorHandler(handler)
    try:
        parser.parseString(s)
    except Parser.ParserException:
        pass
    
    if trace:
        if handler.recognized:
            trace.writeln("Check passed")
        else:
            trace.writeln("Check failed because:")
            trace.indent()
            trace.writeln(str(handler.exc))
            trace.dedent()
    
    return handler.recognized


class Format:
    def __init__(self, **kwargs):
        for k in kwargs.keys():
            if not _allowed_terms.has_key(k):
                raise TypeError("Do not understand keyword argument %r" % \
                                (k, ))
            pass
        for k, v in _allowed_terms.items():
            # Check for fields which aren't set.  Set optional fields to None
            if not kwargs.has_key(k):
                if k in _required_terms:
                    raise TypeError("The keyword argument %r is required" % \
                          (k,))
                kwargs[k] = None

        exp = kwargs["expression"]
        if isinstance(exp, type("")):
            kwargs["expression"] = _FmtUtils.ExpressionLoader(exp)

        exp = kwargs["filter"]
        if isinstance(exp, type("")):
            kwargs["filter"] = _FmtUtils.ExpressionLoader(exp)

        if kwargs["abbrev"] is None:
            kwargs["abbrev"] = kwargs["name"]
        abbrev = kwargs["abbrev"]
        if not check_abbrev(abbrev):
            raise TypeError("abbrev name of %r is not allowed" % \
                            (abbrev,))

        x = kwargs["multirecord"]
        if x == None:
            kwargs["multirecord"] = 1
        elif x not in (0, 1):
            raise TypeError("'multirecord' must be 0 or 1: %s" % repr(x))
        
        self.__dict__.update(kwargs)
        self._parents = []

    def identifyFile(self, infile, trace = 0, debug_level = 0):
        raise NotImplementedError("must be defined in subclass")
    
    def identifyString(self, s, trace = 0, debug_level = 0):
        raise NotImplementedError("must be defined in subclass")
    
    def identify(self, source, trace = 0, debug_level = 0):
        source = saxutils.prepare_input_source(source)
        # Is this correct?  Don't know - don't have Unicode exprerience
        f = source.getCharacterStream() or source.getByteStream()
        try:
            f.tell()
        except (AttributeError, IOError):
            f = ReseekFile.ReseekFile(f)
        return self.identifyFile(f, trace, debug_level)

    def _get_parents_in_depth_order(self):
        # Return a list of self and all the parents, in a depth-first
        # traversal.  There should be no duplicates in the list.
        # This is needed to find the most likely match during
        # best-match resolution in reading/building.
        visited = {}
        format_list = []
        visited[self.name] = 1
        format_list.append(self)
        for parent in self._parents:
            _build_parent_path(parent, visited, format_list)
        return format_list

    def _get_children_in_depth_order(self):
        # Return a list of self and all the parents, in a depth-first
        # traversal.  There should be no duplicates in the list.
        # This is needed to find the most likely match during
        # best-match resolution in writing.
        visited = {}
        format_list = []
        visited[self.name] = 1
        format_list.append(self)
        for child in getattr(self, "_children", ()):
            _build_child_path(child, visited, format_list)
        return format_list

class FormatDef(Format):
    def __init__(self, **kwargs):
        Format.__init__(self, **kwargs)
        assert self.expression is not None, \
               "must have an expression in a bottom-level format definition"
        if self.filter is None:
            # Filter must accept at least a record
            self.filter = self.expression
        self.filter = _FmtUtils.CacheParser(self.filter)
        self._parser_cache = {}
        self._iterator_cache = {}

    def identifyFile(self, infile, trace = 0, debug_level = 0):
        trace = _Trace(trace)
        if self.filter is self.expression:
            trace.writeln("%r: Checking expression" % (self.name,))
        else:
            trace.writeln("%r: Checking filter" % (self.name,))
        if check_parser_file(self.filter, infile, trace, debug_level):
            return self
        return None

    def identifyString(self, s, trace = 0, debug_level = 0):
        trace = _Trace(trace)
        if self.filter is self.expression:
            trace.writeln("%r: Checking expression", (self.name,))
        else:
            trace.writeln("%r: Checking filter", (self.name,))
        if check_parser_string(self.filter, s, trace, debug_level):
            return self
        return None

    def make_parser(self, select_names = None, debug_level = 0):
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
    
    def make_iterator(self, tag = "record", select_names = None,
                      debug_level = 0):
        if select_names is not None:
            select_names = list(select_names)
            select_names.sort()
            key = tuple(select_names), debug_level
        else:
            key = None, debug_level
            
        if not self._iterator_cache.has_key(key):
            exp = self.expression
            if select_names is None:
                exp = Martel.select_names(exp, select_names)
            p = exp.make_iterator(tag, debug_level = debug_level)
            self._parser_cache[key] = p
        return self._parser_cache[key].copy()
        

class FormatGroup(Format):
    def __init__(self, **kwargs):
        Format.__init__(self, **kwargs)
        assert self.expression is None, \
               "must not have an expression in a format grouping"
        if self.filter is not None:
            self.filter = _FmtUtils.CacheParser(self.filter)
        self._children = []

    def identifyFile(self, infile, trace = 0, debug_level = 0):
        # See if the filter test weeds things out
        trace = _Trace(trace)
        trace.writeln("%r: Checking format group filter" % (self.name,))
            
        if not check_parser_file(self.filter, infile, trace, debug_level):
            return None

        for (child, filter) in self._children:
            if filter is not None:
                # Check the associated pre-filter for each format
                trace.writeln("%r: Checking child %r filter" %
                              (self.name, child.name))
                if not check_parser_file(filter, infile, trace, debug_level):
                    continue
            # Ask the child if it knows about this format
            trace.writeln("%r: Asking child %r to identify" %
                          (self.name, child.name))
            trace.indent()
            format = child.identifyFile(infile, trace, debug_level)
            trace.dedent()
            if format is not None:
                trace.writeln("%r: Child %r thinks this is a %r" %
                              (self.name, child.name, format.name))
                return format
        return None
        

    def identifyString(self, s, trace = 0, debug_level = 0):
        # See if the filter test weeds things out
        if not check_parser_string(self.filter, s, trace, debug_level):
            return None
        
        for (child, filter) in self._children:
            if filter is not None:
                # Check the associated pre-filter for each format
                if not check_parser_file(filter, infile, trace, debug_level):
                    continue
            # Ask the child if it knows about this format
            format = child.identifyString(s, trace, debug_level)
            if format is not None:
                return format
        return None
