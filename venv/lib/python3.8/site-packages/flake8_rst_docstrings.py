"""Check Python docstrings validate as reStructuredText (RST).

This is a plugin for the tool flake8 tool for checking Python
soucre code.
"""

import logging
import re
import sys
import textwrap

import tokenize as tk
from itertools import chain, dropwhile

try:
    from StringIO import StringIO
except ImportError:  # Python 3.0 and later
    from io import StringIO
    from io import TextIOWrapper

#####################################
# Start of backported tokenize code #
#####################################

# If possible (python >= 3.2) use tokenize.open to open files, so PEP 263
# encoding markers are interpreted.
try:
    tokenize_open = tk.open
except AttributeError:
    # Fall back on a backport of the encoding aware tokenize open function,
    # which requires we back port tokenize.detect_encoding to implement.
    from codecs import lookup, BOM_UTF8
    from io import open as io_open

    cookie_re = re.compile(r"^[ \t\f]*#.*?coding[:=][ \t]*([-\w.]+)")
    blank_re = re.compile(br"^[ \t\f]*(?:[#\r\n]|$)")

    # I don't think 'blank regular expression' is well named, think
    # it looks for blank line after any Python # comment removed.
    # Key test case of interest is hashbang lines!
    assert blank_re.match(b"\n")
    assert blank_re.match(b"# Comment\n")
    assert blank_re.match(b"#!/usr/bin/python\n")
    assert blank_re.match(b"#!/usr/bin/env python\n")
    assert not blank_re.match(b'"""Welcome\n')
    assert not blank_re.match(b'"""Welcome"""\n')

    def _get_normal_name(orig_enc):
        """Imitates get_normal_name in tokenizer.c (PRIVATE)."""
        # sys.stderr.write("DEBUG: _get_normal_name(%r)\n" % orig_enc)
        # Only care about the first 12 characters.
        enc = orig_enc[:12].lower().replace("_", "-")
        if enc == "utf-8" or enc.startswith("utf-8-"):
            return "utf-8"
        if enc in ("latin-1", "iso-8859-1", "iso-latin-1") or enc.startswith(
            ("latin-1-", "iso-8859-1-", "iso-latin-1-")
        ):
            return "iso-8859-1"
        return orig_enc

    def _find_cookie(line, filename, bom_found):
        """Find encoding string in a line of Python (PRIVATE)."""
        # sys.stderr.write("DEBUG: _find_cookie(%r, %r, %r)\n"
        #                  % (line, filename, bom_found))
        match = cookie_re.match(line)
        if not match:
            return None
        encoding = _get_normal_name(match.group(1))
        try:
            lookup(encoding)
        except LookupError:
            # This behaviour mimics the Python interpreter
            raise SyntaxError(
                "unknown encoding for {!r}: {}".format(filename, encoding)
            )

        if bom_found:
            if encoding != "utf-8":
                # This behaviour mimics the Python interpreter
                raise SyntaxError("encoding problem for {!r}: utf-8".format(filename))
            encoding += "-sig"
        return encoding

    def tokenize_open(filename):
        """Simulate opening a Python file read only with the correct encoding.

        While this was based on the Python 3 standard library function
        tokenize.open in order to backport it to Python 2.7, this proved
        painful.

        Note that because this text will later be fed into ``exex(...)`` we
        would hit SyntaxError encoding declaration in Unicode string, so the
        handle returned has the encoding line masked out!

        Note we don't just remove the line as that would throw off the line
        numbers, it is replaced with a Python comment.
        """
        # sys.stderr.write("DEBUG: tokenize_open(%r)\n" % filename)
        # Will check the first & second lines for an encoding
        # AND REMOVE IT FROM THE TEXT RETURNED
        with io_open(filename, "rb") as handle:
            lines = list(handle)

        # Find the encoding
        first = lines[0] if lines else b""
        second = lines[1] if len(lines) > 1 else b""
        default = "utf-8"
        bom_found = False
        if first.startswith(BOM_UTF8):
            bom_found = True
            first = first[3:]
            default = "utf-8-sig"
        encoding = _find_cookie(first, filename, bom_found)
        if encoding:
            lines[0] = "# original encoding removed\n"
        if not encoding and blank_re.match(first):
            # sys.stderr.write("DEBUG: Trying second line %r\n"
            #                  % second)
            encoding = _find_cookie(second, filename, bom_found)
            if encoding:
                lines[1] = "# original encoding removed\n"
        if not encoding:
            encoding = default

        # sys.stderr.write("DEBUG: tokenize_open using encoding=%r\n"
        #                  % encoding)

        # Apply the encoding, using StringIO as we removed the
        # original encoding to help legacy code using exec.
        # for b in lines:
        #     sys.stderr.write(b"DEBUG: " + b)
        return StringIO("".join(b.decode(encoding) for b in lines))


###################################
# End of backported tokenize code #
###################################

import restructuredtext_lint as rst_lint


__version__ = "0.0.13"


log = logging.getLogger(__name__)

rst_prefix = "RST"
rst_fail_load = 900
rst_fail_parse = 901
rst_fail_all = 902
rst_fail_lint = 903

# Level 1 - info
code_mapping_info = {
    "Possible title underline, too short for the title.": 1,
    "Unexpected possible title overline or transition.": 2,
}

# Level 2 - warning
code_mapping_warning = {
    # XXX ends without a blank line; unexpected unindent:
    "Block quote ends without a blank line; unexpected unindent.": 1,
    "Bullet list ends without a blank line; unexpected unindent.": 2,
    "Definition list ends without a blank line; unexpected unindent.": 3,
    "Enumerated list ends without a blank line; unexpected unindent.": 4,
    "Explicit markup ends without a blank line; unexpected unindent.": 5,
    "Field list ends without a blank line; unexpected unindent.": 6,
    "Literal block ends without a blank line; unexpected unindent.": 7,
    "Option list ends without a blank line; unexpected unindent.": 8,
    # Other:
    "Inline strong start-string without end-string.": 10,
    "Blank line required after table.": 11,
    "Title underline too short.": 12,
    "Inline emphasis start-string without end-string.": 13,
    "Inline literal start-string without end-string.": 14,
    "Inline interpreted text or phrase reference start-string without end-string.": 15,
    "Multiple roles in interpreted text (both prefix and suffix present; only one allowed).": 16,  # noqa: E501
    "Mismatch: both interpreted text role suffix and reference suffix.": 17,
    "Literal block expected; none found.": 18,
    "Inline substitution_reference start-string without end-string.": 19,
}

# Level 3 - error
code_mapping_error = {
    "Unexpected indentation.": 1,
    "Malformed table.": 2,
    # e.g. Unknown directive type "req".
    "Unknown directive type": 3,
    # e.g. Unknown interpreted text role "need".
    "Unknown interpreted text role": 4,
    # e.g. Undefined substitution referenced: "dict".
    "Undefined substitution referenced:": 5,
    # e.g. Unknown target name: "license_txt".
    "Unknown target name:": 6,
}

# Level 4 - severe
code_mapping_severe = {"Unexpected section title.": 1}

code_mappings_by_level = {
    1: code_mapping_info,
    2: code_mapping_warning,
    3: code_mapping_error,
    4: code_mapping_severe,
}


def code_mapping(level, msg, extra_directives, extra_roles, default=99):
    """Return an error code between 0 and 99."""
    try:
        return code_mappings_by_level[level][msg]
    except KeyError:
        pass
    # Following assumes any variable messages take the format
    # of 'Fixed text "variable text".' only:
    # e.g. 'Unknown directive type "req".'
    # ---> 'Unknown directive type'
    # e.g. 'Unknown interpreted text role "need".'
    # ---> 'Unknown interpreted text role'
    if msg.count('"') == 2 and ' "' in msg and msg.endswith('".'):
        txt = msg[: msg.index(' "')]
        value = msg.split('"', 2)[1]
        if txt == "Unknown directive type" and value in extra_directives:
            return 0
        if txt == "Unknown interpreted text role" and value in extra_roles:
            return 0
        return code_mappings_by_level[level].get(txt, default)
    return default


####################################
# Start of code copied from PEP257 #
####################################

# This is the reference implementation of the alogrithm
# in PEP257 for removing the indentation of a docstring,
# which has been placed in the public domain.
#
# This includes the minor change from sys.maxint to
# sys.maxsize for Python 3 compatibility.
#
# https://www.python.org/dev/peps/pep-0257/#handling-docstring-indentation


def trim(docstring):
    """PEP257 docstring indentation trim function."""
    if not docstring:
        return ""
    # Convert tabs to spaces (following the normal Python rules)
    # and split into a list of lines:
    lines = docstring.expandtabs().splitlines()
    # Determine minimum indentation (first line doesn't count):
    indent = sys.maxsize
    for line in lines[1:]:
        stripped = line.lstrip()
        if stripped:
            indent = min(indent, len(line) - len(stripped))
    # Remove indentation (first line is special):
    trimmed = [lines[0].strip()]
    if indent < sys.maxsize:
        for line in lines[1:]:
            trimmed.append(line[indent:].rstrip())
    # Strip off trailing and leading blank lines:
    while trimmed and not trimmed[-1]:
        trimmed.pop()
    while trimmed and not trimmed[0]:
        trimmed.pop(0)
    # Return a single string:
    return "\n".join(trimmed)


##################################
# End of code copied from PEP257 #
##################################


def dequote_docstring(text):
    """Remove the quotes delimiting a docstring."""
    # TODO: Process escaped characters unless raw mode?
    text = text.strip()
    if len(text) > 6 and text[:3] == text[-3:] == '"""':
        # Standard case, """..."""
        return text[3:-3]
    if len(text) > 7 and text[:4] in ('u"""', 'r"""') and text[-3:] == '"""':
        # Unicode, u"""...""", or raw r"""..."""
        return text[4:-3]
    # Other flake8 tools will report atypical quotes:
    if len(text) > 6 and text[:3] == text[-3:] == "'''":
        return text[3:-3]
    if len(text) > 7 and text[:4] in ("u'''", "r'''") and text[-3:] == "'''":
        return text[4:-3]
    if len(text) > 2 and text[0] == text[-1] == '"':
        return text[1:-1]
    if len(text) > 3 and text[:2] in ('u"', 'r"') and text[-1] == '"':
        return text[2:-1]
    if len(text) > 2 and text[0] == text[-1] == "'":
        return text[1:-1]
    if len(text) > 3 and text[:2] in ("u'", "r'") and text[-1] == "'":
        return text[2:-1]
    raise ValueError("Bad quotes!")


##################################################
# Start of code copied from pydocstyle/parser.py #
##################################################


def humanize(string):
    """Make a string human readable."""
    return re.compile(r"(.)([A-Z]+)").sub(r"\1 \2", string).lower()


class Value(object):
    """A generic object with a list of preset fields."""

    def __init__(self, *args):
        """Initialize."""
        if len(self._fields) != len(args):
            raise ValueError(
                "got {} arguments for {} fields for {}: {}".format(
                    len(args), len(self._fields), self.__class__.__name__, self._fields
                )
            )
        vars(self).update(zip(self._fields, args))

    def __hash__(self):
        """Hash."""
        return hash(repr(self))

    def __eq__(self, other):
        """Equality."""
        return other and vars(self) == vars(other)

    def __repr__(self):
        """Representation."""
        kwargs = ", ".join(
            "{}={!r}".format(field, getattr(self, field)) for field in self._fields
        )
        return "{}({})".format(self.__class__.__name__, kwargs)


class Definition(Value):
    """A Python source code definition (could be class, function, etc)."""

    _fields = (
        "name",
        "_source",
        "start",
        "end",
        "decorators",
        "docstring",
        "children",
        "parent",
        "skipped_error_codes",
    )

    _human = property(lambda self: humanize(type(self).__name__))
    kind = property(lambda self: self._human.split()[-1])
    module = property(lambda self: self.parent.module)
    all = property(lambda self: self.module.all)
    _slice = property(lambda self: slice(self.start - 1, self.end))
    is_class = False

    def __iter__(self):
        """Iterate."""
        return chain([self], *self.children)

    @property
    def _publicity(self):
        return {True: "public", False: "private"}[self.is_public]

    @property
    def source(self):
        """Return the source code for the definition."""
        full_src = self._source[self._slice]

        def is_empty_or_comment(line):
            return line.strip() == "" or line.strip().startswith("#")

        filtered_src = dropwhile(is_empty_or_comment, reversed(full_src))
        return "".join(reversed(list(filtered_src)))

    def __str__(self):
        """Definition as a string."""
        out = "in {} {} `{}`".format(self._publicity, self._human, self.name)
        if self.skipped_error_codes:
            out += " (skipping {})".format(self.skipped_error_codes)
        return out


class Module(Definition):
    """A Python source code module."""

    _fields = (
        "name",
        "_source",
        "start",
        "end",
        "decorators",
        "docstring",
        "children",
        "parent",
        "_all",
        "future_imports",
        "skipped_error_codes",
    )
    _nest = staticmethod(lambda s: {"def": Function, "class": Class}[s])
    module = property(lambda self: self)
    all = property(lambda self: self._all)

    @property
    def is_public(self):
        """Is the module public."""
        return not self.name.startswith("_") or self.name.startswith("__")

    def __str__(self):
        """Definition as a string."""
        return "at module level"


class Package(Module):
    """A package is a __init__.py module."""


class Function(Definition):
    """A Python source code function."""

    _nest = staticmethod(lambda s: {"def": NestedFunction, "class": NestedClass}[s])

    @property
    def is_public(self):
        """Return True iff this function should be considered public."""
        if self.all is not None:
            return self.name in self.all
        else:
            return not self.name.startswith("_")

    @property
    def is_test(self):
        """Return True if this function is a test function/method.

        We exclude tests from the imperative mood check, because to phrase
        their docstring in the imperative mood, they would have to start with
        a highly redundant "Test that ...".
        """
        return self.name.startswith("test") or self.name == "runTest"


class NestedFunction(Function):
    """A Python source code nested function."""

    is_public = False


class Method(Function):
    """A Python source code method."""

    @property
    def is_magic(self):
        """Return True iff this method is a magic method (e.g., `__str__`)."""
        return (
            self.name.startswith("__")
            and self.name.endswith("__")
            and self.name not in VARIADIC_MAGIC_METHODS
        )

    @property
    def is_public(self):
        """Return True iff this method should be considered public."""
        # Check if we are a setter/deleter method, and mark as private if so.
        for decorator in self.decorators:
            # Given 'foo', match 'foo.bar' but not 'foobar' or 'sfoo'
            if re.compile(r"^{}\.".format(self.name)).match(decorator.name):
                return False
        name_is_public = (
            not self.name.startswith("_")
            or self.name in VARIADIC_MAGIC_METHODS
            or self.is_magic
        )
        return self.parent.is_public and name_is_public


class Class(Definition):
    """A Python source code class."""

    _nest = staticmethod(lambda s: {"def": Method, "class": NestedClass}[s])
    is_public = Function.is_public
    is_class = True


class NestedClass(Class):
    """A Python source code nested class."""

    @property
    def is_public(self):
        """Return True iff this class should be considered public."""
        return (
            not self.name.startswith("_")
            and self.parent.is_class
            and self.parent.is_public
        )


class Decorator(Value):
    """A decorator for function, method or class."""

    _fields = "name arguments".split()


VARIADIC_MAGIC_METHODS = ("__init__", "__call__", "__new__")


class AllError(Exception):
    """Raised when there is a problem with __all__ when parsing."""

    def __init__(self, message):
        """Initialize the error with a more specific message."""
        Exception.__init__(
            self,
            message
            + textwrap.dedent(
                """
                That means pydocstyle cannot decide which definitions are
                public. Variable __all__ should be present at most once in
                each file, in form
                `__all__ = ('a_public_function', 'APublicClass', ...)`.
                More info on __all__: http://stackoverflow.com/q/44834/. ')
                """
            ),
        )


class TokenStream(object):
    """Token stream."""

    # A logical newline is where a new expression or statement begins. When
    # there is a physical new line, but not a logical one, for example:
    # (x +
    #  y)
    # The token will be tk.NL, not tk.NEWLINE.
    LOGICAL_NEWLINES = {tk.NEWLINE, tk.INDENT, tk.DEDENT}

    def __init__(self, filelike):
        """Initialize."""
        self._generator = tk.generate_tokens(filelike.readline)
        self.current = Token(*next(self._generator, None))
        self.line = self.current.start[0]
        self.log = log
        self.got_logical_newline = True

    def move(self):
        """Move."""
        previous = self.current
        current = self._next_from_generator()
        self.current = None if current is None else Token(*current)
        self.line = self.current.start[0] if self.current else self.line
        self.got_logical_newline = previous.kind in self.LOGICAL_NEWLINES
        return previous

    def _next_from_generator(self):
        try:
            return next(self._generator, None)
        except (SyntaxError, tk.TokenError):
            self.log.warning("error generating tokens", exc_info=True)
            return None

    def __iter__(self):
        """Iterate."""
        while True:
            if self.current is not None:
                yield self.current
            else:
                return
            self.move()


class TokenKind(int):
    """Kind of token."""

    def __repr__(self):
        """Representation."""
        return "tk.{}".format(tk.tok_name[self])


class Token(Value):
    """Token."""

    _fields = "kind value start end source".split()

    def __init__(self, *args):
        """Initialize."""
        super(Token, self).__init__(*args)
        self.kind = TokenKind(self.kind)


class Parser(object):
    """A Python source code parser."""

    def parse(self, filelike, filename):
        """Parse the given file-like object and return its Module object."""
        self.log = log
        self.source = filelike.readlines()
        src = "".join(self.source)
        # This may raise a SyntaxError:
        compile(src, filename, "exec")
        self.stream = TokenStream(StringIO(src))
        self.filename = filename
        self.all = None
        self.future_imports = set()
        self._accumulated_decorators = []
        return self.parse_module()

    # TODO: remove
    def __call__(self, *args, **kwargs):
        """Call the parse method."""
        return self.parse(*args, **kwargs)

    current = property(lambda self: self.stream.current)
    line = property(lambda self: self.stream.line)

    def consume(self, kind):
        """Consume one token and verify it is of the expected kind."""
        next_token = self.stream.move()
        assert next_token.kind == kind

    def leapfrog(self, kind, value=None):
        """Skip tokens in the stream until a certain token kind is reached.

        If `value` is specified, tokens whose values are different will also
        be skipped.
        """
        while self.current is not None:
            if self.current.kind == kind and (
                value is None or self.current.value == value
            ):
                self.consume(kind)
                return
            self.stream.move()

    def parse_docstring(self):
        """Parse a single docstring and return its value."""
        self.log.debug(
            "parsing docstring, token is %r (%s)", self.current.kind, self.current.value
        )
        while self.current.kind in (tk.COMMENT, tk.NEWLINE, tk.NL):
            self.stream.move()
            self.log.debug(
                "parsing docstring, token is %r (%s)",
                self.current.kind,
                self.current.value,
            )
        if self.current.kind == tk.STRING:
            docstring = self.current.value
            self.stream.move()
            return docstring
        return None

    def parse_decorators(self):  # noqa : D401
        """Called after first @ is found.

        Parse decorators into self._accumulated_decorators.
        Continue to do so until encountering the 'def' or 'class' start token.
        """
        name = []
        arguments = []
        at_arguments = False

        while self.current is not None:
            self.log.debug(
                "parsing decorators, current token is %r (%s)",
                self.current.kind,
                self.current.value,
            )
            if self.current.kind == tk.NAME and self.current.value in ["def", "class"]:
                # Done with decorators - found function or class proper
                break
            elif self.current.kind == tk.OP and self.current.value == "@":
                # New decorator found. Store the decorator accumulated so far:
                self._accumulated_decorators.append(
                    Decorator("".join(name), "".join(arguments))
                )
                # Now reset to begin accumulating the new decorator:
                name = []
                arguments = []
                at_arguments = False
            elif self.current.kind == tk.OP and self.current.value == "(":
                at_arguments = True
            elif self.current.kind == tk.OP and self.current.value == ")":
                # Ignore close parenthesis
                pass
            elif self.current.kind == tk.NEWLINE or self.current.kind == tk.NL:
                # Ignore newlines
                pass
            else:
                # Keep accumulating current decorator's name or argument.
                if not at_arguments:
                    name.append(self.current.value)
                else:
                    arguments.append(self.current.value)
            self.stream.move()

        # Add decorator accumulated so far
        self._accumulated_decorators.append(
            Decorator("".join(name), "".join(arguments))
        )

    def parse_definitions(self, class_, all=False):
        """Parse multiple definitions and yield them."""
        while self.current is not None:
            self.log.debug(
                "parsing definition list, current token is %r (%s)",
                self.current.kind,
                self.current.value,
            )
            self.log.debug("got_newline: %s", self.stream.got_logical_newline)
            if all and self.current.value == "__all__":
                self.parse_all()
            elif (
                self.current.kind == tk.OP
                and self.current.value == "@"
                and self.stream.got_logical_newline
            ):
                self.consume(tk.OP)
                self.parse_decorators()
            elif self.current.value in ["def", "class"]:
                yield self.parse_definition(class_._nest(self.current.value))
            elif self.current.kind == tk.INDENT:
                self.consume(tk.INDENT)
                for definition in self.parse_definitions(class_):
                    yield definition
            elif self.current.kind == tk.DEDENT:
                self.consume(tk.DEDENT)
                return
            elif self.current.value == "from":
                self.parse_from_import_statement()
            else:
                self.stream.move()

    def parse_all(self):
        """Parse the __all__ definition in a module."""
        assert self.current.value == "__all__"
        self.consume(tk.NAME)
        if self.current.value != "=":
            raise AllError("Could not evaluate contents of __all__. ")
        self.consume(tk.OP)
        if self.current.value not in "([":
            raise AllError("Could not evaluate contents of __all__. ")
        self.consume(tk.OP)

        self.all = []
        all_content = "("
        while self.current.kind != tk.OP or self.current.value not in ")]":
            if self.current.kind in (tk.NL, tk.COMMENT):
                pass
            elif self.current.kind == tk.STRING or self.current.value == ",":
                all_content += self.current.value
            else:
                raise AllError(
                    "Unexpected token kind in  __all__: {!r}. ".format(
                        self.current.kind
                    )
                )
            self.stream.move()
        self.consume(tk.OP)
        all_content += ")"
        try:
            self.all = eval(all_content, {})
        except BaseException as e:
            raise AllError(
                "Could not evaluate contents of __all__."
                "\bThe value was {}. The exception was:\n{}".format(all_content, e)
            )

    def parse_module(self):
        """Parse a module (and its children) and return a Module object."""
        self.log.debug("parsing module.")
        start = self.line
        docstring = self.parse_docstring()
        children = list(self.parse_definitions(Module, all=True))
        assert self.current is None, self.current
        end = self.line
        cls = Module
        if self.filename.endswith("__init__.py"):
            cls = Package
        module = cls(
            self.filename,
            self.source,
            start,
            end,
            [],
            docstring,
            children,
            None,
            self.all,
            None,
            "",
        )
        for child in module.children:
            child.parent = module
        module.future_imports = self.future_imports
        self.log.debug("finished parsing module.")
        return module

    def parse_definition(self, class_):
        """Parse a definition and return its value in a `class_` object."""
        start = self.line
        self.consume(tk.NAME)
        name = self.current.value
        self.log.debug("parsing %s '%s'", class_.__name__, name)
        self.stream.move()
        if self.current.kind == tk.OP and self.current.value == "(":
            parenthesis_level = 0
            while True:
                if self.current.kind == tk.OP:
                    if self.current.value == "(":
                        parenthesis_level += 1
                    elif self.current.value == ")":
                        parenthesis_level -= 1
                        if parenthesis_level == 0:
                            break
                self.stream.move()
        if self.current.kind != tk.OP or self.current.value != ":":
            self.leapfrog(tk.OP, value=":")
        else:
            self.consume(tk.OP)
        if self.current.kind in (tk.NEWLINE, tk.COMMENT):
            skipped_error_codes = self.parse_skip_comment()
            self.leapfrog(tk.INDENT)
            assert self.current.kind != tk.INDENT
            docstring = self.parse_docstring()
            decorators = self._accumulated_decorators
            self.log.debug("current accumulated decorators: %s", decorators)
            self._accumulated_decorators = []
            self.log.debug("parsing nested definitions.")
            children = list(self.parse_definitions(class_))
            self.log.debug("finished parsing nested definitions for '%s'", name)
            end = self.line - 1
        else:  # one-liner definition
            skipped_error_codes = ""
            docstring = self.parse_docstring()
            decorators = []  # TODO
            children = []
            end = self.line
            self.leapfrog(tk.NEWLINE)
        definition = class_(
            name,
            self.source,
            start,
            end,
            decorators,
            docstring,
            children,
            None,
            skipped_error_codes,
        )
        for child in definition.children:
            child.parent = definition
        self.log.debug(
            "finished parsing %s '%s'. Next token is %r (%s)",
            class_.__name__,
            name,
            self.current.kind,
            self.current.value,
        )
        return definition

    def parse_skip_comment(self):
        """Parse a definition comment for noqa skips."""
        skipped_error_codes = ""
        if self.current.kind == tk.COMMENT:
            if "noqa: " in self.current.value:
                skipped_error_codes = "".join(self.current.value.split("noqa: ")[1:])
            elif self.current.value.startswith("# noqa"):
                skipped_error_codes = "all"
        return skipped_error_codes

    def check_current(self, kind=None, value=None):
        """Verify the current token is of type `kind` and equals `value`."""
        msg = textwrap.dedent(
            """
        Unexpected token at line {self.line}:
        In file: {self.filename}
        Got kind {self.current.kind!r}
        Got value {self.current.value}
        """.format(
                self=self
            )
        )
        kind_valid = self.current.kind == kind if kind else True
        value_valid = self.current.value == value if value else True
        assert kind_valid and value_valid, msg

    def parse_from_import_statement(self):
        """Parse a 'from x import y' statement.

        The purpose is to find __future__ statements.
        """
        self.log.debug("parsing from/import statement.")
        is_future_import = self._parse_from_import_source()
        self._parse_from_import_names(is_future_import)

    def _parse_from_import_source(self):
        """Parse the 'from x import' part in a 'from x import y' statement.

        Return true iff `x` is __future__.
        """
        assert self.current.value == "from", self.current.value
        self.stream.move()
        is_future_import = self.current.value == "__future__"
        self.stream.move()
        while (
            self.current is not None
            and self.current.kind in (tk.DOT, tk.NAME, tk.OP)
            and self.current.value != "import"
        ):
            self.stream.move()
        if self.current is None or self.current.value != "import":
            return False
        self.check_current(value="import")
        assert self.current.value == "import", self.current.value
        self.stream.move()
        return is_future_import

    def _parse_from_import_names(self, is_future_import):
        """Parse the 'y' part in a 'from x import y' statement."""
        if self.current.value == "(":
            self.consume(tk.OP)
            expected_end_kinds = (tk.OP,)
        else:
            expected_end_kinds = (tk.NEWLINE, tk.ENDMARKER)
        while self.current.kind not in expected_end_kinds and not (
            self.current.kind == tk.OP and self.current.value == ";"
        ):
            if self.current.kind != tk.NAME:
                self.stream.move()
                continue
            self.log.debug(
                "parsing import, token is %r (%s)",
                self.current.kind,
                self.current.value,
            )
            if is_future_import:
                self.log.debug("found future import: %s", self.current.value)
                self.future_imports.add(self.current.value)
            self.consume(tk.NAME)
            self.log.debug(
                "parsing import, token is %r (%s)",
                self.current.kind,
                self.current.value,
            )
            if self.current.kind == tk.NAME and self.current.value == "as":
                self.consume(tk.NAME)  # as
                if self.current.kind == tk.NAME:
                    self.consume(tk.NAME)  # new name, irrelevant
            if self.current.value == ",":
                self.consume(tk.OP)
            self.log.debug(
                "parsing import, token is %r (%s)",
                self.current.kind,
                self.current.value,
            )


################################################
# End of code copied from pydocstyle/parser.py #
################################################


parse = Parser()


class reStructuredTextChecker(object):
    """Checker of Python docstrings as reStructuredText."""

    name = "rst-docstrings"
    version = __version__

    STDIN_NAMES = {"stdin", "-", "(none)", None}

    def __init__(self, tree, filename="(none)"):
        """Initialise."""
        self.tree = tree
        self.filename = filename
        try:
            self.load_source()
            self.err = None
        except Exception as err:
            self.source = None
            self.err = err

    @classmethod
    def add_options(cls, parser):
        """Add RST directives and roles options."""
        parser.add_option(
            "--rst-directives",
            metavar="LIST",
            default="",
            parse_from_config=True,
            comma_separated_list=True,
            help="Comma-separated list of additional RST directives.",
        )
        parser.add_option(
            "--rst-roles",
            metavar="LIST",
            default="",
            parse_from_config=True,
            comma_separated_list=True,
            help="Comma-separated list of additional RST roles.",
        )

    @classmethod
    def parse_options(cls, options):
        """Adding black-config option."""
        cls.extra_directives = options.rst_directives
        cls.extra_roles = options.rst_roles

    def run(self):
        """Use docutils to check docstrings are valid RST."""
        # Is there any reason not to call load_source here?
        if self.err is not None:
            assert self.source is None
            msg = "%s%03i %s" % (
                rst_prefix,
                rst_fail_load,
                "Failed to load file: %s" % self.err,
            )
            yield 0, 0, msg, type(self)
            module = []
        try:
            module = parse(StringIO(self.source), self.filename)
        except SyntaxError as err:
            msg = "%s%03i %s" % (
                rst_prefix,
                rst_fail_parse,
                "Failed to parse file: %s" % err,
            )
            yield 0, 0, msg, type(self)
            module = []
        except AllError:
            msg = "%s%03i %s" % (
                rst_prefix,
                rst_fail_all,
                "Failed to parse __all__ entry.",
            )
            yield 0, 0, msg, type(self)
            module = []
        for definition in module:
            if not definition.docstring:
                # People can use flake8-docstrings to report missing
                # docstrings
                continue
            try:
                # Note we use the PEP257 trim algorithm to remove the
                # leading whitespace from each line - this avoids false
                # positive severe error "Unexpected section title."
                unindented = trim(dequote_docstring(definition.docstring))
                # Off load RST validation to reStructuredText-lint
                # which calls docutils internally.
                # TODO: Should we pass the Python filename as filepath?
                rst_errors = list(rst_lint.lint(unindented))
            except Exception as err:
                # e.g. UnicodeDecodeError
                msg = "%s%03i %s" % (
                    rst_prefix,
                    rst_fail_lint,
                    "Failed to lint docstring: %s - %s" % (definition.name, err),
                )
                yield definition.start, 0, msg, type(self)
                continue
            for rst_error in rst_errors:
                # TODO - make this a configuration option?
                if rst_error.level <= 1:
                    continue
                # Levels:
                #
                # 0 - debug   --> we don't receive these
                # 1 - info    --> RST1## codes
                # 2 - warning --> RST2## codes
                # 3 - error   --> RST3## codes
                # 4 - severe  --> RST4## codes
                #
                # Map the string to a unique code:
                msg = rst_error.message.split("\n", 1)[0]
                code = code_mapping(
                    rst_error.level, msg, self.extra_directives, self.extra_roles
                )
                if not code:
                    # We ignored it, e.g. a known Sphinx role
                    continue
                assert 0 < code < 100, code
                code += 100 * rst_error.level
                msg = "%s%03i %s" % (rst_prefix, code, msg)

                # This will return the line number by combining the
                # start of the docstring with the offet within it.
                # We don't know the column number, leaving as zero.
                yield definition.start + rst_error.line, 0, msg, type(self)

    def load_source(self):
        """Load the source for the specified file."""
        if self.filename in self.STDIN_NAMES:
            self.filename = "stdin"
            if sys.version_info[0] < 3:
                self.source = sys.stdin.read()
            else:
                self.source = TextIOWrapper(sys.stdin.buffer, errors="ignore").read()
        else:
            # Could be a Python 2.7 StringIO with no context manager, sigh.
            # with tokenize_open(self.filename) as fd:
            #     self.source = fd.read()
            handle = tokenize_open(self.filename)
            self.source = handle.read()
            handle.close()
