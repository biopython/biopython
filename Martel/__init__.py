# Copyright 2000-2001, Dalke Scientific Software, LLC
# Distributed under the Biopython License Agreement (see the LICENSE file).

"""Martel is a 'regular expressions on steroids' parser generator (DEPRECATED).

A goal of the Biopython project is to reduce the amount of effort needed to
do computational biology. A large part of that work turns out to be parsing
file formats, which lead to the development of Martel, a parser generator
which uses a regular expression as the format description to create a parser
that returns the parse tree using the SAX API common in XML processing.

While intended to be both fast and relatively easy to understand, Martel did
struggle with some very large records (e.g. GenBank files for whole genomes or
chromosomes), and in practice debugging the Martel format specifications for
evolving file formats like GenBank proved non-trivial.

Andrew Dalke is no longer maintaining Martel or Bio.Mindy, and these modules
are now deprecated.  They are no longer used in any of the current Biopython
parsers, and are likely to be removed in a future release of Biopython.
"""
__version__ = "1.50"

import warnings
warnings.warn("Martel and those parts of Biopython depending on it" \
              +" directly, such as Bio.Mindy, are now deprecated, and will" \
              +" be removed in a future release of Biopython.  If you want"\
              +" to continue to use this code, please get in contact with"\
              +" the Biopython developers via the mailing lists to avoid"\
              +" its permanent removal from Biopython", \
              DeprecationWarning)

import Expression
import convert_re
import string
from xml.sax import xmlreader

# This interface is based off of Greg Ewing's Plex parser.

def Str1(s):
    """(s) -> match the literal string"""
    return Expression.Str(s)

def Str(*args):
    """(s1, s2, ...) -> match s1 or s2 or ..."""
    if len(args) == 1:
        return Str1(args[0])
    return Expression.Alt(tuple(map(Str, args)))

def Any(s):
    """(s) -> match any character in s"""
    if len(s) == 1:
        return Expression.Literal(s)
    return Expression.Any(s, 0)

def AnyBut(s):
    """s -> match any character not in s"""
    return Expression.Any(s, 1)

## Untested!
##def AnyChar():
##    """match any character, including newline"""
##    return Expression.Re("(?:.|\n)")


def Seq(*args):
    """exp1, exp2, ... -> match exp1 followed by exp2 followed by ..."""
    # I'm always forgetting and passing strings into this function,
    # so make sure the arguments are expressions
    for arg in args:
        assert isinstance(arg, Expression.Expression), \
               "expecting an Expression, not a %s" % type(arg)
    
    return Expression.Seq(args)

def Alt(*args):
    """exp1, exp2, ... -> match exp1 or (if that fails) match exp2 or ..."""
    # Do some type checking
    for arg in args:
        assert isinstance(arg, Expression.Expression), \
               "expecting an Expression, not a %s" % type(arg)
    return Expression.Alt(args)

def Opt(expr):
    """expr -> match 'expr' 1 or 0 times"""
    assert isinstance(expr, Expression.Expression), \
           "expecting an Expression, not a %s" % type(expr)
    return Expression.MaxRepeat(expr, 0, 1)

def Rep(expr):
    """expr -> match 'expr' as many times as possible, even 0 time"""
    assert isinstance(expr, Expression.Expression), \
           "expecting an Expression, not a %s" % type(expr)
    return Expression.MaxRepeat(expr, 0)

def Rep1(expr):
    """expr -> match 'expr' as many times as possible, but at least once"""
    assert isinstance(expr, Expression.Expression), \
           "expecting an Expression, not a %s" % type(expr)
    return Expression.MaxRepeat(expr, 1)


# These are in Plex, but I don't (yet?) implement them

NoCase = Expression.NoCase

def Case(expr):
    raise NotImplementedError

def Bol():
    raise NotImplementedError

def Eol():
    raise NotImplementedError

def Empty():
    raise NotImplementedError

def Eof():
    raise NotImplementedError
    return Expression.AtEnd()


# Not in Plex, but useful

AnyEol = Expression.AnyEol

def MaxRepeat(expr, min_count, max_count = Expression.MAXREPEAT):
    """expr, min_count, max_count = 65535 -> match between min- and max_count times

    If max_count == 65535 (which is Expression.MAXREPEAT) then there
    is no upper limit.
    """
    assert isinstance(expr, Expression.Expression), \
           "expecting an Expression, not a %s" % type(expr)
    return Expression.MaxRepeat(expr, min_count, max_count)
    
def RepN(expr, count):
    """expr, count -> match the expression 'count' number of time

    This option is handy for named group repeats since you don't have
    to use the name twice; for the min_count and max_count fields.
    """
    return Expression.MaxRepeat(expr, count, count)


def Group(name, expr, attrs = None):
    """name, expr -> use 'name' to describe a successful match of the expression"""
    assert isinstance(expr, Expression.Expression), \
           "expecting an Expression, not a %s" % type(expr)
    return Expression.Group(name, expr, attrs)


def _fix_newlines(s):
    # Replace the characters "\r\n", "\r" and "\n" with \R.
    # Does not affect the substrings '\' + 'n' or '\' + 'n'
    s = string.replace(s, "\r\n", "\n")
    s = string.replace(s, "\r", "\n")
    return string.replace(s, "\n", r"\R")
    

def Re(pattern, fix_newlines = 0):
    """pattern -> the expression tree for the regexp pattern string"""
    if fix_newlines:
        pattern = _fix_newlines(pattern)
    return convert_re.make_expression(pattern)

NullOp = Expression.NullOp
Debug = Expression.Debug

def Assert(expression):
    return Expression.Assert(expression)

def AssertNot(expression):
    return Expression.Assert(expression, invert = 1)

# helper function
def _group(name, exp, attrs):
    if name is None:
        assert not attrs, "Attributes (%s) require a group name" % (attrs,)
        return exp
    return Group(name, exp, attrs)

def Digits(name = None, attrs = None):
    """match one or more decimal digits

    This is the same as (?P<name?attrs>\d+).
    
    If 'name' is not None, the matching text will be put inside a
    group of the given name.  You can optionally include group
    attributes.    
    """
    return _group(name, Re(r"\d+"), attrs)

def Integer(name = None, attrs = None):
    """match an integer (digits w/ optional leading + or - sign)

    If 'name' is not None, the matching text will be put inside a
    group of the given name.  You can optionally include group
    attributes.    
    """
    exp = Re(r"[+-]?\d+")
    return _group(name, exp, attrs)

def Float(name = None, attrs = None):
    """match floating point numbers like 6, 6., -.1, 2.3, +4E-5, ...

    If 'name' is not None, the matching text will be put inside of a
    group of the given name.  You can optionally include group
    attributes.
    """
    exp = Re(r"[+-]?((\d+(\.\d*)?)|\.\d+)([eE][+-]?[0-9]+)?")
    return _group(name, exp, attrs)

def Word(name = None, attrs = None):
    """match a 'word'
 
    A 'word' is defined as '\w+', and \w is [a-zA-Z0-9_].
 
    If 'name' is not None, the matching text will be put inside of a
    group of the given name.  You can optionally include group
    attributes.    
 
    In other words, this is the short way to write (?P<name>\w+).
    """
    exp = Re(r"\w+")
    return _group(name, exp, attrs)
 
def Spaces(name = None, attrs = None):
    """match one or more whitespace (except newline)
 
    "Spaces" is defined as [\\t\\v\\f\\r ]+, which is *not* the same
    as '\\s+'.  (It's missing the '\\n', which is useful since you
    almost never mean for whitespace to go beyond the newline.)
 
    If 'name' is not None, the matching text will be put inside of a
    group of the given name.  You can optionally include group
    attributes.
    """
    exp = Re(r"[\t\v\f ]+")
    return _group(name, exp, attrs)

def Unprintable(name = None, attrs = None):
    """match an unprintable character (characters not in string.printable)

    If 'name' is not None, the matching text will be put inside of a
    group of the given name.  You can optionally include group
    attributes.    
    """
    return _group(name, AnyBut(string.printable), attrs)

def Punctuation(name = None, attrs = None):
    """match a punctuation character (characters in string.punctuation)

    If 'name' is not None, the matching text will be put inside of a
    group of the given name.  You can optionally include group
    attributes.    
    """
    return _group(name, Any(string.punctuation), attrs)


def ToEol(name = None, attrs = None):
    """match everything up to and including the end of line

    If 'name' is not None, the matching text, except for the newline,
    will be put inside a group of the given name.  You can optionally
    include group attributes.    
    """
    if name is None:
        assert not attrs, "Attributes (%s) require a group name" % (attrs,)
        return Re(r"[^\R]*\R")
    else:
        return Group(name, Re(r"[^\R]*"), attrs) + AnyEol()

def UntilEol(name = None, attrs = None):
    """match everything up to but not including the end of line

    If 'name' is not None, the matching text, except for the newline,
    will be put inside a group of the given name.  You can optionally
    include group attributes.    
    """
    if name is None:
        assert not attrs, "Attributes (%s) require a group name" % (attrs,)
        return Re(r"[^\R]*")
    else:
        return Group(name, Re(r"[^\R]*"), attrs)

def SkipLinesUntil(expr):
    """read and ignore lines up to, but excluding, the line matching expr"""
    return Rep(AssertNot(expr) + ToEol())

def SkipLinesTo(expr):
    """read and ignore lines up to and including, the line matching expr"""
    return Rep(AssertNot(expr) + ToEol()) + expr + ToEol()


def ToSep(name = None, sep = None, attrs = None):
    """match all characters up to the given seperator(s)

    This is useful for parsing space, tab, color, or other character
    delimited fields.  There is no default seperator character.
    
    If 'name' is not None, the matching text, except for the seperator
    will be put inside a group of the given name.  You can optionally
    include group attributes.  The seperator character will also be
    consumed.

    Neither "\\r" nor "\\n" may be used as a seperator
    """
    if sep is None:
        # I found it was too easy to make a mistake with a default
        raise TypeError("Must specify a seperator (the 'sep' parameter)")

    assert "\r" not in sep and "\n" not in sep, \
           "cannot use %s as a seperator" % (repr(seperator),)

    exp = Rep(AnyBut(sep + "\r\n"))
    return _group(name, exp, attrs) + Str(sep)

def UntilSep(name = None, sep = None, attrs = None):
    """match all characters up to the given seperators(s)

    This is useful for parsing space, tab, color, or other character
    delimited fields.  There is no default seperator.
    
    If 'name' is not None, the matching text, except for the seperator
    will be put inside a group of the given name.  You can optionally
    include group attributes.  The seperator character will not be
    consumed.

    Neither "\\r" nor "\\n" may be used as a seperator.
    """
    if sep is None:
        # I found it was too easy to make a mistake with a default
        raise TypeError("Must specify a seperator (the 'sep' parameter)")

    assert "\r" not in sep and "\n" not in sep, \
           "cannot use %s as a seperator" % (repr(sep),)

    exp = Rep(AnyBut(sep + "\r\n"))
    return _group(name, exp, attrs)


def DelimitedFields(name = None, sep = None, attrs = None):
    """match 0 or more fields seperated by the given seperator(s)

    This is useful for parsing space, tab, color, or other character
    delimited fields.  There is no default seperator.

    If 'name' is not None, the delimited text, excluding the seperator,
    will be put inside groups of the given name.  You can optionally
    include group attributes.  The seperator character is consumed,
    but not accessible using a group.

    Neither "\\r" nor "\\n" may be used as a seperator.
    The line as a whole is not included in a group.
    """
    if sep is None:
        # I found it was too easy to make a mistake with a default
        raise TypeError("Must specify a sep (via the 'sep' parameter)")

    assert "\r" not in sep and "\n" not in sep, \
           "cannot use %s as a seperator" % (repr(sep),)

    term = _group(name, Rep(AnyBut(sep + "\r\n")), attrs)
    rep = Rep(Any(sep) + term)
    return term + rep + AnyEol()
    
# Allows some optimizations
FastFeature = Expression.FastFeature


# Used when making parsers which read a record at a time
ParseRecords = Expression.ParseRecords
HeaderFooter = Expression.HeaderFooter

# Use this to prune out group names you aren't
# interested in seeing, which reduces the number of method
# calls back to the parser.
def select_names(expression, names):
    # Make a copy so I know I don't share subexpressions which other
    # expressions.
    exp = expression.copy()

    # Use that internal method I told you not to use :)
    exp._select_names(names)

    # Get rid of unnamed groups
    import optimize
    return optimize.optimize_unnamed_groups(exp)

def replace_groups(expr, replacements):
    expr = expr.copy()
    for tagname, replacement_expr in replacements:
        matches = expr._find_groups(tagname)
        for match in matches:
            match.expression = replacement_expr
    return expr

def SimpleRecordFilter(expr, make_reader, reader_args = ()):
    return ParseRecords("dataset", {"format": "*filter*"},
                        Group("record", expr + Rep(ToEol())),
                        make_reader, reader_args)
