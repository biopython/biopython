# Copyright 2000-2001, Dalke Scientific Software, LLC
# Distributed under the Biopython License Agreement (see the LICENSE file).

__version__ = "0.6"

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

def NoCase(expr):
    raise NotImplementedError

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


def Integer(name = None):
    """(name = None) -> match digits, possibly inside of a named group

    If 'name' is not None, the matching text will be put inside a group
    of the given name.
    
    """
    if name is None:
        return Re(r"\d+")
    else:
        return Group(name, Re(r"\d+"))

def SignedInteger(name = None):
    """(name = None) -> match digits with an optional leading minus sign

    If 'name' is not None, the matching text will be put inside of a
    group of the given name.
    """
    if name is None:
        return Re(r"-?\d+")
    else:
        return Group(name, Re("-?\d+"))

def Float(name = None):
    """(name = None) -> match floating point numbers like 6, 6., .1 and 2.3

    If 'name' is not None, the matching text will be put inside of a
    group of the given name.  Can be signed.
    """
    exp = Re(r"[+-]?((\d+(\.(\d+)?)?)|\.\d+)")
    if name is None:
        return exp
    else:
        return Group(name, exp)

def ToEol(name = None):
    """(name = None) -> match everything up to and including the end of line

    If 'name' is not None, the matching text, except for the newline, will
    be put inside a group of the given name.
    
    """
    if name is None:
        return Re(r"[^\R]*\R")
    else:
        return Group(name, Re(r"[^\R]*")) + AnyEol()

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

