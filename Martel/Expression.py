# Copyright 2000-2001, Dalke Scientific Software, LLC
# Distributed under the Biopython License Agreement (see the LICENSE file).

"""Classes for nodes in the Expression tree.

Expression
 |--- Any           - match (or don't match) a set of characters
 |--- AnyEol        - match any newline representation ("\n", "\r" or "\r\n")
 |--- Assert        - used for positive and negative lookahead assertions 
 |--- AtBeginning   - match the beginning of a line
 |--- AtEnd         - match the end of a line
 |--- Debug         - print a debug message
 |--- Dot           - match any character except newline
 |--- Group         - give a group name to an expression
 |--- GroupRef      - match a previously identified expression
 |--- Literal       - match (or don't match) a single character
 |--- MaxRepeat     - greedy repeat of an expression, within min/max bounds
 |--- NullOp        - does nothing (useful as an initial seed)
 |--- PassThrough   - used when overriding 'make_parser'; match its subexp
 |      |--- FastFeature  - keeps information about possibly optional tags
 |      |--- HeaderFooter - files with a header, records and a footer
 |      `--- ParseRecords - parse a record at a time
 |--- Str           - match a given string
 `--- ExpressionList  - expressions containing several subexpressions
        |--- Alt    - subexp1 or subexp2 or subexp3 or ...
        `--- Seq    - subexp1 followed by subexp2 followed by subexp3 ...
"""
import re, string
from xml.sax import xmlreader
import msre_parse  # Modified version of Secret Labs' sre_parse
import Parser

MAXREPEAT = msre_parse.MAXREPEAT

try:
    import IterParser
except SyntaxError:
    IterParser = None

class Expression:
    """Base class for nodes in the Expression tree"""
    def __add__(self, other):
        """returns an Expression to match this Expression then the other one"""
        assert isinstance(other, Expression), "RHS of '+' not an Expression"
        return Seq( (self, other) )
    
    def __or__(self, other):
        """returns an Expression matching this Expression or (if that fails) the other one"""
        assert isinstance(other, Expression), "RHS of '|' not an Expression"
        return Alt( (self, other) )
    
    def group_names(self):
        """the list of group names used by this Expression and its children"""
        return ()

    def _find_groups(self, tag):
        """return a list of all groups matching the given tag"""
        return []

    def features(self):
        """return a list of all features"""
        return []
    
    def _select_names(self, names):
        """internal function used by 'select_names'.

        Don't call this function.  Will likely be removed in future versions.
        """
        # subtrees can be shared so you need to copy first before selecting
        pass
    
    def copy(self):
        """do a deep copy on this Expression tree"""
        raise NotImplementedError
    
    def __str__(self):
        """the corresponding pattern string"""
        raise NotImplementedError
    
    def make_parser(self, debug_level = 0):
        """create a SAX compliant parser for this regexp"""
        import Generate
        tagtable, want_flg, attrlookup = Generate.generate(self, debug_level)
        return Parser.Parser(tagtable, (want_flg, debug_level, attrlookup))

    def make_iterator(self, tag = "record", debug_level = 0):
        """create an iterator for this regexp; the 'tag' defines a record"""
        import Iterator
        return Iterator.Iterator(self.make_parser(debug_level), tag)

    def _modify_leaves(self, func):
        """internal function for manipulating the leaves of an expression

        This really needs to be some sort of visit pattern, but I'm
        not sure the best way to do it.  THIS METHOD MAY CHANGE.
        """
        return func(self)
    
# Any character in a given set: '[abc]'
class Any(Expression):
    def __init__(self, chars, invert = 0):
        """(chars, invert = 0)

        Match any of the characters appearing in the 'chars' string.
        If 'invert' is true, match a character not in the string.
        """
        self.chars = chars
        self.invert = invert

    def copy(self):
        """do a deep copy on this Expression tree"""
        return Any(self.chars, self.invert)

    def __str__(self):
        """the corresponding pattern string"""
        if self.invert:
            return '[^%s]' % _minimize_any_range(self.chars)
        else:
            return '[%s]' % _minimize_any_range(self.chars)

# Lookahead assertions: '(?=...)'
#                       '(?!...)'
class Assert(Expression):
    def __init__(self, expression, invert = 0):
        """(expression, invert = 0)

        A non-consuming assertion using the given expression.
        The default is a positive lookahead, which matches if the expression
        matches at the current position, but does not affect the character
        position.

        If 'invert' is false, this is a negative lookahead assertion,
        and matches if the expression does not match.  Again, the character
        position is not affected.
        """
        self.expression = expression
        self.invert = invert
        
    def copy(self):
        """do a deep copy on this Expression tree"""
        return Assert(self.expression.copy(), self.invert)
    
    def __str__(self):
        """the corresponding pattern string"""
        if self.invert:
            return '(?!%s)' % str(self.expression)
        else:
            return '(?=%s)' % str(self.expression)

    def _modify_leaves(self, func):
        exp = self.expression._modify_leaves(func)
        assert exp is not None
        self.expression = exp
        return self

# At the beginning of the string: '^' in multiline mode
#
# There should be no reason to use this object because the Martel
# grammers all use explicit newlines anyway.

class AtBeginning(Expression):
    """Match the beginning of a line"""
    def copy(self):
        """do a deep copy on this Expression tree"""
        return AtBeginning()
    def __str__(self):
        """the corresponding pattern string"""
        return '^'

# At the end of the string: '$' in multiline mode
#
# There should be no reason to use this object because the Martel
# grammers all use explicit newlines anyway.
class AtEnd(Expression):
    """Match the end of a line"""
    def copy(self):
        """do a deep copy on this Expression tree"""
        return AtEnd()
    def __str__(self):
        """the corresponding pattern string"""
        return '$'

# Print a message when there is a match at this point.
# Helpful for debugging
class Debug(Expression):
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        # There is no pattern for this
        return ""
    def copy(self):
        """do a deep copy on this Expression"""
        return Debug(self.msg)
    
# Any character except newline: '.'
class Dot(Expression):
    """Match any character except newline"""
    def copy(self):
        """do a deep copy on this Expression tree"""
        return Dot()
    def __str__(self):
        """the corresponding pattern string"""
        return '.'

# Read any one of the newline conventions
class AnyEol(Expression):
    """Match a newline ("\n", "\r" or "\r\n")"""
    def copy(self):
        """do a deep copy on this Expression tree"""
        return AnyEol()
    def __str__(self):
        """the corresponding pattern string"""
        return r"(\n|\r\n?)"

# A grouping: '(?P<name>expression)'
#             '(?:expression)'
#             '(expression)'    -- same as (?:...) since I don't track
#                                  or use the group number

# Group names must be valid XML identifiers
def _verify_name(s):
    assert s, "Group name can not be the empty string"
    if not msre_parse.isname(s):
        raise AssertionError, "Illegal character in group name %s" % repr(s)

_fast_quote_lookup = None
def _make_fast_lookup():
    global _fast_quote_lookup
    
    safe = ('ABCDEFGHIJKLMNOPQRSTUVWXYZ'
            'abcdefghijklmnopqrstuvwxyz'
            '0123456789' '_.-')
    lookup = {}
    for c in range(256):
        lookup[chr(c)] = '%%%02X' % c
    for c in safe:
        lookup[c] = c
    
    _fast_quote_lookup = lookup
    
def _quote(s):
    if _fast_quote_lookup is None:
        _make_fast_lookup()
    lookup = _fast_quote_lookup
    terms = []
    if s:
        for c in s:
            terms.append(lookup[c])
    return string.join(terms, "")

def _make_group_pattern(name, expression, attrs):
    if name is None:
        return '(%s)' % str(expression)

    elif attrs:
        # Convert them to the proper URL-encoded form
        terms = []
        for k, v in attrs.items():
            terms.append("%s=%s" % (_quote(k), _quote(v)))
        attrname = name + "?" + string.join(terms, "&")
        return '(?P<%s>%s)' % (attrname, str(expression))
    else:
        return '(?P<%s>%s)' % (name, str(expression))

class Group(Expression):
    def __init__(self, name, expression, attrs = None):
        """(name, expression)

        Create a group named 'name' which matches the given expression
        """
        if name is not None:
            _verify_name(name)
        self.name = name
        self.expression = expression
        if attrs is None:
            attrs = xmlreader.AttributesImpl({})
        elif isinstance(attrs, type({})):
            attrs = xmlreader.AttributesImpl(attrs)
        self.attrs = attrs
            
    def group_names(self):
        """the list of group names used by this Expression and its children"""
        subnames = self.expression.group_names()
        if self.name is not None:
            if self.name in subnames:
                return subnames
            else:
                return (self.name, ) + subnames
        return subnames

    def _find_groups(self, tag):
        """return a list of all groups matching the given tag"""
        x = []
        if self.name == tag:
            x.append(self)
        return x + self.expression._find_groups(tag)

    def features(self):
        """return a list of all features"""
        return self.expression.features()

    def _select_names(self, names):
        """internal function: do not use"""
        if self.name is not None and self.name not in names:
            self.name = None
            self.attrs = xmlreader.AttributesImpl({})
        self.expression._select_names(names)

    def _modify_leaves(self, func):
        exp = self.expression._modify_leaves(func)
        assert exp is not None
        self.expression = exp
        return self

    def copy(self):
        """do a deep copy on this Expression tree"""
        return Group(self.name, self.expression.copy(), self.attrs.copy())

    def __str__(self):
        """the corresponding pattern string"""
        return _make_group_pattern(self.name, self.expression, self.attrs)


# group reference: '(?P<name>.)(?P=name)'
class GroupRef(Expression):
    def __init__(self, name):
        """(name)

        Match the same text previously found by the given named group
        """
        _verify_name(name)
        self.name = name
    def copy(self):
        """do a deep copy on this Expression tree"""
        return GroupRef(self.name)
    def __str__(self):
        """the corresponding pattern string"""
        return "(?P=%s)" % self.name



# A single character: 'a'
# This exists to simplify the sre conversion, but is not directly exposed
# as part of Martel's public expression API
class Literal(Expression):
    def __init__(self, char, invert = 0):
        """(char, invert = 0)

        Match the given character or, if 'invert' is true, match a character
        which is not this character.
        """
        self.char = char
        self.invert = invert

    def copy(self):
        """do a deep copy on this Expression tree"""
        return Literal(self.char, self.invert)

    def __str__(self):
        """the corresponding pattern string"""
        c = escape(self.char)
        if self.invert:
            return '[^%s]' % c
        return c

        
# Greedy repeat: 'a*'
#                'a{3,5}'
#                'a+'
class MaxRepeat(Expression):
    def __init__(self, expression, min_count=0,
                 max_count=MAXREPEAT):
        """(expression, min_count = 0, max_count = MAXREPEAT)

        Match the expression at least 'min_count' times and no more
        than 'max_count' times.  If max_count == MAXREPEAT then
        there is no fixed upper limit.

        min_count and max_count can be strings, in which case they are
        used as "named group repeats."  That is, they are taken to be
        group names and used to find the repeat counts during
        evaluation time.  The current implementation only understands
        named group repeats when min_count == max_count.

        The grouping is greedy.

        WARNING: There is no check to ensure that a match of 0 size is
        repeated indefinitely, as with "(a?)*" against the string "b".
        This will loop forever.

        WARNING: The current implementation does not support
        backtracking in MaxRepeats, so ".*\n" will not match "\n".
        Use a more explicit construct instead, like "[^\n]*\n".

        """
        
        self.expression = expression

        # Do some range checking
        if type(min_count) == type(0):
            assert 0 <= min_count, \
                   "min_count must be non-negative, not %d" % min_count
            if type(max_count) == type(0):
                assert min_count <= max_count, \
                       "min_count (%d) must be <= max_count (%d)" % (min_count, max_count)
        else:
            if type(max_count) == type(0):
                assert max_count >= 0, \
                       "max_count must be non-negative, not %d" % max_count

        self.min_count = min_count
        self.max_count = max_count
        
    def group_names(self):
        """the list of group names used by this Expression and its children"""
        # These are the names created by this Expression, not the names
        # *used* by it.
        return self.expression.group_names()
    def _find_groups(self, tag):
        return self.expression._find_groups(tag)

    def features(self):
        """return a list of all features"""
        return self.expression.features()
    
    def copy(self):
        """do a deep copy on this Expression tree"""
        return MaxRepeat(self.expression.copy(), self.min_count,
                         self.max_count)
    
    def _select_names(self, names):
        """internal function: do not use"""
        self.expression._select_names(names)
        
    def _modify_leaves(self, func):
        exp = self.expression._modify_leaves(func)
        assert exp is not None
        self.expression = exp
        return self
    
    def __str__(self):
        """the corresponding pattern string"""
        min_count = self.min_count
        max_count = self.max_count
        subexp = self.expression

        # If the subexpression is an Alt or Seq, then I need to put
        # them inside their own group, since "(a|b)*" is not the same
        # as "a|b*" and "ab*" is not the same as "(ab)*".
        if isinstance(subexp, ExpressionList):
            need_group = 1
        elif isinstance(subexp, Str):
            # Strings are also special, since it's a special case
            # of Seq( (Literal(s[0]), Literal(s[1]), ... ) )
            need_group = 1
        else:
            need_group = 0

        if need_group:
            s = "(%s)" % str(subexp)
        else:
            s = str(subexp)

        # Find the "extension" put at the end of the expression string
        if type(min_count) == type("") or type(max_count) == type(""):
            # Handle named group repeats
            if min_count == max_count:
                ext = "{%s}" % min_count
            else:
                ext = "{%s,%s}" % (min_count, max_count)
        else:
            # Make things pretty by using the special regexp pattern notation
            if min_count == 0 and max_count == MAXREPEAT:
                ext = "*"
            elif min_count == 1 and max_count == MAXREPEAT:
                ext = "+"
            elif min_count == 0 and max_count == 1:
                ext = "?"
            elif min_count == max_count == 1:
                ext = ""
            elif min_count == max_count:
                ext = "{%d}" % max_count
            elif min_count == 0:
                ext = "{,%d}" % max_count
            elif max_count == MAXREPEAT:
                ext = "{%d,}" % min_count
            else:
                ext = "{%d,%d}" % (min_count, max_count)

        return s + ext

# does nothing
class NullOp(Expression):
    def __init__(self):
        """()

        Doesn't match anything.  This is a null operation.  It's
        useful if you want a valid initial object from which to build,
        as in:

          exp = NullOp()
          for c in string.split(line):
            exp = exp + Str(c)

        (That's contrived -- see Time.py for a real use.)
        """
    def _select_names(self, names):
        pass
    def copy(self):
        return NullOp()
    def __str__(self):
        return ""
    def __add__(self, other):
        return other
    def __or__(self, other):
        raise TypeError("Cannot 'or' a NullOp with anything (only 'and')")


# Match the subexpression.
class PassThrough(Expression):
    def __init__(self, expression):
        """(expression)

        Match the given subexpression.  This class should not be used
        directly.  It is meant for generating specialized parsers which
        read a record at a time.
        """
        self.expression = expression
    def _select_names(self, names):
        self.expression._select_names(names)
    def _modify_leaves(self, func):
        exp = self.expression._modify_leaves(func)
        assert exp is not None
        self.expression = exp
        return self
    def copy(self):
        """do a deep copy on this Expression tree"""
        return PassThrough(self.expression.copy())
    def __str__(self):
        """the corresponding pattern string"""
        return str(self.expression)
    def group_names(self):
        return self.expression.group_names()
    def _find_groups(self, tag):
        return self.expression._find_groups(tag)
    def features(self):
        """return a list of all features"""
        return self.expression.features()

class FastFeature(PassThrough):
    def __init__(self, expression, feature, remove_tags):
        PassThrough.__init__(self, expression)
        self.feature = feature
        self.remove_tags = remove_tags
    def copy(self):
        """do a deep copy on this Expression tree"""
        return FastFeature(self.expression.copy(), self.feature,
                           self.remove_tags[:])
    def features(self):
        return [(self.feature, self.remove_tags)]

class HeaderFooter(PassThrough):
    def __init__(self, format_name, attrs,
                 header_expression, make_header_reader, header_args,
                 record_expression, make_record_reader, record_args,
                 footer_expression, make_footer_reader, footer_args):
        # I added attrs to the parameter list but couldn't make it
        # backwards compatible.  Without this check, it's possible to
        # have the object constructed seemingly okay then have the
        # error appear downstream, making it hard to track down.
        if isinstance(attrs, Expression):
            raise TypeError("Looks like you need an attrs between the format_name and the record_expression")

        
        if header_expression is None:
            assert make_header_reader is None and header_args is None
            exp = MaxRepeat(record_expression, 1)
        else:
            exp = header_expression + MaxRepeat(record_expression, 1)
        if footer_expression is not None:
            exp = exp + footer_expression
        else:
            assert make_footer_reader is None and footer_args is None
        PassThrough.__init__(self, Group(format_name, exp, attrs))

        self.format_name = format_name
        if attrs is None:
            attrs = xmlreader.AttributesImpl({})
        elif isinstance(attrs, type({})):
            attrs = xmlreader.AttributesImpl(attrs)
        self.attrs = attrs
        self.header_expression = header_expression
        self.make_header_reader = make_header_reader
        self.header_args = header_args
        self.record_expression = record_expression
        self.make_record_reader = make_record_reader
        self.record_args = record_args
        self.footer_expression = footer_expression
        self.make_footer_reader = make_footer_reader
        self.footer_args = footer_args

    def copy(self):
        header_exp = self.header_expression
        if header_exp is not None: header_exp = header_exp.copy()

        record_exp = self.record_expression
        if record_exp is not None: record_exp = record_exp.copy()

        footer_exp = self.footer_expression
        if footer_exp is not None: footer_exp = footer_exp.copy()
            
        return HeaderFooter(
            self.format_name, self.attrs.copy(),
            header_exp, self.make_header_reader, self.header_args,
            record_exp, self.make_record_reader, self.record_args,
            footer_exp, self.make_footer_reader, self.footer_args)

    def _modify_leaves(self, func):
        header_exp = self.header_expression
        if header_exp is not None:
            header_exp = header_exp.modify_leaves(func)
            assert header_exp is not None
            self.header_expression = header_exp
        record_exp = self.record_expression
        if record_exp is not None:
            record_exp = record_exp.modify_leaves(func)
            assert record_exp is not None
            self.record_expression = record_exp
        footer_exp = self.footer_expression
        if footer_exp is not None:
            footer_exp = footer_exp.modify_leaves(func)
            assert footer_exp is not None
            self.footer_expression = footer_exp
        return self
            
    def make_parser(self, debug_level = 0):
        import Generate, RecordReader
        want = 0
        if self.header_expression is not None:
            header_tagtable, want_flg, attrlookup = \
                             Generate.generate(self.header_expression,
                                               debug_level = debug_level)
            make_header_reader = self.make_header_reader
            header_args = self.header_args
        else:
            header_tagtable = ()
            want_flg = 0
            attrlookup = {}
            make_header_reader = None,
            header_args = None
            

        record_tagtable, want_flag, tmp_attrlookup = \
                         Generate.generate(self.record_expression,
                                           debug_level = debug_level)
        make_record_reader = self.make_record_reader
        record_args = self.record_args
        attrlookup.update(tmp_attrlookup)
        
        want = want or want_flg

        if self.footer_expression is not None:
            footer_tagtable, want_flag, tmp_attrlookup = \
                             Generate.generate(self.footer_expression,
                                               debug_level = debug_level)
            make_footer_reader = self.make_footer_reader
            footer_args = self.footer_args
            attrlookup.update(tmp_attrlookup)
        else:
            footer_tagtable = ()
            want_flg = 0
            make_footer_reader = None
            footer_args = None
        
        want = want or want_flg

        return Parser.HeaderFooterParser(
            self.format_name, self.attrs,
            make_header_reader, header_args, header_tagtable,
            make_record_reader, record_args, record_tagtable,
            make_footer_reader, footer_args, footer_tagtable,
            (want, debug_level, attrlookup))

    def make_iterator(self, tag, debug_level = 0):
        """create an iterator for this regexp; the 'tag' defines a record"""
        import Iterator
        if tag == self.format_name:
            return self.expression.make_iterator(self, tag)

        if self.header_expression is None:
            header_parser = None
        else:
            header_parser = self.header_expression.make_parser(debug_level)

        assert self.record_expression is not None
        record_parser = self.record_expression.make_parser(debug_level)

        if self.footer_expression is None:
            footer_parser = None
        else:
            footer_parser = self.footer_expression.make_parser(debug_level)

        if isinstance(self.record_expression, Group) and \
           self.record_expression.name == tag and \
           IterParser is not None:
            # There's an optimization for this case
            return IterParser.IterHeaderFooter(
                header_parser, self.make_header_reader, self.header_args,
                record_parser, self.make_record_reader, self.record_args,
                footer_parser, self.make_footer_reader, self.footer_args,
                tag
                )
        
        return Iterator.IteratorHeaderFooter(
            header_parser, self.make_header_reader, self.header_args,
            record_parser, self.make_record_reader, self.record_args,
            footer_parser, self.make_footer_reader, self.footer_args,
            tag
            )

    def group_names(self):
        x = [self.format_name]
        if self.header_expression is not None:
            x.extend(self.header_expression.group_names())
        x.extend(self.expression.group_names())
        if self.footer_expression is not None:
            x.extend(self.footer_expression.group_names())
        return x
    def _find_groups(self, tag):
        assert tag != self.format_name, "can't handle that case"
        x = []
        if self.header_expression is not None:
            x.extend(self.header_expression._find_groups(tag))
        x.extend(self.expression._find_groups(tag))
        if self.footer_expression is not None:
            x.extend(self.footer_expression._find_groups(tag))
        return x
    
    def features(self):
        """return a list of all features"""
        x = []
        if self.header_expression is not None:
            x.extend(self.header_expression.features())
        x.extend(self.expression.features())
        if self.footer_expression is not None:
            x.extend(self.footer_expression.features())
        return x

# Might be useful to allow a minimum record count (likely either 0 or 1)
class ParseRecords(PassThrough):
    def __init__(self, format_name, attrs, record_expression,
                 make_reader, reader_args = ()):
        PassThrough.__init__(self, Group(format_name,
                                         MaxRepeat(record_expression, 1),
                                         attrs))

        # I added attrs to the parameter list but couldn't make it
        # backwards compatible.  Without this check, it's possible to
        # have the object constructed seemingly okay then have the
        # error appear downstream, making it hard to track down.
        if isinstance(attrs, Expression):
            raise TypeError("Looks like you need an attrs between the format_name and the record_expression")
        
        self.format_name = format_name
        if attrs is None:
            attrs = xmlreader.AttributesImpl({})
        elif isinstance(attrs, type({})):
            attrs = xmlreader.AttributesImpl(attrs)
        self.attrs = attrs
        self.record_expression = record_expression
        self.make_reader = make_reader
        self.reader_args = reader_args
    def copy(self):
        """do a deep copy on this Expression tree"""
        return ParseRecords(self.format_name, self.attrs,
                            self.record_expression.copy(),
                            self.make_reader, self.reader_args)
    
    def make_parser(self, debug_level = 0):
        import Generate
        tagtable, want_flg, attrlookup = Generate.generate(
            self.record_expression, debug_level)

        return Parser.RecordParser(self.format_name, self.attrs,
                                   tagtable, (want_flg, debug_level, attrlookup),
                                   self.make_reader, self.reader_args)
    
    def make_iterator(self, tag, debug_level = 0):
        """create an iterator for this regexp; the 'tag' defines a record"""
        import Iterator
        if tag == self.format_name:
            return self.expression.make_iterator(self, tag)

        if isinstance(self.record_expression, Group) and \
           self.record_expression.name == tag and \
           IterParser is not None:
            # There's an optimization for this case
            return IterParser.IterRecords(
                self.record_expression.make_parser(debug_level),
                self.make_reader, self.reader_args, tag)

        return Iterator.IteratorRecords(
            self.record_expression.make_parser(debug_level),
            self.make_reader, self.reader_args, tag)
    
    
    def group_names(self):
        return self.format_name + self.expression.group_names()
    def _find_groups(self, tag):
        assert tag != self.format_name, "can't handle that case"
        return self.expression._find_groups(tag)

    def features(self):
        """return a list of all features"""
        return self.expression.features()

    def _modify_leaves(self, func):
        exp = self.expression.modify_leaves(func)
        assert exp is not None
        return self
    

# A sequence of characters: 'abcdef'
class Str(Expression):
    def __init__(self, s):
        """(s)

        Match the given string exactly (not as a regexp pattern)
        """
        self.string = s
        
    def copy(self):
        """do a deep copy on this Expression tree"""
        return Str(self.string)
    
    def __str__(self):
        """the corresponding pattern string"""
        return escape(self.string)



class ExpressionList(Expression):
    """shares implementation used by 'Expressions with subexpressions'"""
    def group_names(self):
        """the list of group names used by this Expression or its children"""
        names = {}
        for exp in self.expressions:
            for name in exp.group_names():
                names[name] = 1
        return tuple(names.keys())
    def _find_groups(self, tag):
        x = []
        for exp in self.expressions:
            x.extend(exp._find_groups(tag))
        return x
    def features(self):
        """return a list of all features"""
        x = []
        for exp in self.expressions:
            x.extend(exp.features())
        return x
    
    def _select_names(self, names):
        """internal function.  Do not use."""
        for exp in self.expressions:
            exp._select_names(names)
    def copy(self):
        """do a deep copy on this Expression tree"""
        return self.__class__(map(lambda x: x.copy(), self.expressions))
    def _modify_leaves(self, func):
        new_expressions = []
        for exp in self.expressions:
            new_expressions.append(exp._modify_leaves(func))
        assert None not in new_expressions
        self.expressions = tuple(new_expressions)
        return self

# A set of expressions: 'a|b|c'
class Alt(ExpressionList):
    """An Expression tree with a list of alternate matches.

    """
    def __init__(self, expressions):
        """(expressions)

        Match one of a list of alternate expressions.  The expressions are
        tested in their input order.

        For example, Alt( (exp1, exp2, exp3) ) means try to match exp1,
        and if that fails try to match exp2, and if that fails, try to
        match exp3.  If *that* fails, the match failed.

        """
        if isinstance(expressions, type( [] )):
            expressions = tuple(expressions)
        elif isinstance(expressions, Expression):
            raise TypeError("Must pass in a list of expressions, not just a single one (put it inside of ()s")
        self.expressions = expressions

    def __or__(self, other):
        # If the other is also an Alt, I can simplify things by
        # merging together the two lists of subexpressions.
        if isinstance(other, Alt):
            # This is why I convert lists to tuples; I need a
            # homogenous list type for addition.  I chose tuples to
            # help enforce the idea that Expressions should not be
            # changed after they are created.
            return Alt(self.expressions + other.expressions)
        else:
            return Alt(self.expressions + (other,))
        
    def __str__(self):
        """the corresponding pattern string"""
        return string.join(map(str, self.expressions), '|')


# A sequence of expressions: '[ab][cd][ef]'
class Seq(ExpressionList):
    """An Expression matching a set of subexpressions, in sequential order"""
    def __init__(self, expressions):
        """(expressions)

        Match the list of sequential expressions, in order.  Each
        expression starts matching at the point where the previous
        match finished.
        """
        if isinstance(expressions, type( [] )):
            # See 'Alt' for why I'm converting to tuples
            expressions = tuple(expressions)
        elif isinstance(expressions, Expression):
            raise TypeError("Must pass in a list of expressions, not just a single one (put it inside of ()s")
        self.expressions = expressions
        
    def __add__(self, other):
        # Optimize the case of Seq by combining the lists
        if isinstance(other, Seq):
            return Seq(self.expressions + other.expressions)
        return Seq(self.expressions + (other,))
    
    def __str__(self):
        """the corresponding pattern string"""
        # Seq is the lowest priority, so put parens around Alt subexpressions
        patterns = []
        for exp in self.expressions:
            pattern = str(exp)
            if isinstance(exp, Alt):
                patterns.append('(%s)' % pattern)
            else:
                patterns.append(pattern)
        return string.join(patterns, "")


########## Support code for making the pattern string for an expression

# taken from re.escape, except also don't escape " ="
def escape(pattern):
    "Escape all non-alphanumeric characters in pattern."
    result = list(pattern)
    alphanum=string.letters+'_'+string.digits+" ="
    for i in range(len(pattern)):
        char = pattern[i]
        if char not in alphanum:
            if char=='\000': result[i] = '\\000'
            else: result[i] = '\\'+char
    return string.join(result, '')

# Escapes for common control characters
_minimize_escape_chars = {
    "\a": r"\a",
    "\b": r"\b",
    "\n": r"\n",
    "\r": r"\r",
    "\t": r"\t",
    "\f": r"\f",
    "\v": r"\v",
    "[": "\\[",
    "]": "\\]",
    "\\": "\\\\",
    "^": "\^",
    }
    
def _minimize_escape_char(c):
    """(c) -> into an appropriately escaped pattern for the character"""
    x = _minimize_escape_chars.get(c)
    if x is None:
        if ord(c) < 32:
            return '\\' + c
        return c
    return x

def _minimize_escape_range(c1, c2):
    """(c1, c2) -> the pattern for the range bounded by those two characters"""
    # Called when 2 or more successive characters were found in a row.
    # c1 is the first character in the range and c2 is the last.

    if ord(c1) + 1 == ord(c2):
        # Two characters in a row.  Doesn't make sense to use the '-'
        # notation, so
        return _minimize_escape_char(c1) + _minimize_escape_char(c2)

    # Special case for numbers
    if c1 == '0' and c2 == '9':
        return r"\d"

    # Otherwise just use the '-' range.
    return _minimize_escape_char(c1) + '-' + _minimize_escape_char(c2)
    

def _minimize_any_range(s):
    """s -> a string useable inside [] which matches all the characters in s

    For example, passing in "0123456789" returns "\d".

    This code isn't perfect.
    """
    if not(isinstance(s, type(""))):
        s = str(s)
    if not s:
        return s

    # Treat the '-' special since it must occur last.
    # However, this means '!"#....xyz' gets turned into '!-,.-z-'
    has_hyphen = 0
    if '-' in s:
        has_hyphen = 1
        s = string.replace(s, "-", "")

    # Get the ordered list of characters in the string
    chars = list(s)
    chars.sort()
    unique = []
    prev_c = None
    for c in chars:
        if c != prev_c:
            unique.append(c)
            prev_c = c
    
    s = string.join(unique, "")
    t = ""
    prev = None
    prev_pos = 0
    pos = 0

    # Join successive characters which are in ASCII order
    # Eg, "abcdef" gets turned into "a-f"
    for c in unique:
        val = ord(c)
        
        if val - 1 != prev:
            # either beginning of string or non-sequential
            if prev is None:
                # beginning of string
                prev_pos = 0
            else:
                # non-sequential
                # Create the string for the previous range.
                if prev_pos == pos - 1:
                    # If there was one character in the range, use it
                    t = t + _minimize_escape_char(s[prev_pos])
                else:
                    # Two or more characters in a row define a range
                    t = t + _minimize_escape_range(s[prev_pos], s[pos-1])
                prev_pos = pos
        else:
            # Still part of the same sequence, so just advance in the string
            pass

        prev = val
        pos = pos + 1

    # Handle the final sequence block
    if s:
        if prev_pos == pos - 1:
            t = t + _minimize_escape_char(s[prev_pos])
        else:
            t = t + _minimize_escape_range(s[prev_pos], s[pos-1])
    else:
        # Get this case if there was no text except for the hyphen character
        pass

    # Put the hyphen back on the end
    if has_hyphen:
        t = t + '-'

    # Simple fixes for fields that annoy me a lot
    conversions = {
        "\\dA-Z_a-z\xc0-\xd6\xd8-\xf6\xf8-\xff": r"\w",
        }
    t = conversions.get(t, t)
    
    return t

def _make_no_case(node):
    """modify an expression in place to remove case dependencies

    may return a new top-level node
    """
    if isinstance(node, Str):
        x = NullOp()
        s = ""
        for c in node.string:
            up_c = string.upper(c)
            low_c = string.lower(c)
            assert c in (up_c, low_c), "how can this be?"
            if up_c == low_c:
                s = s + c
            else:
                if s:
                    x = x + Str(s)
                    s = ""
                x = x + Any(up_c + low_c)
        if s:
            x = x + Str(s)
        return x
    
    if isinstance(node, Any):
        s = node.chars
        chars = {}
        for c in s:
            chars[c] = 1
        for c in string.upper(s) + string.lower(s):
            if not chars.has_key(c):
                chars[c] = 1
                s = s + c
        return Any(s, node.invert)

    if isinstance(node, Literal):
        c = node.char
        up_c = string.upper(c)
        low_c = string.lower(c)
        if up_c == low_c:
            return node
        return Any(up_c + low_c, node.invert)

    return node

def NoCase(expr):
    """expression -> expression where the text is case insensitive"""
    expr = expr.copy()
    return expr._modify_leaves(_make_no_case)
