# Copyright 2000-2001, Dalke Scientific Software, LLC
# Distributed under the Biopython License Agreement (see the LICENSE file).

"""Classes for nodes in the Expression tree

  Expression
   |--- Any           - match (or don't match) a set of characters
   |--- AnyEol        - match any newline representation ("\n", "\r" or "\r\n")
   |--- Assert        - used for positive and negative lookahead assertions 
   |--- AtBeginning   - match the beginning of a line
   |--- AtEnd         - match the end of a line
   |--- Dot           - match any character except newline
   |--- Group         - give a group name to an expression
   |--- GroupRef      - match a previously identified expression
   |--- Literal       - match (or don't match) a single character
   |--- MaxRepeat     - greedy repeat of an expression, within min/max bounds
   |--- NullOp        - does nothing (useful as an initial seed)
   |--- PassThrough   - used when overriding 'make_parser'; match its subexp
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

    def make_iterator(self, tag, debug_level = 0):
        """create an iterator for this regexp; the 'tag' defines a record"""
        import Iterator
        return Iterator.Iterator(self.make_parser(debug_level), tag)
    
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

    def _select_names(self, names):
        """internal function: do not use"""
        if self.name is not None and self.name not in names:
            self.name = None
        self.expression._select_names(names)

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
    
    def copy(self):
        """do a deep copy on this Expression tree"""
        return MaxRepeat(self.expression.copy(), self.min_count,
                         self.max_count)
    
    def _select_names(self, names):
        """internal function: do not use"""
        self.expression._select_names(names)
        
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
    def group_names(self):
        return ()
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
    def copy(self):
        """do a deep copy on this Expression tree"""
        return PassThrough(self.expression)
    def __str__(self):
        """the corresponding pattern string"""
        return str(self.expression)
    def group_names(self):
        return self.expression.group_names()

class HeaderFooter(PassThrough):
    def __init__(self, format_name,
                 header_expression, make_header_reader, header_args,
                 record_expression, make_record_reader, record_args,
                 footer_expression, make_footer_reader, footer_args):
        if header_expression is None:
            assert make_header_reader is None and header_args is None
            exp = MaxRepeat(record_expression, 1)
        else:
            exp = header_expression + MaxRepeat(record_expression, 1)
        if footer_expression is not None:
            exp = exp + footer_expression
        else:
            assert make_footer_reader is None and footer_args is None
        PassThrough.__init__(self, Group(format_name, exp))

        self.format_name = format_name
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
            self.format_name,
            header_exp, self.make_header_reader, self.header_args,
            record_exp, self.make_record_reader, self.record_args,
            footer_exp, self.make_footer_reader, self.footer_args)

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
            self.format_name,
            make_header_reader, header_args, header_tagtable,
            make_record_reader, record_args, record_tagtable,
            make_footer_reader, footer_args, footer_tagtable,
            (want, debug_level, attrlookup))

    def make_iterator(self, tag, debug_level = 0):
        raise NotImplementedError, "Need an IteratorHeaderFooter"
        return self.expression.make_iterator(self, tag)

    def group_names(self):
        return self.expression.group_names()

# Might be useful to allow a minimum record count (likely either 0 or 1)
class ParseRecords(PassThrough):
    def __init__(self, format_name, record_expression,
                 make_reader, reader_args = ()):
        PassThrough.__init__(self, Group(format_name,
                                         MaxRepeat(record_expression, 1)))
        self.format_name = format_name
        self.record_expression = record_expression
        self.make_reader = make_reader
        self.reader_args = reader_args
    def copy(self):
        """do a deep copy on this Expression tree"""
        return ParseRecords(self.format_name, self.record_expression.copy(),
                            self.make_reader, self.reader_args)
    
    def make_parser(self, debug_level = 0):
        import Generate
        tagtable, want_flg, attrlookup = Generate.generate(
            self.record_expression, debug_level)

        return Parser.RecordParser(self.format_name,
                                   tagtable, (want_flg, debug_level, attrlookup),
                                   self.make_reader, self.reader_args)
    
    def make_iterator(self, tag, debug_level = 0):
        """create an iterator for this regexp; the 'tag' defines a record"""
        import Iterator
        if tag == self.format_name:
            return self.expression.make_iterator(self, tag)
        return Iterator.IteratorRecords(
            self.record_expression.make_parser(debug_level),
            self.make_reader, self.reader_args, tag)
    
    
    def group_names(self):
        return self.expression.group_names()

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
    def _select_names(self, names):
        """internal function.  Do not use."""
        for exp in self.expressions:
            exp._select_names(names)
    def copy(self):
        """do a deep copy on this Expression tree"""
        return self.__class__(map(lambda x: x.copy(), self.expressions))


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
        if type(expressions) == type( [] ):
            expressions = tuple(expressions)
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
        if type(expressions) == type( [] ):
            # See 'Alt' for why I'm converting to tuples
            expressions = tuple(expressions)
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

# taken from re.escape, except also don't escape " "
def escape(pattern):
    "Escape all non-alphanumeric characters in pattern."
    result = list(pattern)
    alphanum=string.letters+'_'+string.digits+" "
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

    This code isn't perfect.  For example, "_ABC...Zabc...z" should
    convert to "\w".
    """
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
    
    return t


