# Copyright 2000-2001, Dalke Scientific Software, LLC
# Distributed under the Biopython License Agreement (see the LICENSE file).

# Create the tables needed for mxTextTools

import string

try:
    from mx import TextTools as TT
except ImportError:
    import TextTools as TT

import msre_parse  # Modified version of Secret Labs' sre_parse
import Expression, convert_re

import Parser

# This was added with mxTextTools 1.2
supports_lookahead = hasattr(TT, "LookAhead")

# NOTE:
#  I use some notation in the comments below, which I'll define here

#   "+1/-1 trick"
#  mxTextTools doesn't like "Call"s which don't consume any
#  characters.  To get around that, I have successful matches return
#  the end position + 1, so that it seems to consume a character.
#  Then I Skip -1 characters to get back to the correct loction.
#    mxTextTools 1.2 and later has a 'LookAhead' command so this
#  trick is not needed for that version

#    ">ignore"
#  Sometimes the tagtables descends into another tagtable.
#  mxTextTools requires there be a tag name for all such descents, for
#  otherwise it cannot create the taglist for the descent.  If there
#  isn't a name, the fake name ">ignore" is used, which is an illegal
#  XML name.  The Parser knows to ignore this tag when generating
#  callbacks.


# With the current implementations, each of the generate_* functions
# return a Parser with a tagtable that starts parsing with the first
# element in the table and succeeds by jumping one past the last
# element.

# This is single threaded!
_generate_count = 0

class GeneratorState:
    def __init__(self, groupref_names, debug_level):
        self.groupref_names = groupref_names
        self.debug_level = debug_level
        self.lookup = {}
        
    def add_group(self, group):
        name = group.name
        attrs = group.attrs
        tag = self.new_group_tag()
        self.lookup[tag] = (name, attrs)
        return tag
    
    def new_group_tag(self):
        global _generate_count
        i = _generate_count
        tag = ">G%d" % i
        i = i + 1
        _generate_count = i
        return tag


# ==============

#  a|b|c|d|...
# Implemented by creating N subtables.  If table i fails, try i+i.  If
# it succeeds, jump to 1 past the end.
#   1. try table1: on fail +1, on success +N+1
#   2. try table2: on fail +1, on success +N
#   3. try table3: on fail +1, on success +N-1
#   ...
#   N. try tableN: on fail +1, on success +2
# N+1. Fail
# N+2. <past the end of table>
#   
# XXX Don't need to create a Table.  Could use a single tagtable by
# merging into one table and retargeting "succeeded" transitions
# (which jumped to the end+1 of the given table) to transition to the
# end of the merged table.
#
def generate_alt(expression, genstate):
    # Get the tagtables for each of the subexpressions
    tables = []
    for expr in expression.expressions:
        tables.append(_generate(expr, genstate))

    # Make a set of branches testing each one
    i = 0
    n = len(tables)
    result = []
    for table in tables:
        result.append( \
            (">ignore", TT.Table, tuple(table), +1, n-i+1)
            )
        i = i + 1
    # Jumps here
    result.append( (None, TT.Fail, TT.Here) )
    return result

#  sequence of successive regexps: abcdef...
#
# Simply catenate the tagtables together, in order.  Works because falling
# off the end of one table means success for that table, so try the next.
def generate_seq(expression, genstate):
    result = []
    if genstate.debug_level == 0:
        for exp in expression.expressions:
            table = _generate(exp, genstate)
            result.extend(table)
    elif genstate.debug_level == 1:
        for exp in expression.expressions:
            table = _generate(exp, genstate)
            result.extend(table)
            result.append( (None, TT.Call, track_position, +1, +1) )
    elif genstate.debug_level == 2:
        for exp in expression.expressions:
            table = _generate(exp, genstate)
            result.extend(table)
            result.append( (None, TT.Call, track_position, +1, +1) )
            table.append( (None, TT.Call, print_info(exp), +1, +1) )
            
    return result

# A literal character, or a character which isn't the given character
def generate_literal(expression, genstate):
    if expression.invert:
        # Can't use "InNot" since it can match EOF
        return [(None, TT.IsIn, convert_re.invert(expression.char))]
                             
    else:
        return [(None, TT.Is, expression.char)]

# A string
def generate_str(expression, genstate):
    return [(None, TT.Word, expression.string)]

# Any character in the given list of characters
def generate_any(expression, genstate):
    if expression.invert:
        # I don't use "IsNotIn" since that allows EOF
        return [(None, TT.IsIn, convert_re.invert(expression.chars))]
                             
    else:
        return [(None, TT.IsIn, expression.chars)]

# Support class for group matches.  Stores the contents of a match so
# it can be used by a future back reference or named group repeat.
class SetGroupValue:
    def __init__(self, name):
        self.name = name
    def __call__(self, taglist, text, l, r, subtags):
        taglist.append( (self.name, l, r, subtags) )
        Parser._match_group[self.name] = text[l:r]

# A 'group', either named or unnamed
def generate_group(expression, genstate):
    tagtable = _generate(expression.expression, genstate)

    name = expression.name
    if name is None:
        assert not expression.attrs, "unnamed group can't have attrs!"
        # Don't really need to do anything, do I?
        return tagtable
        
    if genstate.groupref_names.get(name) != 1:
        if expression.attrs:
            name = genstate.add_group(expression)
        return [(name, TT.Table, tuple(tagtable)) ]
    else:
        # generate the code to save the match information
        if expression.attrs:
            name = genstate.add_group(expression)
        return [(SetGroupValue(name), TT.Table+TT.CallTag,
                                      tuple(tagtable)), ]

# Uglyness to handle named group repeats.

# There are two entry points for this class.  One is 'call', which is
# called by the TextTools 'Call' tag command to test for a match.  The
# other is 'calltag', which is called by the TextTools 'CallTag' tag
# command after a match is found.

# Uses the +1/-1 trick since the match can be of size 0.

# Store the taglist on success as an instance variable so the
# "CallTag" can append the matches.  Again, this is SINGLE-THREADED
# CODE, but the best I can do with the current version of mxTextTools.

class HandleRepeatCount:
    def __init__(self, tagtable, min_count, max_count):
        self.tagtable = tagtable    # The regexp to match
        self.min_count = min_count  # For now, min_count and max_count
        self.max_count = max_count  #   must be the same string

        self.taglist = None

    def _get_ranges(self):
        # Look up the min/max counts and convert to integers
        min_count = self.min_count
        if type(min_count) == type(""):
            min_count = Parser._match_group[min_count] # group must already exist!
            min_count = string.atoi(min_count)  # requires an integer!
        
        max_count = self.max_count
        if type(max_count) == type(""):
            max_count = Parser._match_group[max_count] # group must already exist!
            max_count = string.atoi(max_count)  # requires an integer!

        return min_count, max_count
    
    def call(self, text, x, end):
        # Called by 'TextTools.Call' to detect a match.
        # I do the full match here and store the results for later use.
        # If successful, I return x+1, else return x+0 (the +1/-1 trick)
        min_count, max_count = self._get_ranges()
        assert min_count == max_count, \
               "cannot have different sizes: %s %s" % (min_count, max_count)
        
        tagtable = self.tagtable * min_count
        result, taglist, pos = TT.tag(text, tagtable, x, end)
        if result == 1:
            # Store the taglist for later use
            self.taglist = taglist
            return pos + 1  # +1 because {0} is allowed; Skip -1 later
        else:
            self.taglist = None
            return x
        
    def calltag(self, taglist, text, l, r, subtags):
        # Called by 'CallTag' to say I matched the above 'call'
        assert not subtags, repr(subtags)
        # Append the taglist which we saved from an earlier good match.
        # The '-1' is because of the +1/-1 hack.
        taglist.append( (">ignore", l, r-1, self.taglist) )
        

# These objects forward calls to HandleRepeatCount methods.
# Would like to start the methods directly but methods cannot be pickled.

class _call_calltag:
    def __init__(self, obj):
        self.obj = obj
    def __call__(self, taglist, text, l, r, subtags):
        return self.obj.calltag(taglist, text, l, r, subtags)

class _call_call:
    def __init__(self, obj):
        self.obj = obj
    def __call__(self, text, x, end):
        return self.obj.call(text, x, end)
    
def generate_named_max_repeat(expression, genstate):
    if type(expression.min_count) != type("") or \
       type(expression.max_count) != type(""):
        raise NotImplementedError("Cannot mix numeric and named repeat counts")
    if expression.min_count != expression.max_count:
        raise NotImplementedError("Only a single named repeat count allowed")
    
    tagtable = _generate(expression.expression, genstate)
    counter = HandleRepeatCount(tuple(tagtable),
                                expression.min_count,
                                expression.max_count)

    # You see the +1/-1 trick implemented here
    # Don't use the 'supports_lookahead' here since that complicates matters
    return \
           [(_call_calltag(counter), TT.Call + TT.CallTag,
             _call_call(counter), TT.MatchFail), 
            (None, TT.Skip, -1, TT.MatchFail),]

# It isn't as bad as it looks :)
# Basically, call some other code for named group repeats.
# Everything else is of the form {i,j}.
# Get the subexpression table:
#    generate i copies which must work
#    generate j-i copies which may work, but fail okay.
#    special case when j == 65535, which standard for "unbounded"
def generate_max_repeat(expression, genstate):
    expr = expression.expression
    min_count = expression.min_count
    max_count = expression.max_count

    # Call some other code for named group repeats
    if type(min_count) == type("") or type(max_count) == "":
        return generate_named_max_repeat(expression, genstate)

    assert 0 <= min_count <= max_count, "bad ranges (%sd, %d)" % \
           (min_count, max_count)

    tagtable = _generate(expr, genstate)
    result = []

    # Must repeat at least "i" times.
    for i in range(min_count):
        result.append( (None, TT.SubTable, tuple(tagtable)) )

    # Special case for when the max count means "unbounded"
    if max_count == msre_parse.MAXREPEAT:
        result.append( (None, TT.SubTable, tuple(tagtable),
                                  +1, 0))
    elif min_count == max_count:
        # Special case when i == j
        pass
    else:
        # Generate j-i more tagtables, but allow failures
        offset = max_count - min_count
        for i in range(offset):
            result.append( (">ignore", TT.Table,
                              tuple(tagtable),
                              +offset, +1) )
            offset = offset - 1
    return result

# Doesn't do anything
def generate_null_op(expression, genstate):
    return []

class print_debug:
    """Print debug information"""
    def __init__(self, msg):
        self.msg = msg
    def __call__(self, text, x, end):
        print "Martel:", self.msg
        return x

def generate_debug(expression, genstate):
    return [(None, TT.Call, print_debug(expression.msg), +1, +1)]

# XXX Is this correct?  This is the multiline behaviour which allows
# "^" to match the beginning of a line.
def check_at_beginning(text, x, end):
    if x == 0:
        return x
    if text[x-1] == "\n":
        return x
    return x + 1

# XXX Consider this code broken!
#
# Uses the +1/-1 trick.
def generate_at_beginning(expression, genstate):
    if supports_lookahead:
        return [(None, TT.Call + TT.LookAhead, check_at_beginning, +1,
                 TT.MatchFail),]
    else:
        return [(None, TT.Call, check_at_beginning, +2, +1),
                (None, TT.Skip, -1, TT.MatchFail, TT.MatchFail),
                ]


# XXX Consider this code broken!
#
# XXX If 'check_at_beginning' is correct, then this is wrong since it
# doesn't implement the multiline behaviour.  
def generate_at_end(expression, genstate):
    return [(None, TT.EOF, TT.Here)]


# Match any character except newline (by which I mean just "\012")
def generate_dot(expression, genstate):
    return [(None, TT.IsInSet, TT.invset('\n')), ]

# Match any of the three standard newline conventions
def generate_eol(expression, genstate):
    return [(None, TT.Is, '\n', +1, +3),
            (None, TT.Is, '\r', TT.MatchFail, +1),
            (None, TT.Is, '\n', +1, +1),
            ]

# Used during a negative lookhead test.
# Note the +1/-1 trick.
def check_assert_not(text, x, end, tagtable):
    result, taglist, pos = TT.tag(text, tagtable, x, end)
    if result:
        # This failed
        return x
    # On success, move forward 1, to be removed later
    return x + 1

# Called by the 'Call' tag command when doing negative lookaheads
class CheckAssertNot:
    def __init__(self, tag_words):
        self.tag_words = tag_words
    def __call__(self, text, x, end):
        pos = check_assert_not(text, x, end, self.tag_words)
        return pos


# Used during a positive lookhead test.
# Note the +1/-1 trick.
def check_assert(text, x, end, tag_words):
    result, taglist, pos = TT.tag(text, tag_words, x, end)
    if result:
        # This succeeded, move forward 1, to be removed later
        return x+1
    # failed
    return x

# Called by the 'Call' tag command when doing negative lookaheads
class CheckAssert:
    def __init__(self, tag_words):
        self.tag_words = tag_words
    def __call__(self, text, x, end):
        pos = check_assert(text, x, end, self.tag_words)
        return pos

# Create the tagtable for doing a lookahead assertion.
# Uses the +1/-1 trick.
def generate_assert(expression, genstate):
    tagtable = _generate(expression.expression, genstate)
    if expression.invert:
        func = CheckAssertNot
    else:
        func = CheckAssert
    if supports_lookahead:
        return [
            (None, TT.Call + TT.LookAhead, func(tuple(tagtable)),
             TT.MatchFail),
            ]
    else:
        return [
            (None, TT.Call, func(tuple(tagtable)),
             TT.MatchFail),
            (None, TT.Skip, -1, TT.MatchFail),
            ]

# Create an object which can be called by a 'Call' tag command to
# match the text found by the named group.
# Uses the +1/-1 trick.
class CheckGroupRef:
    def __init__(self, name):
        self.name = name
    def __call__(self, text, x, end):
        match_text = Parser._match_group[self.name]  # group name not yet found
        if text[x:x+len(match_text)] != match_text:
            return x
        return x+len(match_text)+1

# Make the tagtable needed for a named group backreference.
# Uses the +1/-1 trick.
def generate_groupref(expression, genstate):
    # Look up the string from the match group
    # It can be of length 0 or more, so use the +1/-1 trick.
    return [
        (None, TT.Call, CheckGroupRef(expression.name), TT.MatchFail),
        (None, TT.Skip, -1, TT.MatchFail),
        ]

# Used to define parsers which read a record at time.  They contain no
# parse information themselves, but only in their children
def generate_pass_through(expression, genstate):
    return _generate(expression.expression, genstate)

# Mapping from Expression type to function used to generate the tagtable
generate_table = {
    Expression.Alt: generate_alt,
    Expression.Any: generate_any,
    Expression.Assert: generate_assert,
    Expression.AtBeginning: generate_at_beginning,
    Expression.AtEnd: generate_at_end,
    Expression.Debug: generate_debug,
    Expression.Dot: generate_dot,
    Expression.AnyEol: generate_eol,
    Expression.Group: generate_group,
    Expression.GroupRef: generate_groupref,
    Expression.Literal: generate_literal,
    Expression.MaxRepeat: generate_max_repeat,
    Expression.NullOp: generate_null_op,
    Expression.Seq: generate_seq,
    Expression.Str: generate_str,
}

_position = -1
def track_position(text, x, end):
    """store the start position of the farthest successful match

    This value is more useful than mxTextTools' default, which only
    points out the last text region successfully tagged at the top
    level.  This value is the last region successfully tagged
    anywhere.

    Uses a global variable so this is SINGLE THREADED!
    
    """
    
    global _position
    _position = max(x, _position)
    return x

class print_info:
    """Print information after each expression match"""
    def __init__(self, expression):
        s = str(expression)
        if len(s) > 40:
            s = s[:17] + " ... " + s[-17:]
        self.msg = s
    def __call__(self, text, x, end):
        print "Match %s (x=%d): %s" % (repr(text[max(0, x-8):x+8]), x,
                                            repr(self.msg))
        return x

# The internal recursive call
def _generate(expression, genstate):
    try:
        func = generate_table[expression.__class__]
    except KeyError:
        if isinstance(expression, Expression.PassThrough):
            func = generate_pass_through
        else:
            raise AssertionError, \
                  "Unknown Expression object: %s" % repr(expression)
    table = func(expression, genstate)

    if genstate.debug_level == 0 or not table:
        pass
    elif genstate.debug_level == 1:
        table.append( (None, TT.Call, track_position, +1, +1) )
    elif genstate.debug_level == 2:
        table.append( (None, TT.Call, track_position, +1, +1) )
        table.append( (None, TT.Call, print_info(expression), +1, +1) )
    else:
        raise AssertionError, "Unknown debug level: %s" % genstate.debug_level
        
    return table

# Get the tagtable and the want_groupref_names
# Main entry point for record oriented readers
def generate(expression, debug_level = 0):
    """expression -> Parser for the Expression tree"""
    groupref_names = _find_wanted_groupref_names(expression)
    genstate = GeneratorState(groupref_names, debug_level)
    tagtable = _generate(expression, genstate)
    if groupref_names:
        want_groupref_names = 1
    else:
        want_groupref_names = 0
    return tuple(tagtable), want_groupref_names, genstate.lookup

# Get the parser.  Main entry point for everything except record
# oriented readers
def generate_parser(expression, debug_level = 0):
    tagtable, want_groupref_names, attrlookup = generate(expression)
    return Parser.Parser(tagtable, (want_groupref_names, debug_level,
                                    attrlookup))


# This code is ugly.  There are couple possible fixes:
#  1) write a tree iterator which visits each node of the Expression tree.
#       Using it would remove the extra code for doing recursion.
#  2) add a new method to each Expression type to get the names.
#       Using it would remove the many 'if' statements.
#
# There will probably be code in the future to optimize an Expression
# tree.  That code will have tools to manipulate a tree.  I'll wait
# until then to work decide which option to take, or if there's a
# third.

def _find_wanted_groupref_names(expression):
    """expression -> dict of group names wanted by elements of the tree

    The dict is used to during tagtable generation to specify which
    groups need to save their match text.  There's match-time overhead
    for doing that, and the code isn't thread safe, so the intent is
    to save only those groups that are needed.

    The dict value is 1 if the group name is needed, else there is
    no entry in the dict.

    XXX need to make this a method!
    """
    want_names = {}
    if isinstance(expression, Expression.Alt) or \
       isinstance(expression, Expression.Seq):
        for x in expression.expressions:
            want_names.update(_find_wanted_groupref_names(x))

    elif isinstance(expression, Expression.Group) or \
         isinstance(expression, Expression.Assert) or \
         isinstance(expression, Expression.PassThrough):
        want_names.update(_find_wanted_groupref_names(expression.expression))

    elif isinstance(expression, Expression.MaxRepeat):
        if type(expression.min_count) == type(""):
            want_names[expression.min_count] = 1
        if type(expression.max_count) == type(""):
            want_names[expression.max_count] = 1
        want_names.update(_find_wanted_groupref_names(expression.expression))

    elif isinstance(expression, Expression.GroupRef):
        want_names[expression.name] = 1

    elif isinstance(expression, Expression.Literal) or \
         isinstance(expression, Expression.Str) or \
         isinstance(expression, Expression.Any) or \
         isinstance(expression, Expression.AtBeginning) or \
         isinstance(expression, Expression.AtEnd) or \
         isinstance(expression, Expression.Dot) or \
         isinstance(expression, Expression.AnyEol) or \
         isinstance(expression, Expression.Debug) or \
         isinstance(expression, Expression.NullOp):
        pass

    else:
        raise NotImplementedError, "What is a %s?" % repr(expression)

    return want_names
