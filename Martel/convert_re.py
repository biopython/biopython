# Copyright 2000-2001, Dalke Scientific Software, LLC
# Distributed under the Biopython License Agreement (see the LICENSE file).

"""Converts a regular expression pattern string into an Expression tree.

This is not meant to be an externally usable module.

This works by using msre_parse.py to parse the pattern.  The result is
a tree data structure, where the nodes in the tree are tuples.  The
first element of the tuple is the name of the node type.  The format
of the other elements depends on the type.

The conversion routine is pretty simple - convert each msre_parse tuple
node into a Martel Expression node.  It's a recusive implementation.

'msre_parse.py' is a modified version of Secret Labs' 'sre_parse.py'

"""

import string
import msre_parse, Expression

# The msre_parse parser uses a "master pattern object" which keeps
# track of the mapping from group id to group name.  This is okay for
# msre_parse because they only track the last group with a given name.
# I need to get all groups with the same name, so I need a new object
# which stores them in a list that I can use later.

class GroupNames:
    def __init__(self):
        self.flags = 0
        self.open = []
        self.groups = 1
        self.groupdict = {}
    def opengroup(self, name=None):
        gid = self.groups
        self.groups = gid + 1
        if name:
            self.groupdict[name] = self.groupdict.get(name, []) + [gid]
        self.open.append(gid)
        return gid
    def closegroup(self, gid):
        self.open.remove(gid)
    def checkgroup(self, gid):
        return gid < self.groups and gid not in self.open

    def reverse_name(self, id):
        """group number -> group name, or None if there is no name"""
        for key, val in self.groupdict.items():
            if id in val:
                return key
        # Ignore non-named groups



# Convert a 'literal' tuple into a Literal object
def convert_literal(group_names, name, val):
    return Expression.Literal(chr(val), 0)

# Convert a 'not_literal' tuple into a Literal object
def convert_not_literal(group_names, name, val):
    return Expression.Literal(chr(val), 1)

# Convert an 'at_beginning" tuple into an AtBeginning object
# Convert an 'at_end" tuple into an AtEnd object
def convert_at(group_names, name, where):
    if where == "at_beginning":
        return Expression.AtBeginning()
    elif where == "at_end":
        return Expression.AtEnd()
    raise AssertionError("Unknown at name: %s" % repr(where))

# Convert an 'any' tuple into a Dot object
def convert_any(group_names, name, ignore):
    assert ignore is None, "what does it mean when the field is '%s'?" % ignore
    return Expression.Dot()

# Convert an 'assert' tuple into a Assert object, as a positive assertion
def convert_assert(group_names, name, (direction, terms)):
    assert direction == 1, "does not support lookbehind"
    return Expression.Assert(convert_list(group_names, terms), 0)

# Convert an 'assert_not' tuple into a Assert object, as a negative assertion
def convert_assert_not(group_names, name, (direction, terms)):
    assert direction == 1, "does not support lookbehind"
    return Expression.Assert(convert_list(group_names, terms), 1)


# Convert a 'branch' tuple into an Alt object
def convert_branch(group_names, name, (ignore, branches)):
    assert ignore is None, "what is %s?" % repr(ignore)
    results = []
    for branch in branches:
        results.append(convert_list(group_names, branch))
    if len(results) == 1:
        return results[0]
    return Expression.Alt(tuple(results))

# I know, it's only good for ASCII...
def invert(s):
    """s -> a string containing all the characters not present in s"""
    letters = []
    if not(isinstance(s, type(""))):
        s = str(s)
    for c in map(chr, range(256)):
        if c not in s:
            letters.append(c)
    return string.join(letters, "") 

# Map from the msre_parse category names into actual characters.
# I can do this here since I don't worry about non-ASCII character sets.
categories = {
    "category_word": string.letters + "0123456789_",
    "category_digit": string.digits,
    "category_space": "\t\n\v\f\r ",
    "category_newline": "\n\r",
 
    "category_not_word": invert(string.letters + "0123456789_"),
    "category_not_digit": invert(string.digits),
    "category_not_space": invert("\t\n\v\f\r "),
    }

# Convert an 'in' tuple into an Any object
# Pass in the negate flag if given.
def convert_in(group_names, name, terms):
    negate = (terms[0][0] == 'negate')
    s = ""
    for c in terms[negate:]:
        if c[0] == 'literal':
            s = s + chr(c[1])
        elif c[0] == 'range':
            for i in range(c[1][0], c[1][1]+1):
                s = s + chr(i)
        elif c[0] == 'category':
            s = s + categories[c[1]]
        else:
            raise AssertionError("unknown option for 'in': %s" % c[0])
    return Expression.Any(s, negate)


# Convert a 'subpattern' tuple into a Group object
def convert_subpattern(group_names, name, (id, terms)):
    pattern_name = group_names.reverse_name(id)

    # The name in the ?P<group> may contain attr information
    # serialized in a URL-encoded form; if present, deconvolute.
    pos = -1
    attrs = {}
    if pattern_name is not None:
        pos = string.find(pattern_name, "?")
    
    if pos != -1:
        import cgi
        qs = pattern_name[pos+1:]
        if not qs:
            # cgi.parse_qs doesn't like parsing the empty string
            attrs = {}
        else:
            attrs = cgi.parse_qs(pattern_name[pos+1:],
                                 keep_blank_values = 1,
                                 strict_parsing = 1)
        pattern_name = pattern_name[:pos]

        for k, v in attrs.items():
            if len(v) != 1:
                raise AssertionError(
"The attribute name %s was found more than once (%d times) in the tag %s" %
            (repr(k), len(v), repr(pattern_name)))

            attrs[k] = v[0]
    
    return Expression.Group(pattern_name, convert_list(group_names, terms),
                            attrs)

# Convert a 'newline' tuple into an AnyEol object
def convert_newline(group_names, name, ignore):
    assert ignore is None, "what does it mean when field is %s?" % `ignore`
    return Expression.AnyEol()

# Convert a 'max_repeat' tuple into a MaxRepeat object
def convert_max_repeat(group_names, name, (min_count, max_count, terms)):
    return Expression.MaxRepeat(convert_list(group_names, terms),
                                min_count, max_count)

# Convert a 'groupref' tuple into a GroupRef object
def convert_groupref(group_names, name, id):
    assert type(id) != type(0), \
           "Martel cannot use numbered group reference: %d" % id
    # msre_parse returns the list from the GroupNames
    # Map that back to a number.
    pattern_name = group_names.reverse_name(id[0])
    return Expression.GroupRef(pattern_name)

# Map from the tuple typename into the function used to convert the
# tuple into an Expression.

converter_table = {
    "any": convert_any,
    "assert": convert_assert,
    "assert_not": convert_assert_not,
    "at": convert_at,
    "branch": convert_branch,
    "groupref": convert_groupref,
    "in": convert_in,
    "literal": convert_literal,
    "max_repeat":  convert_max_repeat,
    "newline": convert_newline,
    "not_literal": convert_not_literal,
    "subpattern":  convert_subpattern,
    }

# Convert a list of msre_parse tuples into a Seq
def convert_list(group_names, terms):
    # This is always a sequence of terms
    results = []
    for term in terms:
        name = term[0]
        try:
            func = converter_table[name]
        except KeyError:
            raise AssertionError, "Do not understand sre expression %s" % \
                  repr(name)
                  
        results.append( func(*(group_names,) + term) )
    if len(results) == 1:
        return results[0]
    return Expression.Seq(tuple(results))
        
    
# Primary entry point
def make_expression(pattern):
    """pattern -> the Expression tree for the given pattern string"""

    # In the following, the "pattern =" and "x.pattern" are the names
    # used by msre_parse.  They have nothing to do the input pattern.
    
    # Make the msre_parse tuple tree from the string ...
    x = msre_parse.parse(str = pattern, pattern = GroupNames())

    # ... and convert it into an Expression
    return convert_list(x.pattern, x)
