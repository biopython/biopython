# Copyright 2000 by Andrew Dalke.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


# The Prosite patterns are defined at http://www.expasy.ch/txt/prosuser.txt
# 
# The PA  (PAttern) lines  contains the definition of a PROSITE pattern. The
# patterns are described using the following conventions:
# 
#    -  The standard IUPAC one-letter codes for the amino acids are used.
#    -  The symbol `x' is used for a position where any amino acid is accepted.
#    -  Ambiguities are  indicated by  listing the acceptable amino acids for a
#       given position,  between square  parentheses `[  ]'. For example: [ALT]
#       stands for Ala or Leu or Thr.
#    -  Ambiguities are  also indicated  by listing  between a  pair  of  curly
#       brackets `{  }' the  amino acids  that are  not  accepted  at  a  given
#       position. For  example: {AM}  stands for  any amino acid except Ala and
#       Met.
#    -  Each element in a pattern is separated from its neighbor by a `-'.
#    -  Repetition of  an element  of the pattern can be indicated by following
#       that element  with a  numerical value  or  a  numerical  range  between
#       parenthesis. Examples: x(3) corresponds to x-x-x, x(2,4) corresponds to
#       x-x or x-x-x or x-x-x-x.
#    -  When a  pattern is  restricted to  either the  N- or  C-terminal  of  a
#       sequence, that  pattern either starts with a `<' symbol or respectively
#       ends with a `>' symbol.
#    -  A period ends the pattern.
# 
# That boils down to doing these conversions
# 
# [] -> []
# {} -> [^ ]
# -  ->
# () -> {}
# <  -> ^
# >  -> $
# x->X
# . ->

import string, re
import Alphabet, Seq

_prosite_trans = string.maketrans("abcdefghijklmnopqrstuvwxyzX}()<>",
                                  "ABCDEFGHIJKLMNOPQRSTUVW.YZ.]{}^$")

# These fix the known problems in release 16.0 of July 1999.
# To use it, try:
#   pattern = fix_errors.get(pattern, pattern)
fix_errors = {
  "F-[GSTV]-P-R-L-[G>].": "F-[GSTV]-P-R-L-G>.",        # PS00539
  "F-[IVFY]-G-[LM]-M-[G>].": "F-[IVFY]-G-[LM]-M-G>.",  # PS00267
  }

# Both the Prosite pattern and match result act like sequences.
class PrositeAlphabet(Alphabet.Alphabet):
    pass
prosite_alphabet = PrositeAlphabet()

def compile(pattern):
    if not verify_pattern(pattern):
        raise TypeError, "not a legal prosite pattern"
    return Prosite(pattern = pattern)

class Prosite:
    alphabet = prosite_alphabet
    def __init__(self, pattern = None, data = None):
        assert (pattern is None and data is not None) ^ \
               (pattern is not None and data is None), \
               "one and only one of pattern and data can have a value"
        if pattern is not None:
            self.pattern = pattern
        if data is not None:
            self.data = data

    def __repr__(self):
        return "Prosite(%s)" % repr(str(self))
    def __str__(self):
        return string.join(map(str, self.data), "-") + "."
    def __len__(self): return len(self.data)
    def __getitem__(self, i): return self.data[i]
    def __getslice__(self, i, j):
        i = max(i, 0); j = max(j, 0)
        return Prosite(data = self.data[i:j])
    def __getattr__(self, name):
        # Lazy creation of these elements / cache results
        if name == "re":
            self.re = re.compile(prosite_to_re(self.pattern))
            return self.re
        elif name == "grouped_re":
            self.grouped_re = re.compile(prosite_to_grouped_re(self.pattern))
            return self.grouped_re
        elif name == "data":
            self.data = find_terms(self.pattern)
            return self.data
        elif name == "pattern":
            self.pattern = str(self)
            return self.pattern
        raise AttributeError, name

    def tostring(self):
        return str(self)
    
    def search(self, seq, pos=0, endpos=None):
        m = self.grouped_re.search(buffer(seq.data), pos, endpos)
        if m is None:
            return None
        return PrositeMatch(self, seq, m)
    def match(self, seq, pos=0, endpos=None):
        m = self.grouped_re.match(buffer(seq.data), pos, endpos)
        if m is None:
            return None
        return PrositeMatch(self, seq, m)

    # I was thinking about adding sub, subn, findall, etc., but either
    # you just want the string (in which case, use the ".re") or
    # you could be changing to a different alphabet (eg, T->U).


# Elements of a Prosite pattern
class PrositeTerm:
    def __init__(self, letters, ignore, is_begin, is_end, \
                 min_count, max_count):
        self.letters = letters
        self.ignore = ignore
        self.is_begin = is_begin
        self.is_end = is_end
        self.min_count = min_count
        self.max_count = max_count
    def copy(self):
        return PrositeTerm(self.letters, self.ignore, self.is_begin,
                           self.is_end, self.min_count, self.max_count)
    def __str__(self):
        # Convert the term back into Prosite form
        if self.is_begin:
            s = "<"
        else:
            s = ""
        if self.ignore:
            s = s + "{" + self.letters + "}"
        elif len(self.letters) == 1:
            s = s + self.letters
        else:
            s = s + "[" + self.letters + "]"

        if self.min_count == self.max_count:
            if self.min_count == 1:
                pass
            else:
                s = s + "(%d)" % self.min_count
        else:
            s = s + "(%d,%d)" % (self.min_count, self.max_count)
        if self.is_end:
            s = s + ">"
        return s
        
    def no_count_str(self):
        # Convert the term back into Prosite form, without the repeat
        # count fields.
        
        if self.is_begin:
            s = "<"
        else:
            s = ""
        if self.ignore:
            s = s + "{" + self.letters + "}"
        elif len(self.letters) == 1:
            s = s + self.letters
        else:
            s = s + "[" + self.letters + "]"
        return s

# Results of a Prosite match.  Wrapper to the re.MatchObj, but returns
# Seq objects instead of strings.  And lookee - it implements the Seq
# interface too!
class PrositeMatch:
    def __init__(self, prosite, seq, match):
        self.prosite = prosite
        self.seq = seq
        self.match = match
        self.pos = match.pos
        self.endpos = match.pos

        # for Seq.Seq initialization
        self.data = match.group(0)
        self.alphabet = seq.alphabet

    def __repr__(self):
        # XXX this isn't the right way
        return "<PrositeMatch instance at %x>" % id(self)
    def __str__(self):
        return str(self.data)
    def __len__(self): return len(self.data)
    def __getitem__(self, i): return self.data[i]
    def __getslice__(self, i, j):
        i = max(i, 0); j = max(j, 0)
        return Seq(self.data[i:j], self.alphabet)
    
    def mapping(self):
        """return a list of numbers mapping to items of the original pattern

        For example, if the Prosite pattern is "[AP](2)-D." matched against
        "PAD", then the mapping is [1, 1, 2], meaning the first character
        of the match ("P") is from the first Prosite group ("[AP]"), as
        is the second letter ("A").  The 3rd letter ("D") is mapped to
        group 2 of the pattern.
        """

        vals = []
        i = 0
        start = self.start(0)
        try:
            while 1:
                end = self.match.end(i+1)
                while start < end:
                    vals.append(i)
                    start = start + 1
                i = i + 1
        except IndexError:
            pass
        return vals

    def mapped_pattern(self):
        """returns the specific Prosite pattern used to find this sequence

        >>> p = Prosite.compile("[AP](2,3)-D.")
        >>> m = p.search(Seq.Seq("PAD"))
        >>> mapping = m.mapping()
        >>> mapped = m.mapped_pattern()
        >>> print str(m[1]), str(p[mapping[1]]), str(mapped[1])
        P [AP](2,3) [AP]
        >>> print str(mapped)
        [AP]-[AP]-D.
        >>> 

        Note that the original term includes the count, while the
        mapped pattern does the expansion.
        
        """
        return pattern_mapping(self.prosite, self.mapping())

    def start(self, g=0):
        return self.match.start(g)
    def end(self, g=0):
        return self.match.end(g)
    def span(self, g):
        return self.match.span(g)
    def groups(self, default=None):
        result = []
        alphabet = self.alphabet
        for g in self.match.groups(default):
            result.append( Seq.Seq(g, alphabet) )
        return tuple(result)
    def group(self, *groups):
        result = apply(self.match.group, groups)
        if result == ():
            return result
        if len(result) == 1:
            return Seq.Seq(result, self.alphabet)
        retval = []
        for x in result:
            retval.append(Seq.Seq(x, self.alphabet))
        return tuple(retval)

# Use to parse a single element of a pattern
prosite_term_re = re.compile(r"""
(?:
  ([A-Zx])|            # a character OR
  (\[[A-Z]+\])|        # something in []s OR
  (\{[A-Z]+\})         # something in {}s
)(?:\((\d+)(,\d+)?\))? # optional count of the form "(i,j)", ",j" optional
$
""", re.VERBOSE)


def pattern_mapping(prosite, mapping):
    data = []
    for i in mapping:
        x = prosite[i].copy()
        x.min_count = x.max_count = 1
        data.append(x)
    return Prosite(data=data)
    

# This does not verify the pattern is correct - invalid patterns can
# be converted!
def find_terms(pattern):
    if pattern[-1:] != ".":
        raise TypeError, "not a prosite pattern - needs a final '.'"
    pattern = pattern[:-1]
    terms = string.split(pattern, "-")
    result = []
    for term in terms:
        # Starts with a "<"?
        if term[:1] == "<":
            term = term[1:]
            is_begin = 1
        else:
            is_begin = 0
        
        # Ends with a ">"?
        if term[-1:] == ">":
            term = term[:-1]
            is_end = 1
        else:
            is_end = 0

        # Get the elements of the term
        match = prosite_term_re.match(term)
        if match is None:
            raise TypeError, "not a Prosite term (%s)" % repr(term)

        if match.group(1) is not None:
            # Single letter
            ignore = 0
            letters = match.group(1)
        elif match.group(2) is not None:
            # Letters inside of "[]"s
            ignore = 0
            letters = match.group(2)[1:-1]
        elif match.group(3) is not None:
            # Letters inside of "{}"s
            ignore = 1
            letters = match.group(3)[1:-1]
        else:
            raise TypeError, "not a prosite pattern - unknown group?"

        if match.group(4) is not None:
            # there is a minimum number
            min_count = int(match.group(4))
        else:
            # no min, so it's 1
            min_count = 1
        if match.group(5) is not None:
            # there is a maximum number
            max_count = int(match.group(5)[1:])
        else:
            # no max specified, so use the same as the min
            max_count = min_count

        result.append(PrositeTerm(letters, ignore, is_begin,
                                  is_end, min_count, max_count))
    return result


# This does not verify that the pattern is correct - invalid patterns
# can be converted!
def prosite_to_re(pattern):
    """convert a valid Prosite pattern into an re string"""
    s = string.replace(pattern, "{", "[^")
    return string.translate(s, _prosite_trans, "-.")

# This does not verify the pattern is correct - invalid patterns can
# be converted!
def prosite_to_grouped_re(pattern):
    """convert a valid Prosite pattern into an re with groups for each term"""
    s = string.replace(pattern, "{", "[^")
    s = string.translate(s, _prosite_trans, ".")
    if s[:1] == "^":
        s = "^(" + s[1:]
    else:
        s = "(" + s
    if s[-1:] == "$":
        s = s[:-1] + ")$"
    else:
        s = s + ")"
    return string.replace(s, "-", ")(")


prosite_re = re.compile(r"""
^<?                   # starts with an optional "<"
(
  [A-Zx]|             # a character OR
  \[[A-Z]+\]|         # something in []s OR
  \{[A-Z]+\}          # something in {}s
)(\(\d+(,\d+)?\))?    # optional count of the form "(i,j)" (",j" is optional)
(-                    # new terms seperated by a '-'
 (
  [A-Zx]|             # a character OR
  \[[A-Z]+\]|         # something in []s OR
  \{[A-Z]+\}          # something in {}s
 )(\(\d+(,\d+)?\))?   # optional count
)*                    # repeat until done
>?                    # pattern ends with an optional ">"
\.$                   # description ends with a required "."
""", re.VERBOSE)

# This verifies the pattern is correct.
def verify_pattern(pattern):
    """returns 1 if the Prosite pattern is syntactically correct, else 0"""
    return prosite_re.match(pattern) is not None

def _verify_test(infile):
    """verify the patterns from a Prosite file handle"""
    pattern = ""
    while 1:
        line = infile.readline()
        if not line:
            break
        if line[:2] != "PA":
            continue
    
        pattern = pattern + line[5:-1]
        if line[-2] == ".":
            try:
                print "*" * 60
                print pattern
                pattern = fix_errors.get(pattern, pattern)
                p = compile(pattern)
                print prosite_to_re(pattern)
                print repr(p.re)
                print prosite_to_grouped_re(pattern)
                print repr(p.grouped_re)
                terms = str(p)
                if terms != pattern:
                    print "DIFFER", terms, pattern
            except TypeError, msg:
                print "PROBLEM", pattern, msg
            pattern = ""

# Commented out by jchang 4/13/00.
# Specific to Andrew's test environment.
#if __name__ == "__main__":
#    import os
#    infile = os.popen("bzcat /home/dalke/ftps/prosite/prosite.dat.bz2 | grep ^PA")
#    _verify_test(infile)

