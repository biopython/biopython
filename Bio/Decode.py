# Copyright 2002 by Andrew Dalke.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Decode elements from a Std/Martel parsed XML stream (OBSOLETE).

Andrew Dalke is no longer maintaining Martel or Bio.Mindy, and these modules
(and therefore Bio.Decode) have been deprecated.  They are no longer used in
any of the current Biopython parsers, and are likely to be removed in a
future release."""

import warnings
warnings.warn("Martel and those parts of Biopython depending on it" \
              +" directly (such as Bio.Mindy and Bio.Decode) are now" \
              +" deprecated, and will be removed in a future release of"\
              +" Biopython.  If you want to continue to use this code,"\
              +" please get in contact with the Biopython developers via"\
              +" the mailing lists to avoid its permanent removal from"\
              +" Biopython.", \
              DeprecationWarning)

import string
from Bio.Parsers.spark import GenericScanner, GenericParser

def unescape_C(s):
    result = []
    for i in range(len(s)):
        if s[i] != "\\":
            result.append(s[i])
            continue
        c = s[i+1:i+2]
        if c == "x":
            x = s[i+2:i+4]
            if len(x) != 2:
                raise ValueError("invalid \\x escape")
            i = int(x, 16)
            result.append(chr(i))
            continue
        if c in "01234567":
            x = s[i+1:i+4]
            # \octals don't do a length assertion check
            i = int(x, 8)
            result.append(chr(i))
            continue
        result.append(c)
    return "".join(result)

def join_english(fields):
    if not fields:
        return ""
    s = fields[0]
    for field in fields[1:]:
        if s[-1:] == "-" and s[-3:-2] == "-":
            s = s + field
            continue
        if s.find(" ") == -1 and field.find(" ") == -1:
            s = s + field
            continue
        s = s + " " + field
    return (" ".join(s.split())).strip()



def chomp(s, c):
    if s[-1:] == c:
        return s[:-1]
    return s

def lchomp(s, c):
    if s[:1] == c:
        return s[1:]
    return s
    
def chompchomp(s, c):
    if s[:1] == c and s[-1:] == c:
        return s[1:-1]
    return s

def fixspaces(s):
    # s.split breaks down to a list of words
    # " ".join puts them together
    # strip removes leading and trailing spaces
    return " ".join(s.split()).strip()

def join_fixspaces(lines):
    return " ".join((" ".join(lines)).split()).strip()

def tr(s, frm, to):
    table = string.maketrans(frm, to)
    return s.translate(table)

def safe_int(s):
    """converts to int if the number is small, long if it's large"""
    try:
        return int(s)
    except ValueError:
        return long(s)

decode_functions = {
    "chomp": (chomp, str, str),
    "chompchomp": (chompchomp, str, str),
    "chop": (lambda s: s[:-1], str, str),
    "chopchop": (lambda s: s[1:-1], str, str),
    "fixspaces": (fixspaces, str, str),
    "lchomp": (lchomp, str, str),
    "lchop": (lambda s: s[1:], str, str),
    "lower": (lambda s: s.lower(), str, str),
    "lstrip": (lambda s: s.lstrip(), str, str),
    "replace": (lambda s, old, new: s.replace(old, new), str, str),
    "rstrip": (lambda s: s.rstrip(), str, str),
    "str": (str, str, str),
    "strip": (lambda s: s.strip(), str, str),
    "tr": (tr, str, str),
    "unescape.c": (unescape_C, str, str),
    "unescape.doublequote": (lambda s: s.replace('""', '"'), str, str),
    "unescape.singlequote": (lambda s: s.replace("''", "'"), str, str),
    "upper": (lambda s: s.upper(), str, str),

    # List operations
    "join": (lambda lst, s = " ": s.join(lst), list, str),
    "join.english": (join_english, list, str),

    # Integer operations
    "int": (safe_int, [float, str, int], int),
    "int.comma": (lambda s: safe_int(s.replace(",", "")),
                  [float, str, int], int),
    "hex": (hex, str, int),
    "oct": (oct, str, int),
    "add": ((lambda i, j: i+j), int, int),

    # Float operations
    "float": (float, (float, str, int), float),
    
    }

def _fixup_defs():
    # Normalize so the 2nd and 3rd terms are tuples
    for k, v in decode_functions.items():
        f, in_types, out_types = v
        if isinstance(in_types, type([])):
            in_types = tuple(in_types)
        elif not isinstance(in_types, type( () )):
            in_types = (in_types,)

        if isinstance(out_types, type([])):
            out_types = tuple(out_types)
        elif not isinstance(out_types, type( () )):
            out_types = (out_types,)

        decode_functions[k] = (f, in_types, out_types)
_fixup_defs()

class Token:
    def __init__(self, type):
        self.type = type
    def __cmp__(self, other):
        return cmp(self.type, other)
    def __repr__(self):
        return "Token(%r)" % (self.type,)

class ValueToken(Token):
    def __init__(self, type, val):
        Token.__init__(self, type)
        self.val = val
    def __cmp__(self, other):
        return cmp(self.type, other)
    def __repr__(self):
        return "%s(%r)" % (self.__class__.__name__, self.val)
    def __str__(self):
        return str(self.val)

class Integer(ValueToken):
    def __init__(self, val):
        ValueToken.__init__(self, "integer", val)

class Float(ValueToken):
    def __init__(self, val):
        ValueToken.__init__(self, "float", val)

class String(ValueToken):
    def __init__(self, val):
        ValueToken.__init__(self, "string", val)

class FunctionName(ValueToken):
    def __init__(self, val):
        ValueToken.__init__(self, "functionname", val)

class DecodeScanner(GenericScanner):
    def __init__(self):
        GenericScanner.__init__(self)
 
    def tokenize(self, input):
        self.rv = []
        GenericScanner.tokenize(self, input)
        return self.rv

    def t_functionname(self, input):
        r" \w+(\.\w+)*"
        self.rv.append(FunctionName(input))

    def t_pipe(self, input):
        r" \| "
        self.rv.append(Token("pipe"))
        
    def t_open_paren(self, input):
        r" \( "
        self.rv.append(Token("open_paren"))

    def t_close_paren(self, input):
        r" \) "
        self.rv.append(Token("close_paren"))

    def t_comma(self, input):
        r" , "
        self.rv.append(Token("comma"))

    def t_whitespace(self, input):
        r" \s+ "
        pass

    def t_string(self, input):
        r""" "([^"\\]+|\\.)*"|'([^'\\]+|\\.)*' """
        # "'  # emacs cruft
        s = input[1:-1]
        s = unescape_C(s)
        
        self.rv.append(String(s))

    def t_float(self, input):
        r""" [+-]?((\d+(\.\d*)?)|\.\d+)([eE][+-]?[0-9]+)? """
        # See if this is an integer
        try:
            self.rv.append(Integer(safe_int(input)))
        except ValueError:
            self.rv.append(Float(float(input)))

class Function:
    def __init__(self, name, args = ()):
        self.name = name
        self.args = args
    def __str__(self):
        args = self.args
        if not args:
            s = ""
        else:
            s = str(args)[1:-1]
        return "%s(x, %s)" % (self.name, s)
    __repr__ = __str__

class DecodeParser(GenericParser):
    def __init__(self, start = "expression"):
        GenericParser.__init__(self, start)
        self.begin_pos = 0

    def p_expression(self, args):
        """
        expression ::= term
        expression ::= term pipe expression
        """
        if len(args) == 1:
            return [args[0]]
        return [args[0]] + args[2]

    def p_term(self, args):
        """
        term ::= functionname
        term ::= functionname open_paren args close_paren
        """
        if len(args) == 1:
            return Function(args[0].val)
        return Function(args[0].val, tuple([x.val for x in args[2]]))

    def p_args(self, args):
        """
        args ::= arg
        args ::= arg comma args
        """
        if len(args) == 1:
            return [args[0]]
        return [args[0]] + args[2]

    def p_arg(self, args):
        """
        arg ::= string
        arg ::= integer
        arg ::= float
        """
        return args[0]
    
def scan(input):
    scanner = DecodeScanner()
    return scanner.tokenize(input)

def parse(tokens):
    parser = DecodeParser()
    return parser.parse(tokens)

_decoder_cache = {}

class FunctionCall:
    def __init__(self, f, args):
        self.f = f
        self.args = args
    def __call__(self, x):
        return self.f(x, *self.args)

class FunctionCallChain:
    def __init__(self, inner_f, f, args):
        self.inner_f = inner_f
        self.f = f
        self.args = args
    def __call__(self, x):
        return self.f(self.inner_f(x), *self.args)

#### I don't think this is the right way to do things
##class CheckTypes:
##    def __init__(self, f, call_types, return_types):
##        self.f = f
##        self.call_types = call_types
##        self.return_types = return_types
##    def __call__(self, x):
##        if self.call_types is not None:
##            for T in self.call_types:
##                if isinstance(x, T):
##                    break
##            else:
##                raise TypeError(
##                    "Call value %s of type %s, expecting one of %s" %
##                    (x, type(x).__name__,
##                     [T.name for T in self.call_types]))
##        y = self.f(x)

##        if not self.return_types:
##            return y
        
##        for T in self.return_types:
##            if isinstance(y, T):
##                return y
##        raise TypeError("Return value %s of type %s, expecting one of %s" %
##                        (y, type(y).__name__,
##                         [T.name for T in self.return_types]))

def make_decoder(s):
    try:
        return _decoder_cache[s]
    except KeyError:
        pass
    
    functions = parse(scan(s))
    
    f = functions[0]
    fc = decode_functions[f.name][0]
    args = f.args
    if args:
        fc = FunctionCall(fc, args)
    for f in functions[1:]:
        fc = FunctionCallChain(fc, decode_functions[f.name][0], f.args)
    _decoder_cache[s] = fc
    return fc

def _verify_subtypes(subset, total, old_name, new_name):
    for x in subset:
        if x not in total:
            raise TypeError("%s can produce a %r value not accepted by %s" %
                            (old_name, x.__name__, new_name))

_typechecked_decoder_cache = {}
def make_typechecked_decoder(s, input_types = None, output_types = None):
    cache_lookup = (s, input_types, output_types)
    try:
        return _typechecked_decoder_cache[cache_lookup]
    except KeyError:
        pass
    if input_types is not None and not isinstance(input_types, type( () )):
        input_types = (input_types,)
    if output_types is not None and not isinstance(output_types, type( () )):
        output_types = (output_types,)

    functions = parse(scan(s))

    # Make sure the input type(s) are allowed
    f = functions[0]
    fc, in_types, out_types = decode_functions[f.name]
    if input_types is not None:
        for x in input_types:
            if x not in in_types:
                raise TypeError(
                    "the input type includes %r which isn't supported by %s" %
                    (x.__name__, f.name))

    # Do the composition
    old_name = f.name
    input_types = out_types
    args = functions[0].args
    if args:
        fc = FunctionCall(fc, args)
    
    for f in functions[1:]:
        transform_func, in_types, out_types = decode_functions[f.name]
        _verify_subtypes(input_types, in_types, old_name, f.name)
        old_name = f.name
        input_types = out_types
        fc = FunctionCallChain(fc, transform_func, f.args)

    if output_types is not None:
        _verify_subtypes(input_types, output_types, old_name, "the output")
    _typechecked_decoder_cache[cache_lookup] = fc
    return fc
    

def test():
    assert make_decoder("chop")("Andrew") == "Andre"
    assert make_decoder("int")("9") == 9
    assert make_decoder('join(" ")')(["Andrew", "Dalke"]) == \
                                          "Andrew Dalke"
    assert make_decoder('chomp("|")')("|test|") == "|test"
    assert make_decoder('chomp("|")')("|test") == "|test"
    assert make_decoder('chomp("A")|chop')("BA") == ""
    assert make_decoder('chomp("A")|chop')("AB") == "A"
    assert make_decoder('chop|chomp("A")')("AB") == ""
    assert make_decoder('chop|chomp("A")')("BA") == "B"
    assert make_decoder('add(5)')(2) == 7
    assert make_decoder('add(-2)')(5) == 3
    
if __name__ == "__main__":
    test()
