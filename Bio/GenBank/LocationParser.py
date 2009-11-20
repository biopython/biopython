# First pass at a parser for the location fields of a feature table.
# Everything likely to change.
#
# This does NOT cope with the Gap(), Gap(X), or Gap(unkXXX) tokens used
# in CONTIG lines, which are otherwise similar to feature locations.
#
# Based on the DDBJ/EMBL/GenBank Feature Table Definition Version 2.2
# Dec 15 1999 available from EBI, but the documentation is not
# completely internally consistent much less agree with real-life
# examples.  Conflicts resolved to agree with real examples.
#
# This does NOT cope with the Gap(), Gap(X), or Gap(unkXXX) tokens used
# in CONTIG lines, which are otherwise similar to feature locations.
#
# Uses John Aycock's SPARK for parsing
from Bio.Parsers.spark import GenericScanner, GenericParser

class Token:
    def __init__(self, type):
        self.type = type
    def __cmp__(self, other):
        return cmp(self.type, other)
    def __repr__(self):
        return "Tokens(%r)" % (self.type,)

# "38"
class Integer:
    type = "integer"
    def __init__(self, val):
        self.val = val
    def __cmp__(self, other):
        return cmp(self.type, other)
    def __str__(self):
        return str(self.val)
    def __repr__(self):
        return "Integer(%s)" % self.val

# From the BNF definition, this isn't needed.  Does tht mean
# that bases can be refered to with negative numbers?
class UnsignedInteger(Integer):
    type = "unsigned_integer"
    def __repr__(self):
        return "UnsignedInteger(%s)" % self.val

class Symbol:
    type = "symbol"
    def __init__(self, name):
        self.name = name
    def __cmp__(self, other):
        return cmp(self.type, other)
    def __str__(self):
        return str(self.name)
    def __repr__(self):
        return "Symbol(%s)" % repr(self.name)

# ">38"  -- The BNF says ">" is for the lower bound.. seems wrong to me
class LowBound:
    def __init__(self, base):
        self.base = base
    def __repr__(self):
        return "LowBound(%r)" % self.base

# "<38"
class HighBound:
    def __init__(self, base):
        self.base = base
    def __repr__(self):
        return "HighBound(%r)" % self.base

# 12.34
class TwoBound:
    def __init__(self, low, high):
        self.low = low
        self.high = high
    def __repr__(self):
        return "TwoBound(%r, %r)" % (self.low, self.high)

# 12^34
class Between:
    def __init__(self, low, high):
        self.low = low
        self.high = high
    def __repr__(self):
        return "Between(%r, %r)" % (self.low, self.high)

# 12..34
class Range:
    def __init__(self, low, high):
        self.low = low
        self.high = high
    def __repr__(self):
        return "Range(%r, %r)" % (self.low, self.high)

class Function:
    def __init__(self, name, args):
        self.name = name
        self.args = args
    def __repr__(self):
        return "Function(%r, %r)" % (self.name, self.args)

class AbsoluteLocation:
    def __init__(self, path, local_location):
        self.path = path
        self.local_location = local_location
    def __repr__(self):
        return "AbsoluteLocation(%r, %r)" % (self.path, self.local_location)

class Path:
    def __init__(self, database, accession):
        self.database = database
        self.accession = accession
    def __repr__(self):
        return "Path(%r, %r)" % (self.database, self.accession)

class FeatureName:
    def __init__(self, path, label):
        self.path = path
        self.label = label
    def __repr__(self):
        return "FeatureName(%r, %r)" % (self.path, self.label)
    
class LocationScanner(GenericScanner):
    def __init__(self):
        GenericScanner.__init__(self)

    def tokenize(self, input):
        self.rv = []
        GenericScanner.tokenize(self, input)
        return self.rv

    def t_double_colon(self, input):
        r" :: "
        self.rv.append(Token("double_colon"))
    def t_double_dot(self, input):
        r" \.\. "
        self.rv.append(Token("double_dot"))
    def t_dot(self, input):
        r" \.(?!\.) "
        self.rv.append(Token("dot"))
    def t_caret(self, input):
        r" \^ "
        self.rv.append(Token("caret"))
    def t_comma(self, input):
        r" \, "
        self.rv.append(Token("comma"))
    def t_integer(self, input):
        r" -?[0-9]+ "
        self.rv.append(Integer(int(input)))
    def t_unsigned_integer(self, input):
        r" [0-9]+ "
        self.rv.append(UnsignedInteger(int(input)))
    def t_colon(self, input):
        r" :(?!:) "
        self.rv.append(Token("colon"))
    def t_open_paren(self, input):
        r" \( "
        self.rv.append(Token("open_paren"))
    def t_close_paren(self, input):
        r" \) "
        self.rv.append(Token("close_paren"))
    def t_symbol(self, input):
        r" [A-Za-z0-9_'*-][A-Za-z0-9_'*.-]* "
        # Needed an extra '.'
        self.rv.append(Symbol(input))
    def t_less_than(self, input):
        r" < "
        self.rv.append(Token("less_than"))
    def t_greater_than(self, input):
        r" > "
        self.rv.append(Token("greater_than"))

# punctuation .. hmm, isn't needed for location
#        r''' [ !#$%&'()*+,\-./:;<=>?@\[\\\]^_`{|}~] '''

class LocationParser(GenericParser):
    def __init__(self, start='location'):
        GenericParser.__init__(self, start)
        self.begin_pos = 0

    def p_location(self, args):
        """
        location ::= absolute_location
        location ::= feature_name
        location ::= function
        """
        return args[0]
    
    def p_function(self, args):
        """
        function ::= functional_operator open_paren location_list close_paren
        """
        return Function(args[0].name, args[2])
    
    def p_absolute_location(self, args):
        """
        absolute_location ::= local_location
        absolute_location ::= path colon local_location
        """
        if len(args) == 1:
            return AbsoluteLocation(None, args[-1])
        return AbsoluteLocation(args[0], args[-1])
    
    def p_path(self, args):
        """
        path ::= database double_colon primary_accession
        path ::= primary_accession
        """
        if len(args) == 3:
            return Path(args[0], args[2])
        return Path(None, args[0])
    
    def p_feature_name(self, args):
        """
        feature_name ::= path colon feature_label
        feature_name ::= feature_label
        """
        if len(args) == 3:
            return FeatureName(args[0], args[2])
        return FeatureName(None, args[0])

    def p_feature_label(self, args):
        """
        label ::= symbol
        """
        return args[0].name

    def p_local_location(self, args):
        """
        local_location ::= base_position
        local_location ::= between_position
        local_location ::= base_range
        """
        return args[0]
    def p_location_list(self, args):
        """
        location_list ::= location
        location_list ::= location_list comma location
        """
        if len(args) == 1:
            return args
        return args[0] + [args[2]]

    def p_functional_operator(self, args):
        """
        functional_operator ::= symbol
        """
        return args[0]

    def p_base_position(self, args):
        """
        base_position ::= integer
        base_position ::= low_base_bound
        base_position ::= high_base_bound
        base_position ::= two_base_bound
        """
        return args[0]

    def p_low_base_bound(self, args):
        """
        low_base_bound ::= greater_than integer
        """
        return LowBound(args[1])

    def p_high_base_bound(self, args):
        """
        high_base_bound ::= less_than integer
        """
        return HighBound(args[1])

    def p_two_base_bound_1(self, args):
        """
        two_base_bound ::= open_paren base_position dot base_position close_paren
        """
        # main example doesn't have parens but others do.. (?)
        return TwoBound(args[1], args[3])

    def p_two_base_bound_2(self, args):
        """
        two_base_bound ::= base_position dot base_position
        """
        # two_base_bound with no parentheses like 1.6
        return TwoBound(args[0], args[2])
    
    def p_between_position(self, args):
        """
        between_position ::= base_position caret base_position
        """
        return Between(args[0], args[2])

    def p_base_range(self, args):
        """
        base_range ::= base_position double_dot base_position
        base_range ::= function double_dot base_position
        base_range ::= base_position double_dot function
        base_range ::= function double_dot function
        """
        return Range(args[0], args[2])
        
    def p_database(self, args):
        """
        database ::= symbol
        """
        return args[0].name

    def p_primary_accession(self, args):
        """
        primary_accession ::= symbol
        """
        return args[0].name


_cached_scanner = LocationScanner()
def scan(input):
    """Break a location string into a set of tokens"""
    #scanner = LocationScanner()
    #return scanner.tokenize(input)
    return _cached_scanner.tokenize(input)

_cached_parser = LocationParser()
def parse(tokens):
    """Go from a set of tokens to an object representation"""
    #print "I have", tokens
    #parser = LocationParser()
    #return parser.parse(tokens)
    return _cached_parser.parse(tokens)
