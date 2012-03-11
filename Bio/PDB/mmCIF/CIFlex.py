#!/usr/bin/python

import sys  # to get args
import time # to time function
import re  # to supply RE args to lexer
import ply.lex as lex  # lexer
from ply.lex import TOKEN  # to assign complex docstrings to tokens

class CIFlex:
    ### Token source:
    # http://www.iucr.org/resources/cif/spec/version1.1/cifsyntax
    
    ### Single characters
    double_quote = r'"'
    single_quote = r"'"
    pound = r"\#"
    dollar = r"$"
    underscore = r"_"
    semi = r";"
    lbracket = r"["
    rbracket = r"]"
    space = r" "
    tab = r"\t"
    
    ### Character Sets
    # re character class rules
    # for using unescaped special characters:
    # rbracket must be first
    # hyphen must be first or last (we'll assume last)
    # caret must not be first
    # backslash must be escaped
    
    #<OrdinaryChar>  :  { '!' | '%' | '&' | '(' | ')' | '*' | '+' | ',' | '-' | '.' | '/' | '0' | '1' | '2' | '3' | '4' | '5' | '6' | '7' | '8' | '9' | ':' | '<' | '=' | '>' | '?' | '@' | 'A' | 'B' | 'C' | 'D' | 'E' | 'F' | 'G' | 'H' | 'I' | 'J' | 'K' | 'L' | 'M' | 'N' | 'O' | 'P' | 'Q' | 'R' | 'S' | 'T' | 'U' | 'V' | 'W' | 'X' | 'Y' | 'Z' | '\' | '^' | '`' | 'a' | 'b' | 'c' | 'd' | 'e' | 'f' | 'g' | 'h' | 'i' | 'j' | 'k' | 'l' | 'm' | 'n' | 'o' | 'p' | 'q' | 'r' | 's' | 't' | 'u' | 'v' | 'w' | 'x' | 'y' | 'z' | '{' | '|' | '}' | '~' }
    _ordinary_char = r"!%&()*+,./0-9:<=>?@A-Z`a-z{|}~^\\-"
    ordinary_char = r"[" + _ordinary_char + r"]"
    semi_ordinary_char = r"[" + semi + _ordinary_char + r"]"
    
    #<NonBlankChar>  :  <OrdinaryChar> | <double_quote> | '#' | '$' | <single_quote> | '_' |';' | '[' | ']'
    _non_blank_char = rbracket + double_quote + pound + dollar + single_quote + underscore + semi + lbracket + _ordinary_char
    non_blank_char = r"[" + _non_blank_char + r"]"
    
    #<TextLeadChar>  :  <OrdinaryChar> | <double_quote> | '#' | '$' | <single_quote> | '_' | <SP> | <HT> |'[' | ']'
    _text_lead_char = rbracket + double_quote + pound + dollar + single_quote + underscore + space + tab + lbracket + _ordinary_char
    text_lead_char = r"[" + _text_lead_char + r"]"
    
    #<AnyPrintChar>  :  <OrdinaryChar> | <double_quote> | '#' | '$' | <single_quote> | '_' | <SP> | <HT> | ';' | '[' | ']'
    _any_print_char = rbracket + double_quote + pound + dollar + single_quote + underscore + space + tab + semi + lbracket + _ordinary_char
    any_print_char = r"[" + _any_print_char + r"]"

    ### WhiteSpace and Comments
    # (CR or CRLF or LF) hopefully will handle any platform's EOL
    _eol = r"\r\n|\n|\r"
    eol = r"(" + _eol + r")"
    noteol = r"[^" + eol + r"]"
    #<Comments>  :  { '#' {<AnyPrintChar>}* <eol>}+
    comments = r"(" + pound + any_print_char + r"*" + eol + r")+"
    #<TokenizedComments>  :     { <SP> | <HT> | <eol> |}+ <Comments>
    tokenized_comments = r"[" + space + tab + eol + r"]+" + comments
    #<WhiteSpace>  :    { <SP> | <HT> | <eol> | <TokenizedComments>}+
    #whitespace = r"(" + space + r"|" + tab +  r"|" + eol +  r"|" + tokenized_comments + r")+" 
    _whitespace = r" \t" + _eol
    whitespace = r"[" + _whitespace + r"]" 

    ### Character Strings and Text Fields
    #<CharString>  :    <UnquotedString> | <SingleQuotedString> | <DoubleQuotedString>
    # XXX yacc will do this XXX 
    #<eol><UnquotedString>  :   <eol><OrdinaryChar> {<NonBlankChar>}*
    eol_unquoted_string = r"^" + ordinary_char + non_blank_char + r"*"
    #<noteol><UnquotedString>  :    <noteol>{<OrdinaryChar>|';'} {<NonBlankChar>}*
    ### XXX this matches only 2+ char strings because 
    ### the first char matches noteol
    ### and the second char matches the required semi_ordinary_char
    # noteol_unquoted_string = noteol + semi_ordinary_char + non_blank_char + r"*"
    #noteol_unquoted_string = r"[ \t](" + semi_ordinary_char + non_blank_char + r"*)"
    ### XXX provided this is the last type of match, should work
    ### eol_unquoted_string and semi_text_field 
    ### together should enforce the semicolon rule
    noteol_unquoted_string = non_blank_char + r"+"
    #<SingleQuotedString><WhiteSpace>  :    <single_quote>{<AnyPrintChar>}* <single_quote> <WhiteSpace>
    single_quoted_string = single_quote + any_print_char + r"*" + single_quote + whitespace
    #<DoubleQuotedString><WhiteSpace>  :    <double_quote> {<AnyPrintChar>}* <double_quote> <WhiteSpace>
    double_quoted_string = double_quote + any_print_char + r"*" + double_quote + whitespace 

    #<TextField>  :     { <SemiColonTextField> }
    ## XXX yacc would do this; redundant XXX 
    #<eol><SemiColonTextField> : <eol>';' { {<AnyPrintChar>}* <eol>
    #                            {{<TextLeadChar> {<AnyPrintChar>}*}? <eol>}*
    #                            } ';'
    _semi_text_header = r"^" + semi + any_print_char + r"*" + eol
    _semi_text_line = r"(" + text_lead_char + any_print_char + r"*)?" + eol
    _semi_text_end = r"^" + semi + whitespace
    semi_text_field = _semi_text_header + r"(" + _semi_text_line + r")*" + _semi_text_end
    
    ### Numeric Values
    #<Number>  :    {<Integer> | <Float> }
    #<Numeric>  :   { <Number> | <Number> '(' <UnsignedInteger> ')' }
    ## XXX yacc will do this; not sure what second numeric type is XXX
    #<Digit>  :     { '0' | '1' | '2' | '3' | '4' | '5' | '6' | '7' | '8' | '9' }
    #<UnsignedInteger>  :   { <Digit> }+
    unsigned_integer = r"\d+"
    #<Integer>  :   { '+' | '-' }? <UnsignedInteger>
    integer = r"[+-]?" + unsigned_integer
    #<Exponent>  :  { {'e' | 'E' } | {'e' | 'E' } { '+' | '- ' } } <UnsignedInteger>
    exponent = r"e{1,2}[+-]?" + unsigned_integer
    #<Float> : { <Integer><Exponent> |
    #          { {'+'|'-'} ? { {<Digit>} * '.' <UnsignedInteger> } |
    #          { <Digit>} + '.' } } {<Exponent>} ? } }
    mantissa = r"[+-]?(\d*\.\d+|\d+\.?)"
    float_type = mantissa + r"(" + exponent + r")?"

    ### Tags and Values
    #<Value>  :     { '.' | '?' | <Numeric> | <CharString> | <TextField> }
    # XXX yacc will do this XXX
    #<Tag>  :   '_'{ <NonBlankChar>}+
    tag = underscore + non_blank_char + r"+"

    ### Reserved Words
    #<DATA_>  :     {'D' | 'd'} {'A' | 'a'} {'T' | 't'} {'A' | 'a'} '_'
    data = r"data_"
    #<LOOP_>  :     {'L' | 'l'} {'O' | 'o'} {'O' | 'o'} {'P' | 'p'} '_'
    loop = r"loop_"
    #<GLOBAL_>  :   {'G' | 'g'} {'L' | 'l'} {'O' | 'o'} {'B' | 'b'} {'A' | 'a'} {'L' | 'l'} '_'
    global_type = r"global_"
    #<SAVE_>  :     {'S' | 's'} {'A' | 'a'} {'V' | 'v'} {'E' | 'e'} '_'
    save = r"save_"
    #<STOP_>  :     {'S' | 's'} {'T' | 't'} {'O' | 'o'} {'P' | 'p'}'_'
    stop = r"stop_"

    ### Basic Structure of a CIF
    # XXX yacc will do all of this XXX
    #<CIF>  :   <Comments>? <WhiteSpace>? { <DataBlock> { <WhiteSpace> <DataBlock> }* { <WhiteSpace> }? }?
    #<DataBlock>  :     <DataBlockHeading> {<WhiteSpace> { <DataItems> | <SaveFrame>} }*
    #<DataBlockHeading>  :  <DATA_> { <NonBlankChar> }+
    #<SaveFrame>  :     <SaveFrameHeading> { <WhiteSpace> <DataItems> }+ <WhiteSpace> <SAVE_>
    #<SaveFrameHeading>  :  <SAVE_> { <NonBlankChar> }+
    #<DataItems>  :     <Tag> <WhiteSpace> <Value> |
    #                   <LoopHeader><LoopBody>
    #<LoopHeader>  :    <LOOP_> {<WhiteSpace> <Tag>}+
    #<LoopBody>  :  <Value> { <WhiteSpace> <Value> }*

    tokens = (
        "COMMENTS",
        # Reserved words:
        "DATA",
        "LOOP",
        "GLOBAL",
        "SAVE",
        "STOP",
        # Contents:
        "TAG",
        ### VALUE:
        "INAPPLICABLE",
        "UNKNOWN",
        "EOL_UNQUOTED_STRING",
        "NOTEOL_UNQUOTED_STRING",
        "DOUBLE_QUOTED_STRING",
        "SINGLE_QUOTED_STRING",
        "SEMI_TEXT_FIELD",
        "INTEGER",
        "FLOAT",
    )
    
    states = (
        ("data", "exclusive"),
        ("loop", "exclusive"),
        ("semi", "exclusive"),
    )
        
    @TOKEN(comments)
    def t_ANY_COMMENTS(self,t):
        t.lexer.lineno += t.value.count('\n')
        if t.lexer.current_state() == "loop":
            t.lexer.pop_state()
        return None
        
    @TOKEN(data)
    def t_DATA(self,t):
        t.value = t.value[5:]
        t.lexer.push_state("data")
        return t
    
    @TOKEN(loop)
    def t_data_LOOP(self,t):
        if t.lexer.current_state() == "loop":
            warnings.warn("ERROR: Illegal nested loop.", RuntimeWarning)
        t.lexer.push_state("loop")
        return t
    
    @TOKEN(global_type)
    def t_data_GLOBAL(self,t):
        return t
    
    @TOKEN(save)
    def t_data_SAVE(self,t):
        warnings.warn("ERROR: found save frame, this parser is not intended to handle dictionaries.", RuntimeWarning)
        return t
    
    @TOKEN(stop)
    def t_data_STOP(self,t):
        return t
    
    @TOKEN(tag)
    def t_data_loop_TAG(self,t):
        return t
    
    @TOKEN(semi_text_field)
    def t_data_SEMI_TEXT_FIELD(self,t):
        # remove \n by splitting into lines and joining
        t.value = "".join(t.value.splitlines())
        return t

    @TOKEN(double_quoted_string)
    def t_data_loop_DOUBLE_QUOTED_STRING(self,t):
        return t
    
    @TOKEN(single_quoted_string)
    def t_data_loop_SINGLE_QUOTED_STRING(self,t):
        return t

    @TOKEN(integer)
    def t_data_loop_INTEGER(self,t):
        return t

    @TOKEN(float_type)
    def t_data_loop_FLOAT(self,t):
        return t

    def t_data_loop_INAPPLICABLE(self,t):
        r"\."
        return t
         
    def t_data_loop_UNKNOWN(self,t):
        r"\?"
        return t
        
    @TOKEN(eol_unquoted_string)
    def t_data_loop_EOL_UNQUOTED_STRING(self,t):
        return t

    @TOKEN(noteol_unquoted_string)
    def t_data_loop_NOTEOL_UNQUOTED_STRING(self,t):
        return t

    # Ignored characters: spaces and tabs
    t_ANY_ignore  = ' \t'

    # Newline rule to track line numbers
    def t_ANY_newline(self,t):
        r'\n+'
        t.lexer.lineno += len(t.value)

    # Error handling rule
    def t_ANY_error(self,t):
        print "Illegal character '%s'" % t.value[0]
        t.lexer.skipped_lines += t.value.count('\n')
        t.lexer.skip(1)

    ### Public methods 

    def build(self, **kwargs):
        self._lexstart = time.clock()
        # set re.MULTILINE while preserving any user reflags
        re_old = 0
        if "reflags" in kwargs.keys():
            re_old = kwargs["reflags"]
        kwargs["reflags"] = re_old | re.MULTILINE | re.I

        self.lexer = lex.lex(module=self,**kwargs)
        # store number of skipped lines
        self.lexer.skipped_lines = 0
        self.lex_init = time.clock() - self._lexstart

    def test(self, data):
        self.lexer.input(data)
        while True:
            token = self.lexer.token()
            if not token:
                break
            #print token 
        if self.lexer.skipped_lines:
            print "Skipped %s lines" % self.lexer.skipped_lines
        self.lex_end = time.clock() - self._lexstart
        print "Runtime: ", self.lex_end 

if len(sys.argv) == 2:
    filename = sys.argv[1]

    with open(filename) as fh:
        # start new lexer class instance
        m = CIFlex()  
        # build lexer with optional flags (i.e. debug=1 or optimize=1)
        m.build(debug=1)
        # lex file
        m.test(fh.read())


# vim:sw=4:ts=4:expandtab
