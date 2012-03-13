# Created 2012 
# Lenna X. Peterson
# arklenna@gmail.com

import sys  # to get args
import warnings # to warn
import time # to time function
import re  # to supply RE args to lexer
import ply.lex as lex  # lexer
from ply.lex import TOKEN  # to assign complex docstrings to tokens

class CIFlex:
    # <AnyPrintChar> : [^\r\n]
    # <NonBlankChar> : [^ \t\r\n]
    # <eol> : (\r\n|\n|\r)
    
    tokens = (
        "COMMENT",
        # Reserved words:
        "DATA",
        "LOOP",
        "GLOBAL",
        "SAVE",
        "STOP",
        # Underscore begins tag
        "TAG",
        # Values
        "VALUE",
        "SEMI_HEADER",
        "SEMI_LINES",
        #"SEMI_END",
        #"SEMI_TEXT_FIELD",
        "SEMI_ERROR",
        "UNQ_NUM_STRING",
        "INTEGER",
        "UNQ_STRING",
    )
    
    states = (
        ("data", "exclusive"),
        ("loop", "exclusive"),
    )
    
    # Literal hash, any print chars
    def t_ANY_COMMENT(self,t):
        r"""\#[^\r\n]*"""
        return None
    
    # DATA_, non blank chars
    def t_DATA(self,t):
        r"DATA_[^ \t\r\n]+"
        t.lexer.push_state("data")
        t.value = t.value[5:]
        return t
        
    def t_data_LOOP(self,t):
        r"LOOP_"
        return t

    def t_ANY_GLOBAL(self,t):
        r"GLOBAL_"
        return t
        
    # SAVE_, non blank chars
    def t_ANY_SAVE(self,t):
        r"SAVE_[^ \t\r\n]+"
        t.value = t.value[5:]
        return t
        
    def t_ANY_STOP(self,t):
        r"STOP_"
        return t
    
    # _, non blank chars
    def t_data_TAG(self,t):
        r"_[^ \t\r\n]+"
        return t
    
    # line anchor, semi, any print chars
    _semi_text_header = r"^;[^\r\n]*$"
    # line anchor, any non-semi print char, any print chars
    _semi_text_lines = r"((^[^;\r\n][^\r\n]*)?(\r\n|\r|\n))*;[ \t\n]"
    # line anchor, semi, whitespace
    _semi_text_end = r"^;[ \t\r\n]"
    semi_text_field = _semi_text_header + _semi_text_lines + _semi_text_end

    #def t_data_SEMI_END(self,t):
    #    t.value = t.value.rstrip()
    #    return t
    #t_data_SEMI_END.__doc__ = _semi_text_end

    def t_data_SEMI_HEADER(self,t):
        return t
    t_data_SEMI_HEADER.__doc__ = _semi_text_header

    def t_data_SEMI_LINES(self,t):
        return t
    t_data_SEMI_LINES.__doc__ = _semi_text_lines
    
    #def t_data_SEMI_TEXT_FIELD(self,t):
        ## add \n to count
        #t.lexer.lineno += t.value.count('\n')
        ## remove \n by splitting into lines and joining
        #t.value = "".join(t.value.splitlines())
        #t.type = "VALUE"
        #return t
    #t_data_SEMI_TEXT_FIELD.__doc__ = semi_text_field
    
    # line anchor, semi, non blank chars
    # (that hasn't been recognized in a SEMI_TEXT_FIELD)
    def t_data_SEMI_ERROR(self,t):
        r"^;[^ \t\r\n]*"
        warnings.warn("ERROR: found illegal ';', removing", RuntimeWarning)
        t.value = t.value[1:]
        t.type = "VALUE"
        return t

    #def t_data_UNQ_NUM_STRING(self,t):
        #r"\d+[^ \t\r\n]*"
        #return t

    #def t_data_INTEGER(self,t):
        #r"[+-]?\d+"
        ##t.type="VALUE"
        #return t
  
    def t_data_UNQ_STRING(self,t):
        r"[^ \t\r\n\[\]$][^ \t\r\n]*"
        #t.type="VALUE"
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
        self.skipped_lines += t.value.count('\n')
        t.lexer.skip(1)
    
    ### Public methods

    def build(self, **kwargs):
        # set re.MULTILINE while preserving any user reflags
        re_old = 0
        if "reflags" in kwargs.keys():
            re_old = kwargs["reflags"]
        kwargs["reflags"] = re_old | re.MULTILINE | re.I

        self.lexer = lex.lex(module=self,**kwargs)

    def test(self, data):
        self.lexer.input(data)
        while True:
            token = self.lexer.token()
            if not token:
                break
            print token 

if len(sys.argv) == 2:
    filename = sys.argv[1]
    with open(filename) as fh:
        m = CIFlex()
        m.build(debug=1)
        m.test(fh.read())

# vim:sw=4:ts=4:expandtab
