# Created 2012
# Lenna X. Peterson
# arklenna@gmail.com

import sys  # to get args
import warnings  # to warn
import time  # to time function
import re  # to supply RE args to lexer
import ply.lex as lex  # lexer
from ply.lex import TOKEN  # to assign complex docstrings to tokens


class CIFlex:
    '''
    Build PLY lexer for CIF/mmCIF data format and tokenize input.

    Arguments:
    data (default None): bulk input (i.e. fh.read())
    **kwargs: arguments to lexer (i.e. debug=1)
            Consult PLY documentation.

    '''
    # Note that PLY uses docstrings to store RE.
    def __init__(self, data=None, **kwargs):
        '''
        Build lexer with optional args, tokenize if input provided.

        '''
        # combine re.MULTILINE and re.I with user reflags
        re_old = 0
        if "reflags" in kwargs.keys():
            re_old = kwargs["reflags"]
        kwargs["reflags"] = re_old | re.MULTILINE | re.I
        self._kwargs = kwargs
        self.build()
        if data is not None:
            self._data = data
            self.input()

    #### Lexer tokens #####
    # <AnyPrintChar> : [^\r\n]
    # <NonBlankChar> : [^ \t\r\n]
    # <eol> : (\r\n|\n|\r)

    tok_type = (
        None,
        "NAME",  # 1
        "LOOP",  # 2
        "DATA",  # 3
        "SEMICOLONS",  # 4
        "DOUBLEQUOTED",  # 5
        "QUOTED",  # 6
        "SIMPLE",  # 7
    )

    tokens = (
        "COMMENT",
        # Reserved words:
        '3',  # DATA
        '2',  # LOOP
        "GLOBAL",
        "SAVE",
        "STOP",
        '1',  # NAME
        # Value types:
        '4',  # SEMICOLONS
        "SEMI_ERROR",
        '6',  # QUOTED
        '5',  # DOUBLEQUOTED
        '7',  # SIMPLE
    )

    states = (
        ("data", "exclusive"),
        ("loop", "exclusive"),
    )

    # Literal hash, any print chars
    def t_ANY_COMMENT(self, t):
        r"""\#[^\r\n]*"""
        t.lexer.lineno += t.value.count('\n')
        if t.lexer.current_state() == "loop":
            t.lexer.pop_state()
        return None

    # DATA=3
    # DATA_, non blank chars
    def t_3(self, t):
        r"DATA_[^ \t\r\n]+"
        t.lexer.push_state("data")
        t.value = t.value[5:]
        return t

    # LOOP=2
    def t_data_2(self, t):
        r"LOOP_"
        if t.lexer.current_state() == "loop":
            warnings.warn("ERROR: Illegal nested loop.", RuntimeWarning)
        else:
            t.lexer.push_state("loop")
            return t
        return t

    def t_ANY_GLOBAL(self, t):
        r"GLOBAL_"
        return t

    # SAVE_, non blank chars
    def t_ANY_SAVE(self, t):
        r"SAVE_[^ \t\r\n]+"
        t.value = t.value[5:]
        return t

    def t_ANY_STOP(self, t):
        r"STOP_"
        return t

    # NAME=1
    # _, non blank chars
    def t_data_loop_1(self, t):
        r"_[^ \t\r\n]+"
        return t

    # line anchor, semi, any print chars
    _semi_text_header = r"^;[^\r\n]*$"
    # line anchor, any non-semi print char, any print chars, eol,
    #     semi, whitespace
    _semi_text_lines = r"((^[^;\r\n][^\r\n]*)?(\r\n|\r|\n))*;[ \t\r\n]"
    semi_text_field = _semi_text_header + _semi_text_lines

    # SEMICOLONS=4
    # @TOKEN sets docstring of next token definition
    @TOKEN(semi_text_field)
    def t_data_loop_4(self, t):
        # add \n to count
        t.lexer.lineno += t.value.count('\n')
        # remove \n by splitting into lines and joining
        sep = ""
        if " " in t.value.strip():
            sep = " "
        t.value = sep.join(t.value.splitlines())[1:-1]
        return t

    # line anchor, semi, non blank chars
    # (that hasn't been recognized in a SEMI_TEXT_FIELD)
    def t_data_loop_SEMI_ERROR(self, t):
        r"^;[^ \t\r\n]*"
        warnings.warn("ERROR: found illegal ';', removing", RuntimeWarning)
        t.value = t.value[1:]
        t.type = "7"
        return t

    # QUOTED=6
    # cif allows internal quotes w/o trailing whitespace, i.e. 'it's a dog'
    # reluctant * prevents grabbing too much
    # sq, non-newline chars (reluctant!), sq, whitespace
    def t_data_loop_6(self, t):
        r"'[^\r\n]*?'[ \t\r\n]"
        t.value = t.value.rstrip()[1:-1]
        return t

    # DOUBLEQUOTED=5
    # dq, non-newline chars (reluctant!), dq, whitespace
    def t_data_loop_5(self, t):
        r'"[^\r\n]*?"[ \t\r\n]'
        t.value = t.value.rstrip()[1:-1]
        return t

    # SIMPLE=7
    # unquoted strings may not begin with brackets or $
    # [not space brackets $], non blank chars
    def t_data_loop_7(self, t):
        r"[^ \t\r\n\[\]$][^ \t\r\n]*"
        return t

    # Ignored characters: spaces and tabs
    t_ANY_ignore = ' \t'

    # Newline rule to track line numbers
    def t_ANY_newline(self, t):
        r'\n+'
        t.lexer.lineno += len(t.value)

    # Error handling rule
    def t_ANY_error(self, t):
        print "Illegal value '%s'" % t.value
        self.skipped_lines += t.value.count('\n')
        t.lexer.skip(1)

    ##### Public methods #####
    def build(self):
        '''Build the lexer and start timers'''
        self._lexstart = time.clock()
        self.lexer = lex.lex(module=self, **self._kwargs)
        self.skipped_lines = 0
        self.lex_init = time.clock() - self._lexstart
        print "Lexer started", self.lex_init

    def input(self):
        '''Feed the data into the lexer'''
        self.lexer.input(self._data)

    def getToken(self):
        '''Emit token type, value as a 2-tuple'''
        token = self.lexer.token()
        if token:
            return (token.type, token.value)
        return (None, None)

    def _test(self):
        '''Iterate through and print tokens and time'''
        while True:
            tok = self.lexer.token()
            if not tok:
                break
            print self.tok_type[int(tok.type)], tok.value
        self.lex_end = time.clock() - self._lexstart
        print "Lexer runtime:", self.lex_end
        print "Skipped %s lines" % self.skipped_lines


if __name__ == "__main__":
    if len(sys.argv) == 2:
        filename = sys.argv[1]
        with open(filename) as fh:
            filedump = fh.read()
        m = CIFlex(filedump, debug=1)
        m._test()

# vim:sw=4:ts=4:expandtab
