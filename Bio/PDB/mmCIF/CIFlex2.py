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

    Options:
    filename: filename to tokenize (default None)
    **kwargs: arguments to lexer (i.e. debug=1)
            Consult PLY documentation.

    '''
    # Note that PLY uses docstrings to store RE.
    def __init__(self, filename=None, **kwargs):
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
        # Pass input to lexer if provided
        if filename is not None:
            with open(filename) as fh:
                filedump = fh.read()
            self.lexer.input(filedump)

    #### Lexer tokens #####
    # REFERENCE:
    # http://www.iucr.org/resources/cif/spec/version1.1/cifsyntax

    # Tuple indices match C/flex integer tokens
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

    # PLY requires token names to be strings
    tokens = (
        "COMMENT",
        '1',  # NAME
        # Reserved words:
        "GLOBAL",
        "SAVE",
        "STOP",
        '2',  # LOOP
        '3',  # DATA
        # Value types:
        "SEMI_ERROR",
        '4',  # SEMICOLONS
        '5',  # DOUBLEQUOTED
        '6',  # QUOTED
        '7',  # SIMPLE
    )

    ### Lexer states
    # exclusive states only allow tokens assigned to that state
    # Assignments appear in between t_ and _NAME in definition
    # t_ANY_ token will be matched in any state
    states = (
        ("data", "exclusive"),
        ("loop", "exclusive"),
    )

    ### Regex description reference
    # <AnyPrintChar> : [^\r\n]
    # <NonBlankChar> : [^ \t\r\n]
    # <eol> : (\r\n|\n|\r)

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
        # This state is intended to eventually handle multiple data blocks
        t.lexer.push_state("data")
        t.value = t.value[5:]
        return t

    # LOOP=2
    def t_data_2(self, t):
        r"LOOP_"
        if t.lexer.current_state() == "loop":
            warnings.warn("LEX ERROR: Illegal nested loop.", RuntimeWarning)
        else:
            t.lexer.push_state("loop")
            return t
        return t

    #def t_ANY_GLOBAL(self, t):
        r"GLOBAL_"
        warnings.warn("LEX ERROR: Unhandled 'GLOBAL_' found", RuntimeWarning)
        return t

    # SAVE_, non blank chars
    def t_ANY_SAVE(self, t):
        r"SAVE_[^ \t\r\n]*"
        warnings.warn("LEX ERROR: Unhandled 'SAVE_' found", RuntimeWarning)
        if len(t.value) > 5:
            t.value = t.value[5:]
        return t

    def t_ANY_STOP(self, t):
        r"STOP_"
        warnings.warn("LEX ERROR: Unhandled 'STOP_' found", RuntimeWarning)
        return t

    # NAME=1
    # _, non blank chars
    def t_data_loop_1(self, t):
        r"_[^ \t\r\n]+"
        return t

    # SEMICOLONS=4
    # line anchor, semi, any print chars
    _semi_text_header = r"^;[^\r\n]*$"
    # line anchor, any non-semi print char, any print chars, eol,
    #     semi, whitespace
    _semi_text_lines = r"((^[^;\r\n][^\r\n]*)?(\r\n|\r|\n))*;[ \t\r\n]"
    semi_text_field = _semi_text_header + _semi_text_lines

    # @TOKEN sets docstring/re of next token definition
    @TOKEN(semi_text_field)
    def t_data_loop_4(self, t):
        # add \n to count
        t.lexer.lineno += t.value.count('\n')
        # remove \n by splitting into lines and joining
        # break field into list of stripped lines
        lines = [val.strip() for val in t.value.splitlines()]
        # autodetect spaces in content (i.e. sentences vs. gene sequence)
        sep = ""
        whitespace = " "
        if whitespace in t.value.strip():
            sep = whitespace
        # remove semicolons with slice
        t.value = sep.join(lines)[1:-1]
        return t

    # line anchor, semi, non blank chars
    def t_data_loop_SEMI_ERROR(self, t):
        r"^;[^ \t\r\n]*"
        # XXX Removing the semicolon may not be the ideal behavior
        warnings.warn("LEX WARNING: removing illegal ';'", RuntimeWarning)
        t.value = t.value[1:]
        t.type = "7"
        return t

    # DOUBLEQUOTED=5
    # dq, non-newline chars (reluctant!), dq, whitespace
    def t_data_loop_5(self, t):
        r'"[^\r\n]*?"[ \t\r\n]'
        t.value = t.value.rstrip()[1:-1]
        return t

    # QUOTED=6
    # cif allows internal quotes w/o trailing whitespace, i.e. 'it's a dog'
    # reluctant * prevents grabbing too much
    # sq, non-newline chars (reluctant!), sq, whitespace
    def t_data_loop_6(self, t):
        r"'[^\r\n]*?'[ \t\r\n]"
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
        self._lex_init = time.clock() - self._lexstart
        #print "Lexer started", self._lex_init

    # The following 3 classmethods (open_file, close_file, get_token)
    #     emulate the behavior of the C module for MMCIF2Dict
    @classmethod
    def open_file(cls, filename):
        '''Create class instance of this class (emulates C module)'''
        cls._MMCIF2DictInstance = CIFlex(filename)

    @classmethod
    def get_token(cls):
        '''Emit int token type, value as a 2-tuple (emulates C module)'''
        token = cls._MMCIF2DictInstance.lexer.token()
        if token:
            return (int(token.type), token.value)
        # tuple of None avoids 'NoneType not iterable' warning
        return (None, None)

    @classmethod
    def close_file(cls):
        '''Pretend to close file (emulates C module)'''
        pass

    def _test(self):
        '''Iterate through and print tokens and time'''
        while True:
            tok = self.lexer.token()
            if not tok:
                break
            # Get token name from tok_type by index
            print self.tok_type[int(tok.type)], tok.value
        self._lex_end = time.clock() - self._lexstart
        #print "Lexer runtime:", self._lex_end
        print "Lexer skipped %s lines" % self.skipped_lines


if __name__ == "__main__":
    if len(sys.argv) == 2:
        filename = sys.argv[1]
        m = CIFlex(filename, debug=1)
        m._test()

# vim:sw=4:ts=4:expandtab
