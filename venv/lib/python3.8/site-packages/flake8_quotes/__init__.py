import optparse
import tokenize
import warnings

# Polyfill stdin loading/reading lines
# https://gitlab.com/pycqa/flake8-polyfill/blob/1.0.1/src/flake8_polyfill/stdin.py#L52-57
try:
    from flake8.engine import pep8
    stdin_get_value = pep8.stdin_get_value
    readlines = pep8.readlines
except ImportError:
    from flake8 import utils
    import pycodestyle
    stdin_get_value = utils.stdin_get_value
    readlines = pycodestyle.readlines

from flake8_quotes.__about__ import __version__
from flake8_quotes.docstring_detection import get_docstring_tokens


class QuoteChecker(object):
    name = __name__
    version = __version__

    INLINE_QUOTES = {
        # When user wants only single quotes
        '\'': {
            'good_single': '\'',
            'bad_single': '"',
        },
        # When user wants only double quotes
        '"': {
            'good_single': '"',
            'bad_single': '\'',
        },
    }
    # Provide aliases for Windows CLI support
    #   https://github.com/zheller/flake8-quotes/issues/49
    INLINE_QUOTES['single'] = INLINE_QUOTES['\'']
    INLINE_QUOTES['double'] = INLINE_QUOTES['"']

    MULTILINE_QUOTES = {
        '\'': {
            'good_multiline': '\'\'\'',
            'bad_multiline': '"""',
        },
        '"': {
            'good_multiline': '"""',
            'bad_multiline': '\'\'\'',
        },
    }
    # Provide Windows CLI and multi-quote aliases
    MULTILINE_QUOTES['single'] = MULTILINE_QUOTES['\'']
    MULTILINE_QUOTES['double'] = MULTILINE_QUOTES['"']
    MULTILINE_QUOTES['\'\'\''] = MULTILINE_QUOTES['\'']
    MULTILINE_QUOTES['"""'] = MULTILINE_QUOTES['"']

    DOCSTRING_QUOTES = {
        '\'': {
            'good_docstring': '\'\'\'',
            'bad_docstring': '"""',
        },
        '"': {
            'good_docstring': '"""',
            'bad_docstring': '\'\'\'',
        },
    }
    # Provide Windows CLI and docstring-quote aliases
    DOCSTRING_QUOTES['single'] = DOCSTRING_QUOTES['\'']
    DOCSTRING_QUOTES['double'] = DOCSTRING_QUOTES['"']
    DOCSTRING_QUOTES['\'\'\''] = DOCSTRING_QUOTES['\'']
    DOCSTRING_QUOTES['"""'] = DOCSTRING_QUOTES['"']

    def __init__(self, tree, lines=None, filename='(none)'):
        self.filename = filename
        self.lines = lines

    @staticmethod
    def _register_opt(parser, *args, **kwargs):
        """
        Handler to register an option for both Flake8 3.x and 2.x.

        This is based on:
        https://github.com/PyCQA/flake8/blob/3.0.0b2/docs/source/plugin-development/cross-compatibility.rst#option-handling-on-flake8-2-and-3

        It only supports `parse_from_config` from the original function and it
        uses the `Option` object returned to get the string.
        """
        try:
            # Flake8 3.x registration
            parser.add_option(*args, **kwargs)
        except (optparse.OptionError, TypeError):
            # Flake8 2.x registration
            parse_from_config = kwargs.pop('parse_from_config', False)
            option = parser.add_option(*args, **kwargs)
            if parse_from_config:
                parser.config_options.append(option.get_opt_string().lstrip('-'))

    @classmethod
    def add_options(cls, parser):
        cls._register_opt(parser, '--quotes', action='store',
                          parse_from_config=True, type='choice',
                          choices=sorted(cls.INLINE_QUOTES.keys()),
                          help='Deprecated alias for `--inline-quotes`')
        cls._register_opt(parser, '--inline-quotes', default='\'',
                          action='store', parse_from_config=True, type='choice',
                          choices=sorted(cls.INLINE_QUOTES.keys()),
                          help='Quote to expect in all files (default: \')')
        cls._register_opt(parser, '--multiline-quotes', default=None, action='store',
                          parse_from_config=True, type='choice',
                          choices=sorted(cls.MULTILINE_QUOTES.keys()),
                          help='Quote to expect in all files (default: """)')
        cls._register_opt(parser, '--docstring-quotes', default=None, action='store',
                          parse_from_config=True, type='choice',
                          choices=sorted(cls.DOCSTRING_QUOTES.keys()),
                          help='Quote to expect in all files (default: """)')
        cls._register_opt(parser, '--avoid-escape', default=None, action='store_true',
                          parse_from_config=True,
                          help='Avoiding escaping same quotes in inline strings (enabled by default)')
        cls._register_opt(parser, '--no-avoid-escape', dest='avoid_escape', default=None, action='store_false',
                          parse_from_config=False,
                          help='Disable avoiding escaping same quotes in inline strings')

    @classmethod
    def parse_options(cls, options):
        # Define our default config
        # cls.config = {good_single: ', good_multiline: ''', bad_single: ", bad_multiline: """}
        cls.config = {}
        cls.config.update(cls.INLINE_QUOTES['\''])
        cls.config.update(cls.MULTILINE_QUOTES['"""'])
        cls.config.update(cls.DOCSTRING_QUOTES['"""'])

        # If `options.quotes` was specified, then use it
        if hasattr(options, 'quotes') and options.quotes is not None:
            # https://docs.python.org/2/library/warnings.html#warnings.warn
            warnings.warn('flake8-quotes has deprecated `quotes` in favor of `inline-quotes`. '
                          'Please update your configuration')
            cls.config.update(cls.INLINE_QUOTES[options.quotes])
        # Otherwise, use the supported `inline_quotes`
        else:
            # cls.config = {good_single: ', good_multiline: """, bad_single: ", bad_multiline: '''}
            #   -> {good_single: ", good_multiline: """, bad_single: ', bad_multiline: '''}
            cls.config.update(cls.INLINE_QUOTES[options.inline_quotes])

        # If multiline quotes was specified, overload our config with those options
        if hasattr(options, 'multiline_quotes') and options.multiline_quotes is not None:
            # cls.config = {good_single: ', good_multiline: """, bad_single: ", bad_multiline: '''}
            #   -> {good_single: ', good_multiline: ''', bad_single: ", bad_multiline: """}
            cls.config.update(cls.MULTILINE_QUOTES[options.multiline_quotes])

        # If docstring quotes was specified, overload our config with those options
        if hasattr(options, 'docstring_quotes') and options.docstring_quotes is not None:
            cls.config.update(cls.DOCSTRING_QUOTES[options.docstring_quotes])

        # If avoid escaped specified, add to config
        if hasattr(options, 'avoid_escape') and options.avoid_escape is not None:
            cls.config.update({'avoid_escape': options.avoid_escape})
        else:
            cls.config.update({'avoid_escape': True})

    def get_file_contents(self):
        if self.filename in ('stdin', '-', None):
            return stdin_get_value().splitlines(True)
        else:
            if self.lines:
                return self.lines
            else:
                return readlines(self.filename)

    def run(self):
        file_contents = self.get_file_contents()

        noqa_line_numbers = self.get_noqa_lines(file_contents)
        errors = self.get_quotes_errors(file_contents)

        for error in errors:
            if error.get('line') not in noqa_line_numbers:
                yield (error.get('line'), error.get('col'), error.get('message'), type(self))

    def get_noqa_lines(self, file_contents):
        tokens = [Token(t) for t in tokenize.generate_tokens(lambda L=iter(file_contents): next(L))]
        return [token.start_row
                for token in tokens
                if token.type == tokenize.COMMENT and token.string.endswith('noqa')]

    def get_quotes_errors(self, file_contents):
        tokens = [Token(t) for t in tokenize.generate_tokens(lambda L=iter(file_contents): next(L))]
        docstring_tokens = get_docstring_tokens(tokens)

        for token in tokens:

            if token.type != tokenize.STRING:
                # ignore non strings
                continue

            # Remove any prefixes in strings like `u` from `u"foo"`
            # DEV: `last_quote_char` is 1 character, even for multiline strings
            #   `"foo"`   -> `"foo"`
            #   `b"foo"`  -> `"foo"`
            #   `br"foo"` -> `"foo"`
            #   `b"""foo"""` -> `"""foo"""`
            last_quote_char = token.string[-1]
            first_quote_index = token.string.index(last_quote_char)
            prefix = token.string[:first_quote_index].lower()
            unprefixed_string = token.string[first_quote_index:]

            # Determine if our string is multiline-based
            #   "foo"[0] * 3 = " * 3 = """
            #   "foo"[0:3] = "fo
            #   """foo"""[0:3] = """
            is_docstring = token in docstring_tokens
            is_multiline_string = unprefixed_string[0] * 3 == unprefixed_string[0:3]
            start_row, start_col = token.start

            # If our string is a docstring
            # DEV: Docstring quotes must come before multiline quotes as it can as a multiline quote
            if is_docstring:
                if self.config['good_docstring'] in unprefixed_string:
                    continue

                yield {
                    'message': 'Q002 Remove bad quotes from docstring',
                    'line': start_row,
                    'col': start_col,
                }
            # Otherwise if our string is multiline
            elif is_multiline_string:
                # If our string is or containing a known good string, then ignore it
                #   (""")foo""" -> good (continue)
                #   '''foo(""")''' -> good (continue)
                #   (''')foo''' -> possibly bad
                if self.config['good_multiline'] in unprefixed_string:
                    continue

                # Output our error
                yield {
                    'message': 'Q001 Remove bad quotes from multiline string',
                    'line': start_row,
                    'col': start_col,
                }
            # Otherwise (string is inline quote)
            else:
                #   'This is a string'       -> Good
                #   'This is a "string"'     -> Good
                #   'This is a \"string\"'   -> Good
                #   'This is a \'string\''   -> Bad (Q003)  Escaped inner quotes
                #   '"This" is a \'string\'' -> Good        Changing outer quotes would not avoid escaping
                #   "This is a string"       -> Bad (Q000)
                #   "This is a 'string'"     -> Good        Avoids escaped inner quotes
                #   "This is a \"string\""   -> Bad (Q000)
                #   "\"This\" is a 'string'" -> Good
                
                string_contents = unprefixed_string[1:-1]
                
                # If string preferred type, check for escapes
                if last_quote_char == self.config['good_single']:
                    if not self.config['avoid_escape'] or 'r' in prefix:
                        continue
                    if self.config['good_single'] in string_contents and not self.config['bad_single'] in string_contents:
                        yield {
                            'message': 'Q003 Change outer quotes to avoid escaping inner quotes',
                            'line': start_row,
                            'col': start_col,
                        }
                    continue
                
                # If not preferred type, only allow use to avoid escapes.
                if not self.config['good_single'] in string_contents:
                    yield {
                        'message': 'Q000 Remove bad quotes',
                        'line': start_row,
                        'col': start_col,
                    }
        


class Token:
    """Python 2 and 3 compatible token"""
    def __init__(self, token):
        self.token = token

    @property
    def type(self):
        return self.token[0]

    @property
    def string(self):
        return self.token[1]

    @property
    def start(self):
        return self.token[2]

    @property
    def start_row(self):
        return self.token[2][0]

    @property
    def start_col(self):
        return self.token[2][1]
