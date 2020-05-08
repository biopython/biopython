# -*- coding: utf-8 -*-
"""Implementation of pydocstyle integration with Flake8.

pydocstyle docstrings convention needs error code and class parser for be
included as module into flake8
"""

import re
import sys
try:
    import pydocstyle as pep257
    module_name = 'pydocstyle'
except ImportError:
    import pep257
    module_name = 'pep257'

if sys.version_info >= (3, 2):
    from tokenize import open as tokenize_open
else:
    tokenize_open = open

__version__ = '1.5.0'
__all__ = ('pep257Checker',)


class _ContainsAll(object):
    def __contains__(self, code):  # type: (str) -> bool
        return True


class EnvironError(pep257.Error):

    def __init__(self, err):
        super(EnvironError, self).__init__(
            code='D998',
            short_desc='EnvironmentError: ' + str(err),
            context=None,
        )

    @property
    def line(self):
        """Return 0 as line number for EnvironmentError."""
        return 0


class AllError(pep257.Error):

    def __init__(self, err):
        super(AllError, self).__init__(
            code='D999',
            short_desc=str(err).partition('\n')[0],
            context=None,
        )

    @property
    def line(self):
        """pep257.AllError does not contain line number. Return 0 instead."""
        return 0


class pep257Checker(object):
    """Flake8 needs a class to check python file."""

    name = 'flake8-docstrings'
    version = '{}, {}: {}'.format(__version__, module_name, pep257.__version__)

    def __init__(self, tree, filename, lines):
        """Placeholder."""
        self.tree = tree
        self.filename = filename
        self.checker = pep257.ConventionChecker()
        self.source = ''.join(lines)

    @classmethod
    def add_options(cls, parser):
        """Add plugin configuration option to flake8."""
        parser.add_option(
            '--docstring-convention', action='store', parse_from_config=True,
            default="pep257", choices=sorted(pep257.conventions) + ['all'],
            help=(
                "pydocstyle docstring convention, default 'pep257'. "
                "Use the special value 'all' to enable all codes (note: "
                "some codes are conflicting so you'll need to then exclude "
                "those)."
            )
        )
        parser.add_option(
            "--ignore-decorators",
            action="store",
            parse_from_config=True,
            default=None,
            help=(
                "pydocstyle ignore-decorators regular expression, "
                "default None. "
                "Ignore any functions or methods that are decorated by "
                "a function with a name fitting this regular expression. "
                "The default is not ignore any decorated functions. "
            ),
        )

    @classmethod
    def parse_options(cls, options):
        """Parse the configuration options given to flake8."""
        cls.convention = options.docstring_convention
        cls.ignore_decorators = (
            re.compile(options.ignore_decorators) if options.ignore_decorators
            else None
        )

    def _check_source(self):
        try:
            for err in self.checker.check_source(
                self.source,
                self.filename,
                ignore_decorators=self.ignore_decorators,
            ):
                yield err
        except pep257.AllError as err:
            yield AllError(err)
        except EnvironmentError as err:
            yield EnvironError(err)

    def run(self):
        """Use directly check() api from pydocstyle."""
        if self.convention == 'all':
            checked_codes = _ContainsAll()
        else:
            checked_codes = (
                pep257.conventions[self.convention] | {'D998', 'D999'}
            )
        for error in self._check_source():
            if isinstance(error, pep257.Error) and error.code in checked_codes:
                # NOTE(sigmavirus24): Fixes GitLab#3
                message = '%s %s' % (error.code, error.short_desc)
                yield (error.line, 0, message, type(self))
