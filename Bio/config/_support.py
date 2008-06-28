# Copyright 2002 by Jeffrey Chang, Andrew Dalke.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This is based on some older code by Andrew Dalke.

"""Support code for dealing with registries.

Functions:
find_submodules   Find all the modules in a package.
load_module       Load a module and return it.  Raise ImportError if not found.
safe_load_module  Like load_module, but returns None if not found.

make_rate_limited_function   Limit the rate at which a function can run.

make_cached_expression       Caches the make_parser method of expressions.

"""
import sys
import os
import time

from Bio.WWW import RequestLimiter

def find_submodules(modulename):
    """find_submodules(modulename) -> list of module names

    Look inside a package or module and recursively find all the
    modules that exist within it.
    
    """
    # First, figure out where this module exists in the file system.
    module = safe_load_module(modulename)
    if not module:   # This is not a valid python module or package.
        return []
    filename = module.__file__

    # If this is an actual module (rather than a package), then
    # there's no more submodules and we're done.
    if not filename.endswith("__init__.py") and \
       not filename.endswith("__init__.pyc") and \
       not filename.endswith("__init__.pyo"):
        return [modulename]

    # Since it's a package, get a list of all the modules inside it
    # and recurse on those.
    dirname = os.path.dirname(filename)
    submodulenames = {}    # prevent duplicates
    for filename in os.listdir(dirname):
        filename = os.path.splitext(filename)[0]
        if filename == '__init__':
            continue
        elif not filename:
            continue
        name = "%s.%s" % (modulename, filename)
        submodulenames[name] = 1
    submodulenames = submodulenames.keys()
    submodulenames.sort()

    submodules = []
    for name in submodulenames:
        try:
            x = find_submodules(name)
        except ImportError, x:
            raise
            pass   # ignore things that aren't valid modules (e.g. CVS)
        else:
            submodules.extend(x)

    return submodules

def load_module(modulename):
    """load_module(modulename) -> module"""
    try:
        module = __import__(modulename, {}, {}, modulename.split(".")[:-1])
    except SyntaxError, exc:
        raise
    except ImportError, exc:
        raise ImportError("%s during import of %r" % (exc, modulename)), \
              None, sys.exc_info()[2]
    return module

def safe_load_module(modulename):
    """safe_load_module(modulename) -> module or None"""
    try:
        module = load_module(modulename)
    except ImportError, x:
        if str(x).find("during import of") == -1:
            raise
        module = None
    return module

class make_rate_limited_function:
    """make_rate_limited_function(function, delay) -> callable object

    Create a version of function that does not run more often than
    once every delay seconds.

    """
    def __init__(self, function, delay):
        self.fn = function
        self.limiter = RequestLimiter(delay)
    def __call__(self, *args, **keywds):
        self.limiter.wait()
        return self.fn(*args, **keywds)


# Only caches parsers for make_parser, not iterators
class make_cached_expression:
    """make_cached_expression(expression) -> cached expression object"""
    def __init__(self, expression):
        self.expression = expression
        self._parsers = {}   # debug_level -> parser
    def make_parser(self, debug_level=0):
        if self._parsers.get(debug_level) is None:
            parser = self.expression.make_parser(debug_level=debug_level)
            self._parsers[debug_level] = parser
        return self._parsers[debug_level].copy()

