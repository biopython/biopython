# Copyright 2002 by Jeffrey Chang, Andrew Dalke.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This is based on some older code by Andrew Dalke.

"""XXX document

Functions:
find_submodules   Find all the modules in a package.
load_module       Load a module and return it.  Raise ImportError if not found.
safe_load_module  Like load_module, but returns None if not found.

make_rate_limited_function   Limit the rate at which a function can run.
make_timed_function          Limit the amount of time a function can run.

make_cached_expression       Caches the make_parser method of expressions.

"""
import sys
import os
import time

from Bio.WWW import RequestLimiter
from Bio.MultiProc.copen import copen_fn


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

class make_timed_function:
    """make_timed_function(function, timeout[, retval2pickleable_fn][, pickleable2retval_fn]) -> callable object

    Create a version of function that times out if it does not
    complete within timeout seconds.

    Currently, there's an implementation limitation such that function
    must return a pickleable object (or nothing).  If the function
    returns an object that's not pickleable, then please set
    retval2pickleable_fn and pickleable2retval_fn to a pair of
    callbacks to convert the return value of the function to a
    pickleable form.  If this is impossible, then this function should
    not be used.

    """
    def __init__(self, function, timeout,
                 retval2pickleable_fn=None, pickleable2retval_fn=None):
        self.fn = function
        self.timeout = timeout
        self.retval2pickleable_fn = retval2pickleable_fn or (lambda x: x)
        self.pickleable2retval_fn = pickleable2retval_fn or (lambda x: x)
    def _call_fn(self, *args, **keywds):
        retval = self.fn(*args, **keywds)
        return self.retval2pickleable_fn(retval)
    def __call__(self, *args, **keywds):
        end_time = time.time() + self.timeout
        handle = copen_fn(self._call_fn, *args, **keywds)
        while time.time() < end_time:
            if handle.poll():
                break
            time.sleep(0.01)
        else:
            handle.close()
            raise IOError, "timed out"
        pickleable = handle.read()
        return self.pickleable2retval_fn(pickleable)

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

