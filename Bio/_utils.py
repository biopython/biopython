# Copyright 2010 by Eric Talevich. All rights reserved.
# Copyright 2012 by Wibowo Arindrarto. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Common utility functions for various Bio submodules."""

from __future__ import print_function

import os


def iterlen(items):
    """Count the number of items in an iterable.

    If the argument supports len(items), and some iterators do, then
    this returns len(items). Otherwise it will scan over the entries
    in order to count them.

    Exhausts a generator, but doesn't require creating a full list.

    >>> iterlen("abcde")
    5
    >>> iterlen(iter("abcde"))
    5

    """
    try:
        # e.g. Under Python 2, the xrange iterator defines __len__
        return len(items)
    except TypeError:
        for i, x in enumerate(items):
            count = i
        return count + 1


def read_forward(handle):
    """Reads through whitespaces, returns the first non-whitespace line."""
    while True:
        line = handle.readline()
        # if line is empty or line has characters and stripping does not remove
        # them, return the line
        if (not line) or (line and line.strip()):
            return line


def trim_str(string, max_len, concat_char):
    """Truncates the given string for display."""
    if len(string) > max_len:
        return string[:max_len - len(concat_char)] + concat_char
    return string


def getattr_str(obj, attr, fmt=None, fallback='?'):
    """Returns a string of the given object's attribute, defaulting to the
    fallback value if attribute is not present."""
    if hasattr(obj, attr):
        if fmt is not None:
            return fmt % getattr(obj, attr)
        return str(getattr(obj, attr))
    return fallback


def find_test_dir(start_dir=None):
    """Finds the absolute path of Biopython's Tests directory.

    Arguments:
    start_dir -- Initial directory to begin lookup (default to current dir)

    If the directory is not found up the filesystem's root directory, an
    exception will be raised.

    """
    if not start_dir:
        # no callbacks in function signatures!
        # defaults to the current directory
        # (using __file__ would give the installed Biopython)
        start_dir = "."

    target = os.path.abspath(start_dir)
    while True:
        if os.path.isdir(os.path.join(target, "Bio")) \
        and os.path.isdir(os.path.join(target, "Tests")):
            # Good, we're in the Biopython root now
            return os.path.abspath(os.path.join(target, "Tests"))
        # Recurse up the tree
        # TODO - Test this on Windows
        new, tmp = os.path.split(target)
        if target == new:
            # Reached root
            break
        target = new
    raise ValueError("Not within Biopython source tree: %r" % os.path.abspath(start_dir))


def run_doctest(target_dir=None, *args, **kwargs):
    """Runs doctest for the importing module."""
    import doctest

    # default doctest options
    default_kwargs = {
        'optionflags': doctest.ELLIPSIS,
    }
    kwargs.update(default_kwargs)

    cur_dir = os.path.abspath(os.curdir)

    print("Runing doctests...")
    try:
        os.chdir(find_test_dir(target_dir))
        doctest.testmod(*args, **kwargs)
    finally:
        # and revert back to initial directory
        os.chdir(cur_dir)
    print("Done")

if __name__ == "__main__":
    run_doctest()
