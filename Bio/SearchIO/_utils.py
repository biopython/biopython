# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Common SearchIO utility functions."""

import os


def find_test_dir(test_dir='Tests', start_dir=None):
    """Finds the absolute path of Biopython's Tests directory.

    Arguments:
    start_dir -- Initial directory to begin lookup.
    test_dir -- Biopython test directory name.

    If the directory is not found up the filesystem's root directory, an
    exception will be raised.

    """
    # no callbacks in function signatures!
    # defaults to the current _utils directory
    if start_dir is None:
        start_dir = os.path.abspath(__file__)

    # raise error if search goes all the way to root without results
    # to prevent infinite loop
    if start_dir == os.path.dirname(start_dir):
        raise IOError("%r directory not found" % test_dir)

    target_dir = os.path.join(start_dir, test_dir)
    # recurse if the test directory is not present in the current directory
    if not os.path.isdir(target_dir):
        parent = os.path.dirname(start_dir)
        return find_test_dir(test_dir=test_dir, start_dir=parent)

    return target_dir


def run_doctest(*args, **kwargs):
    """Runs doctest for the importing module."""
    import doctest

    # default doctest options
    default_kwargs = {
        'optionflags': doctest.ELLIPSIS,
    }
    kwargs.update(default_kwargs)

    test_dir = find_test_dir()
    cur_dir = os.path.abspath(os.curdir)

    print "Runing doctests..."
    # change to test directory
    os.chdir(test_dir)
    doctest.testmod(*args, **kwargs)
    # and revert back to initial directory
    os.chdir(cur_dir)
    print "Done"
