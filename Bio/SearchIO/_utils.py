# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Common SearchIO utility functions."""

import os


def get_processor(format, mapping):
    """Returns the object to process the given format according to the mapping.

    Arguments:
    format -- Lower case string denoting one of the supported formats.
    mapping -- Dictionary of format and object name mapping.

    """
    # map file format to iterator name
    try:
        obj_info = mapping[format]
    except KeyError:
        # handle the errors with helpful messages
        if format is None:
            raise ValueError("Format required (lower case string)")
        elif not isinstance(format, basestring):
            raise TypeError("Need a string for the file format (lower case)")
        elif format != format.lower():
            raise ValueError("Format string %r should be lower case" %
                    format)
        else:
            raise ValueError("Unknown format %r. Supported formats are "
                    "%r" % (format, "', '".join(mapping.keys())))

    mod_name, obj_name = obj_info
    mod = __import__('Bio.SearchIO.%s' % mod_name, fromlist=[''])

    return getattr(mod, obj_name)


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
