# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Common SearchIO utility functions."""

import os


TEST_DIR = 'Tests'    # Biopython test directory name
MOD_DIR = 'Bio'       # Biopython main library directory name


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


def find_test_dir(start_dir=None):
    """Finds the absolute path of Biopython's Tests directory.

    Arguments:
    start_dir -- Initial directory to begin lookup.

    If the directory is not found up the filesystem's root directory, an
    exception will be raised.

    """
    # no callbacks in function signatures!
    # defaults to the current _utils directory
    if start_dir is None:
        start_dir = os.path.dirname(os.path.abspath(__file__))

    # raise error if search goes all the way to root without results
    # to prevent infinite loop
    # this happens when the Bio directory is not a parent directory of the
    # running module
    if start_dir == os.path.dirname(start_dir):
        raise IOError("Biopython directory not found")

    parent, cur = os.path.split(start_dir)
    if cur == MOD_DIR:
        target_dir = os.path.join(parent, TEST_DIR)
        assert os.path.isdir(target_dir), \
                "Directory %r not found in Biopython's test directory" % \
                TEST_DIR
        return target_dir

    return find_test_dir(start_dir=parent)


def run_doctest(target_dir=None, *args, **kwargs):
    """Runs doctest for the importing module."""
    import doctest

    # default doctest options
    default_kwargs = {
        'optionflags': doctest.ELLIPSIS,
    }
    kwargs.update(default_kwargs)

    test_dir = find_test_dir()
    # set target directory if it's not None
    if target_dir is not None:
        target_dir = os.path.join(test_dir, target_dir)
    else:
        target_dir = test_dir
    cur_dir = os.path.abspath(os.curdir)

    print "Runing doctests..."
    # change to test directory
    os.chdir(target_dir)
    doctest.testmod(*args, **kwargs)
    # and revert back to initial directory
    os.chdir(cur_dir)
    print "Done"
