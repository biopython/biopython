# Copyright 2001 by Gavin E. Crooks.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.



# TODO: This code should be relocated somewhere else.
import sys, os

def findResource(resource_name):
    """Return the location of the given resource.

    It's convenient to store associated resource files (such as text files, or
    image files) within the Python package that uses them. This
    function will search the Python module path for these files, and return
    the absolute location of the file within the file system. 

    For example, suppose the package 'SomePackage' contains a module,
    'SomeModule.py', and an associated resource, 'SomeModule.txt'. Then
    'findResource("SomePackage/SomeModule.txt")' will return the absolute path
    to the resource, which will look something like this:
    '/usr/local/lib/python/site-packages/SomePackage/SomeModule.txt'
    """
    #Adapted from Finding files on the Python path, Mitch Chapman, 2001/03/15,
    #http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/52224

    for dirname in sys.path:
        candidate = os.path.join(dirname, resource_name)
        if os.path.isfile(candidate):
            return candidate
    raise IOError("Can't find file %s" % filename)
















