#!/usr/bin/env python
#
# Copyright 2002 by Michael Hoffman.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""A set of generic bits of code under Bio.GFF (DEPRECATED).

This is part of the "old" Bio.GFF module by Michael Hoffman, which offered
access to a MySQL database holding GFF data loaded by BioPerl. This code has
now been deprecated, and will probably be removed in order to free the Bio.GFF
namespace for a new GFF parser in Biopython (including GFF3 support).
"""

import os
import sys
import tempfile

class AppendableListDictionary(dict):
    """
    a dictionary of lists
    """
    def append_to(self, key, value):
        try:
            dict.__getitem__(self, key).append(value)
        except KeyError:
            self[key] = [value]

class ForgivingDictionary(AppendableListDictionary):
    def __getitem__(self, key):
        try:
            return dict.__getitem__(self, key)
        except KeyError:
            return None

class VerboseDict(dict):
    def __str__(self):
        dict_copy = {}
        for key in self:
            dict_copy[key] = str(self[key])
        return str(dict_copy)

class VerboseList(list):
    def __str__(self):
        return str(map(lambda x: str(x), self))

class TempFile(file):
    def __init__(self, suffix = ".python-temp", keep = 0):
        self.removed = 0
        self.keep = keep
        # XXX: this is a race condition:
        file.__init__(self, tempfile.mktemp(suffix), "w")

    def __del__(self):
        self.remove()

    def remove(self):
        if self.keep == 0:
            if self.removed == 0:
                try:
                    try:
                        self.close()
                        os.remove(self.name)
                    finally:
                        self.removed = 1
                except OSError:
                    pass

class SurrogateNotInitedError(AttributeError):
    pass

class Surrogate(object):
    """
    the data is stored in _data
    """
    def __init__(self, data):
        self._data = data

    def __getattr__(self, name):
        if name == "_data":
            raise SurrogateNotInitedError(name)
        else:
            try:
                return getattr(self._data, name)
            except SurrogateNotInitedError:
                raise SurrogateNotInitedError(name)


def defline_text(defline):
    if defline[0] == ">":
        return defline[1:]
    else:
        return defline

def is_nestable(x):
    """
    Returns 1 if x is a tuple or list (sequence types that can nest)
    Returns 0 otherwise

    >>> is_nestable("string")
    0
    >>> is_nestable((0,))
    1
    >>> is_nestable(range(5))
    1
    """
    return isinstance(x, (tuple, list))

def dump_list(l):
    """
    returns strings of list
    """
    try:
        return '[%s]' % ', '.join(map(str, l))
    except TypeError:
        return str(l)

def reverse_text(text):
    """
    >>> reverse_text('abracadabra')
    'arbadacarba'
    """
    l = list(text)
    l.reverse()
    return ''.join(l)

class ArgsParser(object):
    """
    >>> unparsed_args = ["moocow"]
    >>> args = ArgsParser(unparsed_args, [('infile', 'defaultin'), ('outfile', 'defaultout')])
    >>> args.infile
    'moocow'
    >>> args.outfile
    'defaultout'
    """
    def __init__(self, args, defaults):
        for i, default in enumerate(defaults):
            try:
                self.__dict__[default[0]] = args[i]
                continue
            except TypeError:
                pass
            except IndexError:
                pass
            self.__dict__[default[0]] = default[1]

def all(iterator):
    return [item for item in iterator]

def _test():
    """Run the Bio.GFF.GenericTools module's doctests (PRIVATE).

    This will try and locate the unit tests directory, and run the doctests
    from there in order that the relative paths used in the examples work.
    """
    import doctest
    print "Runing doctests..."
    doctest.testmod()
    print "Done"

if __name__ == "__main__":
    if __debug__:
        _test()
