# Copyright 2002 by Jeffrey Chang, Andrew Dalke.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This is based on some older code by Andrew Dalke.

"""This module implements some base classes used in the Registry
system for Biopython.


Classes:
Registry             Implements a Biopython Registry.
RegisterableObject   Base class for objects in the Registry.
RegisterableGroup    Base class for groups of objects in the Registry.

"""
import re

import _support

_legal_abbrev = re.compile(r"[a-zA-Z][a-zA-Z0-9_]*$")
check_abbrev = _legal_abbrev.match

class Registry:
    """This is a dictionary-like object for storing and retrieving
    objects in a registry.

    Methods:
    register    Add a RegisterableObject into the Registry.
        Dictionary interface:
    __getitem__
    get
    keys
    values
    items

    """
    def __init__(self, name, load_path=None):
        """Registry(name[, load_path])

        Create a new registry.  name is the name of the registry.
        load_path is an optional path (e.g. Bio.config.dbdefs) that
        contains objects for the registry.

        """
        self._name = name
        self._load_path = load_path
        self._name_table, self._abbrev_table = {}, {}
        self._autoloaded = self._autoloading = 0

    def _autoload(self):
        if self._autoloaded or self._autoloading:
            return
        self._autoloading = 1
        self._load(self._load_path)
        self._autoloading = 0
        self._autoloaded = 1

    def _load(self, path):
        if path is None:
            return
        # Get a list of all the modules in that path.
        modulenames = _support.find_submodules(path)
        modulenames = filter(lambda x: not x.startswith("_"), modulenames)
        modulenames.sort()

        # Now load each one of the modules and look for
        # RegisterableObject objects in them.
        
        for name in modulenames:
            module = _support.load_module(name)
            for name, obj in module.__dict__.items():
                if name.startswith("_") or \
                       not isinstance(obj, RegisterableObject):
                    continue
                self.register(obj)

    def register(self, obj):
        """S.register(self, obj)

        Add an object to the registry.  obj must be a
        RegisterableObject object.

        """
        self._autoload()
        name, abbrev = obj.name, obj.abbrev
        abbrev = abbrev or name
        if self._name_table.has_key(name):
            raise ValueError("%r is a duplicate entry" % (name,))
        if self._abbrev_table.has_key(abbrev):
            raise ValueError("%r is a duplicate entry" % (abbrev,))

        self._name_table[name] = obj
        self._abbrev_table[abbrev] = obj

    def __getitem__(self, name):
        self._autoload()
        return self._name_table[name]  # raises KeyError for unknown entries

    def get(self, name, default=None):
        self._autoload()
        return self._name_table.get(name, default)

    def keys(self):
        self._autoload()
        return self._name_table.keys()
    def values(self):
        self._autoload()
        return self._name_table.values()
    def items(self):
        self._autoload()
        return self._name_table.items()

    def __str__(self):
        objs = self.keys()
        objs.sort()
        if not objs:
            return self._name
        obj_str = ', '.join(map(repr, objs))
        return "%s, exporting %s" % (self._name, obj_str)
    __repr__ = __str__

class RegisterableObject:
    """This is a base class for objects that can be added to a registry.

    Members:
    name                The name of the object.
    abbrev              An abbreviation for the name
    doc                 Documentation describing the object.

    """
    def __init__(self, name, abbrev, doc):
        """RegisterableObject(name, abbrev, doc)"""
        self.name = name
        self.abbrev = abbrev or name
        if not check_abbrev(self.abbrev):
            raise ValueError, "abbrev name of %r is not allowed" % self.abbrev
        self.doc = doc

        
class RegisterableGroup(RegisterableObject):
    """This is a base class for a RegisterableObject that groups many
    objects together.

    Methods:
    add         Add an object to the end of the group.
    add_after   Add an object to the group after another object.
    add_before  Add an object to the group before another object.

    """
    def __init__(self, name, abbrev, doc):
        RegisterableObject.__init__(self, name, abbrev, doc)
        self.objs = []
    def add(self, obj, index=None):
        if index is None:
            index = len(self.objs)
        self.objs.insert(index, obj)
    def add_after(self, obj, after):
        for i in range(len(self.objs)):
            if self.objs[i] == after:
                break
        else:
            raise ValueError, "I couldn't find the insertion point"
        self.add(obj, i+1)
    def add_before(self, obj, before):
        for i in range(len(self.objs)):
            if self.objs[i] == before:
                break
        else:
            raise ValueError, "I couldn't find the insertion point"
        self.add(obj, i)
