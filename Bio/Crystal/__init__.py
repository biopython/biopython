# Copyright 2002 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Represent the NDB Atlas structure (a minimal subset of PDB format) (DEPRECATED).

Hetero, Crystal and Chain exist to represent the NDB Atlas structure.  Atlas
is a minimal subset of the PDB format.  Hetero supports a 3 alphanumeric code.
The NDB web interface is located at http://ndbserver.rutgers.edu

Bio.Crystal.Hetero substitute is Bio.PDB.Atom
Bio.Crystal.Chain substitute is Bio.PDB.Chain
Bio.Crystal.Crystal substitute is Bio.PDB.Structure

Using Bio.PDB you can navigate the data as below::

    from Bio.PDB.PDBParser import PDBParser
    parser = PDBParser(PERMISSIVE=1)
    # PDB NDB Only file
    structure = parser.get_structure("001", "001_msd.pbd")
    for model in structure:
        print('Model ',model)
        for chain in model:
            print('Chain ', chain)
            for residue in chain:
                print('Res ', residue)
                for atom in residue:
                    print('Atom ', atom)

Bio.Crystal is self-contained, with the main functionality covered by Bio.PDB.
"""

import copy
import warnings
from functools import reduce

from Bio import BiopythonDeprecationWarning

warnings.warn(
    "Bio.Crystal has been deprecated, and we intend to remove it"
    " in a future release of Biopython. Please use Bio.PDB instead"
    " to parse NDB files.",
    BiopythonDeprecationWarning,
)


class CrystalError(Exception):
    """Class to manage errors."""

    pass


def wrap_line(line):
    """Add end of line at character eighty, to match PDB record standard."""
    output = ""
    for i in range(0, len(line), 80):
        output += "%s\n" % line[i : i + 80]
    return output


def validate_key(key):
    """Check if key is a string and has at least one character."""
    if not isinstance(key, str):
        raise CrystalError("chain requires a string label")
    if len(key) != 1:
        raise CrystalError("chain label should contain one letter")


class Hetero:
    """Class to support the PDB hetero codes.

    Supports only the 3 alphanumeric code.
    The annotation is available from http://xray.bmc.uu.se/hicup/xname.html
    """

    def __init__(self, data):
        """Initialize the class."""
        # Enforce string storage
        if not isinstance(data, str):
            raise CrystalError("Hetero data must be an alphanumeric string")
        if data.isalnum() == 0:
            raise CrystalError("Hetero data must be an alphanumeric string")
        if len(data) > 3:
            raise CrystalError("Hetero data may contain up to 3 characters")
        if len(data) < 1:
            raise CrystalError("Hetero data must not be empty")

        self.data = data[:].lower()

    def __eq__(self, other):
        return self.data == other.data

    def __repr__(self):
        return "%s" % self.data

    def __str__(self):
        return "%s" % self.data

    def __len__(self):
        return len(self.data)


class Chain:
    """Class representing a sequence of Hetero elements."""

    def __init__(self, residues=""):
        """Initialize the class."""
        self.data = []
        if isinstance(residues, str):
            residues = residues.replace("*", " ")
            residues = residues.strip()
            elements = residues.split()
            self.data = [Hetero(x) for x in elements]
        elif isinstance(residues, list):
            for element in residues:
                if not isinstance(element, Hetero):
                    raise CrystalError("Text must be a string")
            for residue in residues:
                self.data.append(residue)
        elif isinstance(residues, Chain):
            for residue in residues:
                self.data.append(residue)
        self.validate()

    def validate(self):
        """Check all data elements are of type Hetero."""
        data = self.data
        for element in data:
            self.validate_element(element)

    def validate_element(self, element):
        """Check the element is of type Hetero."""
        if not isinstance(element, Hetero):
            raise TypeError

    def __str__(self):
        output = ""
        for element in self.data:
            output = output + "%s " % element
        output = output.strip()
        output = wrap_line(output)
        return output

    def __eq__(self, other):
        if len(self.data) != len(other.data):
            return 0
        ok = reduce(
            lambda x, y: x and y, map(lambda x, y: x == y, self.data, other.data)
        )
        return ok

    def __ne__(self, other):
        """Return true iff self is not equal to other."""
        return not self.__eq__(other)

    def __len__(self):
        return len(self.data)

    def __getitem__(self, index):
        if isinstance(index, int):
            return self.data[index]
        elif isinstance(index, slice):
            return self.__class__(self.data[index])
        else:
            raise TypeError

    def __setitem__(self, index, value):
        if isinstance(index, int):
            try:
                self.validate_element(value)
            except TypeError:
                value = Hetero(value.lower())
            self.data[index] = value
        elif isinstance(index, slice):
            if isinstance(value, Chain):
                self.data[index] = value.data
            elif isinstance(value, type(self.data)):
                self.data[index] = value
            elif isinstance(value, str):
                self.data[index] = Chain(value).data
            else:
                raise TypeError
        else:
            raise TypeError

    def __delitem__(self, index):
        del self.data[index]

    def __contains__(self, item):
        try:
            self.validate_element(item)
        except TypeError:
            item = Hetero(item.lower())
        return item in self.data

    def append(self, item):
        """Add Hetero element."""
        try:
            self.validate_element(item)
        except TypeError:
            item = Hetero(item.lower())
        self.data.append(item)

    def insert(self, i, item):
        """Insert Hetero element in position i of the Chain."""
        try:
            self.validate_element(item)
        except TypeError:
            item = Hetero(item.lower())
        self.data.insert(i, item)

    def remove(self, item):
        """Delete Hetero element."""
        item = Hetero(item.lower())
        self.data.remove(item)

    def count(self, item):
        """Return number of elements in the Chain."""
        try:
            self.validate_element(item)
        except TypeError:
            item = Hetero(item.lower())
        return self.data.count(item)

    def index(self, item):
        """Find the index of the item."""
        try:
            self.validate_element(item)
        except TypeError:
            item = Hetero(item.lower())
        return self.data.index(item)

    def __add__(self, other):
        if isinstance(other, Chain):
            return self.__class__(self.data + other.data)
        elif isinstance(other, str):
            return self.__class__(self.data + Chain(other).data)
        else:
            raise TypeError

    def __radd__(self, other):
        if isinstance(other, Chain):
            return self.__class__(other.data + self.data)
        elif isinstance(other, str):
            return self.__class__(Chain(other).data + self.data)
        else:
            raise TypeError

    def __iadd__(self, other):
        if isinstance(other, Chain):
            self.data += other.data
        elif isinstance(other, str):
            self.data += Chain(other).data
        else:
            raise TypeError
        return self


class Crystal:
    """Represents a dictionary of labeled chains from the same structure."""

    def __init__(self, data=None):
        """Initialize the class."""
        # Enforce storage
        if not isinstance(data, dict):
            raise CrystalError("Crystal must be a dictionary")
        if data is None:
            self.data = {}
        else:
            self.data = data
        self.fix()

    def fix(self):
        """Change element of type string to type Chain.

        All elements of Crystal shall be Chain.
        """
        data = self.data
        for key in data:
            element = data[key]
            if isinstance(element, Chain):
                pass
            elif isinstance(element, str):
                data[key] = Chain(element)
            else:
                raise TypeError

    def __repr__(self):
        output = ""
        for key in sorted(self.data):
            output += "%s : %s\n" % (key, self.data[key])
        return output

    def __str__(self):
        output = ""
        for key in sorted(self.data):
            output += "%s : %s\n" % (key, self.data[key])
        return output

    def tostring(self):
        """Return Chains and correspondent Heteros."""
        return self.data

    def __len__(self):
        return len(self.data)

    def __getitem__(self, key):
        return self.data[key]

    def __setitem__(self, key, item):
        if isinstance(item, Chain):
            self.data[key] = item
        elif isinstance(item, str):
            self.data[key] = Chain(item)
        else:
            raise TypeError

    def __delitem__(self, key):
        del self.data[key]

    def clear(self):
        """Empty the data."""
        self.data.clear()

    def copy(self):
        """Copy the Crystal object."""
        return copy.copy(self)

    def keys(self):
        """Return all Chain labels."""
        return self.data.keys()

    def items(self):
        """Return all tuples (Chain label, Hetero)."""
        return self.data.items()

    def values(self):
        """Return all Hetero in the Chains."""
        return self.data.values()

    def __contains__(self, value):
        return value in self.data

    def has_key(self, key):
        """Return true if the Chain Label is in the dictionary."""
        return key in self.data

    def get(self, key, failobj=None):
        """Return Hetero for the given Chain Label."""
        return self.data.get(key, failobj)

    def setdefault(self, key, failobj=None):
        """Return Hetero for the given Chain Label, if Chain Label is not there add it."""
        if key not in self.data:
            self.data[key] = failobj
        return self.data[key]

    def popitem(self):
        """Return and delete a Chain."""
        return self.data.popitem()
