# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Base class for Residue, Chain, Model and Structure classes.

It is a simple container class, with list and dictionary like properties.
"""

from collections import deque
from copy import copy

import numpy as np

from Bio.PDB.PDBExceptions import PDBConstructionException


class Entity:
    """Basic container object for PDB hierarchy.

    Structure, Model, Chain and Residue are subclasses of Entity.
    It deals with storage and lookup.
    """

    def __init__(self, id):
        """Initialize the class."""
        self._id = id
        self.full_id = None
        self.parent = None
        self.child_list = []
        self.child_dict = {}
        # Dictionary that keeps additional properties
        self.xtra = {}

    # Special methods

    def __len__(self):
        """Return the number of children."""
        return len(self.child_list)

    def __getitem__(self, id):
        """Return the child with given id."""
        return self.child_dict[id]

    def __delitem__(self, id):
        """Remove a child."""
        return self.detach_child(id)

    def __contains__(self, id):
        """Check if there is a child element with the given id."""
        return id in self.child_dict

    def __iter__(self):
        """Iterate over children."""
        yield from self.child_list

    # Generic id-based comparison methods considers all parents as well as children
    # Works for all Entities - Atoms have comparable custom operators
    def __eq__(self, other):
        """Test for equality. This compares full_id including the IDs of all parents."""
        if isinstance(other, type(self)):
            if self.parent is None:
                return self.id == other.id
            else:
                return self.full_id[1:] == other.full_id[1:]
        else:
            return NotImplemented

    def __ne__(self, other):
        """Test for inequality."""
        if isinstance(other, type(self)):
            if self.parent is None:
                return self.id != other.id
            else:
                return self.full_id[1:] != other.full_id[1:]
        else:
            return NotImplemented

    def __gt__(self, other):
        """Test greater than."""
        if isinstance(other, type(self)):
            if self.parent is None:
                return self.id > other.id
            else:
                return self.full_id[1:] > other.full_id[1:]
        else:
            return NotImplemented

    def __ge__(self, other):
        """Test greater or equal."""
        if isinstance(other, type(self)):
            if self.parent is None:
                return self.id >= other.id
            else:
                return self.full_id[1:] >= other.full_id[1:]
        else:
            return NotImplemented

    def __lt__(self, other):
        """Test less than."""
        if isinstance(other, type(self)):
            if self.parent is None:
                return self.id < other.id
            else:
                return self.full_id[1:] < other.full_id[1:]
        else:
            return NotImplemented

    def __le__(self, other):
        """Test less or equal."""
        if isinstance(other, type(self)):
            if self.parent is None:
                return self.id <= other.id
            else:
                return self.full_id[1:] <= other.full_id[1:]
        else:
            return NotImplemented

    def __hash__(self):
        """Hash method to allow uniqueness (set)."""
        return hash(self.full_id)

    # Private methods

    def _reset_full_id(self):
        """Reset the full_id (PRIVATE).

        Resets the full_id of this entity and
        recursively of all its children based on their ID.
        """
        for child in self:
            try:
                child._reset_full_id()
            except AttributeError:
                pass  # Atoms do not cache their full ids.
        self.full_id = self._generate_full_id()

    def _generate_full_id(self):
        """Generate full_id (PRIVATE).

        Generate the full_id of the Entity based on its
        Id and the IDs of the parents.
        """
        entity_id = self.get_id()
        parts = [entity_id]
        parent = self.get_parent()
        while parent is not None:
            entity_id = parent.get_id()
            parts.append(entity_id)
            parent = parent.get_parent()
        parts.reverse()
        return tuple(parts)

    # Public methods

    @property
    def id(self):
        """Return identifier."""
        return self._id

    @id.setter
    def id(self, value):
        """Change the id of this entity.

        This will update the child_dict of this entity's parent
        and invalidate all cached full ids involving this entity.

        @raises: ValueError
        """
        if value == self._id:
            return
        if self.parent:
            if value in self.parent.child_dict:
                raise ValueError(
                    f"Cannot change id from `{self._id}` to `{value}`."
                    f" The id `{value}` is already used for a sibling of this entity."
                )
            del self.parent.child_dict[self._id]
            self.parent.child_dict[value] = self

        self._id = value
        self._reset_full_id()

    def get_level(self):
        """Return level in hierarchy.

        A - atom
        R - residue
        C - chain
        M - model
        S - structure
        """
        return self.level

    def set_parent(self, entity):
        """Set the parent Entity object."""
        self.parent = entity
        self._reset_full_id()

    def detach_parent(self):
        """Detach the parent."""
        self.parent = None

    def detach_child(self, id):
        """Remove a child."""
        child = self.child_dict[id]
        child.detach_parent()
        del self.child_dict[id]
        self.child_list.remove(child)

    def add(self, entity):
        """Add a child to the Entity."""
        entity_id = entity.get_id()
        if self.has_id(entity_id):
            raise PDBConstructionException(f"{entity_id} defined twice")
        entity.set_parent(self)
        self.child_list.append(entity)
        self.child_dict[entity_id] = entity

    def insert(self, pos, entity):
        """Add a child to the Entity at a specified position."""
        entity_id = entity.get_id()
        if self.has_id(entity_id):
            raise PDBConstructionException(f"{entity_id} defined twice")
        entity.set_parent(self)
        self.child_list[pos:pos] = [entity]
        self.child_dict[entity_id] = entity

    def get_iterator(self):
        """Return iterator over children."""
        yield from self.child_list

    def get_list(self):
        """Return a copy of the list of children."""
        return copy(self.child_list)

    def has_id(self, id):
        """Check if a child with given id exists."""
        return id in self.child_dict

    def get_parent(self):
        """Return the parent Entity object."""
        return self.parent

    def get_id(self):
        """Return the id."""
        return self.id

    def get_full_id(self):
        """Return the full id.

        The full id is a tuple containing all id's starting from
        the top object (Structure) down to the current object. A full id for
        a Residue object e.g. is something like:

        ("1abc", 0, "A", (" ", 10, "A"))

        This corresponds to:

        Structure with id "1abc"
        Model with id 0
        Chain with id "A"
        Residue with id (" ", 10, "A")

        The Residue id indicates that the residue is not a hetero-residue
        (or a water) because it has a blank hetero field, that its sequence
        identifier is 10 and its insertion code "A".
        """
        if self.full_id is None:
            self.full_id = self._generate_full_id()
        return self.full_id

    def transform(self, rot, tran):
        """Apply rotation and translation to the atomic coordinates.

        :param rot: A right multiplying rotation matrix
        :type rot: 3x3 NumPy array

        :param tran: the translation vector
        :type tran: size 3 NumPy array

        Examples
        --------
        This is an incomplete but illustrative example::

            from numpy import pi, array
            from Bio.PDB.vectors import Vector, rotmat
            rotation = rotmat(pi, Vector(1, 0, 0))
            translation = array((0, 0, 1), 'f')
            entity.transform(rotation, translation)

        """
        for o in self.get_list():
            o.transform(rot, tran)

    def center_of_mass(self, geometric=False):
        """Return the center of mass of the Entity as a numpy array.

        If geometric is True, returns the center of geometry instead.
        """
        # Recursively iterate through children until we get all atom coordinates

        if not len(self):
            raise ValueError(f"{self} does not have children")

        maybe_disordered = {"R", "C"}  # to know when to use get_unpacked_list
        only_atom_level = {"A"}

        entities = deque([self])  # start with [self] to avoid auto-unpacking
        while True:
            e = entities.popleft()
            if e.level in maybe_disordered:
                entities += e.get_unpacked_list()
            else:
                entities += e.child_list

            elevels = {e.level for e in entities}
            if elevels == only_atom_level:
                break  # nothing else to unpack

        coords = np.asarray([a.coord for a in entities], dtype=np.float32)
        if geometric:
            masses = None
        else:
            masses = np.asarray([a.mass for a in entities], dtype=np.float32)

        return np.average(coords, axis=0, weights=masses)

    def copy(self):
        """Copy entity recursively."""
        shallow = copy(self)

        shallow.child_list = []
        shallow.child_dict = {}
        shallow.xtra = copy(self.xtra)

        shallow.detach_parent()

        for child in self.child_list:
            shallow.add(child.copy())
        return shallow


class DisorderedEntityWrapper:
    """Wrapper class to group equivalent Entities.

    This class is a simple wrapper class that groups a number of equivalent
    Entities and forwards all method calls to one of them (the currently selected
    object). DisorderedResidue and DisorderedAtom are subclasses of this class.

    E.g.: A DisorderedAtom object contains a number of Atom objects,
    where each Atom object represents a specific position of a disordered
    atom in the structure.
    """

    def __init__(self, id):
        """Initialize the class."""
        self.id = id
        self.child_dict = {}
        self.selected_child = None
        self.parent = None

    # Special methods

    def __getattr__(self, method):
        """Forward the method call to the selected child."""
        if method == "__setstate__":
            # Avoid issues with recursion when attempting deepcopy
            raise AttributeError
        if not hasattr(self, "selected_child"):
            # Avoid problems with pickling
            # Unpickling goes into infinite loop!
            raise AttributeError
        return getattr(self.selected_child, method)

    def __getitem__(self, id):
        """Return the child with the given id."""
        return self.selected_child[id]

    # XXX Why doesn't this forward to selected_child?
    # (NB: setitem was here before getitem, iter, len, sub)
    def __setitem__(self, id, child):
        """Add a child, associated with a certain id."""
        self.child_dict[id] = child

    def __contains__(self, id):
        """Check if the child has the given id."""
        return id in self.selected_child

    def __iter__(self):
        """Return the number of children."""
        return iter(self.selected_child)

    def __len__(self):
        """Return the number of children."""
        return len(self.selected_child)

    def __sub__(self, other):
        """Subtraction with another object."""
        return self.selected_child - other

    # Sorting
    # Directly compare the selected child
    def __gt__(self, other):
        """Return if child is greater than other."""
        return self.selected_child > other

    def __ge__(self, other):
        """Return if child is greater or equal than other."""
        return self.selected_child >= other

    def __lt__(self, other):
        """Return if child is less than other."""
        return self.selected_child < other

    def __le__(self, other):
        """Return if child is less or equal than other."""
        return self.selected_child <= other

    # Public methods
    def copy(self):
        """Copy disorderd entity recursively."""
        shallow = copy(self)
        shallow.child_dict = {}
        shallow.detach_parent()

        for child in self.disordered_get_list():
            shallow.disordered_add(child.copy())

        return shallow

    def get_id(self):
        """Return the id."""
        return self.id

    def disordered_has_id(self, id):
        """Check if there is an object present associated with this id."""
        return id in self.child_dict

    def detach_parent(self):
        """Detach the parent."""
        self.parent = None
        for child in self.disordered_get_list():
            child.detach_parent()

    def get_parent(self):
        """Return parent."""
        return self.parent

    def set_parent(self, parent):
        """Set the parent for the object and its children."""
        self.parent = parent
        for child in self.disordered_get_list():
            child.set_parent(parent)

    def disordered_select(self, id):
        """Select the object with given id as the currently active object.

        Uncaught method calls are forwarded to the selected child object.
        """
        self.selected_child = self.child_dict[id]

    def disordered_add(self, child):
        """Add disordered entry.

        This is implemented by DisorderedAtom and DisorderedResidue.
        """
        raise NotImplementedError

    def disordered_remove(self, child):
        """Remove disordered entry.

        This is implemented by DisorderedAtom and DisorderedResidue.
        """
        raise NotImplementedError

    def is_disordered(self):
        """Return 2, indicating that this Entity is a collection of Entities."""
        return 2

    def disordered_get_id_list(self):
        """Return a list of id's."""
        # sort id list alphabetically
        return sorted(self.child_dict)

    def disordered_get(self, id=None):
        """Get the child object associated with id.

        If id is None, the currently selected child is returned.
        """
        if id is None:
            return self.selected_child
        return self.child_dict[id]

    def disordered_get_list(self):
        """Return list of children."""
        return list(self.child_dict.values())
