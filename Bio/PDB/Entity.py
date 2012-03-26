# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.  

from copy import copy

from Bio.PDB.PDBExceptions import PDBConstructionException, PDBException

"""Base class for Residue, Chain, Model and Structure classes.

It is a simple container class, with list and dictionary like properties.
"""


class Entity(object):
    """
    Basic container object. Structure, Model, Chain and Residue
    are subclasses of Entity. It deals with storage and lookup.
    """
    def __init__(self, id):
        self.id=id
        self.full_id=None
        self.parent=None
        self.child_list=[]
        self.child_dict={}
        # Dictionary that keeps addictional properties
        self.xtra={}
    
    # Special methods   

    def __len__(self):
        "Return the number of children."
        return len(self.child_list)

    def __getitem__(self, id):
        "Return the child with given id."
        return self.child_dict[id]

    def __delitem__(self, id):
        "Remove a child."
        return self.detach_child(id)

    def __iter__(self):
        "Iterate over children."
        for child in self.child_list:
            yield child

    # Public methods    

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
        "Set the parent Entity object."
        self.parent=entity

    def detach_parent(self):
        "Detach the parent."
        self.parent=None

    def detach_child(self, id):
        "Remove a child."
        child=self.child_dict[id] 
        child.detach_parent()
        del self.child_dict[id]
        self.child_list.remove(child)

    def add(self, entity):
        "Add a child to the Entity."
        entity_id=entity.get_id()
        if self.has_id(entity_id):
            raise PDBConstructionException( \
                "%s defined twice" % str(entity_id))
        entity.set_parent(self)
        self.child_list.append(entity)
        self.child_dict[entity_id]=entity
    
    def insert(self, pos, entity):
        "Add a child to the Entity at a specified position."
        entity_id=entity.get_id()
        if self.has_id(entity_id):
            raise PDBConstructionException( \
                "%s defined twice" % str(entity_id))
        entity.set_parent(self)
        self.child_list[pos:pos] = [entity]
        self.child_dict[entity_id]=entity        

    def get_iterator(self):
        "Return iterator over children."
        for child in self.child_list:
            yield child

    def get_list(self):
        "Return a copy of the list of children."
        return copy(self.child_list)

    def has_id(self, id):
        """True if a child with given id exists."""
        return (id in self.child_dict)

    def get_parent(self):
        "Return the parent Entity object."
        return self.parent

    def get_id(self):
        "Return the id."
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
        (or a water) beacuse it has a blank hetero field, that its sequence 
        identifier is 10 and its insertion code "A".
        """
        if self.full_id==None:
            entity_id=self.get_id()
            l=[entity_id]   
            parent=self.get_parent()
            while not (parent is None):
                entity_id=parent.get_id()
                l.append(entity_id)
                parent=parent.get_parent()
            l.reverse()
            self.full_id=tuple(l)
        return self.full_id

    def transform(self, rot, tran):
        """
        Apply rotation and translation to the atomic coordinates.

        Example:
                >>> rotation=rotmat(pi, Vector(1,0,0))
                >>> translation=array((0,0,1), 'f')
                >>> entity.transform(rotation, translation)

        @param rot: A right multiplying rotation matrix
        @type rot: 3x3 Numeric array

        @param tran: the translation vector
        @type tran: size 3 Numeric array
        """
        for o in self.get_list():
            o.transform(rot, tran)

    def copy(self):
        shallow = copy(self)
        shallow.child_list = copy(self.child_list)
        shallow.child_dict = copy(self.child_dict)
        shallow.xtra = copy(self.xtra)
        shallow.detach_parent()
        for index, child in self.child_dict.items():
            shallow.detach_child(index)
            shallow.add(child.copy())
        return shallow
        

class DisorderedEntityWrapper(object):
    """
    This class is a simple wrapper class that groups a number of equivalent
    Entities and forwards all method calls to one of them (the currently selected 
    object). DisorderedResidue and DisorderedAtom are subclasses of this class.
    
    E.g.: A DisorderedAtom object contains a number of Atom objects,
    where each Atom object represents a specific position of a disordered
    atom in the structure.
    """
    def __init__(self, id):
        self.id=id
        self.child_dict={}
        self.selected_child=None
        self.parent=None    

    # Special methods

    def __getattr__(self, method):
        "Forward the method call to the selected child."
        if not hasattr(self, 'selected_child'):
            # Avoid problems with pickling
            # Unpickling goes into infinite loop!
            raise AttributeError
        return getattr(self.selected_child, method)

    def __getitem__(self, id):
        "Return the child with the given id."
        return self.selected_child[id]

    # XXX Why doesn't this forward to selected_child?
    # (NB: setitem was here before getitem, iter, len, sub)
    def __setitem__(self, id, child):
        "Add a child, associated with a certain id."
        self.child_dict[id]=child

    def __iter__(self):
        "Return the number of children."
        return iter(self.selected_child)

    def __len__(self):
        "Return the number of children."
        return len(self.selected_child)

    def __sub__(self, other):
        """Subtraction with another object."""
        return self.selected_child - other

    # Public methods    

    def get_id(self):
        "Return the id."
        return self.id

    def disordered_has_id(self, id):
        """True if there is an object present associated with this id."""
        return (id in self.child_dict)

    def detach_parent(self):
        "Detach the parent"
        self.parent=None
        for child in self.disordered_get_list():
            child.detach_parent()

    def get_parent(self):
        "Return parent."
        return self.parent

    def set_parent(self, parent):
        "Set the parent for the object and its children."
        self.parent=parent
        for child in self.disordered_get_list():
            child.set_parent(parent)

    def disordered_select(self, id):
        """Select the object with given id as the currently active object.

        Uncaught method calls are forwarded to the selected child object.
        """
        self.selected_child=self.child_dict[id]
    
    def disordered_add(self, child):
        "This is implemented by DisorderedAtom and DisorderedResidue."
        raise NotImplementedError

    def is_disordered(self):
        """
        Return 2, indicating that this Entity is a collection of Entities.
        """
        return 2

    def disordered_get_id_list(self):
        "Return a list of id's."
        l=self.child_dict.keys()
        # sort id list alphabetically
        l.sort()
        return l
        
    def disordered_get(self, id=None):
        """Get the child object associated with id.

        If id is None, the currently selected child is returned.
        """
        if id==None:
            return self.selected_child
        return self.child_dict[id]

    def disordered_get_list(self):
        "Return list of children."
        return self.child_dict.values()
