# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Chain class, used in Structure objects."""

from Bio.PDB.Entity import Entity


class Chain(Entity):
    def __init__(self, id):
        """Initialize the class."""
        self.level = "C"
        Entity.__init__(self, id)

    # Sorting methods: empty chain IDs come last.
    def __gt__(self, other):
        if isinstance(other, Chain):
            if self.id == ' ' and other.id != ' ':
                return 0
            elif self.id != ' ' and other.id == ' ':
                return 1
            else:
                return self.id > other.id
        else:
            return NotImplemented

    def __ge__(self, other):
        if isinstance(other, Chain):
            if self.id == ' ' and other.id != ' ':
                return 0
            elif self.id != ' ' and other.id == ' ':
                return 1
            else:
                return self.id >= other.id
        else:
            return NotImplemented

    def __lt__(self, other):
        if isinstance(other, Chain):
            if self.id == ' ' and other.id != ' ':
                return 0
            elif self.id != ' ' and other.id == ' ':
                return 1
            else:
                return self.id < other.id
        else:
            return NotImplemented

    def __le__(self, other):
        if isinstance(other, Chain):
            if self.id == ' ' and other.id != ' ':
                return 0
            elif self.id != ' ' and other.id == ' ':
                return 1
            else:
                return self.id <= other.id
        else:
            return NotImplemented

    def _translate_id(self, id):
        """Translate sequence identifer to tuple form (PRIVATE).

        A residue id is normally a tuple (hetero flag, sequence identifier,
        insertion code). Since for most residues the hetero flag and the
        insertion code are blank (i.e. " "), you can just use the sequence
        identifier to index a residue in a chain. The _translate_id method
        translates the sequence identifier to the (" ", sequence identifier,
        " ") tuple.

        Arguments:

        - id - int, residue resseq

        """
        if isinstance(id, int):
            id = (' ', id, ' ')
        return id

    def __getitem__(self, id):
        """Return the residue with given id.

        The id of a residue is (hetero flag, sequence identifier, insertion code).
        If id is an int, it is translated to (" ", id, " ") by the _translate_id
        method.

        Arguments:

        - id - (string, int, string) or int

        """
        id = self._translate_id(id)
        return Entity.__getitem__(self, id)

    def __contains__(self, id):
        """Check if a residue with given id is present in this chain.

        Arguments:

        - id - (string, int, string) or int

        """
        id = self._translate_id(id)
        return Entity.__contains__(self, id)

    def __delitem__(self, id):
        """Delete item.

        Arguments:

        - id - (string, int, string) or int

        """
        id = self._translate_id(id)
        return Entity.__delitem__(self, id)

    def __repr__(self):
        return "<Chain id=%s>" % self.get_id()

    # Public methods

    def get_unpacked_list(self):
        """Return a list of undisordered residues.

        Some Residue objects hide several disordered residues
        (DisorderedResidue objects). This method unpacks them,
        ie. it returns a list of simple Residue objects.
        """
        unpacked_list = []
        for residue in self.get_list():
            if residue.is_disordered() == 2:
                for dresidue in residue.disordered_get_list():
                    unpacked_list.append(dresidue)
            else:
                unpacked_list.append(residue)
        return unpacked_list

    def has_id(self, id):
        """Return 1 if a residue with given id is present.

        The id of a residue is (hetero flag, sequence identifier, insertion code).

        If id is an int, it is translated to (" ", id, " ") by the _translate_id
        method.

        Arguments:

        - id - (string, int, string) or int

        """
        id = self._translate_id(id)
        return Entity.has_id(self, id)

    # Public

    def get_residues(self):
        for r in self:
            yield r

    def get_atoms(self):
        for r in self.get_residues():
            for a in r:
                yield a
