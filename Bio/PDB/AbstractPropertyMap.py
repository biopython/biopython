# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Class that maps (chain_id, residue_id) to a residue property."""


class AbstractPropertyMap:
    """Define base class, map holder of residue properties."""

    def __init__(self, property_dict, property_keys, property_list):
        """Initialize the class."""
        self.property_dict = property_dict
        self.property_keys = property_keys
        self.property_list = property_list

    def _translate_id(self, entity_id):
        """Return entity identifier (PRIVATE)."""
        return entity_id

    def __contains__(self, id):
        """Check if the mapping has a property for this residue.

        :param chain_id: chain id
        :type chain_id: char

        :param res_id: residue id
        :type res_id: char

        Examples
        --------
        This is an incomplete but illustrative example::

            if (chain_id, res_id) in apmap:
                res, prop = apmap[(chain_id, res_id)]

        """
        translated_id = self._translate_id(id)
        return translated_id in self.property_dict

    def __getitem__(self, key):
        """Return property for a residue.

        :param chain_id: chain id
        :type chain_id: char

        :param res_id: residue id
        :type res_id: int or (char, int, char)

        :return: some residue property
        :rtype: anything (can be a tuple)
        """
        translated_id = self._translate_id(key)
        return self.property_dict[translated_id]

    def __len__(self):
        """Return number of residues for which the property is available.

        :return: number of residues
        :rtype: int
        """
        return len(self.property_dict)

    def keys(self):
        """Return the list of residues.

        :return: list of residues for which the property was calculated
        :rtype: [(chain_id, res_id), (chain_id, res_id),...]
        """
        return self.property_keys

    def __iter__(self):
        """Iterate over the (entity, property) list.

        Handy alternative to the dictionary-like access.

        :return: iterator

        Examples
        --------
        >>> entity_property_list = [
        ...     ('entity_1', 'property_1'),
        ...     ('entity_2', 'property_2')
        ... ]
        >>> map = AbstractPropertyMap({}, [], entity_property_list)
        >>> for (res, property) in iter(map):
        ...     print(res, property)
        entity_1 property_1
        entity_2 property_2

        """
        for i in range(len(self.property_list)):
            yield self.property_list[i]


class AbstractResiduePropertyMap(AbstractPropertyMap):
    """Define class for residue properties map."""

    def __init__(self, property_dict, property_keys, property_list):
        """Initialize the class."""
        AbstractPropertyMap.__init__(self, property_dict, property_keys, property_list)

    def _translate_id(self, ent_id):
        """Return entity identifier on residue (PRIVATE)."""
        chain_id, res_id = ent_id
        if isinstance(res_id, int):
            ent_id = (chain_id, (" ", res_id, " "))
        return ent_id


class AbstractAtomPropertyMap(AbstractPropertyMap):
    """Define class for atom properties map."""

    def __init__(self, property_dict, property_keys, property_list):
        """Initialize the class."""
        AbstractPropertyMap.__init__(self, property_dict, property_keys, property_list)

    def _translate_id(self, ent_id):
        """Return entity identifier on atoms (PRIVATE)."""
        if len(ent_id) == 4:
            chain_id, res_id, atom_name, icode = ent_id
        else:
            chain_id, res_id, atom_name = ent_id
            icode = None
        if isinstance(res_id, int):
            ent_id = (chain_id, (" ", res_id, " "), atom_name, icode)
        return ent_id
