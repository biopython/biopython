# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Class that maps (chain_id, residue_id) to a residue property."""

from __future__ import print_function


class AbstractPropertyMap(object):
    def __init__(self, property_dict, property_keys, property_list):
        self.property_dict=property_dict
        self.property_keys=property_keys
        self.property_list=property_list

    def _translate_id(self, entity_id):
        return entity_id

    def __contains__(self, id):
        """True if the mapping has a property for this residue.

        Example:
            >>> if (chain_id, res_id) in apmap:
            ...     res, prop = apmap[(chain_id, res_id)]

        @param chain_id: chain id
        @type chain_id: char

        @param res_id: residue id
        @type res_id: char
        """
        translated_id = self._translate_id(id)
        return (translated_id in self.property_dict)

    def __getitem__(self, key):
        """
        Return property for a residue.

        @param chain_id: chain id
        @type chain_id: char

        @param res_id: residue id
        @type res_id: int or (char, int, char)

        @return: some residue property
        @rtype: anything (can be a tuple)
        """
        translated_id=self._translate_id(key)
        return self.property_dict[translated_id]

    def __len__(self):
        """
        Return number of residues for which the property is available.

        @return: number of residues
        @rtype: int
        """
        return len(self.property_dict)

    def has_key(self, id):
        """True if the mapping has a property for this residue.

        (Obsolete; use "id in mapping" instead.)

        Example:

            >>> if apmap.has_key((chain_id, res_id)):
            ...     res, prop = apmap[(chain_id, res_id)]

        Is equivalent to:

            >>> if (chain_id, res_id) in apmap:
            ...     res, prop = apmap[(chain_id, res_id)]

        @param chain_id: chain id
        @type chain_id: char

        @param res_id: residue id
        @type res_id: char
        """
        import warnings
        from Bio import BiopythonDeprecationWarning
        warnings.warn("This function is deprecated; use 'id in mapping' instead", BiopythonDeprecationWarning)
        return (id in self)

    def keys(self):
        """
        Return the list of residues.

        @return: list of residues for which the property was calculated
        @rtype: [(chain_id, res_id), (chain_id, res_id),...]
        """
        return self.property_keys

    def __iter__(self):
        """
        Iterate over the (entity, property) list. Handy alternative to
        the dictionary-like access.

        Example:
            >>> for (res, property) in iter(map):
            ...     print(res, property)

        @return: iterator
        """
        for i in range(0, len(self.property_list)):
            yield self.property_list[i]


class AbstractResiduePropertyMap(AbstractPropertyMap):
    def __init__(self, property_dict, property_keys, property_list):
        AbstractPropertyMap.__init__(self, property_dict, property_keys,
                property_list)

    def _translate_id(self, ent_id):
        chain_id, res_id=ent_id
        if isinstance(res_id, int):
            ent_id=(chain_id, (' ', res_id, ' '))
        return ent_id


class AbstractAtomPropertyMap(AbstractPropertyMap):
    def __init__(self, property_dict, property_keys, property_list):
        AbstractPropertyMap.__init__(self, property_dict, property_keys,
                property_list)

    def _translate_id(self, ent_id):
        if len(ent_id)==4:
            chain_id, res_id, atom_name, icode=ent_id
        else:
            chain_id, res_id, atom_name=ent_id
            icode=None
        if isinstance(res_id, int):
            ent_id=(chain_id, (' ', res_id, ' '), atom_name, icode)
        return ent_id
