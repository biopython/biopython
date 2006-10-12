# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from types import IntType

__doc__="Class that maps (chain_id, residue_id) to a residue property"


class AbstractPropertyMap:
    def __init__(self, property_dict, property_keys, property_list):
        self.property_dict=property_dict
        self.property_keys=property_keys
        self.property_list=property_list

    def _translate_id(self, entity_id):
        return entity_id

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
        """
        Return 1 if the map has a property for this residue, 0 otherwise.

        Example:
            >>> if map.has_key((chain_id, res_id)):
            >>>     res, property=map[(chain_id, res_id)]

        @param chain_id: chain id
        @type chain_id: char 

        @param res_id: residue id
        @type res_id: char 
        """
        translated_id=self._translate_id(id)
        return self.property_dict.has_key(translated_id)

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
            >>>     print res, property

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
        if type(res_id)==IntType:
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
        if type(res_id)==IntType:
            ent_id=(chain_id, (' ', res_id, ' '), atom_name, icode)
        return ent_id



