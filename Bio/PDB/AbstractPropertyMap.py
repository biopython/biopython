
__doc__="Class that maps (chain_id, residue_id) to a residue property"


class AbstractPropertyMap:
    def __init__(self, property_dict, property_keys, property_list):
        self.property_dict=property_dict
        self.property_keys=property_keys
        self.property_list=property_list

    def _translate_res_id(self, res_id):
        if isinstance(res_id, int):
            res_id=(' ', res_id, ' ')
        return res_id

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
        chain_id, res_id=key
        res_id=self._translate_res_id(res_id)
        return self.property_dict[(chain_id, res_id)]

    def __len__(self):
        """
        Return number of residues for which the property is available.

        @return: number of residues
        @rtype: int
        """
        return len(self.property_dict)

    def has_key(self, chain_id, res_id):
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
        res_id=self._translate_res_id(res_id)
        return self.property_dict.has_key((chain_id, res_id))

    def keys(self):
        """
        Return the list of residues.

        @return: list of residues for which the property was calculated
        @rtype: [(chain_id, res_id), (chain_id, res_id),...] 
        """
        return self.property_keys

    def __iter__(self):
        """
        Iterate over the (residue, property) list. Handy alternative to the dictionary-like 
        access.

        Example:
            >>> for (res, property) in iter(map):
            >>>     print res, property

        @return: iterator
        """
        for i in range(0, len(self.property_list)):
            yield self.property_list[i]

