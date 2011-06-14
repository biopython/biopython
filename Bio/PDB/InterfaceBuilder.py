# Copyright (C) 2011, Mikael Trellet (mikael.trellet@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""InterfaceBuilder class, used in Model objects."""

from Bio.PDB.Interface import Interface
from Bio.PDB import Selection
from Bio.PDB import NeighborSearch

class InterfaceBuilder(object):
    """
    Deals with contructing the Interface object. The InterfaceBuilder class is used
    by the Model classe to get an interface from a model.
    """
    def __init__(self, model, id=None, threshold=5.0, include_waters=False, *chains):

        chain_list = []

        # Unpack chain list
        if not chains:
            chain_list = [c.id for c in model]
        else:
             chain_list = self._unpack_chains(chains)
        self.chain_list = sorted(chain_list)

        if not id: # Build name like interface_chainAchainB...chainX from chains
            id = "Interface_%s" %(''.join(self.chain_list))

        self.interface = Interface(id)
        self._build_interface(model, id, threshold, include_waters, *chains)

    def _unpack_chains(self, list_of_tuples):
        """Unpacks a list of tuples into a list of characters"""

        chain_list = []
        chains = set(list_of_tuples)
        for user_chain in chains:
            if user_chain[0] not in chain_list:
                chain_list.append(user_chain[0])
            if user_chain[1] not in chain_list:
                chain_list.append(user_chain[1])

        return chain_list

    def get_interface(self):
        return self.interface

    def _add_residue(self, residue):
        """Adds a residue to an Interface object"""

        self.interface.add(residue)

    def _build_interface(self, model, id, threshold, include_waters=False, *chains):
        """
        Return the interface of a model
        """

        self.threshold=threshold

        # Recover chain list from initial unpacking
        chain_list = self.chain_list

        # Unfold atom list
        atom_list = []
        for c in model:
            if c.id in chain_list:
                atom_list.extend(Selection.unfold_entities(c,'A'))

        # Using of NeighborSearch class in order to get the list of all residues at least than
        # the threshold distance of each others
        ns=NeighborSearch(atom_list)
        pairs=ns.search_all(threshold, 'R')

        if not pairs:
            raise ValueError("No atoms found in the interface")

        # Selection of residues pairs
        # 1. Exclude water contacts
        # 2. Filter same-chain contacts
        # 3. Filter user-defined chain pairs

        uniq_pairs=[]

        for pair in pairs:
             
            pair_resnames = (pair[0].resname, pair[1].resname)
            pair_chains = (pair[0].parent.id, pair[1].parent.id)

            if (not include_waters and 'HOH' in pair_resnames) or (pair_chains[0] == pair_chains[1]):
                continue

            if not (chains and not (pair_chains in chains)):
                uniq_pairs.append(pair)

        # Build the Interface
        # 1. Iterate over the pair list
        # 2. Add residues.

        for resA, resB in uniq_pairs:
            if resA not in self.interface:
                self._add_residue(resA)
            if resB not in self.interface:
                self._add_residue(resB)

        #interface=uniq_pairs
        self.interface.uniq_pairs=uniq_pairs
        # Add neighbors
        # so you can call
        # my_int.neighbors['A'][10] = [list of contacting residues]

#   Public classes
