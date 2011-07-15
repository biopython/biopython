# Copyright (C) 2011, Mikael Trellet (mikael.trellet@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Interface class, used in Structure objects."""

from Bio.PDB.Entity import Entity
from Bio.Data import IUPACData
from Bio.SCOP.Raf import to_one_letter_code


class Interface(Entity):
    """
    The Interface object isn't automatically initialize during a PDB parsing,
    but can be guessed from an existing parsed structure in order to analyse
    the interface between 2 or more chains in a complex.
    """

    def __init__(self, id):
        self.level="I"
        self.id=id
        self.neighbors = {}
        self.uniq_pairs = []

        Entity.__init__(self, id)

    # Override Entity add method
    # Interface doesnt follow strictly
    # other Entity rules.
    #
    # Its childs are residues
    # but it may be useful
    # to list them by chain.
    
    def add(self, entity):
        "Add a child to the Entity."

        entity_id=entity.get_id()
        if not self.has_id(entity_id):

            self.child_list.append(entity)
            if entity.parent.id not in self.child_dict:
                self.child_dict[entity.parent.id] = []
            self.child_dict[entity.parent.id].append(entity)

    def get_chains(self):
        "Get the different chains involved in the Interface object"
        for chain in self.child_dict.keys():
            yield chain

    def set_neighbors(self):
        "Creates residues list of neighbors"
        ## Initializes neighbors dictionnary with interface chains
        for c in self.get_chains():
            self.neighbors[c]={}
            
        for resA, resB in self.uniq_pairs:
        ## Checking for 1st residue (if his chain exist, then if 
        ## it is referenced and finally if his partner is already present)
            if resA not in self.neighbors[resA.parent.id]:
                self.neighbors[resA.parent.id][resA]=[]
                self.neighbors[resA.parent.id][resA].append(resB)
            elif resB not in self.neighbors[resA.parent.id][resA]:
                self.neighbors[resA.parent.id][resA].append(resB)
        ## Checking for 2nd residue
            if resB not in self.neighbors[resB.parent.id]:
                self.neighbors[resB.parent.id][resB]=[]
                self.neighbors[resB.parent.id][resB].append(resB)
            elif resA not in self.neighbors[resB.parent.id][resB]:
                self.neighbors[resB.parent.id][resB].append(resA)
        neighbors=self.neighbors
        return neighbors

    def get_polar_percentage(self):
        "Gets the percentage of polar residues in the interface"
        
        polar_list=getattr(IUPACData, "protein_polarity")
        polar_residues = [r for r in self if to_one_letter_code[r.resname] in polar_list['polar']]
        polar_percentage=float(len(polar_residues))/len(self)
        return polar_percentage

