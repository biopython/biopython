# Copyright (C) 2011, Mikael Trellet (mikael.trellet@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""InterfaceBuilder class, used in Model objects."""

from Bio.PDB.Interface import Interface
from Bio.PDB import Selection
from Bio.PDB import NeighborSearch
from Bio.PDB.NACCESS import NACCESS
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
from Bio.PDB.DSSP import DSSP
import os

class InterfaceBuilder(object):
    """
    Deals with contructing the Interface object. The InterfaceBuilder class is used
    by the Model classe to get an interface from a model.
    """
    def __init__(self, model, id=None, threshold=5.0, rsa_calculation=False, rsa_threshold=1.0, include_waters=False, *chains):

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
        self._build_interface(model, id, threshold, rsa_calculation, rsa_threshold, include_waters, *chains)

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
        
    def _secondary_structure(self, model):
        """Define secondary structure of the interface"""
        
        ss=[]
        #Launch DSSP calculation on the whole model (DSSP limitations)
        dssp=DSSP(model, dssp='/home/mika/dssp/dsspcmbi')
        res_list = [r for r in model.get_residues() if r.id[0] == ' ']
        
        for elt in dssp:
            print 'ELT', elt[1]
            if elt[0] in self.interface:
                ss.append(elt[1])
                
        ss.sort()
        
        helical=(ss.count('H')+ss.count('G')+ss.count('I'))/len(ss)
        strand=(ss.count('E'))/len(ss)
        coil=(len(ss)-helical-strand)/len(ss)
        
        self.interface.secondary_structure=[helical,strand,coil]
        
    def _get_residue_SASA(self, structure):
        """Retrieves residual SASA from a NACCESS-submitted structure object"""

        # Ignore Hydrogens otherwise default NACCESS freaks out
        # Maybe add support for flags in the NACCESS module?
        # From the readme:
        # "By default, the program ignores HETATM records, hydrogens, and waters. If you
        # want these to be considered in the calculation supply a parameter of the form
        # -h, -w and/or -y respectively."

        sasa_l = [res.xtra['EXP_NACCESS']['all_atoms_rel'] for res in structure.get_residues() if res.id[0] == " "]
        sasa = sum(sasa_l)
        return sasa
    
    def _rsa_calculation(self, model, chain_list, rsa_threshold):
        "Uses NACCESS module in order to calculate the Buried Surface Area"
        pairs=[]
        # Create temporary structures to feed NACCESS
        structure_A=Structure("chainA")
        structure_B=Structure("chainB")
        mA = Model(0)
        mB = Model(0)
        mA.add(model[chain_list[0]])
        mB.add(model[chain_list[1]])
        structure_A.add(mA)
        structure_B.add(mB)
        # Calculate SASAs
        nacc_at=NACCESS(model)
        model_values=[]
                
        res_list = [r for r in model.get_residues() if r.id[0] == ' ']
        structure_A_reslist =[r for r in structure_A[0].get_residues() if r.id[0] == ' ']
        structure_B_reslist =[r for r in structure_B[0].get_residues() if r.id[0] == ' ']
        
        for res in res_list:
            model_values.append(float(res.xtra['EXP_NACCESS']['all_atoms_rel']))
            
                
        sas_tot= self._get_residue_SASA(model)
        #print 'Accessible surface area, complex:', sas_tot

        nacc_at=NACCESS(structure_A[0])
        nacc_at=NACCESS(structure_B[0])
        submodel_values=[]
                
        for res in structure_A_reslist:
            if res.id[0]==' ':
                submodel_values.append(float(res.xtra['EXP_NACCESS']['all_atoms_rel']))                
                
        for res in structure_B_reslist:
            if res.id[0]==' ':
                submodel_values.append(float(res.xtra['EXP_NACCESS']['all_atoms_rel']))
        
        count=0        
        for res in res_list:
            if res in structure_A_reslist and ((submodel_values[count] - model_values[count]) > rsa_threshold):
                pairs.append(res)
            elif res in structure_B_reslist and ((submodel_values[count] - model_values[count]) > rsa_threshold):
                pairs.append(res)
            count=count+1
        
        
        sas_A= self._get_residue_SASA(structure_A)
        #print 'Accessible surface aream CHAIN A :', sas_A
        sas_B= self._get_residue_SASA(structure_B)
        #print 'Accessible surface aream CHAIN B :',sas_B
        
        # Calculate BSA
        bsa = sas_A+sas_B-sas_tot
                
        self.interface.accessibility=[bsa, sas_A, sas_B, sas_tot]
        
        return pairs
    
    def _build_interface(self, model, id, threshold, rsa_calculation, rsa_threshold, include_waters=False, *chains):
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
                
        # Accessible surface area calculated for each residue
        # if naccess setup on user computer and rsa_calculation
        # argument is TRUE
        if rsa_calculation and os.system('which naccess') == 0:
            rsa_pairs=self._rsa_calculation(model, chain_list, rsa_threshold)
            
        for res in rsa_pairs:
            if res not in self.interface:
                self._add_residue(res)
        self._secondary_structure(model)
        #interface=uniq_pairs
        self.interface.uniq_pairs=uniq_pairs
        # Add neighbors
        # so you can call
        # my_int.neighbors['A'][10] = [list of contacting residues]

#   Public classes
