# Copyright (C) 2002, Thomas Hamelryck (thamelry@vub.ac.be)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.  

try:
    import numpy
except ImportError:
    from Bio import MissingExternalDependencyError
    raise MissingExternalDependencyError(\
        "Install NumPy if you want to use Bio.PDB.")

from Bio.PDB.PDBExceptions import PDBConstructionException, PDBConstructionWarning

import warnings
def send_pdb_warnings_to_stdout(message, category, filename, lineno,
                                file=None, line=None):
    if category in [PDBConstructionException, PDBConstructionWarning]:
        print message
warnings.resetwarnings()
warnings.showwarning = send_pdb_warnings_to_stdout

def quick_neighbor_search_test():
    #Based on the self test in Bio.PDB.NeighborSearch
    from numpy.random import random
    from Bio.PDB.NeighborSearch import NeighborSearch

    class Atom:
        def __init__(self):
            self.coord=(100*random(3))

        def get_coord(self):
            return self.coord

    for i in range(0, 20):
        al = [Atom() for j in range(100)]
        ns=NeighborSearch(al)
        hits = ns.search_all(5.0)
        assert hits >= 0
    print "Done"    


class TheVoid:
    # Class to hide stderr output
    def write(self, string):
        pass

def run_test():
    from Bio.PDB import PDBParser, PPBuilder, CaPPBuilder


    # first make a PDB parser object
    p=PDBParser(PERMISSIVE=1) 

    # get the structure, call it "example"
    structure=p.get_structure("example", "PDB/a_structure.pdb")

    # now loop over content and print some info
    for model in structure.get_list():
        model_id=model.get_id()
        print "Model %i contains %i chains." % (model_id, len(model))
        for chain in model.get_list():
            chain_id=chain.get_id()
            print "\tChain '%s' contains %i residues." % (chain_id, len(chain))
            for residue in chain.get_list():
                residue_id=residue.get_id()
                hetfield, resseq, icode=residue_id
                print "\t\tResidue ('%s', %i, '%s') contains %i atoms." % (hetfield, resseq, icode, len(residue))
                # check if there is disorder due to a point mutation --- this is rare
                if residue.is_disordered()==2:
                    print "\t\t\tThere is a point mutation present in the crystal at this position."
                    s="\t\t\tResidues at this position are "
                    for resname in residue.disordered_get_id_list():
                        s=s+resname+" "
                    print s[:-1]+"."
                # count the number of disordered atoms
                if residue.is_disordered()==1:
                    disordered_count=0
                    for atom in residue.get_list():
                        if atom.is_disordered():
                            disordered_count=disordered_count+1
                    if disordered_count>0:
                        print "\t\t\tThe residue contains %i disordered atoms." % disordered_count


    print "Polypeptides using C-N"
    ppb=PPBuilder()
    for pp in ppb.build_peptides(structure[1]):
        print pp

    print "Polypeptides using CA-CA"
    ppb=CaPPBuilder()
    for pp in ppb.build_peptides(structure[1]):
        print pp

    print "NeighborSearch test"
    quick_neighbor_search_test()

run_test()

warnings.resetwarnings() #clean up in case this affects other unit tests
