# Copyright (C) 2002, Thomas Hamelryck (thamelry@vub.ac.be)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.  


from Bio.PDB.PDBParser import PDBParser 
from Bio.PDB.Polypeptide import CaPPBuilder, PPBuilder


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
for pp in ppb.build_peptides(structure, 1):
	print pp

print "Polypeptides using CA-CA"
ppb=CaPPBuilder()
for pp in ppb.build_peptides(structure, 1):
	print pp

					
