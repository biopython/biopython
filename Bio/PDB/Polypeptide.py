from Bio.SCOP.Raf import to_one_letter_code
from Bio.PDB.PathFinder import PathFinder
from Bio.PDB.NeighborSearch import NeighborSearch


def build_peptides(structure, model_id=0, aa_only=1, radius=2.0):
	model=structure[model_id]
	atom_list=[]
	for chain in model.get_list():
		for residue in chain.get_list():
			if aa_only:
				resname=residue.get_resname()
				if not to_one_letter_code.has_key(resname):
					continue
			for atom_name in ["C", "N"]:
				if residue.has_id(atom_name):
					a=residue[atom_name]
					if not a.is_disordered():
						atom_list.append(a)
	if len(atom_list)>2:
		ns=NeighborSearch(atom_list, 2)
		neighbors=ns.search_all(radius)
	else:
		return []
	pf=PathFinder()
	for a, b in neighbors:
		na=a.get_id()
		nb=b.get_id()
		residue1=a.get_parent()
		residue2=b.get_parent()
		if na=="C" and nb=="N":
			pf.add_edge(residue1, residue2)
		elif na=="N" and nb=="C":
			pf.add_edge(residue2, residue1)
	return pf.get_paths()

if __name__=="__main__":

	from Bio.PDB.PDBParser import PDBParser

	p=PDBParser()
	s=p.get_structure("scr", "/home/tham/data/pdb/SMALL.pdb")

	for l in build_peptides(s):
		for r in l:
			print r
		print




