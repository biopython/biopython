from UserList import UserList

from Bio.SCOP.Raf import to_one_letter_code
from Bio.PDB.PathFinder import PathFinder
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBExceptions import PDBException


class Polypeptide(UserList):
	def __init__(self, residue_list):
		UserList.__init__(self, residue_list)

	def __repr__(self):
		# get chain, resseq for start residue
		s, m, c, r=self[0].get_full_id()
		if c==" ":
			c=""
		start=c+str(r[1])
		# get chain, resseq for end residue
		s, m, c, r=self[-1].get_full_id()
		if c==" ":
			c=""
		end=c+str(r[1])
		return "<Polypeptide start=%s end=%s length=%s>" % (start, end, len(self))


def build_peptides(structure, model_id=0, aa_only=1, radius=1.8):
	"""
	Extract all polypeptides from a structure.

	Arguments:
	o structure - structure object
	o model_id - polypeptides are extracted from this model (default 0)
	o aa_only - only consider standard AA if this flag is set (default 1) 
	o radius - maximum bond length for a C-N bond (default 1.8)

	Returns a list of lists of residues. Each list of residues
	represents a polypeptide.
	"""
	model=structure[model_id]
	atom_list=[]
	# get all N and C atoms from this model
	for chain in model.get_list():
		for residue in chain.get_list():
			if aa_only:
				# only use standard AA
				resname=residue.get_resname()
				if not to_one_letter_code.has_key(resname):
					# not a standard AA so skip
					continue
			# get N and C atoms from residue
			for atom_name in ["C", "N"]:
				if residue.has_id(atom_name):
					a=residue[atom_name]
					if not a.is_disordered():
						# skip disordered N or C atoms
						atom_list.append(a)
	if len(atom_list)>2:
		ns=NeighborSearch(atom_list)
		# all atom pairs within radius
		# should be (C,N) but can be (N,N) and (C,C) 
		# if the structure is funny
		neighbors=ns.search_all(radius)
	else:
		# not a single (N,C) pair found
		return []
	# Now find connected stretches, ie. polypeptides
	pf=PathFinder()
	for a, b in neighbors:
		# get atom names
		na=a.get_id()
		nb=b.get_id()
		# get residues
		residue1=a.get_parent()
		residue2=b.get_parent()
		# check if it is indeed a peptide bond
		# residue at N-terminus should come first!
		try:
			if na=="C" and nb=="N":
				pf.add_edge(residue1, residue2)
			elif na=="N" and nb=="C":
				pf.add_edge(residue2, residue1)
		except PDBException:
			print "WARNING: some strange peptide bonds are present"
	# return a list of lists of residues
	# each list of residues represents a polypeptide
	polypeptide_list=[]
	# create Polypeptide objects and put them in a list
	for residue_list in pf.get_paths():
		polypeptide=Polypeptide(residue_list)
		polypeptide_list.append(polypeptide)
	return polypeptide_list
		

if __name__=="__main__":

	import sys

	from Bio.PDB.PDBParser import PDBParser

	p=PDBParser(PERMISSIVE=1)

	for file in sys.argv[1:]:
		print "Extracting polypeptides from ", file
		s=p.get_structure("scr", file)
		for pp in build_peptides(s):
			print pp
		print




