import re

from Bio.PDB.PDBIO import PDBIO
from Bio.SVDSuperimposer import SVDSuperimposer


_hydrogen=re.compile("[123 ]*H.*")


class Selector:
	"""
	Only accepts residues with right chainid
	and between start and end. Remove hydrogens.
        Only use model 0.
	"""
	def __init__(self, chainid, start, end):
		self.chainid=chainid
		self.start=start
		self.end=end

	def accept_model(self, model):
		# model - only keep model 0
		if model.get_id()==0:
			return 1
		return 0

	def accept_chain(self, chain):
		if chain.get_id()==self.chainid:
			return 1
		return 0

	def accept_residue(self, residue):
		# residue - between start and end
		hetatm_flag, resseq, icode=residue.get_id()
		if hetatm_flag!=" ":
			# skip HETATMS
			return 0
		if icode!=" ":
			print "WARNING: Icode at ", residue.get_id()
		if self.start<=resseq<=self.end:
			return 1
		return 0

	def accept_atom(self, atom):
		# atoms - get rid of hydrogens
		name=atom.get_id()
		if _hydrogen.match(name):
			return 0
		else:
			return 1


def extract(structure, chainid, start, end, filename):
	"""
	Write out selected portion to filename.
	"""
	sel=Selector(chainid, start, end)
	io=PDBIO()
	io.set_structure(structure)
	io.save(filename, sel)



if __name__=="__main__":

	from Bio.PDB.PDBParser import PDBParser

	import sys

	p=PDBParser()
	s=p.get_structure("scr", sys.argv[1])

	extract(s, "A", 1, 100, "out.pdb")
	
