ATOM_FORMAT_STRING="%s%5i %-4s%c%3s %c%4i%c   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s\n"


class Select:
	"""
	Dummy class.
	This selects which entities will be written out.
	"""
	def accept_model(self, model):
		return 1

	def accept_chain(self, chain):
		return 1

	def accept_residue(self, residue):
		return 1

	def accept_atom(self, atom):
		return 1


class PDBIO:
	def __init__(self):
		pass
	
	# private mathods

	def _get_atom_line(self, atom, hetfield, segid, atom_number, resname, 
		resseq, icode, chain_id, element="  ", charge="  "):
		"""
		Returns an ATOM PDB string.
		"""
		if hetfield!=" ":
			record_type="HETATM"
		else:
			record_type="ATOM  "
		name=atom.get_fullname()
		altloc=atom.get_altloc()
		x, y, z=atom.get_coord()
		bfactor=atom.get_bfactor()
		occupancy=atom.get_occupancy()
		args=(record_type, atom_number, name, altloc, resname, chain_id,
			resseq, icode, x, y, z, occupancy, bfactor, segid,
			element, charge)
		return ATOM_FORMAT_STRING % args

	# Public methods

	def set_structure(self, structure):
		self.structure=structure

	def save(self, filename, select=Select()):
		"""
		select --- selects which entities will be written.
			Should have the following methods:
				accept_model(model)
				accept_chain(chain)
				accept_residue(residue)
				accept_atom(atom)
			These methods should return 1 if the entity
			is to be written out, 0 otherwise.
		"""
		get_atom_line=self._get_atom_line
		fp=open(filename, "w")
		# multiple models?
		if len(self.structure)>1:
			model_flag=1
		else:
			model_flag=0
		for model in self.structure.get_list():
			if not select.accept_model(model):
				continue
			# necessary for ENDMDL 
			# do not write ENDMDL if no residues were written
			# for this model
			model_residues_written=0
			atom_number=1
			if model_flag:
				fp.write("MODEL \n")
			for chain in model.get_list():
				if not select.accept_chain(chain):
					continue
				chain_id=chain.get_id()
				# necessary for TER 
				# do not write TER if no residues were written
				# for this chain
				chain_residues_written=0
				for residue in chain.get_unpacked_list():
					if not select.accept_residue(residue):
						continue
					hetfield, resseq, icode=residue.get_id()
					resname=residue.get_resname()  
					segid=residue.get_segid()
					for atom in residue.get_unpacked_list():
						if select.accept_atom(atom):
							chain_residues_written=1
							model_residues_written=1
							s=get_atom_line(atom, hetfield, segid, atom_number, resname,
								resseq, icode, chain_id)
							fp.write(s)
							atom_number=atom_number+1
				if chain_residues_written:
					fp.write("TER\n")
			if model_flag and model_residues_written:
				fp.write("ENDMDL\n")
		fp.close()

if __name__=="__main__":
	
	from Bio.PDB.PDBParser import PDBParser

	import sys

	p=PDBParser(PERMISSIVE=1)

	s=p.get_structure("test", sys.argv[1])

	io=PDBIO()

	io.set_structure(s)

	io.save("out.pdb")





