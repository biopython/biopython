ATOM_FORMAT_STRING="%s%5i %-4s%c%3s %c%4i%c   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s\n"


class PDBIO:
	def __init__(self):
		pass

	# Private methods

	def _get_atom_line(self, atom, hetfield, segid, atom_number, resname, 
				resseq, icode, chain_id, element="  ", charge="  "):
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

	def save(self, filename):
		fp=open(filename, "w")
		get_atom_line=self._get_atom_line
		# multiple models?
		if len(self.structure)>1:
			model_flag=1
		else:
			model_flag=0
		for model in self.structure.get_list():
			atom_number=1
			if model_flag:
				fp.write("MODEL \n")
			for chain in model.get_list():
				chain_id=chain.get_id()
				for residue in chain.get_unpacked_list():
					hetfield, resseq, icode=residue.get_id()
					resname=residue.get_resname()  
					segid=residue.get_segid()
					for atom in residue.get_unpacked_list():
						s=get_atom_line(atom, hetfield, segid, atom_number, resname,
								resseq, icode, chain_id)
						fp.write(s)
						atom_number=atom_number+1
				fp.write("TER")
			if model_flag:
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





