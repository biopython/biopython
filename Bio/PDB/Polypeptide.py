from Bio.SCOP.Raf import to_one_letter_code
from Bio.PDB.PDBExceptions import PDBException
from Numeric import sum


def dist_sq(p, q):
	"Return squared distance between coordinates."
	d=p-q
	return sum(d*d)


class Polypeptide(list):
	"""
	A polypeptide is simply a list of Residue objects.
	"""
	def __repr__(self):
		start=self[0].get_id()[1]
		end=self[-1].get_id()[1]
		s="<Polypeptide start=%s end=%s>" % (start, end)
		return s


class _PPBuilder:
	"""
	Base class to extract polypeptides.
	It checks if two consecutive residues in a chain 
	are connected. The connectivity test is implemented by a 
	subclass.
	"""
	def __init__(self, radius_sq):
		self.radius_sq=radius_sq

	def _accept(self, residue):
		"Check if the residue is an amino acid."
		resname=residue.get_resname()
		if to_one_letter_code.has_key(resname):
			# not a standard AA so skip
			return 1
		else:
			return 0
	
	def build_peptides(self, structure, model_id=0, aa_only=1):
		"""
		Build and return a list of Polypeptide objects.
		model_id --- only this model is examined
		aa_only --- if 1, the residue name needs to be standard AA
		"""
		is_connected=self._is_connected
		accept=self._accept
		model=structure[model_id]
		pp_list=[]
		for chain in model.get_list():
			prev=None
			pp=None
			for next in chain.get_list():
				if prev:
					if is_connected(prev, next):
						if pp is None:
							pp=Polypeptide()
							pp.append(prev)
							pp_list.append(pp)
						pp.append(next)
					else:
						pp=None
				if aa_only:
					if accept(next):
						prev=next
					else:
						prev=None
						pp=None
				else:
					prev=next
		return pp_list


class CaPPBuilder(_PPBuilder):
	"""
	Use CA--CA distance to find polypeptides.
	"""
	def __init__(self, radius=4.3):
		_PPBuilder.__init__(self, radius*radius)

	def _is_connected(self, prev, next):
		for r in [prev, next]:
			if not r.has_id("CA"):
				return 0
		n_ca=next["CA"]
		p_ca=prev["CA"]
		for a in [n_ca, p_ca]:
			# if a CA is disordered we consider
			# it as not part of a polypeptide 
			# chain
			if a.is_disordered():
				return 0
		nc=n_ca.get_coord()
		pc=p_ca.get_coord()
		if dist_sq(nc, pc)<self.radius_sq:
			return 1
		else:
			return 0


class PPBuilder(_PPBuilder):
	"""
	Use C--N distance to find polypeptides.
	"""
	def __init__(self, radius=1.8):
		_PPBuilder.__init__(self, radius*radius)

	def _is_connected(self, prev, next):
		if not prev.has_id("C"):
			return 0
		if not next.has_id("N"):
			return 0
		test_dist=self._test_dist
		c=prev["C"]
		n=next["N"]
		# Test all disordered atom positions!
		if c.is_disordered():
			clist=c.disordered_get_list()
		else:
			clist=[c]
		if n.is_disordered():
			nlist=n.disordered_get_list()
		else:
			nlist=[n]
		for nn in nlist:
			for cc in clist:
				# To form a peptide bond, N and C must be 
				# within radius and have the same altloc
				# identifier
				n_altloc=nn.get_altloc()
				c_altloc=cc.get_altloc()
				if n_altloc==c_altloc: 
					if test_dist(nn, cc):
						# Select the disordered atoms that
						# are indeed bonded
						if c.is_disordered():
							c.disordered_select(c_altloc)
						if n.is_disordered():
							n.disordered_select(n_altloc)
						return 1
		return 0

	def _test_dist(self, c, n):
		"Return 1 if distance between atoms<radius"
		cc=c.get_coord()
		nc=n.get_coord()
		#print dist_sq(cc, nc), self.radius_sq
		if dist_sq(cc, nc)<self.radius_sq:
			return 1
		else:
			return 0
	

if __name__=="__main__":

	import sys

	from Bio.PDB.PDBParser import PDBParser

	p=PDBParser(PERMISSIVE=1)

	s=p.get_structure("scr", sys.argv[1])

	ppb=PPBuilder()

	print "C-N"
	for pp in ppb.build_peptides(s):
		print pp

	ppb=CaPPBuilder()

	print "CA-CA"
	for pp in ppb.build_peptides(s):
		print pp



