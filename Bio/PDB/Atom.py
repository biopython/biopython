# Copyright (C) 2002, Thomas Hamelryck (thamelry@vub.ac.be)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.           

# My stuff
from Entity import DisorderedEntityWrapper


class Atom:
	def __init__(self, name, coord, bfactor, occupancy, altloc, fullname):
		"""Atom object.

		The Atom object stores atom name (both with and without spaces), 
		coordinates, B factor, occupancy, alternative location specifier
		and (optionally) anisotropic B factor and standard deviations of 
		B factor and positions.
  
		Arguments:
		o name - string, atom name e.g. "CA", note that spaces are normally stripped
		o coord - Numpy array (Float0, size 3), atomic coordinates
		o bfactor - float
		o occupancy - float
		o altloc - string, alternative location specifier for disordered atoms
		o fullname - string, full atom name, including spaces, e.g. " CA ". Normally
		these spaces are stripped from the atom name. 
		"""
		self.level="A"
		# Reference to the residue 
		self.parent=None
		# the atomic data
		self.name=name		# eg. CA, spaces are removed from atom name
		self.fullname=fullname	# e.g. " CA ", spaces included
		self.coord=coord
		self.bfactor=bfactor
		self.occupancy=occupancy
		self.altloc=altloc
		self.full_id=None	# (structure id, model id, chain id, residue id, atom id)
		self.id=name		# id of atom is the atom name (e.g. "CA")
		self.disordered_flag=0
		self.anisou_array=None
		self.siguij_array=None
		self.sigatm_array=None

	# Special methods	

	def __repr__(self):
		return "<Atom %s>" % self.get_id()

	# set methods

	def set_bfactor(self, bfactor):
		self.bfactor=bfactor

	def set_coord(self, coord):
		self.coord=coord

	def set_altloc(self, altloc):
		self.altloc=altloc

	def set_occupancy(self, occupancy):
		self.occupancy=occupancy

	def set_sigatm(self, sigatm_array):
		"""Set standard deviation of atomic parameters.

		The standard deviation of atomic parameters consists
		of 3 positional, 1 B factor and 1 occupancy standard 
		deviation.

		Arguments:
		o sigatm_array - Numpy array (length 5), standard deviations of atomic 
		parameters.
		"""
		self.sigatm_array=sigatm_array

	def set_siguij(self, siguij_array):
		"""Set standard deviations of anisotropic temperature factors.

		Arguments:
		o siguij_array - Numpy array (length 6), standard deviations of 
		anisotropic temperature factors.
		"""
		self.siguij_array=siguij_array

	def set_anisou(self, anisou_array):
		"""Set anisotropic B factor.

		Arguments:
		o anisou_array -  Numpy array (length 6), anisotropic B factor.
		"""
		self.anisou_array=anisou_array


	# Public methods	

	def flag_disorder(self):
		"""Set the disordered flag to 1.

		The disordered flag indicates whether the atom is disordered or not.
		"""
		self.disordered_flag=1

	def is_disordered(self):
		"Return the disordered flag (1 if disordered, 0 otherwise)."
		return self.disordered_flag 

	def set_parent(self, parent):
		"""Set the parent residue.

		Arguments:
		o parent - Residue object
		"""
		self.parent=parent
	
	def detach_parent(self):
		"Remove reference to parent."
		self.parent=None

	def get_sigatm(self):
		"Return standard deviation of atomic parameters."
		return self.sigatm_array

	def get_siguij(self):
		"Return standard deviations of anisotropic temperature factors."
		return self.siguij_array

	def get_anisou(self):
		"Return anisotropic B factor."
		return self.anisou_array

	def get_parent(self):
		"Return parent residue."
		return self.parent

	def destroy(self):
		self.parent=None
		del self.id
		del self.name

	def get_name(self):
		"Return atom name."
		return self.name

	def get_id(self):
		"Return the id of the atom (which is its atom name)."
		return self.id

	def get_full_id(self):
		"""Return the full id of the atom.

		The full id of an atom is the tuple 
		(structure id, model id, chain id, residue id, atom name, altloc).
		"""
		return self.parent.get_full_id()+((self.name, self.altloc),)
	
	def get_coord(self):
		"Return atomic coordinates."
		return self.coord

	def get_bfactor(self):
		"Return B factor."
		return self.bfactor

	def get_occupancy(self):
		"Return occupancy."
		return self.occupancy

	def get_fullname(self):
		"Return the atom name, including leading and trailing spaces."
		return self.fullname

	def get_altloc(self):
		"Return alternative location specifier."
		return self.altloc

	def get_level(self):
		return self.level


class DisorderedAtom(DisorderedEntityWrapper):
	"""
	This class contains all Atom objects that represent the same disordered
	atom. One of these atoms is "selected" and all method calls not caught
	by DisorderedAtom are forwarded to the selected Atom object. In that way, a
	DisorderedAtom behaves exactly like a normal Atom. By default, the selected 
	Atom object represents the Atom object with the highest occupancy, but a 
	different Atom object can be selected by using the disordered_select(altloc) 
	method. 
	"""
	def __init__(self, id):
		"""
		Arguments:
		o id - string, atom name
		"""
		self.last_occupancy=-1
		DisorderedEntityWrapper.__init__(self, id)

	# Special methods

	def __repr__(self):
		return "<Disordered Atom %s>" % self.get_id() 

	def disordered_add(self, atom):
		"Add a disordered atom."
		# Add atom to dict, use altloc as key	
		atom.flag_disorder()
		altloc=atom.get_altloc()
		occupancy=atom.get_occupancy()
		self[altloc]=atom
		if occupancy>self.last_occupancy:
			self.last_occupancy=occupancy
			self.disordered_select(altloc)
