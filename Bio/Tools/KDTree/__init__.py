try:
	import Numeric
except ImportError:
	raise ImportError, "This module requires NumPy"

import _KDTree

class KDTree:
	def __init__(self, dim, bucket_size=1):
		self.dim=dim
		self.kdt=_KDTree.KDTree(dim, bucket_size)
		self.built=0

	def set_coords(self, coords):
		"""Add the coordinates of the points.

		o coords - two dimensional Numpy array of type "f". E.g. if the 
		points have dimensionality D and there are N points, the coords 
		array should be NxD dimensional. 
		"""
		if len(coords.shape)!=2 or coords.shape[1]!=self.dim:
				raise Exception, "Expected a Nx%i Numpy array" % self.dim
		if coords.typecode()!="f":
				raise Exception, "Expected a Numpy array of type float" 
		self.kdt.set_data(coords, coords.shape[0])

	def search(self, center, radius):
		"""Search all points within radius of center.

		o center - one dimensional Numpy array of type "f". E.g. if the 
		points have dimensionality D, the center array should be D 
		dimensional. 
		o radius - float>=0
		"""
		if center.shape!=(self.dim,):
				raise Exception, "Expected a %i-dimensional Numpy array" % self.dim
		if center.typecode()!="f":
				raise Exception, "Expected a Numpy array of type float" 
		self.kdt.search_center_radius(center, radius)

	def neighbor_search(self, radius):
		"""Fixed neighbor search.

		Search all point pairs that are within radius.

		o radius - float (>0)
		"""
		self.kdt.neighbor_search(radius)

	def neighbor_get_indices(self):
		"""Return Fixed neighbor search results.

		Returns a Nx2 dim Numpy array containing
		the indices of the point pairs, where N
		is the number of neighbor pairs.
		"""
		a=self.kdt.neighbor_get_indices()
		if a is None:
			return a
		# return as Nx2 dim Numpy array, where N
		# is number of neighbor pairs.
		a.shape=(-1, 2)
		return a

	def neighbor_get_radii(self):
		a=self.kdt.neighbor_get_radii()
		return a

	def get_radii(self):
		"""Return radii.

		Return the list of distances from center.
		"""
		return self.kdt.get_radii()
	
	def get_indices(self):
		"""Return the list of indices.

		The indices refer to the original coords Numpy array. The
		coordinates with these indices were within radius of center.
		"""
		return self.kdt.get_indices()

	def neighbor_simple_search(self, radius):
		if self.built==1:
			raise Exception, "KD tree is already built."
		self.kdt.neighbor_simple_search(radius)

if __name__=="__main__":

	from RandomArray import *
	from Bio.Tools.KDTree import *

	DIM=3
	N=9000

	kdt=KDTree(3)	

	coords=random((N, DIM)).astype("f")
	
	kdt.set_coords(coords)

	R=0.01

	kdt.neighbor_simple_search(R)
	print "Simple search : ", len(kdt.neighbor_get_radii())
	
	kdt.neighbor_search(R)
	print "KD tree search : ", len(kdt.neighbors_get_radii())



		
