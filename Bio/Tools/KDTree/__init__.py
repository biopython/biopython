try:
	import Numeric
except ImportError:
	raise ImportError, "This module requires NumPy"

import _KDTree

class KDTree:
	def __init__(self, dim):
		self.dim=dim
		self.kdt=_KDTree.KDTree(dim)

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

		o center - two dimensional Numpy array of type "f". E.g. if the 
		points have dimensionality D, the center array should be 1xD 
		dimensional. 
		o radius - float>=0
		"""
		if len(center.shape)!=2 or center.shape[1]!=self.dim:
				raise Exception, "Expected a 1x%i Numpy array" % self.dim
		if center.typecode()!="f":
				raise Exception, "Expected a Numpy array of type float" 
		self.kdt.search_center_radius(center, radius)

	def neighbor_search(self, radius):
		"""Fixed neighbor search.

		Search all point pairs that are within radius.

		o radius - float (>0)
		"""
		self.kdt.neighbor_search(radius)

	def get_neighbors(self):
		"""Return Fixed neighbor search results.

		Returns a Nx2 dim Numpy array containing
		the indices of the point pairs, where N
		is the number of neighbor pairs.
		"""
		a=self.kdt.get_neighbor_indices()
		if a is None:
			return a
		# return as Nx2 dim Numpy array, where N
		# is number of neighbor pairs.
		a.shape=(-1, 2)
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

if __name__=="__main__":

	from RandomArray import *

	DIM=3
	N=9000

	i=1
	while(i):
		kdt=KDTree(3)	
		coords=random((N, DIM)).astype("f")
		#center=random((1, DIM)).astype("f")
		kdt.set_coords(coords)
		#kdt.search(center, 0.1)
		#indices=kdt.get_indices()
		#radii=kdt.get_radii()
		kdt.neighbor_search(0.01)
		indices=kdt.get_neighbors()
		if indices is None:
			l="no"
		else:
			l=len(indices)/2
		print "Cycle ", i,
		print " found ", l, " points."
		i=i+1



		
