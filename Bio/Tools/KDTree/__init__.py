try:
	import Numeric
	from RandomArray import *
except ImportError:
	raise ImportError, "This module requires NumPy"

from Bio.Tools.KDTree import _KDTree 

def _neighbor_test(nr_points, dim, bucket_size, radius):
		kdt=_KDTree.KDTree(dim, bucket_size)
		coords=random((nr_points, dim)).astype("f")
		kdt.set_data(coords, nr_points)
		kdt.neighbor_search(radius)
		l1=len(kdt.neighbor_get_radii())
		kdt.neighbor_simple_search(radius)
		l2=len(kdt.neighbor_get_radii())
		if l1==l2:
			print "Passed."
		else:
			print "Not passed: %i <> %i." % (l1, l2)

def _test(nr_points, dim, bucket_size, radius):
	radius_sq=radius*radius
	kdt=_KDTree.KDTree(dim, bucket_size)
	coords=random((nr_points, dim)).astype("f")
	center=coords[0]
	kdt.set_data(coords, nr_points)
	kdt.search_center_radius(center, radius)
	l1=len(kdt.get_indices())
	l2=0
	for i in range(0, nr_points):
		p=coords[i]
		if _KDTree.KDTREE_dist(p, center, dim)<=radius_sq:
			l2=l2+1
	if l1==l2:
		print "Passed."
	else:
		print "Not passed: %i <> %i." % (l1, l2)

class KDTree:
	def __init__(self, dim, bucket_size=1):
		self.dim=dim
		self.kdt=_KDTree.KDTree(dim, bucket_size)
		self.built=0

	# Set data

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
		self.built=1

	# Fixed radius search for a point

	def search(self, center, radius):
		"""Search all points within radius of center.

		o center - one dimensional Numpy array of type "f". E.g. if the 
		points have dimensionality D, the center array should be D 
		dimensional. 
		o radius - float>=0
		"""
		if not self.built:
				raise Exception, "No points set specified"
		if center.shape!=(self.dim,):
				raise Exception, "Expected a %i-dimensional Numpy array" % self.dim
		if center.typecode()!="f":
				raise Exception, "Expected a Numpy array of type float" 
		self.kdt.search_center_radius(center, radius)

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

	# Fixed radius search for all points


	def neighbor_search(self, radius):
		"""Fixed neighbor search.

		Search all point pairs that are within radius.

		o radius - float (>0)
		"""
		if not self.built:
				raise Exception, "No points set specified"
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

if __name__=="__main__":

	nr_points=1000
	dim=3
	bucket_size=10
	radius=0.05

	while(1):
 		 _neighbor_test(nr_points, dim, bucket_size, radius)
		 _test(nr_points, dim, bucket_size, radius)

