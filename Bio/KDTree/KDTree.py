"""
The KD tree data structure can be used for all kinds of searches that
involve N-dimensional vectors, e.g.  neighbor searches (find all points
within a radius of a given point) or finding all point pairs in a set
that are within a certain radius of each other. See "Computational Geometry: 
Algorithms and Applications" (Mark de Berg, Marc van Kreveld, Mark Overmars, 
Otfried Schwarzkopf). Author: Thomas Hamelryck.
"""

try:
    import Numeric
    from Numeric import sum, sqrt
    from RandomArray import *
except ImportError:
    raise ImportError, "This module requires Numeric (precursor to NumPy)"

import CKDTree 

def _dist(p, q):
    diff=p-q
    return sqrt(sum(diff*diff))

def _neighbor_test(nr_points, dim, bucket_size, radius):
    """ Test all fixed radius neighbor search.

    Test all fixed radius neighbor search using the 
    KD tree C module.

    o nr_points - number of points used in test
    o dim - dimension of coords
    o bucket_size - nr of points per tree node
    o radius - radius of search (typically 0.05 or so) 
    """
    # KD tree search
    kdt=CKDTree.KDTree(dim, bucket_size)
    coords=random((nr_points, dim)).astype("f")
    kdt.set_data(coords, nr_points)
    kdt.neighbor_search(radius)
    r=kdt.neighbor_get_radii()
    if r is None:
        l1=0
    else:
        l1=len(r)
    # now do a slow search to compare results
    kdt.neighbor_simple_search(radius)
    r=kdt.neighbor_get_radii()
    if r is None:
        l2=0
    else:
        l2=len(r)
    if l1==l2:
        print "Passed."
    else:
        print "Not passed: %i <> %i." % (l1, l2)

def _test(nr_points, dim, bucket_size, radius):
    """Test neighbor search.

    Test neighbor search using the KD tree C module.

    o nr_points - number of points used in test
    o dim - dimension of coords
    o bucket_size - nr of points per tree node
    o radius - radius of search (typically 0.05 or so) 
    """
    # kd tree search
    kdt=CKDTree.KDTree(dim, bucket_size)
    coords=random((nr_points, dim)).astype("f")
    center=coords[0]
    kdt.set_data(coords, nr_points)
    kdt.search_center_radius(center, radius)
    r=kdt.get_indices()
    if r is None:
        l1=0
    else:
        l1=len(r)
    l2=0
    # now do a manual search to compare results
    for i in range(0, nr_points):
        p=coords[i]
        if _dist(p, center)<=radius:
            l2=l2+1
    if l1==l2:
        print "Passed."
    else:
        print "Not passed: %i <> %i." % (l1, l2)

class KDTree:
    """
    KD tree implementation (C++, SWIG python wrapper)

    The KD tree data structure can be used for all kinds of searches that
    involve N-dimensional vectors, e.g.  neighbor searches (find all points
    within a radius of a given point) or finding all point pairs in a set
    that are within a certain radius of each other.

    Reference:

    Computational Geometry: Algorithms and Applications
    Second Edition
    Mark de Berg, Marc van Kreveld, Mark Overmars, Otfried Schwarzkopf
    published by Springer-Verlag
    2nd rev. ed. 2000. 
    ISBN: 3-540-65620-0

    The KD tree data structure is described in chapter 5, pg. 99. 

    The following article made clear to me that the nodes should 
    contain more than one point (this leads to dramatic speed 
    improvements for the "all fixed radius neighbor search", see
    below):

    JL Bentley, "Kd trees for semidynamic point sets," in Sixth Annual ACM
    Symposium on Computational Geometry, vol. 91. San Francisco, 1990

    This KD implementation also performs a "all fixed radius neighbor search",
    i.e. it can find all point pairs in a set that are within a certain radius
    of each other. As far as I know the algorithm has not been published.
    """

    def __init__(self, dim, bucket_size=1):
        self.dim=dim
        self.kdt=CKDTree.KDTree(dim, bucket_size)
        self.built=0

    # Set data

    def set_coords(self, coords):
        """Add the coordinates of the points.

        o coords - two dimensional Numeric array of type "f". E.g. if the 
        points have dimensionality D and there are N points, the coords 
        array should be NxD dimensional. 
        """
        if min(coords)<=-1e6 or max(coords)>=1e6:
                raise Exception, "Points should lie between -1e6 and 1e6"
        if len(coords.shape)!=2 or coords.shape[1]!=self.dim:
                raise Exception, "Expected a Nx%i Numeric array" % self.dim
        if coords.typecode()!="f":
                raise Exception, "Expected a Numeric array of type float" 
        self.kdt.set_data(coords, coords.shape[0])
        self.built=1

    # Fixed radius search for a point

    def search(self, center, radius):
        """Search all points within radius of center.

        o center - one dimensional Numeric array of type "f". E.g. if the 
        points have dimensionality D, the center array should be D 
        dimensional. 
        o radius - float>0
        """
        if not self.built:
                raise Exception, "No point set specified"
        if center.shape!=(self.dim,):
                raise Exception, "Expected a %i-dimensional Numeric array" % self.dim
        if center.typecode()!="f":
                raise Exception, "Expected a Numeric array of type float" 
        self.kdt.search_center_radius(center, radius)

    def get_radii(self):
        """Return radii.

        Return the list of distances from center after
        a neighbor search.
        """
        a=self.kdt.get_radii()
        if a is None:
            return []
        return a
    
    def get_indices(self):
        """Return the list of indices.

        Return the list of indices after a neighbor search.
        The indices refer to the original coords Numeric array. The
        coordinates with these indices were within radius of center.

        For an index pair, the first index<second index. 
        """
        a=self.kdt.get_indices()
        if a is None:
            return []
        return a

    # Fixed radius search for all points


    def all_search(self, radius):
        """All fixed neighbor search.

        Search all point pairs that are within radius.

        o radius - float (>0)
        """
        if not self.built:
                raise Exception, "No point set specified"
        self.kdt.neighbor_search(radius)

    def all_get_indices(self):
        """Return All Fixed Neighbor Search results.

        Return a Nx2 dim Numeric array containing
        the indices of the point pairs, where N
        is the number of neighbor pairs.
        """
        a=self.kdt.neighbor_get_indices()
        if a is None:
            return [] 
        # return as Nx2 dim Numeric array, where N
        # is number of neighbor pairs.
        a.shape=(-1, 2)
        return a

    def all_get_radii(self):
        """Return All Fixed Neighbor Search results.

        Return an N-dim array containing the distances
        of all the point pairs, where N is the number 
        of neighbor pairs..
        """
        a=self.kdt.neighbor_get_radii()
        if a is None:
            return [] 
        return a

if __name__=="__main__":

    from RandomArray import *

    nr_points=100000
    dim=3
    bucket_size=10
    query_radius=10

    coords=(200*random((nr_points, dim))).astype("f")

    kdtree=KDTree(dim, bucket_size)

    # enter coords
    kdtree.set_coords(coords)

    # Find all point pairs within radius

    kdtree.all_search(query_radius)

    # get indices & radii of points

    # indices is a list of tuples. Each tuple contains the 
    # two indices of a point pair within query_radius of 
    # each other.
    indices=kdtree.all_get_indices() 
    radii=kdtree.all_get_radii()

    print "Found %i point pairs within radius %f." % (len(indices), query_radius)

    # Do 10 individual queries

    for i in range(0, 10):
        # pick a random center
        center=random(dim).astype("f")
        
        # search neighbors
        kdtree.search(center, query_radius)

        # get indices & radii of points
        indices=kdtree.get_indices()
        radii=kdtree.get_radii()

        x, y, z=center
        print "Found %i points in radius %f around center (%.2f, %.2f, %.2f)." % (len(indices), query_radius, x, y, z)

