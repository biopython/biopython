#TODO - What if numpy was installed after Biopython,
#meaning KDTree's C code hasn't been compiled?
try :
    import numpy
except ImportError :
    from Bio import MissingExternalDependencyError
    raise MissingExternalDependencyError(\
        "Install NumPy if you want to use Bio.KDTree.")

from Bio.KDTree.KDTree import _neighbor_test, _test

nr_points=5000
dim=3
bucket_size=5
radius=0.01

for i in range(0, 10):
	_neighbor_test(nr_points, dim, bucket_size, radius)
	_test(nr_points, dim, bucket_size, radius)
