from Bio.KDTree import _neighbor_test, _test

nr_points=5000
dim=3
bucket_size=5
radius=0.01

for i in range(0, 10):
	_neighbor_test(nr_points, dim, bucket_size, radius)
	_test(nr_points, dim, bucket_size, radius)
