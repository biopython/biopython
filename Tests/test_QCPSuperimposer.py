# Copyright 2017 by Maximilian Greil.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import unittest

try:
    from numpy import array
    from numpy import dot  # missing in old PyPy's micronumpy
    from numpy import array_equal
    from numpy import around 
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.QCPSuperimposer.")

try:
    from Bio.PDB.QCPSuperimposer import QCPSuperimposer
except ImportError:
    from Bio import MissingExternalDependencyError
    raise MissingExternalDependencyError(
        "C module in Bio.QCPSuperimposer not compiled")

class QCSuperimposerTest(unittest.TestCase):

	def setUp(self):
		self.x = array([[51.65, -1.90, 50.07],
		              [50.40, -1.23, 50.65],
		              [50.68, -0.04, 51.54],
		              [50.22, -0.02, 52.85]])

		self.y = array([[51.30, -2.99, 46.54],
		              [51.09, -1.88, 47.58],
		              [52.36, -1.20, 48.03],
		              [52.71, -1.18, 49.38]])

		self.sup = QCPSuperimposer()
		self.sup.set(self.x,self.y)

	#Public methods

	def testSet(self):
		"""Test method set(x,y)"""
		assert(array_equal(self.sup.reference_coords, self.x))
		assert(array_equal(self.sup.coords, self.y))
		assert(self.sup.transformed_coords == None)
		assert(self.sup.rot == None)
		assert(self.sup.tran == None)
		assert(self.sup.rms == None)
		assert(self.sup.init_rms == None)

	def testRun(self):
		"""Test method run()"""
		self.sup.run()
		assert(array_equal(self.sup.reference_coords, self.x))
		assert(array_equal(self.sup.coords, self.y))
		assert(self.sup.transformed_coords == None)
		calc_rot = array([[0.68304939, -0.5227742, -0.51004967],
 					    [0.53664482, 0.83293151, -0.13504605],
 					    [0.49543503, -0.18147239, 0.84947743]])
		assert(array_equal(around(self.sup.rot, decimals=3), around(calc_rot, decimals=3)))
		calc_tran = array([-7.43885125, 36.51522275, 36.81135533])
		assert(array_equal(around(self.sup.tran, decimals=3), around(calc_tran, decimals=3)))
		calc_rms = 0.003
		assert(float(('%.3f' % self.sup.rms)) == calc_rms)
		assert(self.sup.init_rms == None)

	def testGet_transformed(self):
		"""Test method get_transformed()"""
		self.sup.run()
		transformed_coords = array([[49.05456082, -1.23928381, 50.58427566],
 					 [ 50.02204971, -0.39367917, 51.42494181],
 			         [ 51.47738545, -0.57287128, 51.06760944],
 					 [ 52.39602307, -0.98417088, 52.03318837]])
		assert(array_equal(around(self.sup.get_transformed(), decimals=3), around(transformed_coords, decimals=3)))

	def testGet_init_rms(self):
		"""Test method get_init_rms()"""
		x = array([[1.1, 1.2, 1.3],
				  [1.4, 1.5, 1.6],
				  [1.7, 1.8, 1.9]])
		y = array([[1.0, 1.0, 1.0],
				  [1.0, 1.0, 1.0],
				  [1.0, 1.0, 1.0]])
		self.sup.set(x,y)
		assert(self.sup.init_rms == None)
		init_rms = array([0.81, 0.9, 0.98])
		assert(array_equal(around(self.sup.get_init_rms(), decimals=2), around(init_rms, decimals=2)))

	def testGet_rotran(self):
		"""Test method get_rotran()"""
		self.sup.run()
		calc_rot = array([[0.68304939, -0.5227742, -0.51004967],
 					    [0.53664482, 0.83293151, -0.13504605],
 					    [0.49543503, -0.18147239, 0.84947743]])
		calc_tran = array([-7.43885125, 36.51522275, 36.81135533])
		rot, tran = self.sup.get_rotran()
		assert(array_equal(around(rot, decimals=3), around(calc_rot, decimals=3)))
		assert(array_equal(around(tran, decimals=3), around(calc_tran, decimals=3)))
		
	def testGet_rms(self):
		"""Test method get_rms()"""
		self.sup.run()
		calc_rms = 0.003
		assert(float(('%.3f' % self.sup.get_rms())) == calc_rms)

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
