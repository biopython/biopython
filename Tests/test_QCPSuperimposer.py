# Copyright 2017 by Maximilian Greil.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import unittest

try:
    from numpy import array
    from numpy import dot  # missing in old PyPy's micronumpy
    from numpy import around
    from numpy import array_equal
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


class QCPSuperimposerTest(unittest.TestCase):

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
        self.sup.set(self.x, self.y)

    # Public methods

    def test_set(self):
        self.assertTrue(
            array_equal(around(self.sup.reference_coords, decimals=3), around(self.x, decimals=3)))
        self.assertTrue(
            array_equal(around(self.sup.coords, decimals=3), around(self.y, decimals=3)))
        self.assertIsNone(self.sup.transformed_coords)
        self.assertIsNone(self.sup.rot)
        self.assertIsNone(self.sup.tran)
        self.assertIsNone(self.sup.rms)
        self.assertIsNone(self.sup.init_rms)

    def test_run(self):
        self.sup.run()
        self.assertTrue(
            array_equal(around(self.sup.reference_coords, decimals=3), around(self.x, decimals=3)))
        self.assertTrue(
            array_equal(around(self.sup.coords, decimals=3), around(self.y, decimals=3)))
        self.assertIsNone(self.sup.transformed_coords)
        calc_rot = array([[0.68304939, -0.5227742, -0.51004967],
                          [0.53664482, 0.83293151, -0.13504605],
                          [0.49543503, -0.18147239, 0.84947743]])
        self.assertTrue(
            array_equal(around(self.sup.rot, decimals=3), around(calc_rot, decimals=3)))
        calc_tran = array([-7.43885125, 36.51522275, 36.81135533])
        self.assertTrue(
            array_equal(around(self.sup.tran, decimals=3), around(calc_tran, decimals=3)))
        calc_rms = 0.003
        self.assertEqual(float('%.3f' % self.sup.rms), calc_rms)
        self.assertIsNone(self.sup.init_rms)

    def test_get_transformed(self):
        self.sup.run()
        transformed_coords = array([[49.05456082, -1.23928381, 50.58427566],
                                    [50.02204971, -0.39367917, 51.42494181],
                                    [51.47738545, -0.57287128, 51.06760944],
                                    [52.39602307, -0.98417088, 52.03318837]])
        self.assertTrue(
            array_equal(around(self.sup.get_transformed(), decimals=3), around(transformed_coords, decimals=3)))

    def test_get_init_rms(self):
        x = array([[1.1, 1.2, 1.3],
                  [1.4, 1.5, 1.6],
                  [1.7, 1.8, 1.9]])
        y = array([[1.0, 1.0, 1.0],
                  [1.0, 1.0, 1.0],
                  [1.0, 1.0, 1.0]])
        self.sup.set(x, y)
        self.assertIsNone(self.sup.init_rms)
        init_rms = array([0.81, 0.9, 0.98])
        self.assertTrue(
            array_equal(around(self.sup.get_init_rms(), decimals=2), around(init_rms, decimals=2)))

    def test_get_rotran(self):
        self.sup.run()
        calc_rot = array([[0.68304939, -0.5227742, -0.51004967],
                          [0.53664482, 0.83293151, -0.13504605],
                          [0.49543503, -0.18147239, 0.84947743]])
        calc_tran = array([-7.43885125, 36.51522275, 36.81135533])
        rot, tran = self.sup.get_rotran()
        self.assertTrue(
            array_equal(around(rot, decimals=3), around(calc_rot, decimals=3)))
        self.assertTrue(
            array_equal(around(tran, decimals=3), around(calc_tran, decimals=3)))

    def test_get_rms(self):
        self.sup.run()
        calc_rms = 0.003
        self.assertEqual(float('%.3f' % self.sup.get_rms()), calc_rms)

    # Old test from Bio/PDB/QCPSuperimposer/__init__.py

    def test_oldTest(self):
        x = array([[-2.803, -15.373, 24.556],
                   [0.893, -16.062, 25.147],
                   [1.368, -12.371, 25.885],
                   [-1.651, -12.153, 28.177],
                   [-0.440, -15.218, 30.068],
                   [2.551, -13.273, 31.372],
                   [0.105, -11.330, 33.567]])

        y = array([[-14.739, -18.673, 15.040],
                   [-12.473, -15.810, 16.074],
                   [-14.802, -13.307, 14.408],
                   [-17.782, -14.852, 16.171],
                   [-16.124, -14.617, 19.584],
                   [-15.029, -11.037, 18.902],
                   [-18.577, -10.001, 17.996]])

        self.sup.set(x, y)
        self.assertTrue(
            array_equal(around(self.sup.reference_coords, decimals=3), around(x, decimals=3)))
        self.assertTrue(
            array_equal(around(self.sup.coords, decimals=3), around(y, decimals=3)))

        self.sup.run()
        rot = array([[0.72216358, 0.69118937, -0.0271479],
                     [-0.52038257, 0.51700833, -0.67963547],
                     [-0.45572112, 0.50493528, 0.73304748]])
        tran = array([11.68878393, -4.13245037, 6.05208344])
        rms = 0.7191064509622271
        self.assertTrue(
            array_equal(around(self.sup.rot, decimals=3), around(rot, decimals=3)))
        self.assertTrue(
            array_equal(around(self.sup.tran, decimals=3), around(tran, decimals=3)))
        self.assertEqual(float('%.3f' % self.sup.rms), around(rms, decimals=3))

        rms_get = self.sup.get_rms()
        self.assertTrue(
            array_equal(around(rms_get, decimals=3), around(rms, decimals=3)))

        rot_get, tran_get = self.sup.get_rotran()
        self.assertTrue(
            array_equal(around(rot_get, decimals=3), around(rot, decimals=3)))
        self.assertTrue(
            array_equal(around(tran_get, decimals=3), around(tran, decimals=3)))

        y_on_x1 = dot(y, rot) + tran
        y_x_solution = array([[3.90787281, -16.37976032, 30.16808376],
                              [3.58322456, -12.81122728, 28.91874135],
                              [1.3580194, -13.96815768, 26.05958411],
                              [-0.79347336, -15.93647897, 28.48288439],
                              [-1.27379224, -12.9456459, 30.78004989],
                              [-2.03519089, -10.6822696, 27.81728956],
                              [-4.72366028, -13.05646024, 26.54536695]])
        self.assertTrue(
            array_equal(around(y_on_x1, decimals=3), around(y_x_solution, decimals=3)))

        y_on_x2 = self.sup.get_transformed()
        self.assertTrue(
            array_equal(around(y_on_x2, decimals=3), around(y_x_solution, decimals=3)))

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
