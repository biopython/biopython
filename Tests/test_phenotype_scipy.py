# Copyright 2014-2016 Marco Galardini.  All rights reserved.
# Adapted from test_Mymodule.py by Jeff Chang
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

try:
    import numpy
except ImportError:
    from Bio import MissingExternalDependencyError
    raise MissingExternalDependencyError(
        "Install NumPy if you want to use Bio.phenotype.")
try:
    import scipy
except ImportError:
    from Bio import MissingExternalDependencyError
    raise MissingExternalDependencyError(
        "Install SciPy if you want to use fitting in Bio.phenotype.")

import json
import unittest

from Bio import BiopythonExperimentalWarning

import warnings
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)
    from Bio import phenotype

# Example plate files
JSON_PLATE = 'phenotype/Plate.json'


class TestPhenoMicro(unittest.TestCase):

    def test_WellRecord(self):
        """Test basic functionalities of WellRecord objects."""
        with open(JSON_PLATE) as handle:
            p = json.load(handle)

        times = p['measurements']['Hour']
        w = phenotype.phen_micro.WellRecord('A10',
                                            signals=dict([(times[i], p['measurements']['A10'][i])
                                                          for i in range(len(times))]))

        w1 = phenotype.phen_micro.WellRecord('H12',
                                             signals=dict([(times[i], p['measurements']['H12'][i])
                                                           for i in range(len(times))]))

        # self.assertIsInstance(w.plate,
        #                       phenotype.phen_micro.PlateRecord)
        self.assertTrue(isinstance(w.plate, phenotype.phen_micro.PlateRecord))
        self.assertEqual(w.id, 'A10')
        self.assertEqual(len(w), len(times))
        self.assertEqual(len(w), 384)
        self.assertEqual(max(w), (95.75, 217.0))
        self.assertEqual(min(w), (0.0, 37.0))
        self.assertEqual(max(w, key=lambda x: x[1]),
                         (16.75, 313.0))
        self.assertEqual(min(w, key=lambda x: x[1]),
                         (0.25, 29.0))
        self.assertEqual(len(w[:]), 96)
        self.assertEqual(w[1], 29.)
        self.assertEqual(w[12], 272.)
        self.assertEqual(w[1:5], [29., 35., 39., 43.])
        self.assertRaises(ValueError, w.__getitem__, 'a')
        self.assertAlmostEqual(w[1:2:.25][0], 29.)
        self.assertAlmostEqual(w[1.3567], 33.7196)
        self.assertEqual(w.get_raw()[0], (0.0, 37.0))
        self.assertEqual(w.get_raw()[-1], (95.75, 217.0))
        self.assertEqual(w.get_times()[0], 0.0)
        self.assertEqual(w.get_times()[-1], 95.75)
        self.assertEqual(w.get_signals()[0], 37.0)
        self.assertEqual(w.get_signals()[-1], 217.0)
        self.assertEqual(repr(w),
                         "WellRecord('(0.0, 37.0), (0.25, 29.0), (0.5, 32.0)," +
                         " (0.75, 30.0), (1.0, 29.0), ..., (95.75, 217.0)')")
        self.assertEqual(str(w),
                         "Well ID: A10\nTime points: 384\nMinum signal 0.25 at " +
                         "time 29.00\nMaximum signal 16.75 at time " +
                         "313.00\nWellRecord('(0.0, 37.0), (0.25, 29.0), " +
                         "(0.5, 32.0), (0.75, 30.0), " +
                         "(1.0, 29.0), ..., (95.75, 217.0)')")

        w.fit()
        self.assertAlmostEqual(w.area, 20879.5)
        self.assertEqual(w.model, 'gompertz')
        self.assertAlmostEqual(w.lag, 6.0425868725090357, places=5)
        self.assertAlmostEqual(w.plateau, 188.51404344898586, places=5)
        self.assertAlmostEqual(w.slope, 48.190618284831132, places=5)
        self.assertAlmostEqual(w.v, 0.10000000000000001, places=5)
        self.assertAlmostEqual(w.y0, 45.879770069807989, places=5)
        self.assertEqual(w.max, 313.0)
        self.assertEqual(w.min, 29.0)
        self.assertEqual(w.average_height, 217.82552083333334)

        self.assertRaises(TypeError, w.__add__, 'a')

        w2 = w + w1
        self.assertEqual(w2.id, 'A10')
        self.assertEqual(len(w2), len(times))
        self.assertEqual(len(w2), 384)
        self.assertEqual(max(w2), (95.75, 327.0))
        self.assertEqual(min(w2), (0.0, 63.0))
        self.assertEqual(max(w2, key=lambda x: x[1]),
                         (18.25, 357.0))
        self.assertEqual(min(w2, key=lambda x: x[1]),
                         (0.25, 55.0))
        self.assertEqual(w2[1], 71.)
        self.assertEqual(w2[12], 316.)
        self.assertEqual(w2[1:5], [71.0, 88.0, 94.0, 94.0])
        self.assertAlmostEqual(w2[1:2:.25][0], 71.0)
        self.assertAlmostEqual(w2[1.3567], 77.7196)
        self.assertEqual(w2.get_raw()[0], (0.0, 63.0))
        self.assertEqual(w2.get_raw()[-1], (95.75, 327.0))
        self.assertEqual(w2.get_times()[0], 0.0)
        self.assertEqual(w2.get_times()[-1], 95.75)
        self.assertEqual(w2.get_signals()[0], 63.0)
        self.assertEqual(w2.get_signals()[-1], 327.0)

        self.assertRaises(TypeError, w.__sub__, 'a')

        w2 = w - w1
        self.assertEqual(w2.id, 'A10')
        self.assertEqual(len(w2), len(times))
        self.assertEqual(len(w2), 384)
        self.assertEqual(max(w2), (95.75, 107.0))
        self.assertEqual(min(w2), (0.0, 11.0))
        self.assertEqual(max(w2, key=lambda x: x[1]),
                         (15.75, 274.0))
        self.assertEqual(min(w2, key=lambda x: x[1]),
                         (3.25, -20.0))
        self.assertEqual(w2[1], -13.)
        self.assertEqual(w2[12], 228.)
        self.assertEqual(w2[1:5], [-13.0, -18.0, -16.0, -8.0])
        self.assertAlmostEqual(w2[1:2:.25][0], -13.0)
        self.assertAlmostEqual(w2[1.3567], -10.2804)
        self.assertEqual(w2.get_raw()[0], (0.0, 11.0))
        self.assertEqual(w2.get_raw()[-1], (95.75, 107.0))
        self.assertEqual(w2.get_times()[0], 0.0)
        self.assertEqual(w2.get_times()[-1], 95.75)
        self.assertEqual(w2.get_signals()[0], 11.0)
        self.assertEqual(w2.get_signals()[-1], 107.0)

        w[1] = 1

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
