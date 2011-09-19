# Copyright 2004-2008 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# See the Biopython Tutorial for an explanation of the biological
# background of these tests.

import unittest

try:
    import numpy
    from numpy import asarray
    del asarray
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.kNN.")

from Bio import kNN

xs = [[-53, -200.78],
      [117, -267.14],
      [57, -163.47],
      [16, -190.30],
      [11, -220.94],
      [85, -193.94],
      [16, -182.71],
      [15, -180.41],
      [-26, -181.73],
      [58, -259.87],
      [126, -414.53],
      [191, -249.57],
      [113, -265.28],
      [145, -312.99],
      [154, -213.83],
      [147, -380.85],
      [93, -291.13]]

ys = [1,
      1,
      1,
      1,
      1,
      1,
      1,
      1,
      1,
      1,
      0,
      0,
      0,
      0,
      0,
      0,
      0]

class TestKNN(unittest.TestCase):

    def test_calculate_model(self):
        k = 3
        model = kNN.train(xs, ys, k)
        self.assertEqual(model.classes, set([0,1]))
        n = len(xs)
        for i in range(n):
            self.assertAlmostEqual(model.xs[i,0], xs[i][0], places=4)
            self.assertAlmostEqual(model.xs[i,1], xs[i][1], places=4)
            self.assertEqual(model.ys[i], ys[i])
        self.assertEqual(model.k, k)

    def test_classify(self):
        k = 3
        model = kNN.train(xs, ys, k)
        result = kNN.classify(model, [6,-173.143442352])
        self.assertEqual(result, 1)
        result = kNN.classify(model, [309, -271.005880394])
        self.assertEqual(result, 0)

    def test_calculate_probability(self):
        k = 3
        model = kNN.train(xs, ys, k)
        weights = kNN.calculate(model, [6,-173.143442352])
        self.assertAlmostEqual(weights[0], 0.0, places=6)
        self.assertAlmostEqual(weights[1], 3.0, places=6)
        weights = kNN.calculate(model, [309, -271.005880394])
        self.assertAlmostEqual(weights[0], 3.0, places=6)
        self.assertAlmostEqual(weights[1], 0.0, places=6)
        weights = kNN.calculate(model, [117, -267.13999999999999])
        self.assertAlmostEqual(weights[0], 2.0, places=6)
        self.assertAlmostEqual(weights[1], 1.0, places=6)

    def test_model_accuracy(self):
        correct = 0
        k = 3
        model = kNN.train(xs, ys, k)
        predictions = [1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0]
        for i in range(len(predictions)):
            prediction = kNN.classify(model, xs[i])
            self.assertEqual(prediction, predictions[i])
            if prediction==ys[i]:
                correct+=1
        self.assertEqual(correct, 15)

    def test_leave_one_out(self):
        correct = 0
        k = 3
        model = kNN.train(xs, ys, k)
        predictions = [1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1]
        for i in range(len(predictions)):
            model = kNN.train(xs[:i]+xs[i+1:], ys[:i]+ys[i+1:], k)
            prediction = kNN.classify(model, xs[i])
            self.assertEqual(prediction, predictions[i])
            if prediction==ys[i]:
                correct+=1
        self.assertEqual(correct, 13)

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
