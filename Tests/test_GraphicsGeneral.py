# Copyright 2001 by Brad Chapman.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Test for graphics things that don't really deserve there own test module."""


import os
import random
import unittest

from Bio import MissingExternalDependencyError

try:
    # Skip the test if reportlab is not installed
    import reportlab as r

    del r
except ImportError:
    raise MissingExternalDependencyError(
        "Install reportlab if you want to use Bio.Graphics."
    ) from None

# the stuff we're testing
from Bio.Graphics.Comparative import ComparativeScatterPlot


class ComparativeTest(unittest.TestCase):
    """Do tests for modules involved with comparing data."""

    def setUp(self):
        self.min_num_points = 1
        self.max_num_points = 500

        self.min_point_num = 0
        self.max_point_num = 200

    def _make_random_points(self, num_two_d_lists):
        """Make a bunch of random points for testing plots."""
        plot_info = []
        random.seed(num_two_d_lists)  # for reproducibility
        for two_d_list in range(num_two_d_lists):
            cur_list = []
            num_points = random.randrange(self.min_num_points, self.max_num_points)
            for point in range(num_points):
                x_point = random.randrange(self.min_point_num, self.max_point_num)
                y_point = random.randrange(self.min_point_num, self.max_point_num)

                cur_list.append((x_point, y_point))

            plot_info.append(cur_list)
        return plot_info

    def test_simple_scatter_plot_1(self):
        """Test creation of a simple ScatterPlot with one list."""
        compare_plot = ComparativeScatterPlot()
        compare_plot.display_info = self._make_random_points(1)
        output_file = os.path.join(os.getcwd(), "Graphics", "scatter_test_1.pdf")
        compare_plot.draw_to_file(output_file, "Testing Scatter Plots")

    def test_simple_scatter_plot_7(self):
        """Test creation of a simple ScatterPlot with more lists."""
        compare_plot = ComparativeScatterPlot()
        # There are 6 pre-defined colors and symbols, doing more:
        compare_plot.display_info = self._make_random_points(7)
        output_file = os.path.join(os.getcwd(), "Graphics", "scatter_test_7.pdf")
        compare_plot.draw_to_file(output_file, "Testing Scatter Plots")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
