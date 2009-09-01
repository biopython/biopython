# Copyright 2001 by Brad Chapman.  All rights reserved.
# Copyright 2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Test for bitmap output via ReportLab (requires an extra dependency).

The primary purpose of this is to flag if renderPM is missing, but also
to check Bio.Graphics can make a bitmap (e.g. PNG).

The example itself is essentially a repeat from test_GraphicsGeneral.py.
"""

import os
import random
import unittest

from Bio import MissingExternalDependencyError
try:
    # Skip the test if reportlab is not installed
    import reportlab as r
    del r
except:
    raise MissingExternalDependencyError(\
        "Install ReportLab if you want to use Bio.Graphics.")
try:
    # Skip the test if reportlab is not installed
    from reportlab.graphics import renderPM
except:
    raise MissingExternalDependencyError(\
        "Install ReportLab's renderPM module if you want to create "
        "bitmaps with Bio.Graphics.")

# the stuff we're testing
from Bio.Graphics.Comparative import ComparativeScatterPlot


class ComparativeTest(unittest.TestCase):
    """Do tests for modules involved with comparing data."""
    def setUp(self):
        self.min_two_d_lists = 1
        self.max_two_d_lists = 7

        self.min_num_points = 1
        self.max_num_points = 500

        self.min_point_num = 0
        self.max_point_num = 200

    def _make_random_points(self):
        """Make a bunch of random points for testing plots."""
        plot_info = []
        num_two_d_lists = random.randrange(self.min_two_d_lists,
                                           self.max_two_d_lists)

        for two_d_list in range(num_two_d_lists):
            cur_list = []
            num_points = random.randrange(self.min_num_points,
                                          self.max_num_points)
            for point in range(num_points):
                x_point = random.randrange(self.min_point_num,
                                           self.max_point_num)
                y_point = random.randrange(self.min_point_num,
                                           self.max_point_num)

                cur_list.append((x_point, y_point))

            plot_info.append(cur_list)
        return plot_info
                
    def test_simple_scatter_plot(self):
        """Test creation of a simple PNG scatter plot."""
        compare_plot = ComparativeScatterPlot("png")
        compare_plot.display_info = self._make_random_points()

        output_file = os.path.join(os.getcwd(), "Graphics", "scatter_test.png")
        try:
            compare_plot.draw_to_file(output_file, "Testing Scatter Plots")
        # there is a bug in reportlab which occasionally generates an
        # error here.
        except IndexError:
            pass

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
