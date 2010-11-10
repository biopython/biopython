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
from Bio import MissingPythonDependencyError

try:
    # Skip the test if reportlab is not installed
    import reportlab as r
    del r
except:
    raise MissingPythonDependencyError(
        "Install ReportLab if you want to use Bio.Graphics.")
try:
    # Skip the test if reportlab is not installed
    from reportlab.graphics import renderPM
except:
    raise MissingPythonDependencyError(\
        "Install ReportLab's renderPM module if you want to create "
        "bitmaps with Bio.Graphics.")
try:
    # Skip the test if PIL is not installed
    import Image as i
    del i
except:
    raise MissingPythonDependencyError(\
        "Install PIL (Python Imaging Library) if you want to create "
        "bitmaps with Bio.Graphics.")

from reportlab.graphics.renderPM import RenderPMError

# the stuff we're testing
from Bio.Graphics.Comparative import ComparativeScatterPlot

# We're not really using the unittest framework, because we need to
# raise the dependency error BEFORE the invidual tests in order that
# this be skipped by run_tests.py

def real_test():
    min_two_d_lists = 1
    max_two_d_lists = 7

    min_num_points = 1
    max_num_points = 500

    min_point_num = 0
    max_point_num = 200

    plot_info = []
    num_two_d_lists = random.randrange(min_two_d_lists, max_two_d_lists)

    for two_d_list in range(num_two_d_lists):
        cur_list = []
        num_points = random.randrange(min_num_points, max_num_points)
        for point in range(num_points):
            x_point = random.randrange(min_point_num, max_point_num)
            y_point = random.randrange(min_point_num, max_point_num)
            cur_list.append((x_point, y_point))
        plot_info.append(cur_list)
    
    compare_plot = ComparativeScatterPlot("png")
    compare_plot.display_info = plot_info

    output_file = os.path.join(os.getcwd(), "Graphics", "scatter_test.png")
    try:
        compare_plot.draw_to_file(output_file, "Testing Scatter Plots")
    # there is a bug in reportlab which occasionally generates an
    # error here.
    except IndexError:
        pass
    except IOError, err:
        if "encoder zip not available" in str(err):
            raise MissingExternalDependencyError(
                "Check zip encoder installed for PIL and ReportLab renderPM")
        else:
            raise err
    except RenderPMError, err :
        if str(err).startswith("Can't setFont(") :
            #TODO - can we raise the error BEFORE the unit test function
            #is run? That way it can be skipped in run_tests.py
            raise MissingExternalDependencyError(\
                "Check the fonts needed by ReportLab if you want "
                "bitmaps from Bio.Graphics\n" + str(err))
        else :
            raise err
    
    return True

#Run the actual test BEFORE the unittest stuff gets called
real_test()
               
class ComparativeTest(unittest.TestCase):
    """Do tests for modules involved with comparing data."""
    def test_simple_scatter_plot(self):
        """Test creation of a simple PNG scatter plot."""
        #Dummy method to show up via run_tests.py
        pass

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
