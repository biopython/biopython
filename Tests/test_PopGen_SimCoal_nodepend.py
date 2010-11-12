# Copyright 2006 by Tiago Antao <tiagoantao@gmail.com>.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
import unittest
from Bio.PopGen import SimCoal
from Bio.PopGen.SimCoal.Template import generate_simcoal_from_template

#Tests simcoal related code. Note: this case doesn't require simcoal
#test_PopGen_SimCoal tests code that requires simcoal


class TemplateTest(unittest.TestCase):
    def test_template_full(self):
        """Full template creation test
        """
        generate_simcoal_from_template('simple',
            [(1, [('SNP', [24, 0.0005, 0.0])])],
            [('sample_size', [30]),
             ('pop_size', [100])],
            'PopGen')
        #Confirm the files match (ignoring any switch of line endings
        #possible if the input file used a different OS convention)
        handle = open(os.path.join('PopGen', 'simple.par'), "rU")
        old = handle.readlines()
        handle.close()
        handle = open(os.path.join('PopGen', 'simple_100_30.par'))
        new = handle.readlines()
        handle.close()
        assert old==new, "Error - Old:\n%s\n\nNew:\n%s\n" % (old, new)
        #assert(os.stat('PopGen' + os.sep + 'simple.par').st_size ==
        #       os.stat('PopGen' + os.sep + 'simple_100_30.par').st_size)

    def tearDown(self):
        if os.path.isfile(os.path.join('PopGen', 'tmp.par')):
            #This is a temp file create by the Bio.PopGen.SimCoal.Template
            #function generate_simcoal_from_template
            os.remove(os.path.join('PopGen', 'tmp.par'))
        if os.path.isfile(os.path.join('PopGen', 'simple_100_30.par')):
            #This won't exist if the template generation failed:
            os.remove(os.path.join('PopGen', 'simple_100_30.par'))

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
