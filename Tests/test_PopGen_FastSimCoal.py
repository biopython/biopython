# Copyright 2014 by Melissa Gymrek <mgymrek@mit.edu>.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
import shutil
import unittest
import warnings

from Bio import MissingExternalDependencyError
from Bio import BiopythonDeprecationWarning

with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonDeprecationWarning)
    from Bio.PopGen.SimCoal.Controller import FastSimCoalController

# Tests fastsimcoal related code. Note: this case requires fastsimcoal 2.51

found = False
for path in os.environ['PATH'].split(os.pathsep):
    try:
        for filename in os.listdir(path):
            if filename == "fsc252" \
                    or (filename.lower() == "fsc252.exe"):
                found = True
                fastsimcoal_dir = path
    except os.error:
        pass  # Path doesn't exist - correct to pass
if not found:
    raise MissingExternalDependencyError(
        "Install fastsimcoal if you want to use "
        "Bio.PopGen.SimCoal.Controller.FastSimCoalController.")


class AppTest(unittest.TestCase):
    """Tests fastsimcoal execution via biopython.
    """
    def setUp(self):
        self.tidy()

    def tearDown(self):
        self.tidy()

    def tidy(self):
        if not os.path.isdir(os.path.join('PopGen', 'simple')):
            # Unit test must have failed to invoke fastsimcoal251,
            # and thus it never created the directory.
            return
        shutil.rmtree(os.path.join('PopGen', 'simple'))

    def test_fastsimcoal(self):
        """Test fastsimcoal execution.
        """
        ctrl = FastSimCoalController(fastsimcoal_dir=fastsimcoal_dir)
        ctrl.run_fastsimcoal('simple.par', 50, par_dir='PopGen')
        assert os.path.isdir(os.path.join('PopGen', 'simple')), \
            "Output directory not created!"
        assert(len(os.listdir(os.path.join('PopGen', 'simple'))) == 52)

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
