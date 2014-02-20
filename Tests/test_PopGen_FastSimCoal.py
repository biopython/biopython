# Copyright 2014 by Melissa Gymrek <mgymrek@mit.edu>.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
import unittest
from Bio.PopGen.SimCoal.Controller import FastSimCoalController
from Bio import MissingExternalDependencyError

#Tests simcoal related code. Note: this case requires fastsimcoal21

def CheckForExecutable(exe):
    found = False
    for path in os.environ['PATH'].split(os.pathsep):
        try:
            for filename in os.listdir(path):
                if filename == exe \
                or (filename.lower() == exe+".exe"):
                    found = True
                    exe_dir = path
        except os.error:
            pass  # Path doesn't exist - correct to pass
    if not found:
        raise MissingExternalDependencyError(
            "Install %s if you want to use Bio.PopGen.SimCoal."%exe)
    return exe_dir

class AppTest(unittest.TestCase):
    """Tests fastsimcoal execution via biopython.
    """
    def setUp(self):
        self.tidy()

    def tearDown(self):
        self.tidy()

    def tidy(self):
        if not os.path.isdir(os.path.join('PopGen', 'simple')):
            #Unit test must have failed to invoke simcaol,
            #and thus it never created the directory.
            return
        for file in os.listdir(os.path.join('PopGen', 'simple')):
            os.remove(os.sep.join(['PopGen', 'simple', file]))
        os.rmdir(os.path.join('PopGen', 'simple'))

    def test_fastsimcoal(self):
        """Test fastsimcoal execution.
        """
        fastsimcoal_dir = CheckForExecutable("fastsimcoal21")
        print fastsimcoal_dir
        ctrl = FastSimCoalController(fastsimcoal_dir=fastsimcoal_dir)
        ctrl.run_fastsimcoal('simple.par', 50, par_dir = 'PopGen')
        assert os.path.isdir(os.path.join('PopGen', 'simple')), \
               "Output directory not created!"
        assert(len(os.listdir(os.path.join('PopGen', 'simple'))) == 52)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
