# Copyright 2007 by Tiago Antao <tiagoantao@gmail.com>.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
import unittest
from Bio.PopGen import SimCoal
from Bio.PopGen.SimCoal.Controller import SimCoalController, FastSimCoalController
from Bio import MissingExternalDependencyError

#Tests simcoal related code. Note: this case requires simcoal and fastsimcoal
#test_PopGen_SimCoal_nodepend tests code that does not require simcoal or fastsimcoal

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
    """Tests simcoal execution via biopython.
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

    def test_simcoal(self):
        """Test simcoal execution.
        """
        simcoal_dir = CheckForExecutable("simcoal2")
        ctrl = SimCoalController(simcoal_dir)
        ctrl.run_simcoal('simple.par', 50, par_dir = 'PopGen')
        assert os.path.isdir(os.path.join('PopGen', 'simple')), \
               "Output directory not created!"
        assert( len(os.listdir(os.path.join('PopGen', 'simple'))) == 52)
        self.tidy()

    def test_fastsimcoal(self):
        """Test simcoal execution.
        """
        fastsimcoal_dir = CheckForExecutable("fastsimcoal21")
        print fastsimcoal_dir
        ctrl = FastSimCoalController(fastsimcoal_dir=fastsimcoal_dir)
        ctrl.run_fastsimcoal('simple.par', 50, par_dir = 'PopGen')
        assert os.path.isdir(os.path.join('PopGen', 'simple')), \
               "Output directory not created!"
        assert( len(os.listdir(os.path.join('PopGen', 'simple'))) == 52)
        self.tidy()

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
