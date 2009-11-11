# Copyright 2006 by Tiago Antao <tiagoantao@gmail.com>.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import commands
import os
import shutil
import tempfile
import unittest
from Bio.PopGen import FDist
from Bio.PopGen.FDist import Controller
from Bio import MissingExternalDependencyError

#Tests fdist related code. Note: this case requires fdist
#test_PopGen_FDist_nodepend tests code that does not require fdist

found = False
for path in os.environ['PATH'].split(os.pathsep):
    try:
        list = os.listdir(path)
        for file in os.listdir(path):
            if file.startswith('fdist2'):
                found = True
    except os.error:
        pass #Path doesn't exist - correct to pass
if not found:
    raise MissingExternalDependencyError(\
        "Install FDist if you want to use Bio.PopGen.FDist.")


class AppTest(unittest.TestCase):
    """Tests the fdist suite of applications.
    """
    def _copyfile(self, inname, outname):
        shutil.copyfile(
            'PopGen' + os.sep + inname,
            self.dirname + os.sep + outname)

    def setUp(self):
        self.ctrl = Controller.FDistController()
        self.dirname = tempfile.mkdtemp()
        self._copyfile('data_fst_outfile', 'data_fst_outfile')
        self._copyfile('fdist1', 'infile')
        self._copyfile('out.dat', 'out.dat')
        self._copyfile('out.cpl', 'out.cpl')

    def tearDown(self):
        #Not sure how exactly, but its possible the temp directory
        #may not (still) exist.
        if os.path.isdir(self.dirname):
            for file in os.listdir(self.dirname):
                os.remove(self.dirname + os.sep + file)
            os.rmdir(self.dirname)

    def test_datacal(self):
        """Test datacal execution.
        """
        fst, samp_size = self.ctrl.run_datacal(data_dir = self.dirname)
        assert (fst - 0.44 < 0.01)
        assert (samp_size == 11)

    def test_fdist(self):
        """Test fdist execution.
        """
        #The number of simulations in real life should be at least 10000,
        #see the fdist2 documentation.
        fst = self.ctrl.run_fdist(npops = 15, nsamples = 10, fst = 0.1,
                sample_size = 20, mut = 0, num_sims = 100,
                data_dir = self.dirname)
        assert(abs(fst - 0.1) < 0.02) #Stochastic result...

    def test_fdist_force_fst(self):
        """Test fdist execution approximating Fst.
        """
        #The number of simulations in real life should be at least 10000,
        #see the fdist2 documentation.
        fst = self.ctrl.run_fdist_force_fst(npops = 15, nsamples = 10,
                fst = 0.1,
                sample_size = 20, mut = 0, num_sims = 100,
                data_dir = self.dirname)
        assert(abs(fst - 0.09) < 0.05) #Stochastic result...

    def test_cplot(self):
        """Test cplot execution.
        """
        cpl_interval =self.ctrl.run_cplot(data_dir = self.dirname)
        assert(len(cpl_interval) == 8)

    def test_pv(self):
        """Test pv execution.
        """
        pv_data = self.ctrl.run_pv(data_dir = self.dirname)
        assert(len(pv_data) == 4)


if __name__ == "__main__":
    print "Running fdist tests, which might take some time, please wait"
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
