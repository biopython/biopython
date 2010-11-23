# Copyright 2010 by Tiago Antao <tiagoantao@gmail.com>.  All rights reserved.
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

#Tests DFDist related code. Note: this case requires Dfdist (four binaries)
#test_PopGen_FDist_nodepend tests code that does not require fdist2 or Dfdist

wanted = dict()
for path in os.environ['PATH'].split(os.pathsep):
    try:
        list = os.listdir(path)
        for file in os.listdir(path):
            for f in ['Dfdist', 'Ddatacal', 'pv2', 'cplot2']:
                if file == f or file.lower() == f.lower()+".exe":
                    wanted[f] = file
    except os.error:
        pass #Path doesn't exist - correct to pass
if len(wanted) != 4:
    raise MissingExternalDependencyError(\
        "Install Dfdist, Ddatacal, pv2 and cplot2 if you want to use DFDist with Bio.PopGen.FDist.")
del wanted


class AppTest(unittest.TestCase):
    """Tests the Dfdist suite of applications.
    """
    def _copyfile(self, inname, outname):
        shutil.copyfile(
            'PopGen' + os.sep + inname,
            self.dirname + os.sep + outname)

    def setUp(self):
        self.ctrl = Controller.FDistController()
        self.dirname = tempfile.mkdtemp()
        self._copyfile('data_dfst_outfile', 'data_fst_outfile')
        self._copyfile('dfdist1', 'infile')
        self._copyfile('dout.dat', 'out.dat')
        self._copyfile('dout.cpl', 'out.cpl')

    def tearDown(self):
        #Not sure how exactly, but its possible the temp directory
        #may not (still) exist.
        if os.path.isdir(self.dirname):
            for file in os.listdir(self.dirname):
                os.remove(self.dirname + os.sep + file)
            os.rmdir(self.dirname)

    def test_ddatacal(self):
        """Test Ddatacal execution.
        """
        fst, samp_size, loci, pops, F, obs = \
            self.ctrl.run_datacal(data_dir = self.dirname, version=2)
        self.assertTrue(fst - 0.23 < 0.02)
        self.assertEqual(samp_size, 32)
        self.assertEqual(loci, 300)
        self.assertEqual(pops, 2)
        self.assertTrue(F - 0.11 < 0.01)
        self.assertEqual(obs, 300)

    def test_dfdist(self):
        """Test Dfdist execution.
        """
        #The number of simulations in real life should be at least 10000,
        #see the fdist2 documentation.
        fst = self.ctrl.run_fdist(npops = 15, nsamples = 10, fst = 0.1,
                sample_size = 20, mut = 0, num_sims = 100,
                data_dir = self.dirname, is_dominant = True)
        self.assertTrue(abs(fst - 0.1) < 0.025,
                        "Stochastic result, expected %f close to 0.1" % fst)

    def atest_dfdist_force_fst(self):
        """Test dfdist execution approximating Fst.
           THIS IS TOO SLOW
        """
        #The number of simulations in real life should be at least 10000,
        #see the fdist2 documentation.
        fst = self.ctrl.run_fdist_force_fst(npops = 15, nsamples = 10,
                fst = 0.1,
                sample_size = 20, mut = 0, num_sims = 100,
                data_dir = self.dirname, is_dominant=True)
        self.assertTrue(abs(fst - 0.09) < 0.05,
                        "Stochastic result, expected %f close to 0.09" % fst)

    def test_cplot2(self):
        """Test cplot2 execution.
        """
        cpl_interval =self.ctrl.run_cplot(data_dir = self.dirname, version=2)
        self.assertEqual(len(cpl_interval), 300)

    def test_pv2(self):
        """Test pv2 execution.
        """
        pv_data = self.ctrl.run_pv(data_dir = self.dirname, version=2)
        self.assertEqual(len(pv_data), 300)


if __name__ == "__main__":
    print "Running fdist tests, which might take some time, please wait"
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
