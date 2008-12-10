# Copyright 2008 by Bartek Wilczynski.  All rights reserved.
# Adapted from test_Mymodule.py by Jeff Chang
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
import sys
import unittest

from Bio import Motif

def run_tests(argv):
    test_suite = testing_suite()
    runner = unittest.TextTestRunner(sys.stdout, verbosity = 2)
    runner.run(test_suite)

def testing_suite():
    """Generate the suite of tests.
    """
    test_suite = unittest.TestSuite()

    test_loader = unittest.TestLoader()
    test_loader.testMethodPrefix = 't_'
    tests = [MotifTestsBasic]
    
    for test in tests:
        cur_suite = test_loader.loadTestsFromTestCase(test)
        test_suite.addTest(cur_suite)

    return test_suite

class MotifTestsBasic(unittest.TestCase):
    def setUp(self):
        self.ACin = open("Motif/alignace.out")
        self.MEMEin = open("Motif/meme.out")
        self.PFMin = open("Motif/SRF.pfm")
        self.SITESin = open("Motif/Arnt.sites")
        self.TFout = "Motif/tf.out"
        self.FAout = "Motif/fa.out"
        self.PFMout = "Motif/fa.out"
        from Bio.Seq import Seq
        self.m=Motif.Motif()
        self.m.add_instance(Seq("ATATA",self.m.alphabet))
        
    def tearDown(self):
        self.ACin.close()
        self.MEMEin.close()
        self.PFMin.close()
        self.SITESin.close()
        if os.path.exists(self.TFout):
            os.remove(self.TFout)
        if os.path.exists(self.FAout):
            os.remove(self.FAout)

    def t_alignace_parsing(self):
        """Test to be sure that Motif can parse AlignAce output files.
        """
        parser= Motif.AlignAceParser()
        record=parser.parse(self.ACin)
        assert len(record.motifs)==16
        
    def t_meme_parsing(self):
        """Test to be sure that Motif can parse MEME output files.
        """
        parser= Motif.MEMEParser()
        record=parser.parse(self.MEMEin)
        assert len(record.motifs)==1

    def t_pfm_parsing(self):
        """Test to be sure that Motif can parse pfm  files.
        """
        motif= Motif.Motif()
        motif.from_jaspar_pfm(self.PFMin)
        assert motif.length==12

    def t_sites_parsing(self):
        """Test to be sure that Motif can parse sites files.
        """
        motif= Motif.Motif()
        motif.from_jaspar_sites(self.SITESin)
        assert motif.length==6

    def t_FAoutput(self):
        """Ensure that we can write proper FASTA output files.
        """
        output_handle = open(self.FAout, "w")
        output_handle.write(self.m.to_fasta())
        output_handle.close()

    def t_TFoutput(self):
        """Ensure that we can write proper TransFac output files.
        """
        output_handle = open(self.TFout, "w")
        output_handle.write(self.m.to_transfac())
        output_handle.close()

    def t_pfm_output(self):
        """Ensure that we can write proper pfm output files.
        """
        output_handle = open(self.PFMout, "w")
        output_handle.write(self.m.to_transfac())
        output_handle.close()
        
        

if __name__ == "__main__":
    sys.exit(run_tests(sys.argv))
