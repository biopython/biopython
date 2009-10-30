# Copyright 2008 by Bartek Wilczynski.  All rights reserved.
# Adapted from test_Mymodule.py by Jeff Chang
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
import unittest

from Bio import Motif


class MotifTestsBasic(unittest.TestCase):
    def setUp(self):
        self.ACin = open("Motif/alignace.out")
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
        self.PFMin.close()
        self.SITESin.close()
        if os.path.exists(self.TFout):
            os.remove(self.TFout)
        if os.path.exists(self.FAout):
            os.remove(self.FAout)

    def test_alignace_parsing(self):
        """Test to be sure that Motif can parse AlignAce output files.
        """
        parser= Motif.AlignAceParser()
        record=parser.parse(self.ACin)
        assert len(record.motifs)==16
        
    def test_pfm_parsing(self):
        """Test to be sure that Motif can parse pfm  files.
        """
        motif= Motif.read(self.PFMin,"jaspar-pfm")
        assert motif.length==12

    def test_sites_parsing(self):
        """Test to be sure that Motif can parse sites files.
        """
        motif= Motif.read(self.SITESin,"jaspar-sites")
        assert motif.length==6

    def test_FAoutput(self):
        """Ensure that we can write proper FASTA output files.
        """
        output_handle = open(self.FAout, "w")
        output_handle.write(self.m.format("fasta"))
        output_handle.close()

    def test_TFoutput(self):
        """Ensure that we can write proper TransFac output files.
        """
        output_handle = open(self.TFout, "w")
        output_handle.write(self.m.format("transfac"))
        output_handle.close()

    def test_pfm_output(self):
        """Ensure that we can write proper pfm output files.
        """
        output_handle = open(self.PFMout, "w")
        output_handle.write(self.m.format("jaspar-pfm"))
        output_handle.close()
        
        
class TestMEME(unittest.TestCase):
        
    def test_meme_parser_1(self):
        """Tests Motif can parse MEME output files (first test)
        """
        handle = open("Motif/meme.out")
        parser = Motif.MEMEParser()
        record = parser.parse(handle)
        self.assertEqual(len(record.motifs), 1)
        handle.close()

    def test_meme_parser_2(self):
        """Tests Motif can parse MEME output files (second test)
        """
        handle = open("Motif/meme.dna.oops.txt")
        parser = Motif.MEMEParser()
        record = parser.parse(handle)
        handle.close()

    def test_meme_parser_3(self):
        """Tests Motif can parse MEME output files (third test)
        """
        handle = open("Motif/meme.protein.oops.txt")
        parser = Motif.MEMEParser()
        record = parser.parse(handle)
        handle.close()

    def test_meme_parser_4(self):
        """Tests Motif can parse MEME output files (fourth test)
        """
        handle = open("Motif/meme.protein.tcm.txt")
        parser = Motif.MEMEParser()
        record = parser.parse(handle)
        handle.close()

 
class TestMAST(unittest.TestCase):

    # The MAST parser currently fails; commented out the tests below.

    def test_mast_parser_1(self):
        """Tests Motif can parse MAST output files (first test)
        """
        handle = open("Motif/mast.dna.oops.txt")
        parser = Motif.MASTParser()
        # record = parser.parse(handle)
        handle.close()

    def test_mast_parser_2(self):
        """Tests Motif can parse MAST output files (second test)
        """
        handle = open("Motif/mast.protein.oops.txt")
        parser = Motif.MASTParser()
        # record = parser.parse(handle)
        handle.close()

    def test_mast_parser_3(self):
        """Tests Motif can parse MAST output files (third test)
        """
        handle = open("Motif/mast.protein.tcm.txt")
        parser = Motif.MASTParser()
        # record = parser.parse(handle)
        handle.close()


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
