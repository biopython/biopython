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
        """Test if Motif can parse MEME output files (first test)
        """
        from Bio.Alphabet import IUPAC
        handle = open("Motif/meme.out")
        parser = Motif.MEMEParser()
        record = parser.parse(handle)
        self.assertEqual(record.version, '3.5.7')
        self.assertEqual(record.datafile, 'test.fa')
        self.assertEqual(record.alphabet, IUPAC.unambiguous_dna)
        self.assertEqual(len(record.sequence_names), 10)
        self.assertEqual(record.sequence_names[0], 'SEQ1;')
        self.assertEqual(record.sequence_names[1], 'SEQ2;')
        self.assertEqual(record.sequence_names[2], 'SEQ3;')
        self.assertEqual(record.sequence_names[3], 'SEQ4;')
        self.assertEqual(record.sequence_names[4], 'SEQ5;')
        self.assertEqual(record.sequence_names[5], 'SEQ6;')
        self.assertEqual(record.sequence_names[6], 'SEQ7;')
        self.assertEqual(record.sequence_names[7], 'SEQ8;')
        self.assertEqual(record.sequence_names[8], 'SEQ9;')
        self.assertEqual(record.sequence_names[9], 'SEQ10;')
        self.assertEqual(record.command, 'meme test.fa -dna -w 10 -dir /home/bartek/MetaMotif/meme')
        self.assertEqual(len(record.motifs), 1)
        motif = record.motifs[0]
        self.assertEqual(motif.num_occurrences, 10)
        self.assertAlmostEqual(motif.evalue, 1.1e-22)
        self.assertEqual(motif.alphabet, IUPAC.unambiguous_dna)
        self.assertEqual(motif.name, "Motif 1")
        self.assertEqual(len(motif.instances), 10)
        self.assertAlmostEqual(motif.instances[0].pvalue, 8.71e-07)
        self.assertAlmostEqual(motif.instances[1].pvalue, 8.71e-07)
        self.assertAlmostEqual(motif.instances[2].pvalue, 8.71e-07)
        self.assertAlmostEqual(motif.instances[3].pvalue, 8.71e-07)
        self.assertAlmostEqual(motif.instances[4].pvalue, 8.71e-07)
        self.assertAlmostEqual(motif.instances[5].pvalue, 8.71e-07)
        self.assertAlmostEqual(motif.instances[6].pvalue, 8.71e-07)
        self.assertAlmostEqual(motif.instances[7].pvalue, 8.71e-07)
        self.assertAlmostEqual(motif.instances[8].pvalue, 8.71e-07)
        self.assertAlmostEqual(motif.instances[9].pvalue, 8.71e-07)
        self.assertEqual(motif.instances[0].sequence_name, 'SEQ10;')
        self.assertEqual(motif.instances[1].sequence_name, 'SEQ9;')
        self.assertEqual(motif.instances[2].sequence_name, 'SEQ8;')
        self.assertEqual(motif.instances[3].sequence_name, 'SEQ7;')
        self.assertEqual(motif.instances[4].sequence_name, 'SEQ6;')
        self.assertEqual(motif.instances[5].sequence_name, 'SEQ5;')
        self.assertEqual(motif.instances[6].sequence_name, 'SEQ4;')
        self.assertEqual(motif.instances[7].sequence_name, 'SEQ3;')
        self.assertEqual(motif.instances[8].sequence_name, 'SEQ2;')
        self.assertEqual(motif.instances[9].sequence_name, 'SEQ1;')
        self.assertEqual(motif.instances[0].start, 3)
        self.assertEqual(motif.instances[1].start, 93)
        self.assertEqual(motif.instances[2].start, 172)
        self.assertEqual(motif.instances[3].start, 177)
        self.assertEqual(motif.instances[4].start, 105)
        self.assertEqual(motif.instances[5].start, 185)
        self.assertEqual(motif.instances[6].start, 173)
        self.assertEqual(motif.instances[7].start, 112)
        self.assertEqual(motif.instances[8].start, 172)
        self.assertEqual(motif.instances[9].start, 52)
        self.assertEqual(motif.instances[0].strand, '+')
        self.assertEqual(motif.instances[1].strand, '+')
        self.assertEqual(motif.instances[2].strand, '+')
        self.assertEqual(motif.instances[3].strand, '+')
        self.assertEqual(motif.instances[4].strand, '+')
        self.assertEqual(motif.instances[5].strand, '+')
        self.assertEqual(motif.instances[6].strand, '+')
        self.assertEqual(motif.instances[7].strand, '+')
        self.assertEqual(motif.instances[8].strand, '+')
        self.assertEqual(motif.instances[9].strand, '+')
        self.assertEqual(motif.instances[0].length, 10)
        self.assertEqual(motif.instances[1].length, 10)
        self.assertEqual(motif.instances[2].length, 10)
        self.assertEqual(motif.instances[3].length, 10)
        self.assertEqual(motif.instances[4].length, 10)
        self.assertEqual(motif.instances[5].length, 10)
        self.assertEqual(motif.instances[6].length, 10)
        self.assertEqual(motif.instances[7].length, 10)
        self.assertEqual(motif.instances[8].length, 10)
        self.assertEqual(motif.instances[9].length, 10)
        self.assertEqual(motif.instances[0].motif_name, 'Motif 1')
        self.assertEqual(motif.instances[1].motif_name, 'Motif 1')
        self.assertEqual(motif.instances[2].motif_name, 'Motif 1')
        self.assertEqual(motif.instances[3].motif_name, 'Motif 1')
        self.assertEqual(motif.instances[4].motif_name, 'Motif 1')
        self.assertEqual(motif.instances[5].motif_name, 'Motif 1')
        self.assertEqual(motif.instances[6].motif_name, 'Motif 1')
        self.assertEqual(motif.instances[7].motif_name, 'Motif 1')
        self.assertEqual(motif.instances[8].motif_name, 'Motif 1')
        self.assertEqual(motif.instances[9].motif_name, 'Motif 1')
        self.assertEqual(motif.instances[0].alphabet, IUPAC.unambiguous_dna)
        self.assertEqual(motif.instances[1].alphabet, IUPAC.unambiguous_dna)
        self.assertEqual(motif.instances[2].alphabet, IUPAC.unambiguous_dna)
        self.assertEqual(motif.instances[3].alphabet, IUPAC.unambiguous_dna)
        self.assertEqual(motif.instances[4].alphabet, IUPAC.unambiguous_dna)
        self.assertEqual(motif.instances[5].alphabet, IUPAC.unambiguous_dna)
        self.assertEqual(motif.instances[6].alphabet, IUPAC.unambiguous_dna)
        self.assertEqual(motif.instances[7].alphabet, IUPAC.unambiguous_dna)
        self.assertEqual(motif.instances[8].alphabet, IUPAC.unambiguous_dna)
        self.assertEqual(motif.instances[9].alphabet, IUPAC.unambiguous_dna)
        self.assertEqual(motif.instances[0].data, "CTCAATCGTA")
        self.assertEqual(motif.instances[1].data, "CTCAATCGTA")
        self.assertEqual(motif.instances[2].data, "CTCAATCGTA")
        self.assertEqual(motif.instances[3].data, "CTCAATCGTA")
        self.assertEqual(motif.instances[4].data, "CTCAATCGTA")
        self.assertEqual(motif.instances[5].data, "CTCAATCGTA")
        self.assertEqual(motif.instances[6].data, "CTCAATCGTA")
        self.assertEqual(motif.instances[7].data, "CTCAATCGTA")
        self.assertEqual(motif.instances[8].data, "CTCAATCGTA")
        self.assertEqual(motif.instances[9].data, "CTCAATCGTA")
        handle.close()

    def test_meme_parser_2(self):
        """Test if Motif can parse MEME output files (second test)
        """
        from Bio.Alphabet import IUPAC
        handle = open("Motif/meme.dna.oops.txt")
        parser = Motif.MEMEParser()
        record = parser.parse(handle)
        self.assertEqual(record.version, '3.0')
        self.assertEqual(record.datafile, 'INO_up800.s')
        self.assertEqual(record.alphabet, IUPAC.unambiguous_dna)
        self.assertEqual(len(record.sequence_names), 7)
        self.assertEqual(record.sequence_names[0], 'CHO1')
        self.assertEqual(record.sequence_names[1], 'CHO2')
        self.assertEqual(record.sequence_names[2], 'FAS1')
        self.assertEqual(record.sequence_names[3], 'FAS2')
        self.assertEqual(record.sequence_names[4], 'ACC1')
        self.assertEqual(record.sequence_names[5], 'INO1')
        self.assertEqual(record.sequence_names[6], 'OPI3')
        self.assertEqual(record.command, 'meme -mod oops -dna -revcomp -nmotifs 2 -bfile yeast.nc.6.freq INO_up800.s')
        self.assertEqual(len(record.motifs), 2)
        motif = record.motifs[0]
        self.assertEqual(motif.num_occurrences, 7)
        self.assertAlmostEqual(motif.evalue, 0.2)
        self.assertEqual(motif.alphabet, IUPAC.unambiguous_dna)
        self.assertEqual(motif.name, "Motif 1")
        self.assertEqual(len(motif.instances), 7)
        self.assertAlmostEqual(motif.instances[0].pvalue, 1.85e-08)
        self.assertAlmostEqual(motif.instances[1].pvalue, 1.85e-08)
        self.assertAlmostEqual(motif.instances[2].pvalue, 1.52e-07)
        self.assertAlmostEqual(motif.instances[3].pvalue, 2.52e-07)
        self.assertAlmostEqual(motif.instances[4].pvalue, 4.23e-07)
        self.assertAlmostEqual(motif.instances[5].pvalue, 9.43e-07)
        self.assertAlmostEqual(motif.instances[6].pvalue, 3.32e-06)
        self.assertEqual(motif.instances[0].sequence_name, 'INO1')
        self.assertEqual(motif.instances[1].sequence_name, 'FAS1')
        self.assertEqual(motif.instances[2].sequence_name, 'ACC1')
        self.assertEqual(motif.instances[3].sequence_name, 'CHO2')
        self.assertEqual(motif.instances[4].sequence_name, 'CHO1')
        self.assertEqual(motif.instances[5].sequence_name, 'FAS2')
        self.assertEqual(motif.instances[6].sequence_name, 'OPI3')
        self.assertEqual(motif.instances[0].strand, '-')
        self.assertEqual(motif.instances[1].strand, '+')
        self.assertEqual(motif.instances[2].strand, '+')
        self.assertEqual(motif.instances[3].strand, '+')
        self.assertEqual(motif.instances[4].strand, '+')
        self.assertEqual(motif.instances[5].strand, '+')
        self.assertEqual(motif.instances[6].strand, '+')
        self.assertEqual(motif.instances[0].length, 12)
        self.assertEqual(motif.instances[1].length, 12)
        self.assertEqual(motif.instances[2].length, 12)
        self.assertEqual(motif.instances[3].length, 12)
        self.assertEqual(motif.instances[4].length, 12)
        self.assertEqual(motif.instances[5].length, 12)
        self.assertEqual(motif.instances[6].length, 12)
        self.assertEqual(motif.instances[0].start, 620)
        self.assertEqual(motif.instances[1].start,  95)
        self.assertEqual(motif.instances[2].start,  83)
        self.assertEqual(motif.instances[3].start, 354)
        self.assertEqual(motif.instances[4].start, 611)
        self.assertEqual(motif.instances[5].start, 567)
        self.assertEqual(motif.instances[6].start, 340)
        self.assertEqual(motif.instances[0].data, "TTCACATGCCGC")
        self.assertEqual(motif.instances[1].data, "TTCACATGCCGC")
        self.assertEqual(motif.instances[2].data, "TTCACATGGCCC")
        self.assertEqual(motif.instances[3].data, "TTCTCATGCCGC")
        self.assertEqual(motif.instances[4].data, "TTCACACGGCAC")
        self.assertEqual(motif.instances[5].data, "TTCACATGCTAC")
        self.assertEqual(motif.instances[6].data, "TTCAGATCGCTC")
        motif = record.motifs[1]
        self.assertEqual(motif.num_occurrences, 7)
        self.assertAlmostEqual(motif.evalue, 110)
        self.assertEqual(motif.alphabet, IUPAC.unambiguous_dna)
        self.assertEqual(motif.name, "Motif 2")
        self.assertEqual(len(motif.instances), 7)
        self.assertAlmostEqual(motif.instances[0].pvalue, 3.24e-07)
        self.assertAlmostEqual(motif.instances[1].pvalue, 3.24e-07)
        self.assertAlmostEqual(motif.instances[2].pvalue, 3.24e-07)
        self.assertAlmostEqual(motif.instances[3].pvalue, 5.29e-06)
        self.assertAlmostEqual(motif.instances[4].pvalue, 6.25e-06)
        self.assertAlmostEqual(motif.instances[5].pvalue, 8.48e-06)
        self.assertAlmostEqual(motif.instances[6].pvalue, 8.48e-06)
        self.assertEqual(motif.instances[0].sequence_name, 'OPI3')
        self.assertEqual(motif.instances[1].sequence_name, 'ACC1')
        self.assertEqual(motif.instances[2].sequence_name, 'CHO1')
        self.assertEqual(motif.instances[3].sequence_name, 'INO1')
        self.assertEqual(motif.instances[4].sequence_name, 'FAS1')
        self.assertEqual(motif.instances[5].sequence_name, 'FAS2')
        self.assertEqual(motif.instances[6].sequence_name, 'CHO2')
        self.assertEqual(motif.instances[0].strand, '-')
        self.assertEqual(motif.instances[1].strand, '+')
        self.assertEqual(motif.instances[2].strand, '-')
        self.assertEqual(motif.instances[3].strand, '-')
        self.assertEqual(motif.instances[4].strand, '+')
        self.assertEqual(motif.instances[5].strand, '-')
        self.assertEqual(motif.instances[6].strand, '-')
        self.assertEqual(motif.instances[0].length, 10)
        self.assertEqual(motif.instances[1].length, 10)
        self.assertEqual(motif.instances[2].length, 10)
        self.assertEqual(motif.instances[3].length, 10)
        self.assertEqual(motif.instances[4].length, 10)
        self.assertEqual(motif.instances[5].length, 10)
        self.assertEqual(motif.instances[6].length, 10)
        self.assertEqual(motif.instances[0].start, 186)
        self.assertEqual(motif.instances[1].start, 232)
        self.assertEqual(motif.instances[2].start, 559)
        self.assertEqual(motif.instances[3].start, 283)
        self.assertEqual(motif.instances[4].start,  44)
        self.assertEqual(motif.instances[5].start, 185)
        self.assertEqual(motif.instances[6].start, 413)
        self.assertEqual(motif.instances[0].data, "TCTGGCACAG")
        self.assertEqual(motif.instances[1].data, "TCTGGCACAG")
        self.assertEqual(motif.instances[2].data, "TCTGGCACAG")
        self.assertEqual(motif.instances[3].data, "GCGGGCGCAG")
        self.assertEqual(motif.instances[4].data, "GCAGGCACGG")
        self.assertEqual(motif.instances[5].data, "TCTGGCACTC")
        self.assertEqual(motif.instances[6].data, "TCTGGCATCG")
        handle.close()

    def test_meme_parser_3(self):
        """Test if Motif can parse MEME output files (third test)
        """
        handle = open("Motif/meme.protein.oops.txt")
        parser = Motif.MEMEParser()
        record = parser.parse(handle)
        handle.close()

    def test_meme_parser_4(self):
        """Test if Motif can parse MEME output files (fourth test)
        """
        handle = open("Motif/meme.protein.tcm.txt")
        parser = Motif.MEMEParser()
        record = parser.parse(handle)
        handle.close()

 
class TestMAST(unittest.TestCase):

    # The MAST parser currently fails; commented out the tests below.

    def test_mast_parser_1(self):
        """Test if Motif can parse MAST output files (first test)
        """
        handle = open("Motif/mast.dna.oops.txt")
        parser = Motif.MASTParser()
        # record = parser.parse(handle)
        handle.close()

    def test_mast_parser_2(self):
        """Test if Motif can parse MAST output files (second test)
        """
        handle = open("Motif/mast.protein.oops.txt")
        parser = Motif.MASTParser()
        # record = parser.parse(handle)
        handle.close()

    def test_mast_parser_3(self):
        """Test if Motif can parse MAST output files (third test)
        """
        handle = open("Motif/mast.protein.tcm.txt")
        parser = Motif.MASTParser()
        # record = parser.parse(handle)
        handle.close()


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
