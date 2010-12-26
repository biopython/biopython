# Copyright 1999 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


import os
import unittest

from Bio.ExPASy import Enzyme


class TestEnzyme(unittest.TestCase):

    def test_lipoprotein(self):
        "Parsing ENZYME record for lipoprotein lipase (3.1.1.34)"
        filename = os.path.join( 'Enzymes', 'lipoprotein.txt')
        handle = open(filename)
        record = Enzyme.read(handle)
        handle.close()
        self.assertEqual(record["ID"], "3.1.1.34")
        self.assertEqual(record["DE"], "Lipoprotein lipase.")
        self.assertEqual(len(record["AN"]), 3)
        self.assertEqual(record["AN"][0], "Clearing factor lipase.")
        self.assertEqual(record["AN"][1], "Diacylglycerol lipase.")
        self.assertEqual(record["AN"][2], "Diglyceride lipase.")
        self.assertEqual(record["CA"], "Triacylglycerol + H(2)O = diacylglycerol + a carboxylate.")
        self.assertEqual(record["CC"][0], 'Hydrolyzes triacylglycerols in chylomicrons and very low-density lipoproteins (VLDL).')
        self.assertEqual(record["CC"][1], "Also hydrolyzes diacylglycerol.")
        self.assertEqual(record['PR'], ["PDOC00110"])
        self.assertEqual(record["DR"][0], ["P11151", "LIPL_BOVIN"])
        self.assertEqual(record["DR"][1], ["P11153", "LIPL_CAVPO"])
        self.assertEqual(record["DR"][2], ["P11602", "LIPL_CHICK"])
        self.assertEqual(record["DR"][3], ["P55031", "LIPL_FELCA"])
        self.assertEqual(record["DR"][4], ["P06858", "LIPL_HUMAN"])
        self.assertEqual(record["DR"][5], ["P11152", "LIPL_MOUSE"])
        self.assertEqual(record["DR"][6], ["O46647", "LIPL_MUSVI"])
        self.assertEqual(record["DR"][7], ["P49060", "LIPL_PAPAN"])
        self.assertEqual(record["DR"][8], ["P49923", "LIPL_PIG"])
        self.assertEqual(record["DR"][9], ["Q06000", "LIPL_RAT"])
        self.assertEqual(record["DR"][10], ["Q29524", "LIPL_SHEEP"])

    def test_proline(self):
        "Parsing ENZYME record for proline racemase (5.1.1.4)"
        filename = os.path.join( 'Enzymes', 'proline.txt')
        handle = open(filename)
        record = Enzyme.read(handle)
        handle.close()
        self.assertEqual(record["ID"], "5.1.1.4")
        self.assertEqual(record["DE"], "Proline racemase.")
        self.assertEqual(record["CA"], "L-proline = D-proline.")
        self.assertEqual(len(record["DR"]), 9)
        self.assertEqual(record["DR"][0], ["Q17ZY4", "PRAC_CLOD6"])
        self.assertEqual(record["DR"][1], ["A8DEZ8", "PRAC_CLODI"])
        self.assertEqual(record["DR"][2], ["Q4DA80", "PRCMA_TRYCR"])
        self.assertEqual(record["DR"][3], ["Q868H8", "PRCMB_TRYCR"])
        self.assertEqual(record["DR"][4], ["Q3SX04", "PRCM_BOVIN"])
        self.assertEqual(record["DR"][5], ["Q96EM0", "PRCM_HUMAN"])
        self.assertEqual(record["DR"][6], ["Q9CXA2", "PRCM_MOUSE"])
        self.assertEqual(record["DR"][7], ["Q5RC28", "PRCM_PONAB"])
        self.assertEqual(record["DR"][8], ["Q66II5", "PRCM_XENTR"])

    def test_valine(self):
        "Parsing ENZYME record for valine decarboxylase (4.1.1.14)"
        filename = os.path.join( 'Enzymes', 'valine.txt')
        handle = open(filename)
        record = Enzyme.read(handle)
        handle.close()
        self.assertEqual(record["ID"], "4.1.1.14")
        self.assertEqual(record["DE"], "Valine decarboxylase.")
        self.assertEqual(record["CA"], "L-valine = 2-methylpropanamine + CO(2).")
        self.assertEqual(record["CF"], "Pyridoxal 5'-phosphate.")
        self.assertEqual(record["CC"], ["Also acts on L-leucine."])
        self.assertEqual(len(record["DR"]), 0)

    def test_lactate(self):
        "Parsing ENZYME record for lactate racemase (5.1.2.1)"
        filename = os.path.join( 'Enzymes', 'lactate.txt')
        handle = open(filename)
        record = Enzyme.read(handle)
        handle.close()
        self.assertEqual(record["ID"], "5.1.2.1")
        self.assertEqual(record["DE"], "Lactate racemase.")
        self.assertEqual(len(record["AN"]), 3)
        self.assertEqual(record["AN"][0], "Hydroxyacid racemase.")
        self.assertEqual(record["AN"][1], "Lactic acid racemase.")
        self.assertEqual(record["AN"][2], "Lacticoracemase.")
        self.assertEqual(record["CA"], "(S)-lactate = (R)-lactate.")
        self.assertEqual(len(record["DR"]), 0)

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
