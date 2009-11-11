# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import unittest
import os.path
from Bio.SwissProt import KeyWList

class KeyWListTest(unittest.TestCase):

    def test_parse(self):
        "Parsing keywlist.txt"

        filename = os.path.join("SwissProt", "keywlist.txt")
        handle = open(filename)
        records = KeyWList.parse(handle)

        # Testing the first record
        record = records.next()
        self.assertEqual(record["ID"], "2Fe-2S.")
        self.assertEqual(record["AC"], "KW-0001")
        self.assertEqual(record["DE"], "Protein which contains at least one 2Fe-2S iron-sulfur cluster: 2 iron atoms complexed to 2 inorganic sulfides and 4 sulfur atoms of cysteines from the protein.")
        self.assertEqual(record["SY"], "Fe2S2; [2Fe-2S] cluster; [Fe2S2] cluster; Fe2/S2 (inorganic) cluster; Di-mu-sulfido-diiron; 2 iron, 2 sulfur cluster binding.")
        self.assertEqual(len(record["GO"]), 1)
        self.assertEqual(record["GO"], ["GO:0051537; 2 iron, 2 sulfur cluster binding"])
        self.assertEqual(len(record["HI"]), 2)
        self.assertEqual(record["HI"][0], "Ligand: Iron; Iron-sulfur; 2Fe-2S.")
        self.assertEqual(record["HI"][1], "Ligand: Metal-binding; 2Fe-2S.")
        self.assertEqual(record["CA"], "Ligand.")

        # Testing the second record
        record = records.next()
        self.assertEqual(record["IC"], "Molecular function.")
        self.assertEqual(record["AC"], "KW-9992")
        self.assertEqual(record["DE"], "Keywords assigned to proteins due to their particular molecular function.")

        # Testing the third record
        record = records.next()
        self.assertEqual(record["ID"], "Zymogen.")
        self.assertEqual(record["AC"], "KW-0865")
        self.assertEqual(record["DE"], "The enzymatically inactive precursor of mostly proteolytic enzymes.")
        self.assertEqual(record["SY"], "Proenzyme.")
        self.assertEqual(len(record["HI"]), 1)
        self.assertEqual(record["HI"][0], "PTM: Zymogen.")
        self.assertEqual(record["CA"], "PTM.")

        handle.close()

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
