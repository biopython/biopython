# Copyright 2001 by Tarjei Mikkelsen.  All rights reserved.
# Revisions copyright 2007 by Michiel de Hoon. All rights reserved.
# Revisions copyright 2017 by Peter Cock. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests the basic functionality of the KEGG parsers."""

import unittest

from Bio.KEGG import Enzyme
from Bio.KEGG import Compound
from Bio.KEGG import Map
from Bio.Pathway import System


class EnzymeTests(unittest.TestCase):
    """Tests for Bio.KEGG.Enzyme"""

    def test_sample(self):
        with open("KEGG/enzyme.sample") as handle:
            records = list(Enzyme.parse(handle))
        self.assertEqual(len(records), 8)
        self.assertEqual(records[0].entry, "1.1.1.1")
        self.assertEqual(records[-1].entry, "2.7.2.1")

    def test_irregular(self):
        with open("KEGG/enzyme.irregular") as handle:
            records = list(Enzyme.parse(handle))
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0].entry, "1.14.18.1")
        self.assertEqual(records[-1].entry, "3.4.21.50")

    def test_new(self):
        with open("KEGG/enzyme.new") as handle:
            records = list(Enzyme.parse(handle))
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].entry, "6.2.1.25")

    def test_4letter(self):
        with open("KEGG/enzyme.4letter") as handle:
            records = list(Enzyme.parse(handle))
            self.assertEqual(len(records), 1)
        self.assertEqual(records[0].entry, "5.4.2.2")
        self.assertEqual(len(records[0].genes), 3776)
        self.assertEqual(records[0].genes[0],
                         ('HSA', ['5236', '55276']))
        self.assertEqual(records[0].genes[8],
                         ('CSAB', ['103224690', '103246223']))


class CompoundTests(unittest.TestCase):
    """Bio.KEGG.Compound tests."""

    def test_sample(self):
        with open("KEGG/compound.sample") as handle:
            records = list(Compound.parse(handle))
        self.assertEqual(len(records), 8)
        self.assertEqual(records[0].entry, "C00023")

    def test_irregular(self):
        with open("KEGG/compound.irregular") as handle:
            records = list(Compound.parse(handle))
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0].entry, "C01454")


class MapTests(unittest.TestCase):
    """Bio.KEGG.Map tests."""

    def test_map00950(self):
        system = System()
        with open("KEGG/map00950.rea") as handle:
            for reaction in Map.parse(handle):
                system.add_reaction(reaction)
        rxs = system.reactions()
        self.assertEqual(len(rxs), 56)
        # sort the reaction output by the string names, so that the
        # output will be consistent between python versions
        rxs.sort(key=lambda x: str(x))
        self.assertEqual(str(rxs[0]),
                         "(R)-N-Methylcoclaurine + (S)-Coclaurine + NADPH + O2 "
                         "<=> 2'-Norberbamunine + 2 H2O + NADP")
        self.assertEqual(str(rxs[-1]),
                         "Tyramine <=> Dopamine")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
