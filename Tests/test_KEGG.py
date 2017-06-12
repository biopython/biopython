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
        self.assertEqual(records[0].name,
                         ['Alcohol dehydrogenase', 'Aldehyde reductase'])
        self.assertEqual(records[0].pathway,
                         [('PATH', 'MAP00010', 'Glycolysis / Gluconeogenesis'),
                          ('PATH', 'MAP00071', 'Fatty acid metabolism'),
                          ('PATH', 'MAP00120', 'Bile acid biosynthesis'),
                          ('PATH', 'MAP00350', 'Tyrosine metabolism'),
                          ('PATH', 'MAP00561', 'Glycerolipid metabolism')])
        self.assertEqual(records[0].structures,
                         [('PDB', ['1A4U', '1A71', '1A72', '1ADB', '1ADC',
                                   '1ADF', '1ADG', '1AGN', '1AXE', '1AXG',
                                   '1B14', '1B15', '1B16', '1B2L', '1BTO',
                                   '1CDO', '1D1S', '1D1T', '1DDA', '1DEH',
                                   '1E3E', '1E3I', '1E3L', '1EE2', '1HDX',
                                   '1HDY', '1HDZ', '1HET', '1HEU', '1HF3',
                                   '1HLD', '1HSO', '1HSZ', '1HT0', '1HTB',
                                   '1LDE', '1LDY', '1QLH', '1QLJ', '1TEH',
                                   '2OHX', '2OXI', '3BTO', '3HUD', '5ADH',
                                   '6ADH', '7ADH'])])
        self.assertEqual(records[0].dblinks,
                         [('IUBMB Enzyme Nomenclature', ['1.1.1.1']),
                          ('ExPASy - ENZYME nomenclature database', ['1.1.1.1']),
                          ('WIT (What Is There) Metabolic Reconstruction', ['1.1.1.1']),
                          ('BRENDA, the Enzyme Database', ['1.1.1.1']),
                          ('SCOP (Structural Classification of Proteins)', ['1.1.1.1'])])
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
        self.assertEqual(records[0].mass, "")  # Why?
        self.assertEqual(records[0].formula, "Fe")
        self.assertEqual(records[0].name,
                         ['Iron', 'Fe2+', 'Fe(II)', 'Fe3+', 'Fe(III)'])
        self.assertEqual(records[0].pathway,
                         [('PATH', 'MAP00860', 'Porphyrin and chlorophyll metabolism')])
        self.assertEqual(records[0].enzyme[0], ('1.1.3.22', 'C'))
        self.assertEqual(records[0].structures, [])
        self.assertEqual(records[0].dblinks[0], ('CAS', ['7439-89-6']))

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
