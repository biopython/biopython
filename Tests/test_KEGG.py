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
                         ['alcohol dehydrogenase',
                          'aldehyde reductase',
                          'ADH',
                          'alcohol dehydrogenase (NAD)',
                          'aliphatic alcohol dehydrogenase',
                          'ethanol dehydrogenase',
                          'NAD-dependent alcohol dehydrogenase',
                          'NAD-specific aromatic alcohol dehydrogenase',
                          'NADH-alcohol dehydrogenase',
                          'NADH-aldehyde dehydrogenase',
                          'primary alcohol dehydrogenase',
                          'yeast alcohol dehydrogenase'])
        self.assertEqual(records[0].pathway,
                         [('PATH', 'ec00010', 'Glycolysis / Gluconeogenesis'),
                          ('PATH', 'ec00071', 'Fatty acid degradation'),
                          ('PATH', 'ec00260', 'Glycine, serine and threonine metabolism'),
                          ('PATH', 'ec00350', 'Tyrosine metabolism'),
                          ('PATH', 'ec00592', 'alpha-Linolenic acid metabolism'),
                          ('PATH', 'ec00625', 'Chloroalkane and chloroalkene degradation'),
                          ('PATH', 'ec00626', 'Naphthalene degradation'),
                          ('PATH', 'ec00830', 'Retinol metabolism'),
                          ('PATH', 'ec00980', 'Metabolism of xenobiotics by cytochrome P450'),
                          ('PATH', 'ec00982', 'Drug metabolism - cytochrome P450'),
                          ('PATH', 'ec01100', 'Metabolic pathways'),
                          ('PATH', 'ec01110', 'Biosynthesis of secondary metabolites'),
                          ('PATH', 'ec01120', 'Microbial metabolism in diverse environments'),
                          ('PATH', 'ec01130', 'Biosynthesis of antibiotics')])
        self.assertEqual(records[0].dblinks,
                         [('ExplorEnz - The Enzyme Database', ['1.1.1.1']),
                          ('IUBMB Enzyme Nomenclature', ['1.1.1.1']),
                          ('ExPASy - ENZYME nomenclature database', ['1.1.1.1']),
                          ('UM-BBD (Biocatalysis/Biodegradation Database)', ['1.1.1.1']),
                          ('BRENDA, the Enzyme Database', ['1.1.1.1']),
                          ('CAS', ['9031-72-5'])])
        self.assertEqual(records[-1].entry, "2.7.2.1")
        self.assertEqual(str(records[-1]).replace(" ", "").split("\n")[:10],
                         ['ENTRYEC2.7.2.1', 'NAMEacetatekinase', 'acetokinase',
                          'AckA', 'AK', 'acetickinase', 'acetatekinase(phosphorylating)',
                          'CLASSTransferases;', 'Transferringphosphorus-containinggroups;',
                          'Phosphotransferaseswithacarboxygroupasacceptor'])

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
        self.assertEqual(records[0].genes[0],
                         ('HSA', ['5236', '55276']))
        self.assertEqual(records[0].genes[8],
                         ('CSAB', ['103224690', '103246223']))

    def test_exceptions(self):
        with open("KEGG/enzyme.sample") as handle:
            with self.assertRaises(ValueError) as context:
                list(Enzyme.read(handle))
            self.assertTrue("More than one record found in handle" in str(context.exception))
            records = Enzyme.parse(handle)
            for i in range(0, 6):
                next(records)
            self.assertRaises(StopIteration, next, records)


class CompoundTests(unittest.TestCase):
    """Bio.KEGG.Compound tests."""

    def test_sample(self):
        with open("KEGG/compound.sample") as handle:
            records = list(Compound.parse(handle))
        self.assertEqual(len(records), 8)
        self.assertEqual(records[1].entry, "C00017")
        self.assertEqual(records[1].mass, "")  # Why?
        self.assertEqual(records[1].formula, "C2H4NO2R(C2H2NOR)n")
        self.assertEqual(records[1].name,
                         ['Protein'])
        self.assertEqual(records[1].pathway,
                         [('PATH', 'map00450', 'Selenocompound metabolism')])
        self.assertEqual(len(records[1].enzyme), 21)
        self.assertEqual(records[1].enzyme[0], ('2.3.2.6'))
        self.assertEqual(records[1].structures, [])
        self.assertEqual(records[1].dblinks[0], ('PubChem', ['3319']))
        self.assertEqual(str(records[-1]).replace(" ", "").split("\n")[:10],
                         ['ENTRYC01386', 'NAMENH2Mec', '7-Amino-4-methylcoumarin',
                          'FORMULAC10H9NO2', 'DBLINKSCAS:26093-31-2',
                          'PubChem:4580', 'ChEBI:51771', 'ChEMBL:CHEMBL270672',
                          'KNApSAcK:C00048593', 'PDB-CCD:MCM'])

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
