#!/usr/bin/env python
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

__version__ = "$Revision: 1.6 $"

import doctest, unittest
import random
import sys

if 'requires_wise' in sys.modules:
    del sys.modules['requires_wise']
import requires_wise

from Bio.Wise import psw

class TestPSW(unittest.TestCase):
    def test_Alignment_normal(self):
        a = psw.Alignment()

        a.append(psw.ColumnUnit(0, 98, "SEQUENCE"))
        a.append(psw.ColumnUnit(1, 200, "SEQUENCE"))
        a.append(psw.ColumnUnit(0, 98, "INSERT"))
        a.append(psw.ColumnUnit(1, 201, "SEQUENCE"))
        a.append(psw.ColumnUnit(0, 98, "END"))
        a.append(psw.ColumnUnit(1, 201, "END"))

        self.assertEqual(str(a), "[SEQUENCE(98, 200), INSERT(98, 201), END(98, 201)]")

    def test_Alignment_assertions(self):
        a = psw.Alignment()

        self.assertRaises(AssertionError, a.append, psw.ColumnUnit(1, 200, "SEQUENCE"))
        a.append(psw.ColumnUnit(0, 98, "SEQUENCE"))
        self.assertRaises(AssertionError, a.append, psw.ColumnUnit(0, 200, "SEQUENCE"))
        a.append(psw.ColumnUnit(1, 200, "SEQUENCE"))
        self.assertRaises(AssertionError, a.append, psw.ColumnUnit(1, 200, "SEQUENCE"))

    def test_AlignmentColumn_kinds(self):
        ac = psw.AlignmentColumn(psw.ColumnUnit(0, random.randint(0, 9999), "SEQUENCE"))
        ac.append(psw.ColumnUnit(1, random.randint(0, 9999), "INSERT"))
        self.assertEqual(ac.kind, "INSERT")

        ac = psw.AlignmentColumn(psw.ColumnUnit(0, random.randint(0, 9999), "INSERT"))
        ac.append(psw.ColumnUnit(1, random.randint(0, 9999), "SEQUENCE"))
        self.assertEqual(ac.kind, "INSERT")

        ac = psw.AlignmentColumn(psw.ColumnUnit(0, random.randint(0, 9999), "SEQUENCE"))
        ac.append(psw.ColumnUnit(1, random.randint(0, 9999), "SEQUENCE"))
        self.assertEqual(ac.kind, "SEQUENCE")

        ac = psw.AlignmentColumn(psw.ColumnUnit(0, random.randint(0, 9999), "SEQUENCE"))
        ac.append(psw.ColumnUnit(1, random.randint(0, 9999), "END"))
        self.assertEqual(ac.kind, "END")

    def test_AlignmentColumn_repr(self):
        ac = psw.AlignmentColumn(psw.ColumnUnit(0, 34, "SEQUENCE"))
        ac.append(psw.ColumnUnit(1, 55, "END"))
        self.assertEqual(repr(ac), "END(34, 55)")

    def test_AlignmentColumn_full(self):
        ac = psw.AlignmentColumn(psw.ColumnUnit(0, random.randint(0, 9999), "SEQUENCE"))
        ac.append(psw.ColumnUnit(1, random.randint(0, 9999), "END"))
        self.assertRaises(psw.AlignmentColumnFullException, ac.append, psw.ColumnUnit(1, random.randint(0, 9999), "END"))

    def test_AlignmentColumn_assertions(self):
        self.assertRaises(AssertionError, psw.AlignmentColumn, psw.ColumnUnit(1, random.randint(0, 9999), "SEQUENCE"))

        ac = psw.AlignmentColumn(psw.ColumnUnit(0, random.randint(0, 9999), "SEQUENCE"))
        self.assertRaises(AssertionError, ac.append, psw.ColumnUnit(0, random.randint(0, 9999), "SEQUENCE"))


    def test_ColumnUnit(self):
        self.assertEqual(repr(psw.ColumnUnit(0, 33, "SEQUENCE")),
                         "ColumnUnit(unit=0, column=33, SEQUENCE)")

        self.assertEqual(repr(psw.ColumnUnit(1, 33, "INSERT")),
                         "ColumnUnit(unit=1, column=33, INSERT)")

        self.assertEqual(repr(psw.ColumnUnit(1, 33, "END")),
                         "ColumnUnit(unit=1, column=33, END)")

    PARSED = "[SEQUENCE(39, 22), SEQUENCE(40, 23), SEQUENCE(41, 24), SEQUENCE(42, 25), SEQUENCE(43, 26), SEQUENCE(44, 27), END(0, 27)]"

#    def test_align(self):
#        self.assertEqual(repr(psw.align(("Wise/human_114_g01_exons.fna_01", "Wise/human_114_g02_exons.fna_01"), "introns.bla", 23, 5, quiet=True)), self.PARSED)

def run_tests(argv):
    test_suite = testing_suite()
    runner = unittest.TextTestRunner(sys.stdout, verbosity = 2)
    runner.run(test_suite)

def testing_suite():
    """Generate the suite of tests.
    """
    unittest_suite = unittest.TestSuite()

    test_loader = unittest.TestLoader()
    test_loader.testMethodPrefix = 'test_'
    tests = [TestPSW]
    
    for test in tests:
        cur_suite = test_loader.loadTestsFromTestCase(test)
        unittest_suite.addTest(cur_suite)

    doctest_suite = doctest.DocTestSuite(psw)

    big_suite = unittest.TestSuite((unittest_suite, doctest_suite))

    return big_suite

if __name__ == "__main__":
    unittest_suite = unittest.TestLoader().loadTestsFromName("test_psw")
    doctest_suite = doctest.DocTestSuite(psw)
    suite = unittest.TestSuite((unittest_suite, doctest_suite))
    runner = unittest.TextTestRunner(sys.stdout, verbosity = 2)
    runner.run(suite)
