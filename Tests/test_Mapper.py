# Copyright 2012 Lenna X. Peterson <arklenna@gmail.com>

# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for Mapper module."""

from functools import wraps
import unittest
import warnings

from Bio import BiopythonExperimentalWarning

with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonExperimentalWarning)
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqUtils.Mapper import CoordinateMapper, CDSPosition
    from Bio.SeqUtils.Mapper.MapPositions import (
        CDSPositionError,
        GenomePositionError,
        ProteinPositionError,
    )


# Test data

exons = [(5808, 5860), (6757, 6874), (7767, 7912), (13709, 13785)]
cmap = CoordinateMapper(exons)
g_exons = range(7868, 7875)  # len 7
c_exons = range(270, 279)  # len 9
p_exons = range(90, 93)  # len 3
p_exons_trip = [p for p in p_exons for n in range(3)]  # len 9
c_exon_prs = ((270, 272), (273, 275), (276, 278))  # len 3
g_exon_prs = ((7868, 7870), (7871, 7873), (7874, 7876))  # len 3

g_introns = {
    None: (5860, 5861, 6308, 6309, 6755, 6756),
    "HGVS": (5861, 5862, 6309, 6310, 6756, 6757),
}
c_introns = {
    None: ("51+1", "51+2", "51+449", "52-448", "52-2", "52-1"),
    "HGVS": ("52+1", "52+2", "52+449", "53-448", "53-2", "53-1"),
}
c_intron_tups = ((51, 1), (51, 2), (51, 449), (52, -448), (52, -2), (52, -1))

g_outside = {None: (0, 13785, 14000), "HGVS": (1, 13786, 14001)}
c_outside = {None: ("-5808", "+1", "+216"), "HGVS": ("-5808", "*1", "*216")}
c_outside_tups = ((None, -5808), (None, 1), (None, 216))


# Test decorator


def two_dialects(fn):
    @wraps(fn)
    def call_tests(self):
        orig_dialect = self.dialect
        for dialect in (None, "HGVS"):
            self.dialect = dialect
            fn(self)
        self.dialect = orig_dialect

    return call_tests


# Test suites


class TestCDSPosition(unittest.TestCase):
    """Test that CDSPosition works properly."""

    def setUp(self):
        self.dialect = None

    @two_dialects
    def testGoodIntron(self):
        """A CDSPosition should match good intron values."""
        for c_args, c in zip(c_intron_tups, c_introns[self.dialect]):
            actual = CDSPosition.from_anchor(*c_args).to(self.dialect)
            self.assertEqual(actual, c)

    @two_dialects
    def testGoodOutside(self):
        """A CDSPosition should match good outside-CDS values."""
        for c_args, c in zip(c_outside_tups, c_outside[self.dialect]):
            actual = CDSPosition.from_anchor(*c_args).to(self.dialect)
            self.assertEqual(actual, c)

    def testEqual(self):
        """Two CDSPositions should test equal with same args."""
        for args in c_intron_tups:
            CPos = CDSPosition.from_anchor
            self.assertEqual(CPos(*args), CPos(*args))
        self.assertEqual(str(CPos(*args)), str(CPos(*args)))

    @two_dialects
    def testEqualDialects(self):
        """A CDSPosition should be round-trippable."""
        for c_pos in c_outside[self.dialect]:
            c_actual = CDSPosition.from_dialect(self.dialect, c_pos)
            self.assertEqual(c_pos, c_actual.to(self.dialect))

    def testBadAnchor(self):
        """A CDSPosition should fail with negative CDS anchor."""
        self.assertRaises(CDSPositionError, CDSPosition.from_anchor, -1, 1)

    def testBadOffset(self):
        """A CDSPosition should fail with zero offset."""
        self.assertRaises(CDSPositionError, CDSPosition.from_anchor, None, 0)

    def testValidate(self):
        """Setting CDSPosition anchor offset should fail for invalid cases."""
        intron = CDSPosition("22+4")
        self.assertRaises(CDSPositionError, setattr, intron, "offset", 0)
        self.assertRaises(CDSPositionError, setattr, intron, "anchor", -5)
        exon = CDSPosition("60")
        self.assertRaises(CDSPositionError, setattr, exon, "anchor", None)


class TestCoordinateMapper(unittest.TestCase):
    """Test that CoordinateMapper works properly."""

    def setUp(self):
        self.sf = SeqFeature(sum((FeatureLocation(s, e) for s, e in exons)), type="CDS")

    def testListInit(self):
        # CoordinateMapper should take list of exon pairs
        cm = CoordinateMapper(exons)
        # FIXME CompoundLocation does not have __eq__
        self.assertEqual(str(cm.exons), str(self.sf.location))

    def testSeqRecordInit(self):
        # CoordinateMapper should use first CDS feature of SeqRecord
        sr = SeqRecord("", features=[self.sf])
        cm = CoordinateMapper(sr)
        self.assertEqual(str(cm.exons), str(self.sf.location))

    def testSeqFeatureInit(self):
        # CoordinateMapper should take a CDS SeqFeature
        cm = CoordinateMapper(self.sf)
        self.assertEqual(str(cm.exons), str(self.sf.location))

    def testEmptyInit(self):
        # CoordinateMapper should fail with no arguments
        self.assertRaises(Exception, CoordinateMapper)


class Testg2c(unittest.TestCase):
    """Test success of good input and failure of bad input to g2c."""

    def setUp(self):
        self.dialect = None

    # any integer should return
    def testGoodExons(self):
        # g2c should work for exon positions
        for arg, expected in zip(g_exons, c_exons):
            self.assertEqual(cmap.g2c(arg), expected)

    @two_dialects
    def testGoodIntronsStr(self):
        # g2c should work for intron positions (string)
        dia = self.dialect
        for arg, expect in zip(g_introns[dia], c_introns[dia]):
            actual = cmap.g2c(arg, dia).to(dia)
            self.assertEqual(actual, expect)

    @two_dialects
    def testGoodIntrons(self):
        # g2c should work for intron positions (CDSPosition)
        for arg, tup in zip(g_introns[self.dialect], c_intron_tups):
            actual = cmap.g2c(arg, self.dialect)
            expect = CDSPosition.from_anchor(*tup)
            self.assertEqual(actual, expect)
            self.assertEqual(str(actual), str(expect))

    @two_dialects
    def testGoodOutsideStr(self):
        # g2c should work for outside positions (string)
        dia = self.dialect
        for arg, expected in zip(g_outside[dia], c_outside[dia]):
            actual = cmap.g2c(arg, dia).to(dia)
            self.assertEqual(actual, expected)

    def testGoodOutside(self):
        # g2c should work for outside positions (CDSPosition)
        for arg, tup in zip(g_outside[self.dialect], c_outside_tups):
            actual = cmap.g2c(arg)
            expect = CDSPosition.from_anchor(*tup)
            self.assertEqual(actual, expect)
            self.assertEqual(str(actual), str(expect))

    # TODO should it handle non-exact positions?

    @two_dialects
    def testZeroArg(self):
        # g2c should work for 0 in no dialect and fail for 1-index
        args = (0, self.dialect)
        if self.dialect is None:
            cmap.g2c(*args)
        else:
            self.assertRaises(GenomePositionError, cmap.g2c, *args)

    @two_dialects
    def testBadArg(self):
        # g2c should fail on string, float, None, or negative
        bad = ("string", None, 3.14, -5)
        for arg in bad:
            self.assertRaises(Exception, cmap.g2c, arg)


class Testc2p(unittest.TestCase):
    """Test success of good input and failure of bad input to c2p."""

    # integer within length of exon should return correct protein
    def testGoodExons(self):
        # c2p should work for exon positions
        for arg, expect in zip(c_exons, p_exons_trip):
            self.assertEqual(cmap.c2p(arg), expect)

    # FIXME should CDSPosition return None or raise error?
    def testBadNotProtein(self):
        # c2p should fail for CDSPosition or str
        bad = ("string", CDSPosition.from_anchor(7, -5))
        for arg in bad:
            self.assertRaises(ValueError, cmap.c2p, arg)


class Testp2c(unittest.TestCase):
    """Test success of good input and failure of bad input to p2c."""

    # integer within length of protein should return correct range
    def testGoodExons(self):
        # p2c should work for exon positions
        for arg, expect in zip(p_exons, c_exon_prs):
            self.assertEqual(cmap.p2c(arg), expect)

    def testBadTooLarge(self):
        # p2c should fail for positions longer than the max length
        self.assertRaises(ProteinPositionError, cmap.p2c, 500)

    def testBadNegative(self):
        # p2c should fail for negative protein coordinate
        self.assertRaises(ProteinPositionError, cmap.p2c, -50)


class Testc2g(unittest.TestCase):
    """Test success of good input and failure of bad input to c2g."""

    def setUp(self):
        self.dialect = None

    # integer within length of exon should return correct value
    def testGoodExons(self):
        # c2g should work for exon positions
        for arg, expected in zip(c_exons, g_exons):
            self.assertEqual(cmap.c2g(arg), expected)

    # CDSPosition object should return correct value
    @two_dialects
    def testGoodIntrons(self):
        # c2g should work for introns (CDSPosition)
        for c_args, g in zip(c_intron_tups, g_introns[self.dialect]):
            actual = cmap.c2g(CDSPosition.from_anchor(*c_args))
            self.assertEqual(actual.to(self.dialect), g)

    @two_dialects
    def testGoodOutside(self):
        # c2g should work for outside (CDSPosition)
        for c_args, expect in zip(c_outside_tups, g_outside[self.dialect]):
            actual = cmap.c2g(CDSPosition.from_anchor(*c_args))
            self.assertEqual(actual.to(self.dialect), expect)

    # simple string should return correct value
    #    \d*[+-]\d+
    @two_dialects
    def testGoodIntronsStr(self):
        # c2g should work for intron positions (string)
        dia = self.dialect
        # c_str from dialect to *dialect* matches *dialect* g_str
        for c_str, g in zip(c_introns[dia], g_introns[dia]):
            self.assertEqual(cmap.c2g(c_str, dia).to(dia), g)

        # c_str from dialect as *default* matches *default* g_str
        for c_str, g in zip(c_introns[dia], g_introns[None]):
            self.assertEqual(cmap.c2g(c_str, dia), g)

    @two_dialects
    def testGoodOutsideStr(self):
        # c2g should work for outside positions (string)
        dia = self.dialect
        for expected, arg in zip(g_outside[dia], c_outside[dia]):
            self.assertEqual(cmap.c2g(arg, dia).to(dia), expected)

    def testNegativeInt(self):
        # c2g should treat negative integer as pre-CDS
        self.assertEqual(cmap.c2g(-2), 5806)

    # bad positions in form 123+4 where 123 is not end of an exon
    # or where sign is wrong
    @two_dialects
    def testBadNotEnd(self):
        # c2g should fail on invalid introns
        bad_introns = {
            None: (
                "160+5",  # 160 is exonic
                "51-6",  # 51 is end
                "52+10",  # 52 is start
            ),
            "HGVS": ("160+5", "52-6", "53+10"),
        }
        for arg in bad_introns[self.dialect]:
            self.assertRaises(CDSPositionError, cmap.c2g, arg, self.dialect)

    @two_dialects
    def testBadString(self):
        # c2g should fail on arbitrary string
        bad_vals = {None: ["xxx", "4-5'", "*10"], "HGVS": ["xxx", "4-5"]}
        for arg in bad_vals[self.dialect]:
            self.assertRaises(ValueError, cmap.c2g, arg, self.dialect)

    # integers outside length of exon should raise IndexError
    def testBadTooLarge(self):
        # c2g should fail on values beyond the exon
        num = 400
        assert num > len(cmap.exons)
        self.assertRaises(IndexError, cmap.c2g, num)
        self.assertRaises(Exception, cmap.c2g, "%d%d" % (num, -12))


class Testp2g(unittest.TestCase):
    def testGoodStr(self):
        # p2g should return expected ranges
        for p_arg, gpos in zip(p_exons, g_exon_prs):
            self.assertEqual(cmap.p2g(p_arg), gpos)


class Testg2p(unittest.TestCase):
    def testGoodStr(self):
        # g2p should return correct position for exons
        for g_arg, ppos in zip(g_exons, p_exons_trip):
            self.assertEqual(cmap.g2p(g_arg), ppos)


class TestStrand(unittest.TestCase):
    def setUp(self):
        locations = FeatureLocation(24, 30, -1) + FeatureLocation(35, 40, +1)
        feature = SeqFeature(locations, type="CDS")
        self.cm = CoordinateMapper(feature)
        self.g_exons = [29, 28, 27, 26, 25, 24, 35, 36, 37, 38, 39]
        assert list(self.cm.exons) == self.g_exons
        c_len = 11
        assert c_len == len(list(self.cm.exons))
        self.c_exons = range(c_len)

    def testg2c(self):
        # c2g and g2c should work for mixed-strand exons
        for c, g in zip(self.c_exons, self.g_exons):
            self.assertEqual(c, self.cm.g2c(g))
            self.assertEqual(g, self.cm.c2g(c))
            if c < 6:
                self.assertEqual(self.cm.c2g(c).strand, -1)
            else:
                self.assertEqual(self.cm.c2g(c).strand, +1)


if __name__ == "__main__":
    suite0 = unittest.TestLoader().loadTestsFromTestCase(TestCDSPosition)
    suite1 = unittest.TestLoader().loadTestsFromTestCase(TestCoordinateMapper)
    suite2 = unittest.TestLoader().loadTestsFromTestCase(Testc2g)
    suite3 = unittest.TestLoader().loadTestsFromTestCase(Testg2c)
    suite4 = unittest.TestLoader().loadTestsFromTestCase(Testc2p)
    suite5 = unittest.TestLoader().loadTestsFromTestCase(Testp2c)
    suite6 = unittest.TestLoader().loadTestsFromTestCase(Testp2g)
    suite7 = unittest.TestLoader().loadTestsFromTestCase(Testg2p)
    suite8 = unittest.TestLoader().loadTestsFromTestCase(TestStrand)
    all_suites = [
        suite0,
        suite1,
        suite2,
        suite3,
        suite4,
        suite5,
        suite6,
        suite7,
        suite8,
    ]
    all = unittest.TestSuite(all_suites)
    unittest.TextTestRunner(verbosity=3).run(all)
