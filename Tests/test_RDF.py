# Copyright 2014 by Joachim Baran
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Unittests for RDFization."""
from __future__ import print_function

import unittest
import sys
if sys.version_info[0] == 3:
    maketrans = str.maketrans
else:
    from string import maketrans

from Bio.SeqFeature import FeatureLocation, ExactPosition, WithinPosition, \
                           BetweenPosition, BeforePosition, AfterPosition, \
                           OneOfPosition

class RDFizationMethodTests(unittest.TestCase):

    def test_feature_location(self):
        """Check RDFization of a FeatureLocation object."""
        rdf = FeatureLocation(5, 10)._rdfize('http://unit-test/')
        self.assertEqual(rdf,
                         """<http://unit-test/Location5,10,0> a <http://biohackathon.org/resource/faldo#Region> ;
    <http://biohackathon.org/resource/faldo#begin> <http://unit-test/Location5,10,0/5> ;
    <http://biohackathon.org/resource/faldo#end> <http://unit-test/Location5,10,0/10> .

<http://unit-test/Location5,10,0/5> a <http://biohackathon.org/resource/faldo#ExactPosition> ;
    <http://biohackathon.org/resource/faldo#position> 6 .
<http://unit-test/Location5,10,0/10> a <http://biohackathon.org/resource/faldo#ExactPosition> ;
    <http://biohackathon.org/resource/faldo#position> 10 .""")

        rdf = FeatureLocation(5, 10, strand = -1)._rdfize('http://unit-test/')
        self.assertEqual(rdf,
                         """<http://unit-test/Location5,10,-1/5> a <http://biohackathon.org/resource/faldo#ReverseStrandPosition> .
<http://unit-test/Location5,10,-1/10> a <http://biohackathon.org/resource/faldo#ReverseStrandPosition> .
<http://unit-test/Location5,10,-1> a <http://biohackathon.org/resource/faldo#Region> ;
    <http://biohackathon.org/resource/faldo#begin> <http://unit-test/Location5,10,-1/10> ;
    <http://biohackathon.org/resource/faldo#end> <http://unit-test/Location5,10,-1/5> .

<http://unit-test/Location5,10,-1/10> a <http://biohackathon.org/resource/faldo#ExactPosition> ;
    <http://biohackathon.org/resource/faldo#position> 10 .
<http://unit-test/Location5,10,-1/5> a <http://biohackathon.org/resource/faldo#ExactPosition> ;
    <http://biohackathon.org/resource/faldo#position> 6 .""")

    def test_exact_position(self):
        """Check RDFization of a ExactPosition object."""
        rdf = ExactPosition(12)._rdfize('http://unit-test/')
        self.assertEqual(rdf,
                         """<http://unit-test/12> a <http://biohackathon.org/resource/faldo#ExactPosition> ;
    <http://biohackathon.org/resource/faldo#position> 12 .""")

    def test_within_position(self):
        """Check RDFization of a WithinPosition object."""
        rdf = WithinPosition(10, 10, 20)._rdfize('http://unit-test/')
        self.assertEqual(rdf,
                         """<http://unit-test/10-20> a <http://biohackathon.org/resource/faldo#InRangePosition> ;
    <http://biohackathon.org/resource/faldo#begin> <http://unit-test/10-20/begin> ;
    <http://biohackathon.org/resource/faldo#end> <http://unit-test/10-20/end> .

<http://unit-test/10-20/begin> a <http://biohackathon.org/resource/faldo#ExactPosition> ;
    <http://biohackathon.org/resource/faldo#position> 10 .

<http://unit-test/10-20/end> a <http://biohackathon.org/resource/faldo#ExactPosition> ;
    <http://biohackathon.org/resource/faldo#position> 20 .""")

    def test_between_position(self):
        """Check RDFization of a BetweenPosition object."""
        rdf = BetweenPosition(10, 10, 20)._rdfize('http://unit-test/')
        self.assertEqual(rdf,
                         """<http://unit-test/10-20> a <http://biohackathon.org/resource/faldo#InBetweenPosition> ;
    <http://biohackathon.org/resource/faldo#after> <http://unit-test/10-20/begin> ;
    <http://biohackathon.org/resource/faldo#before> <http://unit-test/10-20/end> .

<http://unit-test/10-20/begin> a <http://biohackathon.org/resource/faldo#ExactPosition> ;
    <http://biohackathon.org/resource/faldo#position> 10 .

<http://unit-test/10-20/end> a <http://biohackathon.org/resource/faldo#ExactPosition> ;
    <http://biohackathon.org/resource/faldo#position> 20 .""")

    def test_before_position(self):
        """Check RDFization of a BeforePosition object."""
        rdf = BeforePosition(5)._rdfize('http://unit-test/')
        self.assertEqual(rdf,
                         """<http://unit-test/<5> a <http://biohackathon.org/resource/faldo#InRangePosition> ;
    <http://biohackathon.org/resource/faldo#end> <http://unit-test/<5/end> .

<http://unit-test/<5/end> a <http://biohackathon.org/resource/faldo#ExactPosition> ;
    <http://biohackathon.org/resource/faldo#position> 5 .""")

    def test_after_position(self):
        """Check RDFization of a AfterPosition object."""
        rdf = AfterPosition(5)._rdfize('http://unit-test/')
        self.assertEqual(rdf,
                         """<http://unit-test/>5> a <http://biohackathon.org/resource/faldo#InRangePosition> ;
    <http://biohackathon.org/resource/faldo#begin> <http://unit-test/>5/begin> .

<http://unit-test/>5/begin> a <http://biohackathon.org/resource/faldo#ExactPosition> ;
    <http://biohackathon.org/resource/faldo#position> 5 .""")

    def test_one_of_position(self):
        """Check RDFization of an OneOfPosition object."""
        rdf = OneOfPosition(10,
                            [
                                ExactPosition(9),
                                ExactPosition(10),
                                ExactPosition(11)
                            ])._rdfize('http://unit-test/')
        self.assertEqual(rdf,
                         """<http://unit-test/9,10,11> a <http://biohackathon.org/resource/faldo#OneOfPosition> .
<http://unit-test/9,10,11/0> a <http://biohackathon.org/resource/faldo#ExactPosition> ;
    <http://biohackathon.org/resource/faldo#position> 9 .
<http://unit-test/9,10,11> <http://biohackathon.org/resource/faldo#position> <http://unit-test/9,10,11/0> .
<http://unit-test/9,10,11/1> a <http://biohackathon.org/resource/faldo#ExactPosition> ;
    <http://biohackathon.org/resource/faldo#position> 10 .
<http://unit-test/9,10,11> <http://biohackathon.org/resource/faldo#position> <http://unit-test/9,10,11/1> .
<http://unit-test/9,10,11/2> a <http://biohackathon.org/resource/faldo#ExactPosition> ;
    <http://biohackathon.org/resource/faldo#position> 11 .
<http://unit-test/9,10,11> <http://biohackathon.org/resource/faldo#position> <http://unit-test/9,10,11/2> .""")

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
