#!/usr/bin/env python
"""Tests HotRand.
"""
# standard library
import sys
import random

# local stuff
import requires_internet

# PyUnit
import unittest
from Bio.HotRand import HotRandom
from Bio.HotRand import hex_convert

def run_tests(argv):
    ALL_TESTS = [HexConversionTest, RandomSequenceTest ]

    runner = unittest.TextTestRunner(sys.stdout, verbosity = 2)
    test_loader = unittest.TestLoader()
    test_loader.testMethodPrefix = 't_'

    for test in ALL_TESTS:
        cur_suite = test_loader.loadTestsFromTestCase(test)
        runner.run(cur_suite)

# --- helper classes and functions

def are_lists_equal( a, b ):
    if( len( a ) != len( b ) ):
        return 0
    for j in range( 0, len( a ) ):
        if( a[ j ] != b[ j ] ):
            return 0
    return 1


def are_items_in_range( a, high, low ):
    for j in range( 0, len( a ) ):
        if( a[ j ] > high ):
            print 'a[ %d ] is %d' % ( j , a[ j ] )
            return 0
        if( a[ j ] < low ):
            print 'a[ %d ] is %d' % ( j , a[ j ] )
            return 0
    return 1

# --- the actual test classes

class HexConversionTest(unittest.TestCase):
    """Some examples of conversions from hex strings.
    """

    def t_convert(self):
        """Test conversion of hex strings.
        """
        nums = [ '0000', 'abcd', '1234', '5555', '4321', 'aaaa', 'ffff', '7', '21' ]
        actual = []
        expected = [ 0, 43981, 4660, 21845, 17185, 43690, 65535, 7, 33 ]
        for num in nums:
            actual.append( hex_convert( num ) )
        assert are_lists_equal( actual, expected ), "Did not convert string."

class RandomSequenceTest(unittest.TestCase):
    """Test sequence of random numbers.
    """

    def t_get_random_range(self):
        """Get a sequence of random numbers.
        """
        return
        rand_seq = []
        hot_random = HotRandom()
        for j in range( 0, 200 ):
            rand_num = hot_random.hot_rand( 91, 37 )
            rand_seq.append( rand_num )
        assert are_items_in_range( rand_seq, 91, 37 ) , "Got an out of range number"
        rand_seq = []
        for j in range( 0, 200 ):
            rand_num = hot_random.hot_rand( 19, 0 )
            rand_seq.append( rand_num )
        assert are_items_in_range( rand_seq, 19, 0 ) , "Got an out of range number"
        rand_seq = []
        for j in range( 0, 200 ):
            rand_num = hot_random.hot_rand( 61, 4 )
            rand_seq.append( rand_num )
        assert are_items_in_range( rand_seq, 61, 4 ) , "Got an out of range number"


if __name__ == "__main__":
    sys.exit(run_tests(sys.argv))
