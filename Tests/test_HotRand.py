#!/usr/bin/env python
"""Tests HotRand.
"""
# standard library
import sys

# local stuff
if sys.modules.has_key('requires_internet'):
    del sys.modules['requires_internet']
import requires_internet

# PyUnit
import unittest
from Bio.HotRand import HotRandom

def run_tests(argv):
    runner = unittest.TextTestRunner(sys.stdout, verbosity = 2)
    test_loader = unittest.TestLoader()
    test_loader.testMethodPrefix = 't_'
    cur_suite = test_loader.loadTestsFromTestCase(RandomSequenceTest)
    runner.run(cur_suite)

# --- helper classes and functions

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
