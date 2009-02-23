#!/usr/bin/env python
"""Tests HotRand.
"""
# standard library
import sys
import unittest

# local stuff
import requires_internet
requires_internet.check()

from Bio.HotRand import HotRandom


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

    def test_get_random_range(self):
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
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
