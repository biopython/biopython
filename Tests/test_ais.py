#!/usr/bin/env python
"""Tests Ais.
"""
# standard library
import sys
import random
import string

# local stuff
import requires_internet

# PyUnit
import unittest
from Bio.Seq import Seq
from Bio.Align.Generic import Alignment
from Bio.Align.AlignInfo import SummaryInfo
from Bio.Alphabet import RNAAlphabet
from Bio.Alphabet import Gapped
from Bio.Ais import match_sequence
from Bio.Ais import Lymphocyte
from Bio.Ais import Immune

def run_tests(argv):
    ALL_TESTS = [ MatchSequenceTest, ImmuneTest ]

    runner = unittest.TextTestRunner(sys.stdout, verbosity = 2)
    test_loader = unittest.TestLoader()
    test_loader.testMethodPrefix = 't_'

    for test in ALL_TESTS:
        cur_suite = test_loader.loadTestsFromTestCase(test)
        runner.run(cur_suite)

# --- helper classes and functions


# --- the actual test classes

class MatchSequenceTest(unittest.TestCase):
    """Some examples of subsequence matching.
    """

    def t_match_sequence(self):
        """Test conversion of hex strings.
        """
        params = [ [ 'patch', 'latch', [ 3, 4, 5 ] ], \
                   [ 'strike', 'strand', [ 2, 3, 4  ] ], \
                   [ 'alter', 'water', [ 2, 3, 4 ] ], \
                   [ 'termite', 'germane', [ 2, 3, 4 ] ]
            ]
        actual = []
        expected = [ 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0 ]
        for j in range( len( params ) ):
            paramlist = params[ j ]
            for threshold in paramlist[ 2 ]:
                assert match_sequence( paramlist[ 0 ], paramlist[ 1 ], threshold ) == expected.pop( 0 ), "Did not match correctly."

class ImmuneTest( unittest.TestCase ):
    """
    Test Immune class.
    """

    def setUp( self ):
        return
        alphabet = 'abcdefghijklmnopqrstuvwxyz'
        test_seqs = { \
            'Bovine': 'ADATLSQITNNIDPVGR..IQMRTRRTLRGHLAKIYAMHWGTDSRLLVSASQDGKLIIWD', \
            'Human' : 'ADATLSQITNNIDPVGR..IQMRTRRTLRGHLAKIYAMHWGTDSRLLVSASQDGKLIIWD', \
            'Mouse' : 'HDVELHQVAERVEALGQ..FVMKTRRTLKGHGNKVLCMDWCKDKRRIVSSSQDGKVIVWD', \
            'Yeast' : 'QDASLFQMANKVTSLTKNKINLKPNIVLKGHNNKISDFRWSRDSKRILSASQDGFMLIWD', \
            'Oreni' :'SASRDKSIIMWKLTRDET.NYGIPQRSLKGHSHFVSDVVISSDG.FALSGAWDGTLRLWD' \
            }

        align = Alignment( Gapped( RNAAlphabet() ) )
        for ( name, seq ) in test_seqs.items():
            seq = seq.strip()
            align.add_sequence( name, seq.lower() )
        # self.immune = Immune( align, alphabet, 5 )
        self.immune = Immune( align, alphabet, do_tuning = 0)


    def t_compute_accum_weight( self ):
        return
        accum_weight = self.immune.compute_accum_weight()
        assert accum_weight == 5, 'Incorrect cumulative weight'

    def t_search_accum_weight( self ):
        return
        for j in range( 0, 5 ):
            location = self.immune.search_accum_weight( j + 1 )
            assert j  == location, 'Search failed'
        location = self.immune.search_accum_weight( 0 )
        assert location == 0, 'Search failed'

if __name__ == "__main__":
    sys.exit(run_tests(sys.argv))
