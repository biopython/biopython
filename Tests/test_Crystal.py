# Copyright 2002 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# python unittest framework
import unittest
import copy
import sys

# modules to be tested
from Bio.Crystal import Hetero
from Bio.Crystal import Chain
from Bio.Crystal import Crystal
from Bio.Crystal import CrystalError


class ChainTestCase(unittest.TestCase):

    def setUp( self ):
        self.a = 'C A A C T A G G T C A C U A G G T C A G'
        self.b = 'C T G A C C T A G T G A C C T A G T T G'
        self.c = 'THR LYS LEU ASN GLY MET VAL LEU LEU CYS LYS VAL CYS GLY ASP'
        self.d = 'THR LYS LEU ASN GLY MET VAL LEU LEU CYS LYS VAL CYS GLY ASP '
        self.e = 'TYR LYS LEU ASN GLY MET VAL LEU LEU CYS LYS VAL CYS GLY ASP '
        self.f = 'THR LYS LEU ASN GLY MET VAL LEU LEU CYS LYS VAL CYS GLY SER '
        self.g = 'C A A C T A G G T C A C U A G G T C A T'
        self.h = 'G A A C T A G G T C A C U A G G T C A G'

    def testEquals( self ):
        first = Chain( self.a )
        second = Chain( self.a )
        self.assertEquals( first, second )

        first = Chain( self.b )
        second = Chain( self.b )
        self.assertEquals( first, second )

        first = Chain( self.c )
        second = Chain( self.c )
        self.assertEquals( first, second )

        first = Chain( self.a )
        second = Chain( self.g )
        self.assertNotEquals( first, second )

        first = Chain( self.a )
        second = Chain( self.h )
        self.assertNotEquals( first, second )

        first = Chain( self.c )
        second = Chain( self.e )
        self.assertNotEquals( first, second )

        first = Chain( self.c )
        second = Chain( self.f )
        self.assertNotEquals( first, second )


    def testLen( self ):
        chain = Chain( self.a )
        elements = self.a.strip().split()
        num_elements = len( elements )
        self.assertEquals( len( chain ), num_elements )

        chain = Chain( self.b )
        elements = self.b.strip().split()
        num_elements = len( elements )
        self.assertEquals( len( chain ), num_elements )

        chain = Chain( self.c )
        elements = self.c.strip().split()
        num_elements = len( elements )
        self.assertEquals( len( chain ), num_elements )


    def testAppend( self ):
        chain = Chain( self.a[:] )
        chain.append( 'U' )
        elements = self.a.strip().split()
        num_elements = len( elements )
        last_element = chain.data[ -1 ]
        self.assertEquals( 'u', last_element.data )
        self.assertEquals( len( chain ), num_elements + 1 )

        chain = Chain( self.a[:] )
        chain.append( Hetero( 'A' ) )
        elements = self.a.strip().split()
        num_elements = len( elements )
        last_element = chain.data[ -1 ]
        self.assertEquals( 'a', last_element.data )
        self.assertEquals( len( chain ), num_elements + 1 )

        chain = Chain( self.b[:] )
        chain.append( 't' )
        elements = self.b.strip().split()
        num_elements = len( elements )
        last_element = chain.data[ -1 ]
        self.assertEquals( 't', last_element.data )
        self.assertEquals( len( chain ), num_elements + 1 )

        chain = Chain( self.b[:] )
        chain.append( Hetero( 'C' ) )
        elements = self.b.strip().split()
        num_elements = len( elements )
        last_element = chain.data[ -1 ]
        self.assertEquals( 'c', last_element.data )
        self.assertEquals( len( chain ), num_elements + 1 )

        chain = Chain( self.c[:] )
        chain.append( 'ser' )
        elements = self.c.strip().split()
        num_elements = len( elements )
        last_element = chain.data[ -1 ]
        self.assertEquals( 'ser', last_element.data )
        self.assertEquals( len( chain ), num_elements + 1 )


    def testInsert( self ):
        chain = Chain( self.a[:] )
        i = 4
        chain.insert( i, 'g' )
        elements = self.a.strip().split()
        num_elements = len( elements )
        target_element = chain.data[ i ]
        self.assertEquals( 'g', target_element.data )
        self.assertEquals( len( chain ), num_elements + 1 )

        chain = Chain( self.a[:] )
        i = 0
        chain.insert( i, 't' )
        elements = self.a.strip().split()
        num_elements = len( elements )
        target_element = chain.data[ i ]
        self.assertEquals( 't', target_element.data )
        self.assertEquals( len( chain ), num_elements + 1 )

        chain = Chain( self.b[:] )
        i = 9
        chain.insert( i, Hetero( 'a' ) )
        elements = self.a.strip().split()
        num_elements = len( elements )
        target_element = chain.data[ i ]
        self.assertEquals( 'a', target_element.data )
        self.assertEquals( len( chain ), num_elements + 1 )

        chain = Chain( self.c[:] )
        i = 5
        chain.insert( i, 'gln' )
        elements = self.c.strip().split()
        num_elements = len( elements )
        target_element = chain.data[ i ]
        self.assertEquals( 'gln', target_element.data )
        self.assertEquals( len( chain ), num_elements + 1 )


    def testRemove( self ):

        chain = Chain( self.a[:] )
        elements = self.a.strip().split()
        num_elements = len( elements )
        num_a = chain.data.count( Hetero( 'a' ) )
        chain.remove( 'a' )
        num_a_remaining = chain.data.count( Hetero( 'a' ) )
        self.assertEquals( num_a_remaining, num_a - 1 )
        self.assertEquals( len( chain ), num_elements - 1 )

        chain = Chain( self.b[:] )
        elements = self.b.strip().split()
        num_elements = len( elements )
        num_b = chain.data.count( Hetero( 't' ) )
        chain.remove( 't' )
        num_b_remaining = chain.data.count( Hetero( 't' ) )
        self.assertEquals( num_b_remaining, num_b - 1 )
        self.assertEquals( len( chain ), num_elements - 1 )

        chain = Chain( self.c[:] )
        elements = self.c.strip().split()
        num_elements = len( elements )
        num_leu = chain.data.count( Hetero( 'leu' ) )
        chain.remove( 'leu' )
        num_leu_remaining = chain.data.count( Hetero( 'leu' ) )
        self.assertEquals( num_leu_remaining, num_leu - 1 )
        self.assertEquals( len( chain ), num_elements - 1 )


    def testCount( self ):
        chain = Chain( self.a[:] )
        num_a = chain.data.count( Hetero( 'a' ) )
        self.assertEquals( chain.count( 'a' ), num_a )

        chain = Chain( self.b[:] )
        num_a = chain.data.count( Hetero( 't' ) )
        self.assertEquals( chain.count( 't' ), num_a )

        chain = Chain( self.c[:] )
        num_a = chain.data.count( Hetero( 'leu' ) )
        self.assertEquals( chain.count( 'leu' ), num_a )

        chain = Chain( self.c[:] )
        num_a = chain.data.count( Hetero( 'cys' ) )
        self.assertEquals( chain.count( 'cys' ), num_a )


    def testIndex( self ):
        chain = Chain( self.a[:] )
        index_g = chain.data.index( Hetero( 'g' ) )
        self.assertEquals( chain.index( 'g' ), index_g )

        chain = Chain( self.b[:] )
        index_c = chain.data.index( Hetero( 'c' ) )
        self.assertEquals( chain.index( 'c' ), index_c )

        chain = Chain( self.c[:] )
        index_met = chain.data.index( Hetero( 'met' ) )
        self.assertEquals( chain.index( 'met' ), index_met )

    def testGetItem( self ):
        chain = Chain( self.a[:] )
        element_3 = chain.data[ 3 ]
        self.assertEquals( chain[ 3 ], element_3 )

        chain = Chain( self.a[:] )
        element_0 = chain.data[ 0 ]
        self.assertEquals( chain[ 0 ], element_0 )

        chain = Chain( self.b[:] )
        element_7 = chain.data[ 7 ]
        self.assertEquals( chain[ 7 ], element_7 )

        chain = Chain( self.b[:] )
        last_element = chain.data[ -1 ]
        self.assertEquals( chain[ -1 ], last_element )

        chain = Chain( self.c[:] )
        element_8 = chain.data[ 8 ]
        self.assertEquals( chain[ 8 ], element_8 )

    def testSetItem( self ):
        chain = Chain( self.a[:] )
        chain[ 2 ] = 't'
        element_2 = chain.data[ 2 ]
        self.assertEquals( chain[ 2 ], element_2 )

        chain = Chain( self.a[:] )
        chain[ 0 ] = Hetero( 'U' )
        element_0 = chain.data[ 0 ]
        self.assertEquals( chain[ 0 ], element_0 )

        chain = Chain( self.b[:] )
        chain[ -1 ] = Hetero( 'c' )
        last_element = chain.data[ -1 ]
        self.assertEquals( chain[ -1 ], last_element )

        chain = Chain( self.b[:] )
        chain[ 1 ] = 'a'
        element_1 = chain.data[ 1 ]
        self.assertEquals( chain[ 1 ], element_1 )

        chain = Chain( self.c[:] )
        chain[ 5 ] = 'ser'
        element_5 = chain.data[ 5 ]
        self.assertEquals( chain[ 5 ], element_5 )

    def testDelItem( self ):

        chain = Chain( self.a[:] )
        elements = self.a.strip().split()
        num_elements = len( elements )
        num_t = chain.data.count( Hetero( 't' ) )
        del chain[ 4 ]
        num_t_remaining = chain.data.count( Hetero( 't' ) )
        self.assertEquals( num_t_remaining, num_t - 1 )
        self.assertEquals( len( chain ), num_elements - 1 )

        chain = Chain( self.a[:] )
        elements = self.a.strip().split()
        num_elements = len( elements )
        num_u = chain.data.count( Hetero( 'u' ) )
        del chain[ 12 ]
        num_u_remaining = 0
        self.assertEquals( num_u_remaining, num_u - 1 )
        self.assertEquals( len( chain ), num_elements - 1 )

        chain = Chain( self.b[:] )
        elements = self.b.strip().split()
        num_elements = len( elements )
        num_c = chain.data.count( Hetero( 'c' ) )
        del chain[ 0 ]
        num_c_remaining = chain.data.count( Hetero( 'c' ) )
        self.assertEquals( num_c_remaining, num_c - 1 )
        self.assertEquals( len( chain ), num_elements - 1 )

        chain = Chain( self.b[:] )
        elements = self.b.strip().split()
        num_elements = len( elements )
        num_g = chain.data.count( Hetero( 't' ) )
        del chain[ 6 ]
        num_g_remaining = chain.data.count( Hetero( 't' ) )
        self.assertEquals( num_g_remaining, num_g - 1 )
        self.assertEquals( len( chain ), num_elements - 1 )

        chain = Chain( self.c[:] )
        elements = self.c.strip().split()
        num_elements = len( elements )
        num_thr = chain.data.count( Hetero( 'thr' ) )
        del chain[ 0 ]
        num_thr_remaining = chain.data.count( Hetero( 'thr' ) )
        self.assertEquals( num_thr_remaining, num_thr - 1 )
        self.assertEquals( len( chain ), num_elements - 1 )

    def testGetSlice( self ):
        chain = Chain( self.a[:]  )
        first = 0
        last = len( chain )
        slice = chain[ : ]
        other = chain.data[:]
        self.assertEquals( slice.data, other )


        chain = Chain( self.a[:]  )
        first = 0
        last = 4
        slice = chain[ first : last ]
        other = chain.data[ first : last ]
        self.assertEquals( slice.data, other )

        chain = Chain( self.b[:] )
        first = 2
        last = len( chain )
        slice = chain[ first: last ]
        other = chain.data[ first : last ]
        self.assertEquals( slice.data, other )

        chain = Chain( self.b[:] )
        first = -1
        slice = chain[ first : ]
        other = chain.data[ first:  ]
        self.assertEquals( slice.data, other )

        chain = Chain( self.c[:] )
        first = 3
        last = 7
        slice = chain[ first : last ]
        other = chain.data[ first: last ]
        self.assertEquals( slice.data, other )

        chain = Chain( self.c[:] )
        first = 3
        last = -1
        slice = chain[ first : last ]
        other = chain.data[ first: last ]
        self.assertEquals( slice.data, other )


    def testSetSlice( self ):
        chain = Chain( self.a[:]  )
        slice = 'G T C A G 5NC G C A T G G'
        chain[ : ] = slice[ 4 : 7 ]
        other = Chain( slice[ 4: 7 ]  )
        self.assertEquals( chain, other )


        chain = Chain( self.c[:]  )
        old_chain = Chain( self.c[ : ] )
        slice = 'MET ILE GLU ILE LYS ASP'
        chain[ 2 : 5 ] = slice
        other = Chain( old_chain.data[ :2 ] + Chain( slice ).data + old_chain.data[ 5: ] )
        self.assertEquals( chain, other )

        chain = Chain( self.c[:]  )
        old_chain = Chain( self.c[ : ] )
        slice = 'CYS GLY ALA GLU CYS VAL TYR'
        chain[ 7 : ] = slice
        other = Chain( old_chain.data[ :7 ] + Chain( slice ).data )
        self.assertEquals( chain, other )

        chain = Chain( self.c[:]  )
        old_chain = Chain( self.c[ : ] )
        slice = 'SER ASN GLU TRP ASP '
        chain[ : 9 ] = slice
        other = Chain( Chain( slice ).data + old_chain.data[ 9: ] )
        self.assertEquals( chain, other )

    def testDelSlice( self ):
        chain = Chain( self.c[ : ] )
        old_chain = Chain( self.c[ : ] )
        del chain[ 3 : 8 ]
        other = Chain( old_chain.data[ :3 ] + old_chain.data[ 8: ] )
        self.assertEquals( chain, other )

        chain = Chain( self.c[:]  )
        old_chain = Chain( self.c[ : ] )
        del chain[ :4 ]
        other = Chain( old_chain.data[ 4: ] )
        self.assertEquals( chain, other )

        chain = Chain( self.c[:]  )
        old_chain = Chain( self.c[ : ] )
        del chain[ 9: ]
        other = Chain( old_chain.data[ :9 ] )
        self.assertEquals( chain, other )


    def testContains( self ):
        chain = Chain( self.c[ : ] )
        self.failIf( 'ser' in chain )
        self.failUnless( 'lys' in chain )
        self.failUnless( 'asp' in chain )

    def testAdd( self ):
        texta = 'G U G G U C U G A U G A G G C C'
        textb = 'G G C C G A A A C U C G U A A G A G U C A C C A C'
        targeta = texta + Chain( textb )
        targetb = Chain( texta ) + textb
        targetc = Chain( texta ) + Chain( textb )
        self.assertEquals( targeta, targetc )
        self.assertEquals( targetb, targetc )
        self.assertEquals( targeta, targetb )
        self.assertEquals( len( targeta ), len( Chain( texta ) ) + len( Chain( textb ) ) )

        targetd = Chain( texta )
        targetd += textb

        targete = Chain( texta )
        targete += Chain( textb )
        self.assertEquals( targetd, targetc )
        self.assertEquals( targete, targetb )




class CrystalTestCase(unittest.TestCase):

    def setUp( self ):

        self.crystal = Crystal( { 'a' : 'T T G A C T C T C T T A A', \
                             'b' : Chain( 'G A G A G T C A' ), \
                             'c' : 'T T G A C T C T C T T A A', \
                             'd' : Chain( 'G A G A G T C A' )
                            } )

    def testLen( self ):
        self.assertEquals( len( self.crystal ), len( self.crystal.data ) )

    def testGetItem( self ):
        self.assertEquals( self.crystal[ 'a' ], self.crystal.data[ 'a' ] )

    def testSetItem( self ):
        target = copy.deepcopy( self.crystal )
        e = 'MET ALA LEU THR ASN ALA GLN ILE LEU ALA VAL ILE ASP SER'
        f = 'LEU GLY GLY GLY LEU GLN GLY THR LEU HIS CYS TYR GLU ILE PRO LEU'
        target[ 'e' ] = e
        target[ 'f' ] = Chain( f )
        self.assertEquals( Chain( e ), target[ 'e' ] )
        self.assertEquals( Chain( f ), target[ 'f' ] )

    def testDelItem( self ):
        target = copy.deepcopy( self.crystal )
        del target[ 'b' ]
        self.failIf( target.data.has_key( 'b' ) )
        self.failUnless( target.data.has_key( 'a' ) )
        self.failUnless( target.data.has_key( 'c' ) )

    def testClear( self ):
        target = copy.deepcopy( self.crystal )
        target.clear()
        self.assertEquals( len( target.data ), 0 )

    def testKeys( self ):
        self.assertEquals( self.crystal.keys(), self.crystal.data.keys() )

    def testValues( self ):
        self.assertEquals( self.crystal.values(), self.crystal.data.values() )

    def testItems( self ):
        self.assertEquals( self.crystal.items(), self.crystal.data.items() )

    def testKeys( self ):
        self.assertEquals( self.crystal.keys(), self.crystal.data.keys() )

    def testHasKey( self ):
        self.failUnless( self.crystal.has_key( 'b' ) )
        self.failUnless( self.crystal.has_key( 'c' ) )
        self.failIf( self.crystal.has_key( 'z' ) )


class HeteroTestCase(unittest.TestCase):

    def testInit( self ):
        self.assertRaises( CrystalError, Hetero, 'abcd' )
        self.assertRaises( CrystalError, Hetero, '' )
        self.assertRaises( CrystalError, Hetero, 'A@#' )
        self.assertRaises( CrystalError, Hetero, [] )
        self.assertRaises( CrystalError, Hetero, {} )

    def testLen( self ):
        bru = Hetero( 'bru' )
        self.assertEquals( len( bru ), 3 )
        _14w = Hetero( '14w' )
        self.assertEquals( len( _14w ), 3 )
        a = Hetero( 'a' )
        self.assertEquals( len( a ), 1 )
        ga = Hetero( 'ga' )
        self.assertEquals( len( ga ), 2 )

    def testEquals( self ):
        u = Hetero( 'u' )
        u1 = Hetero( 'u' )
        self.assertEquals( u, u1 )
        self.assertEquals( u, Hetero( 'U' ) )
        self.assertNotEquals( u, Hetero( 'u1' ) )
        self.assertNotEquals( u, Hetero( 'x' ) )
        gna = Hetero( 'gna' )
        self.assertEquals( gna, Hetero( 'gNA' ) )
        self.assertEquals( gna, Hetero( 'GnA' ) )
        self.assertNotEquals( gna, Hetero( 'gnb' ) )
        self.assertNotEquals( gna, Hetero( 'na' ) )

class CrystalTestSuite(unittest.TestSuite):

    def __init__(self):
        unittest.TestSuite.__init__(self,
                                    map(ChainTestCase, (  "testLen",
                                                          "testEquals",
                                                          "testAppend",
                                                          "testInsert",
                                                          "testRemove",
                                                          "testIndex",
                                                          "testCount",
                                                          "testGetItem",
                                                          "testSetItem",
                                                          "testDelItem",
                                                          "testGetSlice",
                                                          "testSetSlice",
                                                          "testDelSlice",
                                                          "testContains",
                                                          "testAdd" )) + \
                                    map(HeteroTestCase, ("testEquals",
                                                             "testInit",
                                                             "testLen")) + \
                                    map(CrystalTestCase, ("testLen",
                                                         "testGetItem",
                                                         "testSetItem",
                                                         "testDelItem",
                                                         "testClear",
                                                         "testKeys",
                                                         "testValues",
                                                         "testItems",
                                                         "testHasKey")))


def run_tests(argv):
    runner = unittest.TextTestRunner(sys.stdout, verbosity = 2)
    runner.run(CrystalTestSuite())

if __name__ == "__main__":
    sys.exit(run_tests(sys.argv))
