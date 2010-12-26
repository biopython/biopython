# Copyright 2002 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# python unittest framework
import unittest
import copy
import sys

# modules to be tested
from Bio.Crystal import Hetero, Chain, Crystal, CrystalError

class ChainTestCase(unittest.TestCase):

    def setUp(self):
        self.a = 'C A A C T A G G T C A C U A G G T C A G'
        self.b = 'C T G A C C T A G T G A C C T A G T T G'
        self.c = 'THR LYS LEU ASN GLY MET VAL LEU LEU CYS LYS VAL CYS GLY ASP'
        self.d = 'THR LYS LEU ASN GLY MET VAL LEU LEU CYS LYS VAL CYS GLY ASP '
        self.e = 'TYR LYS LEU ASN GLY MET VAL LEU LEU CYS LYS VAL CYS GLY ASP '
        self.f = 'THR LYS LEU ASN GLY MET VAL LEU LEU CYS LYS VAL CYS GLY SER '
        self.g = 'C A A C T A G G T C A C U A G G T C A T'
        self.h = 'G A A C T A G G T C A C U A G G T C A G'

    def testEquals(self):
        first = Chain(self.a)
        second = Chain(self.a)
        self.assertEqual(first, second)

        first = Chain(self.b)
        second = Chain(self.b)
        self.assertEqual(first, second)

        first = Chain(self.c)
        second = Chain(self.c)
        self.assertEqual(first, second)

        first = Chain(self.a)
        second = Chain(self.g)
        self.assertNotEqual(first, second)

        first = Chain(self.a)
        second = Chain(self.h)
        self.assertNotEqual(first, second)

        first = Chain(self.c)
        second = Chain(self.e)
        self.assertNotEqual(first, second)

        first = Chain(self.c)
        second = Chain(self.f)
        self.assertNotEqual(first, second)


    def testLen(self):
        chain = Chain(self.a)
        elements = self.a.strip().split()
        num_elements = len(elements)
        self.assertEqual(len(chain), num_elements)

        chain = Chain(self.b)
        elements = self.b.strip().split()
        num_elements = len(elements)
        self.assertEqual(len(chain), num_elements)

        chain = Chain(self.c)
        elements = self.c.strip().split()
        num_elements = len(elements)
        self.assertEqual(len(chain), num_elements)


    def testAppend(self):
        chain = Chain(self.a[:])
        chain.append('U')
        elements = self.a.strip().split()
        num_elements = len(elements)
        last_element = chain.data[ -1 ]
        self.assertEqual('u', last_element.data)
        self.assertEqual(len(chain), num_elements + 1)

        chain = Chain(self.a[:])
        chain.append(Hetero('A'))
        elements = self.a.strip().split()
        num_elements = len(elements)
        last_element = chain.data[ -1 ]
        self.assertEqual('a', last_element.data)
        self.assertEqual(len(chain), num_elements + 1)

        chain = Chain(self.b[:])
        chain.append('t')
        elements = self.b.strip().split()
        num_elements = len(elements)
        last_element = chain.data[ -1 ]
        self.assertEqual('t', last_element.data)
        self.assertEqual(len(chain), num_elements + 1)

        chain = Chain(self.b[:])
        chain.append(Hetero('C'))
        elements = self.b.strip().split()
        num_elements = len(elements)
        last_element = chain.data[ -1 ]
        self.assertEqual('c', last_element.data)
        self.assertEqual(len(chain), num_elements + 1)

        chain = Chain(self.c[:])
        chain.append('ser')
        elements = self.c.strip().split()
        num_elements = len(elements)
        last_element = chain.data[ -1 ]
        self.assertEqual('ser', last_element.data)
        self.assertEqual(len(chain), num_elements + 1)


    def testInsert(self):
        chain = Chain(self.a[:])
        i = 4
        chain.insert(i, 'g')
        elements = self.a.strip().split()
        num_elements = len(elements)
        target_element = chain.data[ i ]
        self.assertEqual('g', target_element.data)
        self.assertEqual(len(chain), num_elements + 1)

        chain = Chain(self.a[:])
        i = 0
        chain.insert(i, 't')
        elements = self.a.strip().split()
        num_elements = len(elements)
        target_element = chain.data[ i ]
        self.assertEqual('t', target_element.data)
        self.assertEqual(len(chain), num_elements + 1)

        chain = Chain(self.b[:])
        i = 9
        chain.insert(i, Hetero('a'))
        elements = self.a.strip().split()
        num_elements = len(elements)
        target_element = chain.data[ i ]
        self.assertEqual('a', target_element.data)
        self.assertEqual(len(chain), num_elements + 1)

        chain = Chain(self.c[:])
        i = 5
        chain.insert(i, 'gln')
        elements = self.c.strip().split()
        num_elements = len(elements)
        target_element = chain.data[ i ]
        self.assertEqual('gln', target_element.data)
        self.assertEqual(len(chain), num_elements + 1)


    def testRemove(self):

        chain = Chain(self.a[:])
        elements = self.a.strip().split()
        num_elements = len(elements)
        num_a = chain.data.count(Hetero('a'))
        chain.remove('a')
        num_a_remaining = chain.data.count(Hetero('a'))
        self.assertEqual(num_a_remaining, num_a - 1)
        self.assertEqual(len(chain), num_elements - 1)

        chain = Chain(self.b[:])
        elements = self.b.strip().split()
        num_elements = len(elements)
        num_b = chain.data.count(Hetero('t'))
        chain.remove('t')
        num_b_remaining = chain.data.count(Hetero('t'))
        self.assertEqual(num_b_remaining, num_b - 1)
        self.assertEqual(len(chain), num_elements - 1)

        chain = Chain(self.c[:])
        elements = self.c.strip().split()
        num_elements = len(elements)
        num_leu = chain.data.count(Hetero('leu'))
        chain.remove('leu')
        num_leu_remaining = chain.data.count(Hetero('leu'))
        self.assertEqual(num_leu_remaining, num_leu - 1)
        self.assertEqual(len(chain), num_elements - 1)


    def testCount(self):
        chain = Chain(self.a[:])
        num_a = chain.data.count(Hetero('a'))
        self.assertEqual(chain.count('a'), num_a)

        chain = Chain(self.b[:])
        num_a = chain.data.count(Hetero('t'))
        self.assertEqual(chain.count('t'), num_a)

        chain = Chain(self.c[:])
        num_a = chain.data.count(Hetero('leu'))
        self.assertEqual(chain.count('leu'), num_a)

        chain = Chain(self.c[:])
        num_a = chain.data.count(Hetero('cys'))
        self.assertEqual(chain.count('cys'), num_a)


    def testIndex(self):
        chain = Chain(self.a[:])
        index_g = chain.data.index(Hetero('g'))
        self.assertEqual(chain.index('g'), index_g)

        chain = Chain(self.b[:])
        index_c = chain.data.index(Hetero('c'))
        self.assertEqual(chain.index('c'), index_c)

        chain = Chain(self.c[:])
        index_met = chain.data.index(Hetero('met'))
        self.assertEqual(chain.index('met'), index_met)

    def testGetItem(self):
        chain = Chain(self.a[:])
        element_3 = chain.data[ 3 ]
        self.assertEqual(chain[ 3 ], element_3)

        chain = Chain(self.a[:])
        element_0 = chain.data[ 0 ]
        self.assertEqual(chain[ 0 ], element_0)

        chain = Chain(self.b[:])
        element_7 = chain.data[ 7 ]
        self.assertEqual(chain[ 7 ], element_7)

        chain = Chain(self.b[:])
        last_element = chain.data[ -1 ]
        self.assertEqual(chain[ -1 ], last_element)

        chain = Chain(self.c[:])
        element_8 = chain.data[ 8 ]
        self.assertEqual(chain[ 8 ], element_8)

    def testSetItem(self):
        chain = Chain(self.a[:])
        chain[ 2 ] = 't'
        element_2 = chain.data[ 2 ]
        self.assertEqual(chain[ 2 ], element_2)

        chain = Chain(self.a[:])
        chain[ 0 ] = Hetero('U')
        element_0 = chain.data[ 0 ]
        self.assertEqual(chain[ 0 ], element_0)

        chain = Chain(self.b[:])
        chain[ -1 ] = Hetero('c')
        last_element = chain.data[ -1 ]
        self.assertEqual(chain[ -1 ], last_element)

        chain = Chain(self.b[:])
        chain[ 1 ] = 'a'
        element_1 = chain.data[ 1 ]
        self.assertEqual(chain[ 1 ], element_1)

        chain = Chain(self.c[:])
        chain[ 5 ] = 'ser'
        element_5 = chain.data[ 5 ]
        self.assertEqual(chain[ 5 ], element_5)

    def testDelItem(self):

        chain = Chain(self.a[:])
        elements = self.a.strip().split()
        num_elements = len(elements)
        num_t = chain.data.count(Hetero('t'))
        del chain[ 4 ]
        num_t_remaining = chain.data.count(Hetero('t'))
        self.assertEqual(num_t_remaining, num_t - 1)
        self.assertEqual(len(chain), num_elements - 1)

        chain = Chain(self.a[:])
        elements = self.a.strip().split()
        num_elements = len(elements)
        num_u = chain.data.count(Hetero('u'))
        del chain[ 12 ]
        num_u_remaining = 0
        self.assertEqual(num_u_remaining, num_u - 1)
        self.assertEqual(len(chain), num_elements - 1)

        chain = Chain(self.b[:])
        elements = self.b.strip().split()
        num_elements = len(elements)
        num_c = chain.data.count(Hetero('c'))
        del chain[ 0 ]
        num_c_remaining = chain.data.count(Hetero('c'))
        self.assertEqual(num_c_remaining, num_c - 1)
        self.assertEqual(len(chain), num_elements - 1)

        chain = Chain(self.b[:])
        elements = self.b.strip().split()
        num_elements = len(elements)
        num_g = chain.data.count(Hetero('t'))
        del chain[ 6 ]
        num_g_remaining = chain.data.count(Hetero('t'))
        self.assertEqual(num_g_remaining, num_g - 1)
        self.assertEqual(len(chain), num_elements - 1)

        chain = Chain(self.c[:])
        elements = self.c.strip().split()
        num_elements = len(elements)
        num_thr = chain.data.count(Hetero('thr'))
        del chain[ 0 ]
        num_thr_remaining = chain.data.count(Hetero('thr'))
        self.assertEqual(num_thr_remaining, num_thr - 1)
        self.assertEqual(len(chain), num_elements - 1)

    def testGetSlice(self):
        chain = Chain(self.a[:] )
        first = 0
        last = len(chain)
        slice = chain[ : ]
        other = chain.data[:]
        self.assertEqual(slice.data, other)


        chain = Chain(self.a[:] )
        first = 0
        last = 4
        slice = chain[ first : last ]
        other = chain.data[ first : last ]
        self.assertEqual(slice.data, other)

        chain = Chain(self.b[:])
        first = 2
        last = len(chain)
        slice = chain[ first: last ]
        other = chain.data[ first : last ]
        self.assertEqual(slice.data, other)

        chain = Chain(self.b[:])
        first = -1
        slice = chain[ first : ]
        other = chain.data[ first:  ]
        self.assertEqual(slice.data, other)

        chain = Chain(self.c[:])
        first = 3
        last = 7
        slice = chain[ first : last ]
        other = chain.data[ first: last ]
        self.assertEqual(slice.data, other)

        chain = Chain(self.c[:])
        first = 3
        last = -1
        slice = chain[ first : last ]
        other = chain.data[ first: last ]
        self.assertEqual(slice.data, other)


    def testSetSlice(self):
        chain = Chain(self.a[:] )
        slice = 'G T C A G 5NC G C A T G G'
        chain[ : ] = slice[ 4 : 7 ]
        other = Chain(slice[ 4: 7 ] )
        self.assertEqual(chain, other)


        chain = Chain(self.c[:] )
        old_chain = Chain(self.c[ : ])
        slice = 'MET ILE GLU ILE LYS ASP'
        chain[ 2 : 5 ] = slice
        other = Chain(old_chain.data[ :2 ] + Chain(slice).data + old_chain.data[ 5: ])
        self.assertEqual(chain, other)

        chain = Chain(self.c[:] )
        old_chain = Chain(self.c[ : ])
        slice = 'CYS GLY ALA GLU CYS VAL TYR'
        chain[ 7 : ] = slice
        other = Chain(old_chain.data[ :7 ] + Chain(slice).data)
        self.assertEqual(chain, other)

        chain = Chain(self.c[:] )
        old_chain = Chain(self.c[ : ])
        slice = 'SER ASN GLU TRP ASP '
        chain[ : 9 ] = slice
        other = Chain(Chain(slice).data + old_chain.data[ 9: ])
        self.assertEqual(chain, other)

    def testDelSlice(self):
        chain = Chain(self.c[ : ])
        old_chain = Chain(self.c[ : ])
        del chain[ 3 : 8 ]
        other = Chain(old_chain.data[ :3 ] + old_chain.data[ 8: ])
        self.assertEqual(chain, other)

        chain = Chain(self.c[:] )
        old_chain = Chain(self.c[ : ])
        del chain[ :4 ]
        other = Chain(old_chain.data[ 4: ])
        self.assertEqual(chain, other)

        chain = Chain(self.c[:] )
        old_chain = Chain(self.c[ : ])
        del chain[ 9: ]
        other = Chain(old_chain.data[ :9 ])
        self.assertEqual(chain, other)


    def testContains(self):
        chain = Chain(self.c[ : ])
        self.assertFalse('ser' in chain)
        self.assertTrue('lys' in chain)
        self.assertTrue('asp' in chain)

    def testAdd(self):
        texta = 'G U G G U C U G A U G A G G C C'
        textb = 'G G C C G A A A C U C G U A A G A G U C A C C A C'
        targeta = texta + Chain(textb)
        targetb = Chain(texta) + textb
        targetc = Chain(texta) + Chain(textb)
        self.assertEqual(targeta, targetc)
        self.assertEqual(targetb, targetc)
        self.assertEqual(targeta, targetb)
        self.assertEqual(len(targeta), len(Chain(texta)) + len(Chain(textb)))

        targetd = Chain(texta)
        targetd += textb

        targete = Chain(texta)
        targete += Chain(textb)
        self.assertEqual(targetd, targetc)
        self.assertEqual(targete, targetb)




class CrystalTestCase(unittest.TestCase):

    def setUp(self):

        self.crystal = Crystal({ 'a' : 'T T G A C T C T C T T A A', \
                             'b' : Chain('G A G A G T C A'), \
                             'c' : 'T T G A C T C T C T T A A', \
                             'd' : Chain('G A G A G T C A')
                            })

    def testLen(self):
        self.assertEqual(len(self.crystal), len(self.crystal.data))

    def testGetItem(self):
        self.assertEqual(self.crystal[ 'a' ], self.crystal.data[ 'a' ])

    def testSetItem(self):
        target = copy.deepcopy(self.crystal)
        e = 'MET ALA LEU THR ASN ALA GLN ILE LEU ALA VAL ILE ASP SER'
        f = 'LEU GLY GLY GLY LEU GLN GLY THR LEU HIS CYS TYR GLU ILE PRO LEU'
        target[ 'e' ] = e
        target[ 'f' ] = Chain(f)
        self.assertEqual(Chain(e), target[ 'e' ])
        self.assertEqual(Chain(f), target[ 'f' ])

    def testDelItem(self):
        target = copy.deepcopy(self.crystal)
        del target[ 'b' ]
        self.assertFalse('b' in target.data)
        self.assertTrue('a' in target.data)
        self.assertTrue('c' in target.data)

    def testClear(self):
        target = copy.deepcopy(self.crystal)
        target.clear()
        self.assertEqual(len(target.data), 0)

    def testKeys(self):
        self.assertEqual(self.crystal.keys(), self.crystal.data.keys())

    def testValues(self):
        self.assertEqual(self.crystal.values(), self.crystal.data.values())

    def testItems(self):
        self.assertEqual(self.crystal.items(), self.crystal.data.items())

    def testKeys(self):
        self.assertEqual(self.crystal.keys(), self.crystal.data.keys())

    def testHasKey(self):
        self.assertTrue('b' in self.crystal)
        self.assertTrue('c' in self.crystal)
        self.assertFalse('z' in self.crystal)


class HeteroTestCase(unittest.TestCase):

    def testInit(self):
        self.assertRaises(CrystalError, Hetero, 'abcd')
        self.assertRaises(CrystalError, Hetero, '')
        self.assertRaises(CrystalError, Hetero, 'A@#')
        self.assertRaises(CrystalError, Hetero, [])
        self.assertRaises(CrystalError, Hetero, {})

    def testLen(self):
        bru = Hetero('bru')
        self.assertEqual(len(bru), 3)
        _14w = Hetero('14w')
        self.assertEqual(len(_14w), 3)
        a = Hetero('a')
        self.assertEqual(len(a), 1)
        ga = Hetero('ga')
        self.assertEqual(len(ga), 2)

    def testEquals(self):
        u = Hetero('u')
        u1 = Hetero('u')
        self.assertEqual(u, u1)
        self.assertEqual(u, Hetero('U'))
        self.assertNotEqual(u, Hetero('u1'))
        self.assertNotEqual(u, Hetero('x'))
        gna = Hetero('gna')
        self.assertEqual(gna, Hetero('gNA'))
        self.assertEqual(gna, Hetero('GnA'))
        self.assertNotEqual(gna, Hetero('gnb'))
        self.assertNotEqual(gna, Hetero('na'))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)

