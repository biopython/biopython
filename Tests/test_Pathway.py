# Copyright 2001 by Tarjei Mikkelsen.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# python unittest framework
import unittest
import sys

# modules to be tested
from Bio.Pathway import *
from Bio.Pathway.Rep.HashSet import *
from Bio.Pathway.Rep.Graph import *
from Bio.Pathway.Rep.MultiGraph import *


class HashSetTestCase(unittest.TestCase):

    def testEquals(self):
        self.assertEqual(HashSet(['a','b']), HashSet(['b','a']), "not equal to similar")
        self.assertEqual(HashSet(), HashSet(), "empty set not equal to similar")
        self.assertNotEqual(HashSet(['a','b','c']), HashSet(), "non-empty equal to empty")
        self.assertNotEqual(HashSet(['a','b']), HashSet(['a','b','c']), "equal to superset")
        self.assertNotEqual(HashSet(['a','b']), HashSet(['a']), "equal to subset")
        
    def testLen(self):
        a = HashSet()
        self.assertEqual(len(a), 0, "incorrect default size")
        a.add('a')
        a.add('b')
        self.assertEqual(len(a), 2, "incorrect size")
        a.remove('b')
        self.assertEqual(len(a), 1, "incorrect size after removal")
        a.add('a')
        self.assertEqual(len(a), 1, "incorrect size after duplicate add")

    def testContains(self):
        n = HashSet()
        self.assertTrue('a' not in n, "element in empty set")
        self.assertTrue(not n.contains('a'), "element in empty set (2)")
        a = HashSet(['a','b','c','d'])
        self.assertTrue('a' in a, "contained element not found")
        self.assertTrue('d' in a, "contained element not found")
        self.assertTrue('e' not in a, "not contained element found")
        self.assertTrue(68 not in a, "not contained element found")
        self.assertTrue(a.contains('a'), "contained element not found (2)")
        self.assertTrue(a.contains('d'), "contained element not found (2)")
        self.assertTrue(not a.contains('e'), "not contained element found (2)")
        self.assertTrue(not a.contains(68), "not contained element found (2)")
        
    def testList(self):
        a = HashSet(['a', 'b', 'c', 'd', 'e'])
        l = a.list()
        l.sort()
        self.assertEqual(l, ['a', 'b', 'c', 'd', 'e'], "incorrect list")
        l = []
        self.assertTrue('e' in a, "set rep exposure")

    def testSetOps(self):
        n = HashSet()
        a = HashSet(['a', 'b', 'c'])
        b = HashSet(['a', 'd', 'e', 'f'])
        c = HashSet(['g', 'h'])
        # union
        self.assertEqual(a.union(b), HashSet(['a','b','c','d','e','f']), "incorrect union")
        self.assertEqual(a.union(n), a, "incorrect union with empty set")
        # intersection
        self.assertEqual(a.intersection(b), HashSet(['a']), "incorrect intersection")
        self.assertEqual(a.intersection(c), HashSet(), "incorrect intersection")
        self.assertEqual(a.intersection(n), HashSet(), "incorrect intersection with empty set")
        # difference
        self.assertEqual(a.difference(b), HashSet(['b','c']), "incorrect difference")
        self.assertEqual(a.difference(c), HashSet(['a','b','c']), "incorrect difference")
        self.assertEqual(b.difference(a), HashSet(['d','e','f']), "incorrect difference")
        # cartesian product
        self.assertEqual(a.cartesian(c), HashSet([('a','g'),('a','h'),
                                                  ('b','g'),('b','h'),
                                                  ('c','g'),('c','h')]),
                         "incorrect cartesian product")
        self.assertEqual(a.cartesian(n), HashSet(), "incorrect cartesian product")


class GraphTestCase(unittest.TestCase):

    def testEquals(self):
        a = Graph(['a','b','c'])
        a.add_edge('a','b','label1')
        a.add_edge('b','c','label1')
        a.add_edge('b','a','label2')
        b = Graph(['a','b','c'])
        self.assertNotEqual(a, b, "equal to similar nodes, no edges")
        b.add_edge('a','b','label1')
        self.assertNotEqual(a, b, "equal to similar nodes, edge subset")
        b.add_edge('b','c','label1')
        b.add_edge('b','a','label2')
        self.assertEqual(a, b, "not equal to similar")
        c = Graph(['a','b','c'])
        c.add_edge('a','b','label2')
        c.add_edge('b','c','label2')
        c.add_edge('b','a','label1')
        self.assertNotEqual(a, c, "equal to similar with different labels")
        self.assertNotEqual(c, Graph(), "equal to empty graph")
        self.assertEqual(Graph(), Graph(), "empty graph not equal to self")

    def testNodes(self):
        a = Graph()
        self.assertEqual(a.nodes(), [], "default graph not empty")
        a.add_node('a')
        self.assertEqual(a.nodes(), ['a'], "one node not added")
        a.add_node('a')
        self.assertEqual(a.nodes(), ['a'], "duplicate node added")
        a.add_node('b')
        l = a.nodes()
        l.sort()
        self.assertEqual(l, ['a', 'b'], "second node not added")

    def testEdges(self):
        a = Graph(['a','b','c','d'])
        a.add_edge('a','b','label1')
        self.assertEqual(a.child_edges('a'), [('b','label1')], "incorrect child edges")
        a.add_edge('b','a','label2')
        self.assertEqual(a.parent_edges('a'), [('b','label2')], "incorrect parent edges")
        a.add_edge('b','c','label3')
        self.assertEqual(a.parent_edges('c'), [('b','label3')], "incorrect parent edges")
        l = a.children('b')
        l.sort()
        self.assertEqual(l, ['a', 'c'], "incorrect children")
        self.assertEqual(a.children('d'), [], "incorrect children for singleton")
        self.assertEqual(a.parents('a'), ['b'], "incorrect parents")

    def testRemoveNode(self):
        a = Graph(['a','b','c','d','e'])
        a.add_edge('a','e','label1')
        a.add_edge('b','e','label1')
        a.add_edge('c','e','label2')
        a.add_edge('d','e','label3')
        a.add_edge('e','d','label4')
        a.add_edge('a','b','label5')
        a.remove_node('e')
        b = Graph(['a','b','c','d'])
        b.add_edge('a','b','label5')
        self.assertEqual(a, b)#, "incorrect node removal")


class MultiGraphTestCase(unittest.TestCase):

    def testEquals(self):
        a = MultiGraph(['a','b','c'])
        a.add_edge('a','b','label1')
        a.add_edge('b','c','label1')
        a.add_edge('b','a','label2')
        b = MultiGraph(['a','b','c'])
        self.assertNotEqual(a, b, "equal to similar nodes, no edges")
        b.add_edge('a','b','label1')
        self.assertNotEqual(a, b, "equal to similar nodes, edge subset")
        b.add_edge('b','c','label1')
        b.add_edge('b','a','label2')
        self.assertEqual(a, b, "not equal to similar")
        c = MultiGraph(['a','b','c'])
        c.add_edge('a','b','label2')
        c.add_edge('b','c','label2')
        c.add_edge('b','a','label1')
        self.assertNotEqual(a, c, "equal to similar with different labels")
        self.assertNotEqual(c, MultiGraph(), "equal to empty graph")
        self.assertEqual(MultiGraph(), MultiGraph(), "empty graph not equal to self")

    def testNodes(self):
        a = MultiGraph()
        self.assertEqual(a.nodes(), [], "default graph not empty")
        a.add_node('a')
        self.assertEqual(a.nodes(), ['a'], "one node not added")
        a.add_node('a')
        self.assertEqual(a.nodes(), ['a'], "duplicate node added")
        a.add_node('b')
        l = a.nodes()
        l.sort()
        self.assertEqual(l, ['a', 'b'], "second node not added")

    def testEdges(self):
        a = MultiGraph(['a','b','c','d'])
        a.add_edge('a','b','label1')
        self.assertEqual(a.child_edges('a'), [('b','label1')], "incorrect child edges")
        a.add_edge('a','b','label2')
        l = a.child_edges('a')
        l.sort()
        self.assertEqual(l, [('b','label1'),('b','label2')], "incorrect child edges")
        a.add_edge('b','a','label2')
        self.assertEqual(a.parent_edges('a'), [('b','label2')], "incorrect parent edges")
        a.add_edge('b','c','label3')
        self.assertEqual(a.parent_edges('c'), [('b','label3')], "incorrect parent edges")
        l = a.children('b')
        l.sort()
        self.assertEqual(l, ['a', 'c'], "incorrect children")
        self.assertEqual(a.children('d'), [], "incorrect children for singleton")
        self.assertEqual(a.parents('a'), ['b'], "incorrect parents")

    def testRemoveNode(self):
        a = MultiGraph(['a','b','c','d','e'])
        a.add_edge('a','e','label1')
        a.add_edge('b','e','label1')
        a.add_edge('c','e','label2')
        a.add_edge('d','e','label3')
        a.add_edge('e','d','label4')
        a.add_edge('a','b','label5')
        a.remove_node('e')
        b = MultiGraph(['a','b','c','d'])
        b.add_edge('a','b','label5')
        self.assertEqual(a, b)#, "incorrect node removal")

        
class ReactionTestCase(unittest.TestCase):

    def setUp(self):
        self.r_empty = Reaction()
        self.r_prod = Reaction({"a":1})
        self.r_dest = Reaction({"a":-1})
        self.r_1 = Reaction({"a":-1, "b":1})
        self.r_1i = Reaction({"a":-1, "b":1, "c":0})
        self.r_2 = Reaction({"b":-1, "c":1})
        self.r_3 = Reaction({"a":-1, "d":2})
        self.r_4 = Reaction({"c":-1, "d":-1, "a":1, "e":2})

    def testEq(self):
        self.assertEqual(self.r_1, self.r_1i, "not equal to similar")
        self.assertNotEquals(self.r_3, self.r_4, "equal to different")
        
    def testRev(self):
        self.assertEqual(self.r_empty.reverse(), self.r_empty, "empty reversed not empty")
        self.assertEqual(self.r_prod.reverse(), self.r_dest,
                          "reversed reaction not equal to similar")
        self.assertEqual(self.r_4.reverse(), Reaction({"c":1, "d":1, "a":-1, "e":-2}),
                         "reversed reaction not equal to similar")
        self.assertEqual(self.r_3.reverse().reverse(), self.r_3,
                          "double reversal not identity")
    

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
