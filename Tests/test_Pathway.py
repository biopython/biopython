# Copyright 2001 by Tarjei Mikkelsen.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# python unittest framework
import unittest
import sys

# modules to be tested
from Bio.Pathway import Reaction
from Bio.Pathway.Rep.Graph import Graph
from Bio.Pathway.Rep.MultiGraph import MultiGraph


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
        self.assertEqual(a.child_edges('a'), [('b','label1')]) #, "incorrect child edges")
        a.add_edge('b','a','label2')
        self.assertEqual(a.parent_edges('a'), [('b','label2')]) #, "incorrect parent edges")
        a.add_edge('b','c','label3')
        self.assertEqual(a.parent_edges('c'), [('b','label3')]) #, "incorrect parent edges")
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
        self.assertEqual(a.child_edges('a'), [('b','label1')]) #, "incorrect child edges")
        a.add_edge('a','b','label2')
        l = a.child_edges('a')
        l.sort()
        self.assertEqual(l, [('b','label1'),('b','label2')]) #, "incorrect child edges")
        a.add_edge('b','a','label2')
        self.assertEqual(a.parent_edges('a'), [('b','label2')]) #, "incorrect parent edges")
        a.add_edge('b','c','label3')
        self.assertEqual(a.parent_edges('c'), [('b','label3')]) #, "incorrect parent edges")
        l = a.children('b')
        l.sort()
        self.assertEqual(l, ['a', 'c'], "incorrect children")
        self.assertEqual(a.children('d'), [], "incorrect children for singleton")
        self.assertEqual(a.parents('a'), ['b'], "incorrect parents")

    def testRemoveNode(self):
        a = MultiGraph(['a','b','c','d','e'])
        a.add_edge('a','e','label1')
        self.assertEqual(repr(a), "<MultiGraph: ('a': ('e', 'label1'))('b': )('c': )('d': )('e': )>")
        a.add_edge('b','e','label1')
        a.add_edge('c','e','label2')
        a.add_edge('d','e','label3')
        a.add_edge('e','d','label4')
        a.add_edge('a','b','label5')
        self.assertEqual(repr(a), "<MultiGraph: ('a': ('b', 'label5'),('e', 'label1'))('b': ('e', 'label1'))('c': ('e', 'label2'))('d': ('e', 'label3'))('e': ('d', 'label4'))>")
        a.remove_node('e')
        self.assertEqual(repr(a), "<MultiGraph: ('a': ('b', 'label5'))('b': )('c': )('d': )>")
        b = MultiGraph(['a','b','c','d'])
        b.add_edge('a','b','label5')
        self.assertEqual(repr(b), "<MultiGraph: ('a': ('b', 'label5'))('b': )('c': )('d': )>")
        self.assertEqual(repr(a), repr(b))
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
        self.assertEqual(self.r_1, self.r_1i) #, "not equal to similar")
        self.assertNotEqual(self.r_3, self.r_4) #, "equal to different")
        
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
