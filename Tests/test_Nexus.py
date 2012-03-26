import os.path
import unittest
import tempfile
import cStringIO
import sys

from Bio.Nexus import Nexus, Trees


class NexusTest1(unittest.TestCase):
    def setUp(self):
        self.testfile_dir = "Nexus"
        self.handle = open(os.path.join(self.testfile_dir, 
            "test_Nexus_input.nex"))

    def tearDown(self):
        self.handle.close()

    def test_WriteToFileName(self):
        """Test writing to a given filename."""
        filename = "Nexus/test_temp.nex"
        if os.path.isfile(filename):
            os.remove(filename)
        n = Nexus.Nexus(self.handle)
        n.write_nexus_data(filename)
        self.assertTrue(os.path.isfile(filename))
        os.remove(filename)

    def test_NexusTest1(self):
        """Test Nexus module"""
        # check data of main nexus file
        n=Nexus.Nexus(self.handle)
        self.assertEqual(os.path.normpath(n.filename),
                         os.path.normpath("Nexus/test_Nexus_input.nex"))
        self.assertEqual(n.ntax, 9)
        self.assertEqual(n.nchar, 48)
        self.assertEqual(n.datatype, "dna")
        self.assertEqual(n.interleave, True)
        self.assertEqual(n.missing, "?")
        self.assertEqual(n.gap, "-")
        self.assertEqual(n.taxlabels, ["t1",
                                       "t2 the name",
                                       "isn'that [a] strange name?",
                                       "one should be punished, for (that)!",
                                       "t5",
                                       "t6",
                                       "t7",
                                       "t8",
                                       "t9"])
        self.assertEqual(n.charlabels, {0: 'a',
                                        1: 'b',
                                        2: 'c',
                                        4: 'f',
                                        9: 'A',
                                        10: 'B',
                                        22: 'x',
                                        23: "y",
                                        29: "1,2,3 can't decide for a name?!",
                                        47: "final"})
        self.assertEqual(n.charsets, {
            "big":        [0, 2, 4, 6],
            "bigchunk":   [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46],
            "byname":     [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29],
            "c1":         [0, 1, 2, 3, 4, 5, 6, 7],
            "c2":         [8, 9, 10, 11, 12, 13, 14, 15],
            "c3":         [16, 17, 18, 19, 20, 21, 22, 23],
            "firsthalf":  [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23],
            "mix":        [0, 1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46],
            "mux":        [0, 1, 4, 7, 8, 10, 13, 16, 17, 18, 19, 20, 21, 22, 23, 25, 28, 31, 34, 37, 40, 43, 46],
            "pos1":       [0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45],
            "pos2":       [1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46],
            "pos3":       [2, 5, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 44, 47],
            "secondhalf": [24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47],
                                     })
        self.assertEqual(n.taxsets, {
            "normal":    ["isn'that [a] strange name?",
                          'one should be punished, for (that)!',
                          't1',
                          't5',
                          't6',
                          't8'],
            "reference": ["isn'that [a] strange name?",
                          'one should be punished, for (that)!',
                          't1',
                          't2 the name',
                          't5',
                          't6'],
            "tbyname1": ["isn'that [a] strange name?",
                          'one should be punished, for (that)!',
                          't1',
                          't2 the name',
                          't5',
                          't6'],
            "tbyname2": ["isn'that [a] strange name?",
                          'one should be punished, for (that)!',
                          't2 the name',
                          't5',
                          't6',
                          't7'],
            "tbyname3": ['t1',
                          't2 the name'],
                                    })
        self.assertEqual(len(n.charpartitions), 2)
        self.assertTrue('codons' in n.charpartitions)
        self.assertTrue('part' in n.charpartitions)
        self.assertEqual(n.charpartitions['codons'],
            {'a': [0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45],
             'b': [1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46],
             'c': [2, 5, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 44, 47],
            })
        self.assertEqual(n.charpartitions['part'],
            {"one":   [0, 1, 2, 3, 4, 5, 6, 7],
             "three": [16, 17, 18, 19, 20, 21, 22, 23],
             "two":   [8, 9, 10, 11, 12, 13, 14, 15],
            })
        self.assertEqual(n.taxpartitions.keys(), ['taxpart'])
        self.assertEqual(n.taxpartitions['taxpart'], 
            {"badnames":  ["isn'that [a] strange name?",
                           'one should be punished, for (that)!',
                           't2 the name'],
             "goodnames": ['t1', 't5', 't6', 't7', 't8', 't9'],
            })

        # now we check excluding characters, deleting taxa,
        # and exporting adjusted sets
        f1=tempfile.NamedTemporaryFile("w+")
        n.write_nexus_data(f1,
                           delete=['t1','t7'],
                           exclude=n.invert(n.charsets['big']))
        f1.seek(0)
        nf1=Nexus.Nexus(f1)
        self.assertEqual(os.path.normpath(nf1.filename),
                         os.path.normpath(f1.name))
        self.assertEqual(nf1.ntax, 7)
        self.assertEqual(nf1.nchar, 4)
        self.assertEqual(nf1.datatype, "dna")
        self.assertEqual(nf1.interleave, False)
        self.assertEqual(nf1.missing, "?")
        self.assertEqual(nf1.gap, "-")
        self.assertEqual(nf1.taxlabels, ["t2 the name",
                                         "isn'that [a] strange name?",
                                         "one should be punished, for (that)!",
                                         "t5",
                                         "t6",
                                         "t8",
                                         "t9"])
        self.assertEqual(nf1.charlabels, {0: 'a', 1: 'c', 2: 'f'})
        self.assertEqual(nf1.charsets, {'big':       [0, 1, 2, 3],
                                        'bigchunk':  [1, 2, 3],
                                        'byname':    [0, 2, 3],
                                        'c1':        [0, 1, 2, 3],
                                        'firsthalf': [0, 1, 2, 3],
                                        'mix':       [0, 2],
                                        'mux':       [0, 2],
                                        'pos1':      [0, 3],
                                        'pos2':      [2],
                                        'pos3':      [1],
                                       })
        self.assertEqual(nf1.taxsets, {
            'normal':    ["isn'that [a] strange name?",
                          'one should be punished, for (that)!',
                          't5',
                          't6',
                          't8'],
            'reference': ["isn'that [a] strange name?",
                          'one should be punished, for (that)!',
                          't2 the name',
                          't5',
                          't6'],
            'tbyname1':  ["isn'that [a] strange name?",
                          'one should be punished, for (that)!',
                          't2 the name',
                          't5',
                          't6'],
            'tbyname2':  ["isn'that [a] strange name?",
                          'one should be punished, for (that)!',
                          't2 the name',
                          't5',
                          't6'],
            'tbyname3': ['t2 the name'],
                                       })
        self.assertEqual(len(nf1.charpartitions), 2)
        self.assertTrue('codons' in nf1.charpartitions)
        self.assertTrue('part' in nf1.charpartitions)
        self.assertEqual(nf1.charpartitions['codons'], {'a': [0, 3],
                                                        'b': [2],
                                                        'c': [1]})
        self.assertEqual(nf1.charpartitions['part'], {'one': [0, 1, 2, 3]})

        self.assertEqual(nf1.taxpartitions.keys(), ['taxpart'])
        self.assertEqual(nf1.taxpartitions['taxpart'],
            {"badnames":  ["isn'that [a] strange name?",
                           'one should be punished, for (that)!',
                           't2 the name'],
             "goodnames": ['t5', 't6', 't8', 't9'],
            })

        f2=tempfile.NamedTemporaryFile("w+")
        n.write_nexus_data(f2,
                           delete=['t2_the_name'],
                           exclude=range(3,40,4))
        f2.seek(0)
        nf2=Nexus.Nexus(f2)
        self.assertEqual(os.path.normpath(nf2.filename),
                         os.path.normpath(f2.name))
        self.assertEqual(nf2.ntax, 9)
        self.assertEqual(nf2.nchar, 38)
        self.assertEqual(nf2.datatype, "dna")
        self.assertEqual(nf2.interleave, False)
        self.assertEqual(nf2.missing, "?")
        self.assertEqual(nf2.gap, "-")
        self.assertEqual(nf2.taxlabels, ["t1", 
                                         "t2 the name", 
                                         "isn'that [a] strange name?", 
                                         "one should be punished, for (that)!", 
                                         "t5", 
                                         "t6", 
                                         "t7", 
                                         "t8", 
                                         "t9"])
        self.assertEqual(nf2.charlabels, {0: "a",
                                          1: "b",
                                          2: "c",
                                          3: "f",
                                          7: "A",
                                          8: "B",
                                          17: "x",
                                          22: "1,2,3 can't decide for a name?!",
                                          37: "final"})
        self.assertEqual(nf2.charsets,
            {"big":        [0, 2, 3, 5],
             "bigchunk":   [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36],
             "byname":     [0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22],
             "c1":         [0, 1, 2, 3, 4, 5],
             "c2":         [6, 7, 8, 9, 10, 11],
             "c3":         [12, 13, 14, 15, 16, 17],
             "firsthalf":  [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17],
             "mix":        [0, 1, 3, 8, 10, 12, 17, 19, 21, 26, 28, 30, 33, 36],
             "mux":        [0, 1, 3, 6, 8, 10, 12, 13, 14, 15, 16, 17, 19, 21, 26, 28, 30, 33, 36],
             "pos1":       [0, 5, 7, 9, 14, 16, 18, 23, 25, 27, 32, 35],
             "pos2":       [1, 3, 8, 10, 12, 17, 19, 21, 26, 28, 30, 33, 36],
             "pos3":       [2, 4, 6, 11, 13, 15, 20, 22, 24, 29, 31, 34, 37],
             "secondhalf": [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37],
            })

        self.assertEqual(nf2.taxsets,
            {"normal":    ["isn'that [a] strange name?",
                           'one should be punished, for (that)!',
                           't1',
                           't5',
                           't6',
                           't8'],
             "reference": ["isn'that [a] strange name?",
                           'one should be punished, for (that)!',
                           't1',
                           't2 the name',
                           't5',
                           't6'],
             "tbyname1":  ["isn'that [a] strange name?",
                           'one should be punished, for (that)!',
                           't1',
                           't2 the name',
                           't5',
                           't6'],
             "tbyname2":  ["isn'that [a] strange name?",
                           'one should be punished, for (that)!',
                           't2 the name',
                           't5',
                           't6',
                           't7'],
             "tbyname3":  ['t1',
                           't2 the name']})
        self.assertEqual(len(nf2.charpartitions), 2)
        self.assertTrue('codons' in nf2.charpartitions)
        self.assertTrue('part' in nf2.charpartitions)
        self.assertEqual(nf2.charpartitions['codons'],
            {"a": [0, 5, 7, 9, 14, 16, 18, 23, 25, 27, 32, 35],
             "b": [1, 3, 8, 10, 12, 17, 19, 21, 26, 28, 30, 33, 36],
             "c": [2, 4, 6, 11, 13, 15, 20, 22, 24, 29, 31, 34, 37],
            })
        self.assertEqual(nf2.charpartitions['part'],
            {"one":   [0, 1, 2, 3, 4, 5],
             "three": [12, 13, 14, 15, 16, 17],
             "two":   [6, 7, 8, 9, 10, 11],
            })
        self.assertEqual(nf2.taxpartitions.keys(), ['taxpart'])
        self.assertEqual(nf2.taxpartitions['taxpart'],
            {"badnames":  ["isn'that [a] strange name?",
                           'one should be punished, for (that)!',
                           't2 the name'],
             "goodnames": ['t1', 't5', 't6', 't7', 't8', 't9'],
            })
        # check the stepmatrix
        self.assertEqual(n.weighted_stepmatrix(name='matrix_test'),
        """\
usertype matrix_test stepmatrix=5
        A        C        G        T        -
[A]     .       2.40     2.57     2.43     2.43     
[C]    2.40      .       2.28     2.12     2.14     
[G]    2.57     2.28      .       2.31     2.31     
[T]    2.43     2.12     2.31      .       2.14     
[-]    2.43     2.14     2.31     2.14      .       
;
""")


    def test_TreeTest1(self):
        """Test Tree module."""
        n=Nexus.Nexus(self.handle)
        t3=n.trees[2]
        t2=n.trees[2]
        t3.root_with_outgroup(['t1','t5'])
        self.assertEqual(str(t3), "tree tree1 = (((((('one should be punished, for (that)!','isn''that [a] strange name?'),'t2 the name'),t8,t9),t6),t7),(t5,t1));")
        self.assertEqual(t3.is_monophyletic(['t8','t9','t6','t7']), -1)
        self.assertEqual(t3.is_monophyletic(['t1','t5']), 13)
        t3.split(parent_id=t3.search_taxon('t9'))
        stdout = sys.stdout
        try:
            sys.stdout = cStringIO.StringIO()
            t3.display()
            if sys.version_info[0] == 3:
                output = sys.stdout.getvalue()
            else:
                sys.stdout.reset()
                output = sys.stdout.read()
        finally:
            sys.stdout = stdout 
        expected = """\
  #                            taxon            prev            succ    brlen blen (sum)  support              comment
  1    'isn''that [a] strange name?'               2              []   100.00     119.84    10.00                    -
  2                                -               4          [3, 1]     0.40      19.84     0.30                    -
  3 'one should be punished, for (that)!'               2              []     0.50      20.34        -                    -
  4                                -               6          [2, 5]     4.00      19.44     3.00                    -
  5                    't2 the name'               4              []     0.30      19.74        -                    -
  6                                -               9       [4, 7, 8]     2.00      15.44     1.00                    -
  7                               t8               6              []     1.20      16.64        -                    -
  8                               t9               6        [17, 18]     3.40      18.84        -                    -
  9                                -              11         [6, 10]     0.44      13.44    33.00                    -
 10                               t6               9              []     1.00      14.44        -                    -
 11                                -              16         [9, 12]    13.00      13.00    12.00                    -
 12                               t7              11              []    99.90     112.90        -                    -
 13                                -              16        [14, 15]     0.00       0.00     0.00                    -
 14                               t5              13              []    99.00      99.00        -                    -
 15                               t1              13              []     0.98       0.98        -                    -
 16                                -            None        [11, 13]     0.00       0.00        -                    -
 17                              t90               8              []     1.00      19.84        -                    -
 18                              t91               8              []     1.00      19.84        -                    -

Root:  16
"""
        self.assertEqual(len(output.split("\n")), len(expected.split("\n")))
        for l1, l2 in zip(output.split("\n"), expected.split("\n")):
            self.assertEqual(l1, l2)
        self.assertEqual(output, expected)
        self.assertEqual(t3.is_compatible(t2,threshold=0.3), [])

    def test_internal_node_labels(self):
        """Handle text labels on internal nodes.
        """
        ts1b = "(Cephalotaxus:125.000000,(Taxus:100.000000,Torreya:100.000000)"\
               "TT1:25.000000)Taxaceae:90.000000;"
        tree = Trees.Tree(ts1b)
        assert self._get_flat_nodes(tree) == [('Taxaceae', 90.0, None, None),
                ('Cephalotaxus', 125.0, None, None), ('TT1', 25.0, None, None),
                ('Taxus', 100.0, None, None), ('Torreya', 100.0, None, None)]

        ts1c = "(Cephalotaxus:125.000000,(Taxus:100.000000,Torreya:100.000000)"\
                "25.000000)90.000000;"
        tree = Trees.Tree(ts1c)
        assert self._get_flat_nodes(tree) == [(None, 90.0, None, None),
                ('Cephalotaxus', 125.0, None, None), (None, 25.0, None, None),
                ('Taxus', 100.0, None, None), ('Torreya', 100.0, None, None)]

        ts2 = "(((t9:0.385832, (t8:0.445135,t4:0.41401)C:0.024032)B:0.041436,"\
          "t6:0.392496)A:0.0291131, t2:0.497673, ((t0:0.301171,"\
          "t7:0.482152)E:0.0268148, ((t5:0.0984167,t3:0.488578)G:0.0349662,"\
          "t1:0.130208)F:0.0318288)D:0.0273876);"
        tree = Trees.Tree(ts2)

        large_ex_handle = open(os.path.join(self.testfile_dir,
            "int_node_labels.nwk"))
        tree = Trees.Tree(large_ex_handle.read())
        large_ex_handle.close()

    def _get_flat_nodes(self, tree):
        cur_nodes = [tree.node(tree.root)]
        nodedata = []
        while len(cur_nodes) > 0:
            new_nodes = []
            for cur_node in cur_nodes:
                nodedata.append((cur_node.data.taxon,
                    cur_node.data.branchlength, cur_node.data.support,
                    cur_node.data.comment))
                new_nodes.extend([tree.node(nid) for nid in
                    cur_node.get_succ()])
            cur_nodes = new_nodes
        return nodedata

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
