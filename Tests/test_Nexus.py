# Copyright 2005 by Iddo Friedberg.  All rights reserved.
# Revisions copyright 2006-2013,2017 by Peter Cock. All rights reserved.
# Revisions copyright 2008 by Frank Kauff. All rights reserved.
# Revisions copyright 2009 by Michiel de Hoon. All rights reserved.
# Revisions copyright 2015 by Joe Cora. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
from __future__ import print_function

import os.path
import unittest
import tempfile
import sys
from Bio._py3k import StringIO
from Bio._py3k import range
from Bio.Align import MultipleSeqAlignment
from Bio.AlignIO.NexusIO import NexusIterator, NexusWriter
from Bio.SeqRecord import SeqRecord
from Bio.Nexus import Nexus, Trees
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import ambiguous_dna
from Bio import SeqIO


class OldSelfTests(unittest.TestCase):
    """Test cases originally in Nexus.py via __main__"""

    def test_trees_and_taxa_block(self):
        """Basic tree file with TREES and TAXA block"""
        nexus1 = Nexus.Nexus()
        nexus1.read('Nexus/bats.nex')

    def test_data_and_codons_block(self):
        """Simple sequence data file with DATA and CODONS block"""
        nexus2 = Nexus.Nexus()
        nexus2.read('Nexus/codonposset.nex')

    def test_data_sets_trees_unknown_block(self):
        """Sequence data file with DATA, SETS, TREES and an unknown block"""
        nexus3 = Nexus.Nexus()
        nexus3.read('Nexus/test_Nexus_input.nex')

    def test_taxa_and_characters_block(self):
        """Taxa and characters multi-state block"""
        nexus4 = Nexus.Nexus()
        nexus4.read('Nexus/vSysLab_Ganaspidium_multistate.nex')

    def test_taxa_and_characters_with_many_codings_one_without_state(self):
        """Taxa and chr blocks, over 9 codings, 1 character without states"""
        nexus5 = Nexus.Nexus()
        nexus5.read('Nexus/vSysLab_Heptascelio_no-states_10+chars.nex')

    def test_taxa_and_characters_with_many_codings_two_without_state(self):
        """Taxa and chr blocks, over 9 codings, 2 character without states"""
        nexus6 = Nexus.Nexus()
        # TODO: Implement continuous datatype:
        # Bio.Nexus.Nexus.NexusError: Unsupported datatype: continuous
        self.assertRaises(Nexus.NexusError,
                          nexus6.read,
                          'Nexus/vSysLab_Oreiscelio_discrete+continuous.nex')


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

    def test_write_with_dups(self):
        # see issue: biopython/Bio/Nexus/Nexus.py _unique_label() eval error #633
        records = [SeqRecord(Seq("ATGCTGCTGAT", alphabet=ambiguous_dna), id="foo") for _ in range(4)]
        out_file = StringIO()
        self.assertEqual(4, SeqIO.write(records, out_file, "nexus"))

    def test_NexusTest1(self):
        """Test Nexus module"""
        # check data of main nexus file
        n = Nexus.Nexus(self.handle)
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
        self.assertEqual(n.charsets,
                         {"big": [0, 2, 4, 6],
                          "bigchunk": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46],
                          "byname": [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29],
                          "c1": [0, 1, 2, 3, 4, 5, 6, 7],
                          "c2": [8, 9, 10, 11, 12, 13, 14, 15],
                          "c3": [16, 17, 18, 19, 20, 21, 22, 23],
                          "firsthalf": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23],
                          "mix": [0, 1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46],
                          "mux": [0, 1, 4, 7, 8, 10, 13, 16, 17, 18, 19, 20, 21, 22, 23, 25, 28, 31, 34, 37, 40, 43, 46],
                          "pos1": [0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45],
                          "pos2": [1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46],
                          "pos3": [2, 5, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 44, 47],
                          "secondhalf": [24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47],
                          })
        self.assertEqual(n.taxsets,
                         {"normal": ["isn'that [a] strange name?",
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
        self.assertIn('codons', n.charpartitions)
        self.assertIn('part', n.charpartitions)
        self.assertEqual(n.charpartitions['codons'],
                         {'a': [0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45],
                          'b': [1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46],
                          'c': [2, 5, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 44, 47],
                          })
        self.assertEqual(n.charpartitions['part'],
                         {"one": [0, 1, 2, 3, 4, 5, 6, 7],
                          "three": [16, 17, 18, 19, 20, 21, 22, 23],
                          "two": [8, 9, 10, 11, 12, 13, 14, 15],
                          })
        self.assertEqual(list(n.taxpartitions), ['taxpart'])
        self.assertEqual(n.taxpartitions['taxpart'],
                         {"badnames": ["isn'that [a] strange name?",
                                       'one should be punished, for (that)!',
                                       't2 the name'],
                          "goodnames": ['t1', 't5', 't6', 't7', 't8', 't9'],
                          })

        # now we check excluding characters, deleting taxa,
        # and exporting adjusted sets
        f1 = tempfile.NamedTemporaryFile("w+")
        n.write_nexus_data(f1,
                           delete=['t1', 't7'],
                           exclude=n.invert(n.charsets['big']))
        f1.seek(0)
        nf1 = Nexus.Nexus(f1)
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
        self.assertEqual(nf1.charsets, {'big': [0, 1, 2, 3],
                                        'bigchunk': [1, 2, 3],
                                        'byname': [0, 2, 3],
                                        'c1': [0, 1, 2, 3],
                                        'firsthalf': [0, 1, 2, 3],
                                        'mix': [0, 2],
                                        'mux': [0, 2],
                                        'pos1': [0, 3],
                                        'pos2': [2],
                                        'pos3': [1],
                                        })
        self.assertEqual(nf1.taxsets,
                         {'normal': ["isn'that [a] strange name?",
                                     'one should be punished, for (that)!',
                                     't5',
                                     't6',
                                     't8'],
                          'reference': ["isn'that [a] strange name?",
                                        'one should be punished, for (that)!',
                                        't2 the name',
                                        't5',
                                        't6'],
                          'tbyname1': ["isn'that [a] strange name?",
                                       'one should be punished, for (that)!',
                                       't2 the name',
                                       't5',
                                       't6'],
                          'tbyname2': ["isn'that [a] strange name?",
                                       'one should be punished, for (that)!',
                                       't2 the name',
                                       't5',
                                       't6'],
                          'tbyname3': ['t2 the name'],
                          })
        self.assertEqual(len(nf1.charpartitions), 2)
        self.assertIn('codons', nf1.charpartitions)
        self.assertIn('part', nf1.charpartitions)
        self.assertEqual(nf1.charpartitions['codons'], {'a': [0, 3],
                                                        'b': [2],
                                                        'c': [1]})
        self.assertEqual(nf1.charpartitions['part'], {'one': [0, 1, 2, 3]})

        self.assertEqual(list(nf1.taxpartitions), ['taxpart'])
        self.assertEqual(nf1.taxpartitions['taxpart'],
                         {"badnames": ["isn'that [a] strange name?",
                                       'one should be punished, for (that)!',
                                       't2 the name'],
                          "goodnames": ['t5', 't6', 't8', 't9'],
                          })

        f2 = tempfile.NamedTemporaryFile("w+")
        n.write_nexus_data(f2,
                           delete=['t2_the_name'],
                           exclude=list(range(3, 40, 4)))
        f2.seek(0)
        nf2 = Nexus.Nexus(f2)
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
                         {"big": [0, 2, 3, 5],
                          "bigchunk": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36],
                          "byname": [0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22],
                          "c1": [0, 1, 2, 3, 4, 5],
                          "c2": [6, 7, 8, 9, 10, 11],
                          "c3": [12, 13, 14, 15, 16, 17],
                          "firsthalf": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17],
                          "mix": [0, 1, 3, 8, 10, 12, 17, 19, 21, 26, 28, 30, 33, 36],
                          "mux": [0, 1, 3, 6, 8, 10, 12, 13, 14, 15, 16, 17, 19, 21, 26, 28, 30, 33, 36],
                          "pos1": [0, 5, 7, 9, 14, 16, 18, 23, 25, 27, 32, 35],
                          "pos2": [1, 3, 8, 10, 12, 17, 19, 21, 26, 28, 30, 33, 36],
                          "pos3": [2, 4, 6, 11, 13, 15, 20, 22, 24, 29, 31, 34, 37],
                          "secondhalf": [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37],
                          })

        self.assertEqual(nf2.taxsets,
                         {"normal": ["isn'that [a] strange name?",
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
                                       't2 the name']})
        self.assertEqual(len(nf2.charpartitions), 2)
        self.assertIn('codons', nf2.charpartitions)
        self.assertIn('part', nf2.charpartitions)
        self.assertEqual(nf2.charpartitions['codons'],
                         {"a": [0, 5, 7, 9, 14, 16, 18, 23, 25, 27, 32, 35],
                          "b": [1, 3, 8, 10, 12, 17, 19, 21, 26, 28, 30, 33, 36],
                          "c": [2, 4, 6, 11, 13, 15, 20, 22, 24, 29, 31, 34, 37],
                          })
        self.assertEqual(nf2.charpartitions['part'],
                         {"one": [0, 1, 2, 3, 4, 5],
                          "three": [12, 13, 14, 15, 16, 17],
                          "two": [6, 7, 8, 9, 10, 11],
                          })
        self.assertEqual(list(nf2.taxpartitions), ['taxpart'])
        self.assertEqual(nf2.taxpartitions['taxpart'],
                         {"badnames": ["isn'that [a] strange name?",
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
""")  # noqa : W291

    def test_write_alignment(self):
        # Default causes no interleave (columns <= 1000)
        records = [SeqRecord(Seq("ATGCTGCTGA" * 90, alphabet=ambiguous_dna), id=_id) for _id in ["foo", "bar", "baz"]]
        a = MultipleSeqAlignment(records, alphabet=ambiguous_dna)

        handle = StringIO()
        NexusWriter(handle).write_alignment(a)
        handle.seek(0)
        data = handle.read()
        self.assertIn("ATGCTGCTGA" * 90, data)

        # Default causes interleave (columns > 1000)
        records = [SeqRecord(Seq("ATGCTGCTGA" * 110, alphabet=ambiguous_dna), id=_id) for _id in ["foo", "bar", "baz"]]
        a = MultipleSeqAlignment(records, alphabet=ambiguous_dna)
        handle = StringIO()
        NexusWriter(handle).write_alignment(a)
        handle.seek(0)
        data = handle.read()
        self.assertNotIn("ATGCTGCTGA" * 90, data)
        self.assertIn("ATGCTGCTGA" * 7, data)

        # Override interleave: True
        records = [SeqRecord(Seq("ATGCTGCTGA" * 9, alphabet=ambiguous_dna), id=_id) for _id in ["foo", "bar", "baz"]]
        a = MultipleSeqAlignment(records, alphabet=ambiguous_dna)
        handle = StringIO()
        NexusWriter(handle).write_alignment(a, interleave=True)
        handle.seek(0)
        data = handle.read()
        self.assertNotIn("ATGCTGCTGA" * 9, data)
        self.assertIn("ATGCTGCTGA" * 7, data)

        # Override interleave: False
        records = [SeqRecord(Seq("ATGCTGCTGA" * 110, alphabet=ambiguous_dna), id=_id) for _id in ["foo", "bar", "baz"]]
        a = MultipleSeqAlignment(records, alphabet=ambiguous_dna)
        handle = StringIO()
        NexusWriter(handle).write_alignment(a, interleave=False)
        handle.seek(0)
        data = handle.read()
        self.assertIn("ATGCTGCTGA" * 110, data)

    def test_TreeTest1(self):
        """Test Tree module."""
        n = Nexus.Nexus(self.handle)
        t3 = n.trees[2]
        t2 = n.trees[2]
        t3.root_with_outgroup(['t1', 't5'])
        self.assertEqual(str(t3), "tree tree1 = (((((('one should be punished, for (that)!','isn''that [a] strange name?'),'t2 the name'),t8,t9),t6),t7),(t5,t1));")
        self.assertEqual(t3.is_monophyletic(['t8', 't9', 't6', 't7']), -1)
        self.assertEqual(t3.is_monophyletic(['t1', 't5']), 13)
        t3.split(parent_id=t3.search_taxon('t9'))
        stdout = sys.stdout
        try:
            sys.stdout = StringIO()
            t3.display()
            output = sys.stdout.getvalue()
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
        self.assertEqual(t3.is_compatible(t2, threshold=0.3), [])

    def test_TreeTest2(self):
        """Handle text labels on internal nodes.
        """
        ts1b = "(Cephalotaxus:125.000000,(Taxus:100.000000,Torreya:100.000000)"\
               "TT1:25.000000)Taxaceae:90.000000;"
        tree = Trees.Tree(ts1b)
        self.assertEqual(self._get_flat_nodes(tree), [('Taxaceae', 90.0, None, None),
                                                      ('Cephalotaxus', 125.0, None, None),
                                                      ('TT1', 25.0, None, None),
                                                      ('Taxus', 100.0, None, None),
                                                      ('Torreya', 100.0, None, None)])
        tree.prune('Torreya')
        self.assertEqual(tree.all_ids(), [0, 1, 3])
        ts1c = ("(Cephalotaxus:125.000000,(Taxus:100.000000,Torreya:100.000000)"
                "25.000000)90.000000;")
        tree = Trees.Tree(ts1c)
        self.assertEqual(self._get_flat_nodes(tree), [(None, 90.0, None, None),
                                                      ('Cephalotaxus', 125.0, None, None),
                                                      (None, 25.0, None, None),
                                                      ('Taxus', 100.0, None, None),
                                                      ('Torreya', 100.0, None, None)])
        self.assertFalse(tree.has_support())
        with self.assertRaises(Exception) as context:
            tree.randomize()
        self.assertTrue("Either numer of taxa or list of taxa must be specified." in str(context.exception))
        tree_rand = Trees.Tree(ts1c)
        tree_rand.randomize(ntax=4)
        self.assertEqual(sorted(tree_rand.get_taxa()), ['taxon1', 'taxon2',
                                                        'taxon3', 'taxon4'])
        tree.branchlength2support()
        tree.convert_absolute_support(2)
        self.assertEqual(self._get_flat_nodes(tree), [(None, 0.0, 90.0, None),
                                                      ('Cephalotaxus', 0.0, 62.5, None),
                                                      (None, 0.0, 12.5, None),
                                                      ('Taxus', 0.0, 50.0, None),
                                                      ('Torreya', 0.0, 50.0, None)])

        ts2 = ("(((t9:0.385832, (t8:0.445135,t4:0.41401)C:0.024032)B:0.041436,"
               "t6:0.392496)A:0.0291131, t2:0.497673, ((t0:0.301171,"
               "t7:0.482152)E:0.0268148, ((t5:0.0984167,t3:0.488578)G:0.0349662,"
               "t1:0.130208)F:0.0318288)D:0.0273876);")
        tree = Trees.Tree(ts2)
        tree.branchlength2support()
        supports = []
        for i in tree.all_ids():
            node = tree.node(i)
            data = node.get_data()
            supports.append(data.support)
        self.assertEqual(supports, [0.0, 0.0291131, 0.041436, 0.385832, 0.024032,
                                    0.445135, 0.41401, 0.392496, 0.497673,
                                    0.0273876, 0.0268148, 0.301171, 0.482152,
                                    0.0318288, 0.0349662, 0.0984167, 0.488578,
                                    0.130208])
        ts3 = ("(((B 9:0.385832, (C 8:0.445135, C4:0.41401)C:0.024032)B:0.041436,"
               "A 6:0.392496)A:0.0291131, t2:0.497673, ((E 0:0.301171,"
               "E 7:0.482152)E:0.0268148, ((G 5:0.0984167,G 3:0.488578)G:0.0349662,"
               "F 1:0.130208)F:0.0318288)D:0.0273876);")
        self.assertFalse(tree.is_identical(Trees.Tree(ts3)))
        tree = Trees.Tree(ts3)
        self.assertTrue(tree.is_bifurcating())
        self.assertTrue(tree.is_bifurcating(1))
        self.assertEqual([tree.distance(0, n) for n in tree.all_ids()], [0.0,
                                                                         0.0291131,
                                                                         0.0705491,
                                                                         0.4563811,
                                                                         0.0945811,
                                                                         0.5397161,
                                                                         0.5085911,
                                                                         0.4216091,
                                                                         0.497673,
                                                                         0.0273876,
                                                                         0.0542024,
                                                                         0.3553734,
                                                                         0.5363544,
                                                                         0.0592164,
                                                                         0.09418259999999999,
                                                                         0.1925993,
                                                                         0.5827606,
                                                                         0.1894244])

        subtree = tree.set_subtree(10)
        self.assertEqual(sorted(list(subtree)), ['E 0', 'E 7'])
        tree.collapse_genera()
        self.assertEqual(tree.all_ids(), [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 14, 17])

    def test_merge_with_support(self):
        """Test merge_with_support and consensus method."""
        ts1 = ("(((B 9:0.385832, (C 8:0.445135, C 4:0.41401)C:0.024032)B:0.041436,"
               "A 6:0.392496)A:0.0291131, t2:0.497673, ((E 0:0.301171,"
               "E 7:0.482152)E:0.0268148, ((G 5:0.0984167,G 3:0.488578)G:0.0349662,"
               "F 1:0.130208)F:0.0318288)D:0.0273876);")

        tbs1 = ("(((B 9:0.385832, (C 8:0.445135, C 4:0.41401)C:0.024032)B:0.041436,"
                "A 6:0.392496)A:0.0291131, t2:0.497673, ((G 5:0.0984167,"
                "G 3:0.488578)E:0.0268148, ((E 0:0.301171, E 7:0.482152)G:0.0349662,"
                "F 1:0.130208)F:0.0318288)D:0.0273876);")

        tbs2 = ("(((B 9:0.385832,A 6:0.392496 C:0.024032)B:0.041436, (C 8:0.445135,"
                "C 4:0.41401))A:0.0291131, t2:0.497673, ((E 0:0.301171, E 7:0.482152)"
                "E:0.0268148, ((G 5:0.0984167,G 3:0.488578)G:0.0349662,F 1:0.130208)"
                "F:0.0318288)D:0.0273876);")

        t1 = Trees.Tree(ts1)
        tb1 = Trees.Tree(tbs1)
        tb2 = Trees.Tree(tbs2)

        t1.branchlength2support()
        tb1.branchlength2support()
        tb2.branchlength2support()

        t1.merge_with_support(bstrees=[tb1, tb2], threshold=0.2)

        supports = []
        for i in t1.all_ids():
            node = t1.node(i)
            data = node.get_data()
            supports.append(data.support)
        self.assertTrue(supports, [0.0, 1.0, 0.04, 1.0, 0.5, 1.0, 1.0, 1.0,
                                   1.0, 1.0, 0.5, 1.0, 1.0, 0.5, 1.0, 1.0, 1.0,
                                   1.0])

    def test_large_newick(self):
        with open(os.path.join(self.testfile_dir, "int_node_labels.nwk")) as large_ex_handle:
            tree = Trees.Tree(large_ex_handle.read())

    def _get_flat_nodes(self, tree):
        cur_nodes = [tree.node(tree.root)]
        nodedata = []
        while len(cur_nodes) > 0:
            new_nodes = []
            for cur_node in cur_nodes:
                nodedata.append((cur_node.data.taxon,
                                 cur_node.data.branchlength,
                                 cur_node.data.support,
                                 cur_node.data.comment))
                new_nodes.extend([tree.node(nid) for nid in
                                  cur_node.get_succ()])
            cur_nodes = new_nodes
        return nodedata

    def test_NexusComments(self):
        """Test the ability to parse nexus comments at internal and leaf nodes
        """
        # A tree with simple comments throughout the tree.
        ts1b = "((12:0.13,19[&comment1]:0.13)[&comment2]:0.1,(20:0.171,11:0.171):0.13)[&comment3];"
        tree = Trees.Tree(ts1b)
        self.assertEqual(self._get_flat_nodes(tree), [(None, 0.0, None, '[&comment3]'),
                                                      (None, 0.1, None, '[&comment2]'),
                                                      (None, 0.13, None, None), ('12', 0.13, None, None),
                                                      ('19', 0.13, None, '[&comment1]'),
                                                      ('20', 0.171, None, None),
                                                      ('11', 0.171, None, None)])

        # A tree with more complex comments throughout the tree.
        # This is typical of the MCC trees produced by `treeannotator` in the beast-mcmc suite of phylogenetic tools
        # The key difference being tested here is the ability to parse internal node comments that include ','.
        ts1b = "(((9[&rate_range={1.3E-5,0.10958320752991428},height_95%_HPD={0.309132419999969,0.3091324199999691},length_range={3.513906814545109E-4,0.4381986285528381},height_median=0.309132419999969,length_95%_HPD={0.003011577063374571,0.08041621647998398}]:0.055354097721950546,5[&rate_range={1.3E-5,0.10958320752991428},height_95%_HPD={0.309132419999969,0.3091324199999691},length_range={3.865051168833178E-5,0.4391594442572986},height_median=0.309132419999969,length_95%_HPD={0.003011577063374571,0.08041621647998398}]:0.055354097721950546)[&height_95%_HPD={0.3110921040545068,0.38690865205576275},length_range={0.09675588357303178,0.4332959544380489},length_95%_HPD={0.16680375169879613,0.36500804261814374}]:0.20039426358269385)[&height_95%_HPD={0.5289500597932948,0.6973881165460601},length_range={0.02586430194846201,0.29509451958008265},length_95%_HPD={0.0840287249314221,0.2411078625957056}]:0.23042678598484334)[&height_95%_HPD={0.7527502510685965,0.821862094763501},height_median=0.8014438411766163,height=0.795965080422763,posterior=1.0,height_range={0.49863013698599995,0.821862094763501},length=0.0];"
        tree = Trees.Tree(ts1b)
        self.assertEqual(self._get_flat_nodes(tree),
                         [(None, 0.0, None, '[&height_95%_HPD={0.7527502510685965,0.821862094763501},height_median=0.8014438411766163,height=0.795965080422763,posterior=1.0,height_range={0.49863013698599995,0.821862094763501},length=0.0]'),
                          (None, 0.23042678598484334, None, '[&height_95%_HPD={0.5289500597932948,0.6973881165460601},length_range={0.02586430194846201,0.29509451958008265},length_95%_HPD={0.0840287249314221,0.2411078625957056}]'),
                          (None, 0.20039426358269385, None, '[&height_95%_HPD={0.3110921040545068,0.38690865205576275},length_range={0.09675588357303178,0.4332959544380489},length_95%_HPD={0.16680375169879613,0.36500804261814374}]'),
                          ('9', 0.055354097721950546, None, '[&rate_range={1.3E-5,0.10958320752991428},height_95%_HPD={0.309132419999969,0.3091324199999691},length_range={3.513906814545109E-4,0.4381986285528381},height_median=0.309132419999969,length_95%_HPD={0.003011577063374571,0.08041621647998398}]'),
                          ('5', 0.055354097721950546, None, '[&rate_range={1.3E-5,0.10958320752991428},height_95%_HPD={0.309132419999969,0.3091324199999691},length_range={3.865051168833178E-5,0.4391594442572986},height_median=0.309132419999969,length_95%_HPD={0.003011577063374571,0.08041621647998398}]')])


class TestSelf(unittest.TestCase):
    def test_repeated_names_no_taxa(self):
        print("Repeated names without a TAXA block")
        handle = StringIO("""#NEXUS
        [TITLE: NoName]
        begin data;
        dimensions ntax=4 nchar=50;
        format interleave datatype=protein   gap=- symbols="FSTNKEYVQMCLAWPHDRIG";
        matrix
        CYS1_DICDI          -----MKVIL LFVLAVFTVF VSS------- --------RG IPPEEQ----
        ALEU_HORVU          MAHARVLLLA LAVLATAAVA VASSSSFADS NPIRPVTDRA ASTLESAVLG
        CATH_HUMAN          ------MWAT LPLLCAGAWL LGV------- -PVCGAAELS VNSLEK----
        CYS1_DICDI          -----MKVIL LFVLAVFTVF VSS------- --------RG IPPEEQ---X
        ;
        end;
        """)  # noqa : W291
        for a in NexusIterator(handle):
            print(a)
            for r in a:
                print("%r %s %s" % (r.seq, r.name, r.id))
        print("Done")

    def test_repeated_names_with_taxa(self):

        print("Repeated names with a TAXA block")
        handle = StringIO("""#NEXUS
        [TITLE: NoName]
        begin taxa
        CYS1_DICDI
        ALEU_HORVU
        CATH_HUMAN
        CYS1_DICDI;
        end;
        begin data;
        dimensions ntax=4 nchar=50;
        format interleave datatype=protein   gap=- symbols="FSTNKEYVQMCLAWPHDRIG";
        matrix
        CYS1_DICDI          -----MKVIL LFVLAVFTVF VSS------- --------RG IPPEEQ----
        ALEU_HORVU          MAHARVLLLA LAVLATAAVA VASSSSFADS NPIRPVTDRA ASTLESAVLG
        CATH_HUMAN          ------MWAT LPLLCAGAWL LGV------- -PVCGAAELS VNSLEK----
        CYS1_DICDI          -----MKVIL LFVLAVFTVF VSS------- --------RG IPPEEQ---X
        ;
        end;
        """)  # noqa : W291
        for a in NexusIterator(handle):
            print(a)
            for r in a:
                print("%r %s %s" % (r.seq, r.name, r.id))
        print("Done")

    def test_empty_file_read(self):
        self.assertEqual([], list(NexusIterator(StringIO())))

    def test_multiple_output(self):
        records = [SeqRecord(Seq("ATGCTGCTGAT", alphabet=ambiguous_dna), id="foo"),
                   SeqRecord(Seq("ATGCTGCAGAT", alphabet=ambiguous_dna), id="bar"),
                   SeqRecord(Seq("ATGCTGCGGAT", alphabet=ambiguous_dna), id="baz")]
        a = MultipleSeqAlignment(records, alphabet=ambiguous_dna)

        handle = StringIO()
        NexusWriter(handle).write_file([a])
        handle.seek(0)
        data = handle.read()
        self.assertTrue(data.startswith("#NEXUS\nbegin data;\n"), data)
        self.assertTrue(data.endswith("end;\n"), data)

        handle = StringIO()
        try:
            NexusWriter(handle).write_file([a, a])
            assert False, "Should have rejected more than one alignment!"
        except ValueError:
            pass


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
