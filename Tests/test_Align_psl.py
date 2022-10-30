# Copyright 2022 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Align.psl module."""
import unittest
from io import StringIO


from Bio.Align import Alignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SimpleLocation, ExactPosition, CompoundLocation
from Bio import SeqIO
from Bio import Align


try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align.psl."
    ) from None


class TestAlign_dna_rna(unittest.TestCase):

    # The PSL file dna_rna.psl was generated using this command:
    # blat -mask=lower hg38.2bit rna.fa dna_rna.psl

    def setUp(self):
        data = {}
        records = SeqIO.parse("Blat/dna.fa", "fasta")
        for record in records:
            name, start_end = record.id.split(":")
            assert name == "chr3"
            start, end = start_end.split("-")
            start = int(start)
            end = int(end)
            sequence = str(record.seq)
            assert len(sequence) == end - start
            data[start] = sequence
        self.dna = data
        records = SeqIO.parse("Blat/rna.fa", "fasta")
        self.rna = {record.id: record.seq for record in records}

    def test_reading(self):
        """Test parsing dna_rna.psl."""
        path = "Blat/dna_rna.psl"
        alignments = Align.parse(path, "psl")
        self.assertEqual(alignments.metadata["psLayout version"], "3")
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 165)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 39)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 5407))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr3")
        self.assertEqual(alignment.query.id, "NR_111921.1")
        self.assertEqual(len(alignment.target.seq), 198295559)
        self.assertEqual(len(alignment.query.seq), 216)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array( [[48663767, 48663813, 48665640, 48665722, 48669098, 48669174],
                              [       0,        46,      46,      128,      128,      204]]),
                # fmt: on
            )
        )
        dna = Seq(self.dna, length=len(alignment.target.seq))
        alignment.target.seq = dna
        alignment.query.seq = self.rna[alignment.query.id]
        self.assertTrue(
            numpy.array_equal(
                alignment.substitutions,
                # fmt: off
# flake8: noqa
            numpy.array([[53.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0., 35.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0., 50.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  0., 27.,  0.,  0.,  0.,  0.],
                         [ 9.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  7.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0., 16.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  0.,  7.,  0.,  0.,  0.,  0.],
                        ])
            )
        )
        self.assertEqual(alignment.substitutions.alphabet, "ACGTacgt")
        matches = sum(
            alignment.substitutions[c, c] for c in alignment.substitutions.alphabet
        )
        repMatches = sum(
            alignment.substitutions[c, c.swapcase()]
            for c in alignment.substitutions.alphabet
        )
        self.assertEqual(matches, alignment.matches)
        self.assertEqual(repMatches, alignment.repMatches)
        self.assertEqual(
            format(alignment, "psl"),
            """\
165	0	39	0	0	0	2	5203	+	NR_111921.1	216	0	204	chr3	198295559	48663767	48669174	3	46,82,76,	0,46,128,	48663767,48665640,48669098,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 175)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 6)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 1711))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr3")
        self.assertEqual(alignment.query.id, "NR_046654.1")
        self.assertEqual(len(alignment.target.seq), 198295559)
        self.assertEqual(len(alignment.query.seq), 181)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[42530895, 42530958, 42532020, 42532095, 42532563, 42532606],
                             [     181,      118,      118,       43,       43,        0]])
                # fmt: on
            )
        )
        dna = Seq(self.dna, length=len(alignment.target))
        alignment.target.seq = dna
        alignment.query.seq = self.rna[alignment.query.id]
        self.assertTrue(
            numpy.array_equal(
                alignment.substitutions,
                # fmt: off
# flake8: noqa
            numpy.array([[36.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0., 40.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0., 57.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  0., 42.,  0.,  0.,  0.,  0.],
                         [ 2.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  1.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  3.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                        ])
            )
        )
        self.assertEqual(alignment.substitutions.alphabet, "ACGTacgt")
        matches = sum(
            alignment.substitutions[c, c] for c in alignment.substitutions.alphabet
        )
        repMatches = sum(
            alignment.substitutions[c, c.swapcase()]
            for c in alignment.substitutions.alphabet
        )
        self.assertEqual(matches, alignment.matches)
        self.assertEqual(repMatches, alignment.repMatches)
        self.assertEqual(
            format(alignment, "psl"),
            """\
175	0	6	0	0	0	2	1530	-	NR_046654.1	181	0	181	chr3	198295559	42530895	42532606	3	63,75,43,	0,63,138,	42530895,42532020,42532563,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 162)
        self.assertEqual(alignment.misMatches, 2)
        self.assertEqual(alignment.repMatches, 39)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 5409))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr3")
        self.assertEqual(alignment.query.id, "NR_111921.1_modified")
        self.assertEqual(len(alignment.target.seq), 198295559)
        self.assertEqual(len(alignment.query.seq), 220)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[48663767, 48663795, 48663796, 48663813, 48665640,
                              48665716, 48665716, 48665722, 48669098, 48669174],
                             [       3,       31,       31,       48,       48,
                                   124,      126,      132,      132,      208]
                            ])
                # fmt: on
            )
        )
        dna = Seq(self.dna, length=len(alignment.target))
        alignment.target.seq = dna
        alignment.query.seq = self.rna[alignment.query.id]
        self.assertTrue(
            numpy.array_equal(
                alignment.substitutions,
                # fmt: off
# flake8: noqa
            numpy.array([[53.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0., 34.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  2., 48.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  0., 27.,  0.,  0.,  0.,  0.],
                         [ 9.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  7.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0., 16.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  0.,  7.,  0.,  0.,  0.,  0.],
                        ]),
            )
        )
        self.assertEqual(alignment.substitutions.alphabet, "ACGTacgt")
        matches = sum(
            alignment.substitutions[c, c] for c in alignment.substitutions.alphabet
        )
        repMatches = sum(
            alignment.substitutions[c, c.swapcase()]
            for c in alignment.substitutions.alphabet
            if c != "X"
        )
        self.assertEqual(matches, alignment.matches)
        self.assertEqual(repMatches, alignment.repMatches)
        self.assertEqual(
            format(alignment, "psl"),
            """\
162	2	39	0	1	2	3	5204	+	NR_111921.1_modified	220	3	208	chr3	198295559	48663767	48669174	5	28,17,76,6,76,	3,31,48,126,132,	48663767,48663796,48665640,48665716,48669098,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 172)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 6)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 1714))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr3")
        self.assertEqual(alignment.query.id, "NR_046654.1_modified")
        self.assertEqual(len(alignment.target.seq), 198295559)
        self.assertEqual(len(alignment.query.seq), 190)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[42530895, 42530922, 42530922, 42530958, 42532020,
                              42532037, 42532039, 42532095, 42532563, 42532606],
                             [     185,      158,      155,      119,      119,
                                   102,      102,       46,       46,        3],
                            ])
                # fmt: on
            )
        )
        dna = Seq(self.dna, length=len(alignment.target))
        alignment.target.seq = dna
        alignment.query.seq = self.rna[alignment.query.id]
        self.assertTrue(
            numpy.array_equal(
                alignment.substitutions,
                # fmt: off
# flake8: noqa
            numpy.array([[34.,  0.,  0.,  1.,  0.,  0.,  0.,  0.],
                         [ 0., 40.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0., 57.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  0., 41.,  0.,  0.,  0.,  0.],
                         [ 2.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  1.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  3.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                        ]),
            )
        )
        self.assertEqual(alignment.substitutions.alphabet, "ACGTacgt")
        matches = sum(
            alignment.substitutions[c, c] for c in alignment.substitutions.alphabet
        )
        repMatches = sum(
            alignment.substitutions[c, c.swapcase()]
            for c in alignment.substitutions.alphabet
            if c != "X"
        )
        self.assertEqual(matches, alignment.matches)
        self.assertEqual(repMatches, alignment.repMatches)
        self.assertEqual(
            format(alignment, "psl"),
            """\
172	1	6	0	1	3	3	1532	-	NR_046654.1_modified	190	3	185	chr3	198295559	42530895	42532606	5	27,36,17,56,43,	5,35,71,88,144,	42530895,42530922,42532020,42532039,42532563,
""",
        )
        self.assertRaises(StopIteration, next, alignments)

    def test_writing(self):
        """Test writing the alignments in dna_rna.psl."""
        path = "Blat/dna_rna.psl"
        with open(path) as stream:
            original_data = stream.read()
        alignments = Align.parse(path, "psl")
        stream = StringIO()
        n = Align.write(alignments, stream, "psl")
        self.assertEqual(n, 4)
        stream.seek(0)
        written_data = stream.read()
        stream.close()
        self.assertEqual(original_data, written_data)
        # Try this again. This time, we first strip the matches, misMatches,
        # repMatches, and nCount attributes from each alignment, and insert the
        # appropriate sequence data in each alignment. The writer will then
        # recalculate the matches, misMatches, repMatches, and nCount values
        # from the sequence data and the alignment, and store those values in
        # the PSL file.
        alignments = []
        for alignment in Align.parse(path, "psl"):
            del alignment.matches
            del alignment.misMatches
            del alignment.repMatches
            del alignment.nCount
            dna = Seq(self.dna, length=len(alignment.target))
            alignment.target.seq = dna
            alignment.query.seq = self.rna[alignment.sequences[1].id]
            alignments.append(alignment)
        stream = StringIO()
        n = Align.write(alignments, stream, "psl", mask="lower")
        self.assertEqual(n, 4)
        stream.seek(0)
        written_data = stream.read()
        stream.close()
        self.assertEqual(original_data, written_data)


class TestAlign_dna(unittest.TestCase):

    queries = {
        record.id: str(record.seq)
        for record in SeqIO.parse("Blat/fasta_34.fa", "fasta")
    }

    def test_reading_psl_34_001(self):
        """Test parsing psl_34_001.psl and pslx_34_001.pslx."""
        self.check_reading_psl_34_001("psl")
        self.check_reading_psl_34_001("pslx")

    def check_reading_psl_34_001(self, fmt):
        """Check parsing psl_34_001.psl or pslx_34_001.pslx."""
        path = "Blat/%s_34_001.%s" % (fmt, fmt)
        alignments = Align.parse(path, "psl")
        self.assertEqual(alignments.metadata["psLayout version"], "3")
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 16)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 16))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 191154276)
        self.assertEqual(len(alignment.query.seq), 33)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[61646095:61646111], "aggtaaactgccttca"
            )
            self.assertEqual(
                alignment.query.seq[11:27], self.queries[alignment.query.id][11:27]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['a', 'g', 'g', 't', 'a', 'a', 'a', 'c', 't', 'g', 'c', 'c', 't',
              't', 'c', 'a'],
             ['a', 'g', 'g', 't', 'a', 'a', 'a', 'c', 't', 'g', 'c', 'c', 't',
              't', 'c', 'a']], dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[61646095, 61646111],
                             [      11,       27]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
16	0	0	0	0	0	0	0	+	hg18_dna	33	11	27	chr4	191154276	61646095	61646111	1	16,	11,	61646095,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 33)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 33))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 33)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[10271783:10271816],
                "atgagcttccaaggtaaactgccttcaagattc",
            )
            self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['a', 't', 'g', 'a', 'g', 'c', 't', 't', 'c', 'c', 'a', 'a', 'g',
              'g', 't', 'a', 'a', 'a', 'c', 't', 'g', 'c', 'c', 't', 't', 'c',
              'a', 'a', 'g', 'a', 't', 't', 'c'],
             ['a', 't', 'g', 'a', 'g', 'c', 't', 't', 'c', 'c', 'a', 'a', 'g',
              'g', 't', 'a', 'a', 'a', 'c', 't', 'g', 'c', 'c', 't', 't', 'c',
              'a', 'a', 'g', 'a', 't', 't', 'c']], dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[10271783, 10271816],
                             [       0,       33]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
33	0	0	0	0	0	0	0	+	hg18_dna	33	0	33	chr1	249250621	10271783	10271816	1	33,	0,	10271783,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 17)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 17))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 33)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[53575980:53575997], "aaggcagtttaccttgg"
            )
            self.assertEqual(
                alignment.query.seq[8:25], self.queries[alignment.query.id][8:25]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['a', 'a', 'g', 'g', 'c', 'a', 'g', 't', 't', 't', 'a', 'c', 'c',
              't', 't', 'g', 'g'],
             ['a', 'a', 'g', 'g', 'c', 'a', 'g', 't', 't', 't', 'a', 'c', 'c',
              't', 't', 'g', 'g']], dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[53575980, 53575997],
                             [      25,        8]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
17	0	0	0	0	0	0	0	-	hg18_dna	33	8	25	chr2	243199373	53575980	53575997	1	17,	8,	53575980,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 38)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr9")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 141213431)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[85737865:85737906],
                "acaaaggggctgggcgcagtggctcacgcctgtaatcccaa",
            )
            self.assertEqual(
                alignment.query.seq[9:50], self.queries[alignment.query.id][9:50]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['a', 'c', 'a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g',
              'g', 'c', 'g', 'c', 'a', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a',
              'c', 'g', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c', 'c',
              'a', 'a'],
             ['a', 'c', 'a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g',
              'g', 'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a',
              'c', 'a', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c', 'c',
              'a', 'a']], dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[85737865, 85737906],
                             [       9,       50]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
38	3	0	0	0	0	0	0	+	hg19_dna	50	9	50	chr9	141213431	85737865	85737906	1	41,	9,	85737865,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 41)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr8")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 146364022)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[95160479:95160520],
                "cacaaaggggctgggcgtggtggctcacacctgtaatccca",
            )
            self.assertEqual(
                alignment.query.seq[8:49], self.queries[alignment.query.id][8:49]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['c', 'a', 'c', 'a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g',
              'g', 'g', 'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c',
              'a', 'c', 'a', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c',
              'c', 'a'],
             ['c', 'a', 'c', 'a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g',
              'g', 'g', 'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c',
              'a', 'c', 'a', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c',
              'c', 'a']], dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[95160479, 95160520],
                             [       8,       49]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
41	0	0	0	0	0	0	0	+	hg19_dna	50	8	49	chr8	146364022	95160479	95160520	1	41,	8,	95160479,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 33)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 51304566)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[42144400:42144436],
                "aaaggggctgggcgtggtagctcatgcctgtaatcc",
            )
            self.assertEqual(
                alignment.query.seq[11:47], self.queries[alignment.query.id][11:47]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g', 'c',
              'g', 't', 'g', 'g', 't', 'a', 'g', 'c', 't', 'c', 'a', 't', 'g',
              'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c'],
             ['a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g', 'c',
              'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a', 'c', 'a',
              'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c']], dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[42144400, 42144436],
                             [      11,       47]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
33	3	0	0	0	0	0	0	+	hg19_dna	50	11	47	chr22	51304566	42144400	42144436	1	36,	11,	42144400,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 43)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 48))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(alignment.target.seq[183925984:183925990], "aaaaat")
            self.assertEqual(
                alignment.target.seq[183925990:183926028],
                "aaaggggctgggcgtggtggctcacgcctgtaatccca",
            )
            self.assertEqual(
                alignment.query.seq[1:7], self.queries[alignment.query.id][1:7]
            )
            self.assertEqual(
                alignment.query.seq[38:49], self.queries[alignment.query.id][38:49]
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[183925984, 183925990, 183925990, 183926028],
                             [        1,         7,        11,        49]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
43	1	0	0	1	4	0	0	+	hg19_dna	50	1	49	chr2	243199373	183925984	183926028	2	6,38,	1,11,	183925984,183925990,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 34)
        self.assertEqual(alignment.misMatches, 2)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 170))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[35483340:35483365], "caaaggggctgggcgtagtggctga"
            )
            self.assertEqual(alignment.target.seq[35483499:35483510], "cacctgtaatc")
            self.assertEqual(
                alignment.query.seq[10:46], self.queries[alignment.query.id][10:46]
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[35483340, 35483365, 35483499, 35483510],
                             [      10,       35,       35,       46]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
34	2	0	0	0	0	1	134	+	hg19_dna	50	10	46	chr19	59128983	35483340	35483510	2	25,11,	10,35,	35483340,35483499,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 39)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 78077248)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[23891310:23891349],
                "caaaggggctgggcgtggtggctcacacctgtaatccca",
            )
            self.assertEqual(
                alignment.query.seq[10:49], self.queries[alignment.query.id][10:49]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['c', 'a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g',
              'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a', 'c',
              'a', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c', 'c', 'a'],
             ['c', 'a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g',
              'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a', 'c',
              'a', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c', 'c', 'a']],
            dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[23891310, 23891349],
                             [      10,       49]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
39	0	0	0	0	0	0	0	+	hg19_dna	50	10	49	chr18	78077248	23891310	23891349	1	39,	10,	23891310,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 27)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 28))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 78077248)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[43252217:43252245], "ggcgtggtggctcacgcctgtaatccca"
            )
            self.assertEqual(
                alignment.query.seq[21:49], self.queries[alignment.query.id][21:49]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['g', 'g', 'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c',
              'a', 'c', 'g', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c',
              'c', 'a'],
             ['g', 'g', 'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c',
              'a', 'c', 'a', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c',
              'c', 'a']], dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[43252217, 43252245],
                             [      21,       49]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
27	1	0	0	0	0	0	0	+	hg19_dna	50	21	49	chr18	78077248	43252217	43252245	1	28,	21,	43252217,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 44)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 54))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 115169878)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(alignment.target.seq[52759147:52759154], "aaaaatt")
            self.assertEqual(
                alignment.target.seq[52759160:52759198],
                "aaaggggctgggcgtggtggctcacgcctgtaatccca",
            )
            self.assertEqual(
                alignment.query.seq[1:8], self.queries[alignment.query.id][1:8]
            )
            self.assertEqual(
                alignment.query.seq[38:49], self.queries[alignment.query.id][38:49]
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[52759147, 52759154, 52759160, 52759160, 52759198],
                             [       1,        8,        8,       11,       49]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
44	1	0	0	1	3	1	6	+	hg19_dna	50	1	49	chr13	115169878	52759147	52759198	2	7,38,	1,11,	52759147,52759160,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 50)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 50))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[1207056:1207106],
                "caaaaattcacaaaggggctgggcgtggtggctcacacctgtaatcccaa",
            )
            self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['c', 'a', 'a', 'a', 'a', 'a', 't', 't', 'c', 'a', 'c', 'a', 'a',
              'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g', 'c', 'g', 't',
              'g', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a', 'c', 'a', 'c', 'c',
              't', 'g', 't', 'a', 'a', 't', 'c', 'c', 'c', 'a', 'a'],
             ['c', 'a', 'a', 'a', 'a', 'a', 't', 't', 'c', 'a', 'c', 'a', 'a',
              'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g', 'c', 'g', 't',
              'g', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a', 'c', 'a', 'c', 'c',
              't', 'g', 't', 'a', 'a', 't', 'c', 'c', 'c', 'a', 'a']],
            dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[1207056, 1207106],
                             [      0,      50]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
50	0	0	0	0	0	0	0	+	hg19_dna	50	0	50	chr1	249250621	1207056	1207106	1	50,	0,	1207056,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 31)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[61700837:61700871],
                "aaaaatgaacaaaggggctgggcgcggtggctca",
            )
            self.assertEqual(
                alignment.query.seq[1:35], self.queries[alignment.query.id][1:35]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['a', 'a', 'a', 'a', 'a', 't', 'g', 'a', 'a', 'c', 'a', 'a', 'a',
              'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g', 'c', 'g', 'c', 'g',
              'g', 't', 'g', 'g', 'c', 't', 'c', 'a'],
             ['a', 'a', 'a', 'a', 'a', 't', 't', 'c', 'a', 'c', 'a', 'a', 'a',
              'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g', 'c', 'g', 't', 'g',
              'g', 't', 'g', 'g', 'c', 't', 'c', 'a']], dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[61700837, 61700871],
                             [       1,       35]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
31	3	0	0	0	0	0	0	+	hg19_dna	50	1	35	chr1	249250621	61700837	61700871	1	34,	1,	61700837,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 28)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 44))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 191154276)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(alignment.target.seq[37558157:37558167], "tgggattaca")
            self.assertEqual(
                alignment.target.seq[37558173:37558191], "accacgcccagccccttt"
            )
            self.assertEqual(
                alignment.query.seq[11:29], self.queries[alignment.query.id][11:29]
            )
            self.assertEqual(
                alignment.query.seq[39:49], self.queries[alignment.query.id][39:49]
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[37558157, 37558167, 37558173, 37558173, 37558191],
                             [      49,       39,       39,       29,       11]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
28	0	0	0	1	10	1	6	-	hg19_dna	50	11	49	chr4	191154276	37558157	37558191	2	10,18,	1,21,	37558157,37558173,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 35)
        self.assertEqual(alignment.misMatches, 2)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 37))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 51304566)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[48997405:48997442],
                "tgggattacaggcgggagccaccacgcccagcccctt",
            )
            self.assertEqual(
                alignment.query.seq[12:49], self.queries[alignment.query.id][12:49]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 'c',
              'g', 'g', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
              'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't'],
             ['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
              'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
              'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't']],
            dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[48997405, 48997442],
                             [      49,       12]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
35	2	0	0	0	0	0	0	-	hg19_dna	50	12	49	chr22	51304566	48997405	48997442	1	37,	1,	48997405,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 35)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[120641740:120641776],
                "tgggattacaggcgtgagccaccacgcccagcccct",
            )
            self.assertEqual(
                alignment.query.seq[13:49], self.queries[alignment.query.id][13:49]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 'c',
              'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
              'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't'],
             ['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
              'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
              'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't']], dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[120641740, 120641776],
                             [       49,        13]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
35	1	0	0	0	0	0	0	-	hg19_dna	50	13	49	chr2	243199373	120641740	120641776	1	36,	1,	120641740,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 39)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[54017130:54017169],
                "tgggattacaggtgtgagccaccacgcccagcccctttg",
            )
            self.assertEqual(
                alignment.query.seq[10:49], self.queries[alignment.query.id][10:49]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
              'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
              'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't', 't', 'g'],
             ['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
              'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
              'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't', 't', 'g']],
            dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[54017130, 54017169],
                             [      49,       10]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
39	0	0	0	0	0	0	0	-	hg19_dna	50	10	49	chr19	59128983	54017130	54017169	1	39,	1,	54017130,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 36)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[553742:553781],
                "tgggatgacaggggtgaggcaccacgcccagcccctttg",
            )
            self.assertEqual(
                alignment.query.seq[10:49], self.queries[alignment.query.id][10:49]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['t', 'g', 'g', 'g', 'a', 't', 'g', 'a', 'c', 'a', 'g', 'g', 'g',
              'g', 't', 'g', 'a', 'g', 'g', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
              'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't', 't', 'g'],
             ['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
              'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
              'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't', 't', 'g']],
            dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[553742, 553781],
                             [    49,     10]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
36	3	0	0	0	0	0	0	-	hg19_dna	50	10	49	chr19	59128983	553742	553781	1	39,	1,	553742,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 33)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 135534747)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[99388555:99388591],
                "tgggattataggcatgagccaccacgcccagcccct",
            )
            self.assertEqual(
                alignment.query.seq[13:49], self.queries[alignment.query.id][13:49]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 't', 'a', 'g', 'g', 'c',
              'a', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
              'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't'],
             ['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
              'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
              'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't']], dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[99388555, 99388591],
                             [      49,       13]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
33	3	0	0	0	0	0	0	-	hg19_dna	50	13	49	chr10	135534747	99388555	99388591	1	36,	1,	99388555,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 24)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 25))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 135534747)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[112178171:112178196], "tgagtcaccacgcccagcccctttg"
            )
            self.assertEqual(
                alignment.query.seq[10:35], self.queries[alignment.query.id][10:35]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['t', 'g', 'a', 'g', 't', 'c', 'a', 'c', 'c', 'a', 'c', 'g', 'c',
              'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't', 't', 'g'],
             ['t', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g', 'c',
              'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't', 't', 'g']],
            dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[112178171, 112178196],
                             [       35,        10]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
24	1	0	0	0	0	0	0	-	hg19_dna	50	10	35	chr10	135534747	112178171	112178196	1	25,	15,	112178171,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 35)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[39368490:39368526],
                "tgggattacaggcgtgagccaccacgcccagcccct",
            )
            self.assertEqual(
                alignment.query.seq[13:49], self.queries[alignment.query.id][13:49]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 'c',
              'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
              'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't'],
             ['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
              'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
              'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't']], dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[39368490, 39368526],
                             [      49,       13]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
35	1	0	0	0	0	0	0	-	hg19_dna	50	13	49	chr1	249250621	39368490	39368526	1	36,	1,	39368490,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 33)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[220325687:220325721],
                "ggattacaggcgtgagccaccacgcccagcccct",
            )
            self.assertEqual(
                alignment.query.seq[13:47], self.queries[alignment.query.id][13:47]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 'c', 'g', 't',
              'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g', 'c', 'c',
              'c', 'a', 'g', 'c', 'c', 'c', 'c', 't'],
             ['g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't', 'g', 't',
              'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g', 'c', 'c',
              'c', 'a', 'g', 'c', 'c', 'c', 'c', 't']], dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[220325687, 220325721],
                             [       47,        13]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
33	1	0	0	0	0	0	0	-	hg19_dna	50	13	47	chr1	249250621	220325687	220325721	1	34,	3,	220325687,
""",
        )
        self.assertRaises(StopIteration, next, alignments)

    def test_writing_psl_34_001(self):
        """Test writing the alignments in psl_34_001.psl."""
        path = "Blat/psl_34_001.psl"
        with open(path) as stream:
            original_data = stream.read()
        alignments = Align.parse(path, "psl")
        stream = StringIO()
        n = Align.write(alignments, stream, "psl")
        self.assertEqual(n, 22)
        stream.seek(0)
        written_data = stream.read()
        stream.close()
        self.assertEqual(original_data, written_data)

    def test_reading_psl_34_002(self):
        """Test parsing psl_34_002.psl and pslx_34_002.pslx."""
        path = "Blat/psl_34_002.psl"
        self.check_reading_psl_34_002(path)
        path = "Blat/pslx_34_002.pslx"
        self.check_reading_psl_34_002(path)

    def check_reading_psl_34_002(self, path):
        """Check parsing psl_34_002.psl or pslx_34_002.pslx."""
        alignments = Align.parse(path, "psl")
        self.assertEqual(alignments.metadata["psLayout version"], "3")
        self.assertRaises(StopIteration, next, alignments)

    def test_writing_psl_34_002(self):
        """Test writing the alignments in psl_34_002.psl."""
        path = "Blat/psl_34_002.psl"
        with open(path) as stream:
            original_data = stream.read()
        alignments = Align.parse(path, "psl")
        stream = StringIO()
        n = Align.write(alignments, stream, "psl")
        self.assertEqual(n, 0)
        stream.seek(0)
        written_data = stream.read()
        stream.close()
        self.assertEqual(original_data, written_data)

    def test_reading_psl_34_003(self):
        """Test parsing psl_34_003.psl and pslx_34_003.pslx."""
        self.check_reading_psl_34_003("psl")
        self.check_reading_psl_34_003("pslx")

    def check_reading_psl_34_003(self, fmt):
        """Check parsing psl_34_003.psl or pslx_34_003.pslx."""
        path = "Blat/%s_34_003.%s" % (fmt, fmt)
        alignments = Align.parse(path, "psl")
        self.assertEqual(alignments.metadata["psLayout version"], "3")
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 16)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 16))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 191154276)
        self.assertEqual(len(alignment.query.seq), 33)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[61646095:61646111], "aggtaaactgccttca"
            )
            self.assertEqual(
                alignment.query.seq[11:27], self.queries[alignment.query.id][11:27]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['a', 'g', 'g', 't', 'a', 'a', 'a', 'c', 't', 'g', 'c', 'c', 't',
              't', 'c', 'a'],
             ['a', 'g', 'g', 't', 'a', 'a', 'a', 'c', 't', 'g', 'c', 'c', 't',
              't', 'c', 'a']], dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[61646095, 61646111],
                             [      11,       27]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
16	0	0	0	0	0	0	0	+	hg18_dna	33	11	27	chr4	191154276	61646095	61646111	1	16,	11,	61646095,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 33)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 33))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 33)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[10271783:10271816],
                "atgagcttccaaggtaaactgccttcaagattc",
            )
            self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['a', 't', 'g', 'a', 'g', 'c', 't', 't', 'c', 'c', 'a', 'a', 'g',
              'g', 't', 'a', 'a', 'a', 'c', 't', 'g', 'c', 'c', 't', 't', 'c',
              'a', 'a', 'g', 'a', 't', 't', 'c'],
             ['a', 't', 'g', 'a', 'g', 'c', 't', 't', 'c', 'c', 'a', 'a', 'g',
              'g', 't', 'a', 'a', 'a', 'c', 't', 'g', 'c', 'c', 't', 't', 'c',
              'a', 'a', 'g', 'a', 't', 't', 'c']], dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[10271783, 10271816],
                             [       0,       33]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
33	0	0	0	0	0	0	0	+	hg18_dna	33	0	33	chr1	249250621	10271783	10271816	1	33,	0,	10271783,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 17)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 17))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 33)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[53575980:53575997], "aaggcagtttaccttgg"
            )
            self.assertEqual(
                alignment.query.seq[8:25], self.queries[alignment.query.id][8:25]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['a', 'a', 'g', 'g', 'c', 'a', 'g', 't', 't', 't', 'a', 'c', 'c',
              't', 't', 'g', 'g'],
             ['a', 'a', 'g', 'g', 'c', 'a', 'g', 't', 't', 't', 'a', 'c', 'c',
              't', 't', 'g', 'g']], dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[53575980, 53575997],
                             [      25,        8]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
17	0	0	0	0	0	0	0	-	hg18_dna	33	8	25	chr2	243199373	53575980	53575997	1	17,	8,	53575980,
""",
        )
        self.assertRaises(StopIteration, next, alignments)

    def test_writing_psl_34_003(self):
        """Test writing the alignments in psl_34_003.psl."""
        path = "Blat/psl_34_003.psl"
        with open(path) as stream:
            original_data = stream.read()
        alignments = Align.parse(path, "psl")
        stream = StringIO()
        n = Align.write(alignments, stream, "psl")
        self.assertEqual(n, 3)
        stream.seek(0)
        written_data = stream.read()
        stream.close()
        self.assertEqual(original_data, written_data)

    def test_reading_psl_34_004(self):
        """Test parsing psl_34_004.psl and pslx_34_004.pslx."""
        self.check_reading_psl_34_004("psl")
        self.check_reading_psl_34_004("pslx")

    def check_reading_psl_34_004(self, fmt):
        """Check parsing psl_34_004.psl or pslx_34_004.pslx."""
        path = "Blat/%s_34_004.%s" % (fmt, fmt)
        alignments = Align.parse(path, "psl")
        self.assertEqual(alignments.metadata["psLayout version"], "3")
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 38)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr9")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 141213431)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[85737865:85737906],
                "acaaaggggctgggcgcagtggctcacgcctgtaatcccaa",
            )
            self.assertEqual(
                alignment.query.seq[9:50], self.queries[alignment.query.id][9:50]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['a', 'c', 'a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g',
              'g', 'c', 'g', 'c', 'a', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a',
              'c', 'g', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c', 'c',
              'a', 'a'],
             ['a', 'c', 'a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g',
              'g', 'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a',
              'c', 'a', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c', 'c',
              'a', 'a']], dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[85737865, 85737906],
                             [       9,       50]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
38	3	0	0	0	0	0	0	+	hg19_dna	50	9	50	chr9	141213431	85737865	85737906	1	41,	9,	85737865,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 41)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr8")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 146364022)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[95160479:95160520],
                "cacaaaggggctgggcgtggtggctcacacctgtaatccca",
            )
            self.assertEqual(
                alignment.query.seq[8:49], self.queries[alignment.query.id][8:49]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['c', 'a', 'c', 'a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g',
              'g', 'g', 'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c',
              'a', 'c', 'a', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c',
              'c', 'a'],
             ['c', 'a', 'c', 'a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g',
              'g', 'g', 'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c',
              'a', 'c', 'a', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c',
              'c', 'a']], dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[95160479, 95160520],
                             [       8,       49]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
41	0	0	0	0	0	0	0	+	hg19_dna	50	8	49	chr8	146364022	95160479	95160520	1	41,	8,	95160479,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 33)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 51304566)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[42144400:42144436],
                "aaaggggctgggcgtggtagctcatgcctgtaatcc",
            )
            self.assertEqual(
                alignment.query.seq[11:47], self.queries[alignment.query.id][11:47]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g', 'c',
              'g', 't', 'g', 'g', 't', 'a', 'g', 'c', 't', 'c', 'a', 't', 'g',
              'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c'],
             ['a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g', 'c',
              'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a', 'c', 'a',
              'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c']], dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[42144400, 42144436],
                             [      11,       47]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
33	3	0	0	0	0	0	0	+	hg19_dna	50	11	47	chr22	51304566	42144400	42144436	1	36,	11,	42144400,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 43)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 48))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(alignment.target.seq[183925984:183925990], "aaaaat")
            self.assertEqual(
                alignment.target.seq[183925990:183926028],
                "aaaggggctgggcgtggtggctcacgcctgtaatccca",
            )
            self.assertEqual(
                alignment.query.seq[1:7], self.queries[alignment.query.id][1:7]
            )
            self.assertEqual(
                alignment.query.seq[38:49], self.queries[alignment.query.id][38:49]
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[183925984, 183925990, 183925990, 183926028],
                             [        1,         7,        11,        49]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
43	1	0	0	1	4	0	0	+	hg19_dna	50	1	49	chr2	243199373	183925984	183926028	2	6,38,	1,11,	183925984,183925990,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 34)
        self.assertEqual(alignment.misMatches, 2)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 170))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[35483340:35483365], "caaaggggctgggcgtagtggctga"
            )
            self.assertEqual(alignment.target.seq[35483499:35483510], "cacctgtaatc")
            self.assertEqual(
                alignment.query.seq[10:46], self.queries[alignment.query.id][10:46]
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[35483340, 35483365, 35483499, 35483510],
                             [      10,       35,       35,       46]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
34	2	0	0	0	0	1	134	+	hg19_dna	50	10	46	chr19	59128983	35483340	35483510	2	25,11,	10,35,	35483340,35483499,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 39)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 78077248)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[23891310:23891349],
                "caaaggggctgggcgtggtggctcacacctgtaatccca",
            )
            self.assertEqual(
                alignment.query.seq[10:49], self.queries[alignment.query.id][10:49]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['c', 'a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g',
              'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a', 'c',
              'a', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c', 'c', 'a'],
             ['c', 'a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g',
              'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a', 'c',
              'a', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c', 'c', 'a']],
            dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[23891310, 23891349],
                             [      10,       49]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
39	0	0	0	0	0	0	0	+	hg19_dna	50	10	49	chr18	78077248	23891310	23891349	1	39,	10,	23891310,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 27)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 28))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 78077248)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[43252217:43252245], "ggcgtggtggctcacgcctgtaatccca"
            )
            self.assertEqual(
                alignment.query.seq[21:49], self.queries[alignment.query.id][21:49]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['g', 'g', 'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c',
              'a', 'c', 'g', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c',
              'c', 'a'],
             ['g', 'g', 'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c',
              'a', 'c', 'a', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c',
              'c', 'a']], dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[43252217, 43252245],
                             [      21,       49]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
27	1	0	0	0	0	0	0	+	hg19_dna	50	21	49	chr18	78077248	43252217	43252245	1	28,	21,	43252217,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 44)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 54))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 115169878)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(alignment.target.seq[52759147:52759154], "aaaaatt")
            self.assertEqual(
                alignment.target.seq[52759160:52759198],
                "aaaggggctgggcgtggtggctcacgcctgtaatccca",
            )
            self.assertEqual(
                alignment.query.seq[1:8], self.queries[alignment.query.id][1:8]
            )
            self.assertEqual(
                alignment.query.seq[38:49], self.queries[alignment.query.id][38:49]
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[52759147, 52759154, 52759160, 52759160, 52759198],
                             [       1,        8,        8,       11,       49]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
44	1	0	0	1	3	1	6	+	hg19_dna	50	1	49	chr13	115169878	52759147	52759198	2	7,38,	1,11,	52759147,52759160,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 50)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 50))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[1207056:1207106],
                "caaaaattcacaaaggggctgggcgtggtggctcacacctgtaatcccaa",
            )
            self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['c', 'a', 'a', 'a', 'a', 'a', 't', 't', 'c', 'a', 'c', 'a', 'a',
              'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g', 'c', 'g', 't',
              'g', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a', 'c', 'a', 'c', 'c',
              't', 'g', 't', 'a', 'a', 't', 'c', 'c', 'c', 'a', 'a'],
             ['c', 'a', 'a', 'a', 'a', 'a', 't', 't', 'c', 'a', 'c', 'a', 'a',
              'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g', 'c', 'g', 't',
              'g', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a', 'c', 'a', 'c', 'c',
              't', 'g', 't', 'a', 'a', 't', 'c', 'c', 'c', 'a', 'a']],
            dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[1207056, 1207106],
                             [      0,      50]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
50	0	0	0	0	0	0	0	+	hg19_dna	50	0	50	chr1	249250621	1207056	1207106	1	50,	0,	1207056,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 31)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[61700837:61700871],
                "aaaaatgaacaaaggggctgggcgcggtggctca",
            )
            self.assertEqual(
                alignment.query.seq[1:35], self.queries[alignment.query.id][1:35]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['a', 'a', 'a', 'a', 'a', 't', 'g', 'a', 'a', 'c', 'a', 'a', 'a',
              'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g', 'c', 'g', 'c', 'g',
              'g', 't', 'g', 'g', 'c', 't', 'c', 'a'],
             ['a', 'a', 'a', 'a', 'a', 't', 't', 'c', 'a', 'c', 'a', 'a', 'a',
              'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g', 'c', 'g', 't', 'g',
              'g', 't', 'g', 'g', 'c', 't', 'c', 'a']], dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[61700837, 61700871],
                             [       1,       35]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
31	3	0	0	0	0	0	0	+	hg19_dna	50	1	35	chr1	249250621	61700837	61700871	1	34,	1,	61700837,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 28)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 44))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 191154276)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(alignment.target.seq[37558157:37558167], "tgggattaca")
            self.assertEqual(
                alignment.target.seq[37558173:37558191], "accacgcccagccccttt"
            )
            self.assertEqual(
                alignment.query.seq[11:29], self.queries[alignment.query.id][11:29]
            )
            self.assertEqual(
                alignment.query.seq[39:49], self.queries[alignment.query.id][39:49]
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[37558157, 37558167, 37558173, 37558173, 37558191],
                             [      49,       39,       39,       29,       11]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
28	0	0	0	1	10	1	6	-	hg19_dna	50	11	49	chr4	191154276	37558157	37558191	2	10,18,	1,21,	37558157,37558173,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 35)
        self.assertEqual(alignment.misMatches, 2)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 37))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 51304566)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[48997405:48997442],
                "tgggattacaggcgggagccaccacgcccagcccctt",
            )
            self.assertEqual(
                alignment.query.seq[12:49], self.queries[alignment.query.id][12:49]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 'c',
              'g', 'g', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
              'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't'],
             ['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
              'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
              'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't']],
            dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[48997405, 48997442],
                             [      49,       12]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
35	2	0	0	0	0	0	0	-	hg19_dna	50	12	49	chr22	51304566	48997405	48997442	1	37,	1,	48997405,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 35)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[120641740:120641776],
                "tgggattacaggcgtgagccaccacgcccagcccct",
            )
            self.assertEqual(
                alignment.query.seq[13:49], self.queries[alignment.query.id][13:49]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 'c',
              'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
              'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't'],
             ['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
              'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
              'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't']], dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[120641740, 120641776],
                             [       49,        13]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
35	1	0	0	0	0	0	0	-	hg19_dna	50	13	49	chr2	243199373	120641740	120641776	1	36,	1,	120641740,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 39)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[54017130:54017169],
                "tgggattacaggtgtgagccaccacgcccagcccctttg",
            )
            self.assertEqual(
                alignment.query.seq[10:49], self.queries[alignment.query.id][10:49]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
              'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
              'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't', 't', 'g'],
             ['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
              'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
              'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't', 't', 'g']],
            dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[54017130, 54017169],
                             [      49,       10]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
39	0	0	0	0	0	0	0	-	hg19_dna	50	10	49	chr19	59128983	54017130	54017169	1	39,	1,	54017130,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 36)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[553742:553781],
                "tgggatgacaggggtgaggcaccacgcccagcccctttg",
            )
            self.assertEqual(
                alignment.query.seq[10:49], self.queries[alignment.query.id][10:49]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['t', 'g', 'g', 'g', 'a', 't', 'g', 'a', 'c', 'a', 'g', 'g', 'g',
              'g', 't', 'g', 'a', 'g', 'g', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
              'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't', 't', 'g'],
             ['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
              'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
              'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't', 't', 'g']],
            dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[553742, 553781],
                             [    49,     10]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
36	3	0	0	0	0	0	0	-	hg19_dna	50	10	49	chr19	59128983	553742	553781	1	39,	1,	553742,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 33)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 135534747)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[99388555:99388591],
                "tgggattataggcatgagccaccacgcccagcccct",
            )
            self.assertEqual(
                alignment.query.seq[13:49], self.queries[alignment.query.id][13:49]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 't', 'a', 'g', 'g', 'c',
              'a', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
              'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't'],
             ['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
              'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
              'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't']], dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[99388555, 99388591],
                             [      49,       13]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
33	3	0	0	0	0	0	0	-	hg19_dna	50	13	49	chr10	135534747	99388555	99388591	1	36,	1,	99388555,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 24)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 25))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 135534747)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[112178171:112178196], "tgagtcaccacgcccagcccctttg"
            )
            self.assertEqual(
                alignment.query.seq[10:35], self.queries[alignment.query.id][10:35]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['t', 'g', 'a', 'g', 't', 'c', 'a', 'c', 'c', 'a', 'c', 'g', 'c',
              'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't', 't', 'g'],
             ['t', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g', 'c',
              'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't', 't', 'g']],
            dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[112178171, 112178196],
                             [       35,        10]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
24	1	0	0	0	0	0	0	-	hg19_dna	50	10	35	chr10	135534747	112178171	112178196	1	25,	15,	112178171,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 35)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[39368490:39368526],
                "tgggattacaggcgtgagccaccacgcccagcccct",
            )
            self.assertEqual(
                alignment.query.seq[13:49], self.queries[alignment.query.id][13:49]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 'c',
              'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
              'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't'],
             ['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
              'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
              'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't']], dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[39368490, 39368526],
                             [      49,       13]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
35	1	0	0	0	0	0	0	-	hg19_dna	50	13	49	chr1	249250621	39368490	39368526	1	36,	1,	39368490,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 33)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[220325687:220325721],
                "ggattacaggcgtgagccaccacgcccagcccct",
            )
            self.assertEqual(
                alignment.query.seq[13:47], self.queries[alignment.query.id][13:47]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 'c', 'g', 't',
              'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g', 'c', 'c',
              'c', 'a', 'g', 'c', 'c', 'c', 'c', 't'],
             ['g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't', 'g', 't',
              'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g', 'c', 'c',
              'c', 'a', 'g', 'c', 'c', 'c', 'c', 't']], dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[220325687, 220325721],
                             [       47,        13]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
33	1	0	0	0	0	0	0	-	hg19_dna	50	13	47	chr1	249250621	220325687	220325721	1	34,	3,	220325687,
""",
        )
        self.assertRaises(StopIteration, next, alignments)

    def test_writing_psl_34_004(self):
        """Test writing the alignments in psl_34_004.psl."""
        path = "Blat/psl_34_004.psl"
        with open(path) as stream:
            original_data = stream.read()
        alignments = Align.parse(path, "psl")
        stream = StringIO()
        n = Align.write(alignments, stream, "psl")
        self.assertEqual(n, 19)
        stream.seek(0)
        written_data = stream.read()
        stream.close()
        self.assertEqual(original_data, written_data)

    def test_reading_psl_34_005(self):
        """Test parsing psl_34_005.psl and pslx_34_005.pslx."""
        self.check_reading_psl_34_005("psl")
        self.check_reading_psl_34_005("pslx")

    def check_reading_psl_34_005(self, fmt):
        """Check parsing psl_34_005.psl or pslx_34_005.pslx."""
        path = "Blat/%s_34_005.%s" % (fmt, fmt)
        alignments = Align.parse(path, "psl")
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 16)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 16))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 191154276)
        self.assertEqual(len(alignment.query.seq), 33)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[61646095:61646111], "aggtaaactgccttca"
            )
            self.assertEqual(
                alignment.query.seq[11:27], self.queries[alignment.query.id][11:27]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['a', 'g', 'g', 't', 'a', 'a', 'a', 'c', 't', 'g', 'c', 'c', 't',
              't', 'c', 'a'],
             ['a', 'g', 'g', 't', 'a', 'a', 'a', 'c', 't', 'g', 'c', 'c', 't',
              't', 'c', 'a']], dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[61646095, 61646111],
                             [      11,       27]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
16	0	0	0	0	0	0	0	+	hg18_dna	33	11	27	chr4	191154276	61646095	61646111	1	16,	11,	61646095,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 33)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 33))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 33)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[10271783:10271816],
                "atgagcttccaaggtaaactgccttcaagattc",
            )
            self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['a', 't', 'g', 'a', 'g', 'c', 't', 't', 'c', 'c', 'a', 'a', 'g',
              'g', 't', 'a', 'a', 'a', 'c', 't', 'g', 'c', 'c', 't', 't', 'c',
              'a', 'a', 'g', 'a', 't', 't', 'c'],
             ['a', 't', 'g', 'a', 'g', 'c', 't', 't', 'c', 'c', 'a', 'a', 'g',
              'g', 't', 'a', 'a', 'a', 'c', 't', 'g', 'c', 'c', 't', 't', 'c',
              'a', 'a', 'g', 'a', 't', 't', 'c']], dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[10271783, 10271816],
                             [       0,       33]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
33	0	0	0	0	0	0	0	+	hg18_dna	33	0	33	chr1	249250621	10271783	10271816	1	33,	0,	10271783,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 17)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 17))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 33)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[53575980:53575997], "aaggcagtttaccttgg"
            )
            self.assertEqual(
                alignment.query.seq[8:25], self.queries[alignment.query.id][8:25]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['a', 'a', 'g', 'g', 'c', 'a', 'g', 't', 't', 't', 'a', 'c', 'c',
              't', 't', 'g', 'g'],
             ['a', 'a', 'g', 'g', 'c', 'a', 'g', 't', 't', 't', 'a', 'c', 'c',
              't', 't', 'g', 'g']], dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[53575980, 53575997],
                             [      25,        8]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
17	0	0	0	0	0	0	0	-	hg18_dna	33	8	25	chr2	243199373	53575980	53575997	1	17,	8,	53575980,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 38)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr9")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 141213431)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[85737865:85737906],
                "acaaaggggctgggcgcagtggctcacgcctgtaatcccaa",
            )
            self.assertEqual(
                alignment.query.seq[9:], self.queries[alignment.query.id][9:]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['a', 'c', 'a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g',
              'g', 'c', 'g', 'c', 'a', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a',
              'c', 'g', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c', 'c',
              'a', 'a'],
             ['a', 'c', 'a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g',
              'g', 'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a',
              'c', 'a', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c', 'c',
              'a', 'a']], dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[85737865, 85737906],
                             [       9,       50]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
38	3	0	0	0	0	0	0	+	hg19_dna	50	9	50	chr9	141213431	85737865	85737906	1	41,	9,	85737865,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 41)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr8")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 146364022)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[95160479:95160520],
                "cacaaaggggctgggcgtggtggctcacacctgtaatccca",
            )
            self.assertEqual(
                alignment.query.seq[8:49], self.queries[alignment.query.id][8:49]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['c', 'a', 'c', 'a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g',
              'g', 'g', 'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c',
              'a', 'c', 'a', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c',
              'c', 'a'],
             ['c', 'a', 'c', 'a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g',
              'g', 'g', 'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c',
              'a', 'c', 'a', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c',
              'c', 'a']], dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[95160479, 95160520],
                             [       8,       49]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
41	0	0	0	0	0	0	0	+	hg19_dna	50	8	49	chr8	146364022	95160479	95160520	1	41,	8,	95160479,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 33)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 51304566)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[42144400:42144436],
                "aaaggggctgggcgtggtagctcatgcctgtaatcc",
            )
            self.assertEqual(
                alignment.query.seq[11:47], self.queries[alignment.query.id][11:47]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g', 'c',
              'g', 't', 'g', 'g', 't', 'a', 'g', 'c', 't', 'c', 'a', 't', 'g',
              'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c'],
             ['a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g', 'c',
              'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a', 'c', 'a',
              'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c']], dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[42144400, 42144436],
                             [      11,       47]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
33	3	0	0	0	0	0	0	+	hg19_dna	50	11	47	chr22	51304566	42144400	42144436	1	36,	11,	42144400,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 43)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 48))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(alignment.target.seq[183925984:183925990], "aaaaat")
            self.assertEqual(
                alignment.target.seq[183925990:183926028],
                "aaaggggctgggcgtggtggctcacgcctgtaatccca",
            )
            self.assertEqual(
                alignment.query.seq[1:7], self.queries[alignment.query.id][1:7]
            )
            self.assertEqual(
                alignment.query.seq[11:49], self.queries[alignment.query.id][11:49]
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[183925984, 183925990, 183925990, 183926028],
                             [        1,         7,        11,        49]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
43	1	0	0	1	4	0	0	+	hg19_dna	50	1	49	chr2	243199373	183925984	183926028	2	6,38,	1,11,	183925984,183925990,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 34)
        self.assertEqual(alignment.misMatches, 2)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 170))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[35483340:35483365], "caaaggggctgggcgtagtggctga"
            )
            self.assertEqual(alignment.target.seq[35483499:35483510], "cacctgtaatc")
            self.assertEqual(
                alignment.query.seq[10:46], self.queries[alignment.query.id][10:46]
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[35483340, 35483365, 35483499, 35483510],
                             [      10,       35,       35,       46]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
34	2	0	0	0	0	1	134	+	hg19_dna	50	10	46	chr19	59128983	35483340	35483510	2	25,11,	10,35,	35483340,35483499,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 39)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 78077248)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[23891310:23891349],
                "caaaggggctgggcgtggtggctcacacctgtaatccca",
            )
            self.assertEqual(
                alignment.query.seq[10:49], self.queries[alignment.query.id][10:49]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['c', 'a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g',
              'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a', 'c',
              'a', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c', 'c', 'a'],
             ['c', 'a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g',
              'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a', 'c',
              'a', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c', 'c', 'a']],
            dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[23891310, 23891349],
                             [      10,       49]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
39	0	0	0	0	0	0	0	+	hg19_dna	50	10	49	chr18	78077248	23891310	23891349	1	39,	10,	23891310,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 27)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 28))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 78077248)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[43252217:43252245], "ggcgtggtggctcacgcctgtaatccca"
            )
            self.assertEqual(
                alignment.query.seq[21:49], self.queries[alignment.query.id][21:49]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['g', 'g', 'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c',
              'a', 'c', 'g', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c',
              'c', 'a'],
             ['g', 'g', 'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c',
              'a', 'c', 'a', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c',
              'c', 'a']], dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[43252217, 43252245],
                             [      21,       49]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
27	1	0	0	0	0	0	0	+	hg19_dna	50	21	49	chr18	78077248	43252217	43252245	1	28,	21,	43252217,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 44)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 54))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 115169878)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(alignment.target.seq[52759147:52759154], "aaaaatt")
            self.assertEqual(
                alignment.target.seq[52759160:52759198],
                "aaaggggctgggcgtggtggctcacgcctgtaatccca",
            )
            self.assertEqual(
                alignment.query.seq[1:8], self.queries[alignment.query.id][1:8]
            )
            self.assertEqual(
                alignment.query.seq[38:49], self.queries[alignment.query.id][38:49]
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[52759147, 52759154, 52759160, 52759160, 52759198],
                             [       1,        8,        8,       11,       49]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
44	1	0	0	1	3	1	6	+	hg19_dna	50	1	49	chr13	115169878	52759147	52759198	2	7,38,	1,11,	52759147,52759160,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 50)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 50))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[1207056:1207106],
                "caaaaattcacaaaggggctgggcgtggtggctcacacctgtaatcccaa",
            )
            self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['c', 'a', 'a', 'a', 'a', 'a', 't', 't', 'c', 'a', 'c', 'a', 'a',
              'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g', 'c', 'g', 't',
              'g', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a', 'c', 'a', 'c', 'c',
              't', 'g', 't', 'a', 'a', 't', 'c', 'c', 'c', 'a', 'a'],
             ['c', 'a', 'a', 'a', 'a', 'a', 't', 't', 'c', 'a', 'c', 'a', 'a',
              'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g', 'c', 'g', 't',
              'g', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a', 'c', 'a', 'c', 'c',
              't', 'g', 't', 'a', 'a', 't', 'c', 'c', 'c', 'a', 'a']],
            dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[1207056, 1207106],
                             [      0,      50]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
50	0	0	0	0	0	0	0	+	hg19_dna	50	0	50	chr1	249250621	1207056	1207106	1	50,	0,	1207056,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 31)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[61700837:61700871],
                "aaaaatgaacaaaggggctgggcgcggtggctca",
            )
            self.assertEqual(
                alignment.query.seq[1:35], self.queries[alignment.query.id][1:35]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['a', 'a', 'a', 'a', 'a', 't', 'g', 'a', 'a', 'c', 'a', 'a', 'a',
              'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g', 'c', 'g', 'c', 'g',
              'g', 't', 'g', 'g', 'c', 't', 'c', 'a'],
             ['a', 'a', 'a', 'a', 'a', 't', 't', 'c', 'a', 'c', 'a', 'a', 'a',
              'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g', 'c', 'g', 't', 'g',
              'g', 't', 'g', 'g', 'c', 't', 'c', 'a']], dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[61700837, 61700871],
                             [       1,       35]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
31	3	0	0	0	0	0	0	+	hg19_dna	50	1	35	chr1	249250621	61700837	61700871	1	34,	1,	61700837,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 28)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 44))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 191154276)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(alignment.target.seq[37558157:37558167], "tgggattaca")
            self.assertEqual(
                alignment.target.seq[37558173:37558191], "accacgcccagccccttt"
            )
            self.assertEqual(
                alignment.query.seq[11:29], self.queries[alignment.query.id][11:29]
            )
            self.assertEqual(
                alignment.query.seq[39:49], self.queries[alignment.query.id][39:49]
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[37558157, 37558167, 37558173, 37558173, 37558191],
                             [      49,       39,       39,       29,       11]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
28	0	0	0	1	10	1	6	-	hg19_dna	50	11	49	chr4	191154276	37558157	37558191	2	10,18,	1,21,	37558157,37558173,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 35)
        self.assertEqual(alignment.misMatches, 2)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 37))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 51304566)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[48997405:48997442],
                "tgggattacaggcgggagccaccacgcccagcccctt",
            )
            self.assertEqual(
                alignment.query.seq[12:49], self.queries[alignment.query.id][12:49]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 'c',
              'g', 'g', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
              'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't'],
             ['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
              'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
              'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't']],
            dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[48997405, 48997442],
                             [      49,       12]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
35	2	0	0	0	0	0	0	-	hg19_dna	50	12	49	chr22	51304566	48997405	48997442	1	37,	1,	48997405,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 35)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[120641740:120641776],
                "tgggattacaggcgtgagccaccacgcccagcccct",
            )
            self.assertEqual(
                alignment.query.seq[13:49], self.queries[alignment.query.id][13:49]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 'c',
              'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
              'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't'],
             ['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
              'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
              'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't']], dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[120641740, 120641776],
                             [       49,        13]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
35	1	0	0	0	0	0	0	-	hg19_dna	50	13	49	chr2	243199373	120641740	120641776	1	36,	1,	120641740,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 39)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[54017130:54017169],
                "tgggattacaggtgtgagccaccacgcccagcccctttg",
            )
            self.assertEqual(
                alignment.query.seq[10:49], self.queries[alignment.query.id][10:49]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
              'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
              'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't', 't', 'g'],
             ['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
              'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
              'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't', 't', 'g']],
            dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[54017130, 54017169],
                             [      49,       10]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
39	0	0	0	0	0	0	0	-	hg19_dna	50	10	49	chr19	59128983	54017130	54017169	1	39,	1,	54017130,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 36)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[553742:553781],
                "tgggatgacaggggtgaggcaccacgcccagcccctttg",
            )
            self.assertEqual(
                alignment.query.seq[10:49], self.queries[alignment.query.id][10:49]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['t', 'g', 'g', 'g', 'a', 't', 'g', 'a', 'c', 'a', 'g', 'g', 'g',
              'g', 't', 'g', 'a', 'g', 'g', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
              'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't', 't', 'g'],
             ['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
              'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
              'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't', 't', 'g']],
            dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[553742, 553781],
                             [    49,     10]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
36	3	0	0	0	0	0	0	-	hg19_dna	50	10	49	chr19	59128983	553742	553781	1	39,	1,	553742,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 33)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 135534747)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[99388555:99388591],
                "tgggattataggcatgagccaccacgcccagcccct",
            )
            self.assertEqual(
                alignment.query.seq[13:49], self.queries[alignment.query.id][13:49]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 't', 'a', 'g', 'g', 'c',
              'a', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
              'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't'],
             ['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
              'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
              'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't']], dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[99388555, 99388591],
                             [      49,       13]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
33	3	0	0	0	0	0	0	-	hg19_dna	50	13	49	chr10	135534747	99388555	99388591	1	36,	1,	99388555,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 24)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 25))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 135534747)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[112178171:112178196], "tgagtcaccacgcccagcccctttg"
            )
            self.assertEqual(
                alignment.query.seq[10:35], self.queries[alignment.query.id][10:35]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['t', 'g', 'a', 'g', 't', 'c', 'a', 'c', 'c', 'a', 'c', 'g', 'c',
              'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't', 't', 'g'],
             ['t', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g', 'c',
              'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't', 't', 'g']],
            dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[112178171, 112178196],
                             [       35,        10]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
24	1	0	0	0	0	0	0	-	hg19_dna	50	10	35	chr10	135534747	112178171	112178196	1	25,	15,	112178171,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 35)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[39368490:39368526],
                "tgggattacaggcgtgagccaccacgcccagcccct",
            )
            self.assertEqual(
                alignment.query.seq[13:49], self.queries[alignment.query.id][13:49]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 'c',
              'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
              'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't'],
             ['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
              'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
              'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't']], dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
                numpy.array([[39368490, 39368526],
                             [      49,       13]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
35	1	0	0	0	0	0	0	-	hg19_dna	50	13	49	chr1	249250621	39368490	39368526	1	36,	1,	39368490,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 33)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[220325687:220325721],
                "ggattacaggcgtgagccaccacgcccagcccct",
            )
            self.assertEqual(
                alignment.query.seq[13:47], self.queries[alignment.query.id][13:47]
            )
            self.assertTrue(
                numpy.array_equal(
                    numpy.array(alignment, "U"),
                    # fmt: off
# flake8: noqa
numpy.array([['g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 'c', 'g', 't',
              'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g', 'c', 'c',
              'c', 'a', 'g', 'c', 'c', 'c', 'c', 't'],
             ['g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't', 'g', 't',
              'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g', 'c', 'c',
              'c', 'a', 'g', 'c', 'c', 'c', 'c', 't']], dtype='U')
                    # fmt: on
                )
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[220325687, 220325721],
                             [       47,        13]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
33	1	0	0	0	0	0	0	-	hg19_dna	50	13	47	chr1	249250621	220325687	220325721	1	34,	3,	220325687,
""",
        )
        self.assertRaises(StopIteration, next, alignments)

    def test_writing_psl_34_005(self):
        """Test writing the alignments in psl_34_005.psl."""
        path = "Blat/psl_34_005.psl"
        with open(path) as stream:
            original_data = stream.read()
        alignments = Align.parse(path, "psl")
        stream = StringIO()
        n = Align.write(alignments, stream, "psl", header=False)
        self.assertEqual(n, 22)
        stream.seek(0)
        written_data = stream.read()
        stream.close()
        self.assertEqual(original_data, written_data)


class TestAlign_dnax_prot(unittest.TestCase):
    def test_reading_psl_35_001(self):
        """Test parsing psl_35_001.psl and pslx_35_001.pslx."""
        self.check_reading_psl_35_001("psl")
        self.check_reading_psl_35_001("pslx")

    def check_reading_psl_35_001(self, fmt):
        """Check parsing psl_35_001.psl or pslx_35_001.pslx."""
        path = "Blat/%s_35_001.%s" % (fmt, fmt)
        alignments = Align.parse(path, "psl")
        self.assertEqual(alignments.metadata["psLayout version"], "3")
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 52)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertEqual(len(alignment.target.seq), 114364328)
        self.assertEqual(len(alignment.query.seq), 230)
        if fmt == "pslx":
            self.assertEqual(
                alignment.query.seq[61:113],
                "YEVFRTEEEEKIKSQGQDVTSSVYFMKQTISNACGTIGLIHAIANNKDKMHF",
            )
            self.assertEqual(len(alignment.target.seq.defined_ranges), 0)
            self.assertEqual(len(alignment.target.features), 1)
            feature = alignment.target.features[0]
            self.assertEqual(
                feature.location,
                SimpleLocation(
                    ExactPosition(75566694), ExactPosition(75566850), strand=+1
                ),
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(len(feature.qualifiers["translation"]), 1)
            self.assertEqual(
                feature.qualifiers["translation"][0],
                "YEVFRTEEEEKIKSQGQDVTSSVYFMKQTISNACGTIGLIHAIANNKDKMHF",
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[75566694, 75566850],
                             [      61,      113]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
52	0	0	0	0	0	0	0	++	CAG33136.1	230	61	113	chr13	114364328	75566694	75566850	1	52,	61,	75566694,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 44)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertEqual(len(alignment.target.seq), 114364328)
        self.assertEqual(len(alignment.query.seq), 230)
        if fmt == "pslx":
            self.assertEqual(
                alignment.query.seq[17:61],
                "QFLKQLGLHPNWQFVDVYGMDPELLSMVPRPVCAVLLLFPITEK",
            )
            self.assertEqual(len(alignment.target.seq.defined_ranges), 0)
            self.assertEqual(len(alignment.target.features), 1)
            feature = alignment.target.features[0]
            self.assertEqual(
                feature.location,
                SimpleLocation(
                    ExactPosition(75560749), ExactPosition(75560881), strand=+1
                ),
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(len(feature.qualifiers["translation"]), 1)
            self.assertEqual(
                feature.qualifiers["translation"][0],
                "QFLKQLGLHPNWQFVDVYGMDPELLSMVPRPVCAVLLLFPITEK",
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[75560749, 75560881],
                             [      17,       61]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
44	0	0	0	0	0	0	0	++	CAG33136.1	230	17	61	chr13	114364328	75560749	75560881	1	44,	17,	75560749,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 44)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertEqual(len(alignment.target.seq), 114364328)
        self.assertEqual(len(alignment.query.seq), 230)
        if fmt == "pslx":
            self.assertEqual(alignment.query.seq[:15], "MEGQRWLPLEANPEV")
            self.assertEqual(
                alignment.query.seq[113:142], "ESGSTLKKFLEESVSMSPEERARYLENYD"
            )
            self.assertEqual(len(alignment.target.seq.defined_ranges), 0)
            self.assertEqual(len(alignment.target.features), 1)
            feature = alignment.target.features[0]
            self.assertEqual(
                feature.location,
                CompoundLocation(
                    [
                        SimpleLocation(
                            ExactPosition(75549820), ExactPosition(75549865), strand=+1
                        ),
                        SimpleLocation(
                            ExactPosition(75567225), ExactPosition(75567312), strand=+1
                        ),
                    ],
                    operator="join",
                ),
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(len(feature.qualifiers["translation"]), 1)
            self.assertEqual(
                feature.qualifiers["translation"][0],
                "MEGQRWLPLEANPEVESGSTLKKFLEESVSMSPEERARYLENYD",
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[75549820, 75549865, 75567225, 75567225, 75567312],
                             [       0,       15,       15,      113,      142]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
44	0	0	0	1	98	1	17360	++	CAG33136.1	230	0	142	chr13	114364328	75549820	75567312	2	15,29,	0,113,	75549820,75567225,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 47)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertEqual(len(alignment.target.seq), 114364328)
        self.assertEqual(len(alignment.query.seq), 230)
        if fmt == "pslx":
            self.assertEqual(
                alignment.query.seq[183:],
                "DGRKPFPINHGETSDETLLEDAIEVCKKFMERDPDELRFNAIALSAA",
            )
            self.assertEqual(len(alignment.target.seq.defined_ranges), 0)
            self.assertEqual(len(alignment.target.features), 1)
            feature = alignment.target.features[0]
            self.assertEqual(
                feature.location,
                CompoundLocation(
                    [
                        SimpleLocation(
                            ExactPosition(75604767), ExactPosition(75604827), strand=+1
                        ),
                        SimpleLocation(
                            ExactPosition(75605728), ExactPosition(75605809), strand=+1
                        ),
                    ],
                    operator="join",
                ),
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(len(feature.qualifiers["translation"]), 1)
            self.assertEqual(
                feature.qualifiers["translation"][0],
                "DGRKPFPINHGETSDETLLEDAIEVCKKFMERDPDELRFNAIALSAA",
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[75604767, 75604827, 75605728, 75605809],
                             [     183,      203,      203,      230]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
47	0	0	0	0	0	1	901	++	CAG33136.1	230	183	230	chr13	114364328	75604767	75605809	2	20,27,	183,203,	75604767,75605728,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 25)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertEqual(len(alignment.target.seq), 114364328)
        self.assertEqual(len(alignment.query.seq), 230)
        if fmt == "pslx":
            self.assertEqual(alignment.query.seq[158:183], "APSIDEKVDLHFIALVHVDGHLYEL")
            self.assertEqual(len(alignment.target.seq.defined_ranges), 0)
            self.assertEqual(len(alignment.target.features), 1)
            feature = alignment.target.features[0]
            self.assertEqual(
                feature.location,
                SimpleLocation(
                    ExactPosition(75594914), ExactPosition(75594989), strand=+1
                ),
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(len(feature.qualifiers["translation"]), 1)
            self.assertEqual(
                feature.qualifiers["translation"][0], "APSIDEKVDLHFIALVHVDGHLYEL"
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[75594914, 75594989],
                             [     158,      183]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
25	0	0	0	0	0	0	0	++	CAG33136.1	230	158	183	chr13	114364328	75594914	75594989	1	25,	158,	75594914,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 16)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertEqual(len(alignment.target.seq), 114364328)
        self.assertEqual(len(alignment.query.seq), 230)
        if fmt == "pslx":
            self.assertEqual(alignment.query.seq[142:158], "AIRVTHETSAHEGQTE")
            self.assertEqual(len(alignment.target.seq.defined_ranges), 0)
            self.assertEqual(len(alignment.target.features), 1)
            feature = alignment.target.features[0]
            self.assertEqual(
                feature.location,
                SimpleLocation(
                    ExactPosition(75569459), ExactPosition(75569507), strand=+1
                ),
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(len(feature.qualifiers["translation"]), 1)
            self.assertEqual(feature.qualifiers["translation"][0], "AIRVTHETSAHEGQTE")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[75569459, 75569507],
                             [     142,      158]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
16	0	0	0	0	0	0	0	++	CAG33136.1	230	142	158	chr13	114364328	75569459	75569507	1	16,	142,	75569459,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 26)
        self.assertEqual(alignment.misMatches, 8)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertEqual(len(alignment.target.seq), 190214555)
        self.assertEqual(len(alignment.query.seq), 230)
        if fmt == "pslx":
            self.assertEqual(
                alignment.query.seq[76:110], "GQDVTSSVYFMKQTISNACGTIGLIHAIANNKDK"
            )
            self.assertEqual(len(alignment.target.seq.defined_ranges), 0)
            self.assertEqual(len(alignment.target.features), 1)
            feature = alignment.target.features[0]
            self.assertEqual(
                feature.location,
                SimpleLocation(
                    ExactPosition(41260685), ExactPosition(41260787), strand=+1
                ),
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(len(feature.qualifiers["translation"]), 1)
            self.assertEqual(
                feature.qualifiers["translation"][0],
                "GQEVSPKVYFMKQTIGNSCGTIGLIHAVANNQDK",
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[41260685, 41260787],
                             [      76,      110]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
26	8	0	0	0	0	0	0	++	CAG33136.1	230	76	110	chr4	190214555	41260685	41260787	1	34,	76,	41260685,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 37)
        self.assertEqual(alignment.misMatches, 26)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertEqual(len(alignment.target.seq), 190214555)
        self.assertEqual(len(alignment.query.seq), 230)
        if fmt == "pslx":
            self.assertEqual(
                alignment.query.seq[17:59], "QFLKQLGLHPNWQFVDVYGMDPELLSMVPRPVCAVLLLFPIT"
            )
            self.assertEqual(alignment.query.seq[162:183], "DEKVDLHFIALVHVDGHLYEL")
            self.assertEqual(len(alignment.target.seq.defined_ranges), 0)
            self.assertEqual(len(alignment.target.features), 1)
            feature = alignment.target.features[0]
            self.assertEqual(
                feature.location,
                CompoundLocation(
                    [
                        SimpleLocation(
                            ExactPosition(41257605), ExactPosition(41257731), strand=+1
                        ),
                        SimpleLocation(
                            ExactPosition(41263227), ExactPosition(41263290), strand=+1
                        ),
                    ],
                    operator="join",
                ),
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(len(feature.qualifiers["translation"]), 1)
            self.assertEqual(
                feature.qualifiers["translation"][0],
                "QVLSRLGVAGQWRFVDVLGLEEESLGSVPAPACALLLLFPLTDDKVNFHFILFNNVDGHLYEL",
            )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[41257605, 41257731, 41263227, 41263227, 41263290],
                             [      17,       59,       59,      162,      183]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
37	26	0	0	1	103	1	5496	++	CAG33136.1	230	17	183	chr4	190214555	41257605	41263290	2	42,21,	17,162,	41257605,41263227,
""",
        )
        self.assertRaises(StopIteration, next, alignments)

    def test_writing_psl_35_001(self):
        """Test writing the alignments in psl_35_001.psl."""
        path = "Blat/psl_35_001.psl"
        with open(path) as stream:
            original_data = stream.read()
        alignments = Align.parse(path, "psl")
        stream = StringIO()
        n = Align.write(alignments, stream, "psl")
        self.assertEqual(n, 8)
        stream.seek(0)
        written_data = stream.read()
        stream.close()
        self.assertEqual(original_data, written_data)
        # Convert the alignment to a protein alignment and insert the
        # appropriate sequence data. Write this alignment in a PSL file;
        # the writer will recalculate the values for matches, misMatches,
        # repMatches, and nCount from the sequence data and the alignment.
        #
        # The alignments were generated using
        # blat -t=dnax -q=prot hg38.2bit CAG33136.1.fasta psl_35_001.psl
        #
        # To save disk space, we extracted the necessary sequence data using
        #
        # twoBitToFa hg38.2bit:chr13:75549820-75605809 stdout
        # twoBitToFa hg38.2bit:chr4:41257605-41263290 stdout
        #
        # and concatenating the results into file hg38.fa. We will use this
        # file below, and create partially defined Seq objects.
        #
        # Load the protein sequence:
        protein = SeqIO.read("Blat/CAG33136.1.fasta", "fasta")
        protein_alignments = []
        alignments = Align.parse(path, "psl")
        for i, alignment in enumerate(alignments):
            records = SeqIO.parse("Blat/hg38.fa", "fasta")
            for record in records:
                name, start_end = record.id.split(":")
                if name == alignment.sequences[0].id:
                    break
            else:
                raise Exception("Failed to find DNA sequence")
            start, end = start_end.split("-")
            start = int(start)
            end = int(end)
            length = len(alignment.sequences[0])
            sequence = str(record.seq)
            dna = Seq({start: sequence}, length=length)
            alignment.sequences[0].seq = dna
            self.assertEqual(alignment.sequences[1].id, protein.id)
            alignment.sequences[1].seq = protein.seq
            # The alignment is on the forward strand of the DNA sequence:
            self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
            # The protein alignment is also in the forward orientation:
            self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
            # Now extract the aligned sequences:
            aligned_dna = ""
            aligned_protein = ""
            for start, end in alignment.aligned[0]:
                aligned_dna += alignment.sequences[0].seq[start:end]
            for start, end in alignment.aligned[1]:
                aligned_protein += alignment.sequences[1].seq[start:end]
            # Translate the aligned DNA sequence:
            aligned_dna = Seq(aligned_dna)
            aligned_dna_translated = Seq(aligned_dna.translate())
            aligned_protein = Seq(aligned_protein)
            # Create a new alignment including the aligned sequences only:
            records = [
                SeqRecord(aligned_dna_translated, id=alignment.sequences[0].id),
                SeqRecord(aligned_protein, id=alignment.sequences[1].id),
            ]
            coordinates = numpy.array(
                [[0, len(aligned_dna_translated)], [0, len(aligned_protein)]]
            )
            protein_alignment = Alignment(records, coordinates)
            protein_alignments.append(protein_alignment)
            if i == 0:
                self.assertEqual(
                    str(protein_alignment),
                    """\
chr13             0 YEVFRTEEEEKIKSQGQDVTSSVYFMKQTISNACGTIGLIHAIANNKDKMHF 52
                  0 |||||||||||||||||||||||||||||||||||||||||||||||||||| 52
CAG33136.         0 YEVFRTEEEEKIKSQGQDVTSSVYFMKQTISNACGTIGLIHAIANNKDKMHF 52
""",
                )
            elif i == 1:
                self.assertEqual(
                    str(protein_alignment),
                    """\
chr13             0 QFLKQLGLHPNWQFVDVYGMDPELLSMVPRPVCAVLLLFPITEK 44
                  0 |||||||||||||||||||||||||||||||||||||||||||| 44
CAG33136.         0 QFLKQLGLHPNWQFVDVYGMDPELLSMVPRPVCAVLLLFPITEK 44
""",
                )
            elif i == 2:
                self.assertEqual(
                    str(protein_alignment),
                    """\
chr13             0 MEGQRWLPLEANPEVESGSTLKKFLEESVSMSPEERARYLENYD 44
                  0 |||||||||||||||||||||||||||||||||||||||||||| 44
CAG33136.         0 MEGQRWLPLEANPEVESGSTLKKFLEESVSMSPEERARYLENYD 44
""",
                )
            elif i == 3:
                self.assertEqual(
                    str(protein_alignment),
                    """\
chr13             0 DGRKPFPINHGETSDETLLEDAIEVCKKFMERDPDELRFNAIALSAA 47
                  0 ||||||||||||||||||||||||||||||||||||||||||||||| 47
CAG33136.         0 DGRKPFPINHGETSDETLLEDAIEVCKKFMERDPDELRFNAIALSAA 47
""",
                )
            elif i == 4:
                self.assertEqual(
                    str(protein_alignment),
                    """\
chr13             0 APSIDEKVDLHFIALVHVDGHLYEL 25
                  0 ||||||||||||||||||||||||| 25
CAG33136.         0 APSIDEKVDLHFIALVHVDGHLYEL 25
""",
                )
            elif i == 5:
                self.assertEqual(
                    str(protein_alignment),
                    """\
chr13             0 AIRVTHETSAHEGQTE 16
                  0 |||||||||||||||| 16
CAG33136.         0 AIRVTHETSAHEGQTE 16
""",
                )
            elif i == 6:
                self.assertEqual(
                    str(protein_alignment),
                    """\
chr4              0 GQEVSPKVYFMKQTIGNSCGTIGLIHAVANNQDK 34
                  0 ||.|...||||||||.|.|||||||||.|||.|| 34
CAG33136.         0 GQDVTSSVYFMKQTISNACGTIGLIHAIANNKDK 34
""",
                )
            elif i == 7:
                self.assertEqual(
                    str(protein_alignment),
                    """\
chr4              0 QVLSRLGVAGQWRFVDVLGLEEESLGSVPAPACALLLLFPLTDDKVNFHFILFNNVDGHL
                  0 |.|..||....|.||||.|...|.|..||.|.||.|||||.||.||..|||....|||||
CAG33136.         0 QFLKQLGLHPNWQFVDVYGMDPELLSMVPRPVCAVLLLFPITDEKVDLHFIALVHVDGHL

chr4             60 YEL 63
                 60 ||| 63
CAG33136.        60 YEL 63
""",
                )
        # Write the protein alignments to a PSL file:
        stream = StringIO()
        n = Align.write(protein_alignments, stream, "psl", wildcard="X")
        self.assertEqual(n, 8)
        # Read the alignments back in:
        alignments = Align.parse(path, "psl")
        stream.seek(0)
        protein_alignments = Align.parse(stream, "psl")
        for alignment, protein_alignment in zip(alignments, protein_alignments):
            # Confirm that the recalculated values for matches, misMatches,
            # repMatches, and nCount are correct:
            self.assertEqual(alignment.matches, protein_alignment.matches)
            self.assertEqual(alignment.misMatches, protein_alignment.misMatches)
            self.assertEqual(alignment.repMatches, protein_alignment.repMatches)
            self.assertEqual(alignment.nCount, protein_alignment.nCount)

    def test_reading_psl_35_002(self):
        """Test parsing psl_35_002.psl."""
        # See below for a description of the file balAcu1.fa.
        # We use this file here so we can check the SeqFeatures.
        records = SeqIO.parse("Blat/balAcu1.fa", "fasta")
        self.dna = {}
        for record in records:
            name, start_end = record.id.split(":")
            start, end = start_end.split("-")
            start = int(start)
            end = int(end)
            sequence = str(record.seq)
            self.dna[name] = Seq({start: sequence}, length=end)
        self.check_reading_psl_35_002("psl")
        self.check_reading_psl_35_002("pslx")

    def check_reading_psl_35_002(self, fmt):
        """Check parsing psl_35_002.psl or pslx_35_002.pslx."""
        path = "Blat/%s_35_002.%s" % (fmt, fmt)
        alignments = Align.parse(path, "psl")
        self.assertEqual(alignments.metadata["psLayout version"], "3")
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 210)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "KI537979")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertEqual(len(alignment.target.seq), 14052872)
        self.assertEqual(len(alignment.query.seq), 230)
        if fmt == "pslx":
            self.assertEqual(
                alignment.query.seq[17:],
                "QFLKQLGLHPNWQFVDVYGMDPELLSMVPRPVCAVLLLFPITEKYEVFRTEEEEKIKSQGQDVTSSVYFMKQTISNACGTIGLIHAIANNKDKMHFESGSTLKKFLEESVSMSPEERARYLENYDAIRVTHETSAHEGQTEAPSIDEKVDLHFIALVHVDGHLYELDGRKPFPINHGETSDETLLEDAIEVCKKFMERDPDELRFNAIALSAA",
            )
            self.assertEqual(len(alignment.target.seq.defined_ranges), 0)
            self.assertEqual(len(alignment.target.features), 1)
            feature = alignment.target.features[0]
            self.assertEqual(
                feature.location,
                CompoundLocation(
                    [
                        SimpleLocation(
                            ExactPosition(9712654), ExactPosition(9712786), strand=+1
                        ),
                        SimpleLocation(
                            ExactPosition(9715941), ExactPosition(9716097), strand=+1
                        ),
                        SimpleLocation(
                            ExactPosition(9716445), ExactPosition(9716532), strand=+1
                        ),
                        SimpleLocation(
                            ExactPosition(9718374), ExactPosition(9718422), strand=+1
                        ),
                        SimpleLocation(
                            ExactPosition(9739264), ExactPosition(9739339), strand=+1
                        ),
                        SimpleLocation(
                            ExactPosition(9743706), ExactPosition(9743766), strand=+1
                        ),
                        SimpleLocation(
                            ExactPosition(9744511), ExactPosition(9744592), strand=+1
                        ),
                    ],
                    operator="join",
                ),
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(len(feature.qualifiers["translation"]), 1)
            self.assertEqual(
                feature.qualifiers["translation"][0],
                "QFLKQLGLHPNWQFVDVYGMDPELLSMVPRPVCAVLLLFPITEKYEIFRTEEEEKIKSQGQDVTSSVYFMKQTISNACGTIGLIHAIANNKDKMHFESGSTLKKFLEESASMSPEERARYLENYDAIRVTHETSAHEGQTEAPNIDEKVDLHFIALVHVDGHLYELDGRKPFPINHGETSDETLLEDAIEVCKKFMERDPDELRFNAIALSAA",
            )
            # confirm that the feature coordinates are correct by extracting
            # the feature sequence from the target sequence and tranlating it.
            cds = feature.extract(self.dna[alignment.target.id]).translate()
            self.assertEqual(feature.qualifiers["translation"][0], cds)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[9712654, 9712786, 9715941, 9716097, 9716445, 9716532, 9718374,
                              9718422, 9739264, 9739339, 9743706, 9743766, 9744511, 9744592],
                             [     17,      61,      61,     113,     113,     142,     142,
                                  158,     158,     183,     183,     203,     203,     230]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
210	3	0	0	0	0	6	31299	++	CAG33136.1	230	17	230	KI537979	14052872	9712654	9744592	7	44,52,29,16,25,20,27,	17,61,113,142,158,183,203,	9712654,9715941,9716445,9718374,9739264,9743706,9744511,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 207)
        self.assertEqual(alignment.misMatches, 22)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "KI538594")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertEqual(len(alignment.target.seq), 7819582)
        self.assertEqual(len(alignment.query.seq), 230)
        if fmt == "pslx":
            self.assertEqual(alignment.query.seq[:20], "MEGQRWLPLEANPEVTNQFL")
            self.assertEqual(
                alignment.query.seq[21:],
                "QLGLHPNWQFVDVYGMDPELLSMVPRPVCAVLLLFPITEKYEVFRTEEEEKIKSQGQDVTSSVYFMKQTISNACGTIGLIHAIANNKDKMHFESGSTLKKFLEESVSMSPEERARYLENYDAIRVTHETSAHEGQTEAPSIDEKVDLHFIALVHVDGHLYELDGRKPFPINHGETSDETLLEDAIEVCKKFMERDPDELRFNAIALSAA",
            )
            self.assertEqual(len(alignment.target.seq.defined_ranges), 0)
            self.assertEqual(len(alignment.target.features), 1)
            feature = alignment.target.features[0]
            self.assertEqual(
                feature.location,
                CompoundLocation(
                    [
                        SimpleLocation(
                            ExactPosition(2103463), ExactPosition(2103523), strand=+1
                        ),
                        SimpleLocation(
                            ExactPosition(2103522), ExactPosition(2104149), strand=+1
                        ),
                    ],
                    operator="join",
                ),
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(len(feature.qualifiers["translation"]), 1)
            self.assertEqual(
                feature.qualifiers["translation"][0],
                "MEGQCWLPLEANPEVTNQLLQLGLHPNWQFVDVYGMDPELLSMVPRPVCAVLLLFPITEKYEVFRTEEEEKIKSQGQNITSSGYFMRQTISSACGTIGLIHAIANNKDKMHFESGSTLKKFLEESASLSPEERAIYLENYDSIRVTHKTSDHEGQTEAQNIDEKVDLHFIALVHVDGHLYELDGWKPFPINHGETSDATLLRDAIEVFKKFRERDPDERRFNVIALSAA",
            )
            # confirm that the feature coordinates are correct by extracting
            # the feature sequence from the target sequence and tranlating it.
            cds = feature.extract(self.dna[alignment.target.id]).translate()
            self.assertEqual(feature.qualifiers["translation"][0], cds)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[2103463, 2103523, 2103522, 2103522, 2104149],
                             [      0,      20,      20,      21,     230]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
207	22	0	0	1	1	1	-1	++	CAG33136.1	230	0	230	KI538594	7819582	2103463	2104149	2	20,209,	0,21,	2103463,2103522,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 204)
        self.assertEqual(alignment.misMatches, 6)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "KI537194")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertEqual(len(alignment.target.seq), 37111980)
        self.assertEqual(len(alignment.query.seq), 230)
        if fmt == "pslx":
            self.assertEqual(
                alignment.query.seq[:183],
                "MEGQRWLPLEANPEVTNQFLKQLGLHPNWQFVDVYGMDPELLSMVPRPVCAVLLLFPITEKYEVFRTEEEEKIKSQGQDVTSSVYFMKQTISNACGTIGLIHAIANNKDKMHFESGSTLKKFLEESVSMSPEERARYLENYDAIRVTHETSAHEGQTEAPSIDEKVDLHFIALVHVDGHLYEL",
            )
            self.assertEqual(alignment.query.seq[203:], "DAIEVCKKFMERDPDELRFNAIALSAA")
            self.assertEqual(len(alignment.target.seq.defined_ranges), 0)
            self.assertEqual(len(alignment.target.features), 1)
            feature = alignment.target.features[0]
            self.assertEqual(
                feature.location,
                CompoundLocation(
                    [
                        SimpleLocation(
                            ExactPosition(20872472), ExactPosition(20873021), strand=-1
                        ),
                        SimpleLocation(
                            ExactPosition(20872390), ExactPosition(20872471), strand=-1
                        ),
                    ],
                    operator="join",
                ),
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(len(feature.qualifiers["translation"]), 1)
            self.assertEqual(
                feature.qualifiers["translation"][0],
                "MESQRWLPLEANPEVTNQFLKQLGLHPNWQCVDVYGMDPELLSMVPRPVCAVLLLFPITEKYEIFRTEEEEKTKSQGQDVTSSVYFMKQTISNACGTIGLIHAIANNKDKMHFESGSTLKKFLEESASMSPEERARYLENYDAIRVTHETSAHEGQTEAPNIDEKVDLHFIALVHVDGHLYELDAIEVCKKFMERDPDELRFNAIALSAA",
            )
            # confirm that the feature coordinates are correct by extracting
            # the feature sequence from the target sequence and tranlating it.
            cds = feature.extract(self.dna[alignment.target.id]).translate()
            self.assertEqual(feature.qualifiers["translation"][0], cds)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[20873021, 20872472, 20872471, 20872471, 20872390],
                             [       0,      183,      183,      203,      230]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
204	6	0	0	1	20	1	1	+-	CAG33136.1	230	0	230	KI537194	37111980	20872390	20873021	2	183,27,	0,203,	16238959,16239509,
""",
        )
        self.assertRaises(StopIteration, next, alignments)

    def test_writing_psl_35_002(self):
        """Test writing the alignments in psl_35_002.psl."""
        path = "Blat/psl_35_002.psl"
        with open(path) as stream:
            original_data = stream.read()
        alignments = Align.parse(path, "psl")
        stream = StringIO()
        n = Align.write(alignments, stream, "psl")
        self.assertEqual(n, 3)
        stream.seek(0)
        written_data = stream.read()
        stream.close()
        self.assertEqual(original_data, written_data)
        # Convert the alignment to a protein alignment and insert the
        # appropriate sequence data. Write this alignment in a PSL file;
        # the writer will recalculate the values for matches, misMatches,
        # repMatches, and nCount from the sequence data and the alignment.
        #
        # The alignments were generated using
        # blat -t=dnax -q=prot balAcu1.2bit CAG33136.1.fasta psl_35_001.psl
        #
        # To save disk space, we extracted the necessary sequence data using
        #
        # twoBitToFa balAcu1.2bit:KI537979:9712654-9744592 stdout
        # twoBitToFa balAcu1.2bit:KI538594:2103463-2104149 stdout
        # twoBitToFa balAcu1.2bit:KI537194:20872390-20873021 stdout
        #
        # and concatenating the results into file balAcu1.fa. We will use this
        # file below, and create partially defined Seq objects.
        #
        # Load the protein sequence:
        protein = SeqIO.read("Blat/CAG33136.1.fasta", "fasta")
        protein_alignments = []
        alignments = Align.parse(path, "psl")
        for i, alignment in enumerate(alignments):
            records = SeqIO.parse("Blat/balAcu1.fa", "fasta")
            for record in records:
                name, start_end = record.id.split(":")
                if name == alignment.sequences[0].id:
                    break
            else:
                raise Exception("Failed to find DNA sequence")
            start, end = start_end.split("-")
            start = int(start)
            end = int(end)
            length = len(alignment.sequences[0])
            sequence = str(record.seq)
            dna = Seq({start: sequence}, length=length)
            alignment.sequences[0].seq = dna
            self.assertEqual(alignment.sequences[1].id, protein.id)
            alignment.sequences[1].seq = protein.seq
            if i == 0 or i == 1:
                # The alignment is on the forward strand of the DNA sequence:
                self.assertLess(
                    alignment.coordinates[0, 0], alignment.coordinates[0, -1]
                )
            elif i == 2:
                # The alignment is on the reverse strand of the DNA sequence:
                self.assertGreater(
                    alignment.coordinates[0, 0], alignment.coordinates[0, -1]
                )
                # so we take the reverse complement:
                alignment.coordinates[0, :] = len(dna) - alignment.coordinates[0, :]
                alignment.sequences[0].seq = dna.reverse_complement()
            # The protein alignment is always in the forward orientation:
            self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
            # Now extract the aligned sequences:
            aligned_dna = ""
            aligned_protein = ""
            for start, end in alignment.aligned[0]:
                aligned_dna += alignment.sequences[0].seq[start:end]
            for start, end in alignment.aligned[1]:
                aligned_protein += alignment.sequences[1].seq[start:end]
            # Translate the aligned DNA sequence:
            aligned_dna = Seq(aligned_dna)
            aligned_dna_translated = Seq(aligned_dna.translate())
            aligned_protein = Seq(aligned_protein)
            # Create a new alignment including the aligned sequences only:
            records = [
                SeqRecord(aligned_dna_translated, id=alignment.sequences[0].id),
                SeqRecord(aligned_protein, id=alignment.sequences[1].id),
            ]
            coordinates = numpy.array(
                [[0, len(aligned_dna_translated)], [0, len(aligned_protein)]]
            )
            protein_alignment = Alignment(records, coordinates)
            protein_alignments.append(protein_alignment)
            if i == 0:
                self.assertEqual(
                    str(protein_alignment),
                    """\
KI537979          0 QFLKQLGLHPNWQFVDVYGMDPELLSMVPRPVCAVLLLFPITEKYEIFRTEEEEKIKSQG
                  0 ||||||||||||||||||||||||||||||||||||||||||||||.|||||||||||||
CAG33136.         0 QFLKQLGLHPNWQFVDVYGMDPELLSMVPRPVCAVLLLFPITEKYEVFRTEEEEKIKSQG

KI537979         60 QDVTSSVYFMKQTISNACGTIGLIHAIANNKDKMHFESGSTLKKFLEESASMSPEERARY
                 60 |||||||||||||||||||||||||||||||||||||||||||||||||.||||||||||
CAG33136.        60 QDVTSSVYFMKQTISNACGTIGLIHAIANNKDKMHFESGSTLKKFLEESVSMSPEERARY

KI537979        120 LENYDAIRVTHETSAHEGQTEAPNIDEKVDLHFIALVHVDGHLYELDGRKPFPINHGETS
                120 |||||||||||||||||||||||.||||||||||||||||||||||||||||||||||||
CAG33136.       120 LENYDAIRVTHETSAHEGQTEAPSIDEKVDLHFIALVHVDGHLYELDGRKPFPINHGETS

KI537979        180 DETLLEDAIEVCKKFMERDPDELRFNAIALSAA 213
                180 ||||||||||||||||||||||||||||||||| 213
CAG33136.       180 DETLLEDAIEVCKKFMERDPDELRFNAIALSAA 213
""",
                )
            elif i == 1:
                self.assertEqual(
                    str(protein_alignment),
                    """\
KI538594          0 MEGQCWLPLEANPEVTNQLLQLGLHPNWQFVDVYGMDPELLSMVPRPVCAVLLLFPITEK
                  0 ||||.|||||||||||||.|||||||||||||||||||||||||||||||||||||||||
CAG33136.         0 MEGQRWLPLEANPEVTNQFLQLGLHPNWQFVDVYGMDPELLSMVPRPVCAVLLLFPITEK

KI538594         60 YEVFRTEEEEKIKSQGQNITSSGYFMRQTISSACGTIGLIHAIANNKDKMHFESGSTLKK
                 60 |||||||||||||||||..|||.|||.||||.||||||||||||||||||||||||||||
CAG33136.        60 YEVFRTEEEEKIKSQGQDVTSSVYFMKQTISNACGTIGLIHAIANNKDKMHFESGSTLKK

KI538594        120 FLEESASLSPEERAIYLENYDSIRVTHKTSDHEGQTEAQNIDEKVDLHFIALVHVDGHLY
                120 |||||.|.||||||.||||||.|||||.||.|||||||..||||||||||||||||||||
CAG33136.       120 FLEESVSMSPEERARYLENYDAIRVTHETSAHEGQTEAPSIDEKVDLHFIALVHVDGHLY

KI538594        180 ELDGWKPFPINHGETSDATLLRDAIEVFKKFRERDPDERRFNVIALSAA 229
                180 ||||.||||||||||||.|||.|||||.|||.||||||.|||.|||||| 229
CAG33136.       180 ELDGRKPFPINHGETSDETLLEDAIEVCKKFMERDPDELRFNAIALSAA 229
""",
                )
            elif i == 2:
                self.assertEqual(
                    str(protein_alignment),
                    """\
KI537194          0 MESQRWLPLEANPEVTNQFLKQLGLHPNWQCVDVYGMDPELLSMVPRPVCAVLLLFPITE
                  0 ||.|||||||||||||||||||||||||||.|||||||||||||||||||||||||||||
CAG33136.         0 MEGQRWLPLEANPEVTNQFLKQLGLHPNWQFVDVYGMDPELLSMVPRPVCAVLLLFPITE

KI537194         60 KYEIFRTEEEEKTKSQGQDVTSSVYFMKQTISNACGTIGLIHAIANNKDKMHFESGSTLK
                 60 |||.||||||||.|||||||||||||||||||||||||||||||||||||||||||||||
CAG33136.        60 KYEVFRTEEEEKIKSQGQDVTSSVYFMKQTISNACGTIGLIHAIANNKDKMHFESGSTLK

KI537194        120 KFLEESASMSPEERARYLENYDAIRVTHETSAHEGQTEAPNIDEKVDLHFIALVHVDGHL
                120 ||||||.|||||||||||||||||||||||||||||||||.|||||||||||||||||||
CAG33136.       120 KFLEESVSMSPEERARYLENYDAIRVTHETSAHEGQTEAPSIDEKVDLHFIALVHVDGHL

KI537194        180 YELDAIEVCKKFMERDPDELRFNAIALSAA 210
                180 |||||||||||||||||||||||||||||| 210
CAG33136.       180 YELDAIEVCKKFMERDPDELRFNAIALSAA 210
""",
                )
        # Write the protein alignments to a PSL file:
        stream = StringIO()
        n = Align.write(protein_alignments, stream, "psl", wildcard="X")
        self.assertEqual(n, 3)
        # Read the alignments back in:
        alignments = Align.parse(path, "psl")
        stream.seek(0)
        protein_alignments = Align.parse(stream, "psl")
        for alignment, protein_alignment in zip(alignments, protein_alignments):
            # Confirm that the recalculated values for matches, misMatches,
            # repMatches, and nCount are correct:
            self.assertEqual(alignment.matches, protein_alignment.matches)
            self.assertEqual(alignment.misMatches, protein_alignment.misMatches)
            self.assertEqual(alignment.repMatches, protein_alignment.repMatches)
            self.assertEqual(alignment.nCount, protein_alignment.nCount)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
