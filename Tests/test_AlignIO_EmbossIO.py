# Copyright 2008-2014 by Peter Cock.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.AlignIO.EmbossIO module."""
import unittest

from io import StringIO

from Bio.AlignIO.EmbossIO import EmbossIterator

# http://emboss.sourceforge.net/docs/themes/alnformats/align.simple
simple_example = """\
########################################
# Program:  alignret
# Rundate:  Wed Jan 16 17:16:13 2002
# Report_file: stdout
########################################
#=======================================
#
# Aligned_sequences: 4
# 1: IXI_234
# 2: IXI_235
# 3: IXI_236
# 4: IXI_237
# Matrix: EBLOSUM62
# Gap_penalty: 10.0
# Extend_penalty: 0.5
#
# Length: 131
# Identity:      95/131 (72.5%)
# Similarity:   127/131 (96.9%)
# Gaps:          25/131 (19.1%)
# Score: 100.0
#
#
#=======================================

IXI_234            1 TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQAT     50
IXI_235            1 TSPASIRPPAGPSSR---------RPSPPGPRRPTGRPCCSAAPRRPQAT     41
IXI_236            1 TSPASIRPPAGPSSRPAMVSSR--RPSPPPPRRPPGRPCCSAAPPRPQAT     48
IXI_237            1 TSPASLRPPAGPSSRPAMVSSRR-RPSPPGPRRPT----CSAAPRRPQAT     45
                     |||||:|||||||||:::::::  |||||:||||:::::|||||:|||||

IXI_234           51 GGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAG    100
IXI_235           42 GGWKTCSGTCTTSTSTRHRGRSGW----------RASRKSMRAACSRSAG     81
IXI_236           49 GGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSR--G     96
IXI_237           46 GGYKTCSGTCTTSTSTRHRGRSGYSARTTTAACLRASRKSMRAACSR--G     93
                     ||:||||||||||||||||||||:::::::::::|||||||||||||  |

IXI_234          101 SRPNRFAPTLMSSCITSTTGPPAWAGDRSHE    131
IXI_235           82 SRPNRFAPTLMSSCITSTTGPPAWAGDRSHE    112
IXI_236           97 SRPPRFAPPLMSSCITSTTGPPPPAGDRSHE    127
IXI_237           94 SRPNRFAPTLMSSCLTSTTGPPAYAGDRSHE    124
                     |||:||||:|||||:|||||||::|||||||


#---------------------------------------
#---------------------------------------

"""  # noqa: E122 not clear to me, why this comes up here

# http://emboss.sourceforge.net/docs/themes/alnformats/align.pair
with open("Emboss/water.txt") as handle:
    pair_example = handle.read()

with open("Emboss/needle.txt") as handle:
    pair_example2 = handle.read()

with open("Emboss/needle_overhang.txt") as handle:
    pair_example3 = handle.read()


class TestEmbossIO(unittest.TestCase):
    def test_pair_example(self):
        alignments = list(EmbossIterator(StringIO(pair_example)))
        self.assertEqual(len(alignments), 1)
        self.assertEqual(len(alignments[0]), 2)
        self.assertEqual([r.id for r in alignments[0]], ["IXI_234", "IXI_235"])

    def test_simple_example(self):
        alignments = list(EmbossIterator(StringIO(simple_example)))
        self.assertEqual(len(alignments), 1)
        self.assertEqual(len(alignments[0]), 4)
        self.assertEqual(
            [r.id for r in alignments[0]], ["IXI_234", "IXI_235", "IXI_236", "IXI_237"]
        )

    def test_pair_plus_simple(self):
        alignments = list(EmbossIterator(StringIO(pair_example + simple_example)))
        self.assertEqual(len(alignments), 2)
        self.assertEqual(len(alignments[0]), 2)
        self.assertEqual(len(alignments[1]), 4)
        self.assertEqual([r.id for r in alignments[0]], ["IXI_234", "IXI_235"])
        self.assertEqual(
            [r.id for r in alignments[1]], ["IXI_234", "IXI_235", "IXI_236", "IXI_237"]
        )

    def test_pair_example2(self):
        alignments = list(EmbossIterator(StringIO(pair_example2)))
        self.assertEqual(len(alignments), 5)
        self.assertEqual(len(alignments[0]), 2)
        self.assertEqual(
            [r.id for r in alignments[0]], ["ref_rec", "gi|94968718|receiver"]
        )
        self.assertEqual(
            [r.id for r in alignments[4]], ["ref_rec", "gi|94970041|receiver"]
        )

    def test_pair_example3(self):
        alignments = list(EmbossIterator(StringIO(pair_example3)))
        self.assertEqual(len(alignments), 1)
        self.assertEqual(len(alignments[0]), 2)
        self.assertEqual([r.id for r in alignments[0]], ["asis", "asis"])


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
