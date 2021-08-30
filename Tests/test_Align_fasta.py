# Copyright 2021 by Michiel de Hoon.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.Align.fasta module."""
import unittest

from Bio.Align.fasta_m8 import AlignmentIterator

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align.emboss."
    ) from None


class TestFasta(unittest.TestCase):

    def test_m8CB(self):
        # Alignment file obtained by running
        # fasta36 -q -m 8CB seq/mgstm1.aa seq/prot_test.lseg
        # in the fasta36 source distribution
        path = "Fasta/output_m8CB.txt"
        with open(path) as stream:
            alignments = AlignmentIterator(stream)
            self.assertEqual(alignments.commandline, "fasta36 -q -m 8CB seq/mgstm1.aa seq/prot_test.lseg")
            alignments = list(alignments)
        return
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertEqual(alignment.annotations["matrix"], "EBLOSUM62")
        self.assertAlmostEqual(alignment.annotations["gap_penalty"], 10.0)
        self.assertAlmostEqual(alignment.annotations["extend_penalty"], 0.5)
        self.assertEqual(alignment.annotations["identity"], 112)
        self.assertEqual(alignment.annotations["similarity"], 112)
        self.assertEqual(alignment.annotations["gaps"], 19)
        self.assertAlmostEqual(alignment.annotations["score"], 591.5)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 131))
        self.assertEqual(alignment.sequences[0].id, "IXI_234")
        self.assertEqual(alignment.sequences[1].id, "IXI_235")
        self.assertEqual(
            alignment.sequences[0].seq,
            "TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "TSPASIRPPAGPSSRRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[0, 15, 24, 74, 84, 131], [0, 15, 15, 65, 65, 112]]),
            )
        )
        self.assertEqual(
            alignment[0],
            "TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertEqual(
            alignment[1],
            "TSPASIRPPAGPSSR---------RPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGW----------RASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            "|||||||||||||||         ||||||||||||||||||||||||||||||||||||||||||||||||||          |||||||||||||||||||||||||||||||||||||||||||||||",
        )

    def test_m8CC(self):
        # Alignment file obtained by running
        # fasta36 -q -m 8CB seq/mgstm1.aa seq/prot_test.lseg
        # in the fasta36 source distribution
        path = "Fasta/output_m8CC.txt"
        with open(path) as stream:
            alignments = AlignmentIterator(stream)
            self.assertEqual(alignments.commandline, "fasta36 -q -m 8CC seq/mgstm1.aa seq/prot_test.lseg")
            alignments = list(alignments)
        return
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertEqual(alignment.annotations["matrix"], "EBLOSUM62")
        self.assertAlmostEqual(alignment.annotations["gap_penalty"], 10.0)
        self.assertAlmostEqual(alignment.annotations["extend_penalty"], 0.5)
        self.assertEqual(alignment.annotations["identity"], 112)
        self.assertEqual(alignment.annotations["similarity"], 112)
        self.assertEqual(alignment.annotations["gaps"], 19)
        self.assertAlmostEqual(alignment.annotations["score"], 591.5)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 131))
        self.assertEqual(alignment.sequences[0].id, "IXI_234")
        self.assertEqual(alignment.sequences[1].id, "IXI_235")
        self.assertEqual(
            alignment.sequences[0].seq,
            "TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "TSPASIRPPAGPSSRRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[0, 15, 24, 74, 84, 131], [0, 15, 15, 65, 65, 112]]),
            )
        )
        self.assertEqual(
            alignment[0],
            "TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertEqual(
            alignment[1],
            "TSPASIRPPAGPSSR---------RPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGW----------RASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            "|||||||||||||||         ||||||||||||||||||||||||||||||||||||||||||||||||||          |||||||||||||||||||||||||||||||||||||||||||||||",
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
