# Copyright 2008-2014 by Peter Cock.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.Align.emboss module."""
import unittest

from io import StringIO

from Bio.Align.emboss import AlignmentIterator

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align.emboss."
    ) from None


with open("Emboss/needle.txt") as handle:
    pair_example2 = handle.read()


with open("Emboss/needle_overhang.txt") as handle:
    pair_example3 = handle.read()


class TestEmboss(unittest.TestCase):
    def test_pair_example(self):
        # Alignment file obtained from EMBOSS:
        # http://emboss.sourceforge.net/docs/themes/alnformats/align.pair
        path = "Emboss/water.txt"
        with open(path) as stream:
            alignments = AlignmentIterator(stream)
            self.assertEqual(alignments.program, "water")
            self.assertEqual(alignments.rundate, "Wed Jan 16 17:23:19 2002")
            self.assertEqual(alignments.report_file, "stdout")
            alignments = list(alignments)
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

    def test_pair_example_nobrief(self):
        # Variation on the alignment file obtained from EMBOSS
        # (http://emboss.sourceforge.net/docs/themes/alnformats/align.pair)
        # if we include 3 sequences to align against, and we use the -nobrief
        # command line option.
        path = "Emboss/needle_nobrief_multiple.pair"
        with open(path) as stream:
            alignments = AlignmentIterator(stream)
            self.assertEqual(alignments.program, "needle")
            self.assertEqual(alignments.rundate, "Fri 23 Jul 2021 22:45:41")
            self.assertEqual(
                alignments.commandline,
                "needle -asequence seqa.fa -bsequence seqb.fa -datafile EBLOSUM62 -gapopen 10 -gapextend 0.5 -nobrief -outfile stdout",
            )
            self.assertEqual(alignments.align_format, "srspair")
            self.assertEqual(alignments.report_file, "stdout")
            alignments = list(alignments)
        self.assertEqual(len(alignments), 3)
        alignment = alignments[0]
        self.assertEqual(alignment.annotations["matrix"], "EBLOSUM62")
        self.assertAlmostEqual(alignment.annotations["gap_penalty"], 10.0)
        self.assertAlmostEqual(alignment.annotations["extend_penalty"], 0.5)
        self.assertEqual(alignment.annotations["identity"], 112)
        self.assertEqual(alignment.annotations["similarity"], 112)
        self.assertEqual(alignment.annotations["gaps"], 19)
        self.assertAlmostEqual(alignment.annotations["score"], 591.5)
        self.assertEqual(alignment.annotations["longest_identity"], "100.00%")
        self.assertEqual(alignment.annotations["longest_similarity"], "100.00%")
        self.assertEqual(alignment.annotations["shortest_identity"], "85.50%")
        self.assertEqual(alignment.annotations["shortest_similarity"], "85.50%")
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
        alignment = alignments[1]
        self.assertEqual(alignment.annotations["matrix"], "EBLOSUM62")
        self.assertAlmostEqual(alignment.annotations["gap_penalty"], 10.0)
        self.assertAlmostEqual(alignment.annotations["extend_penalty"], 0.5)
        self.assertEqual(alignment.annotations["identity"], 120)
        self.assertEqual(alignment.annotations["similarity"], 120)
        self.assertEqual(alignment.annotations["gaps"], 4)
        self.assertAlmostEqual(alignment.annotations["score"], 618.0)
        self.assertEqual(alignment.annotations["longest_identity"], "94.49%")
        self.assertEqual(alignment.annotations["longest_similarity"], "94.49%")
        self.assertEqual(alignment.annotations["shortest_identity"], "91.60%")
        self.assertEqual(alignment.annotations["shortest_similarity"], "91.60%")
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 131))
        self.assertEqual(alignment.sequences[0].id, "IXI_234")
        self.assertEqual(alignment.sequences[1].id, "IXI_236")
        self.assertEqual(
            alignment.sequences[0].seq,
            "TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "TSPASIRPPAGPSSRPAMVSSRRPSPPPPRRPPGRPCCSAAPPRPQATGGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRGSRPPRFAPPLMSSCITSTTGPPPPAGDRSHE",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[0, 22, 24, 97, 99, 131], [0, 22, 22, 95, 95, 127]]),
            )
        )
        self.assertEqual(
            alignment[0],
            "TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertEqual(
            alignment[1],
            "TSPASIRPPAGPSSRPAMVSSR--RPSPPPPRRPPGRPCCSAAPPRPQATGGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSR--GSRPPRFAPPLMSSCITSTTGPPPPAGDRSHE",
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            "||||||||||||||||||||||  |||||.||||.|||||||||.||||||||||||||||||||||||||||||||||||||||||||||||||||  ||||.||||.|||||||||||||..|||||||",
        )
        alignment = alignments[2]
        self.assertEqual(alignment.annotations["matrix"], "EBLOSUM62")
        self.assertAlmostEqual(alignment.annotations["gap_penalty"], 10.0)
        self.assertAlmostEqual(alignment.annotations["extend_penalty"], 0.5)
        self.assertEqual(alignment.annotations["identity"], 119)
        self.assertEqual(alignment.annotations["similarity"], 124)
        self.assertEqual(alignment.annotations["gaps"], 7)
        self.assertAlmostEqual(alignment.annotations["score"], 609.0)
        self.assertEqual(alignment.annotations["longest_identity"], "95.97%")
        self.assertEqual(alignment.annotations["longest_similarity"], "100.00%")
        self.assertEqual(alignment.annotations["shortest_identity"], "90.84%")
        self.assertEqual(alignment.annotations["shortest_similarity"], "94.66%")
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 131))
        self.assertEqual(alignment.sequences[0].id, "IXI_234")
        self.assertEqual(alignment.sequences[1].id, "IXI_237")
        self.assertEqual(
            alignment.sequences[0].seq,
            "TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "TSPASLRPPAGPSSRPAMVSSRRRPSPPGPRRPTCSAAPRRPQATGGYKTCSGTCTTSTSTRHRGRSGYSARTTTAACLRASRKSMRAACSRGSRPNRFAPTLMSSCLTSTTGPPAYAGDRSHE",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [[0, 23, 24, 35, 39, 97, 99, 131], [0, 23, 23, 34, 34, 92, 92, 124]]
                ),
            )
        )
        self.assertEqual(
            alignment[0],
            "TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertEqual(
            alignment[1],
            "TSPASLRPPAGPSSRPAMVSSRR-RPSPPGPRRPT----CSAAPRRPQATGGYKTCSGTCTTSTSTRHRGRSGYSARTTTAACLRASRKSMRAACSR--GSRPNRFAPTLMSSCLTSTTGPPAYAGDRSHE",
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            "|||||:||||||||||||||||| |||||||||||    |||||||||||||:||||||||||||||||||||:|||||||||||||||||||||||  |||||||||||||||:||||||||:|||||||",
        )

    def test_simple_example(self):
        # Alignment file obtained from EMBOSS:
        # http://emboss.sourceforge.net/docs/themes/alnformats/align.simple
        path = "Emboss/alignret.txt"
        with open(path) as stream:
            alignments = AlignmentIterator(stream)
            self.assertEqual(alignments.program, "alignret")
            self.assertEqual(alignments.rundate, "Wed Jan 16 17:16:13 2002")
            self.assertEqual(alignments.report_file, "stdout")
            alignments = list(alignments)
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertEqual(alignment.annotations["matrix"], "EBLOSUM62")
        self.assertAlmostEqual(alignment.annotations["gap_penalty"], 10.0)
        self.assertAlmostEqual(alignment.annotations["extend_penalty"], 0.5)
        self.assertEqual(alignment.annotations["identity"], 95)
        self.assertEqual(alignment.annotations["similarity"], 127)
        self.assertEqual(alignment.annotations["gaps"], 25)
        self.assertAlmostEqual(alignment.annotations["score"], 100.0)
        self.assertEqual(len(alignment), 4)
        self.assertEqual(alignment.shape, (4, 131))
        self.assertEqual(alignment.sequences[0].id, "IXI_234")
        self.assertEqual(alignment.sequences[1].id, "IXI_235")
        self.assertEqual(alignment.sequences[2].id, "IXI_236")
        self.assertEqual(alignment.sequences[3].id, "IXI_237")
        self.assertEqual(
            alignment.sequences[0].seq,
            "TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "TSPASIRPPAGPSSRRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertEqual(
            alignment.sequences[2].seq,
            "TSPASIRPPAGPSSRPAMVSSRRPSPPPPRRPPGRPCCSAAPPRPQATGGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRGSRPPRFAPPLMSSCITSTTGPPPPAGDRSHE",
        )
        self.assertEqual(
            alignment.sequences[3].seq,
            "TSPASLRPPAGPSSRPAMVSSRRRPSPPGPRRPTCSAAPRRPQATGGYKTCSGTCTTSTSTRHRGRSGYSARTTTAACLRASRKSMRAACSRGSRPNRFAPTLMSSCLTSTTGPPAYAGDRSHE",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [
                        [0, 15, 22, 23, 24, 35, 39, 74, 84, 97, 99, 131],
                        [0, 15, 15, 15, 15, 26, 30, 65, 65, 78, 80, 112],
                        [0, 15, 22, 22, 22, 33, 37, 72, 82, 95, 95, 127],
                        [0, 15, 22, 23, 23, 34, 34, 69, 79, 92, 92, 124],
                    ]
                ),
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
            alignment[2],
            "TSPASIRPPAGPSSRPAMVSSR--RPSPPPPRRPPGRPCCSAAPPRPQATGGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSR--GSRPPRFAPPLMSSCITSTTGPPPPAGDRSHE",
        )
        self.assertEqual(
            alignment[3],
            "TSPASLRPPAGPSSRPAMVSSRR-RPSPPGPRRPT----CSAAPRRPQATGGYKTCSGTCTTSTSTRHRGRSGYSARTTTAACLRASRKSMRAACSR--GSRPNRFAPTLMSSCLTSTTGPPAYAGDRSHE",
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            "|||||:|||||||||:::::::  |||||:||||:::::|||||:|||||||:||||||||||||||||||||:::::::::::|||||||||||||  ||||:||||:|||||:|||||||::|||||||",
        )

    def test_pair_example2(self):
        alignments = AlignmentIterator(StringIO(pair_example2))
        self.assertEqual(alignments.program, "needle")
        self.assertEqual(alignments.rundate, "Sun 27 Apr 2007 17:20:35")
        self.assertEqual(
            alignments.commandline,
            "needle [-asequence] Spo0F.faa [-bsequence] paired_r.faa -sformat2 pearson",
        )
        self.assertEqual(alignments.align_format, "srspair")
        self.assertEqual(alignments.report_file, "ref_rec .needle")
        alignments = list(alignments)
        self.assertEqual(len(alignments), 5)
        alignment = alignments[0]
        self.assertEqual(alignment.annotations["matrix"], "EBLOSUM62")
        self.assertAlmostEqual(alignment.annotations["gap_penalty"], 10.0)
        self.assertAlmostEqual(alignment.annotations["extend_penalty"], 0.5)
        self.assertEqual(alignment.annotations["identity"], 32)
        self.assertEqual(alignment.annotations["similarity"], 64)
        self.assertEqual(alignment.annotations["gaps"], 17)
        self.assertAlmostEqual(alignment.annotations["score"], 112.0)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 124))
        self.assertEqual(alignment.sequences[0].id, "ref_rec")
        self.assertEqual(alignment.sequences[1].id, "gi|94968718|receiver")
        self.assertEqual(
            alignment.sequences[0].seq,
            "KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIVTKERPDLVLLDMKIPGMDGIEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHFAKPFDIDEIRDAV",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "VLLADDHALVRRGFRLMLEDDPEIEIVAEAGDGAQAVKLAGELHPRVVVMDCAMPGMSGMDATKQIRTQWPDIAVLMLTMHSEDTWVRLALEAGANGYILKSAIDLDLIQAVRRVANGET",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [
                        [0, 1, 7, 7, 17, 19, 100, 100, 108, 109, 111, 111],
                        [0, 0, 6, 10, 20, 20, 101, 102, 110, 110, 112, 120],
                    ]
                ),
            )
        )
        self.assertEqual(
            alignment[0],
            "KILIVDD----QYGIRILLNEVFNKEGYQTFQAANGLQALDIVTKERPDLVLLDMKIPGMDGIEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHFAK-PFDIDEIRDAV--------",
        )
        self.assertEqual(
            alignment[1],
            "-VLLADDHALVRRGFRLMLED--DPEIEIVAEAGDGAQAVKLAGELHPRVVVMDCAMPGMSGMDATKQIRTQWPDIAVLMLTMHSEDTWVRLALEAGANGYILKSAIDLDLIQ-AVRRVANGET",
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            " :|:.||    :.|.|::|.:  :.|.....:|.:|.||:.:..:..|.:|::|..:|||.|::..|:::....:|.|:::|.:.|...::.:.|.||..:..| ..|:|.|: ||        ",
        )
        alignment = alignments[1]
        self.assertEqual(alignment.annotations["matrix"], "EBLOSUM62")
        self.assertAlmostEqual(alignment.annotations["gap_penalty"], 10.0)
        self.assertAlmostEqual(alignment.annotations["extend_penalty"], 0.5)
        self.assertEqual(alignment.annotations["identity"], 34)
        self.assertEqual(alignment.annotations["similarity"], 58)
        self.assertEqual(alignment.annotations["gaps"], 9)
        self.assertAlmostEqual(alignment.annotations["score"], 154.0)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 119))
        self.assertEqual(alignment.sequences[0].id, "ref_rec")
        self.assertEqual(alignment.sequences[1].id, "gi|94968761|receiver")
        self.assertEqual(
            alignment.sequences[0].seq,
            "KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIVTKERPDLVLLDMKIPGMDGIEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHFAKPFDIDEIRDAV",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "ILIVDDEANTLASLSRAFRLAGHEATVCDNAVRALEIAKSKPFDLILSDVVMPGRDGLTLLEDLKTAGVQAPVVMMSGQAHIEMAVKATRLGALDFLEKPLSTDKLLLTVENALKLKR",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[0, 1, 111, 111], [0, 0, 110, 118]])
            )
        )
        self.assertEqual(
            alignment[0],
            "KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIVTKERPDLVLLDMKIPGMDGIEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHFAKPFDIDEIRDAV--------",
        )
        self.assertEqual(
            alignment[1],
            "-ILIVDDEANTLASLSRAFRLAGHEATVCDNAVRALEIAKSKPFDLILSDVVMPGRDGLTLLEDLKTAGVQAPVVMMSGQAHIEMAVKATRLGALDFLEKPLSTDKLLLTVENALKLKR",
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            " ||||||:......|:..|...|::.....|.::||:|...:..||:|.|:.:||.||:.:|:.:|.......|::|:....::|..::..||||....||...|::...|        ",
        )
        alignment = alignments[2]
        self.assertEqual(alignment.annotations["matrix"], "EBLOSUM62")
        self.assertAlmostEqual(alignment.annotations["gap_penalty"], 10.0)
        self.assertAlmostEqual(alignment.annotations["extend_penalty"], 0.5)
        self.assertEqual(alignment.annotations["identity"], 29)
        self.assertEqual(alignment.annotations["similarity"], 53)
        self.assertEqual(alignment.annotations["gaps"], 9)
        self.assertAlmostEqual(alignment.annotations["score"], 121.0)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 120))
        self.assertEqual(alignment.sequences[0].id, "ref_rec")
        self.assertEqual(alignment.sequences[1].id, "gi|94967506|receiver")
        self.assertEqual(
            alignment.sequences[0].seq,
            "KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIVTKERPDLVLLDMKIPGMDGIEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHFAKPFDIDEIRDAV",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "LHIVVVDDDPGTCVYIESVFAELGHTCKSFVRPEAAEEYILTHPVDLAIVDVYLGSTTGVEVLRRCRVHRPKLYAVIITGQISLEMAARSIAEGAVDYIQKPIDIDALLNIAERALEHKE",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[0, 0, 111, 111], [0, 1, 112, 120]])
            )
        )
        self.assertEqual(
            alignment[0],
            "-KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIVTKERPDLVLLDMKIPGMDGIEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHFAKPFDIDEIRDAV--------",
        )
        self.assertEqual(
            alignment[1],
            "LHIVVVDDDPGTCVYIESVFAELGHTCKSFVRPEAAEEYILTHPVDLAIVDVYLGSTTGVEVLRRCRVHRPKLYAVIITGQISLEMAARSIAEGAVDYIQKPIDIDALLNIAERALEHKE",
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            " .|::|||..|..:.:..||.:.|:..........|.:.:.....||.::|:.:....|:|:|:|.:|....:..:|:|....|:|...|...||:.:..||.|||.:.:..        ",
        )
        alignment = alignments[3]
        self.assertEqual(alignment.annotations["matrix"], "EBLOSUM62")
        self.assertAlmostEqual(alignment.annotations["gap_penalty"], 10.0)
        self.assertAlmostEqual(alignment.annotations["extend_penalty"], 0.5)
        self.assertEqual(alignment.annotations["identity"], 30)
        self.assertEqual(alignment.annotations["similarity"], 64)
        self.assertEqual(alignment.annotations["gaps"], 9)
        self.assertAlmostEqual(alignment.annotations["score"], 126.0)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 118))
        self.assertEqual(alignment.sequences[0].id, "ref_rec")
        self.assertEqual(alignment.sequences[1].id, "gi|94970045|receiver")
        self.assertEqual(
            alignment.sequences[0].seq,
            "KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIVTKERPDLVLLDMKIPGMDGIEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHFAKPFDIDEIRDAV",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "VLLVEDEEALRAAAGDFLETRGYKIMTARDGTEALSMASKFAERIDVLITDLVMPGISGRVLAQELVKIHPETKVMYMSGYDDETVMVNGEIDSSSAFLRKPFRMDALSAKIREVL",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [
                        [0, 1, 41, 41, 82, 83, 98, 98, 105, 105, 111],
                        [0, 0, 40, 42, 83, 83, 98, 99, 106, 110, 116],
                    ]
                ),
            )
        )
        self.assertEqual(
            alignment[0],
            "KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIVTK--ERPDLVLLDMKIPGMDGIEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHF-AKPFDID----EIRDAV",
        )
        self.assertEqual(
            alignment[1],
            "-VLLVEDEEALRAAAGDFLETRGYKIMTARDGTEALSMASKFAERIDVLITDLVMPGISGRVLAQELVKIHPETKVMYMSGYDD-ETVMVNGEIDSSSAFLRKPFRMDALSAKIREVL",
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            " :|:|:|:..:|....:.....||:...|.:|.:||.:.:|  ||.|:::.|:.:||:.|..:.:.:..|....:|:.|:.|.: :.:..:.|:.:.:.| .|||.:|    :||:.:",
        )
        alignment = alignments[4]
        self.assertEqual(alignment.annotations["matrix"], "EBLOSUM62")
        self.assertAlmostEqual(alignment.annotations["gap_penalty"], 10.0)
        self.assertAlmostEqual(alignment.annotations["extend_penalty"], 0.5)
        self.assertEqual(alignment.annotations["identity"], 35)
        self.assertEqual(alignment.annotations["similarity"], 70)
        self.assertEqual(alignment.annotations["gaps"], 18)
        self.assertAlmostEqual(alignment.annotations["score"], 156.5)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 125))
        self.assertEqual(alignment.sequences[0].id, "ref_rec")
        self.assertEqual(alignment.sequences[1].id, "gi|94970041|receiver")
        self.assertEqual(
            alignment.sequences[0].seq,
            "KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIVTKERPDLVLLDMKIPGMDGIEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHFAKPFDIDEIRDAV",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "TVLLVEDEEGVRKLVRGILSRQGYHVLEATSGEEALEIVRESTQKIDMLLSDVVLVGMSGRELSERLRIQMPSLKVIYMSGYTDDAIVRHGVLTESAEFLQKPFTSDSLLRKVRAVLQKRQ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [
                        [0, 39, 39, 88, 92, 99, 99, 111, 111],
                        [0, 39, 41, 90, 90, 97, 101, 113, 121],
                    ]
                ),
            )
        )
        self.assertEqual(
            alignment[0],
            "KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIV--TKERPDLVLLDMKIPGMDGIEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHFA----KPFDIDEIRDAV--------",
        )
        self.assertEqual(
            alignment[1],
            "TVLLVEDEEGVRKLVRGILSRQGYHVLEATSGEEALEIVRESTQKIDMLLSDVVLVGMSGRELSERLRIQMPSLKVIYMSGYTDDAIVRH----GVLTESAEFLQKPFTSDSLLRKVRAVLQKRQ",
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            ".:|:|:|:.|:|.|:..:.:::||...:|.:|.:||:||  :.::.|::|.|:.:.||.|.|:.:|:::...:::||.|:.|.:..:::.    |.||..|    |||..|.:...|        ",
        )

    def test_pair_example3(self):
        alignments = AlignmentIterator(StringIO(pair_example3))
        self.assertEqual(alignments.program, "needle")
        self.assertEqual(alignments.rundate, "Mon 14 Jul 2008 11:45:42")
        self.assertEqual(
            alignments.commandline,
            "needle [-asequence] asis:TGTGGTTAGGTTTGGTTTTATTGGGGGCTTGGTTTGGGCCCACCCCAAATAGGGAGTGGGGGTATGACCTCAGATAGACGAGCTTATTTTAGGGCGGCGACTATAATTATTTCGTTTCCTACAAGGATTAAAGTTTTTTCTTTTACTGTGGGAGGGGGTTTGGTATTAAGAAACGCTAGTCCGGATGTGGCTCTCCATGATACTTATTGTGTAGTAGCTCATTTTCATTATGTTCTTCGAATGGGAGCAGTCATTGGTATTTTTTTGGTTTTTTTTTGAAATTTTTAGGTTATTTAGACCATTTTTTTTTGTTTCGCTAATTAGAATTTTATTAGCCTTTGGTTTTTTTTTATTTTTTGGGGTTAAGACAAGGTGTCGTTGAATTAGTTTAGCAAAATACTGCTTAAGGTAGGCTATAGGATCTACCTTTTATCTTTCTAATCTTTTGTTTTAGTATAATTGGTCTTCGATTCAACAATTTTTAGTCTTCAGTCTTTTTTTTTATTTTGAAAAGGTTTTAACACTCTTGGTTTTGGAGGCTTTGGCTTTCTTCTTACTCTTAGGAGGATGGGCGCTAGAAAGAGTTTTAAGAGGGTGTGAAAGGGGGTTAATAGC [-bsequence] asis:TTATTAATCTTATGGTTTTGCCGTAAAATTTCTTTCTTTATTTTTTATTGTTAGGATTTTGTTGATTTTATTTTTCTCAAGAATTTTTAGGTCAATTAGACCGGCTTATTTTTTTGTCAGTGTTTAAAGTTTTATTAATTTTTGGGGGGGGGGGGAGACGGGGTGTTATCTGAATTAGTTTTTGGGAGTCTCTAGACATCTCATGGGTTGGCCGGGGGCCTGCCGTCTATAGTTCTTATTCCTTTTAAGGGAGTAAGAATTTCGATTCAGCAACTTTAGTTCACAGTCTTTTTTTTTATTAAGAAAGGTTT -filter",
        )
        self.assertEqual(alignments.align_format, "srspair")
        self.assertEqual(alignments.report_file, "stdout")
        alignments = list(alignments)
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertEqual(alignment.annotations["matrix"], "EDNAFULL")
        self.assertAlmostEqual(alignment.annotations["gap_penalty"], 10.0)
        self.assertAlmostEqual(alignment.annotations["extend_penalty"], 0.5)
        self.assertEqual(alignment.annotations["identity"], 210)
        self.assertEqual(alignment.annotations["similarity"], 210)
        self.assertEqual(alignment.annotations["gaps"], 408)
        self.assertAlmostEqual(alignment.annotations["score"], 561.0)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 667))
        self.assertEqual(alignment.sequences[0].id, "asis")
        self.assertEqual(alignment.sequences[1].id, "asis")
        self.assertEqual(
            alignment.sequences[0].seq,
            "TGTGGTTAGGTTTGGTTTTATTGGGGGCTTGGTTTGGGCCCACCCCAAATAGGGAGTGGGGGTATGACCTCAGATAGACGAGCTTATTTTAGGGCGGCGACTATAATTATTTCGTTTCCTACAAGGATTAAAGTTTTTTCTTTTACTGTGGGAGGGGGTTTGGTATTAAGAAACGCTAGTCCGGATGTGGCTCTCCATGATACTTATTGTGTAGTAGCTCATTTTCATTATGTTCTTCGAATGGGAGCAGTCATTGGTATTTTTTTGGTTTTTTTTTGAAATTTTTAGGTTATTTAGACCATTTTTTTTTGTTTCGCTAATTAGAATTTTATTAGCCTTTGGTTTTTTTTTATTTTTTGGGGTTAAGACAAGGTGTCGTTGAATTAGTTTAGCAAAATACTGCTTAAGGTAGGCTATAGGATCTACCTTTTATCTTTCTAATCTTTTGTTTTAGTATAATTGGTCTTCGATTCAACAATTTTTAGTCTTCAGTCTTTTTTTTTATTTTGAAAAGGTTTTAACACTCTTGGTTTTGGAGGCTTTGGCTTTCTTCTTACTCTTAGGAGGATGGGCGCTAGAAAGAGTTTTAAGAGGGTGTGAAAGGGGGTTAATAGC",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "TTATTAATCTTATGGTTTTGCCGTAAAATTTCTTTCTTTATTTTTTATTGTTAGGATTTTGTTGATTTTATTTTTCTCAAGAATTTTTAGGTCAATTAGACCGGCTTATTTTTTTGTCAGTGTTTAAAGTTTTATTAATTTTTGGGGGGGGGGGGAGACGGGGTGTTATCTGAATTAGTTTTTGGGAGTCTCTAGACATCTCATGGGTTGGCCGGGGGCCTGCCGTCTATAGTTCTTATTCCTTTTAAGGGAGTAAGAATTTCGATTCAGCAACTTTAGTTCACAGTCTTTTTTTTTATTAAGAAAGGTTT",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [
                        [
                            0,
                            162,
                            169,
                            201,
                            210,
                            210,
                            220,
                            222,
                            236,
                            240,
                            244,
                            253,
                            277,
                            277,
                            278,
                            279,
                            300,
                            300,
                            310,
                            310,
                            314,
                            320,
                            334,
                            351,
                            357,
                            357,
                            379,
                            379,
                            390,
                            403,
                            405,
                            407,
                            418,
                            418,
                            442,
                            442,
                            447,
                            447,
                            448,
                            452,
                            455,
                            455,
                            460,
                            465,
                            478,
                            479,
                            509,
                            510,
                            518,
                            615,
                        ],
                        [
                            0,
                            0,
                            7,
                            7,
                            16,
                            22,
                            32,
                            32,
                            46,
                            46,
                            50,
                            50,
                            74,
                            80,
                            81,
                            81,
                            102,
                            107,
                            117,
                            119,
                            123,
                            123,
                            137,
                            137,
                            143,
                            147,
                            169,
                            170,
                            181,
                            181,
                            183,
                            183,
                            194,
                            215,
                            239,
                            241,
                            246,
                            250,
                            251,
                            251,
                            254,
                            255,
                            260,
                            260,
                            273,
                            273,
                            303,
                            303,
                            311,
                            311,
                        ],
                    ]
                ),
            )
        )
        self.assertEqual(
            alignment[0],
            "TGTGGTTAGGTTTGGTTTTATTGGGGGCTTGGTTTGGGCCCACCCCAAATAGGGAGTGGGGGTATGACCTCAGATAGACGAGCTTATTTTAGGGCGGCGACTATAATTATTTCGTTTCCTACAAGGATTAAAGTTTTTTCTTTTACTGTGGGAGGGGGTTTGGTATTAAGAAACGCTAGTCCGGATGTGGCTCTCCATGATACTTATTGT------GTAGTAGCTCATTTTCATTATGTTCTTCGAATGGGAGCAGTCATTGGTATTTTTTTGGTTTTTTTTT------GAAATTTTTAGGTTATTTAGACC-----ATTTTTTTTT--GTTTCGCTAATTAGAATTTTATTAGCCTTTGGTTTTTTTTTATTTTT----TGGGGTTAAGACAAGGTGTCGT-TGAATTAGTTTAGCAAAATACTGCTTAAGGTAGGCTATA---------------------GGATCTACCTTTTATCTTTCTAAT--CTTTT----GTTTTAGT-ATAATTGGTCTTCGATTCAACAATTTTTAGTCTTCAGTCTTTTTTTTTATTTTGAAAAGGTTTTAACACTCTTGGTTTTGGAGGCTTTGGCTTTCTTCTTACTCTTAGGAGGATGGGCGCTAGAAAGAGTTTTAAGAGGGTGTGAAAGGGGGTTAATAGC",
        )
        self.assertEqual(
            alignment[1],
            "------------------------------------------------------------------------------------------------------------------------------------------------------------------TTATTAA--------------------------------TCTTATGGTTTTGCCGTAAAATTTC--TTTCTTTATTTTTT----ATTG---------TTAGGATTTTGTTGATTTTATTTTTCTCAAG-AATTTTTAGGTCAATTAGACCGGCTTATTTTTTTGTCAGTGT------TTAAAGTTTTATTA-----------------ATTTTTGGGGGGGGGGGGAGACGGGGTGTTATCTGAATTAGTTT-------------TT--GGGAGTCTCTAGACATCTCATGGGTTGGCCGGGGGCCTGCCGTCTATAGTTCTTATTCCTTTTAAGGG----AGTAAGAAT-----TTCGATTCAGCAA-CTTTAGTTCACAGTCTTTTTTTTTATTAAG-AAAGGTTT-------------------------------------------------------------------------------------------------",
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            "                                                                                                                                                                  .||||||                                .|||||.||      |||..|..||  ||||.||||.||.|    ||.|         ||.|.|||||.|||.||||.||||      | |||||||||||.|.|||||||     ||||||||.|  ||.|      |||.|.||||||||                 ||||||    .||||...||||..|||||..| |||||||||||             ||  ||.||.||.||                     ||..||.||.|.|||..||||.||  |||||    |    ||| |.|||     |||||||||.||| .||||||...|||||||||||||||||..| ||||||||                                                                                                 ",
        )

    def test_local_water2(self):
        """Test parsing a local alignment."""
        path = "Emboss/water2.txt"
        with open(path) as stream:
            alignments = AlignmentIterator(stream)
            self.assertEqual(alignments.program, "water")
            self.assertEqual(alignments.rundate, "Sat Apr 04 2009 22:08:44")
            self.assertEqual(
                alignments.commandline,
                "water -asequence asis:ACACACTCACACACACTTGGTCAGAGATGCTGTGCTTCTTGGAAGCAAGGNCTCAAAGGCAAGGTGCACGCAGAGGGACGTTTGAGTCTGGGATGAAGCATGTNCGTATTATTTATATGATGGAATTTCACGTTTTTATG -bsequence asis:CGTTTGAGTACTGGGATG -gapopen 10 -gapextend 0.5 -filter",
            )
            self.assertEqual(alignments.align_format, "srspair")
            self.assertEqual(alignments.report_file, "stdout")
            alignments = list(alignments)
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertEqual(alignment.annotations["matrix"], "EDNAFULL")
        self.assertAlmostEqual(alignment.annotations["gap_penalty"], 10.0)
        self.assertAlmostEqual(alignment.annotations["extend_penalty"], 0.5)
        self.assertEqual(alignment.annotations["identity"], 17)
        self.assertEqual(alignment.annotations["similarity"], 17)
        self.assertEqual(alignment.annotations["gaps"], 1)
        self.assertAlmostEqual(alignment.annotations["score"], 75.0)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 18))
        self.assertEqual(alignment.sequences[0].id, "asis")
        self.assertEqual(alignment.sequences[1].id, "asis")
        self.assertEqual(
            repr(alignment.sequences[0].seq),
            "Seq({78: 'CGTTTGAGTCTGGGATG'}, length=95)",
        )
        self.assertEqual(alignment.sequences[0].seq[78:95], "CGTTTGAGTCTGGGATG")
        self.assertEqual(alignment.sequences[1].seq[0:18], "CGTTTGAGTACTGGGATG")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[78, 87, 87, 95], [0, 9, 10, 18]])
            )
        )
        self.assertEqual(alignment[0], "CGTTTGAGT-CTGGGATG")
        self.assertEqual(alignment[1], "CGTTTGAGTACTGGGATG")
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"], "||||||||| ||||||||"
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
