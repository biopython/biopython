# Copyright 2021 by Michiel de Hoon.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.Align.msf module."""
import unittest


from Bio.Align.msf import AlignmentIterator

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align.msf."
    ) from None


class TestMSF(unittest.TestCase):
    def test_protein(self):
        path = "msf/W_prot.msf"
        with open(path) as stream:
            alignments = AlignmentIterator(stream)
            alignments = list(alignments)
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertEqual(len(alignment), 11)
        self.assertEqual(alignment.shape, (11, 99))
        self.assertEqual(alignment.sequences[0].id, "W*01:01:01:01")
        self.assertEqual(alignment.sequences[1].id, "W*01:01:01:02")
        self.assertEqual(alignment.sequences[2].id, "W*01:01:01:03")
        self.assertEqual(alignment.sequences[3].id, "W*01:01:01:04")
        self.assertEqual(alignment.sequences[4].id, "W*01:01:01:05")
        self.assertEqual(alignment.sequences[5].id, "W*01:01:01:06")
        self.assertEqual(alignment.sequences[6].id, "W*02:01")
        self.assertEqual(alignment.sequences[7].id, "W*03:01:01:01")
        self.assertEqual(alignment.sequences[8].id, "W*03:01:01:02")
        self.assertEqual(alignment.sequences[9].id, "W*04:01")
        self.assertEqual(alignment.sequences[10].id, "W*05:01")
        self.assertEqual(
            alignment.sequences[0].seq,
            "GLTPFNGYTAATWTRTAVSSVGMNIPYHGASYLVRNQELRSWTAADKAAQMPWRRNRQSCSKPTCREGGRSGSAKSLRMGRRGCSAQNPKDSHDPPPHL",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "GLTPFNGYTAATWTRTAVSSVGMNIPYHGASYLVRNQELRSWTAADKAAQMPWRRNRQSCSKPTCREGGRSGSAKSLRMGRRGCSAQNPKDSHDPPPHL",
        )
        self.assertEqual(
            alignment.sequences[2].seq,
            "GLTPFNGYTAATWTRTAVSSVGMNIPYHGASYLVRNQELRSWTAADKAAQMPWRRNRQSCSKPTCREGGRSGSAKSLRMGRRGCSAQNPKDSHDPPPHL",
        )
        self.assertEqual(
            alignment.sequences[3].seq,
            "GLTPFNGYTAATWTRTAVSSVGMNIPYHGASYLVRNQELRSWTAADKAAQMPWRRNRQSCSKPTCREGGRSGSAKSLRMGRRGCSAQNPKDSHDPPPHL",
        )
        self.assertEqual(
            alignment.sequences[4].seq,
            "GLTPFNGYTAATWTRTAVSSVGMNIPYHGASYLVRNQELRSWTAADKAAQMPWRRNRQSCSKPTCREGGRSGSAKSLRMGRRGCSAQNPKDSHDPPPHL",
        )
        self.assertEqual(
            alignment.sequences[5].seq,
            "GLTPFNGYTAATWTRTAVSSVGMNIPYHGASYLVRNQELRSWTAADKAAQMPWRRNRQSCSKPTCREGGRSGSAKSLRMGRRGCSAQNPKDSHDPPPHL",
        )
        self.assertEqual(
            alignment.sequences[6].seq,
            "GLTPSNGYTAATWTRTAASSVGMNIPYDGASYLVRNQELRSWTAADKAAQMPWRRNMQSCSKPTCREGGRSGSAKSLRMGRRRCTAQNPKRLT",
        )
        self.assertEqual(
            alignment.sequences[7].seq,
            "GLTPSSGYTAATWTRTAVSSVGMNIPYHGASYLVRNQELRSWTAADKAAQMPWRRNRQSCSKPTCREGGRSGSAKSLRMGRRGCSAQNPKRLT",
        )
        self.assertEqual(
            alignment.sequences[8].seq,
            "GLTPSSGYTAATWTRTAVSSVGMNIPYHGASYLVRNQELRSWTAADKAAQMPWRRNRQSCSKPTCREGGRSGSAKSLRMGRRGCSAQNPKRLT",
        )
        self.assertEqual(
            alignment.sequences[9].seq,
            "GLTPSNGYTAATWTRTAASSVGMNIPYDGASYLVRNQELRSWTAADKAAQMPWRRNMQSCSKPTCREGGRSGSAKSLRMGRRGCSAQNPKRLT",
        )
        self.assertEqual(
            alignment.sequences[10].seq,
            "GLTPSSGYTAATWTRTAVSSVGMNIPYHGASYLVRNQELRSWTAADKAAQMPWRRNRQSCSKPTCREGGRSGSAKSLRMGRRGCSAQNPKDSHDPPPHL",
        )

        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [
                        [0, 93, 99],
                        [0, 93, 99],
                        [0, 93, 99],
                        [0, 93, 99],
                        [0, 93, 99],
                        [0, 93, 99],
                        [0, 93, 93],
                        [0, 93, 93],
                        [0, 93, 93],
                        [0, 93, 93],
                        [0, 93, 99],
                    ]
                ),
            )
        )
        self.assertEqual(
            alignment[0],
            "GLTPFNGYTAATWTRTAVSSVGMNIPYHGASYLVRNQELRSWTAADKAAQMPWRRNRQSCSKPTCREGGRSGSAKSLRMGRRGCSAQNPKDSHDPPPHL",
        )
        self.assertEqual(
            alignment[1],
            "GLTPFNGYTAATWTRTAVSSVGMNIPYHGASYLVRNQELRSWTAADKAAQMPWRRNRQSCSKPTCREGGRSGSAKSLRMGRRGCSAQNPKDSHDPPPHL",
        )
        self.assertEqual(
            alignment[2],
            "GLTPFNGYTAATWTRTAVSSVGMNIPYHGASYLVRNQELRSWTAADKAAQMPWRRNRQSCSKPTCREGGRSGSAKSLRMGRRGCSAQNPKDSHDPPPHL",
        )
        self.assertEqual(
            alignment[3],
            "GLTPFNGYTAATWTRTAVSSVGMNIPYHGASYLVRNQELRSWTAADKAAQMPWRRNRQSCSKPTCREGGRSGSAKSLRMGRRGCSAQNPKDSHDPPPHL",
        )
        self.assertEqual(
            alignment[4],
            "GLTPFNGYTAATWTRTAVSSVGMNIPYHGASYLVRNQELRSWTAADKAAQMPWRRNRQSCSKPTCREGGRSGSAKSLRMGRRGCSAQNPKDSHDPPPHL",
        )
        self.assertEqual(
            alignment[5],
            "GLTPFNGYTAATWTRTAVSSVGMNIPYHGASYLVRNQELRSWTAADKAAQMPWRRNRQSCSKPTCREGGRSGSAKSLRMGRRGCSAQNPKDSHDPPPHL",
        )
        self.assertEqual(
            alignment[6],
            "GLTPSNGYTAATWTRTAASSVGMNIPYDGASYLVRNQELRSWTAADKAAQMPWRRNMQSCSKPTCREGGRSGSAKSLRMGRRRCTAQNPKRLT------",
        )
        self.assertEqual(
            alignment[7],
            "GLTPSSGYTAATWTRTAVSSVGMNIPYHGASYLVRNQELRSWTAADKAAQMPWRRNRQSCSKPTCREGGRSGSAKSLRMGRRGCSAQNPKRLT------",
        )
        self.assertEqual(
            alignment[8],
            "GLTPSSGYTAATWTRTAVSSVGMNIPYHGASYLVRNQELRSWTAADKAAQMPWRRNRQSCSKPTCREGGRSGSAKSLRMGRRGCSAQNPKRLT------",
        )
        self.assertEqual(
            alignment[9],
            "GLTPSNGYTAATWTRTAASSVGMNIPYDGASYLVRNQELRSWTAADKAAQMPWRRNMQSCSKPTCREGGRSGSAKSLRMGRRGCSAQNPKRLT------",
        )
        self.assertEqual(
            alignment[10],
            "GLTPSSGYTAATWTRTAVSSVGMNIPYHGASYLVRNQELRSWTAADKAAQMPWRRNRQSCSKPTCREGGRSGSAKSLRMGRRGCSAQNPKDSHDPPPHL",
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
