# Copyright 2021 by Michiel de Hoon.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.Align.msf module."""
import unittest
import warnings


from Bio import BiopythonParserWarning
from Bio import BiopythonExperimentalWarning

with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonExperimentalWarning)
    from Bio.Align.msf import AlignmentIterator


try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align.msf."
    ) from None


class TestMSF(unittest.TestCase):
    def test_protein1(self):
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

    def test_protein2(self):
        path = "msf/DOA_prot.msf"

        with warnings.catch_warnings(record=True) as w:
            with open(path) as stream:
                alignments = AlignmentIterator(stream)
                alignments = list(alignments)
        self.assertEqual(len(w), 1)
        self.assertIsInstance(w[0].message, BiopythonParserWarning)
        self.assertEqual(
            str(w[0].message), "GCG MSF headers said alignment length 62, but found 250"
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertEqual(len(alignment), 12)
        self.assertEqual(alignment.shape, (12, 250))
        self.assertEqual(alignment.sequences[0].id, "DOA*01:01:01")
        self.assertEqual(alignment.sequences[1].id, "DOA*01:01:02:01")
        self.assertEqual(alignment.sequences[2].id, "DOA*01:01:02:02")
        self.assertEqual(alignment.sequences[3].id, "DOA*01:01:02:03")
        self.assertEqual(alignment.sequences[4].id, "DOA*01:01:03")
        self.assertEqual(alignment.sequences[5].id, "DOA*01:01:04:01")
        self.assertEqual(alignment.sequences[6].id, "DOA*01:01:04:02")
        self.assertEqual(alignment.sequences[7].id, "DOA*01:01:05")
        self.assertEqual(alignment.sequences[8].id, "DOA*01:01:06")
        self.assertEqual(alignment.sequences[9].id, "DOA*01:02")
        self.assertEqual(alignment.sequences[10].id, "DOA*01:03")
        self.assertEqual(alignment.sequences[11].id, "DOA*01:04N")
        self.assertEqual(
            alignment.sequences[0].seq,
            "MALRAGLVLGFHTLMTLLSPQEAGATKADHMGSYGPAFYQSYGASGQFTHEFDEEQLFSVDLKKSEAVWRLPEFGDFARFDPQGGLAGIAAIKAHLDILVERSNRSRAINVPPRVTVLPKSRVELGQPNILICIVDNIFPPVINITWLRNGQTVTEGVAQTSFYSQPDHLFRKFHYLPFVPSAEDVYDCQVEHWGLDAPLLRHWELQVPIPPPDAMETLVCALGLAIGLVGFLVGTVLIIMGTYVSSVPR",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "MALRAGLVLGFHTLMTLLSPQEAGATKADHMGSYGPAFYQSYGASGQFTHEFDEEQLFSVDLKKSEAVWRLPEFGDFARFDPQGGLAGIAAIKAHLDILVERSNRSRAINVPPRVTVLPKSRVELGQPNILICIVDNIFPPVINITWLRNGQTVTEGVAQTSFYSQPDHLFRKFHYLPFVPSAEDVYDCQVEHWGLDAPLLRHWELQVPIPPPDAMETLVCALGLAIGLVGFLVGTVLIIMGTYVSSVPR",
        )
        self.assertEqual(
            alignment.sequences[2].seq,
            "MALRAGLVLGFHTLMTLLSPQEAGATKADHMGSYGPAFYQSYGASGQFTHEFDEEQLFSVDLKKSEAVWRLPEFGDFARFDPQGGLAGIAAIKAHLDILVERSNRSRAINVPPRVTVLPKSRVELGQPNILICIVDNIFPPVINITWLRNGQTVTEGVAQTSFYSQPDHLFRKFHYLPFVPSAEDVYDCQVEHWGLDAPLLRHWELQVPIPPPDAMETLVCALGLAIGLVGFLVGTVLIIMGTYVSSVPR",
        )
        self.assertEqual(
            alignment.sequences[3].seq,
            "MALRAGLVLGFHTLMTLLSPQEAGATKADHMGSYGPAFYQSYGASGQFTHEFDEEQLFSVDLKKSEAVWRLPEFGDFARFDPQGGLAGIAAIKAHLDILVERSNRSRAINVPPRVTVLPKSRVELGQPNILICIVDNIFPPVINITWLRNGQTVTEGVAQTSFYSQPDHLFRKFHYLPFVPSAEDVYDCQVEHWGLDAPLLRHWELQVPIPPPDAMETLVCALGLAIGLVGFLVGTVLIIMGTYVSSVPR",
        )
        self.assertEqual(
            alignment.sequences[4].seq,
            "DHMGSYGPAFYQSYGASGQFTHEFDEEQLFSVDLKKSEAVWRLPEFGDFARFDPQGGLAGIAAIKAHLDILVERSNRSRAINVPPRVTVLPKSRVELGQPNILICIVDNIFPPVINITWLRNGQTVTEGVAQTSFYSQPDHLFRKFHYLPFVPSAEDVYDCQVEHWGLDAPLLRHWELQVPIPPPDAMETLVCALGLAIGLVGFLVGTVLIIMGTYVSSVPR",
        )
        self.assertEqual(
            alignment.sequences[5].seq,
            "MALRAGLVLGFHTLMTLLSPQEAGATKADHMGSYGPAFYQSYGASGQFTHEFDEEQLFSVDLKKSEAVWRLPEFGDFARFDPQGGLAGIAAIKAHLDILVERSNRSRAINVPPRVTVLPKSRVELGQPNILICIVDNIFPPVINITWLRNGQTVTEGVAQTSFYSQPDHLFRKFHYLPFVPSAEDVYDCQVEHWGLDAPLLRHWELQVPIPPPDAMETLVCALGLAIGLVGFLVGTVLIIMGTYVSSVPR",
        )
        self.assertEqual(
            alignment.sequences[6].seq,
            "DHMGSYGPAFYQSYGASGQFTHEFDEEQLFSVDLKKSEAVWRLPEFGDFARFDPQGGLAGIAAIKAHLDILVERSNRSRAINVPPRVTVLPKSRVELGQPNILICIVDNIFPPVINITWLRNGQTVTEGVAQTSFYSQPDHLFRKFHYLPFVPSAEDVYDCQVEHWGLDAPLLRHWELQVPIPPPDAMETLVCALGLAIGLVGFLVGTVLIIMGTYVSSVPR",
        )
        self.assertEqual(
            alignment.sequences[7].seq,
            "MALRAGLVLGFHTLMTLLSPQEAGATKADHMGSYGPAFYQSYGASGQFTHEFDEEQLFSVDLKKSEAVWRLPEFGDFARFDPQGGLAGIAAIKAHLDILVERSNRSRAINVPPRVTVLPKSRVELGQPNILICIVDNIFPPVINITWLRNGQTVTEGVAQTSFYSQPDHLFRKFHYLPFVPSAEDVYDCQVEHWGLDAPLLRHWELQVPIPPPDAMETLVCALGLAIGLVGFLVGTVLIIMGTYVSSVPR",
        )
        self.assertEqual(
            alignment.sequences[8].seq,
            "MALRAGLVLGFHTLMTLLSPQEAGATKADHMGSYGPAFYQSYGASGQFTHEFDEEQLFSVDLKKSEAVWRLPEFGDFARFDPQGGLAGIAAIKAHLDILVERSNRSRAINVPPRVTVLPKSRVELGQPNILICIVDNIFPPVINITWLRNGQTVTEGVAQTSFYSQPDHLFRKFHYLPFVPSAEDVYDCQVEHWGLDAPLLRHWELQVPIPPPDAMETLVCALGLAIGLVGFLVGTVLIIMGTYVSSVPR",
        )
        self.assertEqual(
            alignment.sequences[9].seq,
            "MALRAGLVLGFHTLMTLLSPQEAGATKADHMGSYGPAFYQSYGASGQFTHEFDEEQLFSVDLKKSEAVWRLPEFGDFARFDPQGGLAGIAAIKAHLDILVERSNCSRAINVPPRVTVLPKSRVELGQPNILICIVDNIFPPVINITWLRNGQTVTEGVAQTSFYSQPDHLFRKFHYLPFVPSAEDVYDCQVEHWGLDAPLLRHWELQVPIPPPDAMETLVCALGLAIGLVGFLVGTVLIIMGTYVSSVPR",
        )
        self.assertEqual(
            alignment.sequences[10].seq,
            "MALRAGLVLGFHTLMTLLSPQEAGATKADHMGSYGPAFYQSYGASGQFTHEFDEEQLFSVDLKKSEAVWRLPEFGDFARFDPQGGLAGIAAIKAHLDIVVERSNRSRAINVPPRVTVLPKSRVELGQPNILICIVDNIFPPVINITWLRNGQTVTEGVAQTSFYSQPDHLFRKFHYLPFVPSAEDVYDCQVEHWGLDAPLLRHWELQVPIPPPDAMETLVCALGLAIGLVGFLVGTVLIIMGTYVSSVPR",
        )
        self.assertEqual(
            alignment.sequences[11].seq,
            "MALRAGLVLGFHTLMTLLSPQEAGATKADHMGSYGPPSTSLTAPRASSPMNLMRNSCSLWTX",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [
                        [0, 28, 62, 250],
                        [0, 28, 62, 250],
                        [0, 28, 62, 250],
                        [0, 28, 62, 250],
                        [0, 0, 34, 222],
                        [0, 28, 62, 250],
                        [0, 0, 34, 222],
                        [0, 28, 62, 250],
                        [0, 28, 62, 250],
                        [0, 28, 62, 250],
                        [0, 28, 62, 250],
                        [0, 28, 62, 62],
                    ],
                ),
            )
        )
        self.assertEqual(
            alignment[0],
            "MALRAGLVLGFHTLMTLLSPQEAGATKADHMGSYGPAFYQSYGASGQFTHEFDEEQLFSVDLKKSEAVWRLPEFGDFARFDPQGGLAGIAAIKAHLDILVERSNRSRAINVPPRVTVLPKSRVELGQPNILICIVDNIFPPVINITWLRNGQTVTEGVAQTSFYSQPDHLFRKFHYLPFVPSAEDVYDCQVEHWGLDAPLLRHWELQVPIPPPDAMETLVCALGLAIGLVGFLVGTVLIIMGTYVSSVPR",
        )
        self.assertEqual(
            alignment[1],
            "MALRAGLVLGFHTLMTLLSPQEAGATKADHMGSYGPAFYQSYGASGQFTHEFDEEQLFSVDLKKSEAVWRLPEFGDFARFDPQGGLAGIAAIKAHLDILVERSNRSRAINVPPRVTVLPKSRVELGQPNILICIVDNIFPPVINITWLRNGQTVTEGVAQTSFYSQPDHLFRKFHYLPFVPSAEDVYDCQVEHWGLDAPLLRHWELQVPIPPPDAMETLVCALGLAIGLVGFLVGTVLIIMGTYVSSVPR",
        )
        self.assertEqual(
            alignment[2],
            "MALRAGLVLGFHTLMTLLSPQEAGATKADHMGSYGPAFYQSYGASGQFTHEFDEEQLFSVDLKKSEAVWRLPEFGDFARFDPQGGLAGIAAIKAHLDILVERSNRSRAINVPPRVTVLPKSRVELGQPNILICIVDNIFPPVINITWLRNGQTVTEGVAQTSFYSQPDHLFRKFHYLPFVPSAEDVYDCQVEHWGLDAPLLRHWELQVPIPPPDAMETLVCALGLAIGLVGFLVGTVLIIMGTYVSSVPR",
        )
        self.assertEqual(
            alignment[3],
            "MALRAGLVLGFHTLMTLLSPQEAGATKADHMGSYGPAFYQSYGASGQFTHEFDEEQLFSVDLKKSEAVWRLPEFGDFARFDPQGGLAGIAAIKAHLDILVERSNRSRAINVPPRVTVLPKSRVELGQPNILICIVDNIFPPVINITWLRNGQTVTEGVAQTSFYSQPDHLFRKFHYLPFVPSAEDVYDCQVEHWGLDAPLLRHWELQVPIPPPDAMETLVCALGLAIGLVGFLVGTVLIIMGTYVSSVPR",
        )
        self.assertEqual(
            alignment[4],
            "----------------------------DHMGSYGPAFYQSYGASGQFTHEFDEEQLFSVDLKKSEAVWRLPEFGDFARFDPQGGLAGIAAIKAHLDILVERSNRSRAINVPPRVTVLPKSRVELGQPNILICIVDNIFPPVINITWLRNGQTVTEGVAQTSFYSQPDHLFRKFHYLPFVPSAEDVYDCQVEHWGLDAPLLRHWELQVPIPPPDAMETLVCALGLAIGLVGFLVGTVLIIMGTYVSSVPR",
        )
        self.assertEqual(
            alignment[5],
            "MALRAGLVLGFHTLMTLLSPQEAGATKADHMGSYGPAFYQSYGASGQFTHEFDEEQLFSVDLKKSEAVWRLPEFGDFARFDPQGGLAGIAAIKAHLDILVERSNRSRAINVPPRVTVLPKSRVELGQPNILICIVDNIFPPVINITWLRNGQTVTEGVAQTSFYSQPDHLFRKFHYLPFVPSAEDVYDCQVEHWGLDAPLLRHWELQVPIPPPDAMETLVCALGLAIGLVGFLVGTVLIIMGTYVSSVPR",
        )
        self.assertEqual(
            alignment[6],
            "----------------------------DHMGSYGPAFYQSYGASGQFTHEFDEEQLFSVDLKKSEAVWRLPEFGDFARFDPQGGLAGIAAIKAHLDILVERSNRSRAINVPPRVTVLPKSRVELGQPNILICIVDNIFPPVINITWLRNGQTVTEGVAQTSFYSQPDHLFRKFHYLPFVPSAEDVYDCQVEHWGLDAPLLRHWELQVPIPPPDAMETLVCALGLAIGLVGFLVGTVLIIMGTYVSSVPR",
        )
        self.assertEqual(
            alignment[7],
            "MALRAGLVLGFHTLMTLLSPQEAGATKADHMGSYGPAFYQSYGASGQFTHEFDEEQLFSVDLKKSEAVWRLPEFGDFARFDPQGGLAGIAAIKAHLDILVERSNRSRAINVPPRVTVLPKSRVELGQPNILICIVDNIFPPVINITWLRNGQTVTEGVAQTSFYSQPDHLFRKFHYLPFVPSAEDVYDCQVEHWGLDAPLLRHWELQVPIPPPDAMETLVCALGLAIGLVGFLVGTVLIIMGTYVSSVPR",
        )
        self.assertEqual(
            alignment[8],
            "MALRAGLVLGFHTLMTLLSPQEAGATKADHMGSYGPAFYQSYGASGQFTHEFDEEQLFSVDLKKSEAVWRLPEFGDFARFDPQGGLAGIAAIKAHLDILVERSNRSRAINVPPRVTVLPKSRVELGQPNILICIVDNIFPPVINITWLRNGQTVTEGVAQTSFYSQPDHLFRKFHYLPFVPSAEDVYDCQVEHWGLDAPLLRHWELQVPIPPPDAMETLVCALGLAIGLVGFLVGTVLIIMGTYVSSVPR",
        )
        self.assertEqual(
            alignment[9],
            "MALRAGLVLGFHTLMTLLSPQEAGATKADHMGSYGPAFYQSYGASGQFTHEFDEEQLFSVDLKKSEAVWRLPEFGDFARFDPQGGLAGIAAIKAHLDILVERSNCSRAINVPPRVTVLPKSRVELGQPNILICIVDNIFPPVINITWLRNGQTVTEGVAQTSFYSQPDHLFRKFHYLPFVPSAEDVYDCQVEHWGLDAPLLRHWELQVPIPPPDAMETLVCALGLAIGLVGFLVGTVLIIMGTYVSSVPR",
        )
        self.assertEqual(
            alignment[10],
            "MALRAGLVLGFHTLMTLLSPQEAGATKADHMGSYGPAFYQSYGASGQFTHEFDEEQLFSVDLKKSEAVWRLPEFGDFARFDPQGGLAGIAAIKAHLDIVVERSNRSRAINVPPRVTVLPKSRVELGQPNILICIVDNIFPPVINITWLRNGQTVTEGVAQTSFYSQPDHLFRKFHYLPFVPSAEDVYDCQVEHWGLDAPLLRHWELQVPIPPPDAMETLVCALGLAIGLVGFLVGTVLIIMGTYVSSVPR",
        )
        self.assertEqual(
            alignment[11],
            "MALRAGLVLGFHTLMTLLSPQEAGATKADHMGSYGPPSTSLTAPRASSPMNLMRNSCSLWTX--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------",
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
