# Copyright 2008 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Align.stockholm module."""
import unittest
import warnings
from io import StringIO

from Bio import BiopythonExperimentalWarning

with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonExperimentalWarning)
    from Bio.Align import stockholm


try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align.stockholm."
    ) from None


class TestStockholm_reading(unittest.TestCase):
    def test_reading_example(self):
        """Test parsing Pfam record HAT as the docstring example."""
        path = "Stockholm/example.sth"
        alignments = stockholm.AlignmentIterator(path)
        alignment = next(alignments)
        self.assertEqual(alignment.annotations["identifier"], "HAT")
        self.assertEqual(alignment.annotations["accession"], "PF02184.18")
        self.assertEqual(alignment.annotations["definition"], "HAT (Half-A-TPR) repeat")
        self.assertEqual(alignment.annotations["author"], ["SMART;"])
        self.assertEqual(
            alignment.annotations["source of seed"],
            "Alignment kindly provided by SMART",
        )
        self.assertEqual(alignment.annotations["gathering method"], "21.00 21.00;")
        self.assertEqual(alignment.annotations["trusted cutoff"], "21.00 21.00;")
        self.assertEqual(alignment.annotations["noise cutoff"], "20.90 20.90;")
        self.assertEqual(
            alignment.annotations["build method"], "hmmbuild HMM.ann SEED.ann"
        )
        self.assertEqual(
            alignment.annotations["search method"],
            "hmmsearch -Z 57096847 -E 1000 --cpu 4 HMM pfamseq",
        )
        self.assertEqual(alignment.annotations["type"], "Repeat")
        self.assertEqual(alignment.annotations["clan"], "CL0020")
        self.assertEqual(len(alignment.annotations["references"]), 1)
        self.assertEqual(alignment.annotations["references"][0]["number"], 1)
        self.assertEqual(alignment.annotations["references"][0]["medline"], "9478129")
        self.assertEqual(
            alignment.annotations["references"][0]["title"],
            "The HAT helix, a repetitive motif implicated in RNA processing.",
        )
        self.assertEqual(
            alignment.annotations["references"][0]["author"], "Preker PJ, Keller W;"
        )
        self.assertEqual(
            alignment.annotations["references"][0]["location"],
            "Trends Biochem Sci 1998;23:15-16.",
        )
        self.assertEqual(len(alignment.annotations["database references"]), 3)
        self.assertEqual(
            alignment.annotations["database references"][0],
            {"reference": "INTERPRO; IPR003107;"},
        )
        self.assertEqual(
            alignment.annotations["database references"][1],
            {"reference": "SMART; HAT;"},
        )
        self.assertEqual(
            alignment.annotations["database references"][2],
            {"reference": "SO; 0001068; polypeptide_repeat;"},
        )
        self.assertEqual(
            alignment.annotations["comment"],
            "The HAT (Half A TPR) repeat is found in several RNA processing proteins [1].",
        )
        self.assertEqual(len(alignment.sequences), 3)
        self.assertEqual(alignment.sequences[0].annotations["accession"], "P17886.2")
        self.assertEqual(alignment.sequences[1].annotations["accession"], "P87312.1")
        self.assertEqual(alignment.sequences[1].dbxrefs, ["PDB; 3JB9 R; 185-216;"])
        self.assertEqual(alignment.sequences[2].annotations["accession"], "O16376.2")
        self.assertEqual(alignment.sequences[0].id, "CRN_DROME/191-222")
        self.assertEqual(alignment.sequences[1].id, "CLF1_SCHPO/185-216")
        self.assertEqual(alignment.sequences[2].id, "O16376_CAEEL/201-233")
        self.assertEqual(alignment[0], "KEIDRAREIYERFVYVH-PDVKNWIKFARFEES")
        self.assertEqual(alignment[1], "HENERARGIYERFVVVH-PEVTNWLRWARFEEE")
        self.assertEqual(alignment[2], "KEIDRARSVYQRFLHVHGINVQNWIKYAKFEER")
        self.assertEqual(alignment.sequences[0].seq, "KEIDRAREIYERFVYVHPDVKNWIKFARFEES")
        self.assertEqual(alignment.sequences[1].seq, "HENERARGIYERFVVVHPEVTNWLRWARFEEE")
        self.assertEqual(
            alignment.sequences[1].letter_annotations["secondary structure"],
            "--HHHHHHHHHHHHHHS--HHHHHHHHHHHHH",
        )
        self.assertEqual(
            alignment.sequences[2].seq, "KEIDRARSVYQRFLHVHGINVQNWIKYAKFEER"
        )
        self.assertEqual(
            alignment.column_annotations["consensus secondary structure"],
            "--HHHHHHHHHHHHHHS.--HHHHHHHHHHHHH",
        )
        self.assertEqual(
            alignment.column_annotations["consensus sequence"],
            "KEIDRARuIYERFVaVH.P-VpNWIKaARFEEc",
        )

    def check_alignment_globins45(self, alignment):
        """Check the alignment obtained by parsing hmmalign output."""
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
# flake8: noqa
[[0, 0, 0, 1, 5, 17, 19, 46, 47, 50, 51, 55, 78, 78, 79, 80, 146, 147, 148, 153],
 [0, 1, 1, 1, 5, 17, 19, 46, 47, 50, 51, 55, 78, 78, 79, 80, 146, 147, 148, 153],
 [0, 1, 1, 1, 5, 17, 19, 46, 47, 50, 51, 55, 78, 78, 79, 80, 146, 147, 148, 153],
 [0, 1, 1, 1, 5, 17, 19, 46, 47, 50, 51, 55, 78, 78, 79, 80, 146, 147, 148, 153],
 [0, 1, 1, 1, 5, 17, 19, 46, 47, 50, 51, 55, 78, 78, 79, 80, 146, 147, 148, 153],
 [0, 1, 1, 1, 5, 17, 19, 46, 47, 50, 51, 55, 78, 78, 79, 80, 146, 147, 148, 153],
 [0, 1, 1, 1, 1, 13, 15, 42, 43, 46, 46, 50, 73, 73, 74, 75, 141, 142, 143, 148],
 [0, 0, 0, 1, 5, 17, 19, 46, 46, 49, 49, 49, 72, 73, 74, 74, 140, 141, 141, 141],
 [0, 0, 0, 1, 5, 17, 19, 46, 46, 49, 49, 49, 72, 73, 74, 74, 140, 141, 141, 141],
 [0, 0, 0, 1, 5, 17, 19, 46, 46, 49, 49, 49, 72, 73, 74, 74, 140, 141, 141, 141],
 [0, 0, 0, 1, 5, 17, 19, 46, 46, 49, 49, 49, 72, 73, 74, 74, 140, 141, 141, 141],
 [0, 0, 0, 1, 5, 17, 19, 46, 46, 49, 49, 49, 72, 73, 74, 74, 140, 141, 141, 141],
 [0, 0, 0, 1, 5, 17, 19, 46, 46, 49, 49, 49, 72, 73, 74, 74, 140, 141, 141, 141],
 [0, 0, 0, 1, 5, 17, 19, 46, 46, 49, 49, 49, 72, 73, 74, 74, 140, 141, 141, 141],
 [0, 0, 0, 1, 5, 17, 19, 46, 46, 49, 49, 49, 72, 73, 74, 74, 140, 141, 141, 141],
 [0, 0, 0, 1, 5, 17, 19, 46, 46, 49, 49, 49, 72, 73, 74, 74, 140, 141, 141, 141],
 [0, 0, 0, 1, 5, 17, 19, 46, 46, 49, 49, 49, 72, 73, 74, 74, 140, 141, 141, 141],
 [0, 0, 0, 1, 5, 17, 19, 46, 46, 49, 49, 49, 72, 73, 74, 74, 140, 141, 141, 141],
 [0, 0, 0, 1, 5, 17, 19, 46, 46, 49, 49, 49, 72, 73, 74, 74, 140, 141, 141, 141],
 [0, 0, 0, 1, 5, 17, 19, 46, 46, 49, 49, 49, 72, 73, 74, 74, 140, 141, 141, 141],
 [0, 0, 0, 1, 5, 17, 19, 46, 46, 49, 49, 49, 72, 73, 74, 74, 140, 141, 141, 141],
 [0, 0, 0, 1, 5, 17, 19, 46, 46, 49, 49, 49, 72, 73, 74, 74, 140, 141, 141, 141],
 [0, 1, 1, 1, 5, 17, 19, 46, 46, 49, 49, 49, 72, 73, 74, 74, 140, 141, 141, 141],
 [0, 1, 1, 1, 5, 17, 19, 46, 46, 49, 49, 49, 72, 73, 74, 74, 140, 141, 141, 141],
 [0, 1, 1, 1, 5, 17, 19, 46, 46, 49, 49, 49, 72, 73, 74, 74, 140, 141, 141, 141],
 [0, 1, 1, 1, 5, 17, 19, 46, 47, 50, 50, 50, 73, 74, 75, 75, 141, 142, 142, 142],
 [0, 0, 1, 2, 6, 18, 18, 45, 46, 49, 50, 54, 77, 78, 79, 79, 145, 146, 146, 146],
 [0, 0, 1, 2, 6, 18, 18, 45, 46, 49, 50, 54, 77, 78, 79, 79, 145, 146, 146, 146],
 [0, 0, 1, 2, 6, 18, 18, 45, 46, 49, 50, 54, 77, 78, 79, 79, 145, 146, 146, 146],
 [0, 0, 1, 2, 6, 18, 18, 45, 46, 49, 50, 54, 77, 78, 79, 79, 145, 146, 146, 146],
 [0, 0, 1, 2, 6, 18, 18, 45, 46, 49, 50, 54, 77, 78, 79, 79, 145, 146, 146, 146],
 [0, 1, 1, 2, 6, 18, 18, 45, 46, 49, 50, 54, 77, 78, 79, 79, 145, 146, 146, 146],
 [0, 0, 1, 2, 6, 18, 18, 45, 46, 49, 50, 54, 77, 78, 79, 79, 145, 146, 146, 146],
 [0, 0, 1, 2, 6, 18, 18, 45, 46, 49, 50, 54, 77, 78, 79, 79, 145, 146, 146, 146],
 [0, 0, 1, 2, 6, 18, 18, 45, 46, 49, 50, 54, 77, 78, 79, 79, 145, 146, 146, 146],
 [0, 0, 1, 2, 6, 18, 18, 45, 46, 49, 50, 54, 77, 78, 79, 79, 145, 146, 146, 146],
 [0, 0, 1, 2, 6, 18, 18, 45, 46, 49, 50, 54, 77, 78, 79, 79, 145, 146, 146, 146],
 [0, 0, 1, 2, 6, 18, 18, 45, 46, 49, 50, 54, 77, 78, 79, 79, 145, 146, 146, 146],
 [0, 0, 1, 2, 6, 18, 18, 45, 46, 49, 50, 54, 77, 78, 79, 79, 145, 146, 146, 146],
 [0, 1, 1, 2, 6, 18, 18, 45, 46, 49, 50, 54, 77, 78, 79, 79, 145, 146, 146, 146],
 [0, 1, 1, 2, 6, 18, 18, 45, 46, 49, 50, 54, 77, 78, 79, 79, 145, 146, 146, 146],
 [0, 1, 1, 2, 6, 18, 18, 45, 46, 49, 50, 54, 77, 78, 79, 79, 145, 146, 146, 146],
 [0, 1, 1, 2, 6, 18, 18, 45, 46, 49, 50, 54, 77, 78, 79, 79, 145, 145, 146, 146],
 [0, 1, 1, 2, 6, 18, 18, 45, 46, 49, 50, 54, 77, 78, 79, 79, 145, 146, 146, 146],
 [0, 0, 1, 2, 6, 18, 18, 45, 46, 49, 50, 54, 77, 78, 79, 79, 145, 145, 145, 145],
]
                    # fmt: on
                ),
            )
        )
        self.assertEqual(alignment.sequences[0].id, "MYG_ESCGI")
        self.assertEqual(
            alignment.sequences[0].seq,
            "VLSDAEWQLVLNIWAKVEADVAGHGQDILIRLFKGHPETLEKFDKFKHLKTEAEMKASEDLKKHGNTVLTALGGILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISDAIIHVLHSRHPGDFGADAQAAMNKALELFRKDIAAKYKelgfqg",
        )
        self.assertEqual(
            alignment.sequences[0].letter_annotations["posterior probability"],
            "69****************************************************************************99******************************************************************7******",
        )
        self.assertEqual(alignment.sequences[1].id, "MYG_HORSE")
        self.assertEqual(
            alignment.sequences[1].seq,
            "gLSDGEWQQVLNVWGKVEADIAGHGQEVLIRLFTGHPETLEKFDKFKHLKTEAEMKASEDLKKHGTVVLTALGGILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISDAIIHVLHSKHPGNFGADAQGAMTKALELFRNDIAAKYKelgfqg",
        )
        self.assertEqual(
            alignment.sequences[1].letter_annotations["posterior probability"],
            "889***************************************************************************99******************************************************************7******",
        )
        self.assertEqual(alignment.sequences[2].id, "MYG_PROGU")
        self.assertEqual(
            alignment.sequences[2].seq,
            "gLSDGEWQLVLNVWGKVEGDLSGHGQEVLIRLFKGHPETLEKFDKFKHLKAEDEMRASEELKKHGTTVLTALGGILKKKGQHAAELAPLAQSHATKHKIPVKYLEFISEAIIQVLQSKHPGDFGADAQGAMSKALELFRNDIAAKYKelgfqg",
        )
        self.assertEqual(
            alignment.sequences[2].letter_annotations["posterior probability"],
            "889***************************************************************************99******************************************************************7******",
        )
        self.assertEqual(alignment.sequences[3].id, "MYG_SAISC")
        self.assertEqual(
            alignment.sequences[3].seq,
            "gLSDGEWQLVLNIWGKVEADIPSHGQEVLISLFKGHPETLEKFDKFKHLKSEDEMKASEELKKHGTTVLTALGGILKKKGQHEAELKPLAQSHATKHKIPVKYLELISDAIVHVLQKKHPGDFGADAQGAMKKALELFRNDMAAKYKelgfqg",
        )
        self.assertEqual(
            alignment.sequences[3].letter_annotations["posterior probability"],
            "889***************************************************************************99******************************************************************7******",
        )
        self.assertEqual(alignment.sequences[4].id, "MYG_LYCPI")
        self.assertEqual(
            alignment.sequences[4].seq,
            "gLSDGEWQIVLNIWGKVETDLAGHGQEVLIRLFKNHPETLDKFDKFKHLKTEDEMKGSEDLKKHGNTVLTALGGILKKKGHHEAELKPLAQSHATKHKIPVKYLEFISDAIIQVLQNKHSGDFHADTEAAMKKALELFRNDIAAKYKelgfqg",
        )
        self.assertEqual(
            alignment.sequences[4].letter_annotations["posterior probability"],
            "889***************************************************************************99******************************************************************7******",
        )
        self.assertEqual(alignment.sequences[5].id, "MYG_MOUSE")
        self.assertEqual(
            alignment.sequences[5].seq,
            "gLSDGEWQLVLNVWGKVEADLAGHGQEVLIGLFKTHPETLDKFDKFKNLKSEEDMKGSEDLKKHGCTVLTALGTILKKKGQHAAEIQPLAQSHATKHKIPVKYLEFISEIIIEVLKKRHSGDFGADAQGAMSKALELFRNDIAAKYKelgfqg",
        )
        self.assertEqual(
            alignment.sequences[5].letter_annotations["posterior probability"],
            "889***************************************************************************99******************************************************************7******",
        )
        self.assertEqual(alignment.sequences[6].id, "MYG_MUSAN")
        self.assertEqual(
            alignment.sequences[6].seq,
            "vDWEKVNSVWSAVESDLTAIGQNILLRLFEQYPESQNHFPKFKNKSLGELKDTADIKAQADTVLSALGNIVKKKGSHSQPVKALAATHITTHKIPPHYFTKITTIAVDVLSEMYPSEMNAQVQAAFSGAFKIICSDIEKEYKaanfqg",
        )
        self.assertEqual(
            alignment.sequences[6].letter_annotations["posterior probability"],
            "789***************************************987789*************************99****************************************************************997******",
        )
        self.assertEqual(alignment.sequences[7].id, "HBA_AILME")
        self.assertEqual(
            alignment.sequences[7].seq,
            "VLSPADKTNVKATWDKIGGHAGEYGGEALERTFASFPTTKTYFPHFDLSPGSAQVKAHGKKVADALTTAVGHLDDLPGALSALSDLHAHKLRVDPVNFKLLSHCLLVTLASHHPAEFTPAVHASLDKFFSAVSTVLTSKYR",
        )
        self.assertEqual(
            alignment.sequences[7].letter_annotations["posterior probability"],
            "69********************************************9**9***********************9******************************************************************7",
        )
        self.assertEqual(alignment.sequences[8].id, "HBA_PROLO")
        self.assertEqual(
            alignment.sequences[8].seq,
            "VLSPADKANIKATWDKIGGHAGEYGGEALERTFASFPTTKTYFPHFDLSPGSAQVKAHGKKVADALTLAVGHLDDLPGALSALSDLHAYKLRVDPVNFKLLSHCLLVTLACHHPAEFTPAVHASLDKFFTSVSTVLTSKYR",
        )
        self.assertEqual(
            alignment.sequences[8].letter_annotations["posterior probability"],
            "69********************************************9**9***********************9******************************************************************7",
        )
        self.assertEqual(alignment.sequences[9].id, "HBA_PAGLA")
        self.assertEqual(
            alignment.sequences[9].seq,
            "VLSSADKNNIKATWDKIGSHAGEYGAEALERTFISFPTTKTYFPHFDLSHGSAQVKAHGKKVADALTLAVGHLEDLPNALSALSDLHAYKLRVDPVNFKLLSHCLLVTLACHHPAEFTPAVHSALDKFFSAVSTVLTSKYR",
        )
        self.assertEqual(
            alignment.sequences[9].letter_annotations["posterior probability"],
            "69**********************************************************************989*****************************************************************7",
        )
        self.assertEqual(alignment.sequences[10].id, "HBA_MACFA")
        self.assertEqual(
            alignment.sequences[10].seq,
            "VLSPADKTNVKAAWGKVGGHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTLAVGHVDDMPQALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR",
        )
        self.assertEqual(
            alignment.sequences[10].letter_annotations["posterior probability"],
            "69***********************************************************************9******************************************************************7",
        )
        self.assertEqual(alignment.sequences[11].id, "HBA_MACSI")
        self.assertEqual(
            alignment.sequences[11].seq,
            "VLSPADKTNVKDAWGKVGGHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTLAVGHVDDMPQALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR",
        )
        self.assertEqual(
            alignment.sequences[11].letter_annotations["posterior probability"],
            "69***********************************************************************9******************************************************************7",
        )
        self.assertEqual(alignment.sequences[12].id, "HBA_PONPY")
        self.assertEqual(
            alignment.sequences[12].seq,
            "VLSPADKTNVKTAWGKVGAHAGDYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKDHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR",
        )
        self.assertEqual(
            alignment.sequences[12].letter_annotations["posterior probability"],
            "69***********************************************************************9******************************************************************7",
        )
        self.assertEqual(alignment.sequences[13].id, "HBA2_GALCR")
        self.assertEqual(
            alignment.sequences[13].seq,
            "VLSPTDKSNVKAAWEKVGAHAGDYGAEALERMFLSFPTTKTYFPHFDLSHGSTQVKGHGKKVADALTNAVLHVDDMPSALSALSDLHAHKLRVDPVNFKLLRHCLLVTLACHHPAEFTPAVHASLDKFMASVSTVLTSKYR",
        )
        self.assertEqual(
            alignment.sequences[13].letter_annotations["posterior probability"],
            "69***********************************************************************9******************************************************************7",
        )
        self.assertEqual(alignment.sequences[14].id, "HBA_MESAU")
        self.assertEqual(
            alignment.sequences[14].seq,
            "VLSAKDKTNISEAWGKIGGHAGEYGAEALERMFFVYPTTKTYFPHFDVSHGSAQVKGHGKKVADALTNAVGHLDDLPGALSALSDLHAHKLRVDPVNFKLLSHCLLVTLANHHPADFTPAVHASLDKFFASVSTVLTSKYR",
        )
        self.assertEqual(
            alignment.sequences[14].letter_annotations["posterior probability"],
            "69********************************************888************************9******************************************************************7",
        )
        self.assertEqual(alignment.sequences[15].id, "HBA2_BOSMU")
        self.assertEqual(
            alignment.sequences[15].seq,
            "VLSAADKGNVKAAWGKVGGHAAEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGAKVAAALTKAVGHLDDLPGALSELSDLHAHKLRVDPVNFKLLSHSLLVTLASHLPSDFTPAVHASLDKFLANVSTVLTSKYR",
        )
        self.assertEqual(
            alignment.sequences[15].letter_annotations["posterior probability"],
            "69***********************************************************************9******************************************************************7",
        )
        self.assertEqual(alignment.sequences[16].id, "HBA_ERIEU")
        self.assertEqual(
            alignment.sequences[16].seq,
            "VLSATDKANVKTFWGKLGGHGGEYGGEALDRMFQAHPTTKTYFPHFDLNPGSAQVKGHGKKVADALTTAVNNLDDVPGALSALSDLHAHKLRVDPVNFKLLSHCLLVTLALHHPADFTPAVHASLDKFLATVATVLTSKYR",
        )
        self.assertEqual(
            alignment.sequences[16].letter_annotations["posterior probability"],
            "69********************************************9999***********************99*****************************************************************7",
        )
        self.assertEqual(alignment.sequences[17].id, "HBA_FRAPO")
        self.assertEqual(
            alignment.sequences[17].seq,
            "VLSAADKNNVKGIFGKISSHAEDYGAEALERMFITYPSTKTYFPHFDLSHGSAQVKGHGKKVVAALIEAANHIDDIAGTLSKLSDLHAHKLRVDPVNFKLLGQCFLVVVAIHHPSALTPEVHASLDKFLCAVGNVLTAKYR",
        )
        self.assertEqual(
            alignment.sequences[17].letter_annotations["posterior probability"],
            "69***********************************************************************99*****************************************************************7",
        )
        self.assertEqual(alignment.sequences[18].id, "HBA_PHACO")
        self.assertEqual(
            alignment.sequences[18].seq,
            "VLSAADKNNVKGIFTKIAGHAEEYGAEALERMFITYPSTKTYFPHFDLSHGSAQIKGHGKKVVAALIEAVNHIDDITGTLSKLSDLHAHKLRVDPVNFKLLGQCFLVVVAIHHPSALTPEVHASLDKFLCAVGTVLTAKYR",
        )
        self.assertEqual(
            alignment.sequences[18].letter_annotations["posterior probability"],
            "69***********************************************************************99*****************************************************************7",
        )
        self.assertEqual(alignment.sequences[19].id, "HBA_TRIOC")
        self.assertEqual(
            alignment.sequences[19].seq,
            "VLSANDKTNVKTVFTKITGHAEDYGAETLERMFITYPPTKTYFPHFDLHHGSAQIKAHGKKVVGALIEAVNHIDDIAGALSKLSDLHAQKLRVDPVNFKLLGQCFLVVVAIHHPSVLTPEVHASLDKFLCAVGNVLSAKYR",
        )
        self.assertEqual(
            alignment.sequences[19].letter_annotations["posterior probability"],
            "69********************************************999************************99*****************************************************************7",
        )
        self.assertEqual(alignment.sequences[20].id, "HBA_ANSSE")
        self.assertEqual(
            alignment.sequences[20].seq,
            "VLSAADKGNVKTVFGKIGGHAEEYGAETLQRMFQTFPQTKTYFPHFDLQPGSAQIKAHGKKVAAALVEAANHIDDIAGALSKLSDLHAQKLRVDPVNFKFLGHCFLVVLAIHHPSLLTPEVHASMDKFLCAVATVLTAKYR",
        )
        self.assertEqual(
            alignment.sequences[20].letter_annotations["posterior probability"],
            "69********************************************9999***********************99*****************************************************************7",
        )
        self.assertEqual(alignment.sequences[21].id, "HBA_COLLI")
        self.assertEqual(
            alignment.sequences[21].seq,
            "VLSANDKSNVKAVFAKIGGQAGDLGGEALERLFITYPQTKTYFPHFDLSHGSAQIKGHGKKVAEALVEAANHIDDIAGALSKLSDLHAQKLRVDPVNFKLLGHCFLVVVAVHFPSLLTPEVHASLDKFVLAVGTVLTAKYR",
        )
        self.assertEqual(
            alignment.sequences[21].letter_annotations["posterior probability"],
            "69***********************************************************************99*****************************************************************7",
        )
        self.assertEqual(alignment.sequences[22].id, "HBAD_CHLME")
        self.assertEqual(
            alignment.sequences[22].seq,
            "mLTADDKKLLTQLWEKVAGHQEEFGSEALQRMFLTYPQTKTYFPHFDLHPGSEQVRGHGKKVAAALGNAVKSLDNLSQALSELSNLHAYNLRVDPANFKLLAQCFQVVLATHLGKDYSPEMHAAFDKFLSAVAAVLAEKYR",
        )
        self.assertEqual(
            alignment.sequences[22].letter_annotations["posterior probability"],
            "689*******************************************9999******************************************************************************************7",
        )
        self.assertEqual(alignment.sequences[23].id, "HBAD_PASMO")
        self.assertEqual(
            alignment.sequences[23].seq,
            "mLTAEDKKLIQQIWGKLGGAEEEIGADALWRMFHSYPSTKTYFPHFDLSQGSDQIRGHGKKVVAALSNAIKNLDNLSQALSELSNLHAYNLRVDPVNFKFLSQCLQVSLATRLGKEYSPEVHSAVDKFMSAVASVLAEKYR",
        )
        self.assertEqual(
            alignment.sequences[23].letter_annotations["posterior probability"],
            "699*******************************************9**9******************************************************************************************7",
        )
        self.assertEqual(alignment.sequences[24].id, "HBAZ_HORSE")
        self.assertEqual(
            alignment.sequences[24].seq,
            "sLTKAERTMVVSIWGKISMQADAVGTEALQRLFSSYPQTKTYFPHFDLHEGSPQLRAHGSKVAAAVGDAVKSIDNVAGALAKLSELHAYILRVDPVNFKFLSHCLLVTLASRLPADFTADAHAAWDKFLSIVSSVLTEKYR",
        )
        self.assertEqual(
            alignment.sequences[24].letter_annotations["posterior probability"],
            "689*******************************************9999******************************************************************************************7",
        )
        self.assertEqual(alignment.sequences[25].id, "HBA4_SALIR")
        self.assertEqual(
            alignment.sequences[25].seq,
            "sLSAKDKANVKAIWGKILPKSDEIGEQALSRMLVVYPQTKAYFSHWASVAPGSAPVKKHGITIMNQIDDCVGHMDDLFGFLTKLSELHATKLRVDPTNFKILAHNLIVVIAAYFPAEFTPEIHLSVDKFLQQLALALAEKYR",
        )
        self.assertEqual(
            alignment.sequences[25].letter_annotations["posterior probability"],
            "69********************************************77769************************9*****************************************************************7",
        )
        self.assertEqual(alignment.sequences[26].id, "HBB_ORNAN")
        self.assertEqual(
            alignment.sequences[26].seq,
            "VHLSGGEKSAVTNLWGKVNINELGGEALGRLLVVYPWTQRFFEAFGDLSSAGAVMGNPKVKAHGAKVLTSFGDALKNLDDLKGTFAKLSELHCDKLHVDPENFNRLGNVLIVVLARHFSKDFSPEVQAAWQKLVSGVAHALGHKYH",
        )
        self.assertEqual(
            alignment.sequences[26].letter_annotations["posterior probability"],
            "69****************************************************************************9******************************************************************7",
        )
        self.assertEqual(alignment.sequences[27].id, "HBB_TACAC")
        self.assertEqual(
            alignment.sequences[27].seq,
            "VHLSGSEKTAVTNLWGHVNVNELGGEALGRLLVVYPWTQRFFESFGDLSSADAVMGNAKVKAHGAKVLTSFGDALKNLDNLKGTFAKLSELHCDKLHVDPENFNRLGNVLVVVLARHFSKEFTPEAQAAWQKLVSGVSHALAHKYH",
        )
        self.assertEqual(
            alignment.sequences[27].letter_annotations["posterior probability"],
            "69***********************************************************************************************************************************************7",
        )
        self.assertEqual(alignment.sequences[28].id, "HBE_PONPY")
        self.assertEqual(
            alignment.sequences[28].seq,
            "VHFTAEEKAAVTSLWSKMNVEEAGGEALGRLLVVYPWTQRFFDSFGNLSSPSAILGNPKVKAHGKKVLTSFGDAIKNMDNLKTTFAKLSELHCDKLHVDPENFKLLGNVMVIILATHFGKEFTPEVQAAWQKLVSAVAIALAHKYH",
        )
        self.assertEqual(
            alignment.sequences[28].letter_annotations["posterior probability"],
            "5789*********************************************************************************************************************************************7",
        )
        self.assertEqual(alignment.sequences[29].id, "HBB_SPECI")
        self.assertEqual(
            alignment.sequences[29].seq,
            "VHLSDGEKNAISTAWGKVHAAEVGAEALGRLLVVYPWTQRFFDSFGDLSSASAVMGNAKVKAHGKKVIDSFSNGLKHLDNLKGTFASLSELHCDKLHVDPENFKLLGNMIVIVMAHHLGKDFTPEAQAAFQKVVAGVANALAHKYH",
        )
        self.assertEqual(
            alignment.sequences[29].letter_annotations["posterior probability"],
            "69****************99*****************************************************************************************************************************7",
        )
        self.assertEqual(alignment.sequences[30].id, "HBB_SPETO")
        self.assertEqual(
            alignment.sequences[30].seq,
            "VHLTDGEKNAISTAWGKVNAAEIGAEALGRLLVVYPWTQRFFDSFGDLSSASAVMGNAKVKAHGKKVIDSFSNGLKHLDNLKGTFASLSELHCDKLHVDPENFKLLGNMIVIVMAHHLGKDFTPEAQAAFQKVVAGVANALSHKYH",
        )
        self.assertEqual(
            alignment.sequences[30].letter_annotations["posterior probability"],
            "69****************99*****************************************************************************************************************************7",
        )
        self.assertEqual(alignment.sequences[31].id, "HBB_EQUHE")
        self.assertEqual(
            alignment.sequences[31].seq,
            "vQLSGEEKAAVLALWDKVNEEEVGGEALGRLLVVYPWTQRFFDSFGDLSNPAAVMGNPKVKAHGKKVLHSFGEGVHHLDNLKGTFAQLSELHCDKLHVDPENFRLLGNVLVVVLARHFGKDFTPELQASYQKVVAGVANALAHKYH",
        )
        self.assertEqual(
            alignment.sequences[31].letter_annotations["posterior probability"],
            "579***************99*****************************************************************************************************************************7",
        )
        self.assertEqual(alignment.sequences[32].id, "HBB_SUNMU")
        self.assertEqual(
            alignment.sequences[32].seq,
            "VHLSGEEKACVTGLWGKVNEDEVGAEALGRLLVVYPWTQRFFDSFGDLSSASAVMGNPKVKAHGKKVLHSLGEGVANLDNLKGTFAKLSELHCDKLHVDPENFRLLGNVLVVVLASKFGKEFTPPVQAAFQKVVAGVANALAHKYH",
        )
        self.assertEqual(
            alignment.sequences[32].letter_annotations["posterior probability"],
            "69****************99*****************************************************************************************************************************7",
        )
        self.assertEqual(alignment.sequences[33].id, "HBB_CALAR")
        self.assertEqual(
            alignment.sequences[33].seq,
            "VHLTGEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMNNPKVKAHGKKVLGAFSDGLTHLDNLKGTFAHLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPVVQAAYQKVVAGVANALAHKYH",
        )
        self.assertEqual(
            alignment.sequences[33].letter_annotations["posterior probability"],
            "689**********************************************************************************************************************************************7",
        )
        self.assertEqual(alignment.sequences[34].id, "HBB_MANSP")
        self.assertEqual(
            alignment.sequences[34].seq,
            "VHLTPEEKTAVTTLWGKVNVDEVGGEALGRLLVVYPWTQRFFDSFGDLSSPDAVMGNPKVKAHGKKVLGAFSDGLNHLDNLKGTFAQLSELHCDKLHVDPENFKLLGNVLVCVLAHHFGKEFTPQVQAAYQKVVAGVANALAHKYH",
        )
        self.assertEqual(
            alignment.sequences[34].letter_annotations["posterior probability"],
            "69***********************************************************************************************************************************************7",
        )
        self.assertEqual(alignment.sequences[35].id, "HBB_URSMA")
        self.assertEqual(
            alignment.sequences[35].seq,
            "VHLTGEEKSLVTGLWGKVNVDEVGGEALGRLLVVYPWTQRFFDSFGDLSSADAIMNNPKVKAHGKKVLNSFSDGLKNLDNLKGTFAKLSELHCDKLHVDPENFKLLGNVLVCVLAHHFGKEFTPQVQAAYQKVVAGVANALAHKYH",
        )
        self.assertEqual(
            alignment.sequences[35].letter_annotations["posterior probability"],
            "689**********************************************************************************************************************************************7",
        )
        self.assertEqual(alignment.sequences[36].id, "HBB_RABIT")
        self.assertEqual(
            alignment.sequences[36].seq,
            "VHLSSEEKSAVTALWGKVNVEEVGGEALGRLLVVYPWTQRFFESFGDLSSANAVMNNPKVKAHGKKVLAAFSEGLSHLDNLKGTFAKLSELHCDKLHVDPENFRLLGNVLVIVLSHHFGKEFTPQVQAAYQKVVAGVANALAHKYH",
        )
        self.assertEqual(
            alignment.sequences[36].letter_annotations["posterior probability"],
            "69***********************************************************************************************************************************************7",
        )
        self.assertEqual(alignment.sequences[37].id, "HBB_TUPGL")
        self.assertEqual(
            alignment.sequences[37].seq,
            "VHLSGEEKAAVTGLWGKVDLEKVGGQSLGSLLIVYPWTQRFFDSFGDLSSPSAVMSNPKVKAHGKKVLTSFSDGLNHLDNLKGTFAKLSELHCDKLHVDPENFRLLGNVLVRVLACNFGPEFTPQVQAAFQKVVAGVANALAHKYH",
        )
        self.assertEqual(
            alignment.sequences[37].letter_annotations["posterior probability"],
            "69***********************************************************************************************************************************************7",
        )
        self.assertEqual(alignment.sequences[38].id, "HBB_TRIIN")
        self.assertEqual(
            alignment.sequences[38].seq,
            "VHLTPEEKALVIGLWAKVNVKEYGGEALGRLLVVYPWTQRFFEHFGDLSSASAIMNNPKVKAHGEKVFTSFGDGLKHLEDLKGAFAELSELHCDKLHVDPENFRLLGNVLVCVLARHFGKEFSPEAQAAYQKVVAGVANALAHKYH",
        )
        self.assertEqual(
            alignment.sequences[38].letter_annotations["posterior probability"],
            "69***************************************************************************989*****************************************************************7",
        )
        self.assertEqual(alignment.sequences[39].id, "HBB_COLLI")
        self.assertEqual(
            alignment.sequences[39].seq,
            "vHWSAEEKQLITSIWGKVNVADCGAEALARLLIVYPWTQRFFSSFGNLSSATAISGNPNVKAHGKKVLTSFGDAVKNLDNIKGTFAQLSELHCDKLHVDPENFRLLGDILVIILAAHFGKDFTPECQAAWQKLVRVVAHALARKYH",
        )
        self.assertEqual(
            alignment.sequences[39].letter_annotations["posterior probability"],
            "5779*********************************************************************************************************************************************7",
        )
        self.assertEqual(alignment.sequences[40].id, "HBB_LARRI")
        self.assertEqual(
            alignment.sequences[40].seq,
            "vHWSAEEKQLITGLWGKVNVADCGAEALARLLIVYPWTQRFFASFGNLSSPTAINGNPMVRAHGKKVLTSFGEAVKNLDNIKNTFAQLSELHCDKLHVDPENFRLLGDILIIVLAAHFAKDFTPDSQAAWQKLVRVVAHALARKYH",
        )
        self.assertEqual(
            alignment.sequences[40].letter_annotations["posterior probability"],
            "5779*********************************************************************************************************************************************7",
        )
        self.assertEqual(alignment.sequences[41].id, "HBB1_VAREX")
        self.assertEqual(
            alignment.sequences[41].seq,
            "vHWTAEEKQLICSLWGKIDVGLIGGETLAGLLVIYPWTQRQFSHFGNLSSPTAIAGNPRVKAHGKKVLTSFGDAIKNLDNIKDTFAKLSELHCDKLHVDPTNFKLLGNVLVIVLADHHGKEFTPAHHAAYQKLVNVVSHSLARRYH",
        )
        self.assertEqual(
            alignment.sequences[41].letter_annotations["posterior probability"],
            "66799********************************************************************************************************************************************7",
        )
        self.assertEqual(alignment.sequences[42].id, "HBB2_XENTR")
        self.assertEqual(
            alignment.sequences[42].seq,
            "vHWTAEEKATIASVWGKVDIEQDGHDALSRLLVVYPWTQRYFSSFGNLSNVSAVSGNVKVKAHGNKVLSAVGSAIQHLDDVKSHLKGLSKSHAEDLHVDPENFKRLADVLVIVLAAKLGSAFTPQVQAVWEKLNATLVAALSHGYf",
        )
        self.assertEqual(
            alignment.sequences[42].letter_annotations["posterior probability"],
            "66799*************************************************************************99*************************************************************99889",
        )
        self.assertEqual(alignment.sequences[43].id, "HBBL_RANCA")
        self.assertEqual(
            alignment.sequences[43].seq,
            "vHWTAEEKAVINSVWQKVDVEQDGHEALTRLFIVYPWTQRYFSTFGDLSSPAAIAGNPKVHAHGKKILGAIDNAIHNLDDVKGTLHDLSEEHANELHVDPENFRRLGEVLIVVLGAKLGKAFSPQVQHVWEKFIAVLVDALSHSYH",
        )
        self.assertEqual(
            alignment.sequences[43].letter_annotations["posterior probability"],
            "66799*************************************************************************99*****************************************************************7",
        )
        self.assertEqual(alignment.sequences[44].id, "HBB2_TRICR")
        self.assertEqual(
            alignment.sequences[44].seq,
            "VHLTAEDRKEIAAILGKVNVDSLGGQCLARLIVVNPWSRRYFHDFGDLSSCDAICRNPKVLAHGAKVMRSIVEATKHLDNLREYYADLSVTHSLKFYVDPENFKLFSGIVIVCLALTLQTDFSCHKQLAFEKLMKGVSHALGHGY",
        )
        self.assertEqual(
            alignment.sequences[44].letter_annotations["posterior probability"],
            "69*******************************************************************************************************************************************9988",
        )
        self.assertEqual(
            alignment.column_annotations["reference coordinate annotation"],
            ".xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx......",
        )
        self.assertEqual(
            alignment.column_annotations["consensus posterior probability"],
            ".679*****************************************************************************99******************************************************************7......",
        )

    def check_alignment_pfam1(self, alignment):
        """Check the alignment obtained by parsing Pfam record 120_Rick_ant."""
        self.assertEqual(alignment.annotations["identifier"], "120_Rick_ant")
        self.assertEqual(alignment.annotations["accession"], "PF12574.10")
        self.assertEqual(
            alignment.annotations["definition"], "120 KDa Rickettsia surface antigen"
        )
        self.assertEqual(alignment.annotations["author"], ["Gavin OL;"])
        self.assertEqual(alignment.annotations["source of seed"], "Prosite")
        self.assertEqual(alignment.annotations["gathering method"], "25.00 25.00;")
        self.assertEqual(alignment.annotations["trusted cutoff"], "42.00 39.60;")
        self.assertEqual(alignment.annotations["noise cutoff"], "23.60 21.20;")
        self.assertEqual(
            alignment.annotations["build method"], "hmmbuild HMM.ann SEED.ann"
        )
        self.assertEqual(
            alignment.annotations["search method"],
            "hmmsearch -Z 57096847 -E 1000 --cpu 4 HMM pfamseq",
        )
        self.assertEqual(alignment.annotations["type"], "Family")
        self.assertEqual(alignment.annotations["references"][0]["number"], 1)
        self.assertEqual(alignment.annotations["references"][0]["medline"], "8112862")
        self.assertEqual(
            alignment.annotations["references"][0]["title"],
            "Cloning, sequencing, and expression of the gene coding for an antigenic 120-kilodalton protein of Rickettsia conorii.",
        )
        self.assertEqual(
            alignment.annotations["references"][0]["author"], "Schuenke KW, Walker DH;"
        )
        self.assertEqual(
            alignment.annotations["references"][0]["location"],
            "Infect Immun. 1994;62:904-909.",
        )
        self.assertEqual(len(alignment.annotations["database references"]), 2)
        self.assertEqual(
            alignment.annotations["database references"][0],
            {"reference": "INTERPRO; IPR020954;"},
        )
        self.assertEqual(
            alignment.annotations["database references"][1],
            {"reference": "SO; 0100021; polypeptide_conserved_region;"},
        )
        self.assertEqual(
            alignment.annotations["comment"],
            "This domain family is found in bacteria, and is approximately 40 amino acids in length. This family is a Rickettsia surface antigen of 120 KDa which may be used as an antigen for immune response against the bacterial species.",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[0, 8, 8, 229, 231, 235], [0, 8, 13, 234, 234, 238]]),
            )
        )
        self.assertEqual(alignment.sequences[0].id, "SCA4_RICPR/103-337")
        self.assertEqual(
            alignment.sequences[0].seq,
            "LAEQIAKEEDDRKFRAFLSNQDNYALINKAFEDTKTKKNLEKAEIVGYKNVLSTYSVANGYQGGFQPVQWENQVSASDLRSTVVKNDEGEELCTLNETTVKTKDLIVAKQDGTQVQINSYREINFPIKLDKANGSMHLSMVALKADGTKPAKDKAVYFTAHYEEGPNGKPQLKEISSPQPLKFVGTGDDAVAYIEHGGEIYTLAVTRGKYKEMMKEVALNHGQSVALSQTIAEDL",
        )
        self.assertEqual(alignment.sequences[0].annotations["accession"], "Q9ZD49.2")
        self.assertEqual(alignment.sequences[1].id, "H8K5G2_RICAG/113-350")
        self.assertEqual(
            alignment.sequences[1].seq,
            "LAEQKRKEIEEEKEKDKTLSTFFGNPANREFIDKALENPELKKKLESIEIAGYKNVHNTFSAASGYPGGFKPVQWENQVSANDLRATVVKNDAGDELCTLNETTVKTKPFTVAKQDGTQVQISSYREIDFPIKLDKADGSMHLSMVALKADGTKPSKDKAVYFTAHYEEGPNGKPQLKEISSPKPLKFAGTGDDAIAYIEHGGEIYTLAVTRGKYKEMMKEVELNQGQSVDLSQAEDI",
        )
        self.assertEqual(alignment.sequences[1].annotations["accession"], "H8K5G2.1")
        self.assertEqual(
            alignment[0],
            "LAEQIAKE-----EDDRKFRAFLSNQDNYALINKAFEDTKTKKNLEKAEIVGYKNVLSTYSVANGYQGGFQPVQWENQVSASDLRSTVVKNDEGEELCTLNETTVKTKDLIVAKQDGTQVQINSYREINFPIKLDKANGSMHLSMVALKADGTKPAKDKAVYFTAHYEEGPNGKPQLKEISSPQPLKFVGTGDDAVAYIEHGGEIYTLAVTRGKYKEMMKEVALNHGQSVALSQTIAEDL",
        )
        self.assertEqual(
            alignment[1],
            "LAEQKRKEIEEEKEKDKTLSTFFGNPANREFIDKALENPELKKKLESIEIAGYKNVHNTFSAASGYPGGFKPVQWENQVSANDLRATVVKNDAGDELCTLNETTVKTKPFTVAKQDGTQVQISSYREIDFPIKLDKADGSMHLSMVALKADGTKPSKDKAVYFTAHYEEGPNGKPQLKEISSPKPLKFAGTGDDAIAYIEHGGEIYTLAVTRGKYKEMMKEVELNQGQSVDLSQ--AEDI",
        )
        self.assertEqual(
            alignment.column_annotations["consensus sequence"],
            "LAEQhtKE.....EcD+phpsFhuN.sNhthIsKAhEsschKKpLEphEIsGYKNVhsTaSsAsGY.GGFpPVQWENQVSAsDLRuTVVKNDtG-ELCTLNETTVKTKshhVAKQDGTQVQIsSYREIsFPIKLDKAsGSMHLSMVALKADGTKPuKDKAVYFTAHYEEGPNGKPQLKEISSPpPLKFsGTGDDAlAYIEHGGEIYTLAVTRGKYKEMMKEVtLNpGQSVsLSQ..AEDl",
        )

    def check_alignment_pfam2(self, alignment):
        """Check the alignment obtained by parsing Pfam record 7kD_DNA_binding."""
        self.assertEqual(alignment.annotations["identifier"], "7kD_DNA_binding")
        self.assertEqual(alignment.annotations["accession"], "PF02294.20")
        self.assertEqual(alignment.annotations["definition"], "7kD DNA-binding domain")
        self.assertEqual(len(alignment.annotations["author"]), 2)
        self.assertEqual(
            alignment.annotations["author"][0], "Mian N;0000-0003-4284-4749"
        )
        self.assertEqual(
            alignment.annotations["author"][1], "Bateman A;0000-0002-6982-4660"
        )
        self.assertEqual(
            alignment.annotations["source of seed"], "Pfam-B_8148 (release 5.2)"
        )
        self.assertEqual(alignment.annotations["gathering method"], "25.00 25.00;")
        self.assertEqual(alignment.annotations["trusted cutoff"], "26.60 46.20;")
        self.assertEqual(alignment.annotations["noise cutoff"], "23.20 19.20;")
        self.assertEqual(
            alignment.annotations["build method"], "hmmbuild HMM.ann SEED.ann"
        )
        self.assertEqual(
            alignment.annotations["search method"],
            "hmmsearch -Z 57096847 -E 1000 --cpu 4 HMM pfamseq",
        )
        self.assertEqual(alignment.annotations["type"], "Domain")
        self.assertEqual(alignment.annotations["clan"], "CL0049")
        self.assertEqual(len(alignment.annotations["references"]), 1)
        self.assertEqual(alignment.annotations["references"][0]["number"], 1)
        self.assertEqual(alignment.annotations["references"][0]["medline"], "3130377")
        self.assertEqual(
            alignment.annotations["references"][0]["title"],
            "Microsequence analysis of DNA-binding proteins 7a, 7b, and 7e from the archaebacterium Sulfolobus acidocaldarius.",
        )
        self.assertEqual(
            alignment.annotations["references"][0]["author"],
            "Choli T, Wittmann-Liebold B, Reinhardt R;",
        )
        self.assertEqual(
            alignment.annotations["references"][0]["location"],
            "J Biol Chem 1988;263:7087-7093.",
        )
        self.assertEqual(len(alignment.annotations["database references"]), 3)
        self.assertEqual(
            alignment.annotations["database references"][0]["reference"],
            "INTERPRO; IPR003212;",
        )
        self.assertEqual(
            alignment.annotations["database references"][1]["reference"],
            "SCOP; 1sso; fa;",
        )
        self.assertEqual(
            alignment.annotations["database references"][2]["reference"],
            "SO; 0000417; polypeptide_domain;",
        )
        self.assertEqual(
            alignment.annotations["comment"],
            "This family contains members of the hyper-thermophilic archaebacterium  7kD DNA-binding/endoribonuclease P2 family. There are five 7kD DNA-binding proteins, 7a-7e, found as monomers in the cell. Protein 7e shows the  tightest DNA-binding ability.",
        )
        self.assertEqual(alignment.sequences[0].id, "DN7_METS5/4-61")
        self.assertEqual(alignment.sequences[0].annotations["accession"], "A4YEA2.1")
        self.assertEqual(alignment.sequences[1].id, "DN7A_SACS2/3-61")
        self.assertEqual(alignment.sequences[1].annotations["accession"], "P61991.2")
        self.assertEqual(len(alignment.sequences[1].dbxrefs), 4)
        self.assertEqual(alignment.sequences[1].dbxrefs[0], "PDB; 1SSO A; 2-60;")
        self.assertEqual(alignment.sequences[1].dbxrefs[1], "PDB; 1JIC A; 2-60;")
        self.assertEqual(alignment.sequences[1].dbxrefs[2], "PDB; 2CVR A; 2-60;")
        self.assertEqual(alignment.sequences[1].dbxrefs[3], "PDB; 1B4O A; 2-60;")
        self.assertEqual(alignment.sequences[2].id, "DN7E_SULAC/3-60")
        self.assertEqual(alignment.sequences[2].annotations["accession"], "P13125.2")
        self.assertEqual(
            alignment.sequences[0].seq,
            "KIKFKYKGQDLEVDISKVKKVWKVGKMVSFTYDDNGKTGRGAVSEKDAPKELLNMIGK",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "TVKFKYKGEEKQVDISKIKKVWRVGKMISFTYDEGGGKTGRGAVSEKDAPKELLQMLEK",
        )
        self.assertEqual(
            alignment.sequences[2].seq,
            "KVRFKYKGEEKEVDTSKIKKVWRVGKMVSFTYDDNGKTGRGAVSEKDAPKELMDMLAR",
        )
        self.assertEqual(
            alignment[0], "KIKFKYKGQDLEVDISKVKKVWKVGKMVSFTYDD-NGKTGRGAVSEKDAPKELLNMIGK"
        )
        self.assertEqual(
            alignment[1], "TVKFKYKGEEKQVDISKIKKVWRVGKMISFTYDEGGGKTGRGAVSEKDAPKELLQMLEK"
        )
        self.assertEqual(
            alignment[2], "KVRFKYKGEEKEVDTSKIKKVWRVGKMVSFTYDD-NGKTGRGAVSEKDAPKELMDMLAR"
        )
        self.assertEqual(
            alignment.sequences[1].letter_annotations["secondary structure"],
            "EEEEESSSSEEEEETTTEEEEEESSSSEEEEEE-SSSSEEEEEEETTTS-CHHHHHHTT",
        )
        self.assertEqual(
            alignment.column_annotations["consensus secondary structure"],
            "EEEEESSSSEEEEETTTEEEEEESSSSEEEEEE-SSSSEEEEEEETTTS-CHHHHHHTT",
        )
        self.assertEqual(
            alignment.column_annotations["consensus sequence"],
            "KVKFKYKGEEKEVDISKIKKVWRVGKMVSFTYDD.NGKTGRGAVSEKDAPKELLsMLuK",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[0, 34, 34, 58], [0, 34, 35, 59], [0, 34, 34, 58]]),
            )
        )

    def check_alignment_pfam3(self, alignment):
        """Check the alignment obtained by parsing Pfam record 12TM_1."""
        self.assertEqual(alignment.annotations["identifier"], "12TM_1")
        self.assertEqual(alignment.annotations["accession"], "PF09847.11")
        self.assertEqual(
            alignment.annotations["definition"], "Membrane protein of 12 TMs"
        )
        self.assertEqual(alignment.annotations["previous identifier"], "DUF2074;")
        self.assertEqual(len(alignment.annotations["author"]), 3)
        self.assertEqual(alignment.annotations["author"][0], "COGs;")
        self.assertEqual(
            alignment.annotations["author"][1], "Finn RD;0000-0001-8626-2148"
        )
        self.assertEqual(
            alignment.annotations["author"][2], "Sammut SJ;0000-0003-4472-904X"
        )
        self.assertEqual(alignment.annotations["source of seed"], "COGs (COG3368)")
        self.assertEqual(alignment.annotations["gathering method"], "33.20 33.20;")
        self.assertEqual(alignment.annotations["trusted cutoff"], "33.60 33.20;")
        self.assertEqual(
            alignment.annotations["build method"], "hmmbuild HMM.ann SEED.ann"
        )
        self.assertEqual(
            alignment.annotations["search method"],
            "hmmsearch -Z 57096847 -E 1000 --cpu 4 HMM pfamseq",
        )
        self.assertEqual(alignment.annotations["noise cutoff"], "33.10 32.90;")
        self.assertEqual(alignment.annotations["type"], "Family")
        self.assertEqual(alignment.annotations["clan"], "CL0181")
        self.assertEqual(len(alignment.annotations["database references"]), 2)
        self.assertEqual(
            alignment.annotations["database references"][0]["reference"],
            "INTERPRO; IPR018646;",
        )
        self.assertEqual(
            alignment.annotations["database references"][1]["reference"],
            "SO; 0100021; polypeptide_conserved_region;",
        )
        self.assertEqual(
            alignment.annotations["comment"],
            "This family carries twelve transmembrane regions. It does not have any characteristic nucleotide-binding-domains of the GxSGSGKST type. so it may not be an ATP-binding cassette transporter. However, it may well be a transporter of some description.  ABC transporters always have two nucleotide binding domains; this has two unusual conserved sequence-motifs: 'KDhKxhhR' and 'LxxLP'.",
        )
        self.assertEqual(alignment.sequences[0].id, "O29855_ARCFU/39-477")
        self.assertEqual(alignment.sequences[1].id, "O29125_ARCFU/30-435")
        self.assertEqual(alignment.sequences[2].id, "Q8U2D3_PYRFU/39-485")
        self.assertEqual(alignment.sequences[3].id, "Q5JDA6_THEKO/35-482")
        self.assertEqual(alignment.sequences[4].id, "Q97VM1_SACS2/39-451")
        self.assertEqual(alignment.sequences[5].id, "Q9HM06_THEAC/17-497")
        self.assertEqual(alignment.sequences[6].id, "Q6L2L5_PICTO/38-510")
        self.assertEqual(alignment.sequences[0].annotations["accession"], "O29855.1")
        self.assertEqual(alignment.sequences[1].annotations["accession"], "O29125.1")
        self.assertEqual(alignment.sequences[2].annotations["accession"], "Q8U2D3.1")
        self.assertEqual(alignment.sequences[3].annotations["accession"], "Q5JDA6.1")
        self.assertEqual(alignment.sequences[4].annotations["accession"], "Q97VM1.1")
        self.assertEqual(alignment.sequences[5].annotations["accession"], "Q9HM06.1")
        self.assertEqual(alignment.sequences[6].annotations["accession"], "Q6L2L5.1")
        self.assertEqual(
            alignment.sequences[0].seq,
            "WIRYNALLLKIMFTFAALFSVGPAFFDDKVSYASSLLSLFFFFLMFGTAYAHGYFQVDLSYMHTFYSRSDISKVRFYGFFRLFDWPAVIALLSLLVLVGMRNPAGLLPALLGFLAVIMGALSIVILLGKRLGSVQTGRSLRAAFFRIFGLIAWLVSIYGLYLINQLAIYLMTFKNYEAYDSLFPISYGLWISQPFSAKYAALSLFYFALITLLFFYAVRELSKEEIAKHYGSLKGWKIKRRGKMTAMVIKDFKQLFRNPQLFVIALLPIYGALMQLVFYIKLSEVASVLYLQIFLAITVSSFMSLERSSYITALPLTDLEMKFSKILEGLLIYFVSMGIVAAVVIYKGGNLINSLSLFPTGFAVVLVAVQFSRRLTSEPVNVEAVIATLISFFIVLVPAAVGGVAVLILKAPFSSYAFPVSLAETLAVLAVFALLNRRK",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "SLRVQVAKSFFIMTFLGSFLCWVAFISSGLGLSLIFTLSLVFSQIYPAQNIAISASSRVFEPLRYLPVRFSERMLVVFFIDSINILAFATPTIAVLMVKNLYFGLYSLLWIIAAILLGYSMVFLYYALFGVKVRSGFSKSVLAGILFFAVLVFALRRFQEIPDLTPYLTPHLLLLSYAASSATIKLSTGRVWRSILNPEIVEVKGSSRLSSGSPLRAMLIKDFRLILRKNALFPLIVPLVIVMPNVVSIANMPNLSIFIITTISTLSTIDLRIIGNLENVDFLRMLPLSKRGFVMSKACLIFVISFAASLPAGSIAFIVSQNPFYLFMAFAIPAIVSMLSSLIIFWQKGEEIYFPEVGFLKWIGLLLVNFGAVYAVLSPRFILSQPVADIISSVLTLLAMTALFEK",
        )
        self.assertEqual(
            alignment.sequences[2].seq,
            "NIIWGVFLQSVMYLGLGVMVAVSILYSENEVQKAIFFSSYLIIPFILTLYSTSLATAYLLSSKAVEPLKPLPLGNLNFIVSLTLLIENLPAFVFLIPASLALGNSIASLLGLLWICSTILMGHSLALFLQIKFSGIHVGKGSVVKTLVKVAGFLIIAGIYFIVQALMRILEDNIEVIAPIFRKYFIAFPFAASTIYEPYKSLVLLALYTLPFLALYFYDLKRLGEVLEGIKTYGKVATKYKLTVANPVTAMFRKDYRIIFRKNPYLGTFLSPLLMSIYFIYNLAKEGFPVMMTLFSIMGISVLGLVMLDPAFAMDREVFPFLSSLPIKRREYLLGKMLTVSLSPLTFSAILVLLSCAFNGTEALLLIPFLASPFLTSSIGILYVKHKMGNERIELPVLKFYDGIVMLILSMIPFIIVAIPLFLLSVPKGYLVSGAIILVGALILSKL",
        )
        self.assertEqual(
            alignment.sequences[3].seq,
            "DLKKTLLFQTAMYAVFGLMLFPSLKGERDAVLVMASTYAILPFIIAFYATVTNSSYIASLDLFKPLLPLPIKLGGRYMSVLLLLESLPVMAFMVPGAVRIGMVVSATSGLLVLLWSAVGLMLGHVFGLLVYYSFGKTSSGRFADLKSLAKALGVILIFGLFYGFSYFQDYVLQNYTSIKESLGGYEFIYPLSVLSVDRPSFSAPLAGIYIAILGVAYYVLISRLWVRISEGSYTSGRRRRAGGLGVYPPELALMVKDFKTALRNTPVLTGLLVPIVIPIINVAGIFSNPDIGAFGGRLATITFVAALGWVSAVSVETLTKIEVKSFELLLSLPLERGRFLRGKLLTMAAIPSAVGVLALLGLSLKGFSSPIYLPMAVLVPLATCGIALHVYYHGTEGLALPQGGILKSLAVWILNAVVVGIIAGSWYLSYPIALLLTAAIDALLLWSL",
        )
        self.assertEqual(
            alignment.sequences[4].seq,
            "NAVTIKISNIIAYTIATIVSASISLINKDAPFSFIFLDLIILANIFTTGLNVIFFVTNYDLKTFLLSLPLSERDVNRAVFRGIFEFFYYGFLASIVIAPISTYMITSSVLQALMAELEIIFFFSLSFALVMLLGKRIRLGITSALFRIGTSLIWIVFIMLPYGLTFKYVTIPTYILPIFPFGFLNIEGLLISLLYTGLSVFFAYKQSLKFLSFRLNSQYSTKYSIKLRSPLITYLYKDIRGLLRVPQASFLLTIPVFALIFSFFAPVYAIFYTIFMITTSSIMLILLEASGMQLLLSLPAGLRSSYISKLLIILIIYLIDVLIFSFFNRASLSLIMLPSTITSVELSLFISYNNVIKGKGMRLADPLSFIIREIEINSIIGIASILTFFANIYYSLLFSVLSLIMINIVVYKK",
        )
        self.assertEqual(
            alignment.sequences[5].seq,
            "YVNISYGATSFSFIVFSLILVAPSLMEHRIYTLSSVVLLLFVYSLFINISNSLLFFVSVNINHILDPLRILPVDFPDHVIAVSWFIYTGSSSLFAVLPAIFLAAFLLGDPYILVIGLIWSIFSVLLGYIIGSSIFVAFGSRISGKRTRSTNILRNVGRIVFLVFVFAIFEIILYNANIVNGIIPRLPYPYSYFIPIFNIQSTVFFFHGIYMQATGFIISMVYTALASFAFIYVNRKAFYRLLEPTARNQSRVKTQMKAEVRSRPFSFFSKDLKISSRKSQNLVLLIMPLFFVFPTIMSEVLYAPTSKADPIILYNAMVAFVIVTSSFYSILFLVIEGNGISFIKALPLDGNSIIRWKISAPTFIFAVISISTLAAISVKALMGAAFYIIIIVDMMLYFVSSTVYNMNRLYRKIPDTADTVNFYSFGGQIAFITTFAFTGLIVGSADIFSLFLQDLLRLNAYFFFLINTVIGIIVLLFMVFR",
        )
        self.assertEqual(
            alignment.sequences[6].seq,
            "TILLYYISNSLSFLFFSIVLNGIYYVKGNTNDISSFGIILFMYIFVIGIYSSLTYINGISINNLLSPVRSLPIKVNTDVPFLSWFIYTGSSYIFIIIPSLLFYYFLVHNLNTIILGLIYAFAMLLFGFIITAIAFIYSSRKPRAHTSLNNFLRILLIFVFLGFFYLIIYDPNILRAYSIYISSLPVYIKYIAFPLNIDYAVYFHPDIIATFFEYLSSFIILLIFFFIYKKIRSRLFYSLEYSEEVKSTEVTRTKIKRDSISVSFIKKDIKITARKSQNLTYILMPIIFVLPFLFTIISSRQPFLSLMFSILSLSILISSFYPIFTLIIENNGILIINALPINRKDIAKYKAYFSMIFYSIIITVVSIIIMAYKNIFNLYYVFIIPDLILIFYTAMIINLNRLIKKIPKGASTINYYSFGVFPTIVLFIVSGIIFGLLISPGIIISEFLYHSIKMSFIFDIIPDLIIFLIMIKK",
        )
        self.assertEqual(
            alignment[0],
            "WIRYNALLLKIMFTFAALFSVGPAFFDDKVS----YASSLLSLFFFFLMFGTAYAHGYFQVDL---SYMHTFYSRSDISKVRFYGFFRLFDWPAVIALLS-----LLVLVGMRNPAGLLPALLGFLAVIMGALSIVILLGKRLGSVQTGR-SLRAAFFRIFGLIAWLVSIYGLYLINQLAI------YLMTFKNYEAYDSLFP-------ISYGLWISQPFSAKYAALSLFYF-ALITLLFFYAVRELSKE----EIAKHYGSLK-GWKIKRRGKMTAMVIKDFKQLFRNPQLFVIALLPIYGALM------------------QLVFYIKLSEVASVLYLQIFLAITVSSFMSLERSSYITALPLTDLEMKFSKILEGLLIYF-VSMGIVAAVVIYKG-GNLINSLSLFPTGFAVVLVAVQFSRRL-------TSEPVNVE---AVIATLISFFIVLVPAAVGGVAVLILKAPFS---SYAFPVSLAETLAVLAVFALLNRRK",
        )
        self.assertEqual(
            alignment[1],
            "SLRVQVAKSFFIMTFLGSFLCWVAFISSGLG-----LSLIFTLSLVFSQIYPAQNIAISASS----RVFEPLRYLPVRFSERMLVVF-FIDSINILAFAT---PTIAVLMVKNLYFGLYSLLWIIAAILLG-YSMVFLYYALFGVKV--RSGFSKSVL--AGILFFAVLVFAL---------------RRFQEIPDLTPYLTP----------------------HLLLLSY--AASSATIKLSTGRVWRSILNPEIVEVKGSSR----LSSGSPLRAMLIKDFRLILRKN-ALFPLIVPLVIVMPNVVSIANMPN--------LSIFIITTISTLSTIDLRIIGNLENVDF--------LRMLPLSKRGFVMSKACLIFVISFAASLPAGSIAFIVS--QNPFYLFMAFAIPAIVSMLSSLIIFWQ-------KGEEIYFPEV-GFLKWIGLLLVNFGAVYAVLSPRFILSQPVA------DIISSVLTL----LAMTALFEK",
        )
        self.assertEqual(
            alignment[2],
            "NIIWGVFLQSVMYLGLGVMVAVSILYSENEVQKAIFFSSYLIIPFILTLYSTSLATAYLLSS----KAVEPLKPLPLGNLNFIVSLTLLIENLPAFVFLI-----PASLALGNSIASLLGLLWICSTILMG-HSLALFLQIKFSGIHVGKGSVVKTLVKVAGFLI----IAGIYFIVQALMRILEDNIEVIAPIFRKYFIAFP--------FAASTIYEPYKS--LVLLALYT-LPFLALYFYDLKRLGEVL---EGIKTYGKVATKYKLTVANPVTAMFRKDYRIIFRKNPYLGTFLSPLLMSIYFIYNLAKEGFPVM-----MTLFSIMGISVLGLVMLDPAFAMDREVF------PFLSSLPIKRREYLLGKMLTVSLSPLTFSAILVLLSCAFNG-TEALLLIPFLASPFLTSSIGILYVKHKM------GNERIELPVL-KFYDGIVMLILSMIPFIIVAIPLFLLSVPKG------YLVSGAIIL----VGALILSKL",
        )
        self.assertEqual(
            alignment[3],
            "DLKKTLLFQTAMYAVFGLML-FPSLKGERDA-VLVMASTYAILPFIIAFYATVTNSSYIASL----DLFKPLLPLPIKLGGRYMSVLLLLESLPVMAFMV--PGAVRIGMVVSATSGLLVLLWSAVGLMLG-HVFGLLVYYSFGKTSSGRFADLKSLAKALGVIL----IFGLFYGFSYFQDYVLQNYTSIKESLGGYEFIYP--------LSVLSVDRPSFS--APLAGIYI-AILGVAYYVLISRLWVRI--SEGSYTSGRRRRAGGLGVYPPELALMVKDFKTALRNTPVLTGLLVPIVIPIINVAGIFSNPDIGAFGGRLATITFVAALGWVSAVSVETLTKIEVKSF------ELLLSLPLERGRFLRGKLLTMAAIPSAVGV-LALLGLSLKGFSSPIYLPMAVLVPLATCGIALHVYYH--------GTEGLALPQG-GILKSLAVWILNAVVVGIIAG-SWYLSYPIA------LLLTAA-------IDALLLWSL",
        )
        self.assertEqual(
            alignment[4],
            "NAVTIKISNIIAYTIATIVSASISLINKDAP----FSFIFLDLIILANIFTTGLNVIFFVTNY---DLKTFLLSLPLSERDVNRAVFRGIFEFFYYGFLA--SIVIAPISTYMITSSVLQALMAELEIIFF-FSLSFALVMLLGKRI--RLGITSALFRIGTSLIWIVFIMLPYGL-----------TFKYVTIPTYILPIFP--------FGFLNIEG------LLISLLYTGLSVFFAYKQSLKFLSFRL--------NSQYSTKYSIKLRSPLITYLYKDIRGLLRVPQASFLLTIPVFALIFSFFAPV------------YAIFYTIFMITTSSIMLIL---LEASGM------QLLLSLPAGLRSSYISKLLIILIIYL-------IDVLIFSFFNRASLSLIMLPSTITSVELSLFISYNNVI-----KGKGMRLA---DPLSFIIREIEINSIIGIASILTFFANIYYS------LLFSVLSLI----MINIVVYKK",
        )
        self.assertEqual(
            alignment[5],
            "YVNISYGATSFSFIVFSLILVAPSLMEHRIY----TLSSVVLLLFVYSLFINISNSLLFFVSVNINHILDPLRILPVDFPDHVIAVSWFIYTGSSSLFAVLPAIFLAAFLLGDPYILVIGLIWSIFSVLLG-YIIGSSIFVAFGSRISGKRTRSTNILRNVGRIVFLVFVFAIFEIILYNANIV---NGIIPRLPYPYSYFIPIFNIQSTVFFFHGIYMQATG--FIISMVYT-ALASFAFIYVNRKAFYRLLEP-TARNQSRVKTQMKAEVRSRPFSFFSKDLKISSRKSQNLVLLIMPLFFVFPTIMSEVLYAPTSKADPIILYNAMVAFVIVTSSFYSILFLVIEGNGI------SFIKALPLDGNSIIRWKISAPTFIFAVISISTLAAISVKAL-MGAAFYIIIIVDMMLYFVSSTVYNMNRLYRKIPDTADTVNFYSFGGQIAFITTFAFTGLIVGSADIFSLFLQDLLRLNAYFFFLINTVIGI----IVLLFMVFR",
        )
        self.assertEqual(
            alignment[6],
            "TILLYYISNSLSFLFFSIVLNGIYYVKGNTN----DISSFGIILFMYIFVIGIYSSLTYINGISINNLLSPVRSLPIKVNTDVPFLSWFIYTGSSYIFIIIPSLLFYYFLVHNLNTIILGLIYAFAMLLFG-FIITAIAFI-----YSSRKPRAHTSLNNFLRILLIFVFLGFFYLIIYDPNILRAYSIYISSLPVYIKYIAFPLNIDYAVYFHPDIIATFFE--YLSSFIIL-LIFFFIYKKIRSRLFYSL--EYSEEVKSTEVTRTKIKRDSISVSFIKKDIKITARKSQNLTYILMPIIFVLPFLFTIISSRQPFLS----LMFSILSLSILISSFYPIFTLIIENNGI------LIINALPINRKDIAKYKAYFSMIFYSIIITVVSIIIMAYKN-IFNLYYVFIIPDLILIFYTAMIINLNRLIKKIPKGASTINYYSF-GVFPTIVLFIVSGIIFGLLISPGIIISEFLYHSIKMSFIFDIIPDL----IIFLIMIKK",
        )
        self.assertEqual(
            alignment.column_annotations["consensus sequence"],
            "slhhthhhpphhahhhulhlsssuhhscphs....hhSohhhL.Flhshahsshsshhhhss....clhcPLhsLPlp.tschhulhhhI.shsshsFhs....hltshhlhs.hsulLsLLauhhsllhG.aslshhlhhhhGthhsuRhuhspslh+.hGhllhhh.lhulahl..h.........hhh.pl.thh.hlhP........hhh.sI.t...t..hlluhlYh.hhhhhhahhshp+Lhhpl....h.cspuphppthplphtu..huhhhKDh+hhhRps.sLshllhPlhhsl..lhs.h............hhlhhlthh.shSslhl.hhhhlEssuh.......hlpuLPlscpphhhuKhhhhhlI.hhhuh.hshhshhhph.tpshhhlhhlssshhsshluhhhshpp.......su-slph..h.uhlshIshhllshlhhulssh.shhLs..hu......hlloss.hl....lhhLlhhc+",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
[[0, 20, 21, 31, 31, 31, 32, 58, 59, 59, 80, 81, 93, 93, 93, 93, 119, 120, 129, 134, 135, 137, 138, 138, 145, 147, 152, 156, 160, 163, 168, 168, 168, 169, 184, 184, 185, 193, 197, 199, 206, 207, 207, 224, 224, 224, 224, 225, 229, 234, 234, 237, 259, 260, 274, 274, 274, 274, 274, 274, 293, 296, 302, 308, 310, 334, 334, 337, 338, 340, 347, 348, 348, 374, 375, 375, 375, 375, 383, 383, 383, 404, 405, 414, 414, 417, 423, 426, 430, 439],
 [0, 20, 21, 31, 31, 31, 31, 57, 57, 57, 78, 78, 90, 90, 90, 92, 118, 118, 127, 132, 133, 133, 134, 135, 142, 142, 147, 151, 155, 155, 155, 155, 155, 155, 170, 170, 170, 170, 170, 170, 177, 177, 177, 194, 195, 197, 198, 199, 203, 208, 208, 208, 230, 230, 244, 250, 254, 254, 254, 254, 273, 276, 282, 282, 282, 306, 307, 310, 311, 313, 320, 320, 320, 346, 347, 347, 347, 347, 355, 357, 357, 378, 379, 388, 388, 388, 394, 397, 397, 406],
 [0, 20, 21, 31, 32, 35, 36, 62, 62, 62, 83, 84, 96, 96, 96, 96, 122, 122, 131, 136, 137, 139, 140, 141, 148, 150, 155, 155, 159, 162, 167, 170, 173, 174, 189, 189, 189, 197, 201, 201, 208, 209, 209, 226, 227, 227, 227, 228, 232, 237, 238, 241, 263, 264, 278, 284, 288, 291, 291, 291, 310, 313, 319, 319, 321, 345, 346, 349, 350, 352, 359, 360, 360, 386, 387, 388, 388, 388, 396, 398, 398, 419, 420, 429, 429, 429, 435, 438, 438, 447],
 [0, 20, 20, 30, 30, 33, 34, 60, 60, 60, 81, 82, 94, 94, 95, 97, 123, 123, 132, 137, 138, 140, 141, 142, 149, 151, 156, 156, 160, 163, 168, 171, 174, 175, 190, 190, 190, 198, 202, 202, 209, 210, 210, 227, 228, 228, 229, 230, 234, 239, 240, 243, 265, 266, 280, 286, 290, 293, 294, 298, 317, 320, 326, 326, 328, 352, 353, 356, 356, 358, 365, 366, 367, 393, 393, 393, 393, 393, 401, 403, 403, 424, 424, 433, 433, 433, 439, 439, 439, 448],
 [0, 20, 21, 31, 31
 , 31, 32, 58, 59, 59, 80, 81, 93, 93, 94, 96, 122, 122, 131, 136, 137, 137, 138, 139, 146, 148, 153, 157, 161, 164, 164, 164, 164, 165, 180, 180, 180, 188, 188, 188, 195, 196, 197, 214, 215, 215, 215, 215, 215, 220, 221, 224, 246, 247, 261, 267, 267, 267, 267, 267, 286, 286, 292, 292, 294, 318, 318, 318, 318, 318, 325, 326, 327, 353, 354, 355, 356, 356, 364, 364, 364, 385, 386, 395, 395, 395, 401, 404, 404, 413],
 [0, 20, 21, 31, 31, 31, 32, 58, 59, 62, 83, 84, 96, 98, 99, 101, 127, 127, 136, 141, 142, 144, 145, 146, 153, 155, 160, 164, 168, 171, 176, 179, 179, 180, 195, 202, 203, 211, 215, 215, 222, 223, 223, 240, 241, 243, 244, 244, 248, 253, 254, 257, 279, 280, 294, 300, 304, 307, 308, 312, 331, 334, 340, 340, 342, 366, 367, 370, 371, 373, 380, 381, 381, 407, 408, 409, 410, 415, 423, 425, 426, 447, 448, 457, 460, 463, 469, 472, 472, 481],
 [0, 20, 21, 31, 31, 31, 32, 58, 59, 62, 83, 84, 96, 98, 99, 101, 127, 127, 136, 136, 137, 139, 140, 141, 148, 150, 155, 159, 163, 166, 171, 174, 177, 178, 193, 200, 201, 209, 213, 213, 220, 221, 221, 238, 239, 239, 240, 241, 245, 250, 251, 254, 276, 277, 291, 297, 301, 304, 305, 305, 324, 327, 333, 333, 335, 359, 360, 363, 364, 366, 373, 374, 374, 400, 401, 402, 403, 408, 416, 418, 418, 439, 440, 449, 452, 455, 461, 464, 464, 473],
]
                    # fmt: on
                ),
            )
        )

    def check_alignment_pfam4(self, alignment):
        """Check the alignment obtained by parsing Pfam record 3Beta_HSD."""
        self.assertEqual(alignment.annotations["identifier"], "3Beta_HSD")
        self.assertEqual(alignment.annotations["accession"], "PF01073.21")
        self.assertEqual(
            alignment.annotations["definition"],
            "3-beta hydroxysteroid dehydrogenase/isomerase family",
        )
        self.assertEqual(len(alignment.annotations["author"]), 2)
        self.assertEqual(
            alignment.annotations["author"][0], "Finn RD;0000-0001-8626-2148"
        )
        self.assertEqual(
            alignment.annotations["author"][1], "Bateman A;0000-0002-6982-4660"
        )
        self.assertEqual(
            alignment.annotations["source of seed"], "Pfam-B_504 (release 3.0)"
        )
        self.assertEqual(alignment.annotations["gathering method"], "22.00 22.00;")
        self.assertEqual(alignment.annotations["trusted cutoff"], "22.00 22.00;")
        self.assertEqual(alignment.annotations["noise cutoff"], "21.90 21.90;")
        self.assertEqual(
            alignment.annotations["build method"], "hmmbuild HMM.ann SEED.ann"
        )
        self.assertEqual(
            alignment.annotations["search method"],
            "hmmsearch -Z 57096847 -E 1000 --cpu 4 HMM pfamseq",
        )
        self.assertEqual(alignment.annotations["type"], "Family")
        self.assertEqual(len(alignment.annotations["wikipedia"]), 1)
        self.assertEqual(
            alignment.annotations["wikipedia"][0], "3β-Hydroxysteroid_dehydrogenase"
        )
        self.assertEqual(alignment.annotations["clan"], "CL0063")
        self.assertEqual(len(alignment.annotations["references"]), 1)
        self.assertEqual(alignment.annotations["references"][0]["number"], 1)
        self.assertEqual(alignment.annotations["references"][0]["medline"], "1562516")
        self.assertEqual(
            alignment.annotations["references"][0]["title"],
            "Structure and tissue-specific expression of 3 beta-hydroxysteroid dehydrogenase/5-ene-4-ene isomerase genes in human and rat classical and peripheral steroidogenic tissues.",
        )
        self.assertEqual(
            alignment.annotations["references"][0]["author"],
            "Labrie F, Simard J, Luu-The V, Pelletier G, Belanger A, Lachance Y, Zhao HF, Labrie C, Breton N, de Launoit Y, et al",
        )
        self.assertEqual(
            alignment.annotations["references"][0]["location"],
            "J Steroid Biochem Mol Biol 1992;41:421-435.",
        )
        self.assertEqual(len(alignment.annotations["database references"]), 3)
        self.assertEqual(
            alignment.annotations["database references"][0]["reference"],
            "INTERPRO; IPR002225;",
        )
        self.assertEqual(
            alignment.annotations["database references"][1]["reference"],
            "HOMSTRAD; Epimerase;",
        )
        self.assertEqual(
            alignment.annotations["database references"][2]["reference"],
            "SO; 0100021; polypeptide_conserved_region;",
        )
        self.assertEqual(
            alignment.annotations["comment"],
            "The enzyme 3 beta-hydroxysteroid dehydrogenase/5-ene-4-ene isomerase (3 beta-HSD) catalyses the oxidation and isomerisation of 5-ene-3 beta-hydroxypregnene and 5-ene-hydroxyandrostene steroid precursors into the corresponding 4-ene-ketosteroids necessary for the formation of all classes of steroid hormones.",
        )
        self.assertEqual(len(alignment.sequences), 8)
        self.assertEqual(alignment.sequences[0].id, "3BHS_FOWPN/7-287")
        self.assertEqual(alignment.sequences[1].id, "Q98318_MCV1/5-278")
        self.assertEqual(alignment.sequences[2].id, "3BHS_VACCW/5-274")
        self.assertEqual(alignment.sequences[3].id, "3BHS1_HUMAN/7-288")
        self.assertEqual(alignment.sequences[4].id, "3BHS1_MOUSE/7-288")
        self.assertEqual(alignment.sequences[5].id, "O22813_ARATH/15-299")
        self.assertEqual(alignment.sequences[6].id, "HSDD3_ARATH/16-287")
        self.assertEqual(alignment.sequences[7].id, "ERG26_YEAST/8-280")
        self.assertEqual(alignment.sequences[0].annotations["accession"], "Q67477.2")
        self.assertEqual(alignment.sequences[1].annotations["accession"], "Q98318.1")
        self.assertEqual(alignment.sequences[2].annotations["accession"], "P26670.1")
        self.assertEqual(alignment.sequences[3].annotations["accession"], "P14060.2")
        self.assertEqual(alignment.sequences[4].annotations["accession"], "P24815.3")
        self.assertEqual(alignment.sequences[5].annotations["accession"], "O22813.1")
        self.assertEqual(alignment.sequences[6].annotations["accession"], "A9X4U2.2")
        self.assertEqual(alignment.sequences[7].annotations["accession"], "P53199.1")
        self.assertEqual(
            alignment.sequences[0].seq,
            "VVTGGCGFLGRHIINNLILFESSLKEVRVYDIRIDQWLLDLVEKCNIIKIVPVIGDVRNKSTLDEALRSADVVIHIASINDVAGKFTNDSIMDVNINGTKNVVDSCLYNGVRVLVYTSSYSAVGPNFLGDAMIRGNENTYYQSNHKEAYPLSKQLSEKYILEANGTMSNIGLRLCTCALRPLGVFGEYCPVLETLYRRSYKSRKMYKYADDKVFHSRVYAGNVAWMHILAARNMIENGQHSPLCNNVYYCYDTSPTEHYHDFNMHFFNQLGMDLRNTCLPL",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "AVTGGGGFIGSYIVRALLQCERTLIELRVIDVRWGTKVSSRNVNVVYIYCDVCDTARLCAALEGVDVLIHTAGLVDVMGEYSEDEIYRANVHGTHSALSACVCAGVRFVVYTSSMEVVGPNMRAEPFVGDEKTEYESCHQHCYPRSKAEAEELVLSSNGRRVRGGQRMLTCALRPPGVYGEGNQLLLRLAKNYVRMGLHVPRTVCENALQSRVYVGNVAWMHVLAARALQEPDSRLPGNAYFCYDHSPCMDYEAFNVMLLRSFGVELGGPRLPR",
        )
        self.assertEqual(
            alignment.sequences[2].seq,
            "AVTGGAGFLGRYIVKLLISADDVQEIRVIDIVEDPQPITSKVKVINYIQCDINDFDKVREALDGVNLIIHTAALVDVFGKYTDNEIMKVNYYGTQTILAACVDLGIKYLIYTSSMEAIGPNKHGDPFIGHEHTLYDISPGHVYAKSKRMAEQLVMKANNSVIMNGAKLYTCCLRPTGIYGEGDKLTKVFYEQCKQHGNIMYRTVDDDAVHSRVYVGNVAWMHVLAAKYIQYPGSEIKGNAYFCYDYSPSCSYDMFNLLLMKPLGIEQGSR",
        )
        self.assertEqual(
            alignment.sequences[3].seq,
            "LVTGAGGFLGQRIIRLLVKEKELKEIRVLDKAFGPELREEFSKLQNKTKLTVLEGDILDEPFLKRACQDVSVIIHTACIIDVFGVTHRESIMNVNVKGTQLLLEACVQASVPVFIYTSSIEVAGPNSYKEIIQNGHEEEPLENTWPAPYPHSKKLAEKAVLAANGWNLKNGGTLYTCALRPMYIYGEGSRFLSASINEALNNNGILSSVGKFSTVNPVYVGNVAWAHILALRALQDPKKAPSIRGQFYYISDDTPHQSYDNLNYTLSKEFGLRLDSRWSFPL",
        )
        self.assertEqual(
            alignment.sequences[4].seq,
            "LVTGAGGFVGQRIIKMLVQEKELQEVRALDKVFRPETKEEFSKLQTKTKVTVLEGDILDAQCLRRACQGISVVIHTAAVIDVTGVIPRQTILDVNLKGTQNLLEACVQASVPAFIFCSSVDVAGPNSYKKIVLNGHEEQNHESTWSDPYPYSKKMAEKAVLAANGSMLKNGGTLNTCALRPMYIYGERSPFIFNAIIRALKNKGILCVTGKFSIANPVYVENVAWAHILAARGLRDPKKSTSIQGQFYYISDDTPHQSYDDLNYTLSKEWGLRPNASWSLPL",
        )
        self.assertEqual(
            alignment.sequences[5].seq,
            "VVTGGLGFVGAALCLELVRRGARQVRSFDLRHSSPWSDDLKNSGVRCIQGDVTKKQDVDNALDGADCVLHLASYGMSGKEMLRFGRCDEVNINGTCNVLEAAFKHEITRIVYVSTYNVVFGGKEILNGNEGLPYFPLDDHVDAYSRTKSIAEQLVLKSNGRPFKNGGKRMYTCAIRPAAIYGPGEDRHLPRIVTLTKLGLALFKIGEPSVKSDWIYVENLVLAIILASMGLLDDIPGREGQPVAAGQPYFVSDGYPVNTFEFLRPLLKSLDYDLPKCTISVPFAL",
        )
        self.assertEqual(
            alignment.sequences[6].seq,
            "VVLGGRGFIGRSLVSRLLRLGNWTVRVADSGHTLHLDESDSLLEDALSSGRASYHCVDVRDKPQIVKVTEGSYVVFYMGATDLRSHDYFDCYKVIVQGTRNVISACRESGVRKLIYNSTADVVFDGSQPIRDGDESLRRPLKFQSMLTDFKAQAEALIKLANNRDGLLTCALRSSIVFGPGDTEFVPFLVNLAKSGYAKFILGSGENISDFTYSENVSHAHICAVKALDSQMEFVAGKEFFITNLKPVRFWDFVSHIVEGLGYPRPSIKLPV",
        )
        self.assertEqual(
            alignment.sequences[7].seq,
            "LIIGGSGFLGLHLIQQFFDINPKPDIHIFDVRDLPEKLSKQFTFNVDDIKFHKGDLTSPDDMENAINESKANVVVHCASPMHGQNPDIYDIVNVKGTRNVIDMCKKCGVNILVYTSSAGVIFNGQDVHNADETWPIPEVPMDAYNETKAIAEDMVLKANDPSSDFYTVALRPAGIFGPGDRQLVPGLRQVAKLGQSKFQIGDNNNLFDWTYAGNVADAHVLAAQKLLDPKTRTAVSGETFFITNDTPTYFWALARTVWKADGHIDKHVIVLKR",
        )
        self.assertEqual(
            alignment[0],
            "VVTGGCGFLGRHIINNLILFESSLKEVRVYD------IRIDQWLLDLVEKCNII-KIVPVIGDVRNKSTLDEALRSADVVIHIASINDVAG-KFTNDSIMDVNINGTKNVVDSCLYNGVRVLVYTSSYSAVGPNFLGDAMIRGNENTYYQSN--HKEAYPLSKQLSEKYILEANG-TMSNIGLRLCTCALRPLGVFGEYCPVLETLYRRSYKSR-KMYKYADDKVFHSRVYAGNVAWMHILAARNMIENGQ----HSPLCNNVYYCYDTSPTEHYHDFNMHFFNQLGMDLRN-T---CLPL",
        )
        self.assertEqual(
            alignment[1],
            "AVTGGGGFIGSYIVRALLQCERTLIELRVID------VRWGTKV--SSRNVNVV----YIYCDVCDTARLCAALEGVDVLIHTAGLVDVMG-EYSEDEIYRANVHGTHSALSACVCAGVRFVVYTSSMEVVGPNMRAEPFV-GDEKTEYESC--HQHCYPRSKAEAEELVLSSNGRRVRGGQ-RMLTCALRPPGVYGEGNQLLLRLAKNYVRMGLHVPRTVCENALQSRVYVGNVAWMHVLAARALQEP------DSRLPGNAYFCYDHSPCMDYEAFNVMLLRSFGVELGGP----RLPR",
        )
        self.assertEqual(
            alignment[2],
            "AVTGGAGFLGRYIVKLLISADD-VQEIRVID------IVEDPQP--ITSKVKVIN---YIQCDINDFDKVREALDGVNLIIHTAALVDVFG-KYTDNEIMKVNYYGTQTILAACVDLGIKYLIYTSSMEAIGPNKHGDPFI-GHEHTLYDIS--PGHVYAKSKRMAEQLVMKANNSVIMNGA-KLYTCCLRPTGIYGEGDKLTKVFYEQCKQHGNIMYRTVDDDAVHSRVYVGNVAWMHVLAAKYIQYP------GSEIKGNAYFCYDYSPSCSYDMFNLLLMKPLGIEQGSR--------",
        )
        self.assertEqual(
            alignment[3],
            "LVTGAGGFLGQRIIRLLVKEKE-LKEIRVLD------KAFGPELREEFSKLQNKTKLTVLEGDILDEPFLKRACQDVSVIIHTACIIDVFG-VTHRESIMNVNVKGTQLLLEACVQASVPVFIYTSSIEVAGPNSYKEIIQNGHEEEPLENT--WPAPYPHSKKLAEKAVLAANGWNLKNGG-TLYTCALRPMYIYGEGSRFLSASINEALNNN-GILSSVGKFSTVNPVYVGNVAWAHILALRALQDPKK----APSIRGQFYYISDDTPHQSYDNLNYTLSKEFGLRLDSRW---SFPL",
        )
        self.assertEqual(
            alignment[4],
            "LVTGAGGFVGQRIIKMLVQEKE-LQEVRALD------KVFRPETKEEFSKLQTKTKVTVLEGDILDAQCLRRACQGISVVIHTAAVIDVTG-VIPRQTILDVNLKGTQNLLEACVQASVPAFIFCSSVDVAGPNSYKKIVLNGHEEQNHEST--WSDPYPYSKKMAEKAVLAANGSMLKNGG-TLNTCALRPMYIYGERSPFIFNAIIRALKNKGILCVTGKFSI-ANPVYVENVAWAHILAARGLRDPKK----STSIQGQFYYISDDTPHQSYDDLNYTLSKEWGLRPNASW---SLPL",
        )
        self.assertEqual(
            alignment[5],
            "VVTGGLGFVGAALCLELVRRG--ARQVRSFD------LRHSSPWSDDLKNSGVR----CIQGDVTKKQDVDNALDGADCVLHLASYGMSGKEMLRFGRCDEVNINGTCNVLEAAFKHEITRIVYVSTYNVVFG---GKEILNGNEGLPYFPLDDHVDAYSRTKSIAEQLVLKSNGRPFKNGGKRMYTCAIRPAAIYGPGEDRHLPRIVTLTKLGLALFKIGEPSVKSDWIYVENLVLAIILASMGLLDDIPGREGQPVAAGQPYFVSDGYPVN-TFEFLRPLLKSLDYDLPKCTISVPFAL",
        )
        self.assertEqual(
            alignment[6],
            "VVLGGRGFIGRSLVSRLLRLGN--WTVRVADSGHTLHLDESDSLLEDALSSGRAS---YHCVDVRDKPQIVKVTEGSYVVFYM-GATDLRS-HDYFD-CYKVIVQGTRNVISACRESGVRKLIYNSTADVVFD--GSQPIRDGDESLRRPLK--FQSMLTDFKAQAEALIKLANN---RDG---LLTCALRSSIVFGPGDTEFVPFLVNLAKSGYAKFILGSGENISDFTYSENVSHAHICAVKALDSQ------MEFVAGKEFFITNLKPVR-FWDFVSHIVEGLGYPRPS-I---KLPV",
        )
        self.assertEqual(
            alignment[7],
            "LIIGGSGFLGLHLIQQFFDINP-KPDIHIFD------VRDLPEKLSKQFTFNVDDI-KFHKGDLTSPDDMENAINESKANVVVHCASPMHG--QNPDIYDIVNVKGTRNVIDMCKKCGVNILVYTSSAGVIFN---GQDVHNADETWPIPEV--PMDAYNETKAIAEDMVLKAND-----PSSDFYTVALRPAGIFGPGDRQLVPGLRQVAKLGQSKFQIGDNNNLFDWTYAGNVADAHVLAAQKLLDPKT----RTAVSGETFFITNDTPTY-FWALARTVWKADGHIDKHVI---VLKR",
        )
        self.assertEqual(
            alignment.column_annotations["consensus sequence"],
            "lVTGGuGFlGppIlptLlptcp.lpElRVhD......lchssphh-chppsslts...hlpGDlpDpsplccAlcGssVlIHsAulsDVtG.hhsp-pIhcVNlpGTpNlL-AClpsGVphllYTSSh-VlGPN.hucsllsGcEpp.apss..atcsYscSKphAEchVLpANG..h+NGu.cLhTCALRPsuIYGEGsphlhshlppshKpG.thaphucssshpshVYVGNVAWAHILAA+uLp-Ph.....poslsGpsYFloDsoPsppYcsFNhpLhKshGhchsu.h...sLPl",
        )

        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
[[  0,  21,  22,  23,  24,  31,  31,  38,  40,  48,  48,  49,  50, 51,  76,  77,  84,  84,  85,  89,  90, 125, 127, 128, 133, 134, 144, 144, 165, 165, 167, 169, 170, 171, 172, 173, 203, 203, 213, 214, 237, 239, 239, 257, 258, 276, 276, 277, 277, 281],
 [  0,  21,  22,  23,  24,  31,  31,  38,  38,  46,  46,  46,  46, 46,  71,  72,  79,  79,  80,  84,  85, 120, 122, 123, 128, 128, 138, 138, 159, 160, 162, 164, 165, 166, 166, 167, 197, 198, 208, 209, 232, 232, 232, 250, 251, 269, 270, 270, 270, 274],
 [  0,  21,  22,  22,  23,  30,  30,  37,  37,  45,  46,  46,  46, 46,  71,  72,  79,  79,  80,  84,  85, 120, 122, 123, 128, 128, 138, 138, 159, 160, 162, 164, 165, 166, 166, 167, 197, 198, 208, 209, 232, 232, 232, 250, 251, 269, 270, 270, 270, 270],
 [  0,  21,  22,  22,  23,  30,  30,  37,  39,  47,  48,  49,  50, 51,  76,  77,  84,  84,  85,  89,  90, 125, 127, 128, 133, 134, 144, 144, 165, 166, 168, 170, 171, 172, 172, 173, 203, 203, 213, 214, 237, 239, 239, 257, 258, 276, 277, 278, 278, 282],
 [  0,  21,  22,  22,  23,  30,  30,  37,  39,  47,  48,  49,  50, 51,  76,  77,  84,  84,  85,  89,  90, 125, 127, 128, 133, 134, 144, 144, 165, 166, 168, 170, 171, 172, 172, 173, 203, 204, 214, 214, 237, 239, 239, 257, 258, 276, 277, 278, 278, 282],
 [  0,  21,  21,  21,  22,  29,  29,  36,  38,  46,  46,  46,  46, 46,  71,  72,  79,  80,  81,  85,  86, 121, 121, 121, 126, 127, 137, 139, 160, 161, 163, 165, 166, 167, 168, 169, 199, 200, 210, 211, 234, 236, 240, 258, 258, 276, 277, 278, 281, 285],
 [  0,  21,  22,  22,  22,  29,  35,  42,  44,  52,  53,  53,  53, 53,  78,  78,  85,  85,  86,  90,  90, 125, 125, 126, 131, 132, 142, 142, 163, 163, 163, 165, 166, 166, 166, 166, 196, 197, 207, 208, 231, 231, 231, 249, 249, 267, 267, 268, 268, 272],
 [  0,  21,  22,  22,  23,  30,  30,  37,  39,  47,  48,  49,  49, 50,  75,  76,  83,  83,  83,  87,  88, 123, 123, 123, 128, 129, 139, 139, 160, 160, 160, 160, 161, 162, 163, 164, 194, 195, 205, 206, 229, 231, 231, 249, 249, 267, 268, 269, 269, 273]
]
                    # fmt: on
                ),
            )
        )

    def check_alignment_pfam5(self, alignment):
        """Check the alignment obtained by parsing Pfam record ArsP_1."""
        self.assertEqual(alignment.annotations["identifier"], "ArsP_1")
        self.assertEqual(alignment.annotations["accession"], "PF03773.15")
        self.assertEqual(alignment.annotations["definition"], "Predicted permease")
        self.assertEqual(alignment.annotations["previous identifier"], "DUF318;")
        self.assertEqual(
            alignment.annotations["author"], ["Bateman A;0000-0002-6982-4660"]
        )
        self.assertEqual(alignment.annotations["source of seed"], "COG0701")
        self.assertEqual(alignment.annotations["gathering method"], "32.30 32.30;")
        self.assertEqual(alignment.annotations["trusted cutoff"], "32.30 32.40;")
        self.assertEqual(alignment.annotations["noise cutoff"], "32.20 32.20;")
        self.assertEqual(
            alignment.annotations["build method"], "hmmbuild  --handHMM.ann SEED.ann"
        )
        self.assertEqual(
            alignment.annotations["search method"],
            "hmmsearch -Z 57096847 -E 1000 --cpu 4 HMM pfamseq",
        )
        self.assertEqual(alignment.annotations["type"], "Family")
        self.assertEqual(len(alignment.annotations["nested domains"]), 1)
        self.assertEqual(
            alignment.annotations["nested domains"][0]["accession"], "PF04945;"
        )
        self.assertEqual(
            alignment.annotations["nested domains"][0]["location"], "D4GY01.1/189-231"
        )
        self.assertEqual(len(alignment.annotations["database references"]), 3)
        self.assertEqual(
            alignment.annotations["database references"][0]["reference"],
            "INTERPRO; IPR005524;",
        )
        self.assertEqual(
            alignment.annotations["database references"][1]["reference"], "TC; 2.A.119;"
        )
        self.assertEqual(
            alignment.annotations["database references"][2]["reference"],
            "SO; 0100021; polypeptide_conserved_region;",
        )
        self.assertEqual(
            alignment.annotations["comment"],
            "This family of integral membrane proteins are predicted to be permeases of unknown specificity.",
        )
        self.assertEqual(len(alignment.sequences), 11)
        self.assertEqual(alignment.sequences[0].id, "O26980_METTH/26-325")
        self.assertEqual(alignment.sequences[1].id, "O67395_AQUAE/24-315")
        self.assertEqual(alignment.sequences[2].id, "Q9X092_THEMA/31-364")
        self.assertEqual(alignment.sequences[3].id, "O28037_ARCFU/7-346")
        self.assertEqual(alignment.sequences[4].id, "Y584_METJA/16-362")
        self.assertEqual(alignment.sequences[5].id, "Y2963_MYCTU/18-329")
        self.assertEqual(alignment.sequences[6].id, "D4GY01_HALVD/35-380")
        self.assertEqual(alignment.sequences[7].id, "YCGR_BACSU/7-294")
        self.assertEqual(alignment.sequences[8].id, "Q9KCQ1_BACHD/41-335")
        self.assertEqual(alignment.sequences[9].id, "P72867_SYNY3/3-335")
        self.assertEqual(alignment.sequences[10].id, "P73433_SYNY3/6-329")
        self.assertEqual(alignment.sequences[0].annotations["accession"], "O26980.1")
        self.assertEqual(alignment.sequences[1].annotations["accession"], "O67395.1")
        self.assertEqual(alignment.sequences[2].annotations["accession"], "Q9X092.1")
        self.assertEqual(alignment.sequences[3].annotations["accession"], "O28037.1")
        self.assertEqual(alignment.sequences[4].annotations["accession"], "Q58004.3")
        self.assertEqual(alignment.sequences[5].annotations["accession"], "I6YET7.1")
        self.assertEqual(alignment.sequences[6].annotations["accession"], "D4GY01.1")
        self.assertEqual(alignment.sequences[7].annotations["accession"], "P94395.1")
        self.assertEqual(alignment.sequences[8].annotations["accession"], "Q9KCQ1.1")
        self.assertEqual(alignment.sequences[9].annotations["accession"], "P72867.1")
        self.assertEqual(alignment.sequences[10].annotations["accession"], "P73433.1")
        self.assertEqual(
            alignment.sequences[0].seq,
            "HLGSAVNFFIYDTIKIFILLATLIFVISFIRTYIPPNKVKETLEKRHRYTGNFIAALVGIITPFCSCSAVPLFIGFVEAGVPLGATFSFLISSPMINEIAIILLLGLFGWQITAFYILSGFIIAVLGGILIGKLKMETELEDYVYETLEKMRALGVADVELPKPTLRERYVIAKNEMKDILRRVSPYIVIAIAIGGWIHGYLPEDFLLQYAGADNIFAVPMAVIIGVPLYSNAAGTIPLISALIEKGMAAGTALALMMSITALSLPEMIILRKVMKPKLLATFIAILAVSITLTGYIFNL",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "HLAEALHFFVYDTLKIFTLLTVIIFVVSFIRSFFPLEKTREILSKHKAVALPLAAFLGILTPFCSCSAVPMFIGFVEAGIPLGAAFTFLVASPMVNEVALGLLLTLFGVKVAVLYVIFGVIVAIVAGYVIEKLNPRELIADYVFQVKLGQTQIKEMTFKERLEFAKNNVKEILGKIWIYIIIAIGIGGFIHGYVPQDIVERVAKTAGLIAVPLAVLIGIPLYSNAAGILPVIQALIAKGVPLGTALAFMMATTALSFPEFMILKQIMKPKLIAFFAGIVGISIIAVGYLFNF",
        )
        self.assertEqual(
            alignment.sequences[2].seq,
            "ILNGFYLLHEYAREHVLLCLVPAFFIAGTISVMLKKDAVLKLLGPNAKRIISYPVAAISGGILAVCSCTILPLFGGIYKKGAGIGPATTFLFAGPAINIAAIFLTARVLGWDLGLARLIATITAAVLIGLIMEMIYQERGEGGLAFTSDDDQYGVRGIIFFLIQLGFLVTSSLGINQTLKYSLMTLLGISALFMALFGFKRDTVENWLYETWDFAKKILPYLFIGVFFAGVLTRLLPQQVVTALLGSNSFLSNLVASVIGTLMYFATLTEVPIVQALRELGMAKGPTLALLMAGNSLSLPSMIVITKLLGKKKAFTYFGLVVVFSTLFGMIYGV",
        )
        self.assertEqual(
            alignment.sequences[3].seq,
            "LLAGIQALEEYIALHVLTCLVPAFLIAGALMSMMNKAVLINYLGAATSKLKSFPLAIVSSFFLAVCSCTVIPIASGIYKRTNATAPAMIILWVAPATNILAVTYTGAVLGLELALARIVAAISTAFVVGLILFYVFDRKIASQSDSAMPKAGRLVENNALVLFALLVATLLLPNYLGVGKPYIFKVEVFSVLMLVTTVYALKSFSKEDLKYWMLETWFFVKQIIPLLLVGVFIVGVVGEILKATDVVEVYLGGEGVGQSFLAALIGALSYFATMTEAPFVDTLMKLGMGKGPALALLLAGPGLSLPNMLAIGKLFGVKRAAVYIITIVALSTIAGVVYGE",
        )
        self.assertEqual(
            alignment.sequences[4].seq,
            "MINTIIDYLNVNRVLALLMAFLMAGGIASMINKNFIIKYFGSNTPKYISYTVAAVSGSLLAVCSCTILPLFASIYKRGAGIGPATTFLFSGPAINVLAIFYSAALLGWDIGFLRAVFAVVVSILIGLSMEIIFKSHEKKRALRVPKADKISDRPLYQTITFFALQFIMLLVITASPKLFPTLSMPLYDGFLLKHLLFIILGIILAVTTKIWFKDEEIKNWLRESFTLLKIVFPLLIIGVAIAGAIKAIIPPSYIATYVGGNSITANFIASFIGALMYFATLTEVPIIKALMELGMGVGPAMALLLAGPSLSIPTVLTISKVLGKTKALTYLGLVVIFSTICGYIAGI",
        )
        self.assertEqual(
            alignment.sequences[5].seq,
            "IGHALALTASMTWEILWALILGFALSAVVQAVVRRSTIVTLLGDDRPRTLVIATGLGAASSSCSYAAVALARSLFRKGANFTAAMAFEIGSTNLVVELGIILALLMGWQFTAAEFVGGPIMILVLAVLFRLFVGARLIDAAREQAERGLAGSMEGHAAMDMSIKREGSFWRRLLSPPGFTSIAHVFVMEWLAILRDLILGLLIAGAIAAWVPESFWQSFFLANHPAWSAVWGPIIGPIVAIVSFVCSIGNVPLAAVLWNGGISFGGVIAFIFADLLILPILNIYRKYYGARMMLVLLGTFYASMVVAGYLIE",
        )
        self.assertEqual(
            alignment.sequences[6].seq,
            "SAREALSTTAAMAWVTWWALVVGFAIAGGVEAWTSGEEVSELLEGHGPREIGYGSLFGFVSSSCSYSAIATAKNLFKKGGSAAATLGAFMFASTNLVIEIGAVIWILLGWQFLVADILGGFILIGLMAFGFVYLVPDEVVEQARRNVQDEGSETVRDPVCGMEVDPDETEYSVERDGRTFYFCSKSCKESFDPEEANTTVRERATSLSGWKALADKQWKEWGMLWDEIAIGFVFAGLIAGFIPDAVWTSVFSGPTFGLPVYVFWTAVLGAVIGVATFVCSVGNVPFGAVLFSNGLPFGSVLSYIYADLIVPPIVDAYREYYGTTFAAVLSGMIFVAAVLTGVVIHF",
        )
        self.assertEqual(
            alignment.sequences[7].seq,
            "FLQLNSIFISILIEAIPFILIGVILSGIIQMFVSEEMIARIMPKNRFLAVLFGALAGVLFPACECGIIPITRRLLLKGVPLHAGVAFMLTAPIINPIVLFSTYIAFGNRWSVVFYRGGLALAVSLIIGVILSYQFKDNQLLKPDEPGHHHHHHGTLLQKLGGTLRHAIDEFFSVGKYLIIGAFIAAAMQTYVKTSTLLAIGQNDVSSSLVMMGLAFVLSLCSEVDAFIASSFSSTFSLGSLIAFLVFGAMVDIKNLLMMLAAFKKRFVFLLITYIVVIVLAGSLLVKG",
        )
        self.assertEqual(
            alignment.sequences[8].seq,
            "WMNVNTIFLGIVIEAVPFILLGVFVSALIQIYVKEDTIQRYLPKNAYAALLPAAVLGAIFPICECAIVPIVRRLIKKGMPLHVGVVFLVAAPILNPIVAASTYFAFRTDLTVLYARMGLAFILSIVIGGLLYVMFKNSDQLKWTKEELVGRVPVQSDMELKPKMNRLKQTLYHASDEFFLMGKYLIAGAFIAALFQTFLDRNILVTIGSNEWSSTGVMMAFAFILSLCSEADAFVAASFGSTFTTGSLIAFLVYGPMLDLKNTIMLFAFFKSKFVLAFMITVTVVVFLAVMVLQF",
        )
        self.assertEqual(
            alignment.sequences[9].seq,
            "QLHEAFTIFLSLLVEAIPFLTFGVVLSSALLVFSDEKKLIAYIPRNPFLGAIAGSLVGFMFPVCECGNVPVARRFLMQGLPPSVAVAFLLAAPTINPIVIWSTWVAFRDQPGMVVARVVCSLIITVIVSWVFSRQLDAVPLLKPALGRRLAYLTRPEESPTAIACESPLLQSGTFLLGSGNSGQLLKLDEQAVETLLPPIAPSRWEMFTDNIVQELRELGGMLILGSLIAAVIQVFIPREWILLLGQGTISSILAMMLLSVVVSVCSTVDSFFALSFVSTFTSSSLLAFLVFGPMIDVKSIGLLLSVFQRRIVIYLLLLTGQLTFLLSLAHSY",
        )
        self.assertEqual(
            alignment.sequences[10].seq,
            "EFNLFLDLLGSALLLSLPWLLLGIIISSTFLIWTDEQKWVANFPRNRLLSSLVGSALGFLLPLGAFGSVPLVRRLLLQGAPIPLAVSFLVAAPTLNIFAIVRVLSSRQSQYGLIFLCISCSWLMAIVMGLVFSTYRLARQQAEDEGETALLNIPLLRSGALIILQSSMEASPRQGGLVFASGVNPVADFSWRQKLHLFGRNIIEEFQEFGGVLVIGTAIACGIVFFLPQAWLLQWAGLGPVRQTVLMMGWSFILPLGNFSNPDLLAPLGEQLWRGSMVAFLLWGSLFNLQTIGLWLVTLRLRPLSYLVVLVGLSVFLFAMVTNY",
        )
        self.assertEqual(
            alignment[0],
            "HLGSAVNFFIYDTIKIFILLATLIFVISFIRTYIPPNKVKETLE-KRHRYTGNFIAALVGIITPFCSCSAVPLFIGFVEAGVPLGATF-SFLISSPMINEIAIILLLGLFG--WQITAFYILSGFIIAVLGGILIGKLKMETELEDYVYETLE-------------------------KMRALGVADV------------------ELPKPTLR---ERYV--IAKNEMKDILRRVS-------PYIVIAIAIGGWIHGYL-PEDFLLQYA--GADNIF-------AVPMAVIIGVPLYSNAAGTIPLISALIEKGMAAGTALALMMSITALSLPEMIILRKVMKPKLLATFIAILAVSITLTGYIFNL",
        )
        self.assertEqual(
            alignment[1],
            "HLAEALHFFVYDTLKIFTLLTVIIFVVSFIRSFFPLEKTREIL--SKHKAVALPLAAFLGILTPFCSCSAVPMFIGFVEAGIPLGAAF-TFLVASPMVNEVALGLLLTLFG--VKVAVLYVIFGVIVAIVAGYVIEKLNPRELIADYVFQV---------------------------KLGQTQIKEM-------------------TFKERLE---------FAKNNVKEILGKIW-------IYIIIAIGIGGFIHGYV-PQDIVERVA--KTAGLI-------AVPLAVLIGIPLYSNAAGILPVIQALIAKGVPLGTALAFMMATTALSFPEFMILKQIMKPKLIAFFAGIVGISIIAVGYLFNF",
        )
        self.assertEqual(
            alignment[2],
            "ILNGFYLLHEYAREHVLLCLVPAFFIAGTISVMLKKDAVLKLLGPNAKRIISYPVAAISGGILAVCSCTILPLFGGIYKKGAGIGPAT-TFLFAGPAINIAAIFLTARVLG--WDLGLARLIATITAAVLIGLIMEMIYQERGEGGLAFTSDD-----DQYGVRGIIFFLIQLG--FLVTSSLGINQTLKYS----------LMTLLGISALFM---ALFG--FKRDTVENWLYETWDFAKKILPYLFIGVFFAGVLTRLL-PQQVVTALL--GSNSFL-------SNLVASVIGTLMYFATLTEVPIVQALRELGMAKGPTLALLMAGNSLSLPSMIVITKLLGKKKAFTYFGLVVVFSTLFGMIYGV",
        )
        self.assertEqual(
            alignment[3],
            "LLAGIQALEEYIALHVLTCLVPAFLIAGALMSMMNKAVLINYLGAATSKLKSFPLAIVSSFFLAVCSCTVIPIASGIYKRTNATAPAM-IILWVAPATNILAVTYTGAVLG--LELALARIVAAISTAFVVGLILFYVFDRKIASQSDSAMPKAGRLVEN---NALVLFALLVAT-LLLPNYLGVGKPYIFKV--------EVFSVLMLVTTVY---ALKS--FSKEDLKYWMLETWFFVKQIIPLLLVGVFIVGVVGEILKATDVVEVYL--GGEGVG-------QSFLAALIGALSYFATMTEAPFVDTLMKLGMGKGPALALLLAGPGLSLPNMLAIGKLFGVKRAAVYIITIVALSTIAGVVYGE",
        )
        self.assertEqual(
            alignment[4],
            "---MINTIIDYLNVNRVLALLMAFLMAGGIASMINKNFIIKYFGSNTPKYISYTVAAVSGSLLAVCSCTILPLFASIYKRGAGIGPAT-TFLFSGPAINVLAIFYSAALLG--WDIGFLRAVFAVVVSILIGLSMEIIFKSHEKKRALR-VPKADKISDRPLYQTITFFALQFIMLLVITASPKLFPTLSMPLYDGFLLKHLLFIILGIILAVT---TKIW--FKDEEIKNWLRESFTLLKIVFPLLIIGVAIAGAIKAII-PPSYIATYV--GGNSIT-------ANFIASFIGALMYFATLTEVPIIKALMELGMGVGPAMALLLAGPSLSIPTVLTISKVLGKTKALTYLGLVVIFSTICGYIAGI",
        )
        self.assertEqual(
            alignment[5],
            "-IGHALALTASMTWEILWALILGFALSAVVQAVVRRSTIVTLLGDDRPR--TLVIATGLGAASSSCSYAAVALARSLFRKGANFTAAM-AFEIGSTNLVVELGIILALLMG--WQFTAAEFVGGPIMILVLAVLF-RLFVGARLIDAAREQAERGLAGSMEGHAAMDMS---------IKREGSFWRR------------------LLSPPGFT---S-----IAHVFVMEW-LAIL-------RDLILGLLIAGAIAAWV-PESFWQSFFLANHPAWSA----VWGPIIGPIVAIVSFVCSIGNVPLAAVLWNGGISFGGVIAF-IFADLLILPILNIYRKYYGARMMLVLLGTFYASMVVAGYLIE-",
        )
        self.assertEqual(
            alignment[6],
            "SAREALSTTAAMAWVTWWALVVGFAIAGGVEAWTSGEEVSELLEGHGPREIGY--GSLFGFVSSSCSYSAIATAKNLFKKGGSAAATLGAFMFASTNLVIEIGAVIWILLG--WQFLVADILGGFILIGLMAFGFVYLVPDEVVEQARRNVQDEGSETVRDPVCGMEVDPDETE--YSVERDGRTFYFCSKSCKESFDPEEANTTVRERATSLS---GWKA--LADKQWKEW-GMLW-------DEIAIGFVFAGLIAGFI-PDAVWTSVF--SGPTFGLPVYVFWTAVLGAVIGVATFVCSVGNVPFGAVLFSNGLPFGSVLSY-IYADLIVPPIVDAYREYYGTTFAAVLSGMIFVAAVLTGVVIHF",
        )
        self.assertEqual(
            alignment[7],
            "-FLQLNSIFISILIEAIPFILIGVILSGIIQMFVSEEMIARIM--PKNRFLAVLFGALAGVLFPACECGIIPITRRLLLKGVPLHAGV-AFMLTAPIINPIVLFSTYIAFGNRWSVVFYRGGLALAVSLIIGVILSYQFKDNQLLKPD------------------------------EPGHHHHHHG-------------------TLLQKLG---G-----TLRHAIDEF-FSVG-------KYLIIGAFIAAAMQTYV-KTSTLLAI---GQNDVS-------SSLVMMGLAFVLSLCSEVD-AFIASSFSSTFSLGSLIAFLVFGAMVDIKNLLMMLAAFKKRFVFLLITYIVVIVLAGSLLVKG",
        )
        self.assertEqual(
            alignment[8],
            "-WMNVNTIFLGIVIEAVPFILLGVFVSALIQIYVKEDTIQRYL--PKNAYAALLPAAVLGAIFPICECAIVPIVRRLIKKGMPLHVGV-VFLVAAPILNPIVAASTYFAFRTDLTVLYARMGLAFILSIVIGGLLYVMFKNSDQLKWTKEE---------------------------LVGRVPVQSD------------------MELKPKMN---RLKQ--TLYHASDEF-FLMG-------KYLIAGAFIAALFQTFL-DRNILVTI---GSNEWS-------STGVMMAFAFILSLCSEAD-AFVAASFGSTFTTGSLIAFLVYGPMLDLKNTIMLFAFFKSKFVLAFMITVTVVVFLAVMVLQF",
        )
        self.assertEqual(
            alignment[9],
            "QLHEAFTIFLSLLVEAIPFLTFGVVLSSALLVFSDEKKLIAYI--PRNPFLGAIAGSLVGFMFPVCECGNVPVARRFLMQGLPPSVAV-AFLLAAPTINPIVIWSTWVAFRDQPGMVVARVVCSLIITVIVSWVFSRQLDAVPLLKPALGRRLAYLTRPEESPTAIACESPLLQSGTFLLGSGNSGQLLKLD----------EQAVETLLPPIA---PSRWEMFTDNIVQEL-RELG-------GMLILGSLIAAVIQVFI-PREWILLL---GQGTIS-------SILAMMLLSVVVSVCSTVD-SFFALSFVSTFTSSSLLAFLVFGPMIDVKSIGLLLSVFQRRIVIYLLLLTGQLTFLLSLAHSY",
        )
        self.assertEqual(
            alignment[10],
            "EFNLFLDLLGSALLLSLPWLLLGIIISSTFLIWTDEQKWVANF--PRNRLLSSLVGSALGFLLPLGAFGSVPLVRRLLLQGAPIPLAV-SFLVAAPTLNIFAIVRVLSSRQSQYGLIFLCISCSWLMAIVMGLVFSTYRLARQQAEDEGETALLNIPLLRSGALIILQSSMEA-----SPRQGGLVFA------------------SGVNPVADFSWRQKLHLFGRNIIEEF-QEFG-------GVLVIGTAIACGIVFFL-PQAWLLQWA--GLGPVR-------QTVLMMGWSFILPLGNFSN-PDLLAPLGEQLWRGSMVAFLLWGSLFNLQTIGLWLVTLRLRPLSYLVVLVGLSVFLFAMVTNY",
        )
        self.assertEqual(
            alignment.column_annotations["consensus sequence"],
            ".htphhslhhhhslcslhhLlhuhhluuslpsahscpplhchL..s+s+hluhhlAulhGhlhssCSCuslPlspslhc+GsslusAh.sFLluuPslN.lslhhshhlhG..aplshhcllsuhllulllGllhthlh.spthtcsshph.............hh............lhspssltps...................tlhssls...s.h...hs+ptlcEa.hchh.......shLlIGshIAGsIpsal.Ppshlhshh..Gsssls.......ushluslluhlhahsohsshPhlsuLhspGhshGoslAaLlhGshLslPshhhltphhtt+hshshlshlslhshlsGhlhsh",
        )
        self.assertEqual(
            alignment.column_annotations["reference coordinate annotation"],
            "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx.xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx.xxxxxxxxxxxxxxxxxxxxxx..xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx................................................xxxxxxxx...xxxx..xxxxxxxxx.xxxx.......xxxxxxxxxxxxxxxxx.xxxxxxxxx..xxxxxx.......xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
[[  0,   1,   3,  43,  44,  44,  48,  50,  52,  54,  87,  87, 109, 109, 131, 132, 144, 145, 146, 147, 149, 149, 149, 149, 149, 149, 149, 149, 149, 149, 159, 159, 159, 159, 159, 159, 160, 167, 167, 168, 171, 171, 180, 181, 185, 185, 202, 202, 210, 211, 211, 217, 217, 217, 217, 236, 237, 256, 257, 299, 300],
 [  0,   1,   3,  43,  43,  43,  47,  49,  51,  53,  86,  86, 108, 108, 130, 131, 143, 144, 145, 146, 146, 146, 146, 146, 146, 146, 146, 146, 146, 146, 156, 156, 156, 156, 156, 156, 156, 163, 163, 163, 163, 163, 172, 173, 177, 177, 194, 194, 202, 203, 203, 209, 209, 209, 209, 228, 229, 248, 249, 291, 292],
 [  0,   1,   3,  43,  44,  45,  49,  51,  53,  55,  88,  88, 110, 110, 132, 133, 145, 146, 147, 148, 150, 150, 152, 155, 161, 165, 166, 166, 166, 168, 178, 182, 182, 182, 182, 186, 187, 194, 194, 195, 198, 198, 207, 208, 212, 219, 236, 236, 244, 245, 245, 251, 251, 251, 251, 270, 271, 290, 291, 333, 334],
 [  0,   1,   3,  43,  44,  45,  49,  51,  53,  55,  88,  88, 110, 110, 132, 133, 145, 146, 147, 148, 150, 155, 157, 157, 163, 167, 168, 169, 169, 171, 181, 185, 186, 186, 187, 191, 192, 199, 199, 200, 203, 203, 212, 213, 217, 224, 241, 242, 250, 251, 251, 257, 257, 257, 257, 276, 277, 296, 297, 339, 340],
 [  0,   0,   0,  40,  41,  42,  46,  48,  50,  52,  85,  85, 107, 107, 129, 130, 142, 143, 143, 144, 146, 151, 153, 156, 162, 166, 167, 168, 169, 171, 181, 185, 186, 194, 195, 199, 200, 207, 207, 208, 211, 211, 220, 221, 225, 232, 249, 249, 257, 258, 258, 264, 264, 264, 264, 283, 284, 303, 304, 346, 347],
 [  0,   0,   2,  42,  43,  44,  48,  48,  50,  52,  85,  85, 107, 107, 129, 129, 141, 142, 143, 144, 146, 151, 153, 156, 162, 162, 162, 162, 162, 162, 172, 172, 172, 172, 172, 172, 173, 180, 180, 181, 181, 181, 190, 190, 194, 194, 211, 211, 219, 220, 222, 228, 229, 229, 231, 250, 251, 270, 270, 312, 312],
 [  0,   1,   3,  43,  44,  45,  49,  51,  53,  53,  86,  87, 109, 109, 131, 132, 144, 145, 146, 147, 149, 154, 156, 159, 165, 169, 170, 170, 170, 172, 182, 186, 187, 195, 196, 200, 201, 208, 208, 209, 212, 212, 221, 221, 225, 225, 242, 242, 250, 251, 251, 257, 258, 262, 264, 283, 284, 303, 303, 345, 346],
 [  0,   0,   2,  42,  42,  42,  46,  48,  50,  52,  85,  85, 107, 109, 131, 132, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 154, 154, 154, 154, 154, 154, 154, 161, 161, 162, 162, 162, 171, 171, 175, 175, 192, 192, 200, 200, 200, 206, 206, 206, 206, 225, 225, 244, 245, 287, 288],
 [  0,   0,   2,  42,  42,  42,  46,  48,  50,  52,  85,  85, 107, 109, 131, 132, 144, 145, 146, 147, 147, 147, 147, 147, 147, 147, 147, 147, 147, 147, 157, 157, 157, 157, 157, 157, 158, 165, 165, 166, 169, 169, 178, 178, 182, 182, 199, 199, 207, 207, 207, 213, 213, 213, 213, 232, 232, 251, 252, 294, 295],
 [  0,   1,   3,  43,  43,  43,  47,  49,  51,  53,  86,  86, 108, 110, 132, 133, 145, 146, 147, 148, 150, 155, 157, 160, 166, 170, 171, 172, 173, 175, 185, 189, 189, 189, 189, 193, 194, 201, 201, 202, 205, 207, 216, 216, 220, 220, 237, 237, 245, 245, 245, 251, 251, 251, 251, 270, 270, 289, 290, 332, 333],
 [  0,   1,   3,  43,  43,  43,  47,  49,  51,  53,  86,  86, 108, 110, 132, 133, 145, 146, 147, 148, 150, 155, 157, 160, 166, 170, 170, 170, 170, 170, 180, 180, 180, 180, 180, 180, 181, 188, 191, 192, 195, 197, 206, 206, 210, 210, 227, 227, 235, 236, 236, 242, 242, 242, 242, 261, 261, 280, 281, 323, 324]
]
                    # fmt: on
                ),
            )
        )

    def check_alignment_pfam6(self, alignment):
        """Check the alignment obtained by parsing Pfam record COX2_TM."""
        self.assertEqual(alignment.annotations["identifier"], "COX2_TM")
        self.assertEqual(alignment.annotations["accession"], "PF02790.17")
        self.assertEqual(
            alignment.annotations["definition"],
            "Cytochrome C oxidase subunit II, transmembrane domain",
        )
        self.assertEqual(len(alignment.annotations["author"]), 2)
        self.assertEqual(
            alignment.annotations["author"][0], "Sonnhammer ELL;0000-0002-9015-5588"
        )
        self.assertEqual(
            alignment.annotations["author"][1], "Griffiths-Jones SR;0000-0001-6043-807X"
        )
        self.assertEqual(alignment.annotations["source of seed"], "Prosite")
        self.assertEqual(alignment.annotations["gathering method"], "22.80 18.00;")
        self.assertEqual(alignment.annotations["trusted cutoff"], "22.80 21.40;")
        self.assertEqual(alignment.annotations["noise cutoff"], "22.70 17.90;")
        self.assertEqual(
            alignment.annotations["build method"], "hmmbuild HMM.ann SEED.ann"
        )
        self.assertEqual(
            alignment.annotations["search method"],
            "hmmsearch -Z 57096847 -E 1000 --cpu 4 HMM pfamseq",
        )
        self.assertEqual(alignment.annotations["type"], "Family")
        self.assertEqual(
            alignment.annotations["wikipedia"], ["Cytochrome_c_oxidase_subunit_II"]
        )
        self.assertEqual(len(alignment.annotations["references"]), 1)
        self.assertEqual(alignment.annotations["references"][0]["number"], 1)
        self.assertEqual(alignment.annotations["references"][0]["medline"], "8638158")
        self.assertEqual(
            alignment.annotations["references"][0]["title"],
            "The whole structure of the 13-subunit oxidized cytochrome c oxidase at 2.8 A.",
        )
        self.assertEqual(
            alignment.annotations["references"][0]["author"],
            "Tsukihara T, Aoyama H, Yamashita E, Tomizaki T, Yamaguchi H, Shinzawa-Itoh K, Nakashima R, Yaono R, Yoshikawa S;",
        )
        self.assertEqual(
            alignment.annotations["references"][0]["location"],
            "Science 1996;272:1136-1144.",
        )
        self.assertEqual(len(alignment.annotations["database references"]), 5)
        self.assertEqual(
            alignment.annotations["database references"][0]["reference"],
            "INTERPRO; IPR011759;",
        )
        self.assertEqual(
            alignment.annotations["database references"][1]["reference"],
            "PROSITE; PDOC00075;",
        )
        self.assertEqual(
            alignment.annotations["database references"][2]["reference"],
            "SCOP; 1occ; fa;",
        )
        self.assertEqual(
            alignment.annotations["database references"][2]["comment"],
            "This family corresponds to chains b and o.",
        )
        self.assertEqual(
            alignment.annotations["database references"][3]["reference"], "TC; 3.D.4;"
        )
        self.assertEqual(
            alignment.annotations["database references"][4]["reference"],
            "SO; 0100021; polypeptide_conserved_region;",
        )
        self.assertEqual(
            alignment.annotations["comment"],
            "The N-terminal domain of cytochrome C oxidase contains two transmembrane alpha-helices.",
        )
        self.assertEqual(len(alignment.sequences), 11)
        self.assertEqual(alignment.sequences[0].id, "COX2_SCHPO/11-99")
        self.assertEqual(alignment.sequences[1].id, "COX2_CANGA/17-103")
        self.assertEqual(alignment.sequences[2].id, "COX2_NEUCR/14-102")
        self.assertEqual(alignment.sequences[3].id, "H9D0Q0_EMENI/15-102")
        self.assertEqual(alignment.sequences[4].id, "COX2_ARATH/17-103")
        self.assertEqual(alignment.sequences[5].id, "COX2_ANOGA/1-83")
        self.assertEqual(alignment.sequences[6].id, "COX2_CHICK/1-82")
        self.assertEqual(alignment.sequences[7].id, "COX2_SHEEP/1-83")
        self.assertEqual(alignment.sequences[8].id, "COX2_STRPU/1-83")
        self.assertEqual(alignment.sequences[9].id, "COX2_SYNY3/19-111")
        self.assertEqual(alignment.sequences[10].id, "A1BA41_PARDP/42-128")
        self.assertEqual(alignment.sequences[0].annotations["accession"], "P21534.4")
        self.assertEqual(alignment.sequences[1].annotations["accession"], "P43373.2")
        self.assertEqual(alignment.sequences[2].annotations["accession"], "P00411.2")
        self.assertEqual(alignment.sequences[3].annotations["accession"], "H9D0Q0.1")
        self.assertEqual(alignment.sequences[4].annotations["accession"], "P93285.2")
        self.assertEqual(alignment.sequences[5].annotations["accession"], "P34840.1")
        self.assertEqual(alignment.sequences[6].annotations["accession"], "P18944.1")
        self.assertEqual(alignment.sequences[7].annotations["accession"], "O78750.1")
        self.assertEqual(alignment.sequences[8].annotations["accession"], "P15545.1")
        self.assertEqual(alignment.sequences[9].annotations["accession"], "Q06474.2")
        self.assertEqual(alignment.sequences[10].annotations["accession"], "A1BA41.1")
        self.assertEqual(
            alignment.sequences[0].seq,
            "APSSWALYFQDGASPSYLGVTHLNDYLMFYLTFIFIGVIYAICKAVIEYNYNSHPIAAKYTTHGSIVEFIWTLIPALILILVALPSFKL",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "VPTPYGIYFQDSATPNQEGILELHDNIMFYLFIILGLVSWMLFTIVKTYSKNPMAYKYIKHGQTIEIIWTMFPAVILLIIAFPSFIL",
        )
        self.assertEqual(
            alignment.sequences[2].seq,
            "APSPWGIYFQDSATPQMEGLVELHDNIMYYLVVILFGVGWILLSIIRNYISTKSPISHKYLNHGTLIELIWTITPAVILILIAFPSFKL",
        )
        self.assertEqual(
            alignment.sequences[3].seq,
            "PTPWGIFFQDSASPQMEGIEELHNNIMFYLAIILFTVTWMMITIIRNFVAKKSPIAHKYMNHGTLIELIWTITPAFILILIAFPSFKL",
        )
        self.assertEqual(
            alignment.sequences[4].seq,
            "AEPWQLGFQDAATPIMQGIIDLHHDIFFFLILILVFVLWILVRALWHFHYKKNAIPQRIVHGTTIEILWTIFPSIILMFIAIPSFAL",
        )
        self.assertEqual(
            alignment.sequences[5].seq,
            "MATWANLGLQDSSSPLMEQLNFFHDHTLLILTMITILVGYIMGMLSFNKFTNRFLLHGQTIEIIWTVLPAIILMFIAFPSLRL",
        )
        self.assertEqual(
            alignment.sequences[6].seq,
            "MANHSQLGFQDASSPIMEELVEFHDHALMVALAICSLVLYLLTLMLMEKLSSNTVDAQEVELIWTILPAIVLVLLALPSLQI",
        )
        self.assertEqual(
            alignment.sequences[7].seq,
            "MAYPMQLGFQDATSPIMEELLHFHDHTLMIVFLISSLVLYIISLMLTTKLTHTSTMDAQEVETIWTILPAIILIMIALPSLRI",
        )
        self.assertEqual(
            alignment.sequences[8].seq,
            "MGTWAQFGLQDASSPLMEELTYFHDYALIVLTLITILVFYGLVSLLVSSNTNRFFFEGQELETIWTVIPALILILIALPSLQL",
        )
        self.assertEqual(
            alignment.sequences[9].seq,
            "VSLWYGQNHGLMPVAASADAEKVDGIFNYMMTIATGLFLLVEGVLVYCLIRFRRRKDDQTDGPPIEGNVPLEILWTAIPTVIVFTLAVYSFEV",
        )
        self.assertEqual(
            alignment.sequences[10].seq,
            "PVNGGMNFQPASSPLAHDQQWLDHFVLYIITAVTIFVCLLLLICIVRFNRRANPVPARFTHNTPIEVIWTLVPVLILVAIGAFSLPI",
        )
        self.assertEqual(
            alignment[0],
            "APSSWALY---FQDGASPSYLGVTHLNDYLMFYLTFIFIGVIYAICKAVIEYNYNSHPIAAKYTTHGSI-VEFIWTLIPALILILVALPSFKL",
        )
        self.assertEqual(
            alignment[1],
            "VPTPYGIY---FQDSATPNQEGILELHDNIMFYLFIILGLVSWMLFTIVKTY--SKNPMAYKYIKHGQT-IEIIWTMFPAVILLIIAFPSFIL",
        )
        self.assertEqual(
            alignment[2],
            "APSPWGIY---FQDSATPQMEGLVELHDNIMYYLVVILFGVGWILLSIIRNYISTKSPISHKYLNHGTL-IELIWTITPAVILILIAFPSFKL",
        )
        self.assertEqual(
            alignment[3],
            "-PTPWGIF---FQDSASPQMEGIEELHNNIMFYLAIILFTVTWMMITIIRNFVAKKSPIAHKYMNHGTL-IELIWTITPAFILILIAFPSFKL",
        )
        self.assertEqual(
            alignment[4],
            "-AEPWQLG---FQDAATPIMQGIIDLHHDIFFFLILILVFVLWILVRALWHFHYKKNAIPQR-IVHGTT-IEILWTIFPSIILMFIAIPSFAL",
        )
        self.assertEqual(
            alignment[5],
            "MATWANLG---LQDSSSPLMEQLNFFHDHTLLILTMITILVGYIMGMLSFN------KFTNRFLLHGQT-IEIIWTVLPAIILMFIAFPSLRL",
        )
        self.assertEqual(
            alignment[6],
            "MANHSQLG---FQDASSPIMEELVEFHDHALMVALAICSLVLYLLTLMLME------KLS-SNTVDAQE-VELIWTILPAIVLVLLALPSLQI",
        )
        self.assertEqual(
            alignment[7],
            "MAYPMQLG---FQDATSPIMEELLHFHDHTLMIVFLISSLVLYIISLMLTT------KLTHTSTMDAQE-VETIWTILPAIILIMIALPSLRI",
        )
        self.assertEqual(
            alignment[8],
            "MGTWAQFG---LQDASSPLMEELTYFHDYALIVLTLITILVFYGLVSLLVS------SNTNRFFFEGQE-LETIWTVIPALILILIALPSLQL",
        )
        self.assertEqual(
            alignment[9],
            "VSLWYGQNHGLMPVAASADAEKVDGIFNYMMTIATGLFLLVEGVLVYCLIRFRRRKDDQTDGPPIEGNVPLEILWTAIPTVIVFTLAVYSFEV",
        )
        self.assertEqual(
            alignment[10],
            "-PVNGGMN---FQPASSPLAHDQQWLDHFVLYIITAVTIFVCLLLLICIVRFNRRANPVPAR-FTHNTP-IEVIWTLVPVLILVAIGAFSLPI",
        )
        self.assertEqual(
            alignment.column_annotations["consensus sequence"],
            "hssshsls...FQDuuSP.MEtlhclHDahhhhLshIhhhVhalLshhlhpa..ptpslsp+.hhHGph.lElIWTllPAlILlhIAhPShpL",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [
                        [0, 1, 8, 8, 48, 49, 51, 54, 57, 58, 59, 60, 66, 66, 89],
                        [0, 1, 8, 8, 48, 49, 49, 52, 55, 56, 57, 58, 64, 64, 87],
                        [0, 1, 8, 8, 48, 49, 51, 54, 57, 58, 59, 60, 66, 66, 89],
                        [0, 0, 7, 7, 47, 48, 50, 53, 56, 57, 58, 59, 65, 65, 88],
                        [0, 0, 7, 7, 47, 48, 50, 53, 56, 57, 58, 58, 64, 64, 87],
                        [0, 1, 8, 8, 48, 48, 48, 48, 51, 52, 53, 54, 60, 60, 83],
                        [0, 1, 8, 8, 48, 48, 48, 48, 51, 51, 52, 53, 59, 59, 82],
                        [0, 1, 8, 8, 48, 48, 48, 48, 51, 52, 53, 54, 60, 60, 83],
                        [0, 1, 8, 8, 48, 48, 48, 48, 51, 52, 53, 54, 60, 60, 83],
                        [0, 1, 8, 11, 51, 52, 54, 57, 60, 61, 62, 63, 69, 70, 93],
                        [0, 0, 7, 7, 47, 48, 50, 53, 56, 57, 58, 58, 64, 64, 87],
                    ],
                ),
            )
        )

    def check_alignment_pfam7(self, alignment):
        """Check the alignment obtained by parsing Pfam record Alpha_E1_glycop."""
        self.assertEqual(alignment.annotations["identifier"], "Alpha_E1_glycop")
        self.assertEqual(alignment.annotations["accession"], "PF01589.18")
        self.assertEqual(
            alignment.annotations["definition"], "Alphavirus E1 glycoprotein"
        )
        self.assertEqual(
            alignment.annotations["author"], ["Bateman A;0000-0002-6982-4660"]
        )
        self.assertEqual(
            alignment.annotations["source of seed"], "Pfam-B_587 (release 4.1)"
        )
        self.assertEqual(alignment.annotations["gathering method"], "25.00 25.00;")
        self.assertEqual(alignment.annotations["trusted cutoff"], "34.50 33.30;")
        self.assertEqual(alignment.annotations["noise cutoff"], "24.10 23.50;")
        self.assertEqual(
            alignment.annotations["build method"], "hmmbuild HMM.ann SEED.ann"
        )
        self.assertEqual(
            alignment.annotations["search method"],
            "hmmsearch -Z 57096847 -E 1000 --cpu 4 HMM pfamseq",
        )
        self.assertEqual(alignment.annotations["type"], "Family")
        self.assertEqual(alignment.annotations["wikipedia"], ["Alphavirus"])
        self.assertEqual(alignment.annotations["clan"], "CL0543")
        self.assertEqual(len(alignment.annotations["references"]), 2)
        self.assertEqual(
            alignment.annotations["references"][0]["comment"],
            "This paper includes cryoelectron microscopy images.",
        )
        self.assertEqual(alignment.annotations["references"][0]["number"], 1)
        self.assertEqual(alignment.annotations["references"][0]["medline"], "7867069")
        self.assertEqual(
            alignment.annotations["references"][0]["title"],
            "Nucleocapsid and glycoprotein organization in an enveloped virus.",
        )
        self.assertEqual(
            alignment.annotations["references"][0]["author"],
            "Cheng RH, Kuhn RJ, Olson NH, Rossmann MG, Choi HK, Smith TJ, Baker TS;",
        )
        self.assertEqual(
            alignment.annotations["references"][0]["location"], "Cell 1995;80:621-630."
        )
        self.assertEqual(alignment.annotations["references"][1]["number"], 2)
        self.assertEqual(alignment.annotations["references"][1]["medline"], "8995682")
        self.assertEqual(
            alignment.annotations["references"][1]["title"],
            "Role of glycoprotein PE2 in formation and maturation of the Sindbis virus spike.",
        )
        self.assertEqual(
            alignment.annotations["references"][1]["author"],
            "Carleton M, Lee H, Mulvey M, Brown DT;",
        )
        self.assertEqual(
            alignment.annotations["references"][1]["location"],
            "J Virol 1997;71:1558-1566.",
        )
        self.assertEqual(len(alignment.annotations["database references"]), 4)
        self.assertEqual(
            alignment.annotations["database references"][0],
            {"reference": "INTERPRO; IPR002548;"},
        )
        self.assertEqual(
            alignment.annotations["database references"][1],
            {"reference": "TC; 1.A.34;"},
        )
        self.assertEqual(
            alignment.annotations["database references"][2],
            {"reference": "SCOP; 1rer; fa;"},
        )
        self.assertEqual(
            alignment.annotations["database references"][3],
            {"reference": "SO; 0100021; polypeptide_conserved_region;"},
        )
        self.assertEqual(
            alignment.annotations["comment"],
            "E1 forms a heterodimer with E2 Pfam:PF00943.  The virus spikes are made up of 80 trimers of these heterodimers (sindbis virus) [2].",
        )
        self.assertEqual(len(alignment.sequences), 2)
        self.assertEqual(alignment.sequences[0].id, "POLS_SFV/751-1253")
        self.assertEqual(alignment.sequences[1].id, "POLS_CHIKS/744-1247")
        self.assertEqual(alignment.sequences[0].annotations["accession"], "P03315.1")
        self.assertEqual(len(alignment.sequences[0].dbxrefs), 6)
        self.assertEqual(alignment.sequences[0].dbxrefs[0], "PDB; 2V33 A; 292-382;")
        self.assertEqual(alignment.sequences[0].dbxrefs[1], "PDB; 2ALA A; 1-384;")
        self.assertEqual(alignment.sequences[0].dbxrefs[2], "PDB; 1RER A; 1-391;")
        self.assertEqual(alignment.sequences[0].dbxrefs[3], "PDB; 2V33 B; 292-382;")
        self.assertEqual(alignment.sequences[0].dbxrefs[4], "PDB; 1RER C; 1-391;")
        self.assertEqual(alignment.sequences[0].dbxrefs[5], "PDB; 1RER B; 1-391;")
        self.assertEqual(alignment.sequences[1].annotations["accession"], "Q8JUX5.3")
        self.assertEqual(len(alignment.sequences[1].dbxrefs), 1)
        self.assertEqual(alignment.sequences[1].dbxrefs[0], "PDB; 2RSW A; 1-18;")
        self.assertEqual(
            alignment.sequences[0].seq,
            "PRAHAASVAETMAYLWDQNQALFWLEFAAPVACILIITYCLRNVLCCCKSLSFLVLLSLGATARAYEHSTVMPNVVGFPYKAHIERPGYSPLTLQMQVVETSLEPTLNLEYITCEYKTVVPSPYVKCCGASECSTKEKPDYQCKVYTGVYPFMWGGAYCFCDSENTQLSEAYVDRSDVCRHDHASAYKAHTASLKAKVRVMYGNVNQTVDVYVNGDHAVTIGGTQFIFGPLSSAWTPFDNKIVVYKDEVFNQDFPPYGSGQPGRFGDIQSRTVESNDLYANTALKLARPSPGMVHVPYTQTPSGFKYWLKEKGTALNTKAPFGCQIKTNPVRAMNCAVGNIPVSMNLPDSAFTRIVEAPTIIDLTCTVATCTHSSDFGGVLTLTYKTNKNGDCSVHSHSNVATLQEATAKVKTAGKVTLHFSTASASPSFVVSLCSARATCSASCEPPKDHIVPYAASHSNVVFPDMSGTALSWVQKISGGLGAFAIGAILVLVVVTCIGLRR",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "RTAKAATYQEAAVYLWNEQQPLFWLQALIPLAALIVLCNCLRLLPCCCKTLAFLAVMSIGAHTVSAYEHVTVIPNTVGVPYKTLVNRPGYSPMVLEMELLSVTLEPTLSLDYITCEYKTVIPSPYVKCCGTAECKDKNLPDYSCKVFTGVYPFMWGGAYCFCDAENTQLSEAHVEKSESCKTEFASAYRAHTASASAKLRVLYQGNNITVTAYANGDHAVTVKDAKFIVGPMSSAWTPFDNKIVVYKGDVYNMDYPPFGAGRPGQFGDIQSRTPESKDVYANTQLVLQRPAAGTVHVPYSQAPSGFKYWLKERGASLQHTAPFGCQIATNPVRAMNCAVGNMPISIDIPDAAFTRVVDAPSLTDMSCEVPACTHSSDFGGVAIIKYAVSKKGKCAVHSMTNAVTIREAEIEVEGNSQLQISFSTALASAEFRVQVCSTQVHCAAECHPPKDHIVNYPASHTTLGVQDISATAMSWVQKITGGVGLVVAVAALILIVVLCVSFSR",
        )
        self.assertEqual(
            alignment.sequences[0].letter_annotations["secondary structure"],
            "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX-EEEEEEES-SS--EEEEE--TTSS-EEEEEEEEEEEEEEEEEEEEEE--EEEE----EEETSS--------STT-EEEEEES-SSSISCSCS-SSSSS-EEEEEEEEEE-TTGGGS-EEEEEEEEEEEEEEEEEEETTEEEEEEEESSSS-EEEETTEEEEE---S----S--SSEEE-SS-EEE-----GGG--TTSTTSEEBSSTT-S--EE-S--EE----SSSS---EE----HHHHHHHHS-S-GGGT-STT-EEETTTTEEES---SEEEEEESS-TTTS-EETTS--EEEEEEEEEEEBTTSTTEEEEEEEEEESS-EEEEEEESSTTEEESBSEEEE-TT-EEEEEEEESSSS-EEEEEETTEEEEEE---B----------XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
        )
        self.assertEqual(
            alignment.sequences[1].letter_annotations["secondary structure"],
            "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX--HHHHH---STTTTS--XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
        )
        self.assertEqual(
            alignment[0],
            "PRAHAASVAETMAYLWDQNQALFWLEFAAPVACILIITYCLRNVLCCCKSLSFLVLLSLG-ATARAYEHSTVMPNVVGFPYKAHIERPGYSPLTLQMQVVETSLEPTLNLEYITCEYKTVVPSPYVKCCGASECSTKEKPDYQCKVYTGVYPFMWGGAYCFCDSENTQLSEAYVDRSDVCRHDHASAYKAHTASLKAKVRVMYGNVNQTVDVYVNGDHAVTIGGTQFIFGPLSSAWTPFDNKIVVYKDEVFNQDFPPYGSGQPGRFGDIQSRTVESNDLYANTALKLARPSPGMVHVPYTQTPSGFKYWLKEKGTALNTKAPFGCQIKTNPVRAMNCAVGNIPVSMNLPDSAFTRIVEAPTIIDLTCTVATCTHSSDFGGVLTLTYKTNKNGDCSVHSHSNVATLQEATAKVKTAGKVTLHFSTASASPSFVVSLCSARATCSASCEPPKDHIVPYAASHSNVVFPDMSGTALSWVQKISGGLGAFAIGAILVLVVVTCIGLRR",
        )
        self.assertEqual(
            alignment[1],
            "RTAKAATYQEAAVYLWNEQQPLFWLQALIPLAALIVLCNCLRLLPCCCKTLAFLAVMSIGAHTVSAYEHVTVIPNTVGVPYKTLVNRPGYSPMVLEMELLSVTLEPTLSLDYITCEYKTVIPSPYVKCCGTAECKDKNLPDYSCKVFTGVYPFMWGGAYCFCDAENTQLSEAHVEKSESCKTEFASAYRAHTASASAKLRVLYQGNNITVTAYANGDHAVTVKDAKFIVGPMSSAWTPFDNKIVVYKGDVYNMDYPPFGAGRPGQFGDIQSRTPESKDVYANTQLVLQRPAAGTVHVPYSQAPSGFKYWLKERGASLQHTAPFGCQIATNPVRAMNCAVGNMPISIDIPDAAFTRVVDAPSLTDMSCEVPACTHSSDFGGVAIIKYAVSKKGKCAVHSMTNAVTIREAEIEVEGNSQLQISFSTALASAEFRVQVCSTQVHCAAECHPPKDHIVNYPASHTTLGVQDISATAMSWVQKITGGVGLVVAVAALILIVVLCVSFSR",
        )
        self.assertEqual(
            alignment.column_annotations["consensus secondary structure"],
            "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXEEEEEEESXSSXXEEEEEXXTTSSXEEEEEEEEEEEEEEEEEEEEEEXXEEEEXXXXEEETSSXXXXXXXXSTTXEEEEEES-SCCHCCSCSSTTTTS-EEEEEEEEEEXTTGGGSXEEEEEEEEEEEEEEEEEEETTEEEEEEEESSSSXEEEETTEEEEEXXXSXXXXSXXSSEEEXSSXEEEXXXXXGGGXXTTSTTSEEBSSTTXSXXEEXSXXEEXXXXSSSSXXXEEXXXXHHHHHHHHSXSXGGGTXSTTXEEETTTTEEESXXXSEEEEEESSXTTTSXEETTSXXEEEEEEEEEEEBTTSTTEEEEEEEEEESSXEEEEEEESSTTEEESBSEEEEXTTXEEEEEEEESSSSXEEEEEETTEEEEEEXXXBXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
        )
        self.assertEqual(
            alignment.column_annotations["consensus sequence"],
            ".pA+AAohtEshsYLWsppQsLFWLphhhPlAslllls.CLR.l.CCCKoLuFLslhSlG.tTspAYEHsTVhPNsVGhPYKshlpRPGYSPhsLpMpllpsoLEPTLsL-YITCEYKTVlPSPYVKCCGsuECpsKphPDYpCKVaTGVYPFMWGGAYCFCDuENTQLSEAaV-+S-sC+p-aASAY+AHTAShpAKlRVhYtssN.TVssYsNGDHAVTltsspFIhGPhSSAWTPFDNKIVVYKs-VaN.DaPPaGuGpPGpFGDIQSRTsESpDlYANTtLhLtRPusGhVHVPYoQsPSGFKYWLKE+GsuLpppAPFGCQItTNPVRAMNCAVGNhPlShslPDuAFTRlV-APolhDhoCpVssCTHSSDFGGVhhlpYtssKpGcCuVHShoNssTlpEAphcVcssuplplpFSTA.ASspFhVplCSspspCuApCcPPKDHIVsYsASHoslsh.DhSuTAhSWVQKIoGGlGhhshsAhLlLlVVhCluhpR",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[0, 60, 60, 503], [0, 60, 61, 504]])
            )
        )

    def check_alignment_pfam8(self, alignment):
        """Check the alignment obtained by parsing Pfam record Cyclin_N."""
        self.assertEqual(alignment.annotations["identifier"], "Cyclin_N")
        self.assertEqual(alignment.annotations["accession"], "PF00134.25")
        self.assertEqual(
            alignment.annotations["definition"], "Cyclin, N-terminal domain"
        )
        self.assertEqual(alignment.annotations["previous identifier"], "cyclin;")
        self.assertEqual(
            alignment.annotations["author"],
            [
                "Bateman A;0000-0002-6982-4660",
                "Sonnhammer ELL;0000-0002-9015-5588",
                "Griffiths-Jones SR;0000-0001-6043-807X",
            ],
        )
        self.assertEqual(alignment.annotations["source of seed"], "Prosite")
        self.assertEqual(alignment.annotations["gathering method"], "20.50 20.50;")
        self.assertEqual(alignment.annotations["trusted cutoff"], "20.50 20.50;")
        self.assertEqual(alignment.annotations["noise cutoff"], "20.40 20.40;")
        self.assertEqual(
            alignment.annotations["build method"], "hmmbuild HMM.ann SEED.ann"
        )
        self.assertEqual(
            alignment.annotations["search method"],
            "hmmsearch -Z 57096847 -E 1000 --cpu 4 HMM pfamseq",
        )
        self.assertEqual(alignment.annotations["type"], "Domain")
        self.assertEqual(alignment.annotations["wikipedia"], ["Cyclin"])
        self.assertEqual(alignment.annotations["clan"], "CL0065")
        self.assertEqual(
            alignment.annotations["references"][0]["comment"],
            "The cyclins include an internal duplication, which is related to that found in TFIIB and the RB protein.",
        )
        self.assertEqual(alignment.annotations["references"][0]["number"], 1)
        self.assertEqual(alignment.annotations["references"][0]["medline"], "8152925")
        self.assertEqual(
            alignment.annotations["references"][0]["title"],
            "Evidence for a protein domain superfamily shared by the cyclins, TFIIB and RB/p107.",
        )
        self.assertEqual(
            alignment.annotations["references"][0]["author"],
            "Gibson TJ, Thompson JD, Blocker A, Kouzarides T;",
        )
        self.assertEqual(
            alignment.annotations["references"][0]["location"],
            "Nucleic Acids Res 1994;22:946-952.",
        )
        self.assertEqual(alignment.annotations["references"][1]["number"], 2)
        self.assertEqual(alignment.annotations["references"][1]["medline"], "8591034")
        self.assertEqual(
            alignment.annotations["references"][1]["title"],
            "The crystal structure of cyclin A",
        )
        self.assertEqual(
            alignment.annotations["references"][1]["author"],
            "Brown NR, Noble MEM, Endicott JA, Garman EF, Wakatsuki S, Mitchell E, Rasmussen B, Hunt T, Johnson LN;",
        )
        self.assertEqual(
            alignment.annotations["references"][1]["location"],
            "Structure. 1995;3:1235-1247.",
        )
        self.assertEqual(
            alignment.annotations["references"][2]["comment"],
            "Complex of cyclin and cyclin dependent kinase.",
        )
        self.assertEqual(alignment.annotations["references"][2]["number"], 3)
        self.assertEqual(alignment.annotations["references"][2]["medline"], "8756328")
        self.assertEqual(
            alignment.annotations["references"][2]["title"],
            "Structural basis of cyclin-dependant kinase activation by phosphorylation.",
        )
        self.assertEqual(
            alignment.annotations["references"][2]["author"],
            "Russo AA, Jeffrey PD, Pavletich NP;",
        )
        self.assertEqual(
            alignment.annotations["references"][2]["location"],
            "Nat Struct Biol. 1996;3:696-700.",
        )
        self.assertEqual(alignment.annotations["references"][3]["number"], 4)
        self.assertEqual(alignment.annotations["references"][3]["medline"], "2001396")
        self.assertEqual(
            alignment.annotations["references"][3]["title"],
            "Isolation and characterization of a human cDNA encoding uracil-DNA glycosylase.",
        )
        self.assertEqual(
            alignment.annotations["references"][3]["author"], "Muller SJ, Caradonna S;"
        )
        self.assertEqual(
            alignment.annotations["references"][3]["location"],
            "Biochim Biophys Acta 1991;1088:197-207.",
        )
        self.assertEqual(len(alignment.annotations["database references"]), 5)
        self.assertEqual(
            alignment.annotations["database references"][0],
            {"reference": "INTERPRO; IPR006671;"},
        )
        self.assertEqual(
            alignment.annotations["database references"][1],
            {"reference": "PROSITE; PDOC00264;"},
        )
        self.assertEqual(
            alignment.annotations["database references"][2],
            {"reference": "SCOP; 1vin; fa;"},
        )
        self.assertEqual(
            alignment.annotations["database references"][3],
            {"reference": "HOMSTRAD; cyclin;"},
        )
        self.assertEqual(
            alignment.annotations["database references"][4],
            {"reference": "SO; 0000417; polypeptide_domain;"},
        )
        self.assertEqual(
            alignment.annotations["comment"],
            "Cyclins regulate cyclin dependent kinases (CDKs). Swiss:P22674 is a Uracil-DNA glycosylase that is related to other cyclins [4]. Cyclins contain two domains of similar all-alpha fold, of which this family corresponds with the N-terminal domain.",
        )
        self.assertEqual(len(alignment.sequences), 95)
        self.assertEqual(alignment.sequences[0].id, "CCNB3_CAEEL/115-241")
        self.assertEqual(alignment.sequences[1].id, "T1FNQ9_HELRO/256-381")
        self.assertEqual(alignment.sequences[2].id, "CCNB3_DROME/308-433")
        self.assertEqual(alignment.sequences[3].id, "CCNB3_CHICK/149-274")
        self.assertEqual(alignment.sequences[4].id, "CCB31_ARATH/146-272")
        self.assertEqual(alignment.sequences[5].id, "CCNB1_SOYBN/197-322")
        self.assertEqual(alignment.sequences[6].id, "CCB11_ARATH/168-293")
        self.assertEqual(alignment.sequences[7].id, "CCB12_ARATH/185-310")
        self.assertEqual(alignment.sequences[8].id, "M1A5P4_SOLTU/185-310")
        self.assertEqual(alignment.sequences[9].id, "I1M770_SOYBN/186-311")
        self.assertEqual(alignment.sequences[10].id, "CCB14_ARATH/133-258")
        self.assertEqual(alignment.sequences[11].id, "B8A2G9_MAIZE/228-354")
        self.assertEqual(alignment.sequences[12].id, "B4FZZ7_MAIZE/193-318")
        self.assertEqual(alignment.sequences[13].id, "K7V0X7_MAIZE/191-316")
        self.assertEqual(alignment.sequences[14].id, "CCB21_ORYSJ/157-283")
        self.assertEqual(alignment.sequences[15].id, "CCB21_ARATH/174-300")
        self.assertEqual(alignment.sequences[16].id, "Q9XGI1_SOLLC/180-306")
        self.assertEqual(alignment.sequences[17].id, "C4J9B6_MAIZE/173-299")
        self.assertEqual(alignment.sequences[18].id, "CCB24_ARATH/180-307")
        self.assertEqual(alignment.sequences[19].id, "CCNB_DICDI/188-314")
        self.assertEqual(alignment.sequences[20].id, "O45926_CAEEL/34-159")
        self.assertEqual(alignment.sequences[21].id, "CCNB1_CAEEL/81-206")
        self.assertEqual(alignment.sequences[22].id, "P92162_BOMMO/260-387")
        self.assertEqual(alignment.sequences[23].id, "CCNB_DROME/260-387")
        self.assertEqual(alignment.sequences[24].id, "F6QF79_XENTR/139-264")
        self.assertEqual(alignment.sequences[25].id, "CCNB1_MOUSE/170-295")
        self.assertEqual(alignment.sequences[26].id, "Q28HA1_XENTR/132-257")
        self.assertEqual(alignment.sequences[27].id, "CCNB2_HUMAN/137-262")
        self.assertEqual(alignment.sequences[28].id, "CCNB2_CHICK/141-266")
        self.assertEqual(alignment.sequences[29].id, "A0BXX7_PARTE/85-212")
        self.assertEqual(alignment.sequences[30].id, "Q6BFS2_PARTE/82-210")
        self.assertEqual(alignment.sequences[31].id, "CG21_CANAL/208-334")
        self.assertEqual(alignment.sequences[32].id, "CGS5_YEAST/165-294")
        self.assertEqual(alignment.sequences[33].id, "CGS6_YEAST/123-252")
        self.assertEqual(alignment.sequences[34].id, "CG21_YEAST/212-337")
        self.assertEqual(alignment.sequences[35].id, "CG22_YEAST/232-357")
        self.assertEqual(alignment.sequences[36].id, "CG24_CANAL/234-360")
        self.assertEqual(alignment.sequences[37].id, "CG23_YEAST/171-297")
        self.assertEqual(alignment.sequences[38].id, "CG24_YEAST/210-336")
        self.assertEqual(alignment.sequences[39].id, "CG22_SCHPO/139-265")
        self.assertEqual(alignment.sequences[40].id, "M7PCR8_PNEMU/176-302")
        self.assertEqual(alignment.sequences[41].id, "CG21_EMENI/211-337")
        self.assertEqual(alignment.sequences[42].id, "CG23_SCHPO/206-332")
        self.assertEqual(alignment.sequences[43].id, "CG21_SCHPO/166-292")
        self.assertEqual(alignment.sequences[44].id, "REM1_SCHPO/145-271")
        self.assertEqual(alignment.sequences[45].id, "CCNF_MOUSE/282-406")
        self.assertEqual(alignment.sequences[46].id, "I1K835_SOYBN/85-214")
        self.assertEqual(alignment.sequences[47].id, "J3KY10_ORYBR/200-327")
        self.assertEqual(alignment.sequences[48].id, "M4E4A3_BRARP/154-281")
        self.assertEqual(alignment.sequences[49].id, "CCA22_ARATH/175-302")
        self.assertEqual(alignment.sequences[50].id, "Q39878_SOYBN/207-334")
        self.assertEqual(alignment.sequences[51].id, "T1EGC7_HELRO/166-291")
        self.assertEqual(alignment.sequences[52].id, "CCNA1_CAEEL/215-341")
        self.assertEqual(alignment.sequences[53].id, "CCNA_DROME/206-332")
        self.assertEqual(alignment.sequences[54].id, "W4XJF2_STRPU/207-333")
        self.assertEqual(alignment.sequences[55].id, "CCNA2_MOUSE/171-297")
        self.assertEqual(len(alignment.sequences[55].dbxrefs), 8)
        self.assertEqual(alignment.sequences[55].dbxrefs[0], "PDB; 4I3Z D; 181-307;")
        self.assertEqual(alignment.sequences[55].dbxrefs[1], "PDB; 4II5 D; 181-307;")
        self.assertEqual(alignment.sequences[55].dbxrefs[2], "PDB; 4I3Z B; 181-307;")
        self.assertEqual(alignment.sequences[55].dbxrefs[3], "PDB; 4II5 B; 181-307;")
        self.assertEqual(alignment.sequences[55].dbxrefs[4], "PDB; 3QHW B; 181-307;")
        self.assertEqual(alignment.sequences[55].dbxrefs[5], "PDB; 3QHW D; 181-307;")
        self.assertEqual(alignment.sequences[55].dbxrefs[6], "PDB; 3QHR D; 181-307;")
        self.assertEqual(alignment.sequences[55].dbxrefs[7], "PDB; 3QHR B; 181-307;")
        self.assertEqual(
            alignment.sequences[55].letter_annotations["secondary structure"],
            "HHHHHHHHHHHHT---TTGGGG-SS--HHHHHHHHHHHHHHHHHTT--HHHHHHHHHHHHHHHCC----CCCHHHHHHHHHHHHHHHH-SS---HHHHHHHTTTSS-HHHHHHHHHHHHHHHTT---",
        )
        self.assertEqual(alignment.sequences[56].id, "E9QJ66_DANRE/140-266")
        self.assertEqual(alignment.sequences[57].id, "CCNA1_HUMAN/214-340")
        self.assertEqual(alignment.sequences[58].id, "CG12_YEAST/43-195")
        self.assertEqual(alignment.sequences[59].id, "PUC1_SCHPO/99-227")
        self.assertEqual(alignment.sequences[60].id, "CG13_CANAL/44-172")
        self.assertEqual(alignment.sequences[61].id, "CG11_CANAL/44-172")
        self.assertEqual(alignment.sequences[62].id, "CGH2_SHV21/22-148")
        self.assertEqual(len(alignment.sequences[62].dbxrefs), 6)
        self.assertEqual(alignment.sequences[62].dbxrefs[0], "PDB; 1JOW A; 22-148;")
        self.assertEqual(alignment.sequences[62].dbxrefs[1], "PDB; 1BU2 A; 22-148;")
        self.assertEqual(alignment.sequences[62].dbxrefs[2], "PDB; 2EUF A; 22-148;")
        self.assertEqual(alignment.sequences[62].dbxrefs[3], "PDB; 4TTH A; 22-148;")
        self.assertEqual(alignment.sequences[62].dbxrefs[4], "PDB; 1XO2 A; 22-148;")
        self.assertEqual(alignment.sequences[62].dbxrefs[5], "PDB; 2F2C A; 22-148;")
        self.assertEqual(
            alignment.sequences[62].letter_annotations["secondary structure"],
            "HHHHHHHHHHTTS---SSTTTT-SSS-HHHHHHHHHHHHHHHHHTT--TTHHHHHHHHHHHHHHHS---TTTHHHHHHHHHHHHHHHHSSS---HHHHHHTTTTSS-HHHHHHHHHHHHHHTTT---",
        )
        self.assertEqual(alignment.sequences[63].id, "VCYCL_HHV8P/21-147")
        self.assertEqual(alignment.sequences[64].id, "CCND_CAEEL/72-201")
        self.assertEqual(alignment.sequences[65].id, "Q7KUZ5_DROME/153-280")
        self.assertEqual(alignment.sequences[66].id, "CCND1_RAT/26-153")
        self.assertEqual(alignment.sequences[67].id, "CCND2_MOUSE/24-151")
        self.assertEqual(alignment.sequences[68].id, "CCND3_HUMAN/26-153")
        self.assertEqual(len(alignment.sequences[68].dbxrefs), 2)
        self.assertEqual(alignment.sequences[68].dbxrefs[0], "PDB; 3G33 D; 26-153;")
        self.assertEqual(alignment.sequences[68].dbxrefs[1], "PDB; 3G33 B; 26-153;")
        self.assertEqual(
            alignment.sequences[68].letter_annotations["secondary structure"],
            "HHHHHHHHHGGGGS-SS--TTTSTTT--HHHHHHHHHHHHHHHHHTT--TTHHHHHHHHHHHHHHH----GGGHHHHHHHHHHHHHHHH-SS---TTHHHHHTTTSS-HHHHHHHHHHHHHHTTT---",
        )
        self.assertEqual(alignment.sequences[69].id, "Q9VZP3_DROME/42-165")
        self.assertEqual(alignment.sequences[70].id, "SSN8_YEAST/45-176")
        self.assertEqual(alignment.sequences[71].id, "CCC11_ORYSJ/4-144")
        self.assertEqual(alignment.sequences[72].id, "CCNT_DROME/42-176")
        self.assertEqual(alignment.sequences[73].id, "CCT12_ARATH/28-169")
        self.assertEqual(alignment.sequences[74].id, "Q9VE72_DROME/6-144")
        self.assertEqual(alignment.sequences[75].id, "PCL1_YEAST/19-152")
        self.assertEqual(alignment.sequences[76].id, "PCL2_YEAST/18-146")
        self.assertEqual(alignment.sequences[77].id, "PCL9_YEAST/19-146")
        self.assertEqual(alignment.sequences[78].id, "CCU41_ARATH/23-148")
        self.assertEqual(alignment.sequences[79].id, "Q9VKF0_DROME/205-327")
        self.assertEqual(alignment.sequences[80].id, "CCD11_ARATH/50-182")
        self.assertEqual(alignment.sequences[81].id, "CCD21_ARATH/65-197")
        self.assertEqual(alignment.sequences[82].id, "CCD41_ARATH/45-178")
        self.assertEqual(alignment.sequences[83].id, "Q9SMD4_SOLLC/51-182")
        self.assertEqual(alignment.sequences[84].id, "Q9S7H9_SOLLC/61-190")
        self.assertEqual(alignment.sequences[85].id, "CCD33_ARATH/59-186")
        self.assertEqual(alignment.sequences[86].id, "CCD61_ARATH/26-154")
        self.assertEqual(alignment.sequences[87].id, "CCNE_DROME/330-459")
        self.assertEqual(alignment.sequences[88].id, "CCNE2_MOUSE/112-239")
        self.assertEqual(alignment.sequences[89].id, "CCNE1_CHICK/112-239")
        self.assertEqual(alignment.sequences[90].id, "CCNE1_MOUSE/113-240")
        self.assertEqual(alignment.sequences[91].id, "A0A0R4IZF8_DANRE/117-244")
        self.assertEqual(alignment.sequences[92].id, "F6QUN0_XENTR/114-241")
        self.assertEqual(alignment.sequences[93].id, "W4XEA0_STRPU/126-253")
        self.assertEqual(alignment.sequences[94].id, "CCNE_CAEEL/232-360")
        self.assertEqual(alignment.sequences[0].annotations["accession"], "Q10654.3")
        self.assertEqual(alignment.sequences[1].annotations["accession"], "T1FNQ9.1")
        self.assertEqual(alignment.sequences[2].annotations["accession"], "Q9I7I0.1")
        self.assertEqual(alignment.sequences[3].annotations["accession"], "P39963.1")
        self.assertEqual(alignment.sequences[4].annotations["accession"], "Q9SA32.2")
        self.assertEqual(alignment.sequences[5].annotations["accession"], "P25011.1")
        self.assertEqual(alignment.sequences[6].annotations["accession"], "P30183.2")
        self.assertEqual(alignment.sequences[7].annotations["accession"], "Q39067.2")
        self.assertEqual(alignment.sequences[8].annotations["accession"], "M1A5P4.1")
        self.assertEqual(alignment.sequences[9].annotations["accession"], "I1M770.1")
        self.assertEqual(alignment.sequences[10].annotations["accession"], "O48790.1")
        self.assertEqual(alignment.sequences[11].annotations["accession"], "B8A2G9.1")
        self.assertEqual(alignment.sequences[12].annotations["accession"], "B4FZZ7.1")
        self.assertEqual(alignment.sequences[13].annotations["accession"], "K7V0X7.1")
        self.assertEqual(alignment.sequences[14].annotations["accession"], "Q7XSJ6.2")
        self.assertEqual(alignment.sequences[15].annotations["accession"], "Q39068.2")
        self.assertEqual(alignment.sequences[16].annotations["accession"], "Q9XGI1.1")
        self.assertEqual(alignment.sequences[17].annotations["accession"], "C4J9B6.1")
        self.assertEqual(alignment.sequences[18].annotations["accession"], "Q9SFW6.2")
        self.assertEqual(alignment.sequences[19].annotations["accession"], "P42524.1")
        self.assertEqual(alignment.sequences[20].annotations["accession"], "O45926.2")
        self.assertEqual(alignment.sequences[21].annotations["accession"], "Q10653.1")
        self.assertEqual(alignment.sequences[22].annotations["accession"], "P92162.1")
        self.assertEqual(alignment.sequences[23].annotations["accession"], "P20439.2")
        self.assertEqual(alignment.sequences[24].annotations["accession"], "F6QF79.3")
        self.assertEqual(alignment.sequences[25].annotations["accession"], "P24860.3")
        self.assertEqual(alignment.sequences[26].annotations["accession"], "Q28HA1.1")
        self.assertEqual(alignment.sequences[27].annotations["accession"], "O95067.1")
        self.assertEqual(alignment.sequences[28].annotations["accession"], "P29332.1")
        self.assertEqual(alignment.sequences[29].annotations["accession"], "A0BXX7.1")
        self.assertEqual(alignment.sequences[30].annotations["accession"], "Q6BFS2.1")
        self.assertEqual(alignment.sequences[31].annotations["accession"], "Q5ALY0.1")
        self.assertEqual(alignment.sequences[32].annotations["accession"], "P30283.1")
        self.assertEqual(alignment.sequences[33].annotations["accession"], "P32943.2")
        self.assertEqual(alignment.sequences[34].annotations["accession"], "P24868.1")
        self.assertEqual(alignment.sequences[35].annotations["accession"], "P24869.1")
        self.assertEqual(alignment.sequences[36].annotations["accession"], "Q5A0A9.1")
        self.assertEqual(alignment.sequences[37].annotations["accession"], "P24870.3")
        self.assertEqual(alignment.sequences[38].annotations["accession"], "P24871.2")
        self.assertEqual(alignment.sequences[39].annotations["accession"], "P36630.2")
        self.assertEqual(alignment.sequences[40].annotations["accession"], "M7PCR8.1")
        self.assertEqual(alignment.sequences[41].annotations["accession"], "P30284.1")
        self.assertEqual(alignment.sequences[42].annotations["accession"], "P10815.1")
        self.assertEqual(alignment.sequences[43].annotations["accession"], "P24865.2")
        self.assertEqual(alignment.sequences[44].annotations["accession"], "O14332.1")
        self.assertEqual(alignment.sequences[45].annotations["accession"], "P51944.2")
        self.assertEqual(alignment.sequences[46].annotations["accession"], "I1K835.1")
        self.assertEqual(alignment.sequences[47].annotations["accession"], "J3KY10.1")
        self.assertEqual(alignment.sequences[48].annotations["accession"], "M4E4A3.1")
        self.assertEqual(alignment.sequences[49].annotations["accession"], "Q147G5.1")
        self.assertEqual(alignment.sequences[50].annotations["accession"], "Q39878.1")
        self.assertEqual(alignment.sequences[51].annotations["accession"], "T1EGC7.1")
        self.assertEqual(alignment.sequences[52].annotations["accession"], "P34638.2")
        self.assertEqual(alignment.sequences[53].annotations["accession"], "P14785.3")
        self.assertEqual(alignment.sequences[54].annotations["accession"], "W4XJF2.1")
        self.assertEqual(alignment.sequences[55].annotations["accession"], "P51943.2")
        self.assertEqual(alignment.sequences[56].annotations["accession"], "E9QJ66.1")
        self.assertEqual(alignment.sequences[57].annotations["accession"], "P78396.1")
        self.assertEqual(alignment.sequences[58].annotations["accession"], "P20438.2")
        self.assertEqual(alignment.sequences[59].annotations["accession"], "P25009.1")
        self.assertEqual(alignment.sequences[60].annotations["accession"], "Q5A1N6.1")
        self.assertEqual(alignment.sequences[61].annotations["accession"], "Q59YH3.2")
        self.assertEqual(alignment.sequences[62].annotations["accession"], "Q01043.1")
        self.assertEqual(alignment.sequences[63].annotations["accession"], "Q77Q36.1")
        self.assertEqual(alignment.sequences[64].annotations["accession"], "Q9U2M5.1")
        self.assertEqual(alignment.sequences[65].annotations["accession"], "Q7KUZ5.1")
        self.assertEqual(alignment.sequences[66].annotations["accession"], "P39948.1")
        self.assertEqual(alignment.sequences[67].annotations["accession"], "P30280.1")
        self.assertEqual(alignment.sequences[68].annotations["accession"], "P30281.2")
        self.assertEqual(alignment.sequences[69].annotations["accession"], "Q9VZP3.1")
        self.assertEqual(alignment.sequences[70].annotations["accession"], "P47821.1")
        self.assertEqual(alignment.sequences[71].annotations["accession"], "P93411.1")
        self.assertEqual(alignment.sequences[72].annotations["accession"], "O96433.2")
        self.assertEqual(alignment.sequences[73].annotations["accession"], "Q56YF8.2")
        self.assertEqual(alignment.sequences[74].annotations["accession"], "Q9VE72.1")
        self.assertEqual(alignment.sequences[75].annotations["accession"], "P24867.1")
        self.assertEqual(alignment.sequences[76].annotations["accession"], "P25693.2")
        self.assertEqual(alignment.sequences[77].annotations["accession"], "Q12477.1")
        self.assertEqual(alignment.sequences[78].annotations["accession"], "O80513.1")
        self.assertEqual(alignment.sequences[79].annotations["accession"], "Q9VKF0.1")
        self.assertEqual(alignment.sequences[80].annotations["accession"], "P42751.3")
        self.assertEqual(alignment.sequences[81].annotations["accession"], "P42752.3")
        self.assertEqual(alignment.sequences[82].annotations["accession"], "Q8LGA1.2")
        self.assertEqual(alignment.sequences[83].annotations["accession"], "Q9SMD4.1")
        self.assertEqual(alignment.sequences[84].annotations["accession"], "Q9S7H9.1")
        self.assertEqual(alignment.sequences[85].annotations["accession"], "Q9SN11.1")
        self.assertEqual(alignment.sequences[86].annotations["accession"], "Q9ZR04.1")
        self.assertEqual(alignment.sequences[87].annotations["accession"], "P54733.2")
        self.assertEqual(alignment.sequences[88].annotations["accession"], "Q9Z238.1")
        self.assertEqual(alignment.sequences[89].annotations["accession"], "P49707.1")
        self.assertEqual(alignment.sequences[90].annotations["accession"], "Q61457.2")
        self.assertEqual(
            alignment.sequences[91].annotations["accession"], "A0A0R4IZF8.1"
        )
        self.assertEqual(alignment.sequences[92].annotations["accession"], "F6QUN0.2")
        self.assertEqual(alignment.sequences[93].annotations["accession"], "W4XEA0.1")
        self.assertEqual(alignment.sequences[94].annotations["accession"], "O01501.2")
        self.assertEqual(
            alignment.sequences[0].seq,
            "GIFDYYRHREVHFRVRKYLHKHPEVDVKTRAILIDWMVEIQETFELNHETLYNAVKLTDMYLCKTKNVDKNTIQKLACVAIFIAAKYDERSPPLVDDLIYLSGDRFSRDELLAMERELFATVGYDLG",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "DIFDYYRDREVKFRIPDYMFQQTDLTPSMRAILVDWLVEVQQSFELNHETLYMAVKLIDIFSSKVTIKRNKLQLIGAVALNLACKFEERCPPMLDDFVYVCDDAYPRQEFLKMEELVFQAVGFDIG",
        )
        self.assertEqual(
            alignment.sequences[2].seq,
            "DIFNYLKVREAEFPIADYMPRQIHLTTWMRTLLVDWMVEVQETFELNHETLYLAVKIVDLYLCREVINKEKLQLLGAAAFFIACKYDERQPPLIEDFLYICDGAYNHDELVRMERETLRVIKYDLG",
        )
        self.assertEqual(
            alignment.sequences[3].seq,
            "EIFDYMREREEKFLLPDYMEKQSDISRDMRAILVDWMVEVQENFELNHETLYLAVKLVDHYLVEVVSMRDKLQLIGSTAVLIASKFEERCPPCVDDFLYICDDAYKREELIAMETSILRTLNFDIN",
        )
        self.assertEqual(
            alignment.sequences[4].seq,
            "DIYQFYWTAEALNPALGHYLSAHAEVSPVTRGILINWLIEVHFKFDLMHETLYLTMDLLDRYLSQVPIHKNEMQLIGLTALLLASKYEDYWHPRIKDLISISAESYTREQILGMERSMLKQLKFRLN",
        )
        self.assertEqual(
            alignment.sequences[5].seq,
            "DIYKFYKLVENESRPHDYIGSQPEINERMRAILVDWLIDVHTKFELSLETLYLTINIIDRFLAVKTVPRRELQLVGISAMLMASKYEEIWPPEVNDFVCLSDRAYTHEHILTMEKTILNKLEWTLT",
        )
        self.assertEqual(
            alignment.sequences[6].seq,
            "DIYSFYKSVESEWRPRDYMASQPDINEKMRLILVEWLIDVHVRFELNPETFYLTVNILDRFLSVKPVPRKELQLVGLSALLMSAKYEEIWPPQVEDLVDIADHAYSHKQILVMEKTILSTLEWYLT",
        )
        self.assertEqual(
            alignment.sequences[7].seq,
            "DMYSFYKEVEKESQPKMYMHIQTEMNEKMRAILIDWLLEVHIKFELNLETLYLTVNIIDRFLSVKAVPKRELQLVGISALLIASKYEEIWPPQVNDLVYVTDNAYSSRQILVMEKAILGNLEWYLT",
        )
        self.assertEqual(
            alignment.sequences[8].seq,
            "DIYKFYKLTEDENRPCDYMDSQPEINDRVRAILVDWLIEAHKRFELRPESLYLTVNIMDRFLSEETVPRRELQLLCISSMLIACKYEEIWAPEVNDFLTITDNAYVRDQILLMEKVILGKLEWYLT",
        )
        self.assertEqual(
            alignment.sequences[9].seq,
            "DIYKFYKETEEDGCVHDYMGSQPDINAKMRSILVDWLIEVHRKFELMPETLYLTLNIVDRFLSVKAVPRRELQLVGISSMLIASKYEEIWAPEVNDFVCISDNAYVSEQVLMMEKTILRKLEWYLT",
        )
        self.assertEqual(
            alignment.sequences[10].seq,
            "DIFKFYRTVEEEGGIKDYIGSQPEINEKMRSILIDWLVDVHRKFELMPETLYLTINLVDRFLSLTMVHRRELQLLGLGAMLIACKYEEIWAPEVNDFVCISDNAYNRKQVLAMEKSILGQVEWYIT",
        )
        self.assertEqual(
            alignment.sequences[11].seq,
            "DIYRFYKSTEGTCLPLSSYMSSQAEISERMRAILIDWIIEVQYRLTLMPETLYLTVYIIDQYLSMESVPRKELQLVGISAMLIASKYEEIWAPLVKDLMCLCDNAFTRDQILTKEKAILDMLHWNLT",
        )
        self.assertEqual(
            alignment.sequences[12].seq,
            "DIYTFYKIAQHDRRPCDYIDTQVEINPKMRAILAGWIIEVHHKFELMPETLYLTMYIIDQYLSLQPVLRRELQLVGVSAMLIACKYEEIWAPEVNDFILISDSAYSREQILSMEKGILNSLEWNLT",
        )
        self.assertEqual(
            alignment.sequences[13].seq,
            "DIYTFYKTAQHESRPIDYMGNQPELSPRMRSILADWLIESHRRFQLMPETLYLTIYIVDRYLSLQPTPRRELQLVGVAALLIACKYEEIWAPEVNDLIHIADGAFNRSQILAAEKAILNSMEWNLT",
        )
        self.assertEqual(
            alignment.sequences[14].seq,
            "ELYKFYRENEEMSCVQPDYMSSQGDINEKMRAILIDWLIEVHHKFELMDETLFLTVNIVDRFLEKQVVPRKKLQLVGVTAMLLACKYEEVAVPVVEDLVLISDRAYTKGQILEMEKLILNTLQFNMS",
        )
        self.assertEqual(
            alignment.sequences[15].seq,
            "DLYAFYRTMERFSCVPVDYMMQQIDLNEKMRAILIDWLIEVHDKFDLINETLFLTVNLIDRFLSKQNVMRKKLQLVGLVALLLACKYEEVSVPVVEDLVLISDKAYTRNDVLEMEKTMLSTLQFNIS",
        )
        self.assertEqual(
            alignment.sequences[16].seq,
            "DLFANYRTMEVNSCASPYYMAQQADINERMRSILIDWLIEVHHKFELREETLFLTVNLIDRFLEKQGIVRKKLQLVGLVAMLLACKYEEVCAPLVEDLVLISDKAYTRKEVLEMESMMLNTLQFNMS",
        )
        self.assertEqual(
            alignment.sequences[17].seq,
            "EIYRFYRKTEGASCVPTNYMSSQTDINEKMRGILIDWLIEVHYKLELLEETLFLTVNIIDRFLARENVVRKKLQLAGVTAMLLACKYEEVSVPVVEDLILICDRAYTRADILEMERRIVNTLNFNMS",
        )
        self.assertEqual(
            alignment.sequences[18].seq,
            "DIYCFYKKNECRSCVPPNYMENQHDINERMRGILFDWLIEVHYKFELMEETLYLTINLIDRFLAVHQHIARKKLQLVGVTAMLLACKYEEVSVPVVDDLILISDKAYTRTEILDMEKLMANTLQFNFC",
        )
        self.assertEqual(
            alignment.sequences[19].seq,
            "EIFAYYREKEQIDKIDKDYIKNQYHINERMRAILVDWMMAVHVRFKLLSETFFLSVNIVDRYLAKVMIPVTKLQLVGITAILLACKYEEIYSPQIKDFVHTSDDACTHAEVIDMERQILSTLQFHMS",
        )
        self.assertEqual(
            alignment.sequences[20].seq,
            "DIYNYLVHHEKKYVLDDSFINGGNVNSKMRRILVDWLIQVHLRFHLTPETLHLTIFVLDRIIVKNIVSKAEFQLLGVAALFVASKFEDIYLPDILEYEMITDNTFSKKQIMAMEQTILNALNFDLS",
        )
        self.assertEqual(
            alignment.sequences[21].seq,
            "DIYKYLVHHEKKYLLEECFMEGGEPTPKMRRILVDWLVQVHVRFHLTPETLHLTVFILDRMLQKKVTSKADLQLLGISAMFVASKFEEVYLPDIHDYEFITENTYSKKQILAMEQTILNSLNFDLS",
        )
        self.assertEqual(
            alignment.sequences[22].seq,
            "DIYKYLTELEEKYSIEPDHLKKQTVITGKMRATLIDWLVEVQRQFSLVLETFHLTVGIIDRYLQVVPNVQRNQLQLVGVTAMFIASKYEEIYAPDVGDFVYVTDNAYTKSDVFRCERDIMCKLGFCLA",
        )
        self.assertEqual(
            alignment.sequences[23].seq,
            "DIYDYLYQVELEQPIHKDHLAGQKEVSHKMRAVLIDWINEVHLQFHLAAETFQLAVAIIDRYLQVVKDTKRTYLQLVGVTALFIATKYEELFPPAIGDFVFITDDTYTARQIRQMELQIFKAIDCNLS",
        )
        self.assertEqual(
            alignment.sequences[24].seq,
            "DIYCYLRSLENAQAVRQNYLHGQEVTGNMRAILIDWLVQVQMKFRLLQETMFMTVGIIDRFLQDHPVPKNQLQLVGVTAMFLAAKYEEMYPPEIGDFTFVTDHTYTKAQIRDMEMKVLRVLKFAIG",
        )
        self.assertEqual(
            alignment.sequences[25].seq,
            "DIYAYLRQLEEEQSVRPKYLQGREVTGNMRAILIDWLIQVQMKFRLLQETMYMTVSIIDRFMQNSCVPKKMLQLVGVTAMFIASKYEEMYPPEIGDFAFVTNNTYTKHQIRQMEMKILRVLNFSLG",
        )
        self.assertEqual(
            alignment.sequences[26].seq,
            "DIYNYLKQLEVQQSVRPCYLEGKEINERMRAILVDWIVQVHSRFQLLQETLYMGIAIMDRFLQVQPVSRSKLQLVGVTSLLVASKYEEMYTPEVADFVYITDNAYTASQIREMEMIILRVLNFDLG",
        )
        self.assertEqual(
            alignment.sequences[27].seq,
            "DIYQYLRQLEVLQSINPHFLDGRDINGRMRAILVDWLVQVHSKFRLLQETLYMCVGIMDRFLQVQPVSRKKLQLVGITALLLASKYEEMFSPNIEDFVYITDNAYTSSQIREMETLILKELKFELG",
        )
        self.assertEqual(
            alignment.sequences[28].seq,
            "DIYLYLRQLELQQSVRPHYLDGKTINGRMRAILVDWLVQVHSRFQLLQETLYMCVAVMDRFLQSHPVPRKRLQLVGVTALLLASKYEEMYSPDIADFVYITDNAYNSAEVREMEITILKELNFDLG",
        )
        self.assertEqual(
            alignment.sequences[29].seq,
            "EILQHLLIEENKYTINQYMTPEQQPDINIKMRAILVDWLIDVHAKFKLKDETLYITISLIDRYLALAQVTRMRLQLVGVAALFIACKYEEIYPPALKDFVYITDNAYVKSDVLEMEGLMLQALNFNIC",
        )
        self.assertEqual(
            alignment.sequences[30].seq,
            "EIYTYLLTQEEKYLVSNNYMNEQQQPDLNARMRAILLDWLIDVHLKFKLRDETLYVTTYLIDRFLNFKTTTRQQLQLVGVASLFIACKYEEIYPPDLKDFVYITDNAYTKQDVLEMEGQILQTLDFSIT",
        )
        self.assertEqual(
            alignment.sequences[31].seq,
            "EIFSYYYELETRMLPDPQYLFKQTLLKPRMRSILVDWLVEMHLKFKLLPESLFLAVNVMDRFMSVEVVQIDKLQLLATAALFTAAKYEEVFSPSVKNYAYFTDGSYTPEEVVQAEKYMLTILNFDLN",
        )
        self.assertEqual(
            alignment.sequences[32].seq,
            "EIFAFLYRRELETLPSHNYLLDKTSKYYLRPSMRTILVDWLVEVHEKFQCYPETLFLSINLMDRFLAKNKVTMNKLQLLAVTSLFIAAKFEEVNLPKLAEYAYITDGAASKNDIKNAEMFMLTSLEFNIG",
        )
        self.assertEqual(
            alignment.sequences[33].seq,
            "SIFSHLYEKEIQMLPTHNYLMDTQSPYHLKSSMRALLIDWLVEVHEKFHCLPETLFLAINLLDRFLSQNVVKLNKLQLLCITCLFIACKFEEVKLPKITNFAYVTDGAATVEGIRKAELFVLSSLGYNIS",
        )
        self.assertEqual(
            alignment.sequences[34].seq,
            "DIFDYLHHLEIITLPNKANLYKHKNIKQNRDILVNWIIKIHNKFGLLPETLYLAINIMDRFLCEEVVQLNRLQLVGTSCLFIASKYEEIYSPSIKHFAYETDGACSVEDIKEGERFILEKLDFQIS",
        )
        self.assertEqual(
            alignment.sequences[35].seq,
            "DIFEYLHQLEVITLPKKEDLYQHRNIHQNRDILVNWLVKIHNKFGLLPETLYLAINIMDRFLGKELVQLDKLQLVGTSCLFIASKYEEVYSPSIKHFASETDGACTEDEIKEGEKFILKTLKFNLN",
        )
        self.assertEqual(
            alignment.sequences[36].seq,
            "EIFNYLHELENKFTPDPNYMDFQDDLKWEMRAVLIDWVVQVHARFNLFSETLYLTVNYIDRFLSKRRVSLSRFQLVGAVALFIAAKYEEINCPTVQEIAYMADNAYSIDEFLKAERFMIDVLEFDLG",
        )
        self.assertEqual(
            alignment.sequences[37].seq,
            "EIFEYMRKLEDLYKPNPYYMDKQPELRWSFRSTLIDWIVQVHEKFQLLPETLYLCINIIDRYLCKEVVPVNKFQLVGAASLFIAAKYEEINCPTIKDFVYMSENCYSRNDLLDAERTILNGLEFELG",
        )
        self.assertEqual(
            alignment.sequences[38].seq,
            "DIFYYLRELEVKYRPNPYYMQNQVELTWPFRRTMIDWLVQLHFRFQLLPETLYLTINIVDRFLSKKTVTLNRFQLVGVSALFIAAKFEEINCPTLDDLVYMLENTYTRDDIIRAEQYMIDTLEFEIG",
        )
        self.assertEqual(
            alignment.sequences[39].seq,
            "EIFEYIRKLDLKCLPNPKYMDQQKELTWKMREILNEWLVEIHSNFCLMPETLYLAVNIIDRFLSRRSCSLSKFQLTGITALLIASKYEEVMCPSIQNFVYMTDGAFTVEDVCVAERYMLNVLNFDLS",
        )
        self.assertEqual(
            alignment.sequences[40].seq,
            "EIMSYMRELEVLTLPLPDYMDRQKELQWKMRGILVDWLIEVHAKFRLLPETLFLSVNIIDRFLSLRVCSLPKLQLVGITALFIAAKYEEVMCPSIQNFMYMADGGYTNEEILKAEQYVLQVLGYDMS",
        )
        self.assertEqual(
            alignment.sequences[41].seq,
            "EIFDYLRELEMETLPNPDYIDHQPDLEWKMRGILVDWLIEVHTRFRLLPETLFLAVNIIDRFLSAEVVALDRLQLVGVAAMFIASKYEEVLSPHVANFSHVADETFSDKEILDAERHILATLEYNMS",
        )
        self.assertEqual(
            alignment.sequences[42].seq,
            "DIFEYLNELEIETMPSPTYMDRQKELAWKMRGILTDWLIEVHSRFRLLPETLFLAVNIIDRFLSLRVCSLNKLQLVGIAALFIASKYEEVMCPSVQNFVYMADGGYDEEEILQAERYILRVLEFNLA",
        )
        self.assertEqual(
            alignment.sequences[43].seq,
            "EIFHYMQSLERKLAPPPNYMSVQQEIDWVTRHMLVDWIVQVQIHFRLLPETLFLAVNLIDRFLSIKVVSLQKVQLVGLSALLIACKYEEIHPPSIYNFAHVVQGIFTVDEIIRAERYMLMLLDFDIS",
        )
        self.assertEqual(
            alignment.sequences[44].seq,
            "EILSHMEKLEIRFMPDYRHMSAQPYYVTEMRASVINWIVGVHTCINLLPESLFLSINVLDRFLSLQNVPASKMKLCGATALFIACKYEEIHPPTVKDLEIVLEGEWIGEDICGMEKYMLMVLQYQLG",
        )
        self.assertEqual(
            alignment.sequences[45].seq,
            "SVCQLFQASQAVNKQQIFSVQKGLSDTMRYILIDWLVEVATMKDFTSLCLHLTVECVDRYLRRRLVPRYKLQLLGIACMVICTRFISKEILTIREAVWLTDNTYKYEDLVRVMGEIISALEGKIR",
        )
        self.assertEqual(
            alignment.sequences[46].seq,
            "DIHGYLREMEMQNKRRPMVDYIEKVQKIVTPTMRAILVDWLVEVAVEYKLLSDTLHLSVSYIDRFLSVNPVSKSRLQLLGVSSMLIAAKYEEMDPPGVDEFCSITDHTYDKTEVVKMEADILKSLKFEMG",
        )
        self.assertEqual(
            alignment.sequences[47].seq,
            "DIYMHLREAETRKRPSTDFMETIQKDVNPSMRAILIDWLVEVAEEYRLVPDTLYLTVNYIDRYLSGNEINRQRLQLLGVACMLIAAKYEEICAPQVEEFCYITDNTYFRDEVLEMEASVLNYLKFEMT",
        )
        self.assertEqual(
            alignment.sequences[48].seq,
            "DIYNHLRAAEAKKQPAVDYMATVQKDVNSTMRGILVDWLVEVSEEYRLVPETLYLTVNYIDRYLSGNVISRQKLQLLGVACMMIAAKYEEVCAPQVEEFCYITDNTYLKDEVLDMESAVLNYLKFEMS",
        )
        self.assertEqual(
            alignment.sequences[49].seq,
            "DIYDNIHVAELQQRPLANYMELVQRDIDPDMRKILIDWLVEVSDDYKLVPDTLYLTVNLIDRFLSNSYIERQRLQLLGVSCMLIASKYEELSAPGVEEFCFITANTYTRPEVLSMEIQILNFVHFRLS",
        )
        self.assertEqual(
            alignment.sequences[50].seq,
            "DIYSNIRVTELQRKPLTNYMDKLQKDINPSMRGILVDWLVEVSEEYKLVPDTLYLTVNLIDRYLSTRLIQKQKLQLLGVTCMLIASKYEEMCAPRVEEFCFITDNTYTKEEVLKMEREVLNLVHFQLS",
        )
        self.assertEqual(
            alignment.sequences[51].seq,
            "DILTYGKEAEQRYMAKANYMERQSDINHSMRSILVDWLVEVADEYKLKRETFFLAVNYIDRFLSMMSVIRCRLQLLGAAAMFIAAKYEEIYPPDVAEFVYITDDTYTMKQVLQMEQAILKTLNFLV",
        )
        self.assertEqual(
            alignment.sequences[52].seq,
            "DIIKYMLHRQTKNRASHECFDIQSQVNEEMRTILIDWFSDVVKEYNFQKETFHLAVSLVDRALSMFNIDKMRFQLVGTTSMMIAVKYEEIFPPEIEDFALITDNTYRVPDILLMERFLLGKFDFVVA",
        )
        self.assertEqual(
            alignment.sequences[53].seq,
            "DILEYFRESEKKHRPKPLYMRRQKDISHNMRSILIDWLVEVSEEYKLDTETLYLSVFYLDRFLSQMAVVRSKLQLVGTAAMYIAAKYEEIYPPEVGEFVFLTDDSYTKAQVLRMEQVILKILSFDLC",
        )
        self.assertEqual(
            alignment.sequences[54].seq,
            "EIYQYLKTAESKHRPKHGYMRKQPDITNSMRCILVDWLVEVSEEYRLHNETLYLAAAFIDRFLSQMSVLRAKLQLVGTASMFVASKYEEIYPPDVKEFVYITDDTYSIKQVLRMEHLILKVLSFDLA",
        )
        self.assertEqual(
            alignment.sequences[55].seq,
            "DIHTYLREMEVKCKPKVGYMKRQPDITNSMRAILVDWLVEVGEEYKLQNETLHLAVNYIDRFLSSMSVLRGKLQLVGTAAMLLASKFEEIYPPEVAEFVYITDDTYSKKQVLRMEHLVLKVLAFDLA",
        )
        self.assertEqual(
            alignment.sequences[56].seq,
            "DIHRYLRECEVKYRPKPGYMRKQPDITNCMRVILVDWLVEVGEEYKLCSETLYLAVNYLDRFLSCMSVLRGKLQLVGTAAILLAAKYEEVYPPEVDEFVYITDDTYTKKQLLRMEQHLLRVLAFDMT",
        )
        self.assertEqual(
            alignment.sequences[57].seq,
            "EIYQYLREAEIRHRPKAHYMKKQPDITEGMRTILVDWLVEVGEEYKLRAETLYLAVNFLDRFLSCMSVLRGKLQLVGTAAMLLASKYEEIYPPEVDEFVYITDDTYTKRQLLKMEHLLLKVLAFDLT",
        )
        self.assertEqual(
            alignment.sequences[58].seq,
            "EISTNVIAQSCKFKPNPKLIDQQPEMNPVETRSNIITFLFELSVVTRVTNGIFFHSVRLYDRYCSKRIVLRDQAKLVVATCLWLAAKTWGGCNHIINNVVIPTGGRFYGPNPRARIPRLSELVHYCGDGQVFDESMFLQMERHILDTLNWNIY",
        )
        self.assertEqual(
            alignment.sequences[59].seq,
            "DIIHHLITREKNFLLNVHLSNQQPELRWSMRPALVNFIVEIHNGFDLSIDTLPLSISLMDSYVSRRVVYCKHIQLVACVCLWIASKFHETEDRVPLLQELKLACKNIYAEDLFIRMERHILDTLDWDIS",
        )
        self.assertEqual(
            alignment.sequences[60].seq,
            "EMLHHLLSVEAKTLPNLSLIEQQPEIKLGMRPLLLDFLMEVITILSLSRSTFPLTVNLIDRYCSTRIVKKQHYQLLGLTSLWISCKNLDSKFKVPTLNDLRKICVDSYYKELFVEMEKHILKSLEWVVN",
        )
        self.assertEqual(
            alignment.sequences[61].seq,
            "DIVNTLSQLESLTLVNPAMIDLQPEIQWFMRPFLLDFLIELHSSFKLQPTTLFLCLNIIDRYCAKRIVFKRHYQLVGCTALWIASKYEDKKSRVPTLKELTIMCRNAYDEEMFVQMEMHILSTLDWSIG",
        )
        self.assertEqual(
            alignment.sequences[62].seq,
            "RVLNNLKLRELLLPKFTSLWEIQTEVTVDNRTILLTWMHLLCESFELDKSVFPLSVSILDRYLCKKQGTKKTLQKIGAACVLIGSKIRTVKPMTVSKLTYLSCDCFTNLELINQEKDILEALKWDTE",
        )
        self.assertEqual(
            alignment.sequences[63].seq,
            "IFYNILEIEPRFLTSDSVFGTFQQSLTSHMRKLLGTWMFSVCQEYNLEPNVVALALNLLDRLLLIKQVSKEHFQKTGSACLLVASKLRSLTPISTSSLCYAAADSFSRQELIDQEKELLEKLAWRTE",
        )
        self.assertEqual(
            alignment.sequences[64].seq,
            "DMRAFYNCMEYEEALQPNYHYFTGVQENITPFHREQAIDWIYDVAKEENCDGDVFLLAVSLIDRFMSVQNILKHDIQMIAGVALFIASKLKAPHPMTASKIAYYSDNSCPIDMILQWELLIVTTLQWETE",
        )
        self.assertEqual(
            alignment.sequences[65].seq,
            "LENFLKVEEKHHKIPDTYFSIQKDITPPMRKIVAEWMMEVCAEENCQEEVVLLALNYMDRFLSSKSVRKTQLQILAAACLLLASKLREPSCRALSVDLLVVYTDNSIYKDDLIKWELYVLSRLGWDLS",
        )
        self.assertEqual(
            alignment.sequences[66].seq,
            "RVLRAMLKTEETCAPSVSYFKCVQREIVPSMRKIVATWMLEVCEEQKCEEEVFPLAMNYLDRFLSLEPLKKSRLQLLGATCMFVASKMKETIPLTAEKLCIYTDNSIRPEELLQMELLLVNKLKWNLA",
        )
        self.assertEqual(
            alignment.sequences[67].seq,
            "RVLQNLLTIEERYLPQCSYFKCVQKDIQPYMRRMVATWMLEVCEEQKCEEEVFPLAMNYLDRFLAGVPTPKTHLQLLGAVCMFLASKLKETIPLTAEKLCIYTDNSVKPQELLEWELVVLGKLKWNLA",
        )
        self.assertEqual(
            alignment.sequences[68].seq,
            "RVLQSLLRLEERYVPRASYFQCVQREIKPHMRKMLAYWMLEVCEEQRCEEEVFPLAMNYLDRYLSCVPTRKAQLQLLGAVCMLLASKLRETTPLTIEKLCIYTDHAVSPRQLRDWEVLVLGKLKWDLA",
        )
        self.assertEqual(
            alignment.sequences[69].seq,
            "DIFLTMREQELSRRPLFYLSPQLNERRRMLQLLKLATSAHKLSRCALHLAVYYMDRFVDYYKIRPDKLLLVAITCLHIAAQIENTDAFIPRYSEMNRLVKNAYTAFEYKAVERKILCFLNFELI",
        )
        self.assertEqual(
            alignment.sequences[70].seq,
            "DSKQNGIEQSITKNIPITHRDLHYDKDYNLRIYCYFLIMKLGRRLNIRQYALATAHIYLSRFLIKASVREINLYMLVTTCVYLACKVEECPQYIRTLVSEARTLWPEFIPPDPTKVTEFEFYLLEELESYLI",
        )
        self.assertEqual(
            alignment.sequences[71].seq,
            "NFWTSSHCKQLLDQEDVDKVPQADSDRGITLEEFRLVKIHMSFHIWRLAQQVKVRQRVIATAVTYFRRVYTRKSMTEYDPRLVAPTCLYLASKVEESTVQARLLVFYIKKMCASDEKYRFEIKDILEMEMKLLEALDYYLV",
        )
        self.assertEqual(
            alignment.sequences[72].seq,
            "DKIWYFSNDQLANSPSRRCGIKGDDELQYRQMTAYLIQEMGQRLQVSQLCINTAIVYMHRFYAFHSFTHFHRNSMASASLFLAAKVEEQPRKLEHVIRAANKCLPPTTEQNYAELAQELVFNENVLLQTLGFDVA",
        )
        self.assertEqual(
            alignment.sequences[73].seq,
            "IIPWFFSREEIERNSPSRRDGIDLKTETRLRDSYCTFLEILGERLKVPQVTIATAIFFCHRFFLRQSHAKNDRQTIATVCMLLAGKVEETPVTLEDVIIASYERIHKKDLAGAQRKEVYDQQKELVLIGEELVLSTLNFDLC",
        )
        self.assertEqual(
            alignment.sequences[74].seq,
            "DVMSMQQHVELNKAQTMKPIDYRKMNKPGVVPMYIFECAAKLKMKPLTAACAAIVFHRFFREVKASDYDEFLIAAGSLYLAGKIKEDESVKIRDVINVAYCTLNRGNDPVDLNDEYWSMRDAIVQAELLITRTLCFDLN",
        )
        self.assertEqual(
            alignment.sequences[75].seq,
            "DIIKFLTDTTLRVVPSSNYPTPPGSPGEKHLTRLPSLMTFITRLVRYTNVYTPTLLTAACYLNKLKRILPRDATGLPSTIHRIFLACLILSAKFHNDSSPLNKHWARYTDGLFTLEDINLMERQLLQLLNWDLR",
        )
        self.assertEqual(
            alignment.sequences[76].seq,
            "EMVQYLASTTASIIKIKKTNSMIDIALPAPPLTKFINRLIKHSNVQTPTLMATSVYLAKLRSIIPSNVYGIETTRHRIFLGCLILAAKTLNDSSPLNKHWAEYTDGLLILREVNTIERELLEYFDWDVT",
        )
        self.assertEqual(
            alignment.sequences[77].seq,
            "EMIQFLATSTASIIKIRENNNPIQGCRPPDLSIFIKNVVIQSNVQTPTLMATSVYLNKLKSVIPKNVYGINTTRHRIFLGCLILAAKTLNDSSPWNKHWTTYTEGLLRIREVNTIERELLEYLNWDVR",
        )
        self.assertEqual(
            alignment.sequences[78].seq,
            "RVAESNDLTRRVATQSQRVSVFHGLSRPTITIQSYLERIFKYANCSPSCFVVAYVYLDRFTHRQPSLPINSFNVHRLLITSVMVAAKFLDDLYYNNAYYAKVGGISTKEMNFLELDFLFGLGFELN",
        )
        self.assertEqual(
            alignment.sequences[79].seq,
            "DIFDEKLHPLTHDQVPDNYDTHNPEHRQIYKFVRTLFNAAQLTAECAIITLVYLERLLTYAELDVGPCNWKRMVLGAILLASKVWDDQAVWNVDYCQILKDITVEDMNELERQFLELLQFNIN",
        )
        self.assertEqual(
            alignment.sequences[80].seq,
            "DSIACFIEDERHFVPGHDYLSRFQTRSLDASAREDSVAWILKVQAYYNFQPLTAYLAVNYMDRFLYARRLPETSGWPMQLLAVACLSLAAKMEEILVPSLFDFQVAGVKYLFEAKTIKRMELLVLSVLDWRLR",
        )
        self.assertEqual(
            alignment.sequences[81].seq,
            "DRIKEMLVREIEFCPGTDYVKRLLSGDLDLSVRNQALDWILKVCAHYHFGHLCICLSMNYLDRFLTSYELPKDKDWAAQLLAVSCLSLASKMEETDVPHIVDLQVEDPKFVFEAKTIKRMELLVVTTLNWRLQ",
        )
        self.assertEqual(
            alignment.sequences[82].seq,
            "EIIMEMVEKEKQHLPSDDYIKRLRSGDLDLNVGRRDALNWIWKACEVHQFGPLCFCLAMNYLDRFLSVHDLPSGKGWILQLLAVACLSLAAKIEETEVPMLIDLQVGDPQFVFEAKSVQRMELLVLNKLKWRLR",
        )
        self.assertEqual(
            alignment.sequences[83].seq,
            "EELTSLFSKETEYEISYNVLEKNQSFISSRRESVEWILKTTAYYSFSAQTGFLAVNYFDRFLLFSFNQSLNHKPWMNQLVAVTCLSLAAKVEETDVPLLLDLQVEESGFLFESKTIQRMEMLILSTLKWKMN",
        )
        self.assertEqual(
            alignment.sequences[84].seq,
            "DELATLLSKENEFHLGFQSLISDGSLMGARKEALDWMLRVIAYYGFTATTAVLAVNYFDRFVSGWCFQKDKPWMSQLAAVACLSIAAKVEETQVPLLLDLQVADSRFVFEAKTIQRMELLVLSTLKWKMN",
        )
        self.assertEqual(
            alignment.sequences[85].seq,
            "DELSTLISKQEPCLYDEILDDEFLVLCREKALDWIFKVKSHYGFNSLTALLAVNYFDRFITSRKFQTDKPWMSQLTALACLSLAAKVEEIRVPFLLDFQVEEARYVFEAKTIQRMELLVLSTLDWRMH",
        )
        self.assertEqual(
            alignment.sequences[86].seq,
            "TLPHSLFLVEFQHMPSSHYFHSLKSSAFLLSNRNQAISSITQYSRKFDDPSLTYLAVNYLDRFLSSEDMPQSKPWILKLISLSCVSLSAKMRKPDMSVSDLPVEGEFFDAQMIERMENVILGALKWRMR",
        )
        self.assertEqual(
            alignment.sequences[87].seq,
            "DVWRLMCHRDEQDSRLRSISMLEQHPGLQPRMRAILLDWLIEVCEVYKLHRETFYLAVDYLDRYLHVAHKVQKTHLQLIGITCLFVAAKVEEIYPPKIGEFAYVTDGACTERDILNHEKILLQALDWDIS",
        )
        self.assertEqual(
            alignment.sequences[88].seq,
            "EVWQNMLQKENRYVHDKHFQVLHSDLEPQMRSILLDWLLEVCEVYTLHRETFYLAQDFFDRFMLTQKDVNKNMLQLIGITSLFIASKLEEIYAPKLQEFAYVTDGACSEVDILKMELNILKALKWELC",
        )
        self.assertEqual(
            alignment.sequences[89].seq,
            "DVWKNMINKEETYVRDKLYMQRHPLLQPKMRTILLDWLMEVCEVYKLYRETFYLAQDFFDRFMATQQNVVKTLLQLIGISSLFIAAKLEEIYPPKLHQFAYVTDGACTEDEILSMELIIMKALNWNLN",
        )
        self.assertEqual(
            alignment.sequences[90].seq,
            "EVWRIMLNKEKTYLRDEHFLQRHPLLQARMRAVLLDWLMEVCEVYKLHRETFYLAQDFFDRYMASQHNIIKTLLQLIGISALFIASKLEEIYPPKLHQFAYVTDGACSGDEILTMELMMMKALKWRLS",
        )
        self.assertEqual(
            alignment.sequences[91].seq,
            "EVWNNLLGKDKLYLRDTRVMERHPNLQPKMRAILLDWLMEVCEVYKLHRETFYLGQDYFDRFMATQENVLKTTLQLIGISCLFIAAKMEEIYPPKVHQFAYVTDGACTEDDILSMEIIIMKELNWSLS",
        )
        self.assertEqual(
            alignment.sequences[92].seq,
            "DVWRNMLNKDRTYLRDKNFFQKHPQLQPNMRAILLDWLMEVCEVYKLHRETFYLGQDFFDRFMATQKNVIKSRLQLIGITSLFIAAKLEEIYPPKLHQFAFITDGACTEDEITSMELIIMKDLDWCLS",
        )
        self.assertEqual(
            alignment.sequences[93].seq,
            "EVWTIMTRKEALCPRKHDCLKSHPSLGERMRAILLDWLIEVCEVYRLHRESFYLAADFVDRYLAAKENVPKTKLQLIGITSLFVAAKLEEIYPPKLHEFAYVTDGACTDDQILDQELIMLMTLNWDLT",
        )
        self.assertEqual(
            alignment.sequences[94].seq,
            "KVWSLMVKRDEIPRATRFLLGNHPDMDDEKRRILIDWMMEVCESEKLHRETFHLAVDYVDRYLESSNVECSTDNFQLVGTAALFIAAKYEEIYPPKCIDFAHLTDSAFTCDNIRTMEVLIVKYIGWSLG",
        )
        self.assertEqual(
            alignment[0],
            "GIFDYYRHREV---HFRVRKYL--HKHPE---VDV-KTRAILIDW---MVEIQETFELNHETLYNAVKLTDMYLCKTK-NVDKN------TIQKLACVAIFIAAKY-----------------------DERS--PPLVDDLIYLS--------------GD--RFSRDELLAMERELFATVGYDLG",
        )
        self.assertEqual(
            alignment[1],
            "DIFDYYRDREV---KFRIPDYM--FQQTD---LTP-SMRAILVDW---LVEVQQSFELNHETLYMAVKLIDIFSS-KV-TIKRN------KLQLIGAVALNLACKF-----------------------EERC--PPMLDDFVYVC--------------DD--AYPRQEFLKMEELVFQAVGFDIG",
        )
        self.assertEqual(
            alignment[2],
            "DIFNYLKVREA---EFPIADYM--PRQIH---LTT-WMRTLLVDW---MVEVQETFELNHETLYLAVKIVDLYLC-RE-VINKE------KLQLLGAAAFFIACKY-----------------------DERQ--PPLIEDFLYIC--------------DG--AYNHDELVRMERETLRVIKYDLG",
        )
        self.assertEqual(
            alignment[3],
            "EIFDYMREREE---KFLLPDYM--EKQSD---ISR-DMRAILVDW---MVEVQENFELNHETLYLAVKLVDHYLV-EV-VSMRD------KLQLIGSTAVLIASKF-----------------------EERC--PPCVDDFLYIC--------------DD--AYKREELIAMETSILRTLNFDIN",
        )
        self.assertEqual(
            alignment[4],
            "DIYQFYWTAEA--LNPALGHYL--SAHAE---VSP-VTRGILINW---LIEVHFKFDLMHETLYLTMDLLDRYLS-QV-PIHKN------EMQLIGLTALLLASKY-----------------------EDYW--HPRIKDLISIS--------------AE--SYTREQILGMERSMLKQLKFRLN",
        )
        self.assertEqual(
            alignment[5],
            "DIYKFYKLVEN--ESRP-HDYI--GSQPE---INE-RMRAILVDW---LIDVHTKFELSLETLYLTINIIDRFLA-VK-TVPRR------ELQLVGISAMLMASKY-----------------------EEIW--PPEVNDFVCLS--------------DR--AYTHEHILTMEKTILNKLEWTLT",
        )
        self.assertEqual(
            alignment[6],
            "DIYSFYKSVES--EWRP-RDYM--ASQPD---INE-KMRLILVEW---LIDVHVRFELNPETFYLTVNILDRFLS-VK-PVPRK------ELQLVGLSALLMSAKY-----------------------EEIW--PPQVEDLVDIA--------------DH--AYSHKQILVMEKTILSTLEWYLT",
        )
        self.assertEqual(
            alignment[7],
            "DMYSFYKEVEK--ESQP-KMYM--HIQTE---MNE-KMRAILIDW---LLEVHIKFELNLETLYLTVNIIDRFLS-VK-AVPKR------ELQLVGISALLIASKY-----------------------EEIW--PPQVNDLVYVT--------------DN--AYSSRQILVMEKAILGNLEWYLT",
        )
        self.assertEqual(
            alignment[8],
            "DIYKFYKLTED--ENRP-CDYM--DSQPE---IND-RVRAILVDW---LIEAHKRFELRPESLYLTVNIMDRFLS-EE-TVPRR------ELQLLCISSMLIACKY-----------------------EEIW--APEVNDFLTIT--------------DN--AYVRDQILLMEKVILGKLEWYLT",
        )
        self.assertEqual(
            alignment[9],
            "DIYKFYKETEE--DGCV-HDYM--GSQPD---INA-KMRSILVDW---LIEVHRKFELMPETLYLTLNIVDRFLS-VK-AVPRR------ELQLVGISSMLIASKY-----------------------EEIW--APEVNDFVCIS--------------DN--AYVSEQVLMMEKTILRKLEWYLT",
        )
        self.assertEqual(
            alignment[10],
            "DIFKFYRTVEE--EGGI-KDYI--GSQPE---INE-KMRSILIDW---LVDVHRKFELMPETLYLTINLVDRFLS-LT-MVHRR------ELQLLGLGAMLIACKY-----------------------EEIW--APEVNDFVCIS--------------DN--AYNRKQVLAMEKSILGQVEWYIT",
        )
        self.assertEqual(
            alignment[11],
            "DIYRFYKSTEG--TCLPLSSYM--SSQAE---ISE-RMRAILIDW---IIEVQYRLTLMPETLYLTVYIIDQYLS-ME-SVPRK------ELQLVGISAMLIASKY-----------------------EEIW--APLVKDLMCLC--------------DN--AFTRDQILTKEKAILDMLHWNLT",
        )
        self.assertEqual(
            alignment[12],
            "DIYTFYKIAQH--DRRP-CDYI--DTQVE---INP-KMRAILAGW---IIEVHHKFELMPETLYLTMYIIDQYLS-LQ-PVLRR------ELQLVGVSAMLIACKY-----------------------EEIW--APEVNDFILIS--------------DS--AYSREQILSMEKGILNSLEWNLT",
        )
        self.assertEqual(
            alignment[13],
            "DIYTFYKTAQH--ESRP-IDYM--GNQPE---LSP-RMRSILADW---LIESHRRFQLMPETLYLTIYIVDRYLS-LQ-PTPRR------ELQLVGVAALLIACKY-----------------------EEIW--APEVNDLIHIA--------------DG--AFNRSQILAAEKAILNSMEWNLT",
        )
        self.assertEqual(
            alignment[14],
            "ELYKFYRENEE--MSCVQPDYM--SSQGD---INE-KMRAILIDW---LIEVHHKFELMDETLFLTVNIVDRFLE-KQ-VVPRK------KLQLVGVTAMLLACKY-----------------------EEVA--VPVVEDLVLIS--------------DR--AYTKGQILEMEKLILNTLQFNMS",
        )
        self.assertEqual(
            alignment[15],
            "DLYAFYRTMER--FSCVPVDYM--MQQID---LNE-KMRAILIDW---LIEVHDKFDLINETLFLTVNLIDRFLS-KQ-NVMRK------KLQLVGLVALLLACKY-----------------------EEVS--VPVVEDLVLIS--------------DK--AYTRNDVLEMEKTMLSTLQFNIS",
        )
        self.assertEqual(
            alignment[16],
            "DLFANYRTMEV--NSCASPYYM--AQQAD---INE-RMRSILIDW---LIEVHHKFELREETLFLTVNLIDRFLE-KQ-GIVRK------KLQLVGLVAMLLACKY-----------------------EEVC--APLVEDLVLIS--------------DK--AYTRKEVLEMESMMLNTLQFNMS",
        )
        self.assertEqual(
            alignment[17],
            "EIYRFYRKTEG--ASCVPTNYM--SSQTD---INE-KMRGILIDW---LIEVHYKLELLEETLFLTVNIIDRFLA-RE-NVVRK------KLQLAGVTAMLLACKY-----------------------EEVS--VPVVEDLILIC--------------DR--AYTRADILEMERRIVNTLNFNMS",
        )
        self.assertEqual(
            alignment[18],
            "DIYCFYKKNEC--RSCVPPNYM--ENQHD---INE-RMRGILFDW---LIEVHYKFELMEETLYLTINLIDRFLAVHQ-HIARK------KLQLVGVTAMLLACKY-----------------------EEVS--VPVVDDLILIS--------------DK--AYTRTEILDMEKLMANTLQFNFC",
        )
        self.assertEqual(
            alignment[19],
            "EIFAYYREKEQ--IDKIDKDYI--KNQYH---INE-RMRAILVDW---MMAVHVRFKLLSETFFLSVNIVDRYLA-KV-MIPVT------KLQLVGITAILLACKY-----------------------EEIY--SPQIKDFVHTS--------------DD--ACTHAEVIDMERQILSTLQFHMS",
        )
        self.assertEqual(
            alignment[20],
            "DIYNYLVHHEK--KYVLDDSFI------NGGNVNS-KMRRILVDW---LIQVHLRFHLTPETLHLTIFVLDRIIV-KN-IVSKA------EFQLLGVAALFVASKF-----------------------EDIY--LPDILEYEMIT--------------DN--TFSKKQIMAMEQTILNALNFDLS",
        )
        self.assertEqual(
            alignment[21],
            "DIYKYLVHHEK--KYLLEECFM------EGGEPTP-KMRRILVDW---LVQVHVRFHLTPETLHLTVFILDRMLQ-KK-VTSKA------DLQLLGISAMFVASKF-----------------------EEVY--LPDIHDYEFIT--------------EN--TYSKKQILAMEQTILNSLNFDLS",
        )
        self.assertEqual(
            alignment[22],
            "DIYKYLTELEE--KYSIEPDHL--KKQTV---ITG-KMRATLIDW---LVEVQRQFSLVLETFHLTVGIIDRYLQVVP-NVQRN------QLQLVGVTAMFIASKY-----------------------EEIY--APDVGDFVYVT--------------DN--AYTKSDVFRCERDIMCKLGFCLA",
        )
        self.assertEqual(
            alignment[23],
            "DIYDYLYQVEL--EQPIHKDHL--AGQKE---VSH-KMRAVLIDW---INEVHLQFHLAAETFQLAVAIIDRYLQVVK-DTKRT------YLQLVGVTALFIATKY-----------------------EELF--PPAIGDFVFIT--------------DD--TYTARQIRQMELQIFKAIDCNLS",
        )
        self.assertEqual(
            alignment[24],
            "DIYCYLRSLEN--AQAVRQNYL--HG-QE---VTG-NMRAILIDW---LVQVQMKFRLLQETMFMTVGIIDRFLQ-DH-PVPKN------QLQLVGVTAMFLAAKY-----------------------EEMY--PPEIGDFTFVT--------------DH--TYTKAQIRDMEMKVLRVLKFAIG",
        )
        self.assertEqual(
            alignment[25],
            "DIYAYLRQLEE--EQSVRPKYL--QG-RE---VTG-NMRAILIDW---LIQVQMKFRLLQETMYMTVSIIDRFMQ-NS-CVPKK------MLQLVGVTAMFIASKY-----------------------EEMY--PPEIGDFAFVT--------------NN--TYTKHQIRQMEMKILRVLNFSLG",
        )
        self.assertEqual(
            alignment[26],
            "DIYNYLKQLEV--QQSVRPCYL--EG-KE---INE-RMRAILVDW---IVQVHSRFQLLQETLYMGIAIMDRFLQ-VQ-PVSRS------KLQLVGVTSLLVASKY-----------------------EEMY--TPEVADFVYIT--------------DN--AYTASQIREMEMIILRVLNFDLG",
        )
        self.assertEqual(
            alignment[27],
            "DIYQYLRQLEV--LQSINPHFL--DG-RD---ING-RMRAILVDW---LVQVHSKFRLLQETLYMCVGIMDRFLQ-VQ-PVSRK------KLQLVGITALLLASKY-----------------------EEMF--SPNIEDFVYIT--------------DN--AYTSSQIREMETLILKELKFELG",
        )
        self.assertEqual(
            alignment[28],
            "DIYLYLRQLEL--QQSVRPHYL--DG-KT---ING-RMRAILVDW---LVQVHSRFQLLQETLYMCVAVMDRFLQ-SH-PVPRK------RLQLVGVTALLLASKY-----------------------EEMY--SPDIADFVYIT--------------DN--AYNSAEVREMEITILKELNFDLG",
        )
        self.assertEqual(
            alignment[29],
            "EILQHLLIEEN--KYTI-NQYMTPEQQPD---INI-KMRAILVDW---LIDVHAKFKLKDETLYITISLIDRYLA-LA-QVTRM------RLQLVGVAALFIACKY-----------------------EEIY--PPALKDFVYIT--------------DN--AYVKSDVLEMEGLMLQALNFNIC",
        )
        self.assertEqual(
            alignment[30],
            "EIYTYLLTQEE--KYLVSNNYMNEQQQPD---LNA-RMRAILLDW---LIDVHLKFKLRDETLYVTTYLIDRFLN-FK-TTTRQ------QLQLVGVASLFIACKY-----------------------EEIY--PPDLKDFVYIT--------------DN--AYTKQDVLEMEGQILQTLDFSIT",
        )
        self.assertEqual(
            alignment[31],
            "EIFSYYYELET--RMLPDPQYL--FKQTL---LKP-RMRSILVDW---LVEMHLKFKLLPESLFLAVNVMDRFMS-VE-VVQID------KLQLLATAALFTAAKY-----------------------EEVF--SPSVKNYAYFT--------------DG--SYTPEEVVQAEKYMLTILNFDLN",
        )
        self.assertEqual(
            alignment[32],
            "EIFAFLYRREL--ETLPSHNYL--LDKTSKYYLRP-SMRTILVDW---LVEVHEKFQCYPETLFLSINLMDRFLA-KN-KVTMN------KLQLLAVTSLFIAAKF-----------------------EEVN--LPKLAEYAYIT--------------DG--AASKNDIKNAEMFMLTSLEFNIG",
        )
        self.assertEqual(
            alignment[33],
            "SIFSHLYEKEI--QMLPTHNYL--MDTQSPYHLKS-SMRALLIDW---LVEVHEKFHCLPETLFLAINLLDRFLS-QN-VVKLN------KLQLLCITCLFIACKF-----------------------EEVK--LPKITNFAYVT--------------DG--AATVEGIRKAELFVLSSLGYNIS",
        )
        self.assertEqual(
            alignment[34],
            "DIFDYLHHLEI--ITLPNKANL--YKHKN---IK--QNRDILVNW---IIKIHNKFGLLPETLYLAINIMDRFLC-EE-VVQLN------RLQLVGTSCLFIASKY-----------------------EEIY--SPSIKHFAYET--------------DG--ACSVEDIKEGERFILEKLDFQIS",
        )
        self.assertEqual(
            alignment[35],
            "DIFEYLHQLEV--ITLPKKEDL--YQHRN---IH--QNRDILVNW---LVKIHNKFGLLPETLYLAINIMDRFLG-KE-LVQLD------KLQLVGTSCLFIASKY-----------------------EEVY--SPSIKHFASET--------------DG--ACTEDEIKEGEKFILKTLKFNLN",
        )
        self.assertEqual(
            alignment[36],
            "EIFNYLHELEN--KFTPDPNYM--DFQDD---LKW-EMRAVLIDW---VVQVHARFNLFSETLYLTVNYIDRFLS-KR-RVSLS------RFQLVGAVALFIAAKY-----------------------EEIN--CPTVQEIAYMA--------------DN--AYSIDEFLKAERFMIDVLEFDLG",
        )
        self.assertEqual(
            alignment[37],
            "EIFEYMRKLED--LYKPNPYYM--DKQPE---LRW-SFRSTLIDW---IVQVHEKFQLLPETLYLCINIIDRYLC-KE-VVPVN------KFQLVGAASLFIAAKY-----------------------EEIN--CPTIKDFVYMS--------------EN--CYSRNDLLDAERTILNGLEFELG",
        )
        self.assertEqual(
            alignment[38],
            "DIFYYLRELEV--KYRPNPYYM--QNQVE---LTW-PFRRTMIDW---LVQLHFRFQLLPETLYLTINIVDRFLS-KK-TVTLN------RFQLVGVSALFIAAKF-----------------------EEIN--CPTLDDLVYML--------------EN--TYTRDDIIRAEQYMIDTLEFEIG",
        )
        self.assertEqual(
            alignment[39],
            "EIFEYIRKLDL--KCLPNPKYM--DQQKE---LTW-KMREILNEW---LVEIHSNFCLMPETLYLAVNIIDRFLS-RR-SCSLS------KFQLTGITALLIASKY-----------------------EEVM--CPSIQNFVYMT--------------DG--AFTVEDVCVAERYMLNVLNFDLS",
        )
        self.assertEqual(
            alignment[40],
            "EIMSYMRELEV--LTLPLPDYM--DRQKE---LQW-KMRGILVDW---LIEVHAKFRLLPETLFLSVNIIDRFLS-LR-VCSLP------KLQLVGITALFIAAKY-----------------------EEVM--CPSIQNFMYMA--------------DG--GYTNEEILKAEQYVLQVLGYDMS",
        )
        self.assertEqual(
            alignment[41],
            "EIFDYLRELEM--ETLPNPDYI--DHQPD---LEW-KMRGILVDW---LIEVHTRFRLLPETLFLAVNIIDRFLS-AE-VVALD------RLQLVGVAAMFIASKY-----------------------EEVL--SPHVANFSHVA--------------DE--TFSDKEILDAERHILATLEYNMS",
        )
        self.assertEqual(
            alignment[42],
            "DIFEYLNELEI--ETMPSPTYM--DRQKE---LAW-KMRGILTDW---LIEVHSRFRLLPETLFLAVNIIDRFLS-LR-VCSLN------KLQLVGIAALFIASKY-----------------------EEVM--CPSVQNFVYMA--------------DG--GYDEEEILQAERYILRVLEFNLA",
        )
        self.assertEqual(
            alignment[43],
            "EIFHYMQSLER--KLAPPPNYM--SVQQE---IDW-VTRHMLVDW---IVQVQIHFRLLPETLFLAVNLIDRFLS-IK-VVSLQ------KVQLVGLSALLIACKY-----------------------EEIH--PPSIYNFAHVV--------------QG--IFTVDEIIRAERYMLMLLDFDIS",
        )
        self.assertEqual(
            alignment[44],
            "EILSHMEKLEI--RFMPDYRHM--SAQPY---YVT-EMRASVINW---IVGVHTCINLLPESLFLSINVLDRFLS-LQ-NVPAS------KMKLCGATALFIACKY-----------------------EEIH--PPTVKDLEIVL--------------EG--EWIGEDICGMEKYMLMVLQYQLG",
        )
        self.assertEqual(
            alignment[45],
            "SVCQLFQASQA----VNKQQIF--SVQKG---LSD-TMRYILIDW---LVEVATMKDFTSLCLHLTVECVDRYLR-RR-LVPRY------KLQLLGIACMVICTRFI---------------------SKEIL----TIREAVWLT--------------DN--TYKYEDLVRVMGEIISALEGKIR",
        )
        self.assertEqual(
            alignment[46],
            "DIHGYLREMEMQNKRRPMVDYI-EKVQKI---VTP-TMRAILVDW---LVEVAVEYKLLSDTLHLSVSYIDRFLS-VN-PVSKS------RLQLLGVSSMLIAAKY-----------------------EEMD--PPGVDEFCSIT--------------DH--TYDKTEVVKMEADILKSLKFEMG",
        )
        self.assertEqual(
            alignment[47],
            "DIYMHLREAET--RKRPSTDFM-ETIQKD---VNP-SMRAILIDW---LVEVAEEYRLVPDTLYLTVNYIDRYLS-GN-EINRQ------RLQLLGVACMLIAAKY-----------------------EEIC--APQVEEFCYIT--------------DN--TYFRDEVLEMEASVLNYLKFEMT",
        )
        self.assertEqual(
            alignment[48],
            "DIYNHLRAAEA--KKQPAVDYM-ATVQKD---VNS-TMRGILVDW---LVEVSEEYRLVPETLYLTVNYIDRYLS-GN-VISRQ------KLQLLGVACMMIAAKY-----------------------EEVC--APQVEEFCYIT--------------DN--TYLKDEVLDMESAVLNYLKFEMS",
        )
        self.assertEqual(
            alignment[49],
            "DIYDNIHVAEL--QQRPLANYM-ELVQRD---IDP-DMRKILIDW---LVEVSDDYKLVPDTLYLTVNLIDRFLS-NS-YIERQ------RLQLLGVSCMLIASKY-----------------------EELS--APGVEEFCFIT--------------AN--TYTRPEVLSMEIQILNFVHFRLS",
        )
        self.assertEqual(
            alignment[50],
            "DIYSNIRVTEL--QRKPLTNYM-DKLQKD---INP-SMRGILVDW---LVEVSEEYKLVPDTLYLTVNLIDRYLS-TR-LIQKQ------KLQLLGVTCMLIASKY-----------------------EEMC--APRVEEFCFIT--------------DN--TYTKEEVLKMEREVLNLVHFQLS",
        )
        self.assertEqual(
            alignment[51],
            "DILTYGKEAEQ--RYMAKANYM--ERQSD---INH-SMRSILVDW---LVEVADEYKLKRETFFLAVNYIDRFLS-MM-SVIRC------RLQLLGAAAMFIAAKY-----------------------EEIY--PPDVAEFVYIT--------------DD--TYTMKQVLQMEQAILKTLNF-LV",
        )
        self.assertEqual(
            alignment[52],
            "DIIKYMLHRQT--KNRASHECF--DIQSQ---VNE-EMRTILIDW---FSDVVKEYNFQKETFHLAVSLVDRALS-MF-NIDKM------RFQLVGTTSMMIAVKY-----------------------EEIF--PPEIEDFALIT--------------DN--TYRVPDILLMERFLLGKFDFVVA",
        )
        self.assertEqual(
            alignment[53],
            "DILEYFRESEK--KHRPKPLYM--RRQKD---ISH-NMRSILIDW---LVEVSEEYKLDTETLYLSVFYLDRFLS-QM-AVVRS------KLQLVGTAAMYIAAKY-----------------------EEIY--PPEVGEFVFLT--------------DD--SYTKAQVLRMEQVILKILSFDLC",
        )
        self.assertEqual(
            alignment[54],
            "EIYQYLKTAES--KHRPKHGYM--RKQPD---ITN-SMRCILVDW---LVEVSEEYRLHNETLYLAAAFIDRFLS-QM-SVLRA------KLQLVGTASMFVASKY-----------------------EEIY--PPDVKEFVYIT--------------DD--TYSIKQVLRMEHLILKVLSFDLA",
        )
        self.assertEqual(
            alignment[55],
            "DIHTYLREMEV--KCKPKVGYM--KRQPD---ITN-SMRAILVDW---LVEVGEEYKLQNETLHLAVNYIDRFLS-SM-SVLRG------KLQLVGTAAMLLASKF-----------------------EEIY--PPEVAEFVYIT--------------DD--TYSKKQVLRMEHLVLKVLAFDLA",
        )
        self.assertEqual(
            alignment[56],
            "DIHRYLRECEV--KYRPKPGYM--RKQPD---ITN-CMRVILVDW---LVEVGEEYKLCSETLYLAVNYLDRFLS-CM-SVLRG------KLQLVGTAAILLAAKY-----------------------EEVY--PPEVDEFVYIT--------------DD--TYTKKQLLRMEQHLLRVLAFDMT",
        )
        self.assertEqual(
            alignment[57],
            "EIYQYLREAEI--RHRPKAHYM--KKQPD---ITE-GMRTILVDW---LVEVGEEYKLRAETLYLAVNFLDRFLS-CM-SVLRG------KLQLVGTAAMLLASKY-----------------------EEIY--PPEVDEFVYIT--------------DD--TYTKRQLLKMEHLLLKVLAFDLT",
        )
        self.assertEqual(
            alignment[58],
            "EISTNVIAQSC--KFKPNPKLI--DQQPE---MNPVETRSNIITF---LFELSVVTRVTNGIFFHSVRLYDRYCS-KR-IVLRD------QAKLVVATCLWLAAKTWGGCNHIINNVVIPTGGRFYGPNPRAR--IPRLSELVHYC--------------GDGQVFDESMFLQMERHILDTLNWNIY",
        )
        self.assertEqual(
            alignment[59],
            "DIIHHLITREK--NFLLNVHLS--NQQPE---LRW-SMRPALVNF---IVEIHNGFDLSIDTLPLSISLMDSYVS-RR-VVYCK------HIQLVACVCLWIASKF-----------------------HETEDRVPLLQELKLAC--------------KN--IYAEDLFIRMERHILDTLDWDIS",
        )
        self.assertEqual(
            alignment[60],
            "EMLHHLLSVEA--KTLPNLSLI--EQQPE---IKL-GMRPLLLDF---LMEVITILSLSRSTFPLTVNLIDRYCS-TR-IVKKQ------HYQLLGLTSLWISCKN-----------------------LDSKFKVPTLNDLRKIC--------------VD--SYYKELFVEMEKHILKSLEWVVN",
        )
        self.assertEqual(
            alignment[61],
            "DIVNTLSQLES--LTLVNPAMI--DLQPE---IQW-FMRPFLLDF---LIELHSSFKLQPTTLFLCLNIIDRYCA-KR-IVFKR------HYQLVGCTALWIASKY-----------------------EDKKSRVPTLKELTIMC--------------RN--AYDEEMFVQMEMHILSTLDWSIG",
        )
        self.assertEqual(
            alignment[62],
            "RVLNNLKLREL---LLPKFTSL-WEIQTE---VTV-DNRTILLTW---MHLLCESFELDKSVFPLSVSILDRYLC-KK-QGTKK------TLQKIGAACVLIGSKI-----------------------RTVK--PMTVSKLTYLS--------------CD--CFTNLELINQEKDILEALKWDTE",
        )
        self.assertEqual(
            alignment[63],
            "-IFYNILEIEP--RFLTSDSVFGTFQQS----LTS-HMRKLLGTW---MFSVCQEYNLEPNVVALALNLLDRLLL-IK-QVSKE------HFQKTGSACLLVASKL-----------------------RSLT--PISTSSLCYAA--------------AD--SFSRQELIDQEKELLEKLAWRTE",
        )
        self.assertEqual(
            alignment[64],
            "DMRAFYNCMEYEEALQPNYHYF-TGVQEN---ITP-FHREQAIDW---IYDVAKEENCDGDVFLLAVSLIDRFMS-VQ-NILKH------DIQMIAGVALFIASKL-----------------------KAPH--PMTASKIAYYS--------------DN--SCPIDMILQWELLIVTTLQWETE",
        )
        self.assertEqual(
            alignment[65],
            "--LENFLKVEEKHHKIPDTYF---SIQKD---ITP-PMRKIVAEW---MMEVCAEENCQEEVVLLALNYMDRFLS-SK-SVRKT------QLQILAAACLLLASKL-----------------------REPSCRALSVDLLVVYT--------------DN--SIYKDDLIKWELYVLSRLGWDLS",
        )
        self.assertEqual(
            alignment[66],
            "RVLRAMLKTEE--TCAPSVSYF-KCVQRE---IVP-SMRKIVATW---MLEVCEEQKCEEEVFPLAMNYLDRFLS-LE-PLKKS------RLQLLGATCMFVASKM-----------------------KETI--PLTAEKLCIYT--------------DN--SIRPEELLQMELLLVNKLKWNLA",
        )
        self.assertEqual(
            alignment[67],
            "RVLQNLLTIEE--RYLPQCSYF-KCVQKD---IQP-YMRRMVATW---MLEVCEEQKCEEEVFPLAMNYLDRFLA-GV-PTPKT------HLQLLGAVCMFLASKL-----------------------KETI--PLTAEKLCIYT--------------DN--SVKPQELLEWELVVLGKLKWNLA",
        )
        self.assertEqual(
            alignment[68],
            "RVLQSLLRLEE--RYVPRASYF-QCVQRE---IKP-HMRKMLAYW---MLEVCEEQRCEEEVFPLAMNYLDRYLS-CV-PTRKA------QLQLLGAVCMLLASKL-----------------------RETT--PLTIEKLCIYT--------------DH--AVSPRQLRDWEVLVLGKLKWDLA",
        )
        self.assertEqual(
            alignment[69],
            "DIFLTMREQEL-------------SRRPLFYLSPQLNERRRMLQL---LKLATSAHKLSRCALHLAVYYMDRFVD-YY-KIRPD------KLLLVAITCLHIAAQI-----------------------ENTDAFIPRYSEMNRLV--------------KN--AYTAFEYKAVERKILCFLNFELI",
        )
        self.assertEqual(
            alignment[70],
            "DSKQNGIEQSITKNIPITHRDLHYDKDYN--------LRIYCYFL---IMKLGRRLNIRQYALATAHIYLSRFLI-KA-SVREI------NLYMLVTTCVYLACKV-----------------------EEC---PQYIRTLVSEART----------LWPEFIPPDPTKVTEFEFYLLEELESYLI",
        )
        self.assertEqual(
            alignment[71],
            "----NFWTSSHCKQLLDQEDVDKVPQADSDRGITLEEFRLVKIHMSFHIWRLAQQVKVRQRVIATAVTYFRRVYT-RK-SMTEY------DPRLVAPTCLYLASKV-----------------------EES---TVQARLLVFYIKKM--------CASDEKYRFEIKDILEMEMKLLEALDYYLV",
        )
        self.assertEqual(
            alignment[72],
            "DKIWYFSNDQL-ANSPSRRCGIKGDDELQ--------YRQMTAYL---IQEMGQRLQVSQLCINTAIVYMHRFYA-FH-SFTHF------HRNSMASASLFLAAKV-----------------------EEQ---PRKLEHVIRAANKCL------PPTTEQNYAELAQELVFNENVLLQTLGFDVA",
        )
        self.assertEqual(
            alignment[73],
            "IIPWFFSREEIERNSPSRRDGIDLKTETR--------LRDSYCTF---LEILGERLKVPQVTIATAIFFCHRFFL-RQ-SHAKN------DRQTIATVCMLLAGKV-----------------------EET---PVTLEDVIIASYERIHKKDLAGAQRKEVYDQQKELVLIGEELVLSTLNFDLC",
        )
        self.assertEqual(
            alignment[74],
            "DVMSMQQHVELNKAQTMKPIDYRKMNKPG-----------VVPMY---IFECAAKLKMKPLTAACAAIVFHRFFR----EVKASD----YDEFLIAAGSLYLAGKI-----------------------KEDE--SVKIRDVINVAYCTLNRGNDPVDLNDEYWSM-RDAIVQAELLITRTLCFDLN",
        )
        self.assertEqual(
            alignment[75],
            "DIIKFLTDTTL--RVVPSSNYPTPPGSPG---EKHLTRLPSLMTF---ITRLVRYTNVYTPTLLTAACYLNKLKR----ILPRDATGLPSTIHRIFLACLILSAKF-----------------------HNDS--SPLNKHWARYT--------------DG--LFTLEDINLMERQLLQLLNWDLR",
        )
        self.assertEqual(
            alignment[76],
            "EMVQYLASTTASIIKIKKTNSMIDIALPA----------PPLTKF---INRLIKHSNVQTPTLMATSVYLAKLRS----IIPSNVYGIETTRHRIFLGCLILAAKT-----------------------LNDS--SPLNKHWAEYT--------------DG--LLILREVNTIERELLEYFDWDVT",
        )
        self.assertEqual(
            alignment[77],
            "EMIQFLATSTASIIKIRENNNPIQGCRP-----------PDLSIF---IKNVVIQSNVQTPTLMATSVYLNKLKS----VIPKNVYGINTTRHRIFLGCLILAAKT-----------------------LNDS--SPWNKHWTTYT--------------EG--LLRIREVNTIERELLEYLNWDVR",
        )
        self.assertEqual(
            alignment[78],
            "RVAESNDLTRRVATQSQRVSVFHGLSRPT----------ITIQSY---LERIFKYANCSPSCFVVAYVYLDRFTH-RQPSLPINS----FNVHRLLITSVMVAAKF--------------------------------LDDLYYNNAYY-------AKVG----GISTKEMNFLELDFLFGLGFELN",
        )
        self.assertEqual(
            alignment[79],
            "DIFD------------EKLHPLTHDQVPDNYDTHNPEHRQ-IYKF---VRTLFNAAQLTAECAIITLVYLERLLTYAELDVGPC------NWKRMVLGAILLASKV--------------------------------WDDQAVWNVDYC------QILK----DITVEDMNELERQFLELLQFNIN",
        )
        self.assertEqual(
            alignment[80],
            "DSIACFIEDER--HFVPGHDYLSRFQTRS---LDA-SAREDSVAW---ILKVQAYYNFQPLTAYLAVNYMDRFLY-AR-RLPETS---GWPMQLLAVACLSLAAKM-----------------------EEIL--VPSLFDFQVA---------------GVKYLFEAKTIKRMELLVLSVLDWRLR",
        )
        self.assertEqual(
            alignment[81],
            "DRIKEMLVREI--EFCPGTDYVKRLLSGD---LDL-SVRNQALDW---ILKVCAHYHFGHLCICLSMNYLDRFLT-SY-ELPKDK---DWAAQLLAVSCLSLASKM-----------------------EETD--VPHIVDLQVE---------------DPKFVFEAKTIKRMELLVVTTLNWRLQ",
        )
        self.assertEqual(
            alignment[82],
            "EIIMEMVEKEK--QHLPSDDYIKRLRSGD---LDLNVGRRDALNW---IWKACEVHQFGPLCFCLAMNYLDRFLS-VH-DLPSGK---GWILQLLAVACLSLAAKI-----------------------EETE--VPMLIDLQVG---------------DPQFVFEAKSVQRMELLVLNKLKWRLR",
        )
        self.assertEqual(
            alignment[83],
            "EELTSLFSKET--EYEISYNVLEK----N---QSFISSRRESVEW---ILKTTAYYSFSAQTGFLAVNYFDRFLL--F-SFNQSLNHKPWMNQLVAVTCLSLAAKV-----------------------EETD--VPLLLDLQVE---------------ESGFLFESKTIQRMEMLILSTLKWKMN",
        )
        self.assertEqual(
            alignment[84],
            "DELATLLSKEN--EFHLGFQSLIS----D---GSLMGARKEALDW---MLRVIAYYGFTATTAVLAVNYFDRFVS-GW-CFQKDK---PWMSQLAAVACLSIAAKV-----------------------EETQ--VPLLLDLQVA---------------DSRFVFEAKTIQRMELLVLSTLKWKMN",
        )
        self.assertEqual(
            alignment[85],
            "DELSTLISKQE--------PCLYDEILDD---EFLVLCREKALDW---IFKVKSHYGFNSLTALLAVNYFDRFIT-SR-KFQTDK---PWMSQLTALACLSLAAKV-----------------------EEIR--VPFLLDFQVE---------------EARYVFEAKTIQRMELLVLSTLDWRMH",
        )
        self.assertEqual(
            alignment[86],
            "TLPHSLFLVEF--QHMPSSHYFHSLKSSA---FLL-SNRNQAISS---ITQYSRKFD-DPSLTYLAVNYLDRFLS-SE-DMPQSK---PWILKLISLSCVSLSAKM-----------------------RKPD---MSVSDLPVE---------------GE--FFDAQMIERMENVILGALKWRMR",
        )
        self.assertEqual(
            alignment[87],
            "DVWRLMCHRDEQDSRLRSISML--EQHPG---LQP-RMRAILLDW---LIEVCEVYKLHRETFYLAVDYLDRYLHVAH-KVQKT------HLQLIGITCLFVAAKV-----------------------EEIY--PPKIGEFAYVT--------------DG--ACTERDILNHEKILLQALDWDIS",
        )
        self.assertEqual(
            alignment[88],
            "EVWQNMLQKEN--RYVHDKHFQ--VLHSD---LEP-QMRSILLDW---LLEVCEVYTLHRETFYLAQDFFDRFMLTQK-DVNKN------MLQLIGITSLFIASKL-----------------------EEIY--APKLQEFAYVT--------------DG--ACSEVDILKMELNILKALKWELC",
        )
        self.assertEqual(
            alignment[89],
            "DVWKNMINKEE--TYVRDKLYM--QRHPL---LQP-KMRTILLDW---LMEVCEVYKLYRETFYLAQDFFDRFMATQQ-NVVKT------LLQLIGISSLFIAAKL-----------------------EEIY--PPKLHQFAYVT--------------DG--ACTEDEILSMELIIMKALNWNLN",
        )
        self.assertEqual(
            alignment[90],
            "EVWRIMLNKEK--TYLRDEHFL--QRHPL---LQA-RMRAVLLDW---LMEVCEVYKLHRETFYLAQDFFDRYMASQH-NIIKT------LLQLIGISALFIASKL-----------------------EEIY--PPKLHQFAYVT--------------DG--ACSGDEILTMELMMMKALKWRLS",
        )
        self.assertEqual(
            alignment[91],
            "EVWNNLLGKDK--LYLRDTRVM--ERHPN---LQP-KMRAILLDW---LMEVCEVYKLHRETFYLGQDYFDRFMATQE-NVLKT------TLQLIGISCLFIAAKM-----------------------EEIY--PPKVHQFAYVT--------------DG--ACTEDDILSMEIIIMKELNWSLS",
        )
        self.assertEqual(
            alignment[92],
            "DVWRNMLNKDR--TYLRDKNFF--QKHPQ---LQP-NMRAILLDW---LMEVCEVYKLHRETFYLGQDFFDRFMATQK-NVIKS------RLQLIGITSLFIAAKL-----------------------EEIY--PPKLHQFAFIT--------------DG--ACTEDEITSMELIIMKDLDWCLS",
        )
        self.assertEqual(
            alignment[93],
            "EVWTIMTRKEA--LCPRKHDCL--KSHPS---LGE-RMRAILLDW---LIEVCEVYRLHRESFYLAADFVDRYLAAKE-NVPKT------KLQLIGITSLFVAAKL-----------------------EEIY--PPKLHEFAYVT--------------DG--ACTDDQILDQELIMLMTLNWDLT",
        )
        self.assertEqual(
            alignment[94],
            "KVWSLMVKRDE--IPRATRFLL--GNHPD---MDD-EKRRILIDW---MMEVCESEKLHRETFHLAVDYVDRYLESSNVECSTD------NFQLVGTAALFIAAKY-----------------------EEIY--PPKCIDFAHLT--------------DS--AFTCDNIRTMEVLIVKYIGWSLG",
        )
        self.assertEqual(
            alignment.column_annotations["consensus secondary structure"],
            "HHHHHHHHHHC..HTS-STTCT.TTCTSS...S-H.HHHHHHHHH...HHHHHHHTT--TTHHHHHHHHHHHHHH.HS.---CC......CHHHHHHHHHHHHHHH.......................HSSS..---HHHHHHHT..............TT..SS-HHHHHHHHHHHHHHTTT---",
        )
        self.assertEqual(
            alignment.column_annotations["consensus sequence"],
            "-IapahpptEt..phhsp.sYh..ppps-...ls..pMRsILlDW...Ll-VppcacLhtETLaLulshlDRFLu.tp.sls+s......cLQLlGlsulhlAuKa.......................EElh..sPplp-hshlo..............Ds..saopcpllpMEphlLpsLpasls",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
[[0, 1, 2, 4, 11, 11, 11, 11, 12, 13, 14, 15, 16, 18, 19, 19, 19, 21, 22, 23, 24, 24, 26, 27, 27, 28, 30, 31, 32, 36, 36, 45, 46, 63, 64, 65, 66, 66, 71, 71, 71, 71, 71, 87, 87, 87, 87, 90, 91, 91, 92, 93, 94, 101, 102, 102, 102, 102, 102, 102, 102, 102, 104, 104, 106, 107, 124, 125, 127],
 [0, 1, 2, 4, 11, 11, 11, 11, 12, 13, 14, 15, 16, 18, 19, 19, 19, 21, 22, 23, 24, 24, 26, 27, 27, 28, 30, 31, 32, 36, 36, 45, 46, 63, 63, 64, 65, 65, 70, 70, 70, 70, 70, 86, 86, 86, 86, 89, 90, 90, 91, 92, 93, 100, 101, 101, 101, 101, 101, 101, 101, 101, 103, 103, 105, 106, 123, 124, 126],
 [0, 1, 2, 4, 11, 11, 11, 11, 12, 13, 14, 15, 16, 18, 19, 19, 19, 21, 22, 23, 24, 24, 26, 27, 27, 28, 30, 31, 32, 36, 36, 45, 46, 63, 63, 64, 65, 65, 70, 70, 70, 70, 70, 86, 86, 86, 86, 89, 90, 90, 91, 92, 93, 100, 101, 101, 101, 101, 101, 101, 101, 101, 103, 103, 105, 106, 123, 124, 126],
 [0, 1, 2, 4, 11, 11, 11, 11, 12, 13, 14, 15, 16, 18, 19, 19, 19, 21, 22, 23, 24, 24, 26, 27, 27, 28, 30, 31, 32, 36, 36, 45, 46, 63, 63, 64, 65, 65, 70, 70, 70, 70, 70, 86, 86, 86, 86, 89, 90, 90, 91, 92, 93, 100, 101, 101, 101, 101, 101, 101, 101, 101, 103, 103, 105, 106, 123, 124, 126],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 64, 65, 66, 66, 71, 71, 71, 71, 71, 87, 87, 87, 87, 90, 91, 91, 92, 93, 94, 101, 102, 102, 102, 102, 102, 102, 102, 102, 104, 104, 106, 107, 124, 125, 127],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 15, 16, 18, 19, 19, 19, 21, 22, 23, 24, 24, 26, 27, 27, 28, 30, 31, 32, 36, 36, 45, 46, 63, 63, 64, 65, 65, 70, 70, 70, 70, 70, 86, 86, 86, 86, 89, 90, 90, 91, 92, 93, 100, 101, 101, 101, 101, 101, 101, 101, 101, 103, 103, 105, 106, 123, 124, 126],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 15, 16, 18, 19, 19, 19, 21, 22, 23, 24, 24, 26, 27, 27, 28, 30, 31, 32, 36, 36, 45, 46, 63, 63, 64, 65, 65, 70, 70, 70, 70, 70, 86, 86, 86, 86, 89, 90, 90, 91, 92, 93, 100, 101, 101, 101, 101, 101, 101, 101, 101, 103, 103, 105, 106, 123, 124, 126],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 15, 16, 18, 19, 19, 19, 21, 22, 23, 24, 24, 26, 27, 27, 28, 30, 31, 32, 36, 36, 45, 46, 63, 63, 64, 65, 65, 70, 70, 70, 70, 70, 86, 86, 86, 86, 89, 90, 90, 91, 92, 93, 100, 101, 101, 101, 101, 101, 101, 101, 101, 103, 103, 105, 106, 123, 124, 126],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 15, 16, 18, 19, 19, 19, 21, 22, 23, 24, 24, 26, 27, 27, 28, 30, 31, 32, 36, 36, 45, 46, 63, 63, 64, 65, 65, 70, 70, 70, 70, 70, 86, 86, 86, 86, 89, 90, 90, 91, 92, 93, 100, 101, 101, 101, 101, 101, 101, 101, 101, 103, 103, 105, 106, 123, 124, 126],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 15, 16, 18, 19, 19, 19, 21, 22, 23, 24, 24, 26, 27, 27, 28, 30, 31, 32, 36, 36, 45, 46, 63, 63, 64, 65, 65, 70, 70, 70, 70, 70, 86, 86, 86, 86, 89, 90, 90, 91, 92, 93, 100, 101, 101, 101, 101, 101, 101, 101, 101, 103, 103, 105, 106, 123, 124, 126],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 15, 16, 18, 19, 19, 19, 21, 22, 23, 24, 24, 26, 27, 27, 28, 30, 31, 32, 36, 36, 45, 46, 63, 63, 64, 65, 65, 70, 70, 70, 70, 70, 86, 86, 86, 86, 89, 90, 90, 91, 92, 93, 100, 101, 101, 101, 101, 101, 101, 101, 101, 103, 103, 105, 106, 123, 124, 126],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 64, 65, 66, 66, 71, 71, 71, 71, 71, 87, 87, 87, 87, 90, 91, 91, 92, 93, 94, 101, 102, 102, 102, 102, 102, 102, 102, 102, 104, 104, 106, 107, 124, 125, 127],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 15, 16, 18, 19, 19, 19, 21, 22, 23, 24, 24, 26, 27, 27, 28, 30, 31, 32, 36, 36, 45, 46, 63, 63, 64, 65, 65, 70, 70, 70, 70, 70, 86, 86, 86, 86, 89, 90, 90, 91, 92, 93, 100, 101, 101, 101, 101, 101, 101, 101, 101, 103, 103, 105, 106, 123, 124, 126],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 15, 16, 18, 19, 19, 19, 21, 22, 23, 24, 24, 26, 27, 27, 28, 30, 31, 32, 36, 36, 45, 46, 63, 63, 64, 65, 65, 70, 70, 70, 70, 70, 86, 86, 86, 86, 89, 90, 90, 91, 92, 93, 100, 101, 101, 101, 101, 101, 101, 101, 101, 103, 103, 105, 106, 123, 124, 126],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 64, 65, 66, 66, 71, 71, 71, 71, 71, 87, 87, 87, 87, 90, 91, 91, 92, 93, 94, 101, 102, 102, 102, 102, 102, 102, 102, 102, 104, 104, 106, 107, 124, 125, 127],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 64, 65, 66, 66, 71, 71, 71, 71, 71, 87, 87, 87, 87, 90, 91, 91, 92, 93, 94, 101, 102, 102, 102, 102, 102, 102, 102, 102, 104, 104, 106, 107, 124, 125, 127],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 64, 65, 66, 66, 71, 71, 71, 71, 71, 87, 87, 87, 87, 90, 91, 91, 92, 93, 94, 101, 102, 102, 102, 102, 102, 102, 102, 102, 104, 104, 106, 107, 124, 125, 127],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 64, 65, 66, 66, 71, 71, 71, 71, 71, 87, 87, 87, 87, 90, 91, 91, 92, 93, 94, 101, 102, 102, 102, 102, 102, 102, 102, 102, 104, 104, 106, 107, 124, 125, 127],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 65, 66, 67, 67, 72, 72, 72, 72, 72, 88, 88, 88, 88, 91, 92, 92, 93, 94, 95, 102, 103, 103, 103, 103, 103, 103, 103, 103, 105, 105, 107, 108, 125, 126, 128],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 64, 65, 66, 66, 71, 71, 71, 71, 71, 87, 87, 87, 87, 90, 91, 91, 92, 93, 94, 101, 102, 102, 102, 102, 102, 102, 102, 102, 104, 104, 106, 107, 124, 125, 127],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 20, 20, 20, 21, 24, 26, 27, 27, 28, 30, 31, 32, 36, 36, 45, 46, 63, 63, 64, 65, 65, 70, 70, 70, 70, 70, 86, 86, 86, 86, 89, 90, 90, 91, 92, 93, 100, 101, 101, 101, 101, 101, 101, 101, 101, 103, 103, 105, 106, 123, 124, 126],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 20, 20, 20, 21, 24, 26, 27, 27, 28, 30, 31, 32, 36, 36, 45, 46, 63, 63, 64, 65, 65, 70, 70, 70, 70, 70, 86, 86, 86, 86, 89, 90, 90, 91, 92, 93, 100, 101, 101, 101, 101, 101, 101, 101, 101, 103, 103, 105, 106, 123, 124, 126],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 65, 66, 67, 67, 72, 72, 72, 72, 72, 88, 88, 88, 88, 91, 92, 92, 93, 94, 95, 102, 103, 103, 103, 103, 103, 103, 103, 103, 105, 105, 107, 108, 125, 126, 128],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 65, 66, 67, 67, 72, 72, 72, 72, 72, 88, 88, 88, 88, 91, 92, 92, 93, 94, 95, 102, 103, 103, 103, 103, 103, 103, 103, 103, 105, 105, 107, 108, 125, 126, 128],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 22, 23, 24, 24, 26, 27, 27, 28, 30, 31, 32, 36, 36, 45, 46, 63, 63, 64, 65, 65, 70, 70, 70, 70, 70, 86, 86, 86, 86, 89, 90, 90, 91, 92, 93, 100, 101, 101, 101, 101, 101, 101, 101, 101, 103, 103, 105, 106, 123, 124, 126],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 22, 23, 24, 24, 26, 27, 27, 28, 30, 31, 32, 36, 36, 45, 46, 63, 63, 64, 65, 65, 70, 70, 70, 70, 70, 86, 86, 86, 86, 89, 90, 90, 91, 92, 93, 100, 101, 101, 101, 101, 101, 101, 101, 101, 103, 103, 105, 106, 123, 124, 126],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 22, 23, 24, 24, 26, 27, 27, 28, 30, 31, 32, 36, 36, 45, 46, 63, 63, 64, 65, 65, 70, 70, 70, 70, 70, 86, 86, 86, 86, 89, 90, 90, 91, 92, 93, 100, 101, 101, 101, 101, 101, 101, 101, 101, 103, 103, 105, 106, 123, 124, 126],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 22, 23, 24, 24, 26, 27, 27, 28, 30, 31, 32, 36, 36, 45, 46, 63, 63, 64, 65, 65, 70, 70, 70, 70, 70, 86, 86, 86, 86, 89, 90, 90, 91, 92, 93, 100, 101, 101, 101, 101, 101, 101, 101, 101, 103, 103, 105, 106, 123, 124, 126],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 22, 23, 24, 24, 26, 27, 27, 28, 30, 31, 32, 36, 36, 45, 46, 63, 63, 64, 65, 65, 70, 70, 70, 70, 70, 86, 86, 86, 86, 89, 90, 90, 91, 92, 93, 100, 101, 101, 101, 101, 101, 101, 101, 101, 103, 103, 105, 106, 123, 124, 126],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 15, 16, 18, 19, 20, 21, 23, 24, 25, 26, 26, 28, 29, 29, 30, 32, 33, 34, 38, 38, 47, 48, 65, 65, 66, 67, 67, 72, 72, 72, 72, 72, 88, 88, 88, 88, 91, 92, 92, 93, 94, 95, 102, 103, 103, 103, 103, 103, 103, 103, 103, 105, 105, 107, 108, 125, 126, 128],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22, 24, 25, 26, 27, 27, 29, 30, 30, 31, 33, 34, 35, 39, 39, 48, 49, 66, 66, 67, 68, 68, 73, 73, 73, 73, 73, 89, 89, 89, 89, 92, 93, 93, 94, 95, 96, 103, 104, 104, 104, 104, 104, 104, 104, 104, 106, 106, 108, 109, 126, 127, 129],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 64, 65, 66, 66, 71, 71, 71, 71, 71, 87, 87, 87, 87, 90, 91, 91, 92, 93, 94, 101, 102, 102, 102, 102, 102, 102, 102, 102, 104, 104, 106, 107, 124, 125, 127],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 28, 30, 31, 31, 32, 34, 35, 36, 40, 40, 49, 50, 67, 67, 68, 69, 69, 74, 74, 74, 74, 74, 90, 90, 90, 90, 93, 94, 94, 95, 96, 97, 104, 105, 105, 105, 105, 105, 105, 105, 105, 107, 107, 109, 110, 127, 128, 130],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 28, 30, 31, 31, 32, 34, 35, 36, 40, 40, 49, 50, 67, 67, 68, 69, 69, 74, 74, 74, 74, 74, 90, 90, 90, 90, 93, 94, 94, 95, 96, 97, 104, 105, 105, 105, 105, 105, 105, 105, 105, 107, 107, 109, 110, 127, 128, 130],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 27, 27, 28, 30, 31, 32, 36, 36, 45, 46, 63, 63, 64, 65, 65, 70, 70, 70, 70, 70, 86, 86, 86, 86, 89, 90, 90, 91, 92, 93, 100, 101, 101, 101, 101, 101, 101, 101, 101, 103, 103, 105, 106, 123, 124, 126],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 27, 27, 28, 30, 31, 32, 36, 36, 45, 46, 63, 63, 64, 65, 65, 70, 70, 70, 70, 70, 86, 86, 86, 86, 89, 90, 90, 91, 92, 93, 100, 101, 101, 101, 101, 101, 101, 101, 101, 103, 103, 105, 106, 123, 124, 126],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 64, 65, 66, 66, 71, 71, 71, 71, 71, 87, 87, 87, 87, 90, 91, 91, 92, 93, 94, 101, 102, 102, 102, 102, 102, 102, 102, 102, 104, 104, 106, 107, 124, 125, 127],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 64, 65, 66, 66, 71, 71, 71, 71, 71, 87, 87, 87, 87, 90, 91, 91, 92, 93, 94, 101, 102, 102, 102, 102, 102, 102, 102, 102, 104, 104, 106, 107, 124, 125, 127],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 64, 65, 66, 66, 71, 71, 71, 71, 71, 87, 87, 87, 87, 90, 91, 91, 92, 93, 94, 101, 102, 102, 102, 102, 102, 102, 102, 102, 104, 104, 106, 107, 124, 125, 127],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 64, 65, 66, 66, 71, 71, 71, 71, 71, 87, 87, 87, 87, 90, 91, 91, 92, 93, 94, 101, 102, 102, 102, 102, 102, 102, 102, 102, 104, 104, 106, 107, 124, 125, 127],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 64, 65, 66, 66, 71, 71, 71, 71, 71, 87, 87, 87, 87, 90, 91, 91, 92, 93, 94, 101, 102, 102, 102, 102, 102, 102, 102, 102, 104, 104, 106, 107, 124, 125, 127],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 64, 65, 66, 66, 71, 71, 71, 71, 71, 87, 87, 87, 87, 90, 91, 91, 92, 93, 94, 101, 102, 102, 102, 102, 102, 102, 102, 102, 104, 104, 106, 107, 124, 125, 127],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 64, 65, 66, 66, 71, 71, 71, 71, 71, 87, 87, 87, 87, 90, 91, 91, 92, 93, 94, 101, 102, 102, 102, 102, 102, 102, 102, 102, 104, 104, 106, 107, 124, 125, 127],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 64, 65, 66, 66, 71, 71, 71, 71, 71, 87, 87, 87, 87, 90, 91, 91, 92, 93, 94, 101, 102, 102, 102, 102, 102, 102, 102, 102, 104, 104, 106, 107, 124, 125, 127],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 64, 65, 66, 66, 71, 71, 71, 71, 71, 87, 87, 87, 87, 90, 91, 91, 92, 93, 94, 101, 102, 102, 102, 102, 102, 102, 102, 102, 104, 104, 106, 107, 124, 125, 127],
 [0, 1, 2, 4, 11, 11, 11, 11, 11, 12, 13, 14, 15, 17, 18, 18, 18, 20, 21, 22, 23, 23, 25, 26, 26, 27, 29, 30, 31, 35, 35, 44, 45, 62, 62, 63, 64, 64, 69, 69, 69, 69, 69, 85, 86, 86, 87, 90, 91, 91, 91, 91, 92, 99, 100, 100, 100, 100, 100, 100, 100, 100, 102, 102, 104, 105, 122, 123, 125],
 [0, 1, 2, 4, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21, 22, 22, 23, 25, 26, 27, 28, 28, 30, 31, 31, 32, 34, 35, 36, 40, 40, 49, 50, 67, 67, 68, 69, 69, 74, 74, 74, 74, 74, 90, 90, 90, 90, 93, 94, 94, 95, 96, 97, 104, 105, 105, 105, 105, 105, 105, 105, 105, 107, 107, 109, 110, 127, 128, 130],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 21, 23, 24, 25, 26, 26, 28, 29, 29, 30, 32, 33, 34, 38, 38, 47, 48, 65, 65, 66, 67, 67, 72, 72, 72, 72, 72, 88, 88, 88, 88, 91, 92, 92, 93, 94, 95, 102, 103, 103, 103, 103, 103, 103, 103, 103, 105, 105, 107, 108, 125, 126, 128],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 21, 23, 24, 25, 26, 26, 28, 29, 29, 30, 32, 33, 34, 38, 38, 47, 48, 65, 65, 66, 67, 67, 72, 72, 72, 72, 72, 88, 88, 88, 88, 91, 92, 92, 93, 94, 95, 102, 103, 103, 103, 103, 103, 103, 103, 103, 105, 105, 107, 108, 125, 126, 128],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 21, 23, 24, 25, 26, 26, 28, 29, 29, 30, 32, 33, 34, 38, 38, 47, 48, 65, 65, 66, 67, 67, 72, 72, 72, 72, 72, 88, 88, 88, 88, 91, 92, 92, 93, 94, 95, 102, 103, 103, 103, 103, 103, 103, 103, 103, 105, 105, 107, 108, 125, 126, 128],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 21, 23, 24, 25, 26, 26, 28, 29, 29, 30, 32, 33, 34, 38, 38, 47, 48, 65, 65, 66, 67, 67, 72, 72, 72, 72, 72, 88, 88, 88, 88, 91, 92, 92, 93, 94, 95, 102, 103, 103, 103, 103, 103, 103, 103, 103, 105, 105, 107, 108, 125, 126, 128],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 64, 65, 66, 66, 71, 71, 71, 71, 71, 87, 87, 87, 87, 90, 91, 91, 92, 93, 94, 101, 102, 102, 102, 102, 102, 102, 102, 102, 104, 104, 106, 107, 124, 124, 126],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 64, 65, 66, 66, 71, 71, 71, 71, 71, 87, 87, 87, 87, 90, 91, 91, 92, 93, 94, 101, 102, 102, 102, 102, 102, 102, 102, 102, 104, 104, 106, 107, 124, 125, 127],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 64, 65, 66, 66, 71, 71, 71, 71, 71, 87, 87, 87, 87, 90, 91, 91, 92, 93, 94, 101, 102, 102, 102, 102, 102, 102, 102, 102, 104, 104, 106, 107, 124, 125, 127],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 64, 65, 66, 66, 71, 71, 71, 71, 71, 87, 87, 87, 87, 90, 91, 91, 92, 93, 94, 101, 102, 102, 102, 102, 102, 102, 102, 102, 104, 104, 106, 107, 124, 125, 127],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 64, 65, 66, 66, 71, 71, 71, 71, 71, 87, 87, 87, 87, 90, 91, 91, 92, 93, 94, 101, 102, 102, 102, 102, 102, 102, 102, 102, 104, 104, 106, 107, 124, 125, 127],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 64, 65, 66, 66, 71, 71, 71, 71, 71, 87, 87, 87, 87, 90, 91, 91, 92, 93, 94, 101, 102, 102, 102, 102, 102, 102, 102, 102, 104, 104, 106, 107, 124, 125, 127],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 64, 65, 66, 66, 71, 71, 71, 71, 71, 87, 87, 87, 87, 90, 91, 91, 92, 93, 94, 101, 102, 102, 102, 102, 102, 102, 102, 102, 104, 104, 106, 107, 124, 125, 127],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 29, 30, 32, 33, 34, 38, 38, 47, 48, 65, 65, 66, 67, 67, 72, 72, 72, 72, 72, 88, 89, 110, 111, 114, 115, 115, 116, 117, 118, 125, 126, 126, 126, 126, 126, 126, 126, 126, 128, 130, 132, 133, 150, 151, 153],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 64, 65, 66, 66, 71, 71, 71, 71, 71, 87, 87, 87, 87, 90, 91, 93, 94, 95, 96, 103, 104, 104, 104, 104, 104, 104, 104, 104, 106, 106, 108, 109, 126, 127, 129],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 64, 65, 66, 66, 71, 71, 71, 71, 71, 87, 87, 87, 87, 90, 91, 93, 94, 95, 96, 103, 104, 104, 104, 104, 104, 104, 104, 104, 106, 106, 108, 109, 126, 127, 129],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 64, 65, 66, 66, 71, 71, 71, 71, 71, 87, 87, 87, 87, 90, 91, 93, 94, 95, 96, 103, 104, 104, 104, 104, 104, 104, 104, 104, 106, 106, 108, 109, 126, 127, 129],
 [0, 1, 2, 4, 11, 11, 11, 11, 12, 13, 14, 15, 16, 18, 19, 19, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 64, 65, 66, 66, 71, 71, 71, 71, 71, 87, 87, 87, 87, 90, 91, 91, 92, 93, 94, 101, 102, 102, 102, 102, 102, 102, 102, 102, 104, 104, 106, 107, 124, 125, 127],
 [0, 0, 1, 3, 10, 10, 10, 11, 12, 13, 14, 15, 16, 18, 19, 20, 21, 23, 24, 25, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 64, 65, 66, 66, 71, 71, 71, 71, 71, 87, 87, 87, 87, 90, 91, 91, 92, 93, 94, 101, 102, 102, 102, 102, 102, 102, 102, 102, 104, 104, 106, 107, 124, 125, 127],
 [0, 1, 2, 4, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21, 22, 22, 23, 25, 26, 27, 28, 28, 30, 31, 31, 32, 34, 35, 36, 40, 40, 49, 50, 67, 67, 68, 69, 69, 74, 74, 74, 74, 74, 90, 90, 90, 90, 93, 94, 94, 95, 96, 97, 104, 105, 105, 105, 105, 105, 105, 105, 105, 107, 107, 109, 110, 127, 128, 130],
 [0, 0, 0, 2, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 19, 19, 19, 21, 22, 23, 24, 24, 26, 27, 27, 28, 30, 31, 32, 36, 36, 45, 46, 63, 63, 64, 65, 65, 70, 70, 70, 70, 70, 86, 86, 86, 86, 89, 90, 92, 93, 94, 95, 102, 103, 103, 103, 103, 103, 103, 103, 103, 105, 105, 107, 108, 125, 126, 128],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 21, 23, 24, 25, 26, 26, 28, 29, 29, 30, 32, 33, 34, 38, 38, 47, 48, 65, 65, 66, 67, 67, 72, 72, 72, 72, 72, 88, 88, 88, 88, 91, 92, 92, 93, 94, 95, 102, 103, 103, 103, 103, 103, 103, 103, 103, 105, 105, 107, 108, 125, 126, 128],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 21, 23, 24, 25, 26, 26, 28, 29, 29, 30, 32, 33, 34, 38, 38, 47, 48, 65, 65, 66, 67, 67, 72, 72, 72, 72, 72, 88, 88, 88, 88, 91, 92, 92, 93, 94, 95, 102, 103, 103, 103, 103, 103, 103, 103, 103, 105, 105, 107, 108, 125, 126, 128],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 21, 23, 24, 25, 26, 26, 28, 29, 29, 30, 32, 33, 34, 38, 38, 47, 48, 65, 65, 66, 67, 67, 72, 72, 72, 72, 72, 88, 88, 88, 88, 91, 92, 92, 93, 94, 95, 102, 103, 103, 103, 103, 103, 103, 103, 103, 105, 105, 107, 108, 125, 126, 128],
 [0, 1, 2, 4, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 13, 14, 15, 16, 19, 21, 22, 23, 24, 26, 27, 28, 32, 32, 41, 42, 59, 59, 60, 61, 61, 66, 66, 66, 66, 66, 82, 82, 82, 82, 85, 86, 88, 89, 90, 91, 98, 99, 99, 99, 99, 99, 99, 99, 99, 101, 101, 103, 104, 121, 122, 124],
 [0, 1, 2, 4, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21, 22, 23, 24, 26, 27, 28, 29, 29, 29, 29, 29, 29, 31, 32, 33, 37, 37, 46, 47, 64, 64, 65, 66, 66, 71, 71, 71, 71, 71, 87, 87, 87, 87, 90, 90, 90, 91, 92, 93, 100, 101, 103, 103, 103, 103, 103, 103, 105, 107, 109, 111, 112, 129, 130, 132],
 [0, 0, 0, 0, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 22, 23, 24, 25, 28, 30, 31, 32, 33, 35, 36, 37, 41, 44, 53, 54, 71, 71, 72, 73, 73, 78, 78, 78, 78, 78, 94, 94, 94, 94, 97, 97, 97, 98, 99, 100, 107, 108, 110, 111, 111, 111, 111, 112, 114, 116, 118, 120, 121, 138, 139, 141],
 [0, 1, 2, 4, 11, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 22, 23, 25, 26, 27, 28, 28, 28, 28, 28, 28, 30, 31, 32, 36, 36, 45, 46, 63, 63, 64, 65, 65, 70, 70, 70, 70, 70, 86, 86, 86, 86, 89, 89, 89, 90, 91, 92, 99, 100, 102, 103, 104, 104, 105, 106, 108, 110, 112, 114, 115, 132, 133, 135],
 [0, 1, 2, 4, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21, 22, 23, 24, 26, 27, 28, 29, 29, 29, 29, 29, 29, 31, 32, 33, 37, 37, 46, 47, 64, 64, 65, 66, 66, 71, 71, 71, 71, 71, 87, 87, 87, 87, 90, 90, 90, 91, 92, 93, 100, 101, 103, 104, 105, 111, 112, 113, 115, 117, 119, 121, 122, 139, 140, 142],
 [0, 1, 2, 4, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21, 22, 23, 24, 26, 27, 28, 29, 29, 29, 29, 29, 29, 29, 29, 30, 34, 34, 43, 44, 61, 61, 61, 61, 61, 66, 67, 67, 67, 68, 84, 84, 84, 84, 87, 88, 88, 89, 90, 91, 98, 99, 101, 102, 103, 109, 110, 111, 113, 115, 117, 119, 119, 136, 137, 139],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22, 24, 25, 26, 27, 27, 29, 30, 31, 32, 34, 35, 36, 40, 40, 49, 50, 67, 67, 67, 67, 67, 72, 73, 76, 77, 78, 94, 94, 94, 94, 97, 98, 98, 99, 100, 101, 108, 109, 109, 109, 109, 109, 109, 109, 109, 111, 111, 113, 114, 131, 132, 134],
 [0, 1, 2, 4, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21, 22, 23, 24, 26, 27, 28, 29, 29, 29, 29, 29, 29, 29, 30, 31, 35, 35, 44, 45, 62, 62, 62, 62, 62, 67, 68, 71, 72, 73, 89, 89, 89, 89, 92, 93, 93, 94, 95, 96, 103, 104, 104, 104, 104, 104, 104, 104, 104, 106, 106, 108, 109, 126, 127, 129],
 [0, 1, 2, 4, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21, 22, 23, 24, 26, 27, 28, 28, 28, 28, 28, 28, 28, 28, 29, 30, 34, 34, 43, 44, 61, 61, 61, 61, 61, 66, 67, 70, 71, 72, 88, 88, 88, 88, 91, 92, 92, 93, 94, 95, 102, 103, 103, 103, 103, 103, 103, 103, 103, 105, 105, 107, 108, 125, 126, 128],
 [0, 1, 2, 4, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21, 22, 23, 24, 26, 27, 28, 29, 29, 29, 29, 29, 29, 29, 30, 31, 35, 35, 44, 45, 62, 62, 63, 64, 65, 70, 71, 71, 71, 72, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 95, 96, 98, 99, 99, 99, 100, 101, 103, 103, 103, 105, 106, 123, 124, 126],
 [0, 1, 2, 4, 4, 4, 4, 4, 4, 4, 5, 6, 7, 9, 10, 11, 12, 14, 15, 16, 17, 20, 22, 23, 24, 25, 27, 28, 28, 32, 32, 41, 42, 59, 60, 61, 62, 63, 68, 68, 68, 68, 68, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 91, 92, 94, 95, 96, 96, 97, 98, 100, 100, 100, 102, 103, 120, 121, 123],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22, 24, 25, 26, 27, 27, 29, 30, 30, 31, 33, 34, 35, 39, 39, 48, 49, 66, 66, 67, 68, 68, 73, 74, 74, 75, 76, 92, 92, 92, 92, 95, 96, 96, 97, 98, 99, 106, 106, 106, 106, 106, 106, 106, 106, 106, 108, 110, 112, 113, 130, 131, 133],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22, 24, 25, 26, 27, 27, 29, 30, 30, 31, 33, 34, 35, 39, 39, 48, 49, 66, 66, 67, 68, 68, 73, 74, 74, 75, 76, 92, 92, 92, 92, 95, 96, 96, 97, 98, 99, 106, 106, 106, 106, 106, 106, 106, 106, 106, 108, 110, 112, 113, 130, 131, 133],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22, 24, 25, 26, 27, 27, 29, 30, 31, 32, 34, 35, 36, 40, 40, 49, 50, 67, 67, 68, 69, 69, 74, 75, 75, 76, 77, 93, 93, 93, 93, 96, 97, 97, 98, 99, 100, 107, 107, 107, 107, 107, 107, 107, 107, 107, 109, 111, 113, 114, 131, 132, 134],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22, 22, 22, 22, 23, 23, 25, 26, 27, 28, 30, 31, 32, 36, 36, 45, 46, 63, 63, 63, 64, 64, 69, 70, 73, 74, 75, 91, 91, 91, 91, 94, 95, 95, 96, 97, 98, 105, 105, 105, 105, 105, 105, 105, 105, 105, 107, 109, 111, 112, 129, 130, 132],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22, 22, 22, 22, 23, 23, 25, 26, 27, 28, 30, 31, 32, 36, 36, 45, 46, 63, 63, 64, 65, 65, 70, 71, 71, 72, 73, 89, 89, 89, 89, 92, 93, 93, 94, 95, 96, 103, 103, 103, 103, 103, 103, 103, 103, 103, 105, 107, 109, 110, 127, 128, 130],
 [0, 1, 2, 4, 11, 11, 11, 11, 11, 11, 11, 11, 11, 13, 14, 15, 16, 18, 19, 20, 21, 21, 23, 24, 25, 26, 28, 29, 30, 34, 34, 43, 44, 61, 61, 62, 63, 63, 68, 69, 69, 70, 71, 87, 87, 87, 87, 90, 91, 91, 92, 93, 94, 101, 101, 101, 101, 101, 101, 101, 101, 101, 103, 105, 107, 108, 125, 126, 128],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22, 24, 25, 26, 27, 27, 29, 30, 30, 31, 33, 34, 35, 39, 39, 48, 48, 65, 65, 66, 67, 67, 72, 73, 73, 74, 75, 91, 91, 91, 91, 94, 95, 95, 95, 96, 97, 104, 104, 104, 104, 104, 104, 104, 104, 104, 106, 106, 108, 109, 126, 127, 129],
 [0, 1, 2, 4, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21, 22, 22, 22, 24, 25, 26, 27, 27, 29, 30, 30, 31, 33, 34, 35, 39, 39, 48, 49, 66, 67, 68, 69, 69, 74, 74, 74, 74, 74, 90, 90, 90, 90, 93, 94, 94, 95, 96, 97, 104, 105, 105, 105, 105, 105, 105, 105, 105, 107, 107, 109, 110, 127, 128, 130],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 65, 66, 67, 67, 72, 72, 72, 72, 72, 88, 88, 88, 88, 91, 92, 92, 93, 94, 95, 102, 103, 103, 103, 103, 103, 103, 103, 103, 105, 105, 107, 108, 125, 126, 128],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 65, 66, 67, 67, 72, 72, 72, 72, 72, 88, 88, 88, 88, 91, 92, 92, 93, 94, 95, 102, 103, 103, 103, 103, 103, 103, 103, 103, 105, 105, 107, 108, 125, 126, 128],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 65, 66, 67, 67, 72, 72, 72, 72, 72, 88, 88, 88, 88, 91, 92, 92, 93, 94, 95, 102, 103, 103, 103, 103, 103, 103, 103, 103, 105, 105, 107, 108, 125, 126, 128],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 65, 66, 67, 67, 72, 72, 72, 72, 72, 88, 88, 88, 88, 91, 92, 92, 93, 94, 95, 102, 103, 103, 103, 103, 103, 103, 103, 103, 105, 105, 107, 108, 125, 126, 128],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 65, 66, 67, 67, 72, 72, 72, 72, 72, 88, 88, 88, 88, 91, 92, 92, 93, 94, 95, 102, 103, 103, 103, 103, 103, 103, 103, 103, 105, 105, 107, 108, 125, 126, 128],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 65, 66, 67, 67, 72, 72, 72, 72, 72, 88, 88, 88, 88, 91, 92, 92, 93, 94, 95, 102, 103, 103, 103, 103, 103, 103, 103, 103, 105, 105, 107, 108, 125, 126, 128],
 [0, 1, 2, 4, 11, 11, 11, 12, 13, 14, 15, 16, 17, 19, 20, 20, 20, 22, 23, 24, 25, 25, 27, 28, 28, 29, 31, 32, 33, 37, 37, 46, 47, 64, 65, 66, 67, 68, 73, 73, 73, 73, 73, 89, 89, 89, 89, 92, 93, 93, 94, 95, 96, 103, 104, 104, 104, 104, 104, 104, 104, 104, 106, 106, 108, 109, 126, 127, 129],
]
                    # fmt: on
                ),
            )
        )

    def check_alignment_pfam9(self, alignment):
        """Check the alignment obtained by parsing Pfam record SH3_11."""
        self.assertEqual(alignment.annotations["identifier"], "SH3_11")
        self.assertEqual(alignment.annotations["accession"], "PF18103.3")
        self.assertEqual(
            alignment.annotations["definition"],
            "Retroviral integrase C-terminal SH3 domain",
        )
        self.assertEqual(
            alignment.annotations["author"], ["El-Gebali S;0000-0003-1378-5495"]
        )
        self.assertEqual(alignment.annotations["source of seed"], "ECOD:EUF00899")
        self.assertEqual(alignment.annotations["gathering method"], "25.00 25.00;")
        self.assertEqual(alignment.annotations["trusted cutoff"], "25.20 78.10;")
        self.assertEqual(alignment.annotations["noise cutoff"], "24.90 24.20;")
        self.assertEqual(
            alignment.annotations["build method"], "hmmbuild HMM.ann SEED.ann"
        )
        self.assertEqual(
            alignment.annotations["search method"],
            "hmmsearch -Z 57096847 -E 1000 --cpu 4 HMM pfamseq",
        )
        self.assertEqual(alignment.annotations["type"], "Domain")
        self.assertEqual(
            alignment.annotations["wikipedia"], ["Integrase", "SH3_domain"]
        )
        self.assertEqual(alignment.annotations["clan"], "CL0010")
        self.assertEqual(len(alignment.annotations["references"]), 1)
        self.assertEqual(alignment.annotations["references"][0]["number"], 1)
        self.assertEqual(alignment.annotations["references"][0]["medline"], "20118915")
        self.assertEqual(
            alignment.annotations["references"][0]["title"],
            "Retroviral intasome assembly and inhibition of DNA strand transfer.",
        )
        self.assertEqual(
            alignment.annotations["references"][0]["author"],
            "Hare S, Gupta SS, Valkov E, Engelman A, Cherepanov P;",
        )
        self.assertEqual(
            alignment.annotations["references"][0]["location"],
            "Nature. 2010;464:232-236.",
        )
        self.assertEqual(len(alignment.annotations["database references"]), 2)
        self.assertEqual(
            alignment.annotations["database references"][0]["reference"],
            "INTERPRO; IPR040903;",
        )
        self.assertEqual(
            alignment.annotations["database references"][1]["reference"],
            "SO; 0000417; polypeptide_domain;",
        )
        self.assertEqual(
            alignment.annotations["comment"],
            "This is the carboxy-terminal domain (CTD) found in retroviral integrase, an essential retroviral enzyme that binds both termini of linear viral DNA and inserts them into a host cell chromosome. The CTD adopts an SH3-like fold. Each CTD makes contact with the phosphodiester backbone of both viral DNA molecules, essentially crosslinking the structure [1].",
        )
        self.assertEqual(len(alignment.sequences), 1)
        self.assertEqual(alignment.sequences[0].id, "POL_SFVCP/1064-1126")
        self.assertEqual(alignment.sequences[0].annotations["accession"], "Q87040.1")
        self.assertEqual(
            alignment.sequences[0].seq,
            "RSWSPVVGQLVQERVARPASLRPRWHKPSTVLEVLNPRTVVILDHLGNNRTVSIDNLKPTSHQ",
        )
        self.assertEqual(
            alignment.column_annotations["consensus sequence"],
            "RSWSPVVGQLVQERVARPASLRPRWHKPSTVLEVLNPRTVVILDHLGNNRTVSIDNLKPTSHQ",
        )
        self.assertTrue(
            numpy.array_equal(alignment.coordinates, numpy.array([[0, 63]]))
        )

    def check_alignment_rfam1(self, alignment):
        """Check the alignment obtained by parsing Rfam record BTnc005."""
        self.assertEqual(alignment.annotations["accession"], "RF04178")
        self.assertEqual(alignment.annotations["identifier"], "BTnc005")
        self.assertEqual(
            alignment.annotations["definition"], "Bacteroides sRNA BTnc005"
        )
        self.assertEqual(
            alignment.annotations["author"],
            [
                "Prezza, G",
                "Ryan, D",
                "Mädler, G",
                "Barquist, L; 0000-0003-4732-2667",
                "Westermann, A",
            ],
        )
        self.assertEqual(
            alignment.annotations["source of seed"], "Published; PMID:32678091;"
        )
        self.assertEqual(
            alignment.annotations["source of structure"], "Published; PMID:32678091;"
        )
        self.assertEqual(alignment.annotations["gathering method"], "174.80")
        self.assertEqual(alignment.annotations["trusted cutoff"], "179.30")
        self.assertEqual(alignment.annotations["noise cutoff"], "174.30")
        self.assertEqual(alignment.annotations["type"], "Gene; sRNA;")
        self.assertEqual(alignment.annotations["build method"], "cmbuild -F CM SEED")
        self.assertEqual(
            alignment.annotations["calibration method"], "cmcalibrate --mpi CM"
        )
        self.assertEqual(
            alignment.annotations["search method"],
            "cmsearch --cpu 4 --verbose --nohmmonly -T 30.00 -Z 742849.287494 CM SEQDB",
        )
        self.assertEqual(len(alignment.annotations["database references"]), 1)
        self.assertEqual(
            alignment.annotations["database references"][0],
            {"reference": "SO; 0000655; ncRNA;"},
        )
        self.assertEqual(len(alignment.annotations["references"]), 1)
        self.assertEqual(alignment.annotations["references"][0]["number"], 1)
        self.assertEqual(alignment.annotations["references"][0]["medline"], "32678091")
        self.assertEqual(
            alignment.annotations["references"][0]["title"],
            "A high-resolution transcriptome map identifies small RNA regulation of metabolism in the gut microbe Bacteroides thetaiotaomicron.",
        )
        self.assertEqual(
            alignment.annotations["references"][0]["author"],
            "Ryan D, Jenniches L, Reichardt S, Barquist L, Westermann AJ",
        )
        self.assertEqual(
            alignment.annotations["references"][0]["location"],
            "Nat Commun. 2020;11:3557.",
        )
        self.assertEqual(
            alignment.annotations["comment"],
            "An uncharacterized small RNA discovered in Bacteroides thetaiotaomicron [1]",
        )
        self.assertEqual(
            alignment.annotations["wikipedia"], ["Bacteroides_thetaiotaomicron_sRNA"]
        )
        self.assertEqual(len(alignment.sequences), 3)
        self.assertEqual(alignment.sequences[0].id, "AE015928.1/72774-72978")
        self.assertEqual(alignment.sequences[1].id, "CP000139.1/2819055-2819247")
        self.assertEqual(alignment.sequences[2].id, "FP929033.1/4930704-4930908")
        self.assertEqual(
            alignment.sequences[0].seq,
            "GUAAGUAAAAGUGUAACAGGAAGAAAGUUGCAGCAUAUAUGCGGUGAAUUAUGCGGUGUCAUAGGAAUUGAGGAUUUAUGUAAGAUGCUGAUAAUGAGUAAGGAACCUUAAAGUUAAUCGUUCCCUGUCUCUCCGCAGAACCUACUGGACAAAACAGGACAGUAAGUGGACAAAAACCUACAAAUCAGCGAUUUGUAGGUUUUUU",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "AAAAGUAAGAGUGUAACAGGAAGAAAGUUGCAGCAUAUACGCGGUGAAUUAUUCGGUGUCAUAGGAGUAGAGUCUUUUGGUAAGAUGCUGAUAAUGAGUAGGGGAGAUGAAAGUUAAUCGUUCCCUGUCUCUCCGCUGGAAAGAAUUGCAAAACAAAGAAAAUCCCUGUAAAUUAAUACUUUACGGGGAUUUU",
        )
        self.assertEqual(
            alignment.sequences[2].seq,
            "GUAAGUAAAAGUGUAACAGGAAGAAAGUUGCAGCAUAUAUGCGGUGAAUUAUGCGGUGUCAUAGGAAUUGAGGAUUUAUGUAAGAUGCUGAUAAUGAGUAAGGAACCUUAAAGUUAAUCGUUCCCUGUCUCUCCGCUGAACUAUCCGGACAAAACCGGGCAAUGAACAGUCAAAUCCCACAAAUUCAAUGAUUUGUGGGACUUUU",
        )
        self.assertEqual(
            alignment[0],
            "GUAAGUAAAAGUGUAACAGGAAGAAAGUUGCAGCAUAUAUGCGGUGAAUUAUGCGGUGUCAUAGGAAUUGAGGAUUUAUGUAAGAUGCUGAUAAUGAGUAAGGAACCUUAAAGUUAAUCGUUCCCUGUCUCUCCGCAGAACCUACUGGACAAAACAGGACAGUAAGUGGACAAAAACCUACAAAUCAGC-GAUUUGUAGGUUUUUU",
        )
        self.assertEqual(
            alignment[1],
            "AAAAGUAAGAGUGUAACAGGAAGAAAGUUGCAGCAUAUACGCGGUGAAUUAUUCGGUGUCAUAGGAGUAGAGUCUUUUGGUAAGAUGCUGAUAAUGAGUAGGGGAGAUGAAAGUUAAUCGUUCCCUGUCUCUCCGCUGG---------AAAGAAUUGCAAAACAA--AGA-AAAUCCCUGUAAAUUAAU-ACUUUACGGGGAUUUU",
        )
        self.assertEqual(
            alignment[2],
            "GUAAGUAAAAGUGUAACAGGAAGAAAGUUGCAGCAUAUAUGCGGUGAAUUAUGCGGUGUCAUAGGAAUUGAGGAUUUAUGUAAGAUGCUGAUAAUGAGUAAGGAACCUUAAAGUUAAUCGUUCCCUGUCUCUCCGCUGAACUAUCCGGACAAAACCGGGCAAUGAACAGUCAAA-UCCCACAAAUUCAAUGAUUUGUGGGACUUUU",
        )
        self.assertEqual(
            alignment.column_annotations["consensus secondary structure"],
            ":::::::::::<<<<<<_________>>>>>>,,,,,,,,((((,,,<<<<<-<<<<<<<----<<<_______>>>------>>>>>>>>>>>><<<<<-<<<<_______________>>>>->>->>>,))))---------------------------------------<<<<<<<<<<<____>>>>>>>>>>>:::::",
        )
        self.assertEqual(
            alignment.column_annotations["reference coordinate annotation"],
            "guAAGUAAaAGuGuaaCAGGAAGAAAGuugCaGCAUAUAuGCGGUGAauuaugCgGuguCAUAGgaaUuGAGgauuuauGUAAGaugCuGauaauGaGuaaGGaaccUuAAAGUUAAUCGuuCCCugUCuCUCCGCuGaACuaaCuGGAcAaAAcuGgacAauaAauaGaCAAAacCCcgcaaaucaau.gauuugcgGGguUUUU",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [
                        [0, 139, 148, 165, 167, 170, 171, 174, 175, 189, 189, 205],
                        [0, 139, 139, 156, 156, 159, 159, 162, 163, 177, 177, 193],
                        [0, 139, 148, 165, 167, 170, 171, 174, 174, 188, 189, 205],
                    ]
                ),
            )
        )

    def check_alignment_rfam2(self, alignment):
        """Check the alignment obtained by parsing Rfam record SraC_RyeA."""
        self.assertEqual(alignment.annotations["accession"], "RF00101")
        self.assertEqual(alignment.annotations["identifier"], "SraC_RyeA")
        self.assertEqual(alignment.annotations["definition"], "SraC/RyeA RNA")
        self.assertEqual(
            alignment.annotations["author"], ["Bateman A; 0000-0002-6982-4660"]
        )
        self.assertEqual(alignment.annotations["source of seed"], "Bateman A")
        self.assertEqual(
            alignment.annotations["source of structure"], "Predicted; PFOLD"
        )
        self.assertEqual(alignment.annotations["gathering method"], "37.00")
        self.assertEqual(alignment.annotations["trusted cutoff"], "37.20")
        self.assertEqual(alignment.annotations["noise cutoff"], "36.80")
        self.assertEqual(alignment.annotations["type"], "Gene; sRNA;")
        self.assertEqual(alignment.annotations["build method"], "cmbuild -F CM SEED")
        self.assertEqual(
            alignment.annotations["calibration method"], "cmcalibrate --mpi CM"
        )
        self.assertEqual(
            alignment.annotations["search method"],
            "cmsearch --cpu 4 --verbose --nohmmonly -E 1000 -Z 549862.597050 CM SEQDB",
        )
        self.assertEqual(alignment.annotations["clan"], "CL00105")
        self.assertEqual(len(alignment.annotations["database references"]), 1)
        self.assertEqual(
            alignment.annotations["database references"][0],
            {"reference": "SO; 0000655; ncRNA;"},
        )
        self.assertEqual(len(alignment.annotations["references"]), 2)
        self.assertEqual(alignment.annotations["references"][0]["number"], 1)
        self.assertEqual(alignment.annotations["references"][0]["medline"], "11448770")
        self.assertEqual(
            alignment.annotations["references"][0]["title"],
            "Novel small RNA-encoding genes in the intergenic regions of Escherichia coli.",
        )
        self.assertEqual(
            alignment.annotations["references"][0]["author"],
            "Argaman L, Hershberg R, Vogel J, Bejerano G, Wagner EG, Margalit H, Altuvia S",
        )
        self.assertEqual(
            alignment.annotations["references"][0]["location"],
            "Curr Biol 2001;11:941-950.",
        )
        self.assertEqual(alignment.annotations["references"][1]["number"], 2)
        self.assertEqual(alignment.annotations["references"][1]["medline"], "11445539")
        self.assertEqual(
            alignment.annotations["references"][1]["title"],
            "Identification of novel small RNAs using comparative genomics and microarrays.",
        )
        self.assertEqual(
            alignment.annotations["references"][1]["author"],
            "Wassarman KM, Repoila F, Rosenow C, Storz G, Gottesman S",
        )
        self.assertEqual(
            alignment.annotations["references"][1]["location"],
            "Genes Dev 2001;15:1637-1651.",
        )
        self.assertEqual(
            alignment.annotations["comment"],
            "This RNA was discovered in E. coli during a large scale screens [1-2]. The function of this RNA is unknown. This RNA overlaps RFAM:RF00111 on the opposite strand suggesting that the two may act in a concerted manner.",
        )
        self.assertEqual(alignment.annotations["wikipedia"], ["SraC/RyeA_RNA"])
        self.assertEqual(len(alignment.sequences), 13)
        self.assertEqual(alignment.sequences[0].id, "AL627272.1/127686-127830")
        self.assertEqual(alignment.sequences[1].id, "CP000653.1/2601613-2601758")
        self.assertEqual(alignment.sequences[2].id, "AE017042.1/1756200-1756347")
        self.assertEqual(alignment.sequences[3].id, "CP000034.1/1046100-1046244")
        self.assertEqual(alignment.sequences[4].id, "CP000647.1/2580976-2581120")
        self.assertEqual(alignment.sequences[5].id, "AM286415.1/1991675-1991530")
        self.assertEqual(alignment.sequences[6].id, "CU928145.2/2074283-2074427")
        self.assertEqual(alignment.sequences[7].id, "CP000970.1/1336993-1336849")
        self.assertEqual(alignment.sequences[8].id, "AM933172.1/1226335-1226191")
        self.assertEqual(alignment.sequences[9].id, "AALD02000029.1/37435-37580")
        self.assertEqual(alignment.sequences[10].id, "AALC02000009.1/70496-70641")
        self.assertEqual(alignment.sequences[11].id, "AALF02000003.1/121616-121765")
        self.assertEqual(alignment.sequences[12].id, "AALE02000013.1/38-183")
        self.assertEqual(
            alignment.sequences[0].seq,
            "AAUUAAAAAAAGACCGAAUACGAUUCCUGUAUUCGGUCCAGGGAAAUGGCUCUUGGGAGAGAGCCGUGCGCUAAAAGUUGGCAUUAAUGCAGGCUAAAUCGCCUUGCUCUUUAAGAAUAGAUGACGACGCCAGGUUUUCCAGUUU",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "AGAUAAAAAGAGACCGAAUACGAUUCCUGUAUUCGGUCCAGGGAAAUGGCUCUUGGGAUAGAGCCGUGCGCUAAAAGUUGGCAUUAAUGCAGGCUAACUACGCCUGACACUCUAAGAAUAGAUGACGACGCCAGGUUUUCCAGUCC",
        )
        self.assertEqual(
            alignment.sequences[2].seq,
            "AAUUAAAAAAAGACCGAAUACGAUUCCUGAUAUUCGGUCUAGGGAAAUGGCUCUUGGGAGAGAGCCGUGCGCUAAAAGUUGGUAUUAAUGUAGGCUUAUUCAGCCGCACUUCUUAAGCGUAGCCGAGUACCGACAUUUCGCCAACCUU",
        )
        self.assertEqual(
            alignment.sequences[3].seq,
            "AAUGAAAAAAAGACCGAAUACGAUUCCUGUAUUCGGUCCAGGGAAAUGGCUCUUGGGAGAGAGCCGUGCGCUAAAAGUUGGCAUUAAUGCAGGCUUAGUUGCCUUGCCCUUUAAGAAUAGAUGACGACGCCAGGUUUUCCAGUUU",
        )
        self.assertEqual(
            alignment.sequences[4].seq,
            "AGAUAAAAAGAGACCGAAUACGAUUCCUGUAUUCGGUCCAGGGAAAUGGCUCUUGGGAGAGAGCCGUGCGCUAAAAGUUGGCAUUAAUGCAGGCUAAGUCGCCUUGCACUAUAAGAAUAGUUUAACGCGUCAGCUUUUCCAGUCC",
        )
        self.assertEqual(
            alignment.sequences[5].seq,
            "AAUUAAAAAGAGACCGAAUACGAUUCCUAUAUUCGGUCUAGGGAAAUGGCUCUUGGGAGAGAGCCGUGCGCUAAAAGUUGGCAUUAACGUAGGCUUGUUCAGCCAUACUCUUUAAGAGUAGUCGAGGUCAUGUGUUUCGCCAACUU",
        )
        self.assertEqual(
            alignment.sequences[6].seq,
            "AGAUAAAAAGAGACCGAACACGAUUCCUGUAUUCGGUCCAGGGAAAUGGCUCUUGGGAGAGAGCCGUGCGCUAAAAGUUGGCAUUAAUGCAGGCUUAGUUGCCUUGCCCUUUAAGAAUAGAUGACGACGCCAGGUUUUCCAGUUU",
        )
        self.assertEqual(
            alignment.sequences[7].seq,
            "AGAUAAAAAGAGACCGAACACGAUUCCUGUAUUCGGUCCAGGGAAAUGGCUCUUGGGAGAGAGCCGUGCGCUAAAAGUUGGCAUUAAUGCAGGCUUAGUUGCCUUGCUCUUUAAGAAUAGAUGACGACGCCAGAUUUUCCAGUUU",
        )
        self.assertEqual(
            alignment.sequences[8].seq,
            "AAAUAAAAAGAGACCGAAUACGAUUCCUGUAUUCGGUCCAGGGAAAUGGCUCUUGGGAGAGAGCCGUGCGCUAAAAGUUGGCAUUAAUGCAGGCUAAAUCGCCUUGCCCUUUAAGAAUAGAUGACGACGCCAGGUUUUCCAGUUU",
        )
        self.assertEqual(
            alignment.sequences[9].seq,
            "ACUUAAAAAGAGACCGAAUACGAUUCCUAUAUUCGGUCUAGGGAAAUGGCUCUUGGGAGAGAGCCGUGCGCUAAAAGUUGGCAUUAACGUAGGCUUGUUCAGCCGUACUACUUAAGCGUAGUCGAGUACAUGUGUUUCGCCAACUU",
        )
        self.assertEqual(
            alignment.sequences[10].seq,
            "AGAUAAAAAAAGACCGAAUACGAUUCCUAUAUUCGGUCUAGGGAAAUGGCUCUUGGGAGAGAGCCGUGCGCUAAAAGUUGGCAUUAACGUAGGCUUAUUCAGCCGUACUCCUUAAGCGUAGUCGAGUACAUGUAUUUAGCCAACUU",
        )
        self.assertEqual(
            alignment.sequences[11].seq,
            "AAUUAAAAAAAGACCGAAUACGAUUCCUAUAUUCGGUCUAGGGAAAUGGCUCUUGGGACAGAGCCGUGCGCUAAAAGUUGGCAUUAAUUAACGUAGGCUUAUUCAGCCGUACUCCUUAAGCGUAGUCGAGUACAUGUGUUUAGCCAACUU",
        )
        self.assertEqual(
            alignment.sequences[12].seq,
            "AGUUAAAAAAAGACCGAAUACGAUUCCUAUAUUCGGUCUAGGGAAAGGGCUCUUGGGAGAGAGCCGUGCGCUAAAAGUUGGCAUUAACGUAGGCUUAUUCAGCCGUACUCCUUAAGCGUAGUCGAGUACAUGUGUUUCGCCAACUU",
        )
        self.assertEqual(
            alignment[0],
            "AAUUAAAAAAAGACCGAAUACGAUUCCUG-UAUUCGGUCCAGGGAAAUGGCUCUUGGGAGAGAGCCGUGCGCUAAAAGUUGGCAUUAAU----GCAGGCUAAAU-C-GCCUUGCUCUUUAAGAAUAGAUGACGA-CGCCAGGUUUUCCAGUUU",
        )
        self.assertEqual(
            alignment[1],
            "AGAUAAAAAGAGACCGAAUACGAUUCCUG-UAUUCGGUCCAGGGAAAUGGCUCUUGGGAUAGAGCCGUGCGCUAAAAGUUGGCAUUAAU----GCAGGCUAACUAC-GCCUGACACUCUAAGAAUAGAUGACGA-CGCCAGGUUUUCCAGUCC",
        )
        self.assertEqual(
            alignment[2],
            "AAUUAAAAAAAGACCGAAUACGAUUCCUGAUAUUCGGUCUAGGGAAAUGGCUCUUGGGAGAGAGCCGUGCGCUAAAAGUUGGUAUUAAU----GUAGGCUUAUU-CAGCCGCACUUCUUAAGCGUAGCCGAGUACCGACAUUUCGCCAACCUU",
        )
        self.assertEqual(
            alignment[3],
            "AAUGAAAAAAAGACCGAAUACGAUUCCUG-UAUUCGGUCCAGGGAAAUGGCUCUUGGGAGAGAGCCGUGCGCUAAAAGUUGGCAUUAAU----GCAGGCUUAGU-U-GCCUUGCCCUUUAAGAAUAGAUGACGA-CGCCAGGUUUUCCAGUUU",
        )
        self.assertEqual(
            alignment[4],
            "AGAUAAAAAGAGACCGAAUACGAUUCCUG-UAUUCGGUCCAGGGAAAUGGCUCUUGGGAGAGAGCCGUGCGCUAAAAGUUGGCAUUAAU----GCAGGCUAAGU-C-GCCUUGCACUAUAAGAAUAGUUUAACG-CGUCAGCUUUUCCAGUCC",
        )
        self.assertEqual(
            alignment[5],
            "AAUUAAAAAGAGACCGAAUACGAUUCCUA-UAUUCGGUCUAGGGAAAUGGCUCUUGGGAGAGAGCCGUGCGCUAAAAGUUGGCAUUAAC----GUAGGCUUGUU-CAGCCAUACUCUUUAAGAGUAGUCGAGGU-CAUGUGUUUCGCCAACUU",
        )
        self.assertEqual(
            alignment[6],
            "AGAUAAAAAGAGACCGAACACGAUUCCUG-UAUUCGGUCCAGGGAAAUGGCUCUUGGGAGAGAGCCGUGCGCUAAAAGUUGGCAUUAAU----GCAGGCUUAGU-U-GCCUUGCCCUUUAAGAAUAGAUGACGA-CGCCAGGUUUUCCAGUUU",
        )
        self.assertEqual(
            alignment[7],
            "AGAUAAAAAGAGACCGAACACGAUUCCUG-UAUUCGGUCCAGGGAAAUGGCUCUUGGGAGAGAGCCGUGCGCUAAAAGUUGGCAUUAAU----GCAGGCUUAGU-U-GCCUUGCUCUUUAAGAAUAGAUGACGA-CGCCAGAUUUUCCAGUUU",
        )
        self.assertEqual(
            alignment[8],
            "AAAUAAAAAGAGACCGAAUACGAUUCCUG-UAUUCGGUCCAGGGAAAUGGCUCUUGGGAGAGAGCCGUGCGCUAAAAGUUGGCAUUAAU----GCAGGCUAAAU-C-GCCUUGCCCUUUAAGAAUAGAUGACGA-CGCCAGGUUUUCCAGUUU",
        )
        self.assertEqual(
            alignment[9],
            "ACUUAAAAAGAGACCGAAUACGAUUCCUA-UAUUCGGUCUAGGGAAAUGGCUCUUGGGAGAGAGCCGUGCGCUAAAAGUUGGCAUUAAC----GUAGGCUUGUU-CAGCCGUACUACUUAAGCGUAGUCGAGUA-CAUGUGUUUCGCCAACUU",
        )
        self.assertEqual(
            alignment[10],
            "AGAUAAAAAAAGACCGAAUACGAUUCCUA-UAUUCGGUCUAGGGAAAUGGCUCUUGGGAGAGAGCCGUGCGCUAAAAGUUGGCAUUAAC----GUAGGCUUAUU-CAGCCGUACUCCUUAAGCGUAGUCGAGUA-CAUGUAUUUAGCCAACUU",
        )
        self.assertEqual(
            alignment[11],
            "AAUUAAAAAAAGACCGAAUACGAUUCCUA-UAUUCGGUCUAGGGAAAUGGCUCUUGGGACAGAGCCGUGCGCUAAAAGUUGGCAUUAAUUAACGUAGGCUUAUU-CAGCCGUACUCCUUAAGCGUAGUCGAGUA-CAUGUGUUUAGCCAACUU",
        )
        self.assertEqual(
            alignment[12],
            "AGUUAAAAAAAGACCGAAUACGAUUCCUA-UAUUCGGUCUAGGGAAAGGGCUCUUGGGAGAGAGCCGUGCGCUAAAAGUUGGCAUUAAC----GUAGGCUUAUU-CAGCCGUACUCCUUAAGCGUAGUCGAGUA-CAUGUGUUUCGCCAACUU",
        )
        self.assertEqual(
            alignment.column_annotations["consensus secondary structure"],
            ":::::::::::<<<<<<<<<<_______>.>>>>>>>>>,,,,,,,<<<<<<<<______>>>>>>>>,,<<_________>>,,,,,,....,,,<<<_____._.>>>::::::::::::::::::::::::.::::::::::::::::::",
        )
        self.assertEqual(
            alignment.column_annotations["reference coordinate annotation"],
            "AauUAAAAAaAGaCCGaauacGAUUCCUg.uauuCGGuCuAGGGAAauGGCuCuUGGGAgaGaGCCguGCGCUAAAAGUUGGCAUUAAu....GuAGGCUuAuU.c.GCCuuaCucuUUAAGaaUAGuuGAguA.CgucaguUUuuCcAauUU",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [
                        [0, 29, 29, 88, 88, 99, 99, 100, 100, 127, 127, 145],
                        [0, 29, 29, 88, 88, 99, 100, 101, 101, 128, 128, 146],
                        [0, 29, 30, 89, 89, 100, 100, 101, 102, 129, 130, 148],
                        [0, 29, 29, 88, 88, 99, 99, 100, 100, 127, 127, 145],
                        [0, 29, 29, 88, 88, 99, 99, 100, 100, 127, 127, 145],
                        [0, 29, 29, 88, 88, 99, 99, 100, 101, 128, 128, 146],
                        [0, 29, 29, 88, 88, 99, 99, 100, 100, 127, 127, 145],
                        [0, 29, 29, 88, 88, 99, 99, 100, 100, 127, 127, 145],
                        [0, 29, 29, 88, 88, 99, 99, 100, 100, 127, 127, 145],
                        [0, 29, 29, 88, 88, 99, 99, 100, 101, 128, 128, 146],
                        [0, 29, 29, 88, 88, 99, 99, 100, 101, 128, 128, 146],
                        [0, 29, 29, 88, 92, 103, 103, 104, 105, 132, 132, 150],
                        [0, 29, 29, 88, 88, 99, 99, 100, 101, 128, 128, 146],
                    ]
                ),
            )
        )

    def check_alignment_rfam3(self, alignment):
        """Check the alignment obtained by parsing Rfam record McaS."""
        self.assertEqual(alignment.annotations["accession"], "RF00115")
        self.assertEqual(alignment.annotations["identifier"], "McaS")
        self.assertEqual(alignment.annotations["previous identifier"], "IS061;")
        self.assertEqual(alignment.annotations["definition"], "McaS/IsrA RNA")
        self.assertEqual(
            alignment.annotations["author"], ["Argasinska J; 0000-0003-2678-2824"]
        )
        self.assertEqual(alignment.annotations["source of seed"], "Argasinska J")
        self.assertEqual(
            alignment.annotations["source of structure"], "Predicted; 22289118"
        )
        self.assertEqual(alignment.annotations["gathering method"], "42.00")
        self.assertEqual(alignment.annotations["trusted cutoff"], "42.10")
        self.assertEqual(alignment.annotations["noise cutoff"], "38.70")
        self.assertEqual(alignment.annotations["type"], "Gene; sRNA;")
        self.assertEqual(alignment.annotations["build method"], "cmbuild -n -F CM SEED")
        self.assertEqual(
            alignment.annotations["calibration method"], "cmcalibrate --mpi CM"
        )
        self.assertEqual(
            alignment.annotations["search method"],
            "cmsearch --cpu 4 --verbose --nohmmonly -T 30.00 -Z 604040.189692 --mxsize 128 CM SEQDB",
        )
        self.assertEqual(alignment.annotations["clan"], "CL00106")
        self.assertEqual(len(alignment.annotations["database references"]), 3)
        self.assertEqual(
            alignment.annotations["database references"][0],
            {"reference": "SO; 0001263; ncRNA_gene;"},
        )
        self.assertEqual(
            alignment.annotations["database references"][1],
            {"reference": "GO; 0005515; protein binding;"},
        )
        self.assertEqual(
            alignment.annotations["database references"][2],
            {"reference": "GO; 0006417; regulation of translation;"},
        )
        self.assertEqual(len(alignment.annotations["references"]), 4)
        self.assertEqual(alignment.annotations["references"][0]["number"], 1)
        self.assertEqual(alignment.annotations["references"][0]["medline"], "12069726")
        self.assertEqual(
            alignment.annotations["references"][0]["title"],
            "A bioinformatics based approach to discover small RNA genes in the Escherichia coli genome.",
        )
        self.assertEqual(
            alignment.annotations["references"][0]["author"],
            "Chen S, Lesnik EA, Hall TA, Sampath R, Griffey RH, Ecker DJ, Blyn LB",
        )
        self.assertEqual(
            alignment.annotations["references"][0]["location"],
            "Biosystems 2002;65:157-177.",
        )
        self.assertEqual(alignment.annotations["references"][1]["number"], 2)
        self.assertEqual(alignment.annotations["references"][1]["medline"], "23666921")
        self.assertEqual(
            alignment.annotations["references"][1]["title"],
            "Dual function of the McaS small RNA in controlling biofilm formation.",
        )
        self.assertEqual(
            alignment.annotations["references"][1]["author"],
            "Jorgensen MG, Thomason MK, Havelund J, Valentin-Hansen P, Storz G",
        )
        self.assertEqual(
            alignment.annotations["references"][1]["location"],
            "Genes Dev. 2013;27:1132-1145.",
        )
        self.assertEqual(alignment.annotations["references"][2]["number"], 3)
        self.assertEqual(alignment.annotations["references"][2]["medline"], "22289118")
        self.assertEqual(
            alignment.annotations["references"][2]["title"],
            "A small RNA that regulates motility and biofilm formation in response to changes  in nutrient availability in Escherichia coli.",
        )
        self.assertEqual(
            alignment.annotations["references"][2]["author"],
            "Thomason MK, Fontaine F, De Lay N, Storz G",
        )
        self.assertEqual(
            alignment.annotations["references"][2]["location"],
            "Mol Microbiol. 2012;84:17-35.",
        )
        self.assertEqual(alignment.annotations["references"][3]["number"], 4)
        self.assertEqual(alignment.annotations["references"][3]["medline"], "26609136")
        self.assertEqual(
            alignment.annotations["references"][3]["title"],
            "Ribonucleoprotein particles of bacterial small non-coding RNA IsrA (IS61 or McaS) and its interaction with RNA polymerase core may link transcription to mRNA fate.",
        )
        self.assertEqual(
            alignment.annotations["references"][3]["author"],
            "van Nues RW, Castro-Roa D, Yuzenkova Y, Zenkin N",
        )
        self.assertEqual(
            alignment.annotations["references"][3]["location"],
            "Nucleic Acids Res. 2016;44:2577-2592.",
        )
        self.assertEqual(
            alignment.annotations["comment"],
            "This family consists of several bacterial RNA genes which are found between the abgR and ydaL genes in Escherichia coli and Shigella flexneri.[1] It was discovered using a computational screen of the E. coli genome.[1] Subsequent characterisation of ISO61 region has revealed that the reverse strand is actually a CsrA binding ncRNA called McaS and that it has a role in biofilm formation control.[2] Furthermore, it has been shown that McaS(IsrA) exists as a ribonucleoprotein particles (sRNPs), which involve a defined set of proteins including Hfq, S1, CsrA, ProQ and PNPase.[4]",
        )
        self.assertEqual(alignment.annotations["wikipedia"], ["IS061_RNA"])
        self.assertEqual(len(alignment.sequences), 4)
        self.assertEqual(alignment.sequences[0].id, "CP000036.1/1703842-1703937")
        self.assertEqual(alignment.sequences[1].id, "U00096.3/1405751-1405656")
        self.assertEqual(alignment.sequences[2].id, "CP000034.1/1309299-1309204")
        self.assertEqual(alignment.sequences[3].id, "CP011132.1/1732716-1732810")
        self.assertEqual(
            alignment.sequences[0].seq,
            "ACCGGCGCAGAGGAGACAAUGCCGGACUUAAGACGCGGAUGCACUGCUGUGUGUACUGUAGAGUCUGGCGGAUGUCGACAGACUCUAUUUUUUUAU",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "ACCGGCGCAGAGGAGACAAUGCCGGAUUUAAGACGCGGAUGCACUGCUGUGUGUACUGUAGAGUCUGGCGGAUGUCGACAGACUCUAUUUUUUUAU",
        )
        self.assertEqual(
            alignment.sequences[2].seq,
            "ACCGGUCACCAGGACCCCAGGCCGGAUUUAAGACGAGGAUGCACUGCUGUGUGUACUGUAGAGUCUGGCGGAUGUCGACAGGCUCUAUUUUUUUAU",
        )
        self.assertEqual(
            alignment.sequences[3].seq,
            "ACCCGCCACACGGAAUAAUAACGGGAACACAUGAAGGAUAAACUGCUGUGUGUACUGUAGAGUCUGGCGGAUGUCGACAGGCUCUAUUUUUUUAU",
        )
        self.assertEqual(
            alignment[0],
            "ACCGGCGCAGAGGAGACAAUGCCGGACUUAAGACGCGGAUGCACUGCUGUGUGUACUGUAGAGUCUGGCGGAUGUCGACAGACUCUAUUUUUUUAU",
        )
        self.assertEqual(
            alignment[1],
            "ACCGGCGCAGAGGAGACAAUGCCGGAUUUAAGACGCGGAUGCACUGCUGUGUGUACUGUAGAGUCUGGCGGAUGUCGACAGACUCUAUUUUUUUAU",
        )
        self.assertEqual(
            alignment[2],
            "ACCGGUCACCAGGACCCCAGGCCGGAUUUAAGACGAGGAUGCACUGCUGUGUGUACUGUAGAGUCUGGCGGAUGUCGACAGGCUCUAUUUUUUUAU",
        )
        self.assertEqual(
            alignment[3],
            "ACCCGCCACACGGAAUAAUAACGGGAACACAUG-AAGGAUAAACUGCUGUGUGUACUGUAGAGUCUGGCGGAUGUCGACAGGCUCUAUUUUUUUAU",
        )
        self.assertEqual(
            alignment.column_annotations["consensus secondary structure"],
            ":<<<<<<____________>>>>>>,,,,,,,<<<<<<________>>>>>>-----<<<<<<<<<<-<<_____>>->>>>>>>>>>::::::::",
        )
        self.assertEqual(
            alignment.column_annotations["reference coordinate annotation"],
            "ACCgGccaaaaGGAaacaaggCcGGAuuuaAgaCgcgGAUgcACUGCugcGuGUACUguaGaGuCuGGCGGAUGUCGACaGaCuCuauUUUUUUAU",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [[0, 33, 34, 96], [0, 33, 34, 96], [0, 33, 34, 96], [0, 33, 33, 95]]
                ),
            )
        )

    def check_alignment_rfam4(self, alignment):
        """Check the alignment obtained by parsing Rfam record IRES_KSHV."""
        self.assertEqual(alignment.annotations["accession"], "RF00511")
        self.assertEqual(alignment.annotations["identifier"], "IRES_KSHV")
        self.assertEqual(
            alignment.annotations["definition"],
            "Kaposi's sarcoma-associated herpesvirus internal ribosome entry site",
        )
        self.assertEqual(
            alignment.annotations["author"], ["Moxon SJ; 0000-0003-4644-1816"]
        )
        self.assertEqual(
            alignment.annotations["source of seed"], "Published; 11160685, INFERNAL"
        )
        self.assertEqual(
            alignment.annotations["source of structure"], "Published; PMID:11160685"
        )
        self.assertEqual(alignment.annotations["gathering method"], "100.00")
        self.assertEqual(alignment.annotations["trusted cutoff"], "317.10")
        self.assertEqual(alignment.annotations["noise cutoff"], "30.10")
        self.assertEqual(alignment.annotations["type"], "Cis-reg; IRES;")
        self.assertEqual(alignment.annotations["build method"], "cmbuild -F CM SEED")
        self.assertEqual(
            alignment.annotations["calibration method"], "cmcalibrate --mpi CM"
        )
        self.assertEqual(
            alignment.annotations["search method"],
            "cmsearch --cpu 4 --verbose --nohmmonly -E 1000 -Z 549862.597050 CM SEQDB",
        )
        self.assertEqual(len(alignment.annotations["database references"]), 2)
        self.assertEqual(
            alignment.annotations["database references"][0],
            {"reference": "SO; 0000243; internal_ribosome_entry_site;"},
        )
        self.assertEqual(
            alignment.annotations["database references"][1],
            {"reference": "GO; 0043022; ribosome binding;"},
        )
        self.assertEqual(alignment.annotations["references"][0]["number"], 1)
        self.assertEqual(alignment.annotations["references"][0]["medline"], "11160685")
        self.assertEqual(
            alignment.annotations["references"][0]["title"],
            "Kaposi's sarcoma-associated herpesvirus vCyclin open reading frame contains an internal ribosome entry site.",
        )
        self.assertEqual(
            alignment.annotations["references"][0]["author"], "Bieleski L, Talbot SJ"
        )
        self.assertEqual(
            alignment.annotations["references"][0]["location"],
            "J Virol 2001;75:1864-1869.",
        )
        self.assertEqual(alignment.annotations["references"][1]["number"], 2)
        self.assertEqual(alignment.annotations["references"][1]["medline"], "14993645")
        self.assertEqual(
            alignment.annotations["references"][1]["title"],
            "A polypyrimidine tract facilitates the expression of Kaposi's sarcoma-associated herpesvirus vFLIP through an internal ribosome entry site.",
        )
        self.assertEqual(
            alignment.annotations["references"][1]["author"],
            "Bieleski L, Hindley C, Talbot SJ",
        )
        self.assertEqual(
            alignment.annotations["references"][1]["location"],
            "J Gen Virol 2004;85:615-620.",
        )
        self.assertEqual(
            alignment.annotations["comment"],
            "This family represents the Kaposi's sarcoma-associated herpesvirus (KSHV) internal ribosome entry site (IRES) present in the vCyclin gene. The vCyclin and vFLIP coding sequences are present on a bicistronic transcript and it is thought the IRES may initiate translation of vFLIP from this bicistronic transcript [1,2].",
        )
        self.assertEqual(
            alignment.annotations["wikipedia"],
            [
                "Kaposi's_sarcoma-associated_herpesvirus_internal_ribosome_entry_site_(IRES)"
            ],
        )
        self.assertEqual(len(alignment.sequences), 5)
        self.assertEqual(alignment.sequences[0].id, "AF148805.2/123462-123215")
        self.assertEqual(alignment.sequences[1].id, "U40667.1/2005-2252")
        self.assertEqual(alignment.sequences[2].id, "U79416.1/354-601")
        self.assertEqual(alignment.sequences[3].id, "U93872.2/123729-123482")
        self.assertEqual(alignment.sequences[4].id, "U75698.1/123214-122967")
        self.assertEqual(
            alignment.sequences[0].seq,
            "UUGCUAUGCCGCGGCAGACUCCUUUUCCCGCCAAGAACUUAUAGACCAGGAGAAAGAACUCCUUGAGAAGUUGGCGUGGCGAACAGAGGCAGUCUUAGCGACGGACGUCACUUCCUUCUUGUUACUUAAAUUGCUGGGGGGCUCCCAACACCUGGACUUUUGGCACCACGAGGUCAACACCCUGAUUACAAAAGCCUUAGUUGACCCAAAGACUGGCUCAUUGCCCGCCUCUAUUAUCAGCGCUGCAG",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "UUGCUAUGCCGCGGCAGACUCCUUUUCCCGCCAAGAACUUAUAGACCAGGAGAAAGAACUCCUUGAGAAGUUGGCGUGGCGAACAGAGGCAGUCUUAGCGACGGACGUCACUUCCUUCUUGUUACUUAAAUUGGUGGGGGGCUCCCAACACCUGGACUUUUGGCACCACGAGGUCAACACCCUGAUUACAAAAGCCUUAGUUGACCCAAAGACUGGGUCAUUGCCCGCCUCUAUUAUCAGCGCUGCAG",
        )
        self.assertEqual(
            alignment.sequences[2].seq,
            "UUGCUAUGCCGCGGCAGACUCCUUUUCCCGCCAAGAACUUAUAGACCAGGAGAAAGAACUCCUUGAGAAGUUGGCGUGGCGAACAGAGGCAGUCUUAGCGACGGACGUAACUUCCUUCUUGUUACUUAAAUUGCUGGGGGGCUCCCAACACCUGGACUUUUGGCACCACGAGGUCGACACCCUGAUUACAAAAGCCUUAGUUGACCCAAAGACUGGCUCAUUGCCCGCCUCUAUUAUCAGCGCUGCAG",
        )
        self.assertEqual(
            alignment.sequences[3].seq,
            "UUGCUAUGCCGCGGCAGACUCCUUUUCCCGCCAAGAACUUAUAGACCAGGAGAAAGAACUCCUUGAGAAGUUGGCGUGGCGAACAGAGGCAGUCUUAGCGACGGACGUAACUUCCUUCUUGUUACUUAAAUUGCUGGGGGGCUCCCAACACCUGGACUUUUGGCACCACGAGGUCAACACCCUGAUUACAAAAGCCUUAGUUGACCCAAAGACUGGCUCAUUGCCCGCCUCUAUUAUCAGCGCUGCAG",
        )
        self.assertEqual(
            alignment.sequences[4].seq,
            "UUGCUAUGCCGCGGCAGACUCCUUUUCCCGCCAAGAACUUAUAGACCAGGAGAAAGAACUCCUUGAGAAGUUGGCGUGGCGAACAGAGGCAGUCUUAGCGACGGACGUCACUUCCUUCUUGUUACUUAAAUUGCUGGGGGGCUCCCAACACCUGGACUUUUGGCACCACGAGGUCAACACCCUGAUUACAAAAGCCUUAGUUGACCCAAAGACUGGCUCAUUGCCCGCCUCUAUUAUCAGCGCUGCAG",
        )
        self.assertEqual(
            alignment[0],
            "UUGCUAUGCCGCGGCAGACUCCUUUUCCCGCCAAGAACUUAUAGACCAGGAGAAAGAACUCCUUGAGAAGUUGGCGUGGCGAACAGAGGCAGUCUUAGCGACGGACGUCACUUCCUUCUUGUUACUUAAAUUGCUGGGGGGCUCCCAACACCUGGACUUUUGGCACCACGAGGUCAACACCCUGAUUACAAAAGCCUUAGUUGACCCAAAGACUGGCUCAUUGCCCGCCUCUAUUAUCAGCGCUGCAG",
        )
        self.assertEqual(
            alignment[1],
            "UUGCUAUGCCGCGGCAGACUCCUUUUCCCGCCAAGAACUUAUAGACCAGGAGAAAGAACUCCUUGAGAAGUUGGCGUGGCGAACAGAGGCAGUCUUAGCGACGGACGUCACUUCCUUCUUGUUACUUAAAUUGGUGGGGGGCUCCCAACACCUGGACUUUUGGCACCACGAGGUCAACACCCUGAUUACAAAAGCCUUAGUUGACCCAAAGACUGGGUCAUUGCCCGCCUCUAUUAUCAGCGCUGCAG",
        )
        self.assertEqual(
            alignment[2],
            "UUGCUAUGCCGCGGCAGACUCCUUUUCCCGCCAAGAACUUAUAGACCAGGAGAAAGAACUCCUUGAGAAGUUGGCGUGGCGAACAGAGGCAGUCUUAGCGACGGACGUAACUUCCUUCUUGUUACUUAAAUUGCUGGGGGGCUCCCAACACCUGGACUUUUGGCACCACGAGGUCGACACCCUGAUUACAAAAGCCUUAGUUGACCCAAAGACUGGCUCAUUGCCCGCCUCUAUUAUCAGCGCUGCAG",
        )
        self.assertEqual(
            alignment[3],
            "UUGCUAUGCCGCGGCAGACUCCUUUUCCCGCCAAGAACUUAUAGACCAGGAGAAAGAACUCCUUGAGAAGUUGGCGUGGCGAACAGAGGCAGUCUUAGCGACGGACGUAACUUCCUUCUUGUUACUUAAAUUGCUGGGGGGCUCCCAACACCUGGACUUUUGGCACCACGAGGUCAACACCCUGAUUACAAAAGCCUUAGUUGACCCAAAGACUGGCUCAUUGCCCGCCUCUAUUAUCAGCGCUGCAG",
        )
        self.assertEqual(
            alignment[4],
            "UUGCUAUGCCGCGGCAGACUCCUUUUCCCGCCAAGAACUUAUAGACCAGGAGAAAGAACUCCUUGAGAAGUUGGCGUGGCGAACAGAGGCAGUCUUAGCGACGGACGUCACUUCCUUCUUGUUACUUAAAUUGCUGGGGGGCUCCCAACACCUGGACUUUUGGCACCACGAGGUCAACACCCUGAUUACAAAAGCCUUAGUUGACCCAAAGACUGGCUCAUUGCCCGCCUCUAUUAUCAGCGCUGCAG",
        )
        self.assertEqual(
            alignment.column_annotations["consensus secondary structure"],
            ":::::::::<<<---<<<<<<<<<<---<<<<<--<<<<<-------<<<<<______>>>>>---->>>>>---->>>>>--->>>>>->>>>>-->>>,,<<<_______>>>-----------------((((((((((,,,<<<<----<<-<<<<<<------<<-<<<____>>>->>---->>>>>>>>--->>>>,,,,,,,,,,,<<<_____>>>,))))))-----)))):::::::",
        )
        self.assertEqual(
            alignment.column_annotations["reference coordinate annotation"],
            "UUGCUAUGCCGCGGCaGaCuCCucuUCCCGCCaAGaaCuuAUAGACCaGGaGAAAGAACuCCuUGAGaaGuuGGCGuGGCGAACagaGGCaGuCuUAGCGACGGaCGUaACUuCCUUCUUGUUACUUAAAUUGcuGgGgGGCUCCCaaCACCUGGACuuuuGGCACCACgAGGuCaACaCCCcGAUUACaaaaGCCUUAGuuGACCCAAAGACUGGcUCAUUgCCCGCCcCcAUUAUCagCGCUGCAG",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[0, 248], [0, 248], [0, 248], [0, 248], [0, 248]]),
            )
        )

    def check_alignment_rfam5(self, alignment):
        """Check the alignment obtained by parsing Rfam record BMV3_UPD-PK3."""
        self.assertEqual(alignment.annotations["accession"], "RF01113")
        self.assertEqual(alignment.annotations["identifier"], "BMV3_UPD-PK3")
        self.assertEqual(
            alignment.annotations["definition"],
            "Pseudoknot of upstream pseudoknot domain (UPD) of the 3'UTR",
        )
        self.assertEqual(
            alignment.annotations["author"], ["Wilkinson A; 0000-0001-7406-0151"]
        )
        self.assertEqual(alignment.annotations["source of seed"], "Pseudobase")
        self.assertEqual(alignment.annotations["source of structure"], "Pseudobase")
        self.assertEqual(alignment.annotations["gathering method"], "36.00")
        self.assertEqual(alignment.annotations["trusted cutoff"], "37.00")
        self.assertEqual(alignment.annotations["noise cutoff"], "33.30")
        self.assertEqual(alignment.annotations["type"], "Cis-reg;")
        self.assertEqual(alignment.annotations["build method"], "cmbuild -F CM SEED")
        self.assertEqual(
            alignment.annotations["calibration method"], "cmcalibrate --mpi CM"
        )
        self.assertEqual(
            alignment.annotations["search method"],
            "cmsearch --cpu 4 --verbose --nohmmonly -T 28.00 -Z 549862.597050 CM SEQDB",
        )
        self.assertEqual(len(alignment.annotations["database references"]), 4)
        self.assertEqual(
            alignment.annotations["database references"][0],
            {"reference": "PKBASE; PKB00156;"},
        )
        self.assertEqual(
            alignment.annotations["database references"][1],
            {"reference": "SO; 0005836; regulatory_region;"},
        )
        self.assertEqual(
            alignment.annotations["database references"][2],
            {"reference": "GO; 1904973; positive regulation of viral translation;"},
        )
        self.assertEqual(
            alignment.annotations["database references"][3],
            {"reference": "GO; 0046782; regulation of viral transcription;"},
        )
        self.assertEqual(len(alignment.annotations["references"]), 1)
        self.assertEqual(alignment.annotations["references"][0]["number"], 1)
        self.assertEqual(alignment.annotations["references"][0]["medline"], "7684465")
        self.assertEqual(
            alignment.annotations["references"][0]["title"],
            "Contributions of the brome mosaic virus RNA-3 3'-nontranslated region to replication and translation.",
        )
        self.assertEqual(
            alignment.annotations["references"][0]["author"],
            "Lahser FC, Marsh LE, Hall TC",
        )
        self.assertEqual(
            alignment.annotations["references"][0]["location"],
            "J Virol. 1993;67:3295-3303.",
        )
        self.assertEqual(alignment.annotations["wikipedia"], ["UPSK_RNA"])
        self.assertEqual(
            alignment.annotations["**"], "seedtax: Viruses; unclassified sequences"
        )
        self.assertEqual(len(alignment.sequences), 2)
        self.assertEqual(alignment.sequences[0].id, "X01678.1/2648-2670")
        self.assertEqual(alignment.sequences[1].id, "X58459.1/659-681")
        self.assertEqual(alignment.sequences[0].seq, "ACUUUGGCUAAGUUUAAAAGCUU")
        self.assertEqual(alignment.sequences[1].seq, "ACUUUGGCUAAGGUUAAAAGCUU")
        self.assertEqual(alignment[0], "ACUUUGGCUAAGUUUAAAAGCUU")
        self.assertEqual(alignment[1], "ACUUUGGCUAAGGUUAAAAGCUU")
        self.assertEqual(
            alignment.column_annotations["consensus secondary structure"],
            ":<<<_AAAA>>>::::::aaaa:",
        )
        self.assertEqual(
            alignment.column_annotations["reference coordinate annotation"],
            "ACUUUGGCUAAGuUUAAAAGCUU",
        )
        self.assertTrue(
            numpy.array_equal(alignment.coordinates, numpy.array([[0, 23], [0, 23]]))
        )

    def check_alignment_cath1(self, alignment):
        """Check the alignment obtained by parsing CATH record 3.30.160.60/FF/004774."""
        self.assertEqual(alignment.annotations["identifier"], "3.30.160.60/FF/004774")
        self.assertEqual(alignment.annotations["definition"], "Uncharacterized protein")
        self.assertEqual(alignment.annotations["accession"], "3.30.160.60/FF/004774")
        self.assertEqual(alignment.annotations["type"], "FunFam")
        self.assertEqual(len(alignment.annotations["database references"]), 2)
        self.assertEqual(
            alignment.annotations["database references"][0], {"reference": "CATH: v4.3"}
        )
        self.assertEqual(
            alignment.annotations["database references"][1],
            {"reference": "DOPS: 0.000"},
        )
        self.assertEqual(len(alignment.sequences), 1)
        self.assertEqual(alignment.sequences[0].annotations["accession"], "L7MZX4")
        self.assertEqual(
            alignment.sequences[0].annotations["organism"], "Anolis carolinensis"
        )
        self.assertEqual(alignment.sequences[0].description, "Uncharacterized protein")
        self.assertEqual(len(alignment.sequences[0].dbxrefs), 1)
        self.assertEqual(
            alignment.sequences[0].dbxrefs[0],
            "ORG; Eukaryota; Metazoa; Chordata; Craniata; Sarcopterygii; Lepidosauria; Squamata; Iguania; Dactyloidae; Anolis; Anolis carolinensis;",
        )
        self.assertEqual(alignment.sequences[0].id, "L7MZX4/382-398")
        self.assertEqual(alignment.sequences[0].seq, "GEKPYECLECGKRFTAR")
        self.assertEqual(alignment[0], "GEKPYECLECGKRFTAR")
        self.assertEqual(
            alignment.column_annotations["consensus score"], "00000000000000000"
        )
        self.assertEqual(
            alignment.column_annotations["consensus score 70"], "_________________"
        )
        self.assertEqual(
            alignment.column_annotations["consensus score 80"], "_________________"
        )
        self.assertEqual(
            alignment.column_annotations["consensus score 90"], "_________________"
        )
        self.assertTrue(
            numpy.array_equal(alignment.coordinates, numpy.array([[0, 17]]))
        )

    def check_alignment_cath2(self, alignment):
        """Check the alignment obtained by parsing CATH record 2.105.10.10/FF/000002."""
        self.assertEqual(alignment.annotations["identifier"], "2.105.10.10/FF/000002")
        self.assertEqual(alignment.annotations["definition"], "Adsorption protein P2")
        self.assertEqual(alignment.annotations["accession"], "2.105.10.10/FF/000002")
        self.assertEqual(alignment.annotations["type"], "FunFam")
        self.assertEqual(len(alignment.annotations["database references"]), 2)
        self.assertEqual(
            alignment.annotations["database references"][0], {"reference": "CATH: v4.3"}
        )
        self.assertEqual(
            alignment.annotations["database references"][1],
            {"reference": "DOPS: 0.000"},
        )
        self.assertEqual(alignment.sequences[0].annotations["accession"], "P27378")
        self.assertEqual(alignment.sequences[0].annotations["organism"], "")
        self.assertEqual(alignment.sequences[0].description, "Adsorption protein P2")
        self.assertEqual(len(alignment.sequences[0].dbxrefs), 2)
        self.assertEqual(alignment.sequences[0].dbxrefs[0], "ORG;")
        self.assertEqual(alignment.sequences[0].dbxrefs[1], "GO; GO:0019012;")
        self.assertEqual(len(alignment.sequences), 1)
        self.assertEqual(alignment.sequences[0].id, "P27378/2-64")
        self.assertEqual(
            alignment.sequences[0].seq,
            "ANFNVPKLGVFPVAAVFDIDNVPEDSSATGSRWLPSIYQGGNYWGGGPQALHAQVSNFDSSNR",
        )
        self.assertEqual(
            alignment[0],
            "ANFNVPKLGVFPVAAVFDIDNVPEDSSATGSRWLPSIYQGGNYWGGGPQALHAQVSNFDSSNR",
        )
        self.assertEqual(
            alignment.column_annotations["consensus score"],
            "000000000000000000000000000000000000000000000000000000000000000",
        )
        self.assertEqual(
            alignment.column_annotations["consensus score 70"],
            "_______________________________________________________________",
        )
        self.assertEqual(
            alignment.column_annotations["consensus score 80"],
            "_______________________________________________________________",
        )
        self.assertEqual(
            alignment.column_annotations["consensus score 90"],
            "_______________________________________________________________",
        )
        self.assertTrue(
            numpy.array_equal(alignment.coordinates, numpy.array([[0, 63]]))
        )

    def check_alignment_cath3(self, alignment):
        """Check the alignment obtained by parsing CATH record 1.10.275.10/FF/000026."""
        self.assertEqual(alignment.annotations["identifier"], "1.10.275.10/FF/000026")
        self.assertEqual(alignment.annotations["definition"], "Adenylosuccinate lyase")
        self.assertEqual(alignment.annotations["accession"], "1.10.275.10/FF/000026")
        self.assertEqual(alignment.annotations["type"], "FunFam")
        self.assertEqual(len(alignment.annotations["database references"]), 2)
        self.assertEqual(
            alignment.annotations["database references"][0], {"reference": "CATH: v4.3"}
        )
        self.assertEqual(
            alignment.annotations["database references"][1],
            {"reference": "DOPS: 0.000"},
        )
        self.assertEqual(alignment.sequences[0].annotations["accession"], "Q9X0I0")
        self.assertEqual(
            alignment.sequences[0].annotations["organism"], "Thermotoga maritima MSB8"
        )
        self.assertEqual(alignment.sequences[0].description, "Adenylosuccinate lyase")
        self.assertEqual(len(alignment.sequences[0].dbxrefs), 3)
        self.assertEqual(alignment.sequences[0].dbxrefs[0], "CATH; 1c3c; B:2-92;")
        self.assertEqual(
            alignment.sequences[0].dbxrefs[1],
            "ORG; Bacteria; Thermotogae; Thermotogae; Thermotogales; Thermotogaceae; Thermotoga; Thermotoga maritima;",
        )
        self.assertEqual(alignment.sequences[0].dbxrefs[2], "EC; 4.3.2.2;")
        self.assertEqual(alignment.sequences[1].annotations["accession"], "Q9X0I0")
        self.assertEqual(
            alignment.sequences[1].annotations["organism"], "Thermotoga maritima MSB8"
        )
        self.assertEqual(alignment.sequences[1].description, "Adenylosuccinate lyase")
        self.assertEqual(len(alignment.sequences[1].dbxrefs), 3)
        self.assertEqual(alignment.sequences[1].dbxrefs[0], "CATH; 1c3c; A:2-92;")
        self.assertEqual(
            alignment.sequences[1].dbxrefs[1],
            "ORG; Bacteria; Thermotogae; Thermotogae; Thermotogales; Thermotogaceae; Thermotoga; Thermotoga maritima;",
        )
        self.assertEqual(alignment.sequences[1].dbxrefs[2], "EC; 4.3.2.2;")
        self.assertEqual(alignment.sequences[2].annotations["accession"], "Q9X0I0")
        self.assertEqual(
            alignment.sequences[2].annotations["organism"], "Thermotoga maritima MSB8"
        )
        self.assertEqual(alignment.sequences[2].description, "Adenylosuccinate lyase")
        self.assertEqual(len(alignment.sequences[2].dbxrefs), 2)
        self.assertEqual(
            alignment.sequences[2].dbxrefs[0],
            "ORG; Bacteria; Thermotogae; Thermotogae; Thermotogales; Thermotogaceae; Thermotoga; Thermotoga maritima;",
        )
        self.assertEqual(alignment.sequences[2].dbxrefs[1], "EC; 4.3.2.2;")
        self.assertEqual(alignment.sequences[3].annotations["accession"], "G4FEQ2")
        self.assertEqual(
            alignment.sequences[3].annotations["organism"], "Thermotoga maritima MSB8"
        )
        self.assertEqual(alignment.sequences[3].description, "Adenylosuccinate lyase")
        self.assertEqual(len(alignment.sequences[3].dbxrefs), 2)
        self.assertEqual(
            alignment.sequences[3].dbxrefs[0],
            "ORG; Bacteria; Thermotogae; Thermotogae; Thermotogales; Thermotogaceae; Thermotoga; Thermotoga maritima;",
        )
        self.assertEqual(alignment.sequences[3].dbxrefs[1], "EC; 4.3.2.2;")
        self.assertEqual(len(alignment.sequences), 4)
        self.assertEqual(alignment.sequences[0].id, "1c3cB01/1-91")
        self.assertEqual(alignment.sequences[1].id, "1c3cA01/1-91")
        self.assertEqual(alignment.sequences[2].id, "Q9X0I0/2-92")
        self.assertEqual(alignment.sequences[3].id, "G4FEQ2/2-92")
        self.assertEqual(
            alignment.sequences[0].seq,
            "VERYSLSPMKDLWTEEAKYRRWLEVELAVTRAYEELGMIPKGVTERIRNNAKIDVELFKKIEEKTNHDVVAFVEGIGSMIGEDSRFFHYGL",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "VERYSLSPMKDLWTEEAKYRRWLEVELAVTRAYEELGMIPKGVTERIRNNAKIDVELFKKIEEKTNHDVVAFVEGIGSMIGEDSRFFHYGL",
        )
        self.assertEqual(
            alignment.sequences[2].seq,
            "VERYSLSPMKDLWTEEAKYRRWLEVELAVTRAYEELGMIPKGVTERIRNNAKIDVELFKKIEEKTNHDVVAFVEGIGSMIGEDSRFFHYGL",
        )
        self.assertEqual(
            alignment.sequences[3].seq,
            "VERYSLSPMKDLWTEEAKYRRWLEVELAVTRAYEELGMIPKGVTERIRNNAKIDVELFKKIEEKTNHDVVAFVEGIGSMIGEDSRFFHYGL",
        )
        self.assertEqual(
            alignment[0],
            "VERYSLSPMKDLWTEEAKYRRWLEVELAVTRAYEELGMIPKGVTERIRNNAKIDVELFKKIEEKTNHDVVAFVEGIGSMIGEDSRFFHYGL",
        )
        self.assertEqual(
            alignment[1],
            "VERYSLSPMKDLWTEEAKYRRWLEVELAVTRAYEELGMIPKGVTERIRNNAKIDVELFKKIEEKTNHDVVAFVEGIGSMIGEDSRFFHYGL",
        )
        self.assertEqual(
            alignment[2],
            "VERYSLSPMKDLWTEEAKYRRWLEVELAVTRAYEELGMIPKGVTERIRNNAKIDVELFKKIEEKTNHDVVAFVEGIGSMIGEDSRFFHYGL",
        )
        self.assertEqual(
            alignment[3],
            "VERYSLSPMKDLWTEEAKYRRWLEVELAVTRAYEELGMIPKGVTERIRNNAKIDVELFKKIEEKTNHDVVAFVEGIGSMIGEDSRFFHYGL",
        )
        self.assertEqual(
            alignment.sequences[0].letter_annotations["Catalytic Site Atlas"],
            "__________________________________________________________________0________________________",
        )
        self.assertEqual(
            alignment.column_annotations["consensus score"],
            "0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",
        )
        self.assertEqual(
            alignment.column_annotations["consensus score 70"],
            "___________________________________________________________________________________________",
        )
        self.assertEqual(
            alignment.column_annotations["consensus score 80"],
            "___________________________________________________________________________________________",
        )
        self.assertEqual(
            alignment.column_annotations["consensus score 90"],
            "___________________________________________________________________________________________",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[0, 91], [0, 91], [0, 91], [0, 91]])
            )
        )

    def test_reading_writing_alignments_globins45(self):
        """Test parsing hmmalign output."""
        # File generated by running
        # hmmalign -o globins45.ali globins4.hmm globins45.fa
        # in the HMMER 3.3.2 tutorial
        path = "Stockholm/globins45.ali"
        alignments = stockholm.AlignmentIterator(path)
        alignment = next(alignments)
        self.assertRaises(StopIteration, next, alignments)
        self.check_alignment_globins45(alignment)
        stream = StringIO()
        writer = stockholm.AlignmentWriter(stream)
        alignments = [alignment]
        n = writer.write_file(alignments, mincount=1, maxcount=1)
        self.assertEqual(n, 1)
        stream.seek(0)
        alignments = stockholm.AlignmentIterator(stream)
        alignment = next(alignments)
        stream.close()
        self.check_alignment_globins45(alignment)

    def test_reading_writing_alignments_pfam1(self):
        """Test parsing Pfam record 120_Rick_ant."""
        path = "Stockholm/pfam1.seed.txt"
        alignments = stockholm.AlignmentIterator(path)
        alignment = next(alignments)
        self.assertRaises(StopIteration, next, alignments)
        self.check_alignment_pfam1(alignment)
        stream = StringIO()
        writer = stockholm.AlignmentWriter(stream)
        alignments = [alignment]
        n = writer.write_file(alignments, mincount=1, maxcount=1)
        self.assertEqual(n, 1)
        stream.seek(0)
        alignments = stockholm.AlignmentIterator(stream)
        alignment = next(alignments)
        stream.close()
        self.check_alignment_pfam1(alignment)

    def test_reading_writing_alignments_pfam2(self):
        """Test parsing Pfam record 7kD_DNA_binding."""
        path = "Stockholm/pfam2.seed.txt"
        alignments = stockholm.AlignmentIterator(path)
        alignment = next(alignments)
        self.assertRaises(StopIteration, next, alignments)
        self.check_alignment_pfam2(alignment)
        stream = StringIO()
        writer = stockholm.AlignmentWriter(stream)
        alignments = [alignment]
        n = writer.write_file(alignments, mincount=1, maxcount=1)
        self.assertEqual(n, 1)
        stream.seek(0)
        alignments = stockholm.AlignmentIterator(stream)
        alignment = next(alignments)
        stream.close()
        self.check_alignment_pfam2(alignment)

    def test_reading_writing_alignments_pfam3(self):
        """Test parsing Pfam record 12TM_1."""
        path = "Stockholm/pfam3.seed.txt"
        alignments = stockholm.AlignmentIterator(path)
        alignment = next(alignments)
        self.assertRaises(StopIteration, next, alignments)
        self.check_alignment_pfam3(alignment)
        stream = StringIO()
        writer = stockholm.AlignmentWriter(stream)
        alignments = [alignment]
        n = writer.write_file(alignments, mincount=1, maxcount=1)
        self.assertEqual(n, 1)
        stream.seek(0)
        alignments = stockholm.AlignmentIterator(stream)
        alignment = next(alignments)
        stream.close()
        self.check_alignment_pfam3(alignment)

    def test_reading_writing_alignments_pfam4(self):
        """Test parsing Pfam record 3Beta_HSD."""
        path = "Stockholm/pfam4.seed.txt"
        with open(path, encoding="UTF-8") as stream:
            # encoding depends on locale by default
            alignments = stockholm.AlignmentIterator(stream)
            alignment = next(alignments)
            self.assertRaises(StopIteration, next, alignments)
        self.check_alignment_pfam4(alignment)
        stream = StringIO()
        writer = stockholm.AlignmentWriter(stream)
        alignments = [alignment]
        n = writer.write_file(alignments, mincount=1, maxcount=1)
        self.assertEqual(n, 1)
        stream.seek(0)
        alignments = stockholm.AlignmentIterator(stream)
        alignment = next(alignments)
        stream.close()
        self.check_alignment_pfam4(alignment)

    def test_reading_writing_alignments_pfam5(self):
        """Test parsing Pfam record ArsP_1."""
        path = "Stockholm/pfam5.seed.txt"
        alignments = stockholm.AlignmentIterator(path)
        alignment = next(alignments)
        self.assertRaises(StopIteration, next, alignments)
        self.check_alignment_pfam5(alignment)
        stream = StringIO()
        writer = stockholm.AlignmentWriter(stream)
        alignments = [alignment]
        n = writer.write_file(alignments, mincount=1, maxcount=1)
        self.assertEqual(n, 1)
        stream.seek(0)
        alignments = stockholm.AlignmentIterator(stream)
        alignment = next(alignments)
        stream.close()
        self.check_alignment_pfam5(alignment)

    def test_reading_writing_alignments_pfam6(self):
        """Test parsing Pfam record COX2_TM."""
        path = "Stockholm/pfam6.seed.txt"
        alignments = stockholm.AlignmentIterator(path)
        alignment = next(alignments)
        self.assertRaises(StopIteration, next, alignments)
        self.check_alignment_pfam6(alignment)
        stream = StringIO()
        writer = stockholm.AlignmentWriter(stream)
        alignments = [alignment]
        n = writer.write_file(alignments, mincount=1, maxcount=1)
        self.assertEqual(n, 1)
        stream.seek(0)
        alignments = stockholm.AlignmentIterator(stream)
        alignment = next(alignments)
        stream.close()
        self.check_alignment_pfam6(alignment)

    def test_reading_writing_alignments_pfam7(self):
        """Test parsing Pfam record Alpha_E1_glycop."""
        path = "Stockholm/pfam7.seed.txt"
        alignments = stockholm.AlignmentIterator(path)
        alignment = next(alignments)
        self.assertRaises(StopIteration, next, alignments)
        self.check_alignment_pfam7(alignment)
        stream = StringIO()
        writer = stockholm.AlignmentWriter(stream)
        alignments = [alignment]
        n = writer.write_file(alignments, mincount=1, maxcount=1)
        self.assertEqual(n, 1)
        stream.seek(0)
        alignments = stockholm.AlignmentIterator(stream)
        alignment = next(alignments)
        stream.close()
        self.check_alignment_pfam7(alignment)

    def test_reading_writing_alignments_pfam8(self):
        """Test parsing Pfam record Cyclin_N."""
        path = "Stockholm/pfam8.seed.txt"
        alignments = stockholm.AlignmentIterator(path)
        alignment = next(alignments)
        self.assertRaises(StopIteration, next, alignments)
        self.check_alignment_pfam8(alignment)
        stream = StringIO()
        writer = stockholm.AlignmentWriter(stream)
        alignments = [alignment]
        n = writer.write_file(alignments, mincount=1, maxcount=1)
        self.assertEqual(n, 1)
        stream.seek(0)
        alignments = stockholm.AlignmentIterator(stream)
        alignment = next(alignments)
        stream.close()
        self.check_alignment_pfam8(alignment)

    def test_reading_writing_alignments_pfam9(self):
        """Test parsing Pfam record SH3_11."""
        path = "Stockholm/pfam9.seed.txt"
        alignments = stockholm.AlignmentIterator(path)
        alignment = next(alignments)
        self.assertRaises(StopIteration, next, alignments)
        self.check_alignment_pfam9(alignment)
        stream = StringIO()
        writer = stockholm.AlignmentWriter(stream)
        alignments = [alignment]
        n = writer.write_file(alignments, mincount=1, maxcount=1)
        self.assertEqual(n, 1)
        stream.seek(0)
        alignments = stockholm.AlignmentIterator(stream)
        alignment = next(alignments)
        stream.close()
        self.check_alignment_pfam9(alignment)

    def test_reading_writing_alignments_rfam1(self):
        """Test parsing Rfam record BTnc005."""
        path = "Stockholm/rfam1.seed.txt"
        with open(path, encoding="UTF-8") as stream:
            alignments = stockholm.AlignmentIterator(stream)
            alignment = next(alignments)
            self.assertRaises(StopIteration, next, alignments)
        self.check_alignment_rfam1(alignment)
        stream = StringIO()
        writer = stockholm.AlignmentWriter(stream)
        alignments = [alignment]
        n = writer.write_file(alignments, mincount=1, maxcount=1)
        self.assertEqual(n, 1)
        stream.seek(0)
        alignments = stockholm.AlignmentIterator(stream)
        alignment = next(alignments)
        stream.close()
        self.check_alignment_rfam1(alignment)

    def test_reading_writing_alignments_rfam2(self):
        """Test parsing Rfam record SraC_RyeA."""
        path = "Stockholm/rfam2.seed.txt"
        alignments = stockholm.AlignmentIterator(path)
        alignment = next(alignments)
        self.assertRaises(StopIteration, next, alignments)
        self.check_alignment_rfam2(alignment)
        stream = StringIO()
        writer = stockholm.AlignmentWriter(stream)
        alignments = [alignment]
        n = writer.write_file(alignments, mincount=1, maxcount=1)
        self.assertEqual(n, 1)
        stream.seek(0)
        alignments = stockholm.AlignmentIterator(stream)
        alignment = next(alignments)
        stream.close()
        self.check_alignment_rfam2(alignment)

    def test_reading_writing_alignments_rfam3(self):
        """Test parsing Rfam record McaS."""
        path = "Stockholm/rfam3.seed.txt"
        alignments = stockholm.AlignmentIterator(path)
        alignment = next(alignments)
        self.assertRaises(StopIteration, next, alignments)
        self.check_alignment_rfam3(alignment)
        stream = StringIO()
        writer = stockholm.AlignmentWriter(stream)
        alignments = [alignment]
        n = writer.write_file(alignments, mincount=1, maxcount=1)
        self.assertEqual(n, 1)
        stream.seek(0)
        alignments = stockholm.AlignmentIterator(stream)
        alignment = next(alignments)
        stream.close()
        self.check_alignment_rfam3(alignment)

    def test_reading_writing_alignments_rfam4(self):
        """Test parsing Rfam record IRES_KSHV."""
        path = "Stockholm/rfam4.seed.txt"
        alignments = stockholm.AlignmentIterator(path)
        alignment = next(alignments)
        self.assertRaises(StopIteration, next, alignments)
        self.check_alignment_rfam4(alignment)
        stream = StringIO()
        writer = stockholm.AlignmentWriter(stream)
        alignments = [alignment]
        n = writer.write_file(alignments, mincount=1, maxcount=1)
        self.assertEqual(n, 1)
        stream.seek(0)
        alignments = stockholm.AlignmentIterator(stream)
        alignment = next(alignments)
        stream.close()
        self.check_alignment_rfam4(alignment)

    def test_reading_writing_alignments_rfam5(self):
        """Test parsing Rfam record BMV3_UPD-PK3."""
        path = "Stockholm/rfam5.seed.txt"
        alignments = stockholm.AlignmentIterator(path)
        alignment = next(alignments)
        self.assertRaises(StopIteration, next, alignments)
        self.check_alignment_rfam5(alignment)
        stream = StringIO()
        writer = stockholm.AlignmentWriter(stream)
        alignments = [alignment]
        n = writer.write_file(alignments, mincount=1, maxcount=1)
        self.assertEqual(n, 1)
        stream.seek(0)
        alignments = stockholm.AlignmentIterator(stream)
        alignment = next(alignments)
        stream.close()
        self.check_alignment_rfam5(alignment)

    def test_reading_alignments_cath1(self):
        """Test parsing CATH record 3.30.160.60/FF/004774."""
        path = "Stockholm/cath1.sth"
        alignments = stockholm.AlignmentIterator(path)
        alignment = next(alignments)
        self.assertRaises(StopIteration, next, alignments)
        self.check_alignment_cath1(alignment)
        stream = StringIO()
        writer = stockholm.AlignmentWriter(stream)
        alignments = [alignment]
        n = writer.write_file(alignments, mincount=1, maxcount=1)
        self.assertEqual(n, 1)
        stream.seek(0)
        alignments = stockholm.AlignmentIterator(stream)
        alignment = next(alignments)
        stream.close()
        self.check_alignment_cath1(alignment)

    def test_reading_alignments_cath2(self):
        """Test parsing CATH record 2.105.10.10/FF/000002."""
        path = "Stockholm/cath2.sth"
        alignments = stockholm.AlignmentIterator(path)
        alignment = next(alignments)
        self.assertRaises(StopIteration, next, alignments)
        self.check_alignment_cath2(alignment)
        stream = StringIO()
        writer = stockholm.AlignmentWriter(stream)
        alignments = [alignment]
        n = writer.write_file(alignments, mincount=1, maxcount=1)
        self.assertEqual(n, 1)
        stream.seek(0)
        alignments = stockholm.AlignmentIterator(stream)
        alignment = next(alignments)
        stream.close()
        self.check_alignment_cath2(alignment)

    def test_reading_alignments_cath3(self):
        """Test parsing CATH record 1.10.275.10/FF/000026."""
        path = "Stockholm/cath3.sth"
        alignments = stockholm.AlignmentIterator(path)
        alignment = next(alignments)
        self.assertRaises(StopIteration, next, alignments)
        self.check_alignment_cath3(alignment)
        stream = StringIO()
        writer = stockholm.AlignmentWriter(stream)
        alignments = [alignment]
        n = writer.write_file(alignments, mincount=1, maxcount=1)
        self.assertEqual(n, 1)
        stream.seek(0)
        alignments = stockholm.AlignmentIterator(stream)
        alignment = next(alignments)
        stream.close()
        self.check_alignment_cath3(alignment)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
