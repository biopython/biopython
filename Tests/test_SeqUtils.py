# Copyright 2003 by Iddo Friedberg.  All rights reserved.
# Copyright 2007-2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for SeqUtils module."""
import os
import unittest

from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils import GC_skew
from Bio.SeqUtils import seq1
from Bio.SeqUtils import seq3
from Bio.SeqUtils.CheckSum import crc32
from Bio.SeqUtils.CheckSum import crc64
from Bio.SeqUtils.CheckSum import gcg
from Bio.SeqUtils.CheckSum import seguid
from Bio.SeqUtils import CodonAdaptationIndex
from Bio.SeqUtils.lcc import lcc_mult
from Bio.SeqUtils.lcc import lcc_simp


class SeqUtilsTests(unittest.TestCase):
    # Example of crc64 collision from Sebastian Bassi using the
    # immunoglobulin lambda light chain variable region from Homo sapiens
    # Both sequences share the same CRC64 checksum: 44CAAD88706CC153
    str_light_chain_one = (
        "QSALTQPASVSGSPGQSITISCTGTSSDVGSYNLVSWYQQHPGKAPKLMIYEGSKRPSGV"
        "SNRFSGSKSGNTASLTISGLQAEDEADYYCSSYAGSSTLVFGGGTKLTVL"
    )
    str_light_chain_two = (
        "QSALTQPASVSGSPGQSITISCTGTSSDVGSYNLVSWYQQHPGKAPKLMIYEGSKRPSGV"
        "SNRFSGSKSGNTASLTISGLQAEDEADYYCCSYAGSSTWVFGGGTKLTVL"
    )

    def test_codon_adaptation_index_initialization(self):
        """Test Codon Adaptation Index (CAI) initialization from sequences."""
        # We need CDS sequences to count the codon usage...
        dna_filename = "GenBank/NC_005816.gb"
        record = SeqIO.read(dna_filename, "genbank")
        records = []
        for feature in record.features:
            if feature.type == "CDS" and len(feature.location.parts) == 1:
                start = feature.location.start
                end = feature.location.end
                table = int(feature.qualifiers["transl_table"][0])
                if feature.location.strand == -1:
                    seq = record.seq[start:end].reverse_complement()
                else:
                    seq = record.seq[start:end]
                # Double check we have the CDS sequence expected
                # TODO - Use any cds_start option if/when added to deal with the met
                a = "M" + seq[3:].translate(table)
                b = feature.qualifiers["translation"][0] + "*"
                self.assertEqual(a, b)
                records.append(
                    SeqRecord(
                        seq,
                        id=feature.qualifiers["protein_id"][0],
                        description=feature.qualifiers["product"][0],
                    )
                )

        cai = CodonAdaptationIndex(records)
        # Now check codon usage index (CAI) using this species
        self.assertEqual(
            record.annotations["source"], "Yersinia pestis biovar Microtus str. 91001"
        )
        value = cai.calculate("ATGCGTATCGATCGCGATACGATTAGGCGGATG")
        self.assertAlmostEqual(value, 0.70246, places=5)
        optimized_sequence = cai.optimize(
            "ATGCGTATCGATCGCGATACGATTAGGCGGATG", strict=False
        )
        optimized_value = cai.calculate(optimized_sequence)
        self.assertEqual(optimized_value, 1.0)
        aa_initial = Seq("ATGCGTATCGATCGCGATACGATTAGGCGGATG").translate()
        aa_optimized = optimized_sequence.translate()
        self.assertEqual(aa_initial, aa_optimized)
        with self.assertRaises(KeyError):
            cai.optimize("CAU", "protein", strict=False)
        self.maxDiff = None
        self.assertEqual(
            str(cai),
            """\
AAA	1.000
AAC	0.385
AAG	0.344
AAT	1.000
ACA	1.000
ACC	0.553
ACG	0.319
ACT	0.447
AGA	0.595
AGC	0.967
AGG	0.297
AGT	1.000
ATA	0.581
ATC	0.930
ATG	1.000
ATT	1.000
CAA	0.381
CAC	0.581
CAG	1.000
CAT	1.000
CCA	0.500
CCC	0.500
CCG	1.000
CCT	0.767
CGA	0.568
CGC	0.919
CGG	0.514
CGT	1.000
CTA	0.106
CTC	0.379
CTG	1.000
CTT	0.424
GAA	1.000
GAC	0.633
GAG	0.506
GAT	1.000
GCA	1.000
GCC	0.617
GCG	0.532
GCT	0.809
GGA	1.000
GGC	0.525
GGG	0.575
GGT	0.950
GTA	0.500
GTC	0.618
GTG	0.971
GTT	1.000
TAA	1.000
TAC	0.434
TAG	0.062
TAT	1.000
TCA	1.000
TCC	0.533
TCG	0.233
TCT	0.967
TGA	0.250
TGC	1.000
TGG	1.000
TGT	0.750
TTA	0.455
TTC	1.000
TTG	0.212
TTT	0.886
""",
        )

    def test_codon_adaptation_index_calculation(self):
        """Test Codon Adaptation Index (CAI) calculation for an mRNA."""
        cai = CodonAdaptationIndex([])
        # Use the Codon Adaption Index for E. coli, precalculated by
        # Sharp and Li (Nucleic Acids Res. 1987 Feb 11;15(3):1281-95), Table 1.
        cai["TTT"] = 0.296  # Phe
        cai["TTC"] = 1.000  # Phe
        cai["TTA"] = 0.020  # Leu
        cai["TTG"] = 0.020  # Leu
        cai["CTT"] = 0.042  # Leu
        cai["CTC"] = 0.037  # Leu
        cai["CTA"] = 0.007  # Leu
        cai["CTG"] = 1.000  # Leu
        cai["ATT"] = 0.185  # Ile
        cai["ATC"] = 1.000  # Ile
        cai["ATA"] = 0.003  # Ile
        cai["ATG"] = 1.000  # Met
        cai["GTT"] = 1.000  # Val
        cai["GTC"] = 0.066  # Val
        cai["GTA"] = 0.495  # Val
        cai["GTG"] = 0.221  # Val
        cai["TAT"] = 0.239  # Tyr
        cai["TAC"] = 1.000  # Tyr
        cai["CAT"] = 0.291  # His
        cai["CAC"] = 1.000  # His
        cai["CAA"] = 0.124  # Gln
        cai["CAG"] = 1.000  # Gln
        cai["AAT"] = 0.051  # Asn
        cai["AAC"] = 1.000  # Asn
        cai["AAA"] = 1.000  # Lys
        cai["AAG"] = 0.253  # Lys
        cai["GAT"] = 0.434  # Asp
        cai["GAC"] = 1.000  # Asp
        cai["GAA"] = 1.000  # Glu
        cai["GAG"] = 0.259  # Glu
        cai["TCT"] = 1.000  # Ser
        cai["TCC"] = 0.744  # Ser
        cai["TCA"] = 0.077  # Ser
        cai["TCG"] = 0.017  # Ser
        cai["CCT"] = 0.070  # Pro
        cai["CCC"] = 0.012  # Pro
        cai["CCA"] = 0.135  # Pro
        cai["CCG"] = 1.000  # Pro
        cai["ACT"] = 0.965  # Thr
        cai["ACC"] = 1.000  # Thr
        cai["ACA"] = 0.076  # Thr
        cai["ACG"] = 0.099  # Thr
        cai["GCT"] = 1.000  # Ala
        cai["GCC"] = 0.122  # Ala
        cai["GCA"] = 0.586  # Ala
        cai["GCG"] = 0.424  # Ala
        cai["TGT"] = 0.500  # Cys
        cai["TGC"] = 1.000  # Cys
        cai["TGG"] = 1.000  # Trp
        cai["CGT"] = 1.000  # Arg
        cai["CGC"] = 0.356  # Arg
        cai["CGA"] = 0.004  # Arg
        cai["CGG"] = 0.004  # Arg
        cai["AGT"] = 0.085  # Ser
        cai["AGC"] = 0.410  # Ser
        cai["AGA"] = 0.004  # Arg
        cai["AGG"] = 0.002  # Arg
        cai["GGT"] = 1.000  # Gly
        cai["GGC"] = 0.724  # Gly
        cai["GGA"] = 0.010  # Gly
        cai["GGG"] = 0.019  # Gly
        # Now calculate the CAI for the genes listed in Table 2 of
        # Sharp and Li (Nucleic Acids Res. 1987 Feb 11;15(3):1281-95).
        rpsU = Seq(
            "CCGGTAATTAAAGTACGTGAAAACGAGCCGTTCGACGTAGCTCTGCGTCGCTTCAAGCGTTCCTGCGAAAAAGCAGGTGTTCTGGCGGAAGTTCGTCGTCGTGAGTTCTATGAAAAACCGACTACCGAACGTAAGCGCGCTAAAGCTTCTGCAGTGAAACGTCACGCGAAGAAACTGGCTCGCGAAAACGCACGCCGCACTCGTCTGTAC"
        )
        self.assertAlmostEqual(cai.calculate(rpsU), 0.726, places=3)
        rpoD = Seq(
            "ATGGAGCAAAACCCGCAGTCACAGCTGAAACTTCTTGTCACCCGTGGTAAGGAGCAAGGCTATCTGACCTATGCCGAGGTCAATGACCATCTGCCGGAAGATATCGTCGATTCAGATCAGATCGAAGACATCATCCAAATGATCAACGACATGGGCATTCAGGTGATGGAAGAAGCACCGGATGCCGATGATCTGATGCTGGCTGAAAACACCGCGGACGAAGATGCTGCCGAAGCCGCCGCGCAGGTGCTTTCCAGCGTGGAATCTGAAATCGGGCGCACGACTGACCCGGTACGCATGTACATGCGTGAAATGGGCACCGTTGAACTGTTGACCCGCGAAGGCGAAATTGACATCGCTAAGCGTATTGAAGACGGGATCAACCAGGTTCAATGCTCCGTTGCTGAATATCCGGAAGCGATCACCTATCTGCTGGAACAGTACGATCGTGTTGAAGCAGAAGAAGCGCGTCTGTCCGATCTGATCACCGGCTTTGTTGACCCGAACGCAGAAGAAGATCTGGCACCTACCGCCACTCACGTCGGTTCTGAGCTTTCCCAGGAAGATCTGGACGATGACGAAGATGAAGACGAAGAAGATGGCGATGACGACAGCGCCGATGATGACAACAGCATCGACCCGGAACTGGCTCGCGAAAAATTTGCGGAACTACGCGCTCAGTACGTTGTAACGCGTGACACCATCAAAGCGAAAGGTCGCAGTCACGCTACCGCTCAGGAAGAGATCCTGAAACTGTCTGAAGTATTCAAACAGTTCCGCCTGGTGCCGAAGCAGTTTGACTACCTGGTCAACAGCATGCGCGTCATGATGGACCGCGTTCGTACGCAAGAACGTCTGATCATGAAGCTCTGCGTTGAGCAGTGCAAAATGCCGAAGAAAAACTTCATTACCCTGTTTACCGGCAACGAAACCAGCGATACCTGGTTCAACGCGGCAATTGCGATGAACAAGCCGTGGTCGGAAAAACTGCACGATGTCTCTGAAGAAGTGCATCGCGCCCTGCAAAAACTGCAGCAGATTGAAGAAGAAACCGGCCTGACCATCGAGCAGGTTAAAGATATCAACCGTCGTATGTCCATCGGTGAAGCGAAAGCCCGCCGTGCGAAGAAAGAGATGGTTGAAGCGAACTTACGTCTGGTTATTTCTATCGCTAAGAAATACACCAACCGTGGCTTGCAGTTCCTTGACCTGATTCAGGAAGGCAACATCGGTCTGATGAAAGCGGTTGATAAATTCGAATACCGCCGTGGTTACAAGTTCTCCACCTACGCAACCTGGTGGATCCGTCAGGCGATCACCCGCTCTATCGCGGATCAGGCGCGCACCATCCGTATTCCGGTGCATATGATTGAGACCATCAACAAGCTCAACCGTATTTCTCGCCAGATGCTGCAAGAGATGGGCCGTGAACCGACGCCGGAAGAACTGGCTGAACGTATGCTGATGCCGGAAGACAAGATCCGCAAAGTGCTGAAGATCGCCAAAGAGCCAATCTCCATGGAAACGCCGATCGGTGATGATGAAGATTCGCATCTGGGGGATTTCATCGAGGATACCACCCTCGAGCTGCCGCTGGATTCTGCGACCACCGAAAGCCTGCGTGCGGCAACGCACGACGTGCTGGCTGGCCTGACCGCGCGTGAAGCAAAAGTTCTGCGTATGCGTTTCGGTATCGATATGAACACCGACTACACGCTGGAAGAAGTGGGTAAACAGTTCGACGTTACCCGCGAACGTATCCGTCAGATCGAAGCGAAGGCGCTGCGCAAACTGCGTCACCCGAGCCGTTCTGAAGTGCTGCGTAGCTTCCTGGACGAT"
        )
        self.assertAlmostEqual(cai.calculate(rpoD), 0.582, places=2)
        dnaG = "ATGGCTGGACGAATCCCACGCGTATTCATTAATGATCTGCTGGCACGCACTGACATCGTCGATCTGATCGATGCCCGTGTGAAGCTGAAAAAGCAGGGCAAGAATTTCCACGCGTGTTGTCCATTCCACAACGAGAAAACCCCGTCCTTCACCGTTAACGGTGAGAAACAGTTTTACCACTGCTTTGGATGTGGCGCGCACGGCAACGCGATCGACTTCCTGATGAACTACGACAAGCTCGAGTTCGTCGAAACGGTCGAAGAGCTGGCAGCAATGCACAATCTTGAAGTGCCATTTGAAGCAGGCAGCGGCCCCAGCCAGATCGAGCGCCATCAGAGGCAAACGCTTTATCAGTTGATGGACGGTCTGAATACGTTTTACCAACAATCTTTACAACAACCTGTTGCCACGTCTGCGCGCCAGTATCTGGAAAAACGCGGATTAAGCCACGAGGTTATCGCTCGCTTTGCGATTGGTTTTGCGCCCCCCGGCTGGGACAACGTCCTGAAGCGGTTTGGCGGCAATCCAGAAAATCGCCAGTCATTGATTGATGCGGGGATGTTGGTCACTAACGATCAGGGACGCAGTTACGATCGTTTCCGCGAGCGGGTGATGTTCCCCATTCGCGATAAACGCGGTCGGGTGATTGGTTTTGGCGGGCGCGTGCTGGGCAACGATACCCCCAAATACCTGAACTCGCCGGAAACAGACATTTTCCATAAAGGCCGCCAGCTTTACGGTCTTTATGAAGCGCAGCAGGATAACGCTGAACCCAATCGTCTGCTTGTGGTCGAAGGCTATATGGACGTGGTGGCGCTGGCGCAATACGGCATTAATTACGCCGTTGCGTCGTTAGGTACGTCAACCACCGCCGATCACATACAACTGTTGTTCCGCGCGACCAACAATGTCATTTGCTGTTATGACGGCGACCGTGCAGGCCGCGATGCCGCCTGGCGAGCGCTGGAAACGGCGCTGCCTTACATGACAGACGGCCGTCAGCTACGCTTTATGTTTTTGCCTGATGGCGAAGACCCTGACACGCTAGTACGAAAAGAAGGTAAAGAAGCGTTTGAAGCGCGGATGGAGCAGGCGATGCCACTCTCCGCATTTCTGTTTAACAGTCTGATGCCGCAAGTTGATCTGAGTACCCCTGACGGGCGCGCACGTTTGAGTACGCTGGCACTACCATTGATATCGCAAGTGCCGGGCGAAACGCTGCGAATATATCTTCGTCAGGAATTAGGCAACAAATTAGGCATACTTGATGACAGCCAGCTTGAACGATTAATGCCAAAAGCGGCAGAGAGCGGCGTTTCTCGCCCTGTTCCGCAGCTAAAACGCACGACCATGCGTATACTTATAGGGTTGCTGGTGCAAAATCCAGAATTAGCGACGTTGGTCCCGCCGCTTGAGAATCTGGATGAAAATAAGCTCCCTGGACTTGGCTTATTCAGAGAACTGGTCAACACTTGTCTCTCCCAGCCAGGTCTGACCACCGGGCAACTTTTAGAGCACTATCGTGGTACAAATAATGCTGCCACCCTTGAAAAACTGTCGATGTGGGACGATATAGCAGATAAGAATATTGCTGAGCAAACCTTCACCGACTCACTCAACCATATGTTTGATTCGCTGCTTGAACTGCGCCAGGAAGAGTTAATCGCTCGTGAGCGCACGCATGGTTTAAGCAACGAAGAACGCCTGGAGCTCTGGACATTAAACCAGGAGCTGGCGAAAAAG"
        self.assertAlmostEqual(cai.calculate(dnaG), 0.271, places=3)
        lacI = "GTGAAACCAGTAACGTTATACGATGTCGCAGAGTATGCCGGTGTCTCTTATCAGACCGTTTCCCGCGTGGTGAACCAGGCCAGCCACGTTTCTGCGAAAACGCGGGAAAAAGTGGAAGCGGCGATGGCGGAGCTGAATTACATTCCCAACCGCGTGGCACAACAACTGGCGGGCAAACAGTCGTTGCTGATTGGCGTTGCCACCTCCAGTCTGGCCCTGCACGCGCCGTCGCAAATTGTCGCGGCGATTAAATCTCGCGCCGATCAACTGGGTGCCAGCGTGGTGGTGTCGATGGTAGAACGAAGCGGCGTCGAAGCCTGTAAAGCGGCGGTGCACAATCTTCTCGCGCAACGCGTCAGTGGGCTGATCATTAACTATCCGCTGGATGACCAGGATGCCATTGCTGTGGAAGCTGCCTGCACTAATGTTCCGGCGTTATTTCTTGATGTCTCTGACCAGACACCCATCAACAGTATTATTTTCTCCCATGAAGACGGTACGCGACTGGGCGTGGAGCATCTGGTCGCATTGGGTCACCAGCAAATCGCGCTGTTAGCGGGCCCATTAAGTTCTGTCTCGGCGCGTCTGCGTCTGGCTGGCTGGCATAAATATCTCACTCGCAATCAAATTCAGCCGATAGCGGAACGGGAAGGCGACTGGAGTGCCATGTCCGGTTTTCAACAAACCATGCAAATGCTGAATGAGGGCATCGTTCCCACTGCGATGCTGGTTGCCAACGATCAGATGGCGCTGGGCGCAATGCGCGCCATTACCGAGTCCGGGCTGCGCGTTGGTGCGGATATCTCGGTAGTGGGATACGACGATACCGAAGACAGCTCATGTTATATCCCGCCGTTAACCACCATCAAACAGGATTTTCGCCTGCTGGGGCAAACCAGCGTGGACCGCTTGCTGCAACTCTCTCAGGGCCAGGCGGTGAAGGGCAATCAGCTGTTGCCCGTCTCACTGGTGAAAAGAAAAACCACCCTGGCGCCCAATACGCAAACCGCCTCTCCCCGCGCGTTGGCCGATTCATTAATGCAGCTGGCACGACAGGTTTCCCGACTGGAAAGCGGGCAG"
        self.assertAlmostEqual(cai.calculate(lacI), 0.296, places=2)
        trpR = "ATGGCCCAACAATCACCCTATTCAGCAGCGATGGCAGAACAGCGTCACCAGGAGTGGTTACGTTTTGTCGACCTGCTTAAGAATGCCTACCAAAACGATCTCCATTTACCGTTGTTAAACCTGATGCTGACGCCAGATGAGCGCGAAGCGTTGGGGACTCGCGTGCGTATTGTCGAAGAGCTGTTGCGCGGCGAAATGAGCCAGCGTGAGTTAAAAAATGAACTCGGCGCAGGCATCGCGACGATTACGCGTGGATCTAACAGCCTGAAAGCCGCGCCCGTCGAGCTGCGCCAGTGGCTGGAAGAGGTGTTGCTGAAAAGCGAT"
        self.assertAlmostEqual(cai.calculate(trpR), 0.267, places=2)
        lpp = "ATGAAAGCTACTAAACTGGTACTGGGCGCGGTAATCCTGGGTTCTACTCTGCTGGCAGGTTGCTCCAGCAACGCTAAAATCGATCAGCTGTCTTCTGACGTTCAGACTCTGAACGCTAAAGTTGACCAGCTGAGCAACGACGTGAACGCAATGCGTTCCGACGTTCAGGCTGCTAAAGATGACGCAGCTCGTGCTAACCAGCGTCTGGACAACATGGCTACTAAATACCGCAAG"
        self.assertAlmostEqual(cai.calculate(lpp), 0.849, places=3)

    def test_crc_checksum_collision(self):
        # Explicit testing of crc64 collision:
        self.assertNotEqual(self.str_light_chain_one, self.str_light_chain_two)
        self.assertNotEqual(
            crc32(self.str_light_chain_one), crc32(self.str_light_chain_two)
        )
        self.assertEqual(
            crc64(self.str_light_chain_one), crc64(self.str_light_chain_two)
        )
        self.assertNotEqual(
            gcg(self.str_light_chain_one), gcg(self.str_light_chain_two)
        )
        self.assertNotEqual(
            seguid(self.str_light_chain_one), seguid(self.str_light_chain_two)
        )

    def seq_checksums(
        self,
        seq_str,
        exp_crc32,
        exp_crc64,
        exp_gcg,
        exp_seguid,
        exp_simple_LCC,
        exp_window_LCC,
    ):
        for s in [seq_str, Seq(seq_str), MutableSeq(seq_str)]:
            self.assertEqual(exp_crc32, crc32(s))
            self.assertEqual(exp_crc64, crc64(s))
            self.assertEqual(exp_gcg, gcg(s))
            self.assertEqual(exp_seguid, seguid(s))
            self.assertAlmostEqual(exp_simple_LCC, lcc_simp(s), places=4)
            values = lcc_mult(s, 20)
            self.assertEqual(len(exp_window_LCC), len(values), values)
            for value1, value2 in zip(exp_window_LCC, values):
                self.assertAlmostEqual(value1, value2, places=2)

    def test_checksum1(self):
        self.seq_checksums(
            self.str_light_chain_one,
            2994980265,
            "CRC-44CAAD88706CC153",
            9729,
            "BpBeDdcNUYNsdk46JoJdw7Pd3BI",
            0.5160,
            (
                0.4982,
                0.4794,
                0.4794,
                0.4794,
                0.3241,
                0.2160,
                0.1764,
                0.1764,
                0.1764,
                0.1764,
                0.2657,
                0.2948,
                0.1287,
            ),
        )

    def test_checksum2(self):
        self.seq_checksums(
            self.str_light_chain_two,
            802105214,
            "CRC-44CAAD88706CC153",
            9647,
            "X5XEaayob1nZLOc7eVT9qyczarY",
            0.5343,
            (
                0.4982,
                0.4794,
                0.4794,
                0.4794,
                0.3241,
                0.2160,
                0.1764,
                0.1764,
                0.1764,
                0.1764,
                0.2657,
                0.2948,
                0.1287,
            ),
        )

    def test_checksum3(self):
        self.seq_checksums(
            "ATGCGTATCGATCGCGATACGATTAGGCGGAT",
            817679856,
            "CRC-6234FF451DC6DFC6",
            7959,
            "8WCUbVjBgiRmM10gfR7XJNjbwnE",
            0.9886,
            (
                1.00,
                0.9927,
                0.9927,
                1.00,
                0.9927,
                0.9854,
                0.9927,
                0.9927,
                0.9927,
                0.9794,
                0.9794,
                0.9794,
                0.9794,
            ),
        )

    def test_gc_fraction(self):
        """Tests gc_fraction function."""
        self.assertAlmostEqual(gc_fraction("", "ignore"), 0, places=3)
        self.assertAlmostEqual(gc_fraction("", "weighted"), 0, places=3)
        self.assertAlmostEqual(gc_fraction("", "remove"), 0, places=3)

        seq = "ACGGGCTACCGTATAGGCAAGAGATGATGCCC"
        self.assertAlmostEqual(gc_fraction(seq, "ignore"), 0.5625, places=3)
        self.assertAlmostEqual(gc_fraction(seq, "weighted"), 0.5625, places=3)
        self.assertAlmostEqual(gc_fraction(seq, "remove"), 0.5625, places=3)

        seq = "ACTGSSSS"
        self.assertAlmostEqual(gc_fraction(seq, "ignore"), 0.75, places=3)
        self.assertAlmostEqual(gc_fraction(seq, "weighted"), 0.75, places=3)
        self.assertAlmostEqual(gc_fraction(seq, "remove"), 0.75, places=3)

        # Test ambiguous nucleotide behaviour

        seq = "CCTGNN"
        self.assertAlmostEqual(gc_fraction(seq, "ignore"), 0.5, places=3)
        self.assertAlmostEqual(gc_fraction(seq, "weighted"), 0.667, places=3)
        self.assertAlmostEqual(gc_fraction(seq, "remove"), 0.75, places=3)

        seq = "GDVV"
        self.assertAlmostEqual(gc_fraction(seq, "ignore"), 0.25, places=3)
        self.assertAlmostEqual(gc_fraction(seq, "weighted"), 0.6667, places=3)
        self.assertAlmostEqual(gc_fraction(seq, "remove"), 1.00, places=3)

        with self.assertRaises(ValueError):
            gc_fraction(seq, "other string")

    def test_GC_skew(self):
        s = "A" * 50
        seq = Seq(s)
        record = SeqRecord(seq)
        self.assertEqual(GC_skew(s)[0], 0)
        self.assertEqual(GC_skew(seq)[0], 0)
        self.assertEqual(GC_skew(record)[0], 0)

    def test_seq1_seq3(self):
        s3 = "MetAlaTyrtrpcysthrLYSLEUILEGlYPrOGlNaSnaLapRoTyRLySSeRHisTrpLysThr"
        s1 = "MAYWCTKLIGPQNAPYKSHWKT"
        self.assertEqual(seq1(s3), s1)
        self.assertEqual(seq3(s1).upper(), s3.upper())
        self.assertEqual(seq1(seq3(s1)), s1)
        self.assertEqual(seq3(seq1(s3)).upper(), s3.upper())

    def test_lcc_simp(self):
        s = "ACGATAGC"
        seq = Seq(s)
        record = SeqRecord(seq)
        self.assertAlmostEqual(lcc_simp(s), 0.9528, places=4)
        self.assertAlmostEqual(lcc_simp(seq), 0.9528, places=4)
        self.assertAlmostEqual(lcc_simp(record), 0.9528, places=4)

    def test_lcc_mult(self):
        s = "ACGATAGC"
        seq = Seq(s)
        record = SeqRecord(seq)
        llc_lst = lcc_mult(s, len(s))
        self.assertEqual(len(llc_lst), 1)
        self.assertAlmostEqual(llc_lst[0], 0.9528, places=4)
        llc_lst = lcc_mult(seq, len(seq))
        self.assertEqual(len(llc_lst), 1)
        self.assertAlmostEqual(llc_lst[0], 0.9528, places=4)
        llc_lst = lcc_mult(record, len(record))
        self.assertEqual(len(llc_lst), 1)
        self.assertAlmostEqual(llc_lst[0], 0.9528, places=4)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
