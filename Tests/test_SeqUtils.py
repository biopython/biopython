# Copyright 2003 by Iddo Friedberg.  All rights reserved.
# Copyright 2007-2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


"""Tests for SeqUtils module."""

import os
import unittest

from Bio import SeqIO
from Bio.Alphabet import single_letter_alphabet
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC, seq1, seq3, GC_skew
from Bio.SeqUtils.lcc import lcc_simp, lcc_mult
from Bio.SeqUtils.CheckSum import crc32, crc64, gcg, seguid
from Bio.SeqUtils.CodonUsage import CodonAdaptationIndex


def u_crc32(seq):
    # NOTE - On Python 2 crc32 could return a signed int, but on Python 3 it is
    # always unsigned
    # Docs suggest should use crc32(x) & 0xffffffff for consistency.
    return crc32(seq) & 0xffffffff


class SeqUtilsTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Example of crc64 collision from Sebastian Bassi using the
        # immunoglobulin lambda light chain variable region from Homo sapiens
        # Both sequences share the same CRC64 checksum: 44CAAD88706CC153
        cls.str_light_chain_one = ("QSALTQPASVSGSPGQSITISCTGTSSDVGSYNLVSWYQQ"
                                   "HPGKAPKLMIYEGSKRPSGVSNRFSGSKSGNTASLTISGL"
                                   "QAEDEADYYCSSYAGSSTLVFGGGTKLTVL")
        cls.str_light_chain_two = ("QSALTQPASVSGSPGQSITISCTGTSSDVGSYNLVSWYQQ"
                                   "HPGKAPKLMIYEGSKRPSGVSNRFSGSKSGNTASLTISGL"
                                   "QAEDEADYYCCSYAGSSTWVFGGGTKLTVL")
        X = CodonAdaptationIndex()
        path = os.path.join("CodonUsage", "HighlyExpressedGenes.txt")
        X.generate_index(path)
        cls.X = X

    def test_codon_usage_ecoli(self):
        """Test Codon Adaptation Index (CAI) using default E. coli data."""
        CAI = CodonAdaptationIndex()
        value = CAI.cai_for_gene("ATGCGTATCGATCGCGATACGATTAGGCGGATG")
        self.assertAlmostEqual(value, 0.09978, places=5)

    def test_codon_usage_custom(self):
        """Test Codon Adaptation Index (CAI) using FASTA file for background."""
        # We need a FASTA file of CDS sequences to count the codon usage...
        dna_fasta_filename = "fasta.tmp"
        dna_genbank_filename = "GenBank/NC_005816.gb"
        record = SeqIO.read(dna_genbank_filename, "genbank")
        records = []
        for feature in record.features:
            if feature.type == "CDS" and len(feature.location.parts) == 1:
                start = feature.location.start.position
                end = feature.location.end.position
                table = int(feature.qualifiers["transl_table"][0])
                if feature.strand == -1:
                    seq = record.seq[start:end].reverse_complement()
                else:
                    seq = record.seq[start:end]
                # Double check we have the CDS sequence expected
                # TODO - Use any cds_start option if/when added to deal with the met
                a = "M" + str(seq[3:].translate(table))
                b = feature.qualifiers["translation"][0] + "*"
                self.assertEqual(a, b, "%r vs %r" % (a, b))
                records.append(SeqRecord(seq, id=feature.qualifiers["protein_id"][0],
                                         description=feature.qualifiers["product"][0]))

        with open(dna_fasta_filename, "w") as handle:
            SeqIO.write(records, handle, "fasta")

        CAI = CodonAdaptationIndex()
        # Note - this needs a FASTA file which containing non-ambiguous DNA coding
        # sequences - which should each be a whole number of codons.
        CAI.generate_index(dna_fasta_filename)
        # Now check codon usage index (CAI) using this species
        self.assertEqual(record.annotations["source"],
                         "Yersinia pestis biovar Microtus str. 91001")
        value = CAI.cai_for_gene("ATGCGTATCGATCGCGATACGATTAGGCGGATG")
        self.assertAlmostEqual(value, 0.67213, places=5)
        os.remove(dna_fasta_filename)

    def test_crc_checksum_collision(self):
        # Explicit testing of crc64 collision:
        self.assertNotEqual(self.str_light_chain_one, self.str_light_chain_two)
        self.assertNotEqual(crc32(self.str_light_chain_one), crc32(self.str_light_chain_two))
        self.assertEqual(crc64(self.str_light_chain_one), crc64(self.str_light_chain_two))
        self.assertNotEqual(gcg(self.str_light_chain_one), gcg(self.str_light_chain_two))
        self.assertNotEqual(seguid(self.str_light_chain_one), seguid(self.str_light_chain_two))

    def seq_checksums(self, seq_str, exp_crc32, exp_crc64, exp_gcg, exp_seguid,
                      exp_simple_LCC, exp_window_LCC):
        for s in [seq_str,
                  Seq(seq_str, single_letter_alphabet),
                  MutableSeq(seq_str, single_letter_alphabet)]:
            self.assertEqual(exp_crc32, u_crc32(s))
            self.assertEqual(exp_crc64, crc64(s))
            self.assertEqual(exp_gcg, gcg(s))
            self.assertEqual(exp_seguid, seguid(s))
            self.assertAlmostEqual(exp_simple_LCC, lcc_simp(s), places=2)
            values = lcc_mult(s, 20)
            self.assertEqual(len(exp_window_LCC), len(values))
            for value1, value2 in zip(exp_window_LCC, values):
                self.assertAlmostEqual(value1, value2, places=2)

    def test_checksum1(self):
        self.seq_checksums(self.str_light_chain_one,
                           2994980265,
                           "CRC-44CAAD88706CC153",
                           9729,
                           "BpBeDdcNUYNsdk46JoJdw7Pd3BI",
                           1.03,
                           (0.00, 1.00, 0.96, 0.96, 0.96, 0.65, 0.43, 0.35, 0.35, 0.35, 0.35, 0.53, 0.59, 0.26))

    def test_checksum2(self):
        self.seq_checksums(self.str_light_chain_two,
                           802105214,
                           "CRC-44CAAD88706CC153",
                           9647,
                           "X5XEaayob1nZLOc7eVT9qyczarY",
                           1.07,
                           (0.00, 1.00, 0.96, 0.96, 0.96, 0.65, 0.43, 0.35, 0.35, 0.35, 0.35, 0.53, 0.59, 0.26))

    def test_checksum3(self):
        self.seq_checksums("ATGCGTATCGATCGCGATACGATTAGGCGGAT",
                           817679856,
                           "CRC-6234FF451DC6DFC6",
                           7959,
                           "8WCUbVjBgiRmM10gfR7XJNjbwnE",
                           1.98,
                           (0.00, 2.00, 1.99, 1.99, 2.00, 1.99, 1.97, 1.99, 1.99, 1.99, 1.96, 1.96, 1.96, 1.96))

    def test_GC(self):
        seq = "ACGGGCTACCGTATAGGCAAGAGATGATGCCC"
        self.assertEqual(GC(seq), 56.25)

    def test_GC_skew(self):
        seq = "A" * 50
        self.assertEqual(GC_skew(seq)[0], 0)

    def test_seq1_seq3(self):
        s3 = "MetAlaTyrtrpcysthrLYSLEUILEGlYPrOGlNaSnaLapRoTyRLySSeRHisTrpLysThr"
        s1 = "MAYWCTKLIGPQNAPYKSHWKT"
        self.assertEqual(seq1(s3), s1)
        self.assertEqual(seq3(s1).upper(), s3.upper())
        self.assertEqual(seq1(seq3(s1)), s1)
        self.assertEqual(seq3(seq1(s3)).upper(), s3.upper())

    def test_codon_adaptation_index(self):
        X = self.X
        cai = X.cai_for_gene("ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA")
        self.assertAlmostEqual(cai, 0.6723, places=3)

    def test_index(self):
        X = self.X
        self.assertEqual(len(X.index), 64)
        self.assertAlmostEqual(X.index["AAA"], 1.000, places=3)
        self.assertAlmostEqual(X.index["AAC"], 1.000, places=3)
        self.assertAlmostEqual(X.index["AAG"], 0.219, places=3)
        self.assertAlmostEqual(X.index["AAT"], 0.293, places=3)
        self.assertAlmostEqual(X.index["ACA"], 0.110, places=3)
        self.assertAlmostEqual(X.index["ACC"], 1.000, places=3)
        self.assertAlmostEqual(X.index["ACG"], 0.204, places=3)
        self.assertAlmostEqual(X.index["ACT"], 0.517, places=3)
        self.assertAlmostEqual(X.index["AGA"], 0.018, places=3)
        self.assertAlmostEqual(X.index["AGC"], 0.762, places=3)
        self.assertAlmostEqual(X.index["AGG"], 0.006, places=3)
        self.assertAlmostEqual(X.index["AGT"], 0.195, places=3)
        self.assertAlmostEqual(X.index["ATA"], 0.015, places=3)
        self.assertAlmostEqual(X.index["ATC"], 1.000, places=3)
        self.assertAlmostEqual(X.index["ATG"], 1.000, places=3)
        self.assertAlmostEqual(X.index["ATT"], 0.490, places=3)
        self.assertAlmostEqual(X.index["CAA"], 0.259, places=3)
        self.assertAlmostEqual(X.index["CAC"], 1.000, places=3)
        self.assertAlmostEqual(X.index["CAG"], 1.000, places=3)
        self.assertAlmostEqual(X.index["CAT"], 0.416, places=3)
        self.assertAlmostEqual(X.index["CCA"], 0.247, places=3)
        self.assertAlmostEqual(X.index["CCC"], 0.040, places=3)
        self.assertAlmostEqual(X.index["CCG"], 1.000, places=3)
        self.assertAlmostEqual(X.index["CCT"], 0.161, places=3)
        self.assertAlmostEqual(X.index["CGA"], 0.023, places=3)
        self.assertAlmostEqual(X.index["CGC"], 0.531, places=3)
        self.assertAlmostEqual(X.index["CGG"], 0.014, places=3)
        self.assertAlmostEqual(X.index["CGT"], 1.000, places=3)
        self.assertAlmostEqual(X.index["CTA"], 0.017, places=3)
        self.assertAlmostEqual(X.index["CTC"], 0.100, places=3)
        self.assertAlmostEqual(X.index["CTG"], 1.000, places=3)
        self.assertAlmostEqual(X.index["CTT"], 0.085, places=3)
        self.assertAlmostEqual(X.index["GAA"], 1.000, places=3)
        self.assertAlmostEqual(X.index["GAC"], 1.000, places=3)
        self.assertAlmostEqual(X.index["GAG"], 0.308, places=3)
        self.assertAlmostEqual(X.index["GAT"], 0.886, places=3)
        self.assertAlmostEqual(X.index["GCA"], 0.794, places=3)
        self.assertAlmostEqual(X.index["GCC"], 0.538, places=3)
        self.assertAlmostEqual(X.index["GCG"], 0.937, places=3)
        self.assertAlmostEqual(X.index["GCT"], 1.000, places=3)
        self.assertAlmostEqual(X.index["GGA"], 0.056, places=3)
        self.assertAlmostEqual(X.index["GGC"], 0.892, places=3)
        self.assertAlmostEqual(X.index["GGG"], 0.103, places=3)
        self.assertAlmostEqual(X.index["GGT"], 1.000, places=3)
        self.assertAlmostEqual(X.index["GTA"], 0.465, places=3)
        self.assertAlmostEqual(X.index["GTC"], 0.297, places=3)
        self.assertAlmostEqual(X.index["GTG"], 0.618, places=3)
        self.assertAlmostEqual(X.index["GTT"], 1.000, places=3)
        self.assertAlmostEqual(X.index["TAA"], 1.000, places=3)
        self.assertAlmostEqual(X.index["TAC"], 1.000, places=3)
        self.assertAlmostEqual(X.index["TAG"], 0.012, places=3)
        self.assertAlmostEqual(X.index["TAT"], 0.606, places=3)
        self.assertAlmostEqual(X.index["TCA"], 0.221, places=3)
        self.assertAlmostEqual(X.index["TCC"], 0.785, places=3)
        self.assertAlmostEqual(X.index["TCG"], 0.240, places=3)
        self.assertAlmostEqual(X.index["TCT"], 1.000, places=3)
        self.assertAlmostEqual(X.index["TGA"], 0.081, places=3)
        self.assertAlmostEqual(X.index["TGC"], 1.000, places=3)
        self.assertAlmostEqual(X.index["TGG"], 1.000, places=3)
        self.assertAlmostEqual(X.index["TGT"], 0.721, places=3)
        self.assertAlmostEqual(X.index["TTA"], 0.059, places=3)
        self.assertAlmostEqual(X.index["TTC"], 1.000, places=3)
        self.assertAlmostEqual(X.index["TTG"], 0.072, places=3)
        self.assertAlmostEqual(X.index["TTT"], 0.457, places=3)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
