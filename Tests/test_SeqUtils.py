# Copyright 2007-2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
import unittest

from Bio import SeqIO
from Bio.Alphabet import single_letter_alphabet
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC, quick_FASTA_reader, seq1, seq3
from Bio.SeqUtils.lcc import lcc_simp, lcc_mult
from Bio.SeqUtils.CheckSum import crc32, crc64, gcg, seguid
from Bio.SeqUtils.CodonUsage import CodonAdaptationIndex


def u_crc32(seq):
    #NOTE - On Python 2 crc32 could return a signed int, but on Python 3 it is
    #always unsigned
    #Docs suggest should use crc32(x) & 0xffffffff for consistency.
    return crc32(seq) & 0xffffffff


def simple_LCC(s):
    #Avoid cross platforms with printing floats by doing conversion explicitly
    return "%0.2f" % lcc_simp(s)


def windowed_LCC(s):
    return ", ".join("%0.2f" % v for v in lcc_mult(s, 20))


class SeqUtilsTests(unittest.TestCase):

    def setUp(self):
        # Example of crc64 collision from Sebastian Bassi using the
        # immunoglobulin lambda light chain variable region from Homo sapiens
        # Both sequences share the same CRC64 checksum: 44CAAD88706CC153
        self.str_light_chain_one = "QSALTQPASVSGSPGQSITISCTGTSSDVGSYNLVSWYQQHPGK" \
                        + "APKLMIYEGSKRPSGVSNRFSGSKSGNTASLTISGLQAEDEADY" \
                        + "YCSSYAGSSTLVFGGGTKLTVL"
        self.str_light_chain_two = "QSALTQPASVSGSPGQSITISCTGTSSDVGSYNLVSWYQQHPGK" \
                        + "APKLMIYEGSKRPSGVSNRFSGSKSGNTASLTISGLQAEDEADY" \
                        + "YCCSYAGSSTWVFGGGTKLTVL"

    def test_quick_fasta_reader(self):
        dna_fasta_filename = "Fasta/f002"

        tuple_records = quick_FASTA_reader(dna_fasta_filename)
        self.assertEqual(len(tuple_records), 3)
        seq_records = list(SeqIO.parse(dna_fasta_filename, "fasta"))
        self.assertEqual(len(seq_records), 3)
        for tuple_record, seq_record in zip(tuple_records, seq_records):
            self.assertEqual(tuple_record, (seq_record.description, str(seq_record.seq)))

    def test_codon_usage_ecoli(self):
        """Test Codon Adaptation Index (CAI) using default E. coli data."""
        CAI = CodonAdaptationIndex()
        self.assertEqual("%0.5f" % CAI.cai_for_gene("ATGCGTATCGATCGCGATACGATTAGGCGGATG"),
                         "0.09978")

    def test_codon_usage_custom(self):
        """Test Codon Adaptation Index (CAI) using FASTA file for background."""
        #We need a FASTA file of CDS sequences to count the codon usage...
        dna_fasta_filename = "fasta.tmp"
        dna_genbank_filename = "GenBank/NC_005816.gb"
        record = SeqIO.read(dna_genbank_filename, "genbank")
        records = []
        for feature in record.features:
            if feature.type == "CDS" and not feature.sub_features:
                start = feature.location.start.position
                end = feature.location.end.position
                table = int(feature.qualifiers["transl_table"][0])
                if feature.strand == -1:
                    seq = record.seq[start:end].reverse_complement()
                else:
                    seq = record.seq[start:end]
                #Double check we have the CDS sequence expected
                #TODO - Use any cds_start option if/when added to deal with the met
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
        self.assertEqual("%0.5f" % CAI.cai_for_gene("ATGCGTATCGATCGCGATACGATTAGGCGGATG"),
                         "0.67213")
        os.remove(dna_fasta_filename)

    def test_crc_checksum_collision(self):
        #Explicit testing of crc64 collision:
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
            self.assertEqual(exp_simple_LCC, simple_LCC(s))
            self.assertEqual(exp_window_LCC, windowed_LCC(s))

    def test_checksum1(self):
        self.seq_checksums(self.str_light_chain_one,
                           2994980265,
                           "CRC-44CAAD88706CC153",
                           9729,
                           "BpBeDdcNUYNsdk46JoJdw7Pd3BI",
                           "1.03",
                           "0.00, 1.00, 0.96, 0.96, 0.96, 0.65, 0.43, 0.35, 0.35, 0.35, 0.35, 0.53, 0.59, 0.26")

    def test_checksum2(self):
        self.seq_checksums(self.str_light_chain_two,
                           802105214,
                           "CRC-44CAAD88706CC153",
                           9647,
                           "X5XEaayob1nZLOc7eVT9qyczarY",
                           "1.07",
                           "0.00, 1.00, 0.96, 0.96, 0.96, 0.65, 0.43, 0.35, 0.35, 0.35, 0.35, 0.53, 0.59, 0.26")

    def test_checksum3(self):
        self.seq_checksums("ATGCGTATCGATCGCGATACGATTAGGCGGAT",
                           817679856,
                           "CRC-6234FF451DC6DFC6",
                           7959,
                           "8WCUbVjBgiRmM10gfR7XJNjbwnE",
                           "1.98",
                           "0.00, 2.00, 1.99, 1.99, 2.00, 1.99, 1.97, 1.99, 1.99, 1.99, 1.96, 1.96, 1.96, 1.96")

    def test_GC(self):
        seq = "ACGGGCTACCGTATAGGCAAGAGATGATGCCC"
        self.assertEqual(GC(seq), 56.25)

    def test_seq1_seq3(self):
        s3 = "MetAlaTyrtrpcysthrLYSLEUILEGlYPrOGlNaSnaLapRoTyRLySSeRHisTrpLysThr"
        s1 = "MAYWCTKLIGPQNAPYKSHWKT"
        self.assertEqual(seq1(s3), s1)
        self.assertEqual(seq3(s1).upper(), s3.upper())
        self.assertEqual(seq1(seq3(s1)), s1)
        self.assertEqual(seq3(seq1(s3)).upper(), s3.upper())


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
