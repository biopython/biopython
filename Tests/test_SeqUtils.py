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


import warnings
from Bio import BiopythonDeprecationWarning

with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonDeprecationWarning)
    from Bio.SeqUtils.CodonUsage import CodonAdaptationIndex as OldCodonAdaptationIndex


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

    def test_codon_usage_ecoli(self):
        """Test Codon Adaptation Index (CAI) using default E. coli data."""
        CAI = OldCodonAdaptationIndex()
        value = CAI.cai_for_gene("ATGCGTATCGATCGCGATACGATTAGGCGGATG")
        self.assertAlmostEqual(value, 0.09978, places=5)
        self.assertEqual(
            str(CAI),
            """\
AAA	1.000
AAC	1.000
AAG	0.253
AAT	0.051
ACA	0.076
ACC	1.000
ACG	0.099
ACT	0.965
AGA	0.004
AGC	0.410
AGG	0.002
AGT	0.085
ATA	0.003
ATC	1.000
ATG	1.000
ATT	0.185
CAA	0.124
CAC	1.000
CAG	1.000
CAT	0.291
CCA	0.135
CCC	0.012
CCG	1.000
CCT	0.070
CGA	0.004
CGC	0.356
CGG	0.004
CGT	1.000
CTA	0.007
CTC	0.037
CTG	1.000
CTT	0.042
GAA	1.000
GAC	1.000
GAG	0.259
GAT	0.434
GCA	0.586
GCC	0.122
GCG	0.424
GCT	1.000
GGA	0.010
GGC	0.724
GGG	0.019
GGT	1.000
GTA	0.495
GTC	0.066
GTG	0.221
GTT	1.000
TAC	1.000
TAT	0.239
TCA	0.077
TCC	0.744
TCG	0.017
TCT	1.000
TGC	1.000
TGG	1.000
TGT	0.500
TTA	0.020
TTC	1.000
TTG	0.020
TTT	0.296
""",
        )

    def test_codon_usage_custom_old(self):
        """Test Codon Adaptation Index (CAI) using FASTA file for background."""
        # We need a FASTA file of CDS sequences to count the codon usage...
        dna_fasta_filename = "fasta.tmp"
        dna_genbank_filename = "GenBank/NC_005816.gb"
        record = SeqIO.read(dna_genbank_filename, "genbank")
        records = []
        for feature in record.features:
            if feature.type == "CDS" and len(feature.location.parts) == 1:
                start = feature.location.start
                end = feature.location.end
                table = int(feature.qualifiers["transl_table"][0])
                if feature.strand == -1:
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

        with open(dna_fasta_filename, "w") as handle:
            SeqIO.write(records, handle, "fasta")

        CAI = OldCodonAdaptationIndex()
        # Note - this needs a FASTA file which containing non-ambiguous DNA coding
        # sequences - which should each be a whole number of codons.
        CAI.generate_index(dna_fasta_filename)
        # Now check codon usage index (CAI) using this species
        self.assertEqual(
            record.annotations["source"], "Yersinia pestis biovar Microtus str. 91001"
        )
        value = CAI.cai_for_gene("ATGCGTATCGATCGCGATACGATTAGGCGGATG")
        self.assertAlmostEqual(value, 0.67213, places=5)
        self.assertEqual(
            str(CAI),
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
TAG	0.000
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
        os.remove(dna_fasta_filename)

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
                if feature.strand == -1:
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

    def test_codon_adaptation_index(self):
        X = OldCodonAdaptationIndex()
        path = os.path.join("CodonUsage", "HighlyExpressedGenes.txt")
        X.generate_index(path)
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
        cai = X.cai_for_gene(
            "ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA"
        )
        self.assertAlmostEqual(cai, 0.6723, places=3)

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
