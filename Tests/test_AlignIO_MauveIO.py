# Copyright 2006-2014 by Peter Cock.  All rights reserved.
# Revisions copyright 2011 Brandon Invergo. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.AlignIO.MauveIO module."""
import os
import unittest

from io import StringIO

from Bio import SeqIO
from Bio.AlignIO.MauveIO import MauveIterator
from Bio.AlignIO.MauveIO import MauveWriter


class TestMauveIO(unittest.TestCase):
    MAUVE_TEST_DATA_DIR = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), "Mauve"
    )
    SIMPLE_XMFA = os.path.join(MAUVE_TEST_DATA_DIR, "simple.xmfa")
    SIMPLE_FA = os.path.join(MAUVE_TEST_DATA_DIR, "simple.fa")

    def test_one(self):
        ids = []
        with open(self.SIMPLE_XMFA) as handle:
            for alignment in MauveIterator(handle):
                for record in alignment:
                    ids.append(record.id)

        self.assertEqual(
            ids,
            [
                "1/0-5670",
                "2/0-5670",
                "1/5670-9940",
                "2/7140-11410",
                "1/9940-14910",
                "2/5670-7140",
                "2/11410-12880",
            ],
        )

        expected = """ATTCGCACAT AAGAATGTAC CTTGCTGTAA TTTATACTCA
            GCAGGTGGTG CAGACATCAT AACAAAAGAA GACTCTTGTT GTACTAGATA TTGTGTAGCA
            TCACGACCAC ACACACATGG AATGGAAACA CCTGTCTTAA GATTATCATA AGATAGAGTA
            CCCATATACA TCACAGCTTC TACACCCGTT AAGGTAGTAG TTTTCTGACC ACAATGTTTA
            CACACCACAT TAAGAACTCG CTTTGCAGAT TCCAAATTAG CATGCTGTAG AAGATGGGTC
            ATAGTTTCTC TGACATCACC AAGCTCGCCA ACAGTTTTAT TACTGTAAGC GAGTATGAGT
            GCACAAAAGT TAGCAGCATC ACCAGCACGG GCTCTATAAT AAGCCTCTTG AAGTGCTGGT
            GCATTGAATT TGACTTCAAG CTGTTGAAGT GCTAATAAAA CACTAGACAA ATAACAATTG
            TTATCAGCCC ATTTAATTGA AGTTAAACCA CCAACTTGAG GAAATTTCCA TTTCTTTGTG
            TGGTTTAAAG CAGACATGTA CCTACCAAGA AAACTCTCAT CAAGAGTATG GTAGTACTCG
            AAAGCTTCAC TACGTAGTGT GTCATCACTA GGTAGTACAA AGAAAGTCTT ACCCTCATGA
            TTTACATGAG GTTTAATTTT TGTAACATCA GCACCATCCA AGTATGTTGG ACCAAACTGC
            TGTCCATATG TCATAGACAT ATCCACAAGC TGTGTGTGGA GATTAGTGTT GTCCACAGTT
            GTGAACACTT TTATAGTCTT AACCTCCCGC AGGGATAAGA GACTCTTTAG TTTGTCAAGT
            GAAAGAACCT CACCGTCAAG ATGAAACTCG ACGGGGCTCT CCAGAGTGTG GTACACAATT
            TTGTCACCAC GCTTAAGAAA TTCAACACCT AACTCTGTAC GCTGTCCTGA ATAGGACCAA
            TCTCTGTAAG AGCCAGCCAA AGAAACTGTT TCTACAAAGT GCTCCTCAGA TGTCTTTGAT
            GACGAAGTGA GGTATCCATT ATATGTAGTA ACAGCATCTG GTGATGATAC TGACACTACG
            GCAGGAGCTT TAAGAGAACG CATACAGCGC GCAGCCTCTT CAAGATTAAA ACCATGTGTC
            ACATAACCAA TTGGCATTGT GACAAGCGGC TCATTTAGAG AGTTCAGCTT CGTAATAATA
            GAAGCTACAG GCTCTTTACT AGTATAAAAG AAGAATCGGA CACCATAGTC AACGATGCCC
            TCTTGAATTT TAATTCCTTT ATACTTACGT TGGATGGTTG CCATTATGGC TCTAACATCC
            ATGCATATAG GCATTAATTT TCTTGTCTCT TCAGCATGAG CAAGCATTTC TCTCAAATTC
            CAGGATACAG TTCCTAGAAT CTCTTCCTTA GCATTAGGTG CTTCTGAAGG TAGTACATAA
            AATGCAGATT TGCATTTCTT AAGAGCAGTC TTAGCTTCCT CAAGTGTATA """
        self.assertEqual(
            str(record.seq).replace("-", ""),
            expected.replace(" ", "").replace("\n", ""),
        )

    def test_sequence_positions(self):
        with open(self.SIMPLE_FA) as handle:
            seqs = list(SeqIO.parse(handle, "fasta"))

        with open(self.SIMPLE_XMFA) as handle:
            aln_list = list(MauveIterator(handle))

        for aln in aln_list:
            for record in aln:
                if not record.seq.startswith("-"):
                    expected = record.seq[0:10]
                    # seqs 0, 1 are ids 1, 2
                    actual = seqs[int(record.name) - 1].seq
                    # Slice out portion mentioned in file
                    actual = actual[
                        record.annotations["start"] : record.annotations["end"]
                    ]

                    if record.annotations["strand"] < 0:
                        actual = actual.reverse_complement()
                    # Slice first 10 chars for comparison, don't want to
                    # get any '-'s by accident
                    actual = actual[0:10]
                    if len(actual) == 0:
                        # We can't test sequences which don't provide
                        # proper annotation start/end/strand information
                        continue
                    self.assertEqual(expected, actual)

    def test_write_read(self):
        with open(self.SIMPLE_XMFA) as handle:
            aln_list = list(MauveIterator(handle))

        handle = StringIO()
        MauveWriter(handle).write_file(aln_list)
        handle.seek(0)
        aln_list_out = list(MauveIterator(handle))

        for a1, a2 in zip(aln_list, aln_list_out):
            self.assertEqual(len(a1), len(a2))
            for r1, r2 in zip(a1, a2):
                self.assertEqual(r1.id, r2.id)
                self.assertEqual(r1.seq, r2.seq)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
