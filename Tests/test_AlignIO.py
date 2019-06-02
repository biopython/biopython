# Copyright 2008 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for AlignIO module."""


import unittest


from Bio._py3k import StringIO
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo, MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

test_write_read_alignment_formats = sorted(AlignIO._FormatToWriter)
test_write_read_align_with_seq_count = test_write_read_alignment_formats \
                                     + ["fasta", "tab"]


def str_summary(text, max_len=40):
    """Cuts *text* if too long for summaries."""
    if len(text) <= max_len:
        return text
    else:
        return text[:max_len - 4] + "..." + text[-3:]


def simple_alignment_comparison(alignments, alignments2, format):
    assert len(alignments) == len(alignments2)
    for a1, a2 in zip(alignments, alignments2):
        assert a1.get_alignment_length() == a2.get_alignment_length()
        assert len(a1) == len(a2)
        for r1, r2 in zip(a1, a2):
            # Check the bare minimum (ID and sequence) as
            # many formats can't store more than that.

            # Check the sequence
            assert str(r1.seq) == str(r2.seq)

            # Beware of different quirks and limitations in the
            # valid character sets and the identifier lengths!
            if format in ["phylip", "phylip-sequential"]:
                assert r1.id.replace("[", "").replace("]", "")[:10] == r2.id, \
                       "'%s' vs '%s'" % (r1.id, r2.id)
            elif format == "phylip-relaxed":
                assert r1.id.replace(" ", "").replace(':', '|') == r2.id, \
                        "'%s' vs '%s'" % (r1.id, r2.id)
            elif format == "clustal":
                assert r1.id.replace(" ", "_")[:30] == r2.id, \
                       "'%s' vs '%s'" % (r1.id, r2.id)
            elif format in ["stockholm", "maf"]:
                assert r1.id.replace(" ", "_") == r2.id, \
                       "'%s' vs '%s'" % (r1.id, r2.id)
            elif format == "fasta":
                assert r1.id.split()[0] == r2.id
            else:
                assert r1.id == r2.id, \
                       "'%s' vs '%s'" % (r1.id, r2.id)


class TestAlignIO_exceptions(unittest.TestCase):

    t_formats = list(AlignIO._FormatToWriter) + list(SeqIO._FormatToWriter)

    def test_phylip_reject_duplicate(self):
        """Check that writing duplicated IDs after truncation fails for PHYLIP."""
        handle = StringIO()
        sequences = [SeqRecord(Seq('AAAA'), id='longsequencename1'),
                     SeqRecord(Seq('AAAA'), id='longsequencename2'),
                     SeqRecord(Seq('AAAA'), id='other_sequence'), ]
        alignment = MultipleSeqAlignment(sequences)
        with self.assertRaises(ValueError) as cm:
            AlignIO.write(alignment, handle, 'phylip')
        self.assertEqual("Repeated name 'longsequen' (originally 'longsequencename2'), possibly due to truncation", str(cm.exception))

    def test_parsing_empty_files(self):
        """Check that parsing an empty file returns an empty list."""
        for t_format in AlignIO._FormatToIterator:
            handle = StringIO()
            alignments = list(AlignIO.parse(handle, t_format))
            self.assertEqual(alignments, [])

    def test_writing_empty_files(self):
        """Check that writers can cope with no alignments."""
        for t_format in self.t_formats:
            handle = StringIO()
            number = AlignIO.write([], handle, t_format)
            self.assertEqual(number, 0)

    def test_writing_not_alignments(self):
        """Check that writers reject records that are not alignments."""
        path = "Clustalw/opuntia.aln"
        records = list(AlignIO.read(path, "clustal"))
        for t_format in self.t_formats:
            handle = StringIO()
            self.assertRaises(Exception, AlignIO.write, [records], handle, t_format)


class TestAlignIO_reading(unittest.TestCase):

    def check_simple_write_read(self, alignments, indent=" "):
        for format in test_write_read_align_with_seq_count:
            records_per_alignment = len(alignments[0])
            for a in alignments:
                if records_per_alignment != len(a):
                    records_per_alignment = None
            # Can we expect this format to work?
            if not records_per_alignment \
                    and format not in test_write_read_alignment_formats:
                continue

            # Going to write to a handle...
            handle = StringIO()

            if format == 'nexus':
                with self.assertRaises(ValueError) as cm:
                    c = AlignIO.write(alignments, handle=handle, format=format)
                self.assertEqual("We can only write one Alignment to a Nexus file.", str(cm.exception))
                continue
            c = AlignIO.write(alignments, handle=handle, format=format)
            self.assertEqual(c, len(alignments))

            # First, try with the seq_count
            if records_per_alignment:
                handle.flush()
                handle.seek(0)
                alignments2 = list(AlignIO.parse(
                    handle=handle,
                    format=format,
                    seq_count=records_per_alignment))
                simple_alignment_comparison(alignments, alignments2, format)

            if format in test_write_read_alignment_formats:
                # Don't need the seq_count
                handle.flush()
                handle.seek(0)
                alignments2 = list(AlignIO.parse(handle=handle, format=format))
                simple_alignment_comparison(alignments, alignments2, format)

            # Try writing just one Alignment (not a list)
            handle = StringIO()
            AlignIO.write(alignments[0:1], handle, format)
            self.assertEqual(handle.getvalue(), alignments[0].format(format))

    def check_parse1(self, path, fmt, length, m=None):
        # Try using the iterator with a for loop and a handle
        with open(path, "r") as handle:
            alignments = list(AlignIO.parse(handle, format=fmt))
            self.assertEqual(len(alignments), length)
        if m is not None:
            for alignment in alignments:
                self.assertEqual(len(alignment), m)
        return alignments

    def check_parse2(self, path, fmt, length):
        # Try using the iterator with a for loop and a filename not handle
        counter = 0
        for record in AlignIO.parse(path, format=fmt):
            counter += 1
        self.assertEqual(counter, length)

    def check_parse3(self, path, fmt, length):
        # Try using the iterator with the next() method
        counter = 0
        alignments = AlignIO.parse(path, format=fmt)
        while True:
            try:
                alignment = next(alignments)
            except StopIteration:
                break
            self.assertIsNotNone(alignment)
            counter += 1
        self.assertEqual(counter, length)

    def check_parse4(self, path, fmt, length):
        # Try a mixture of next() and list
        counter = 0
        alignments = AlignIO.parse(path, format=fmt)
        alignment = next(alignments)
        counter = 1
        counter += len(list(alignments))
        self.assertEqual(counter, length)

    def check_parse5(self, path, fmt, length):
        # Try a mixture of next() and for loop
        alignments = AlignIO.parse(path, format=fmt)
        alignment = next(alignments)
        counter = 1
        for alignment in alignments:
            counter += 1
        self.assertEqual(counter, length)

    def check_read(self, path, fmt, m, k):
        # Check Bio.AlignIO.read(...)
        with open(path) as handle:
            alignment = AlignIO.read(handle, format=fmt)
        self.assertIsInstance(alignment, MultipleSeqAlignment)
        self.assertEqual(len(alignment), m)
        self.assertEqual(alignment.get_alignment_length(), k)
        return alignment

    def check_read_fails(self, path, fmt):
        with open(path) as handle:
            self.assertRaises(ValueError, AlignIO.read, handle, format=fmt)

    def check_summary(self, alignment, text, column_annotations={}):
        lines = []
        for record in alignment:
            line = "  %s %s" % (str_summary(str(record.seq)), record.id)
            lines.append(line)
        self.assertEqual(lines, text)
        self.assertEqual(alignment.column_annotations, column_annotations)

    def check_summary2(self, alignment, text):
        vertical_threshold = 5
        lines = []
        alignment_len = alignment.get_alignment_length()
        rec_count = len(alignment)
        if rec_count < vertical_threshold:
            # Show each sequence row horizontally
            for record in alignment:
                line = "  %s %s" % (str_summary(str(record.seq)), record.id)
                lines.append(line)
            for key, value in alignment.column_annotations.items():
                if isinstance(value, str):
                    line = "  %s %s" % (str_summary(str(value)), key)
                    lines.append(line)
        else:
            # Show each sequence row vertically
            for i in range(min(5, alignment_len)):
                line = "  " + str_summary(alignment[:, i]) + " alignment column %i" % i
                lines.append(line)
            if alignment_len > 5:
                i = alignment_len - 1
                line = "  " + str_summary("|" * rec_count) + " ..."
                lines.append(line)
                line = "  " + str_summary(alignment[:, i]) + " alignment column %i" % i
                lines.append(line)
        self.assertEqual(lines, text)

    def test_reading_alignments_clustal1(self):
        path = 'Clustalw/cw02.aln'
        self.check_parse1(path, "clustal", 1, 2)
        self.check_parse2(path, "clustal", 1)
        self.check_parse3(path, "clustal", 1)
        self.check_parse4(path, "clustal", 1)
        self.check_parse5(path, "clustal", 1)
        alignment = self.check_read(path, "clustal", 2, 601)
        self.check_summary(alignment, [
"  MENSDSNDKGSDQSAAQRRSQMDRLDREEAFYQFVN...SVV gi|4959044|gb|AAD34209.1|AF069",
"  ---------MSPQTETKASVGFKAGVKEYKLTYYTP...--- gi|671626|emb|CAA85685.1|"],
                                      {'clustal_consensus': '          * *: ::    :.   :*  :  :. : . :*  ::   .:   **                  **:...   *.*** ..          .:*   * *: .* :*        : :* .*                   *::.  .    .:: :*..*  :* .*   .. .  :    .  :    *. .:: : .      .* .  :  *.:     ..::   * .  ::  :  .*.    :.    :. .  .  .* **.*..  :..  *.. .    . ::*                         :.: .*:    :     * ::   ***  . * :. .  .  :  *: .:: :::   ..   . : :   ::  *    *  : .. :.* . ::.  :: * :  :   * *   :..  * ..  * :**                             .  .:. ..   :*.  ..: :. .  .:* * :   : * .             ..*:.  .**   *.*... :  ::   :* .*  ::* : :.  :.    :   '})
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="clustal", seq_count=2))), 3)
        handle.close()

    def test_reading_alignments_clustal2(self):
        path = 'Clustalw/opuntia.aln'
        self.check_parse1(path, "clustal", 1, 7)
        self.check_parse2(path, "clustal", 1)
        self.check_parse3(path, "clustal", 1)
        self.check_parse4(path, "clustal", 1)
        self.check_parse5(path, "clustal", 1)
        alignment = self.check_read(path, "clustal", 7, 156)
        self.check_summary2(alignment, [
"  TTTTTTT alignment column 0",
"  AAAAAAA alignment column 1",
"  TTTTTTT alignment column 2",
"  AAAAAAA alignment column 3",
"  CCCCCCC alignment column 4",
"  ||||||| ...",
"  AAAAAAA alignment column 155"])

        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="clustal", seq_count=7))), 3)
        handle.close()

    def test_reading_alignments_clustal3(self):
        path = 'Clustalw/hedgehog.aln'
        self.check_parse1(path, "clustal", 1, 5)
        self.check_parse2(path, "clustal", 1)
        self.check_parse3(path, "clustal", 1)
        self.check_parse4(path, "clustal", 1)
        self.check_parse5(path, "clustal", 1)
        alignment = self.check_read(path, "clustal", 5, 447)
        self.check_summary2(alignment, [
"  M---- alignment column 0",
"  F---- alignment column 1",
"  N---- alignment column 2",
"  L---- alignment column 3",
"  V---- alignment column 4",
"  ||||| ...",
"  ---SS alignment column 446"])

        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="clustal", seq_count=5))), 3)
        handle.close()

    def test_reading_alignments_clustal4(self):
        path = 'Clustalw/odd_consensus.aln'
        self.check_parse1(path, "clustal", 1, 2)
        self.check_parse2(path, "clustal", 1)
        self.check_parse3(path, "clustal", 1)
        self.check_parse4(path, "clustal", 1)
        self.check_parse5(path, "clustal", 1)
        alignment = self.check_read(path, "clustal", 2, 687)
        self.check_summary(alignment, [
"  ------------------------------------...TAG AT3G20900.1-CDS",
"  ATGAACAAAGTAGCGAGGAAGAACAAAACATCAGGT...TAG AT3G20900.1-SEQ"],
                                      {'clustal_consensus': '                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            *       *  *** ***** *   *  **      *******************************************************************************************************************************************************************************'})
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="clustal", seq_count=2))), 3)
        handle.close()

    def test_reading_alignments_clustal5(self):
        path = 'Clustalw/protein.aln'
        self.check_parse1(path, "clustal", 1, 20)
        self.check_parse2(path, "clustal", 1)
        self.check_parse3(path, "clustal", 1)
        self.check_parse4(path, "clustal", 1)
        self.check_parse5(path, "clustal", 1)
        alignment = self.check_read(path, "clustal", 20, 411)
        self.check_summary2(alignment, [
"  -M------------------ alignment column 0",
"  -T------------------ alignment column 1",
"  -V------------------ alignment column 2",
"  -L-----------------M alignment column 3",
"  -E---------------MMS alignment column 4",
"  |||||||||||||||||||| ...",
"  -------------------T alignment column 410"])

        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="clustal", seq_count=20))), 3)
        handle.close()

    def test_reading_alignments_clustal6(self):
        path = 'Clustalw/promals3d.aln'
        self.check_parse1(path, "clustal", 1, 20)
        self.check_parse2(path, "clustal", 1)
        self.check_parse3(path, "clustal", 1)
        self.check_parse4(path, "clustal", 1)
        self.check_parse5(path, "clustal", 1)
        alignment = self.check_read(path, "clustal", 20, 414)
        self.check_summary2(alignment, [
"  MMMMMMMMMMMMMMMM-M-- alignment column 0",
"  -----------------T-- alignment column 1",
"  -----------------V-- alignment column 2",
"  -----------------L-- alignment column 3",
"  -S---------------E-- alignment column 4",
"  |||||||||||||||||||| ...",
"  -T------------------ alignment column 413"])

        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="clustal", seq_count=20))), 3)
        handle.close()

    def test_reading_alignments_fasta(self):
        path = 'GFF/multi.fna'  # Trivial nucleotide alignment
        self.check_parse1(path, "fasta", 1, 3)
        self.check_parse2(path, "fasta", 1)
        self.check_parse3(path, "fasta", 1)
        self.check_parse4(path, "fasta", 1)
        self.check_parse5(path, "fasta", 1)
        alignment = self.check_read(path, "fasta", 3, 8)
        self.check_summary(alignment, [
"  ACGTCGCG test1",
"  GGGGCCCC test2",
"  AAACACAC test3"])

        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))

        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="fasta", seq_count=3))), 3)
        handle.close()

    def test_reading_alignments_nexus1(self):
        path = 'Nexus/test_Nexus_input.nex'
        self.check_parse1(path, "nexus", 1, 9)
        self.check_parse2(path, "nexus", 1)
        self.check_parse3(path, "nexus", 1)
        self.check_parse4(path, "nexus", 1)
        self.check_parse5(path, "nexus", 1)
        alignment = self.check_read(path, "nexus", 9, 48)
        self.check_summary2(alignment, [
"  AAAAAAAAc alignment column 0",
"  -----c?tc alignment column 1",
"  CCCCCCCCc alignment column 2",
"  --c-?a-tc alignment column 3",
"  GGGGGGGGc alignment column 4",
"  ||||||||| ...",
"  tt--?ag?c alignment column 47"])

        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()

    def test_reading_alignments_nexus2(self):
        path = 'Nexus/codonposset.nex'
        self.check_parse1(path, "nexus", 1, 2)
        self.check_parse2(path, "nexus", 1)
        self.check_parse3(path, "nexus", 1)
        self.check_parse4(path, "nexus", 1)
        self.check_parse5(path, "nexus", 1)
        alignment = self.check_read(path, "nexus", 2, 22)
        self.check_summary(alignment, [
"  AAAAAGGCATTGTGGTGGGAAT Aegotheles",
"  ?????????TTGTGGTGGGAAT Aerodramus"])

        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()

    def test_reading_alignments_stockholm1(self):
        path = 'Stockholm/simple.sth'
        self.check_parse1(path, "stockholm", 1, 2)
        self.check_parse2(path, "stockholm", 1)
        self.check_parse3(path, "stockholm", 1)
        self.check_parse4(path, "stockholm", 1)
        self.check_parse5(path, "stockholm", 1)
        alignment = self.check_read(path, "stockholm", 2, 104)
        self.check_summary(alignment, [
"  UUAAUCGAGCUCAACACUCUUCGUAUAUCCUC-UCA...UGU AP001509.1",
"  AAAAUUGAAUAUCGUUUUACUUGUUUAU-GUCGUGA...GAU AE007476.1"],
                                      {'secondary_structure': '.................<<<<<<<<...<<<<<<<........>>>>>>>........<<<<<<<.......>>>>>>>..>>>>>>>>...............'})
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="stockholm", seq_count=2))), 3)
        handle.close()

    def test_reading_alignments_stockholm2(self):
        path = 'Stockholm/funny.sth'
        self.check_parse1(path, "stockholm", 1, 6)
        self.check_parse2(path, "stockholm", 1)
        self.check_parse3(path, "stockholm", 1)
        self.check_parse4(path, "stockholm", 1)
        self.check_parse5(path, "stockholm", 1)
        alignment = self.check_read(path, "stockholm", 6, 43)
        self.check_summary2(alignment, [
"  MMMEEE alignment column 0",
"  TQIVVV alignment column 1",
"  CHEMMM alignment column 2",
"  RVALLL alignment column 3",
"  ASDTTT alignment column 4",
"  |||||| ...",
"  SYSEEE alignment column 42"])

        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="stockholm", seq_count=6))), 3)
        handle.close()

    def test_reading_alignments_phylip1(self):
        path = 'Phylip/reference_dna.phy'
        self.check_parse1(path, "phylip", 1, 6)
        self.check_parse2(path, "phylip", 1)
        self.check_parse3(path, "phylip", 1)
        self.check_parse4(path, "phylip", 1)
        self.check_parse5(path, "phylip", 1)
        alignment = self.check_read(path, "phylip", 6, 13)
        self.check_summary2(alignment, [
"  CCTTCG alignment column 0",
"  GGAAAG alignment column 1",
"  ATAAAC alignment column 2",
"  TTTTAA alignment column 3",
"  GAGGAG alignment column 4",
"  |||||| ...",
"  CTTTTC alignment column 12"])

        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="phylip", seq_count=6))), 3)
        handle.close()

    def test_reading_alignments_phylip2(self):
        path = 'Phylip/reference_dna2.phy'
        self.check_parse1(path, "phylip", 1, 6)
        self.check_parse2(path, "phylip", 1)
        self.check_parse3(path, "phylip", 1)
        self.check_parse4(path, "phylip", 1)
        self.check_parse5(path, "phylip", 1)
        alignment = self.check_read(path, "phylip", 6, 39)
        self.check_summary2(alignment, [
"  CCTTCG alignment column 0",
"  GGAAAG alignment column 1",
"  ATAAAC alignment column 2",
"  TTTTAA alignment column 3",
"  GAGGAG alignment column 4",
"  |||||| ...",
"  CTTTTC alignment column 38"])
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="phylip", seq_count=6))), 3)
        handle.close()

    def test_reading_alignments_phylip3(self):
        path = 'Phylip/hennigian.phy'
        self.check_parse1(path, "phylip", 1, 10)
        self.check_parse2(path, "phylip", 1)
        self.check_parse3(path, "phylip", 1)
        self.check_parse4(path, "phylip", 1)
        self.check_parse5(path, "phylip", 1)
        alignment = self.check_read(path, "phylip", 10, 40)
        self.check_summary2(alignment, [
"  CCCCCAAAAA alignment column 0",
"  AAAAACCCCC alignment column 1",
"  CCCAAAAAAA alignment column 2",
"  AAACCAAAAA alignment column 3",
"  CCAAAAAAAA alignment column 4",
"  |||||||||| ...",
"  AAAAAAAAAA alignment column 39"])
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="phylip", seq_count=10))), 3)
        handle.close()

    def test_reading_alignments_phylip4(self):
        path = 'Phylip/horses.phy'
        self.check_parse1(path, "phylip", 1, 10)
        self.check_parse2(path, "phylip", 1)
        self.check_parse3(path, "phylip", 1)
        self.check_parse4(path, "phylip", 1)
        self.check_parse5(path, "phylip", 1)
        alignment = self.check_read(path, "phylip", 10, 40)
        self.check_summary2(alignment, [
"  AACCCCCCCC alignment column 0",
"  AAAACCCCCC alignment column 1",
"  AAAAAAAAAC alignment column 2",
"  ACAAAAAAAA alignment column 3",
"  ACACCCCCCC alignment column 4",
"  |||||||||| ...",
"  AAAAAAAAAA alignment column 39"])
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="phylip", seq_count=10))), 3)
        handle.close()

    def test_reading_alignments_phylip5(self):
        path = 'Phylip/random.phy'
        self.check_parse1(path, "phylip", 1, 10)
        self.check_parse2(path, "phylip", 1)
        self.check_parse3(path, "phylip", 1)
        self.check_parse4(path, "phylip", 1)
        self.check_parse5(path, "phylip", 1)
        alignment = self.check_read(path, "phylip", 10, 40)
        self.check_summary2(alignment, [
"  CAAAACAAAC alignment column 0",
"  AACAACCACC alignment column 1",
"  CAAAACAAAA alignment column 2",
"  ACAACACACA alignment column 3",
"  CCAAAACCAA alignment column 4",
"  |||||||||| ...",
"  AAAAAAAAAA alignment column 39"])
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="phylip", seq_count=10))), 3)
        handle.close()

    def test_reading_alignments_phylip6(self):
        path = 'Phylip/interlaced.phy'
        self.check_parse1(path, "phylip", 1, 3)
        self.check_parse2(path, "phylip", 1)
        self.check_parse3(path, "phylip", 1)
        self.check_parse4(path, "phylip", 1)
        self.check_parse5(path, "phylip", 1)
        alignment = self.check_read(path, "phylip", 3, 384)
        self.check_summary(alignment, [
"  -----MKVILLFVLAVFTVFVSS-------------...I-- CYS1_DICDI",
"  MAHARVLLLALAVLATAAVAVASSSSFADSNPIRPV...VAA ALEU_HORVU",
"  ------MWATLPLLCAGAWLLGV--------PVCGA...PLV CATH_HUMAN"])
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="phylip", seq_count=3))), 3)
        handle.close()

    def test_reading_alignments_phylip7(self):
        path = 'Phylip/interlaced2.phy'
        self.check_parse1(path, "phylip", 1, 4)
        self.check_parse2(path, "phylip", 1)
        self.check_parse3(path, "phylip", 1)
        self.check_parse4(path, "phylip", 1)
        self.check_parse5(path, "phylip", 1)
        alignment = self.check_read(path, "phylip", 4, 131)
        self.check_summary(alignment, [
"  TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTG...SHE IXI_234",
"  TSPASIRPPAGPSSR---------RPSPPGPRRPTG...SHE IXI_235",
"  TSPASIRPPAGPSSRPAMVSSR--RPSPPPPRRPPG...SHE IXI_236",
"  TSPASLRPPAGPSSRPAMVSSRR-RPSPPGPRRPT-...SHE IXI_237"])
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="phylip", seq_count=4))), 3)
        handle.close()

    def test_reading_alignments_phylip8(self):
        path = 'ExtendedPhylip/primates.phyx'
        self.check_parse1(path, "phylip-relaxed", 1, 12)
        self.check_parse2(path, "phylip-relaxed", 1)
        self.check_parse3(path, "phylip-relaxed", 1)
        self.check_parse4(path, "phylip-relaxed", 1)
        self.check_parse5(path, "phylip-relaxed", 1)
        alignment = self.check_read(path, "phylip-relaxed", 12, 898)
        self.check_summary2(alignment, [
"  AAAAAAAAAAAA alignment column 0",
"  AAAAAAAAAAAA alignment column 1",
"  GGGGGGGGGGGG alignment column 2",
"  TCCCCCCCCCCC alignment column 3",
"  TTTTTTTTTTTT alignment column 4",
"  |||||||||||| ...",
"  TTTTTTTTTTTT alignment column 897"])

        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="phylip-relaxed", seq_count=12))), 3)
        handle.close()

    def test_reading_alignments_phylip9(self):
        path = 'Phylip/sequential.phy'
        self.check_parse1(path, "phylip-sequential", 1, 3)
        self.check_parse2(path, "phylip-sequential", 1)
        self.check_parse3(path, "phylip-sequential", 1)
        self.check_parse4(path, "phylip-sequential", 1)
        self.check_parse5(path, "phylip-sequential", 1)
        alignment = self.check_read(path, "phylip-sequential", 3, 384)
        self.check_summary(alignment, [
"  -----MKVILLFVLAVFTVFVSS-------------...I-- CYS1_DICDI",
"  MAHARVLLLALAVLATAAVAVASSSSFADSNPIRPV...VAA ALEU_HORVU",
"  ------MWATLPLLCAGAWLLGV--------PVCGA...PLV CATH_HUMAN"])

        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="phylip-sequential", seq_count=3))), 3)
        handle.close()

    def test_reading_alignments_phylip10(self):
        path = 'Phylip/sequential2.phy'
        self.check_parse1(path, "phylip-sequential", 1, 4)
        self.check_parse2(path, "phylip-sequential", 1)
        self.check_parse3(path, "phylip-sequential", 1)
        self.check_parse4(path, "phylip-sequential", 1)
        self.check_parse5(path, "phylip-sequential", 1)
        alignment = self.check_read(path, "phylip-sequential", 4, 131)
        self.check_summary(alignment, [
"  TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTG...SHE IXI_234",
"  TSPASIRPPAGPSSR---------RPSPPGPRRPTG...SHE IXI_235",
"  TSPASIRPPAGPSSRPAMVSSR--RPSPPPPRRPPG...SHE IXI_236",
"  TSPASLRPPAGPSSRPAMVSSRR-RPSPPGPRRPT-...SHE IXI_237"])

        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="phylip-sequential", seq_count=4))), 3)
        handle.close()

    def test_reading_alignments_emboss1(self):
        path = 'Emboss/alignret.txt'
        self.check_parse1(path, "emboss", 1, 4)
        self.check_parse2(path, "emboss", 1)
        self.check_parse3(path, "emboss", 1)
        self.check_parse4(path, "emboss", 1)
        self.check_parse5(path, "emboss", 1)
        alignment = self.check_read(path, "emboss", 4, 131)
        self.check_summary(alignment, [
"  TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTG...SHE IXI_234",
"  TSPASIRPPAGPSSR---------RPSPPGPRRPTG...SHE IXI_235",
"  TSPASIRPPAGPSSRPAMVSSR--RPSPPPPRRPPG...SHE IXI_236",
"  TSPASLRPPAGPSSRPAMVSSRR-RPSPPGPRRPT-...SHE IXI_237"])
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))

    def test_reading_alignments_emboss2(self):
        path = 'Emboss/needle.txt'
        alignments = self.check_parse1(path, "emboss", 5, 2)
        self.check_parse2(path, "emboss", 5)
        self.check_parse3(path, "emboss", 5)
        self.check_parse4(path, "emboss", 5)
        self.check_parse5(path, "emboss", 5)
        self.check_read_fails(path, "emboss")
        # Show the alignment
        self.assertEqual(alignments[0].get_alignment_length(), 124)
        self.assertEqual(alignments[1].get_alignment_length(), 119)
        self.assertEqual(alignments[2].get_alignment_length(), 120)
        self.assertEqual(alignments[3].get_alignment_length(), 118)
        self.assertEqual(alignments[4].get_alignment_length(), 125)
        self.check_summary(alignments[0], [
"  KILIVDD----QYGIRILLNEVFNKEGYQTFQAANG...--- ref_rec",
"  -VLLADDHALVRRGFRLMLED--DPEIEIVAEAGDG...GET gi|94968718|receiver"])
        self.check_summary(alignments[1], [
"  KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQAL...--- ref_rec",
"  -ILIVDDEANTLASLSRAFRLAGHEATVCDNAVRAL...LKR gi|94968761|receiver"])
        self.check_summary(alignments[2], [
"  -KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQA...--- ref_rec",
"  LHIVVVDDDPGTCVYIESVFAELGHTCKSFVRPEAA...HKE gi|94967506|receiver"])
        self.check_summary(alignments[3], [
"  KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQAL...DAV ref_rec",
"  -VLLVEDEEALRAAAGDFLETRGYKIMTARDGTEAL...EVL gi|94970045|receiver"])
        self.check_summary(alignments[4], [
"  KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQAL...--- ref_rec",
"  TVLLVEDEEGVRKLVRGILSRQGYHVLEATSGEEAL...KRQ gi|94970041|receiver"])
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignments[0])
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        # Some alignment file formats have magic characters which mean
        # use the letter in this position in the first sequence.
        # They should all have been converted by the parser, but if
        # not reversing the record order might expose an error.  Maybe.
        alignments.reverse()
        self.check_simple_write_read(alignments)

    def test_reading_alignments_emboss3(self):
        path = 'Emboss/needle_asis.txt'
        self.check_parse1(path, "emboss", 1, 2)
        self.check_parse2(path, "emboss", 1)
        self.check_parse3(path, "emboss", 1)
        self.check_parse4(path, "emboss", 1)
        self.check_parse5(path, "emboss", 1)
        alignment = self.check_read(path, "emboss", 2, 3653)
        self.check_summary(alignment, [
"  TATTTTTTGGATTTTTTTCTAGATTTTCTAGGTTAT...GAA asis",
"  ------------------------------------...GAA asis"])
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))

    def test_reading_alignments_emboss4(self):
        path = 'Emboss/water.txt'
        self.check_parse1(path, "emboss", 1, 2)
        self.check_parse2(path, "emboss", 1)
        self.check_parse3(path, "emboss", 1)
        self.check_parse4(path, "emboss", 1)
        self.check_parse5(path, "emboss", 1)
        alignment = self.check_read(path, "emboss", 2, 131)
        self.check_summary(alignment, [
"  TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTG...SHE IXI_234",
"  TSPASIRPPAGPSSR---------RPSPPGPRRPTG...SHE IXI_235"])

        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))

    def test_reading_alignments_emboss5(self):
        path = 'Emboss/water2.txt'
        self.check_parse1(path, "emboss", 1, 2)
        self.check_parse2(path, "emboss", 1)
        self.check_parse3(path, "emboss", 1)
        self.check_parse4(path, "emboss", 1)
        self.check_parse5(path, "emboss", 1)
        alignment = self.check_read(path, "emboss", 2, 18)
        self.check_summary(alignment, [
"  CGTTTGAGT-CTGGGATG asis",
"  CGTTTGAGTACTGGGATG asis"])
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))

    def test_reading_alignments_emboss6(self):
        path = 'Emboss/matcher_simple.txt'
        self.check_parse1(path, "emboss", 1, 2)
        self.check_parse2(path, "emboss", 1)
        self.check_parse3(path, "emboss", 1)
        self.check_parse4(path, "emboss", 1)
        self.check_parse5(path, "emboss", 1)
        alignment = self.check_read(path, "emboss", 2, 16)
        self.check_summary(alignment, [
"  GPPPQSPDENRAGESS AF069992_1",
"  GVPPEEAGAAVAAESS CAA85685.1"])
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))

    def test_reading_alignments_emboss7(self):
        path = 'Emboss/matcher_pair.txt'
        alignments = self.check_parse1(path, "emboss", 5, 2)
        self.check_parse2(path, "emboss", 5)
        self.check_parse3(path, "emboss", 5)
        self.check_parse4(path, "emboss", 5)
        self.check_parse5(path, "emboss", 5)
        self.check_read_fails(path, "emboss")
        self.assertEquals(alignments[0].get_alignment_length(), 145)
        self.assertEquals(alignments[1].get_alignment_length(), 13)
        self.assertEquals(alignments[2].get_alignment_length(), 18)
        self.assertEquals(alignments[3].get_alignment_length(), 10)
        self.assertEquals(alignments[4].get_alignment_length(), 10)
        self.check_summary(alignments[0], [
"  LSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFP...SKY HBA_HUMAN",
"  LTPEEKSAVTALWGKV--NVDEVGGEALGRLLVVYP...HKY HBB_HUMAN"])
        self.check_summary(alignments[1], [
"  KKVADALTNAVAH HBA_HUMAN",
"  QKVVAGVANALAH HBB_HUMAN"])
        self.check_summary(alignments[2], [
"  KLRVDPVNFKLLSHCLLV HBA_HUMAN",
"  KVNVDEVGGEALGRLLVV HBB_HUMAN"])
        self.check_summary(alignments[3], [
"  LSALSDLHAH HBA_HUMAN",
"  LGAFSDGLAH HBB_HUMAN"])
        self.check_summary(alignments[4], [
"  VKAAWGKVGA HBA_HUMAN",
"  VQAAYQKVVA HBB_HUMAN"])
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignments[0])
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        # Some alignment file formats have magic characters which mean
        # use the letter in this position in the first sequence.
        # They should all have been converted by the parser, but if
        # not reversing the record order might expose an error.  Maybe.
        alignments.reverse()
        self.check_simple_write_read(alignments)

    def test_reading_alignments_emboss8(self):
        path = 'Emboss/emboss_pair_aln_full_blank_line.txt'
        self.check_parse1(path, "emboss", 1, 2)
        self.check_parse2(path, "emboss", 1)
        self.check_parse3(path, "emboss", 1)
        self.check_parse4(path, "emboss", 1)
        self.check_parse5(path, "emboss", 1)
        alignment = self.check_read(path, "emboss", 2, 1450)
        self.check_summary(alignment, [
"  GGCAGGTGCATAGCTTGAGCCTAGGAGTTCAAGTCC...AAA hg38_chrX_131691529_131830643_47210_48660",
"  G--------------------------TTCAAGGCC...AAA mm10_chrX_50555743_50635321_27140_27743"])
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))

    def test_reading_alignments_fasta_m10_1(self):
        path = 'Fasta/output001.m10'
        alignments = self.check_parse1(path, "fasta-m10", 4, 2)
        self.check_parse2(path, "fasta-m10", 4)
        self.check_parse3(path, "fasta-m10", 4)
        self.check_parse4(path, "fasta-m10", 4)
        self.check_parse5(path, "fasta-m10", 4)
        self.check_read_fails(path, "fasta-m10")
        self.assertEqual(alignments[0].get_alignment_length(), 108)
        self.check_summary(alignments[0], [
"  SGSNT-RRRAISRPVRLTAEED---QEIRKRAAECG...LSR gi|10955263|ref|NP_052604.1|",
"  AGSGAPRRRGSGLASRISEQSEALLQEAAKHAAEFG...LSR gi|152973457|ref|YP_001338508.1|"])
        self.assertEqual(alignments[1].get_alignment_length(), 64)
        self.check_summary(alignments[1], [
"  AAECGKTVSGFLRAAALGKKVNSLTDDRVLKEV-MR...AIT gi|10955263|ref|NP_052604.1|",
"  ASRQGCTVGG--KMDSVQDKASDKDKERVMKNINIM...TLT gi|152973588|ref|YP_001338639.1|"])
        self.assertEqual(alignments[2].get_alignment_length(), 38)
        self.check_summary(alignments[2], [
"  MKKDKKYQIEAIKNKDKTLFIVYATDIYSPSEFFSKIE gi|10955264|ref|NP_052605.1|",
"  IKKDLGVSFLKLKNREKTLIVDALKKKYPVAELLSVLQ gi|152973462|ref|YP_001338513.1|"])
        self.assertEqual(alignments[3].get_alignment_length(), 43)
        self.check_summary(alignments[3], [
"  SELHSKLPKSIDKIHEDIKKQLSC-SLIMKKIDVEM...TYC gi|10955265|ref|NP_052606.1|",
"  SRINSDVARRIPGIHRDPKDRLSSLKQVEEALDMLI...EYC gi|152973545|ref|YP_001338596.1|"])
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignments[0])
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        # Some alignment file formats have magic characters which mean
        # use the letter in this position in the first sequence.
        # They should all have been converted by the parser, but if
        # not reversing the record order might expose an error.  Maybe.
        alignments.reverse()
        self.check_simple_write_read(alignments)

    def test_reading_alignments_fasta_m10_2(self):
        path = 'Fasta/output002.m10'
        alignments = self.check_parse1(path, "fasta-m10", 6, 2)
        self.check_parse2(path, "fasta-m10", 6)
        self.check_parse3(path, "fasta-m10", 6)
        self.check_parse4(path, "fasta-m10", 6)
        self.check_parse5(path, "fasta-m10", 6)
        self.check_read_fails(path, "fasta-m10")
        self.assertEqual(alignments[0].get_alignment_length(), 88)
        self.check_summary(alignments[0], [
"  SGSNTRRRAISRPVR--LTAEEDQEIRKRAAECG-K...AEV gi|10955263|ref|NP_052604.1|",
"  SQRSTRRKPENQPTRVILFNKPYDVLPQFTDEAGRK...VQV gi|162139799|ref|NP_309634.2|"])
        self.assertEqual(alignments[1].get_alignment_length(), 53)
        self.check_summary(alignments[1], [
"  EIRKRAAECGKTVSGFLRAAA-LGKKV----NSLTD...KKL gi|10955263|ref|NP_052604.1|",
"  EIKPRGTSKGEAIAAFMQEAPFIGRTPVFLGDDLTD...VKI gi|15831859|ref|NP_310632.1|"])
        self.assertEqual(alignments[2].get_alignment_length(), 92)
        self.check_summary(alignments[2], [
"  SEFFSKIESDLKKKKSKGDVFFDLIIPNG-----GK...ATS gi|10955264|ref|NP_052605.1|",
"  TELNSELAKAMKVDAQRG-AFVSQVLPNSSAAKAGI...QSS gi|15829419|ref|NP_308192.1|"])
        self.assertEqual(alignments[5].get_alignment_length(), 157)
        self.check_summary(alignments[5], [
"  QYIMTTSNGDRVRAKIYKRGSIQFQGKYLQIASLIN...REI gi|10955265|ref|NP_052606.1|",
"  EFIRLLSDHDQFEKDQISELTVAANALKLEVAK--N...KKV gi|15833861|ref|NP_312634.1|"])
        summary = AlignInfo.SummaryInfo(alignments[0])
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        # Some alignment file formats have magic characters which mean
        # use the letter in this position in the first sequence.
        # They should all have been converted by the parser, but if
        # not reversing the record order might expose an error.  Maybe.
        alignments.reverse()
        self.check_simple_write_read(alignments)

    def test_reading_alignments_fasta_m10_3(self):
        path = 'Fasta/output003.m10'
        alignments = self.check_parse1(path, "fasta-m10", 3, 2)
        self.check_parse2(path, "fasta-m10", 3)
        self.check_parse3(path, "fasta-m10", 3)
        self.check_parse4(path, "fasta-m10", 3)
        self.check_parse5(path, "fasta-m10", 3)
        self.check_read_fails(path, "fasta-m10")
        self.assertEqual(alignments[0].get_alignment_length(), 55)
        self.check_summary(alignments[0], [
"  ISISNNKDQYEELQKEQGERDLKTVDQLVRIAAAGG...IAA gi|152973837|ref|YP_001338874.1|",
"  VRLTAEEDQ--EIRKRAAECG-KTVSGFLRAAALGK...LGA gi|10955263|ref|NP_052604.1|"])
        self.assertEqual(alignments[1].get_alignment_length(), 22)
        self.check_summary(alignments[1], [
"  DDAEHLFRTLSSR-LDALQDGN gi|152973840|ref|YP_001338877.1|",
"  DDRANLFEFLSEEGITITEDNN gi|10955265|ref|NP_052606.1|"])
        self.assertEqual(alignments[2].get_alignment_length(), 63)
        self.check_summary(alignments[2], [
"  VFGSFEQPKGEHLSGQVSEQ--RDTAFADQNEQVIR...QAM gi|152973841|ref|YP_001338878.1|",
"  VYTSFN---GEKFSSYTLNKVTKTDEYNDLSELSAS...KGI gi|10955264|ref|NP_052605.1|"])
        summary = AlignInfo.SummaryInfo(alignments[0])
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        # Some alignment file formats have magic characters which mean
        # use the letter in this position in the first sequence.
        # They should all have been converted by the parser, but if
        # not reversing the record order might expose an error.  Maybe.
        alignments.reverse()
        self.check_simple_write_read(alignments)

    def test_reading_alignments_fasta_m10_4(self):
        path = 'Fasta/output004.m10'
        self.check_parse1(path, "fasta-m10", 1, 2)
        self.check_parse2(path, "fasta-m10", 1)
        self.check_parse3(path, "fasta-m10", 1)
        self.check_parse4(path, "fasta-m10", 1)
        self.check_parse5(path, "fasta-m10", 1)
        alignment = self.check_read(path, "fasta-m10", 2, 102)
        self.check_summary(alignment, [
"  AAAAAAGATAAAAAATATCAAATAGAAGCAATAAAA...TCA ref|NC_002127.1|:c1351-971",
"  AGAGAAAATAAAACAAGTAATAAAATATTAATGGAA...ACA ref|NC_002695.1|:1970775-1971404"])
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))

    def test_reading_alignments_fasta_m10_5(self):
        path = 'Fasta/output005.m10'
        self.check_parse1(path, "fasta-m10", 1, 2)
        self.check_parse2(path, "fasta-m10", 1)
        self.check_parse3(path, "fasta-m10", 1)
        self.check_parse4(path, "fasta-m10", 1)
        self.check_parse5(path, "fasta-m10", 1)
        alignment = self.check_read(path, "fasta-m10", 2, 110)
        self.check_summary(alignment, [
"  IKNKDKTLFIVYAT-DIYSPSEFFSKIESDLKKKKS...LSK gi|10955264|ref|NP_052605.1|",
"  IKDELPVAFCSWASLDLECEVKYINDVTSLYAKDWM...MSE gi|10955282|ref|NP_052623.1|"])
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))

    def test_reading_alignments_fasta_m10_6(self):
        path = 'Fasta/output006.m10'
        self.check_parse1(path, "fasta-m10", 1, 2)
        self.check_parse2(path, "fasta-m10", 1)
        self.check_parse3(path, "fasta-m10", 1)
        self.check_parse4(path, "fasta-m10", 1)
        self.check_parse5(path, "fasta-m10", 1)
        alignment = self.check_read(path, "fasta-m10", 2, 131)
        self.check_summary(alignment, [
"  GCAACGCTTCAAGAACTGGAATTAGGAACCGTGACA...CAT query",
"  GCAACGCTTCAAGAACTGGAATTAGGAACCGTGACA...CAT gi|116660610|gb|EG558221.1|EG558221"])
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))

    def test_reading_alignments_fasta_m10_7(self):
        path = 'Fasta/output007.m10'
        alignments = self.check_parse1(path, "fasta-m10", 9, 2)
        self.check_parse2(path, "fasta-m10", 9)
        self.check_parse3(path, "fasta-m10", 9)
        self.check_parse4(path, "fasta-m10", 9)
        self.check_parse5(path, "fasta-m10", 9)
        self.check_read_fails(path, "fasta-m10")
        self.assertEqual(alignments[0].get_alignment_length(), 108)
        self.check_summary(alignments[0], [
"  SGSNT-RRRAISRPVRLTAEED---QEIRKRAAECG...LSR gi|10955263|ref|NP_052604.1|",
"  AGSGAPRRRGSGLASRISEQSEALLQEAAKHAAEFG...LSR gi|152973457|ref|YP_001338508.1|"])
        self.assertEqual(alignments[1].get_alignment_length(), 64)
        self.check_summary(alignments[1], [
"  AAECGKTVSGFLRAAALGKKVNSLTDDRVLKEV-MR...AIT gi|10955263|ref|NP_052604.1|",
"  ASRQGCTVGG--KMDSVQDKASDKDKERVMKNINIM...TLT gi|152973588|ref|YP_001338639.1|"])
        self.assertEqual(alignments[2].get_alignment_length(), 45)
        self.check_summary(alignments[2], [
"  EIRKRAAECGKTVSGFLRAAA-----LGKKVNSLTD...VMR gi|10955263|ref|NP_052604.1|",
"  ELVKLIADMGISVRALLRKNVEPYEELGLEEDKFTD...MLQ gi|152973480|ref|YP_001338531.1|"])
        self.assertEqual(alignments[8].get_alignment_length(), 64)
        self.check_summary(alignments[8], [
"  ISGTYKGIDFLIKLMPSGGNTTIGRASGQNNTYFDE...FSD gi|10955265|ref|NP_052606.1|",
"  IDGVITAFD-LRTGMNISKDKVVAQIQGMDPVW---...YPD gi|152973505|ref|YP_001338556.1|"])
        summary = AlignInfo.SummaryInfo(alignments[0])
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        # Some alignment file formats have magic characters which mean
        # use the letter in this position in the first sequence.
        # They should all have been converted by the parser, but if
        # not reversing the record order might expose an error.  Maybe.
        alignments.reverse()
        self.check_simple_write_read(alignments)

    def test_reading_alignments_fasta_m10_8(self):
        path = 'Fasta/output008.m10'
        alignments = self.check_parse1(path, "fasta-m10", 12, 2)
        self.check_parse2(path, "fasta-m10", 12)
        self.check_parse3(path, "fasta-m10", 12)
        self.check_parse4(path, "fasta-m10", 12)
        self.check_parse5(path, "fasta-m10", 12)
        self.check_read_fails(path, "fasta-m10")
        self.assertEqual(alignments[0].get_alignment_length(), 65)
        self.check_summary(alignments[0], [
"  LQHRHPHQQQQQQQQQQQQQQQQQQQQQQQQQQQH-...QML sp|Q9NSY1|BMP2K_HUMAN",
"  IPHQLPHALRHRPAQEAAHASQLHPAQPGCGQPLHG...GLL gi|283855822|gb|GQ290312.1|"])
        self.assertEqual(alignments[1].get_alignment_length(), 201)
        self.check_summary(alignments[1], [
"  GPEIL---LGQ-GPPQQPPQQHRVLQQLQQGDWRLQ...NRS sp|Q9NSY1|BMP2K_HUMAN",
"  GPELLRALLQQNGCGTQPLRVPTVLPG*AMAVLHAG...QKS gi|57163782|ref|NM_001009242.1|"])
        self.assertEqual(alignments[2].get_alignment_length(), 348)
        self.check_summary(alignments[2], [
"  MNGTEGPNFYVPFSNATGVVRSPFEYPQYYLAEPWQ...APA sp|P08100|OPSD_HUMAN",
"  MNGTEGPNFYVPFSNKTGVVRSPFEYPQYYLAEPWQ...APA gi|57163782|ref|NM_001009242.1|"])
        self.assertEqual(alignments[11].get_alignment_length(), 31)
        self.check_summary(alignments[11], [
"  AQQQESATTQKAEKEVTRMVIIMVIAFLICW sp|P08100|OPSD_HUMAN",
"  SQQIRNATTMMMTMRVTSFSAFWVVADSCCW gi|283855822|gb|GQ290312.1|"])
        summary = AlignInfo.SummaryInfo(alignments[0])
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        # Some alignment file formats have magic characters which mean
        # use the letter in this position in the first sequence.
        # They should all have been converted by the parser, but if
        # not reversing the record order might expose an error.  Maybe.
        alignments.reverse()
        self.check_simple_write_read(alignments)

    def test_reading_alignments_ig(self):
        path = 'IntelliGenetics/VIF_mase-pro.txt'
        self.check_parse1(path, "ig", 1, 16)
        self.check_parse2(path, "ig", 1)
        self.check_parse3(path, "ig", 1)
        self.check_parse4(path, "ig", 1)
        self.check_parse5(path, "ig", 1)
        alignment = self.check_read(path, "ig", 16, 298)
        self.check_summary2(alignment, [
"  MMMMMMMMMMMMMMMM alignment column 0",
"  EEEEEEETEEEENEEE alignment column 1",
"  NNNNNNNAEEEEQRKK alignment column 2",
"  --------DEEEEE-- alignment column 3",
"  --------KKKKKK-- alignment column 4",
"  |||||||||||||||| ...",
"  HHHHHHH-AAAAL-R- alignment column 297"])
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="ig", seq_count=16))), 3)
        handle.close()

    def test_reading_alignments_pir(self):
        path = 'NBRF/clustalw.pir'
        self.check_parse1(path, "pir", 1, 2)
        self.check_parse2(path, "pir", 1)
        self.check_parse3(path, "pir", 1)
        self.check_parse4(path, "pir", 1)
        self.check_parse5(path, "pir", 1)
        alignment = self.check_read(path, "pir", 2, 2527)
        self.check_summary(alignment, [
"  ------------------------------------...--- 804Angiostrongylus_cantonensis",
"  ------------------------------------...--- 815Parelaphostrongylus_odocoil"])
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        try:
            info_content = summary.information_content()
        except ValueError as err:
            if str(err) != "Error in alphabet: not Nucleotide or Protein, supply expected frequencies":
                raise err

        with open(path, "r") as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(len(list(AlignIO.parse(handle=handle, format="pir", seq_count=2))), 3)
        handle.close()

    def test_reading_alignments_maf1(self):
        path = 'MAF/humor.maf'
        alignments = self.check_parse1(path, "maf", 2, 3)
        self.check_parse2(path, "maf", 2)
        self.check_parse3(path, "maf", 2)
        self.check_parse4(path, "maf", 2)
        self.check_parse5(path, "maf", 2)
        self.check_read_fails(path, "maf")
        self.assertEqual(alignments[0].get_alignment_length(), 5486)
        self.check_summary(alignments[0], [
"  gcacagcctttactccctgactgcgtttatattctg...CCG NM_006987",
"  gcacagcctttactccctgactgcgtttatattctg...TTG mm3",
"  gcacagcctttactccctgactgcgtttatattctg...CCG rn3"])
        self.assertEqual(alignments[1].get_alignment_length(), 5753)
        self.check_summary(alignments[1], [
"  tttgtccatgttggtcaggctggtctcgaactcccc...GGT NM_018289",
"  tttgtccatgttggtcaggctggtctcgaactcccc...GGT mm3",
"  tttgtccatgttggtcaggctggtctcgaactcccc...GGT rn3"])
        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignments[1])
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        alignments.reverse()
        self.check_simple_write_read(alignments)

    def test_reading_alignments_maf2(self):
        path = "MAF/bug2453.maf"
        alignments = self.check_parse1(path, "maf", 3)
        self.check_parse2(path, "maf", 3)
        self.check_parse3(path, "maf", 3)
        self.check_parse4(path, "maf", 3)
        self.check_parse5(path, "maf", 3)
        self.check_read_fails(path, "maf")
        self.assertEqual(len(alignments[0]), 5)
        self.assertEqual(alignments[0].get_alignment_length(), 42)
        self.check_summary2(alignments[0], [
"  AAA-- alignment column 0",
"  AAAAA alignment column 1",
"  AAAAA alignment column 2",
"  ---T- alignment column 3",
"  GGGGG alignment column 4",
"  ||||| ...",
"  GGGGG alignment column 41"])
        self.assertEqual(len(alignments[1]), 5)
        self.assertEqual(alignments[1].get_alignment_length(), 6)
        self.check_summary2(alignments[1], [
"  TTTTt alignment column 0",
"  AAAAa alignment column 1",
"  AAAAa alignment column 2",
"  AAAAg alignment column 3",
"  GGGGg alignment column 4",
"  ||||| ...",
"  AAAAa alignment column 5"])
        self.assertEqual(len(alignments[2]), 4)
        self.assertEqual(alignments[2].get_alignment_length(), 13)
        self.check_summary(alignments[2], [
"  gcagctgaaaaca hg16.chr7",
"  gcagctgaaaaca panTro1.chr6",
"  gcagctgaaaaca baboon",
"  ACAGCTGAAAATA mm4.chr6"])
        summary = AlignInfo.SummaryInfo(alignments[1])
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        alignments.reverse()
        self.check_simple_write_read(alignments)

    def test_reading_alignments_maf3(self):
        path = "MAF/ucsc_test.maf"
        alignments = self.check_parse1(path, "maf", 3)
        self.check_parse2(path, "maf", 3)
        self.check_parse3(path, "maf", 3)
        self.check_parse4(path, "maf", 3)
        self.check_parse5(path, "maf", 3)
        self.check_read_fails(path, "maf")
        self.assertEqual(len(alignments[0]), 5)
        self.assertEqual(alignments[0].get_alignment_length(), 42)
        self.check_summary2(alignments[0], [
"  AAA-- alignment column 0",
"  AAAAA alignment column 1",
"  AAAAA alignment column 2",
"  ---T- alignment column 3",
"  GGGGG alignment column 4",
"  ||||| ...",
"  GGGGG alignment column 41"])
        self.assertEqual(len(alignments[1]), 5)
        self.assertEqual(alignments[1].get_alignment_length(), 6)
        self.check_summary2(alignments[1], [
"  TTTTt alignment column 0",
"  AAAAa alignment column 1",
"  AAAAa alignment column 2",
"  AAAAg alignment column 3",
"  GGGGg alignment column 4",
"  ||||| ...",
"  AAAAa alignment column 5"])
        self.assertEqual(len(alignments[2]), 4)
        self.assertEqual(alignments[2].get_alignment_length(), 13)
        self.check_summary2(alignments[2], [
"  gcagctgaaaaca hg16.chr7",
"  gcagctgaaaaca panTro1.chr6",
"  gcagctgaaaaca baboon",
"  ACAGCTGAAAATA mm4.chr6"])
        summary = AlignInfo.SummaryInfo(alignments[2])
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        try:
            info_content = summary.information_content()
        except ValueError as err:
            if str(err) != "Error in alphabet: not Nucleotide or Protein, supply expected frequencies":
                raise err
        alignments.reverse()
        self.check_simple_write_read(alignments)

    def test_reading_alignments_maf4(self):
        path = "MAF/ucsc_mm9_chr10.maf"
        alignments = self.check_parse1(path, "maf", 48)
        self.check_parse2(path, "maf", 48)
        self.check_parse3(path, "maf", 48)
        self.check_parse4(path, "maf", 48)
        self.check_parse5(path, "maf", 48)
        self.check_read_fails(path, "maf")
        self.assertEqual(len(alignments[0]), 2)
        self.assertEqual(alignments[0].get_alignment_length(), 164)
        self.check_summary(alignments[0], [
"  TCATAGGTATTTATTTTTAAATATGGTTTGCTTTAT...GTT mm9.chr10",
"  TCACAGATATTTACTATTAAATATGGTTTGTTATAT...GTT oryCun1.scaffold_133159"])
        self.assertEqual(len(alignments[1]), 4)
        self.assertEqual(alignments[1].get_alignment_length(), 466)
        self.check_summary(alignments[1], [
"  AGTCTTTCCAATGGGACCTGTGAGTCCTAACTATGC...CTG mm9.chr10",
"  AGTCTTCATAAGTGGAAATATAAGTTTTAATTATTC...TTC ponAbe2.chr6",
"  AGTCTTCATAAGTGGAAATATAAGTTTTAATTATTC...TTC panTro2.chr6",
"  AGTCTTCATAAGTGGAAATATAAGTTTTAATTATTC...TTC hg18.chr6"])
        self.assertEqual(len(alignments[2]), 5)
        self.assertEqual(alignments[2].get_alignment_length(), 127)
        self.check_summary2(alignments[2], [
"  TTTTT alignment column 0",
"  GGGGG alignment column 1",
"  GGGGG alignment column 2",
"  GGGGG alignment column 3",
"  TTTTC alignment column 4",
"  ||||| ...",
"  CCCCC alignment column 126"])
        self.assertEqual(len(alignments[47]), 6)
        self.assertEqual(alignments[47].get_alignment_length(), 46)
        self.check_summary2(alignments[47], [
"  TTTTTT alignment column 0",
"  GGGGGG alignment column 1",
"  TTTTTT alignment column 2",
"  TTTTTT alignment column 3",
"  TGGGAT alignment column 4",
"  |||||| ...",
"  tTTTT- alignment column 45"])
        summary = AlignInfo.SummaryInfo(alignments[47])
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))
        alignments.reverse()
        self.check_simple_write_read(alignments)

    def test_reading_alignments_mauve(self):
        path = 'Mauve/simple.xmfa'
        alignments = self.check_parse1(path, "mauve", 5)
        self.check_parse2(path, "mauve", 5)
        self.check_parse3(path, "mauve", 5)
        self.check_parse4(path, "mauve", 5)
        self.check_parse5(path, "mauve", 5)
        self.check_read_fails(path, "mauve")
        self.assertEqual(len(alignments[0]), 2)
        self.assertEqual(alignments[0].get_alignment_length(), 5670)
        self.check_summary(alignments[0], [
"  ATATTAGGTTTTTACCTACCCAGGAAAAGCCAACCA...AAT 1/0-5670",
"  ATATTAGGTTTTTACCTACCCAGGAAAAGCCAACCA...AAT 2/0-5670"])
        self.assertEqual(len(alignments[1]), 2)
        self.assertEqual(alignments[1].get_alignment_length(), 4420)
        self.check_summary(alignments[1], [
"  GAACATCAGCACCTGAGTTGCTAAAGTCATTTAGAG...CTC 1/5670-9940",
"  GAACATCAGCACCTGAGTTGCTAAAGTCATTTAGAG...CTC 2/7140-11410"])
        self.assertEqual(len(alignments[2]), 1)
        self.assertEqual(alignments[2].get_alignment_length(), 4970)
        self.check_summary(alignments[2], [
"  TCTACCAACCACCACAGACATCAATCACTTCTGCTG...GAC 1/9940-14910"])
        self.assertEqual(len(alignments[3]), 1)
        self.assertEqual(alignments[3].get_alignment_length(), 1470)
        self.assertEqual(len(alignments[4]), 1)
        self.assertEqual(alignments[4].get_alignment_length(), 1470)
        self.check_summary(alignments[4], [
"  ATTCGCACATAAGAATGTACCTTGCTGTAATTTATA...ATA 2/11410-12880"])

        # Check AlignInfo.SummaryInfo likes the alignment
        summary = AlignInfo.SummaryInfo(alignments[4])
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        with self.assertRaises(ValueError) as cm:
            info_content = summary.information_content()
        self.assertEqual("Error in alphabet: not Nucleotide or Protein, supply expected frequencies", str(cm.exception))

        # Some alignment file formats have magic characters which mean
        # use the letter in this position in the first sequence.
        # They should all have been converted by the parser, but if
        # not reversing the record order might expose an error.  Maybe.
        alignments.reverse()
        self.check_simple_write_read(alignments)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
