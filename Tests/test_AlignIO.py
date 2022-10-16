# Copyright 2008 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for AlignIO module."""
import unittest
import warnings

from io import StringIO

from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import AlignInfo
from Bio.Align import MultipleSeqAlignment
from Bio.Data import IUPACData
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

test_write_read_alignment_formats = sorted(AlignIO._FormatToWriter)
test_write_read_align_with_seq_count = test_write_read_alignment_formats + [
    "fasta",
    "tab",
]


class TestAlignIO_exceptions(unittest.TestCase):

    t_formats = list(AlignIO._FormatToWriter) + list(SeqIO._FormatToWriter)

    def test_phylip_reject_duplicate(self):
        """Check that writing duplicated IDs after truncation fails for PHYLIP."""
        handle = StringIO()
        sequences = [
            SeqRecord(Seq("AAAA"), id="longsequencename1"),
            SeqRecord(Seq("AAAA"), id="longsequencename2"),
            SeqRecord(Seq("AAAA"), id="other_sequence"),
        ]
        alignment = MultipleSeqAlignment(sequences)
        with self.assertRaises(ValueError) as cm:
            AlignIO.write(alignment, handle, "phylip")
        self.assertEqual(
            "Repeated name 'longsequen' (originally 'longsequencename2'), possibly due to truncation",
            str(cm.exception),
        )

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
    def simple_alignment_comparison(self, alignments, alignments2, fmt):
        self.assertEqual(len(alignments), len(alignments2))
        for a1, a2 in zip(alignments, alignments2):
            self.assertEqual(a1.get_alignment_length(), a2.get_alignment_length())
            self.assertEqual(len(a1), len(a2))
            for r1, r2 in zip(a1, a2):
                # Check the bare minimum (ID and sequence) as
                # many formats can't store more than that.
                # Check the sequence:
                self.assertEqual(r1.seq, r2.seq)
                # Beware of different quirks and limitations in the
                # valid character sets and the identifier lengths!
                if fmt in ["phylip", "phylip-sequential"]:
                    id1 = r1.id.replace("[", "").replace("]", "")[:10]
                elif fmt == "phylip-relaxed":
                    id1 = r1.id.replace(" ", "").replace(":", "|")
                elif fmt == "clustal":
                    id1 = r1.id.replace(" ", "_")[:30]
                elif fmt in ["stockholm", "maf"]:
                    id1 = r1.id.replace(" ", "_")
                elif fmt == "fasta":
                    id1 = r1.id.split()[0]
                else:
                    id1 = r1.id
                id2 = r2.id
                self.assertEqual(id1, id2)

    def check_reverse_write_read(self, alignments, indent=" "):
        alignments.reverse()
        for fmt in test_write_read_align_with_seq_count:
            records_per_alignment = len(alignments[0])
            for a in alignments:
                if records_per_alignment != len(a):
                    records_per_alignment = None
            # Can we expect this format to work?
            if (
                not records_per_alignment
                and fmt not in test_write_read_alignment_formats
            ):
                continue

            # Going to write to a handle...
            handle = StringIO()

            if fmt == "nexus":
                with self.assertRaises(ValueError) as cm:
                    c = AlignIO.write(alignments, handle=handle, format=fmt)
                self.assertEqual(
                    "We can only write one Alignment to a Nexus file.",
                    str(cm.exception),
                )
                continue
            c = AlignIO.write(alignments, handle=handle, format=fmt)
            self.assertEqual(c, len(alignments))

            # First, try with the seq_count
            if records_per_alignment:
                handle.flush()
                handle.seek(0)
                alignments2 = list(
                    AlignIO.parse(
                        handle=handle, format=fmt, seq_count=records_per_alignment
                    )
                )
                self.simple_alignment_comparison(alignments, alignments2, fmt)

            if fmt in test_write_read_alignment_formats:
                # Don't need the seq_count
                handle.flush()
                handle.seek(0)
                alignments2 = list(AlignIO.parse(handle=handle, format=fmt))
                self.simple_alignment_comparison(alignments, alignments2, fmt)

            # Try writing just one Alignment (not a list)
            handle = StringIO()
            AlignIO.write(alignments[0:1], handle, fmt)
            self.assertEqual(handle.getvalue(), format(alignments[0], fmt))

    def check_iterator_for_loop_handle(self, path, fmt, length, m=None):
        # Try using the iterator with a for loop and a handle
        with open(path) as handle:
            alignments = list(AlignIO.parse(handle, format=fmt))
            self.assertEqual(len(alignments), length)
        if m is not None:
            for alignment in alignments:
                self.assertEqual(len(alignment), m)
        return alignments

    def check_iterator_for_loop_filename(self, path, fmt, length):
        # Try using the iterator with a for loop and a filename not handle
        counter = 0
        for record in AlignIO.parse(path, format=fmt):
            counter += 1
        self.assertEqual(counter, length)

    def check_iterator_next(self, path, fmt, length):
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

    def check_iterator_next_and_list(self, path, fmt, length):
        # Try a mixture of next() and list
        counter = 0
        alignments = AlignIO.parse(path, format=fmt)
        alignment = next(alignments)
        counter = 1
        counter += len(list(alignments))
        self.assertEqual(counter, length)

    def check_iterator_next_for_loop(self, path, fmt, length):
        # Try a mixture of next() and for loop
        alignments = AlignIO.parse(path, format=fmt)
        alignment = next(alignments)
        counter = 1
        for alignment in alignments:
            counter += 1
        self.assertEqual(counter, length)

    def check_write_three_times_and_read(self, path, fmt, m):
        with open(path) as handle:
            data = handle.read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        self.assertEqual(
            len(list(AlignIO.parse(handle=handle, format=fmt, seq_count=m))), 3
        )
        handle.close()

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

    def check_alignment_rows(self, alignment, sequences, column_annotations=None):
        max_len = 40
        items = []
        for record in alignment:
            name = record.id
            sequence = record.seq
            if len(sequence) > max_len:
                sequence = sequence[: max_len - 6] + "..." + sequence[-3:]
            item = (name, sequence)
            items.append(item)
        self.assertEqual(sequences, sorted(items))
        if column_annotations is None:
            self.assertEqual(alignment.column_annotations, {})
        else:
            self.assertEqual(alignment.column_annotations, column_annotations)

    def check_alignment_columns(self, alignment, columns):
        alignment_len = alignment.get_alignment_length()
        # Compare each sequence column
        for index in range(min(5, alignment_len)):
            self.assertEqual(alignment[:, index], columns[index])
        if alignment_len > 5:
            self.assertEqual(alignment[:, -1], columns[-1])

    def check_summary_simple(self, alignment):
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()

    def check_summary(self, alignment, molecule_type):
        # Check AlignInfo.SummaryInfo likes the alignment; smoke test only
        if molecule_type == "DNA":
            letters = IUPACData.unambiguous_dna_letters
            all_letters = IUPACData.ambiguous_dna_letters
        elif molecule_type == "RNA":
            letters = IUPACData.unambiguous_rna_letters
            all_letters = IUPACData.ambiguous_rna_letters
        elif molecule_type == "protein":
            letters = IUPACData.protein_letters
            all_letters = IUPACData.protein_letters
        else:
            raise ValueError(f"Unknown molecule type '{molecule_type}'")
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary(skip_chars=None, letters=letters)
        e_freq = 1.0 / len(letters)
        all_letters = all_letters.upper() + all_letters.lower()
        e_freq_table = dict.fromkeys(all_letters, e_freq)
        info_content = summary.information_content(
            e_freq_table=e_freq_table, chars_to_ignore=["N", "X"]
        )

    def check_summary_pir(self, alignment):
        letters = IUPACData.unambiguous_dna_letters
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()
        # gap_consensus = summary.gap_consensus()
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary(skip_chars=None, letters=letters)
        e_freq = 1.0 / len(letters)
        all_letters = letters.upper() + letters.lower()
        e_freq_table = dict.fromkeys(all_letters, e_freq)
        info_content = summary.information_content(e_freq_table=e_freq_table)

    def test_reading_alignments_clustal1(self):
        path = "Clustalw/clustalw.aln"
        self.check_iterator_for_loop_handle(path, "clustal", 1, 2)
        self.check_iterator_for_loop_filename(path, "clustal", 1)
        self.check_iterator_next(path, "clustal", 1)
        self.check_iterator_next_and_list(path, "clustal", 1)
        self.check_iterator_next_for_loop(path, "clustal", 1)
        self.check_write_three_times_and_read(path, "clustal", 2)
        alignment = self.check_read(path, "clustal", 2, 601)
        self.check_alignment_rows(
            alignment,
            [
                (
                    "gi|4959044|gb|AAD34209.1|AF069",
                    "MENSDSNDKGSDQSAAQRRSQMDRLDREEAFYQF...SVV",
                ),
                (
                    "gi|671626|emb|CAA85685.1|",
                    "---------MSPQTETKASVGFKAGVKEYKLTYY...---",
                ),
            ],
            {
                "clustal_consensus": "          * *: ::    :.   :*  :  :. : . :*  ::   .:   **                  **:...   *.*** ..          .:*   * *: .* :*        : :* .*                   *::.  .    .:: :*..*  :* .*   .. .  :    .  :    *. .:: : .      .* .  :  *.:     ..::   * .  ::  :  .*.    :.    :. .  .  .* **.*..  :..  *.. .    . ::*                         :.: .*:    :     * ::   ***  . * :. .  .  :  *: .:: :::   ..   . : :   ::  *    *  : .. :.* . ::.  :: * :  :   * *   :..  * ..  * :**                             .  .:. ..   :*.  ..: :. .  .:* * :   : * .             ..*:.  .**   *.*... :  ::   :* .*  ::* : :.  :.    :   "
            },
        )
        self.check_summary(alignment, "protein")

    def test_reading_alignments_clustal2(self):
        path = "Clustalw/opuntia.aln"
        self.check_iterator_for_loop_handle(path, "clustal", 1, 7)
        self.check_iterator_for_loop_filename(path, "clustal", 1)
        self.check_iterator_next(path, "clustal", 1)
        self.check_iterator_next_and_list(path, "clustal", 1)
        self.check_iterator_next_for_loop(path, "clustal", 1)
        self.check_write_three_times_and_read(path, "clustal", 7)
        alignment = self.check_read(path, "clustal", 7, 156)
        self.check_alignment_columns(
            alignment,
            ["TTTTTTT", "AAAAAAA", "TTTTTTT", "AAAAAAA", "CCCCCCC", "AAAAAAA"],
        )
        self.check_summary(alignment, "DNA")

    def test_reading_alignments_clustal3(self):
        path = "Clustalw/hedgehog.aln"
        self.check_iterator_for_loop_handle(path, "clustal", 1, 5)
        self.check_iterator_for_loop_filename(path, "clustal", 1)
        self.check_iterator_next(path, "clustal", 1)
        self.check_iterator_next_and_list(path, "clustal", 1)
        self.check_iterator_next_for_loop(path, "clustal", 1)
        self.check_write_three_times_and_read(path, "clustal", 5)
        alignment = self.check_read(path, "clustal", 5, 447)
        self.check_alignment_columns(
            alignment, ["M----", "F----", "N----", "L----", "V----", "---SS"]
        )
        self.check_summary(alignment, "protein")

    def test_reading_alignments_clustal4(self):
        path = "Clustalw/odd_consensus.aln"
        self.check_iterator_for_loop_handle(path, "clustal", 1, 2)
        self.check_iterator_for_loop_filename(path, "clustal", 1)
        self.check_iterator_next(path, "clustal", 1)
        self.check_iterator_next_and_list(path, "clustal", 1)
        self.check_iterator_next_for_loop(path, "clustal", 1)
        self.check_write_three_times_and_read(path, "clustal", 2)
        alignment = self.check_read(path, "clustal", 2, 687)
        self.check_alignment_rows(
            alignment,
            [
                ("AT3G20900.1-CDS", "----------------------------------...TAG"),
                ("AT3G20900.1-SEQ", "ATGAACAAAGTAGCGAGGAAGAACAAAACATCAG...TAG"),
            ],
            {
                "clustal_consensus": "                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            *       *  *** ***** *   *  **      *******************************************************************************************************************************************************************************"
            },
        )
        self.check_summary(alignment, "DNA")

    def test_reading_alignments_clustal5(self):
        path = "Clustalw/protein.aln"
        self.check_iterator_for_loop_handle(path, "clustal", 1, 20)
        self.check_iterator_for_loop_filename(path, "clustal", 1)
        self.check_iterator_next(path, "clustal", 1)
        self.check_iterator_next_and_list(path, "clustal", 1)
        self.check_iterator_next_for_loop(path, "clustal", 1)
        self.check_write_three_times_and_read(path, "clustal", 20)
        alignment = self.check_read(path, "clustal", 20, 411)
        self.check_alignment_columns(
            alignment,
            [
                "-M------------------",
                "-T------------------",
                "-V------------------",
                "-L-----------------M",
                "-E---------------MMS",
                "-------------------T",
            ],
        )
        self.check_summary(alignment, "protein")

    def test_reading_alignments_clustal6(self):
        path = "Clustalw/promals3d.aln"
        self.check_iterator_for_loop_handle(path, "clustal", 1, 20)
        self.check_iterator_for_loop_filename(path, "clustal", 1)
        self.check_iterator_next(path, "clustal", 1)
        self.check_iterator_next_and_list(path, "clustal", 1)
        self.check_iterator_next_for_loop(path, "clustal", 1)
        self.check_write_three_times_and_read(path, "clustal", 20)
        alignment = self.check_read(path, "clustal", 20, 414)
        self.check_alignment_columns(
            alignment,
            [
                "MMMMMMMMMMMMMMMM-M--",
                "-----------------T--",
                "-----------------V--",
                "-----------------L--",
                "-S---------------E--",
                "-T------------------",
            ],
        )
        self.check_summary(alignment, "protein")

    def test_reading_alignments_fasta(self):
        path = "GFF/multi.fna"  # Trivial nucleotide alignment
        self.check_iterator_for_loop_handle(path, "fasta", 1, 3)
        self.check_iterator_for_loop_filename(path, "fasta", 1)
        self.check_iterator_next(path, "fasta", 1)
        self.check_iterator_next_and_list(path, "fasta", 1)
        self.check_iterator_next_for_loop(path, "fasta", 1)
        self.check_write_three_times_and_read(path, "fasta", 3)
        alignment = self.check_read(path, "fasta", 3, 8)
        self.check_alignment_rows(
            alignment,
            [("test1", "ACGTCGCG"), ("test2", "GGGGCCCC"), ("test3", "AAACACAC")],
        )
        self.check_summary(alignment, "DNA")

    def test_reading_alignments_nexus1(self):
        path = "Nexus/test_Nexus_input.nex"
        self.check_iterator_for_loop_handle(path, "nexus", 1, 9)
        self.check_iterator_for_loop_filename(path, "nexus", 1)
        self.check_iterator_next(path, "nexus", 1)
        self.check_iterator_next_and_list(path, "nexus", 1)
        self.check_iterator_next_for_loop(path, "nexus", 1)
        alignment = self.check_read(path, "nexus", 9, 48)
        self.check_alignment_columns(
            alignment,
            [
                "AAAAAAAAc",
                "-----c?tc",
                "CCCCCCCCc",
                "--c-?a-tc",
                "GGGGGGGGc",
                "tt--?ag?c",
            ],
        )
        self.check_summary_simple(alignment)

    def test_reading_alignments_nexus2(self):
        path = "Nexus/codonposset.nex"
        self.check_iterator_for_loop_handle(path, "nexus", 1, 2)
        self.check_iterator_for_loop_filename(path, "nexus", 1)
        self.check_iterator_next(path, "nexus", 1)
        self.check_iterator_next_and_list(path, "nexus", 1)
        self.check_iterator_next_for_loop(path, "nexus", 1)
        alignment = self.check_read(path, "nexus", 2, 22)
        self.check_alignment_rows(
            alignment,
            [
                ("Aegotheles", "AAAAAGGCATTGTGGTGGGAAT"),
                ("Aerodramus", "?????????TTGTGGTGGGAAT"),
            ],
        )
        self.check_summary_simple(alignment)

    def test_reading_alignments_msf1(self):
        path = "msf/DOA_prot.msf"
        with self.assertRaisesRegex(
            ValueError,
            "GCG MSF header said alignment length 62, "
            "but 11 of 12 sequences said Len: 250",
        ):
            AlignIO.read(path, "msf")

    def test_reading_alignments_msf2(self):
        path = "msf/W_prot.msf"
        with warnings.catch_warnings(record=True) as w:
            self.check_iterator_for_loop_handle(path, "msf", 1, 11)
            self.check_iterator_for_loop_filename(path, "msf", 1)
            self.check_iterator_next(path, "msf", 1)
            self.check_iterator_next_and_list(path, "msf", 1)
            self.check_iterator_next_for_loop(path, "msf", 1)
            alignment = self.check_read(path, "msf", 11, 99)
        warning_msgs = {str(_.message) for _ in w}
        self.assertIn(
            "One of more alignment sequences were truncated and have been gap padded",
            warning_msgs,
        )
        self.check_alignment_columns(
            alignment,
            [
                "GGGGGGGGGGG",
                "LLLLLLLLLLL",
                "TTTTTTTTTTT",
                "PPPPPPPPPPP",
                "FFFFFFSSSSS",
                # ...
                "LLLLLL----L",
            ],
        )
        self.check_summary_simple(alignment)

    def test_reading_alignments_stockholm1(self):
        path = "Stockholm/simple.sth"
        self.check_iterator_for_loop_handle(path, "stockholm", 1, 2)
        self.check_iterator_for_loop_filename(path, "stockholm", 1)
        self.check_iterator_next(path, "stockholm", 1)
        self.check_iterator_next_and_list(path, "stockholm", 1)
        self.check_iterator_next_for_loop(path, "stockholm", 1)
        self.check_write_three_times_and_read(path, "stockholm", 2)
        alignment = self.check_read(path, "stockholm", 2, 104)
        self.check_alignment_rows(
            alignment,
            [
                ("AE007476.1", "AAAAUUGAAUAUCGUUUUACUUGUUUAU-GUCGU...GAU"),
                ("AP001509.1", "UUAAUCGAGCUCAACACUCUUCGUAUAUCCUC-U...UGU"),
            ],
            {
                "secondary_structure": ".................<<<<<<<<...<<<<<<<........>>>>>>>........<<<<<<<.......>>>>>>>..>>>>>>>>..............."
            },
        )
        self.check_summary(alignment, "RNA")

    def test_reading_alignments_stockholm2(self):
        path = "Stockholm/funny.sth"
        self.check_iterator_for_loop_handle(path, "stockholm", 1, 6)
        self.check_iterator_for_loop_filename(path, "stockholm", 1)
        self.check_iterator_next(path, "stockholm", 1)
        self.check_iterator_next_and_list(path, "stockholm", 1)
        self.check_iterator_next_for_loop(path, "stockholm", 1)
        self.check_write_three_times_and_read(path, "stockholm", 6)
        alignment = self.check_read(path, "stockholm", 6, 43)
        self.check_alignment_columns(
            alignment, ["MMMEEE", "TQIVVV", "CHEMMM", "RVALLL", "ASDTTT", "SYSEEE"]
        )
        self.check_summary(alignment, "protein")

    def test_reading_alignments_phylip1(self):
        path = "Phylip/reference_dna.phy"
        self.check_iterator_for_loop_handle(path, "phylip", 1, 6)
        self.check_iterator_for_loop_filename(path, "phylip", 1)
        self.check_iterator_next(path, "phylip", 1)
        self.check_iterator_next_and_list(path, "phylip", 1)
        self.check_iterator_next_for_loop(path, "phylip", 1)
        self.check_write_three_times_and_read(path, "phylip", 6)
        alignment = self.check_read(path, "phylip", 6, 13)
        self.check_alignment_columns(
            alignment, ["CCTTCG", "GGAAAG", "ATAAAC", "TTTTAA", "GAGGAG", "CTTTTC"]
        )
        self.check_summary(alignment, "DNA")

    def test_reading_alignments_phylip2(self):
        path = "Phylip/reference_dna2.phy"
        self.check_iterator_for_loop_handle(path, "phylip", 1, 6)
        self.check_iterator_for_loop_filename(path, "phylip", 1)
        self.check_iterator_next(path, "phylip", 1)
        self.check_iterator_next_and_list(path, "phylip", 1)
        self.check_iterator_next_for_loop(path, "phylip", 1)
        self.check_write_three_times_and_read(path, "phylip", 6)
        alignment = self.check_read(path, "phylip", 6, 39)
        self.check_alignment_columns(
            alignment, ["CCTTCG", "GGAAAG", "ATAAAC", "TTTTAA", "GAGGAG", "CTTTTC"]
        )
        self.check_summary(alignment, "DNA")

    def test_reading_alignments_phylip3(self):
        path = "Phylip/hennigian.phy"
        self.check_iterator_for_loop_handle(path, "phylip", 1, 10)
        self.check_iterator_for_loop_filename(path, "phylip", 1)
        self.check_iterator_next(path, "phylip", 1)
        self.check_iterator_next_and_list(path, "phylip", 1)
        self.check_iterator_next_for_loop(path, "phylip", 1)
        self.check_write_three_times_and_read(path, "phylip", 10)
        alignment = self.check_read(path, "phylip", 10, 40)
        self.check_alignment_columns(
            alignment,
            [
                "CCCCCAAAAA",
                "AAAAACCCCC",
                "CCCAAAAAAA",
                "AAACCAAAAA",
                "CCAAAAAAAA",
                "AAAAAAAAAA",
            ],
        )
        self.check_summary(alignment, "DNA")

    def test_reading_alignments_phylip4(self):
        path = "Phylip/horses.phy"
        self.check_iterator_for_loop_handle(path, "phylip", 1, 10)
        self.check_iterator_for_loop_filename(path, "phylip", 1)
        self.check_iterator_next(path, "phylip", 1)
        self.check_iterator_next_and_list(path, "phylip", 1)
        self.check_iterator_next_for_loop(path, "phylip", 1)
        self.check_write_three_times_and_read(path, "phylip", 10)
        alignment = self.check_read(path, "phylip", 10, 40)
        self.check_alignment_columns(
            alignment,
            [
                "AACCCCCCCC",
                "AAAACCCCCC",
                "AAAAAAAAAC",
                "ACAAAAAAAA",
                "ACACCCCCCC",
                "AAAAAAAAAA",
            ],
        )
        self.check_summary(alignment, "DNA")

    def test_reading_alignments_phylip5(self):
        path = "Phylip/random.phy"
        self.check_iterator_for_loop_handle(path, "phylip", 1, 10)
        self.check_iterator_for_loop_filename(path, "phylip", 1)
        self.check_iterator_next(path, "phylip", 1)
        self.check_iterator_next_and_list(path, "phylip", 1)
        self.check_iterator_next_for_loop(path, "phylip", 1)
        self.check_write_three_times_and_read(path, "phylip", 10)
        alignment = self.check_read(path, "phylip", 10, 40)
        self.check_alignment_columns(
            alignment,
            [
                "CAAAACAAAC",
                "AACAACCACC",
                "CAAAACAAAA",
                "ACAACACACA",
                "CCAAAACCAA",
                "AAAAAAAAAA",
            ],
        )
        self.check_summary(alignment, "DNA")

    def test_reading_alignments_phylip6(self):
        path = "Phylip/interlaced.phy"
        self.check_iterator_for_loop_handle(path, "phylip", 1, 3)
        self.check_iterator_for_loop_filename(path, "phylip", 1)
        self.check_iterator_next(path, "phylip", 1)
        self.check_iterator_next_and_list(path, "phylip", 1)
        self.check_iterator_next_for_loop(path, "phylip", 1)
        self.check_write_three_times_and_read(path, "phylip", 3)
        alignment = self.check_read(path, "phylip", 3, 384)
        self.check_alignment_rows(
            alignment,
            [
                ("ALEU_HORVU", "MAHARVLLLALAVLATAAVAVASSSSFADSNPIR...VAA"),
                ("CATH_HUMAN", "------MWATLPLLCAGAWLLGV--------PVC...PLV"),
                ("CYS1_DICDI", "-----MKVILLFVLAVFTVFVSS-----------...I--"),
            ],
        )
        self.check_summary(alignment, "protein")

    def test_reading_alignments_phylip7(self):
        path = "Phylip/interlaced2.phy"
        self.check_iterator_for_loop_handle(path, "phylip", 1, 4)
        self.check_iterator_for_loop_filename(path, "phylip", 1)
        self.check_iterator_next(path, "phylip", 1)
        self.check_iterator_next_and_list(path, "phylip", 1)
        self.check_iterator_next_for_loop(path, "phylip", 1)
        self.check_write_three_times_and_read(path, "phylip", 4)
        alignment = self.check_read(path, "phylip", 4, 131)
        self.check_alignment_rows(
            alignment,
            [
                ("IXI_234", "TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRP...SHE"),
                ("IXI_235", "TSPASIRPPAGPSSR---------RPSPPGPRRP...SHE"),
                ("IXI_236", "TSPASIRPPAGPSSRPAMVSSR--RPSPPPPRRP...SHE"),
                ("IXI_237", "TSPASLRPPAGPSSRPAMVSSRR-RPSPPGPRRP...SHE"),
            ],
        )
        self.check_summary(alignment, "protein")

    def test_reading_alignments_phylip8(self):
        path = "ExtendedPhylip/primates.phyx"
        self.check_iterator_for_loop_handle(path, "phylip-relaxed", 1, 12)
        self.check_iterator_for_loop_filename(path, "phylip-relaxed", 1)
        self.check_iterator_next(path, "phylip-relaxed", 1)
        self.check_iterator_next_and_list(path, "phylip-relaxed", 1)
        self.check_iterator_next_for_loop(path, "phylip-relaxed", 1)
        self.check_write_three_times_and_read(path, "phylip-relaxed", 12)
        alignment = self.check_read(path, "phylip-relaxed", 12, 898)
        self.check_alignment_columns(
            alignment,
            [
                "AAAAAAAAAAAA",
                "AAAAAAAAAAAA",
                "GGGGGGGGGGGG",
                "TCCCCCCCCCCC",
                "TTTTTTTTTTTT",
                "TTTTTTTTTTTT",
            ],
        )
        self.check_summary(alignment, "DNA")

    def test_reading_alignments_phylip9(self):
        path = "Phylip/sequential.phy"
        self.check_iterator_for_loop_handle(path, "phylip-sequential", 1, 3)
        self.check_iterator_for_loop_filename(path, "phylip-sequential", 1)
        self.check_iterator_next(path, "phylip-sequential", 1)
        self.check_iterator_next_and_list(path, "phylip-sequential", 1)
        self.check_iterator_next_for_loop(path, "phylip-sequential", 1)
        self.check_write_three_times_and_read(path, "phylip-sequential", 3)
        alignment = self.check_read(path, "phylip-sequential", 3, 384)
        self.check_alignment_rows(
            alignment,
            [
                ("ALEU_HORVU", "MAHARVLLLALAVLATAAVAVASSSSFADSNPIR...VAA"),
                ("CATH_HUMAN", "------MWATLPLLCAGAWLLGV--------PVC...PLV"),
                ("CYS1_DICDI", "-----MKVILLFVLAVFTVFVSS-----------...I--"),
            ],
        )
        self.check_summary(alignment, "protein")

    def test_reading_alignments_phylip10(self):
        path = "Phylip/sequential2.phy"
        self.check_iterator_for_loop_handle(path, "phylip-sequential", 1, 4)
        self.check_iterator_for_loop_filename(path, "phylip-sequential", 1)
        self.check_iterator_next(path, "phylip-sequential", 1)
        self.check_iterator_next_and_list(path, "phylip-sequential", 1)
        self.check_iterator_next_for_loop(path, "phylip-sequential", 1)
        self.check_write_three_times_and_read(path, "phylip-sequential", 4)
        alignment = self.check_read(path, "phylip-sequential", 4, 131)
        self.check_alignment_rows(
            alignment,
            [
                ("IXI_234", "TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRP...SHE"),
                ("IXI_235", "TSPASIRPPAGPSSR---------RPSPPGPRRP...SHE"),
                ("IXI_236", "TSPASIRPPAGPSSRPAMVSSR--RPSPPPPRRP...SHE"),
                ("IXI_237", "TSPASLRPPAGPSSRPAMVSSRR-RPSPPGPRRP...SHE"),
            ],
        )
        self.check_summary(alignment, "protein")

    def test_reading_alignments_emboss1(self):
        path = "Emboss/alignret.txt"
        self.check_iterator_for_loop_handle(path, "emboss", 1, 4)
        self.check_iterator_for_loop_filename(path, "emboss", 1)
        self.check_iterator_next(path, "emboss", 1)
        self.check_iterator_next_and_list(path, "emboss", 1)
        self.check_iterator_next_for_loop(path, "emboss", 1)
        alignment = self.check_read(path, "emboss", 4, 131)
        self.check_alignment_rows(
            alignment,
            [
                ("IXI_234", "TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRP...SHE"),
                ("IXI_235", "TSPASIRPPAGPSSR---------RPSPPGPRRP...SHE"),
                ("IXI_236", "TSPASIRPPAGPSSRPAMVSSR--RPSPPPPRRP...SHE"),
                ("IXI_237", "TSPASLRPPAGPSSRPAMVSSRR-RPSPPGPRRP...SHE"),
            ],
        )
        self.check_summary(alignment, "protein")

    def test_reading_alignments_emboss2(self):
        path = "Emboss/needle.txt"
        alignments = self.check_iterator_for_loop_handle(path, "emboss", 5, 2)
        self.check_iterator_for_loop_filename(path, "emboss", 5)
        self.check_iterator_next(path, "emboss", 5)
        self.check_iterator_next_and_list(path, "emboss", 5)
        self.check_iterator_next_for_loop(path, "emboss", 5)
        self.check_read_fails(path, "emboss")
        # Show the alignment
        self.assertEqual(alignments[0].get_alignment_length(), 124)
        self.assertEqual(alignments[1].get_alignment_length(), 119)
        self.assertEqual(alignments[2].get_alignment_length(), 120)
        self.assertEqual(alignments[3].get_alignment_length(), 118)
        self.assertEqual(alignments[4].get_alignment_length(), 125)
        self.check_alignment_rows(
            alignments[0],
            [
                ("gi|94968718|receiver", "-VLLADDHALVRRGFRLMLED--DPEIEIVAEAG...GET"),
                ("ref_rec", "KILIVDD----QYGIRILLNEVFNKEGYQTFQAA...---"),
            ],
        )
        self.check_alignment_rows(
            alignments[1],
            [
                ("gi|94968761|receiver", "-ILIVDDEANTLASLSRAFRLAGHEATVCDNAVR...LKR"),
                ("ref_rec", "KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQ...---"),
            ],
        )
        self.check_alignment_rows(
            alignments[2],
            [
                ("gi|94967506|receiver", "LHIVVVDDDPGTCVYIESVFAELGHTCKSFVRPE...HKE"),
                ("ref_rec", "-KILIVDDQYGIRILLNEVFNKEGYQTFQAANGL...---"),
            ],
        )
        self.check_alignment_rows(
            alignments[3],
            [
                ("gi|94970045|receiver", "-VLLVEDEEALRAAAGDFLETRGYKIMTARDGTE...EVL"),
                ("ref_rec", "KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQ...DAV"),
            ],
        )
        self.check_alignment_rows(
            alignments[4],
            [
                ("gi|94970041|receiver", "TVLLVEDEEGVRKLVRGILSRQGYHVLEATSGEE...KRQ"),
                ("ref_rec", "KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQ...---"),
            ],
        )
        self.check_summary(alignments[0], "protein")
        self.check_reverse_write_read(alignments)

    def test_reading_alignments_emboss3(self):
        path = "Emboss/needle_asis.txt"
        self.check_iterator_for_loop_handle(path, "emboss", 1, 2)
        self.check_iterator_for_loop_filename(path, "emboss", 1)
        self.check_iterator_next(path, "emboss", 1)
        self.check_iterator_next_and_list(path, "emboss", 1)
        self.check_iterator_next_for_loop(path, "emboss", 1)
        alignment = self.check_read(path, "emboss", 2, 3653)
        self.check_alignment_rows(
            alignment,
            [
                ("asis", "----------------------------------...GAA"),
                ("asis", "TATTTTTTGGATTTTTTTCTAGATTTTCTAGGTT...GAA"),
            ],
        )
        self.check_summary(alignment, "DNA")

    def test_reading_alignments_emboss4(self):
        path = "Emboss/water.txt"
        self.check_iterator_for_loop_handle(path, "emboss", 1, 2)
        self.check_iterator_for_loop_filename(path, "emboss", 1)
        self.check_iterator_next(path, "emboss", 1)
        self.check_iterator_next_and_list(path, "emboss", 1)
        self.check_iterator_next_for_loop(path, "emboss", 1)
        alignment = self.check_read(path, "emboss", 2, 131)
        self.check_alignment_rows(
            alignment,
            [
                ("IXI_234", "TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRP...SHE"),
                ("IXI_235", "TSPASIRPPAGPSSR---------RPSPPGPRRP...SHE"),
            ],
        )
        self.check_summary(alignment, "protein")

    def test_reading_alignments_emboss5(self):
        path = "Emboss/water2.txt"
        self.check_iterator_for_loop_handle(path, "emboss", 1, 2)
        self.check_iterator_for_loop_filename(path, "emboss", 1)
        self.check_iterator_next(path, "emboss", 1)
        self.check_iterator_next_and_list(path, "emboss", 1)
        self.check_iterator_next_for_loop(path, "emboss", 1)
        alignment = self.check_read(path, "emboss", 2, 18)
        self.check_alignment_rows(
            alignment, [("asis", "CGTTTGAGT-CTGGGATG"), ("asis", "CGTTTGAGTACTGGGATG")]
        )
        self.check_summary(alignment, "DNA")

    def test_reading_alignments_emboss6(self):
        path = "Emboss/matcher_simple.txt"
        self.check_iterator_for_loop_handle(path, "emboss", 1, 2)
        self.check_iterator_for_loop_filename(path, "emboss", 1)
        self.check_iterator_next(path, "emboss", 1)
        self.check_iterator_next_and_list(path, "emboss", 1)
        self.check_iterator_next_for_loop(path, "emboss", 1)
        alignment = self.check_read(path, "emboss", 2, 16)
        self.check_alignment_rows(
            alignment,
            [("AF069992_1", "GPPPQSPDENRAGESS"), ("CAA85685.1", "GVPPEEAGAAVAAESS")],
        )
        self.check_summary(alignment, "protein")

    def test_reading_alignments_emboss7(self):
        path = "Emboss/matcher_pair.txt"
        alignments = self.check_iterator_for_loop_handle(path, "emboss", 5, 2)
        self.check_iterator_for_loop_filename(path, "emboss", 5)
        self.check_iterator_next(path, "emboss", 5)
        self.check_iterator_next_and_list(path, "emboss", 5)
        self.check_iterator_next_for_loop(path, "emboss", 5)
        self.check_read_fails(path, "emboss")
        self.assertEqual(alignments[0].get_alignment_length(), 145)
        self.assertEqual(alignments[1].get_alignment_length(), 13)
        self.assertEqual(alignments[2].get_alignment_length(), 18)
        self.assertEqual(alignments[3].get_alignment_length(), 10)
        self.assertEqual(alignments[4].get_alignment_length(), 10)
        self.check_alignment_rows(
            alignments[0],
            [
                ("HBA_HUMAN", "LSPADKTNVKAAWGKVGAHAGEYGAEALERMFLS...SKY"),
                ("HBB_HUMAN", "LTPEEKSAVTALWGKV--NVDEVGGEALGRLLVV...HKY"),
            ],
        )
        self.check_alignment_rows(
            alignments[1],
            [("HBA_HUMAN", "KKVADALTNAVAH"), ("HBB_HUMAN", "QKVVAGVANALAH")],
        )
        self.check_alignment_rows(
            alignments[2],
            [("HBA_HUMAN", "KLRVDPVNFKLLSHCLLV"), ("HBB_HUMAN", "KVNVDEVGGEALGRLLVV")],
        )
        self.check_alignment_rows(
            alignments[3], [("HBA_HUMAN", "LSALSDLHAH"), ("HBB_HUMAN", "LGAFSDGLAH")]
        )
        self.check_alignment_rows(
            alignments[4], [("HBA_HUMAN", "VKAAWGKVGA"), ("HBB_HUMAN", "VQAAYQKVVA")]
        )
        self.check_summary(alignments[0], "protein")
        self.check_reverse_write_read(alignments)

    def test_reading_alignments_emboss8(self):
        path = "Emboss/emboss_pair_aln_full_blank_line.txt"
        self.check_iterator_for_loop_handle(path, "emboss", 1, 2)
        self.check_iterator_for_loop_filename(path, "emboss", 1)
        self.check_iterator_next(path, "emboss", 1)
        self.check_iterator_next_and_list(path, "emboss", 1)
        self.check_iterator_next_for_loop(path, "emboss", 1)
        alignment = self.check_read(path, "emboss", 2, 1450)
        self.check_alignment_rows(
            alignment,
            [
                (
                    "hg38_chrX_131691529_131830643_47210_48660",
                    "GGCAGGTGCATAGCTTGAGCCTAGGAGTTCAAGT...AAA",
                ),
                (
                    "mm10_chrX_50555743_50635321_27140_27743",
                    "G--------------------------TTCAAGG...AAA",
                ),
            ],
        )
        self.check_summary(alignment, "DNA")

    def test_reading_alignments_fasta_m10_1(self):
        path = "Fasta/output001.m10"
        alignments = self.check_iterator_for_loop_handle(path, "fasta-m10", 4, 2)
        self.check_iterator_for_loop_filename(path, "fasta-m10", 4)
        self.check_iterator_next(path, "fasta-m10", 4)
        self.check_iterator_next_and_list(path, "fasta-m10", 4)
        self.check_iterator_next_for_loop(path, "fasta-m10", 4)
        self.check_read_fails(path, "fasta-m10")
        self.assertEqual(alignments[0].get_alignment_length(), 108)
        self.check_alignment_rows(
            alignments[0],
            [
                (
                    "gi|10955263|ref|NP_052604.1|",
                    "SGSNT-RRRAISRPVRLTAEED---QEIRKRAAE...LSR",
                ),
                (
                    "gi|152973457|ref|YP_001338508.1|",
                    "AGSGAPRRRGSGLASRISEQSEALLQEAAKHAAE...LSR",
                ),
            ],
        )
        self.assertEqual(alignments[1].get_alignment_length(), 64)
        self.check_alignment_rows(
            alignments[1],
            [
                (
                    "gi|10955263|ref|NP_052604.1|",
                    "AAECGKTVSGFLRAAALGKKVNSLTDDRVLKEV-...AIT",
                ),
                (
                    "gi|152973588|ref|YP_001338639.1|",
                    "ASRQGCTVGG--KMDSVQDKASDKDKERVMKNIN...TLT",
                ),
            ],
        )
        self.assertEqual(alignments[2].get_alignment_length(), 38)
        self.check_alignment_rows(
            alignments[2],
            [
                (
                    "gi|10955264|ref|NP_052605.1|",
                    "MKKDKKYQIEAIKNKDKTLFIVYATDIYSPSEFFSKIE",
                ),
                (
                    "gi|152973462|ref|YP_001338513.1|",
                    "IKKDLGVSFLKLKNREKTLIVDALKKKYPVAELLSVLQ",
                ),
            ],
        )
        self.assertEqual(alignments[3].get_alignment_length(), 43)
        self.check_alignment_rows(
            alignments[3],
            [
                (
                    "gi|10955265|ref|NP_052606.1|",
                    "SELHSKLPKSIDKIHEDIKKQLSC-SLIMKKIDV...TYC",
                ),
                (
                    "gi|152973545|ref|YP_001338596.1|",
                    "SRINSDVARRIPGIHRDPKDRLSSLKQVEEALDM...EYC",
                ),
            ],
        )
        self.check_summary(alignments[0], "protein")
        self.check_reverse_write_read(alignments)

    def test_reading_alignments_fasta_m10_2(self):
        path = "Fasta/output002.m10"
        alignments = self.check_iterator_for_loop_handle(path, "fasta-m10", 6, 2)
        self.check_iterator_for_loop_filename(path, "fasta-m10", 6)
        self.check_iterator_next(path, "fasta-m10", 6)
        self.check_iterator_next_and_list(path, "fasta-m10", 6)
        self.check_iterator_next_for_loop(path, "fasta-m10", 6)
        self.check_read_fails(path, "fasta-m10")
        self.assertEqual(alignments[0].get_alignment_length(), 88)
        self.check_alignment_rows(
            alignments[0],
            [
                (
                    "gi|10955263|ref|NP_052604.1|",
                    "SGSNTRRRAISRPVR--LTAEEDQEIRKRAAECG...AEV",
                ),
                (
                    "gi|162139799|ref|NP_309634.2|",
                    "SQRSTRRKPENQPTRVILFNKPYDVLPQFTDEAG...VQV",
                ),
            ],
        )
        self.assertEqual(alignments[1].get_alignment_length(), 53)
        self.check_alignment_rows(
            alignments[1],
            [
                (
                    "gi|10955263|ref|NP_052604.1|",
                    "EIRKRAAECGKTVSGFLRAAA-LGKKV----NSL...KKL",
                ),
                (
                    "gi|15831859|ref|NP_310632.1|",
                    "EIKPRGTSKGEAIAAFMQEAPFIGRTPVFLGDDL...VKI",
                ),
            ],
        )
        self.assertEqual(alignments[2].get_alignment_length(), 92)
        self.check_alignment_rows(
            alignments[2],
            [
                (
                    "gi|10955264|ref|NP_052605.1|",
                    "SEFFSKIESDLKKKKSKGDVFFDLIIPNG-----...ATS",
                ),
                (
                    "gi|15829419|ref|NP_308192.1|",
                    "TELNSELAKAMKVDAQRG-AFVSQVLPNSSAAKA...QSS",
                ),
            ],
        )
        self.assertEqual(alignments[5].get_alignment_length(), 157)
        self.check_alignment_rows(
            alignments[5],
            [
                (
                    "gi|10955265|ref|NP_052606.1|",
                    "QYIMTTSNGDRVRAKIYKRGSIQFQGKYLQIASL...REI",
                ),
                (
                    "gi|15833861|ref|NP_312634.1|",
                    "EFIRLLSDHDQFEKDQISELTVAANALKLEVAK-...KKV",
                ),
            ],
        )
        self.check_summary(alignments[0], "protein")
        self.check_reverse_write_read(alignments)

    def test_reading_alignments_fasta_m10_3(self):
        path = "Fasta/output003.m10"
        alignments = self.check_iterator_for_loop_handle(path, "fasta-m10", 3, 2)
        self.check_iterator_for_loop_filename(path, "fasta-m10", 3)
        self.check_iterator_next(path, "fasta-m10", 3)
        self.check_iterator_next_and_list(path, "fasta-m10", 3)
        self.check_iterator_next_for_loop(path, "fasta-m10", 3)
        self.check_read_fails(path, "fasta-m10")
        self.assertEqual(alignments[0].get_alignment_length(), 55)
        self.check_alignment_rows(
            alignments[0],
            [
                (
                    "gi|10955263|ref|NP_052604.1|",
                    "VRLTAEEDQ--EIRKRAAECG-KTVSGFLRAAAL...LGA",
                ),
                (
                    "gi|152973837|ref|YP_001338874.1|",
                    "ISISNNKDQYEELQKEQGERDLKTVDQLVRIAAA...IAA",
                ),
            ],
        )
        self.assertEqual(alignments[1].get_alignment_length(), 22)
        self.check_alignment_rows(
            alignments[1],
            [
                ("gi|10955265|ref|NP_052606.1|", "DDRANLFEFLSEEGITITEDNN"),
                ("gi|152973840|ref|YP_001338877.1|", "DDAEHLFRTLSSR-LDALQDGN"),
            ],
        )
        self.assertEqual(alignments[2].get_alignment_length(), 63)
        self.check_alignment_rows(
            alignments[2],
            [
                (
                    "gi|10955264|ref|NP_052605.1|",
                    "VYTSFN---GEKFSSYTLNKVTKTDEYNDLSELS...KGI",
                ),
                (
                    "gi|152973841|ref|YP_001338878.1|",
                    "VFGSFEQPKGEHLSGQVSEQ--RDTAFADQNEQV...QAM",
                ),
            ],
        )
        self.check_summary(alignments[0], "protein")
        self.check_reverse_write_read(alignments)

    def test_reading_alignments_fasta_m10_4(self):
        path = "Fasta/output004.m10"
        self.check_iterator_for_loop_handle(path, "fasta-m10", 1, 2)
        self.check_iterator_for_loop_filename(path, "fasta-m10", 1)
        self.check_iterator_next(path, "fasta-m10", 1)
        self.check_iterator_next_and_list(path, "fasta-m10", 1)
        self.check_iterator_next_for_loop(path, "fasta-m10", 1)
        alignment = self.check_read(path, "fasta-m10", 2, 102)
        self.check_alignment_rows(
            alignment,
            [
                (
                    "ref|NC_002127.1|:c1351-971",
                    "AAAAAAGATAAAAAATATCAAATAGAAGCAATAA...TCA",
                ),
                (
                    "ref|NC_002695.1|:1970775-1971404",
                    "AGAGAAAATAAAACAAGTAATAAAATATTAATGG...ACA",
                ),
            ],
        )
        self.check_summary(alignment, "DNA")

    def test_reading_alignments_fasta_m10_5(self):
        path = "Fasta/output005.m10"
        self.check_iterator_for_loop_handle(path, "fasta-m10", 1, 2)
        self.check_iterator_for_loop_filename(path, "fasta-m10", 1)
        self.check_iterator_next(path, "fasta-m10", 1)
        self.check_iterator_next_and_list(path, "fasta-m10", 1)
        self.check_iterator_next_for_loop(path, "fasta-m10", 1)
        alignment = self.check_read(path, "fasta-m10", 2, 110)
        self.check_alignment_rows(
            alignment,
            [
                (
                    "gi|10955264|ref|NP_052605.1|",
                    "IKNKDKTLFIVYAT-DIYSPSEFFSKIESDLKKK...LSK",
                ),
                (
                    "gi|10955282|ref|NP_052623.1|",
                    "IKDELPVAFCSWASLDLECEVKYINDVTSLYAKD...MSE",
                ),
            ],
        )
        self.check_summary(alignment, "protein")

    def test_reading_alignments_fasta_m10_6(self):
        path = "Fasta/output006.m10"
        self.check_iterator_for_loop_handle(path, "fasta-m10", 1, 2)
        self.check_iterator_for_loop_filename(path, "fasta-m10", 1)
        self.check_iterator_next(path, "fasta-m10", 1)
        self.check_iterator_next_and_list(path, "fasta-m10", 1)
        self.check_iterator_next_for_loop(path, "fasta-m10", 1)
        alignment = self.check_read(path, "fasta-m10", 2, 131)
        self.check_alignment_rows(
            alignment,
            [
                (
                    "gi|116660610|gb|EG558221.1|EG558221",
                    "GCAACGCTTCAAGAACTGGAATTAGGAACCGTGA...CAT",
                ),
                ("query", "GCAACGCTTCAAGAACTGGAATTAGGAACCGTGA...CAT"),
            ],
        )
        self.check_summary(alignment, "DNA")

    def test_reading_alignments_fasta_m10_7(self):
        path = "Fasta/output007.m10"
        alignments = self.check_iterator_for_loop_handle(path, "fasta-m10", 9, 2)
        self.check_iterator_for_loop_filename(path, "fasta-m10", 9)
        self.check_iterator_next(path, "fasta-m10", 9)
        self.check_iterator_next_and_list(path, "fasta-m10", 9)
        self.check_iterator_next_for_loop(path, "fasta-m10", 9)
        self.check_read_fails(path, "fasta-m10")
        self.assertEqual(alignments[0].get_alignment_length(), 108)
        self.check_alignment_rows(
            alignments[0],
            [
                (
                    "gi|10955263|ref|NP_052604.1|",
                    "SGSNT-RRRAISRPVRLTAEED---QEIRKRAAE...LSR",
                ),
                (
                    "gi|152973457|ref|YP_001338508.1|",
                    "AGSGAPRRRGSGLASRISEQSEALLQEAAKHAAE...LSR",
                ),
            ],
        )
        self.assertEqual(alignments[1].get_alignment_length(), 64)
        self.check_alignment_rows(
            alignments[1],
            [
                (
                    "gi|10955263|ref|NP_052604.1|",
                    "AAECGKTVSGFLRAAALGKKVNSLTDDRVLKEV-...AIT",
                ),
                (
                    "gi|152973588|ref|YP_001338639.1|",
                    "ASRQGCTVGG--KMDSVQDKASDKDKERVMKNIN...TLT",
                ),
            ],
        )
        self.assertEqual(alignments[2].get_alignment_length(), 45)
        self.check_alignment_rows(
            alignments[2],
            [
                (
                    "gi|10955263|ref|NP_052604.1|",
                    "EIRKRAAECGKTVSGFLRAAA-----LGKKVNSL...VMR",
                ),
                (
                    "gi|152973480|ref|YP_001338531.1|",
                    "ELVKLIADMGISVRALLRKNVEPYEELGLEEDKF...MLQ",
                ),
            ],
        )
        self.assertEqual(alignments[8].get_alignment_length(), 64)
        self.check_alignment_rows(
            alignments[8],
            [
                (
                    "gi|10955265|ref|NP_052606.1|",
                    "ISGTYKGIDFLIKLMPSGGNTTIGRASGQNNTYF...FSD",
                ),
                (
                    "gi|152973505|ref|YP_001338556.1|",
                    "IDGVITAFD-LRTGMNISKDKVVAQIQGMDPVW-...YPD",
                ),
            ],
        )
        self.check_summary(alignments[0], "protein")
        self.check_reverse_write_read(alignments)

    def test_reading_alignments_fasta_m10_8(self):
        path = "Fasta/output008.m10"
        alignments = self.check_iterator_for_loop_handle(path, "fasta-m10", 12, 2)
        self.check_iterator_for_loop_filename(path, "fasta-m10", 12)
        self.check_iterator_next(path, "fasta-m10", 12)
        self.check_iterator_next_and_list(path, "fasta-m10", 12)
        self.check_iterator_next_for_loop(path, "fasta-m10", 12)
        self.check_read_fails(path, "fasta-m10")
        self.assertEqual(alignments[0].get_alignment_length(), 65)
        self.check_alignment_rows(
            alignments[0],
            [
                (
                    "gi|283855822|gb|GQ290312.1|",
                    "IPHQLPHALRHRPAQEAAHASQLHPAQPGCGQPL...GLL",
                ),
                ("sp|Q9NSY1|BMP2K_HUMAN", "LQHRHPHQQQQQQQQQQQQQQQQQQQQQQQQQQQ...QML"),
            ],
        )
        self.assertEqual(alignments[1].get_alignment_length(), 201)
        self.check_alignment_rows(
            alignments[1],
            [
                (
                    "gi|57163782|ref|NM_001009242.1|",
                    "GPELLRALLQQNGCGTQPLRVPTVLPG*AMAVLH...QKS",
                ),
                ("sp|Q9NSY1|BMP2K_HUMAN", "GPEIL---LGQ-GPPQQPPQQHRVLQQLQQGDWR...NRS"),
            ],
        )
        self.assertEqual(alignments[2].get_alignment_length(), 348)
        self.check_alignment_rows(
            alignments[2],
            [
                (
                    "gi|57163782|ref|NM_001009242.1|",
                    "MNGTEGPNFYVPFSNKTGVVRSPFEYPQYYLAEP...APA",
                ),
                ("sp|P08100|OPSD_HUMAN", "MNGTEGPNFYVPFSNATGVVRSPFEYPQYYLAEP...APA"),
            ],
        )
        self.assertEqual(alignments[11].get_alignment_length(), 31)
        self.check_alignment_rows(
            alignments[11],
            [
                ("gi|283855822|gb|GQ290312.1|", "SQQIRNATTMMMTMRVTSFSAFWVVADSCCW"),
                ("sp|P08100|OPSD_HUMAN", "AQQQESATTQKAEKEVTRMVIIMVIAFLICW"),
            ],
        )
        self.check_summary(alignments[0], "protein")
        self.check_reverse_write_read(alignments)

    def test_reading_alignments_ig(self):
        path = "IntelliGenetics/VIF_mase-pro.txt"
        self.check_iterator_for_loop_handle(path, "ig", 1, 16)
        self.check_iterator_for_loop_filename(path, "ig", 1)
        self.check_iterator_next(path, "ig", 1)
        self.check_iterator_next_and_list(path, "ig", 1)
        self.check_iterator_next_for_loop(path, "ig", 1)
        self.check_write_three_times_and_read(path, "ig", 16)
        alignment = self.check_read(path, "ig", 16, 298)
        self.check_alignment_columns(
            alignment,
            [
                "MMMMMMMMMMMMMMMM",
                "EEEEEEETEEEENEEE",
                "NNNNNNNAEEEEQRKK",
                "--------DEEEEE--",
                "--------KKKKKK--",
                "HHHHHHH-AAAAL-R-",
            ],
        )
        self.check_summary(alignment, "protein")

    def test_reading_alignments_pir(self):
        path = "NBRF/clustalw.pir"
        self.check_iterator_for_loop_handle(path, "pir", 1, 2)
        self.check_iterator_for_loop_filename(path, "pir", 1)
        self.check_iterator_next(path, "pir", 1)
        self.check_iterator_next_and_list(path, "pir", 1)
        self.check_iterator_next_for_loop(path, "pir", 1)
        self.check_write_three_times_and_read(path, "pir", 2)
        alignment = self.check_read(path, "pir", 2, 2527)
        self.check_alignment_rows(
            alignment,
            [
                (
                    "804Angiostrongylus_cantonensis",
                    "----------------------------------...---",
                ),
                (
                    "815Parelaphostrongylus_odocoil",
                    "----------------------------------...---",
                ),
            ],
        )
        self.check_summary_pir(alignment)

    def test_reading_alignments_maf1(self):
        path = "MAF/humor.maf"
        alignments = self.check_iterator_for_loop_handle(path, "maf", 2, 3)
        self.check_iterator_for_loop_filename(path, "maf", 2)
        self.check_iterator_next(path, "maf", 2)
        self.check_iterator_next_and_list(path, "maf", 2)
        self.check_iterator_next_for_loop(path, "maf", 2)
        self.check_read_fails(path, "maf")
        self.assertEqual(alignments[0].get_alignment_length(), 5486)
        self.check_alignment_rows(
            alignments[0],
            [
                ("NM_006987", "gcacagcctttactccctgactgcgtttatattc...CCG"),
                ("mm3", "gcacagcctttactccctgactgcgtttatattc...TTG"),
                ("rn3", "gcacagcctttactccctgactgcgtttatattc...CCG"),
            ],
        )
        self.assertEqual(alignments[1].get_alignment_length(), 5753)
        self.check_alignment_rows(
            alignments[1],
            [
                ("NM_018289", "tttgtccatgttggtcaggctggtctcgaactcc...GGT"),
                ("mm3", "tttgtccatgttggtcaggctggtctcgaactcc...GGT"),
                ("rn3", "tttgtccatgttggtcaggctggtctcgaactcc...GGT"),
            ],
        )
        self.check_summary(alignments[1], "DNA")
        self.check_reverse_write_read(alignments)

    def test_reading_alignments_maf2(self):
        path = "MAF/bug2453.maf"
        alignments = self.check_iterator_for_loop_handle(path, "maf", 3)
        self.check_iterator_for_loop_filename(path, "maf", 3)
        self.check_iterator_next(path, "maf", 3)
        self.check_iterator_next_and_list(path, "maf", 3)
        self.check_iterator_next_for_loop(path, "maf", 3)
        self.check_read_fails(path, "maf")
        self.assertEqual(len(alignments[0]), 5)
        self.assertEqual(alignments[0].get_alignment_length(), 42)
        self.check_alignment_columns(
            alignments[0], ["AAA--", "AAAAA", "AAAAA", "---T-", "GGGGG", "GGGGG"]
        )
        self.assertEqual(len(alignments[1]), 5)
        self.assertEqual(alignments[1].get_alignment_length(), 6)
        self.check_alignment_columns(
            alignments[1], ["TTTTt", "AAAAa", "AAAAa", "AAAAg", "GGGGg", "AAAAa"]
        )
        self.assertEqual(len(alignments[2]), 4)
        self.assertEqual(alignments[2].get_alignment_length(), 13)
        self.check_alignment_rows(
            alignments[2],
            [
                ("baboon", "gcagctgaaaaca"),
                ("hg16.chr7", "gcagctgaaaaca"),
                ("mm4.chr6", "ACAGCTGAAAATA"),
                ("panTro1.chr6", "gcagctgaaaaca"),
            ],
        )
        self.check_summary(alignments[1], "DNA")
        self.check_reverse_write_read(alignments)

    def test_reading_alignments_maf3(self):
        path = "MAF/ucsc_test.maf"
        alignments = self.check_iterator_for_loop_handle(path, "maf", 3)
        self.check_iterator_for_loop_filename(path, "maf", 3)
        self.check_iterator_next(path, "maf", 3)
        self.check_iterator_next_and_list(path, "maf", 3)
        self.check_iterator_next_for_loop(path, "maf", 3)
        self.check_read_fails(path, "maf")
        self.assertEqual(len(alignments[0]), 5)
        self.assertEqual(alignments[0].get_alignment_length(), 42)
        self.check_alignment_columns(
            alignments[0], ["AAA--", "AAAAA", "AAAAA", "---T-", "GGGGG", "GGGGG"]
        )
        self.assertEqual(len(alignments[1]), 5)
        self.assertEqual(alignments[1].get_alignment_length(), 6)
        self.check_alignment_columns(
            alignments[1], ["TTTTt", "AAAAa", "AAAAa", "AAAAg", "GGGGg", "AAAAa"]
        )
        self.assertEqual(len(alignments[2]), 4)
        self.assertEqual(alignments[2].get_alignment_length(), 13)
        self.check_alignment_columns(
            alignments[2], ["gggA", "cccC", "aaaA", "gggG", "cccC", "aaaA"]
        )
        self.check_summary(alignments[2], "DNA")
        self.check_reverse_write_read(alignments)

    def test_reading_alignments_maf4(self):
        path = "MAF/ucsc_mm9_chr10.maf"
        alignments = self.check_iterator_for_loop_handle(path, "maf", 48)
        self.check_iterator_for_loop_filename(path, "maf", 48)
        self.check_iterator_next(path, "maf", 48)
        self.check_iterator_next_and_list(path, "maf", 48)
        self.check_iterator_next_for_loop(path, "maf", 48)
        self.check_read_fails(path, "maf")
        self.assertEqual(len(alignments[0]), 2)
        self.assertEqual(alignments[0].get_alignment_length(), 164)
        self.check_alignment_rows(
            alignments[0],
            [
                ("mm9.chr10", "TCATAGGTATTTATTTTTAAATATGGTTTGCTTT...GTT"),
                ("oryCun1.scaffold_133159", "TCACAGATATTTACTATTAAATATGGTTTGTTAT...GTT"),
            ],
        )
        self.assertEqual(len(alignments[1]), 4)
        self.assertEqual(alignments[1].get_alignment_length(), 466)
        self.check_alignment_rows(
            alignments[1],
            [
                ("hg18.chr6", "AGTCTTCATAAGTGGAAATATAAGTTTTAATTAT...TTC"),
                ("mm9.chr10", "AGTCTTTCCAATGGGACCTGTGAGTCCTAACTAT...CTG"),
                ("panTro2.chr6", "AGTCTTCATAAGTGGAAATATAAGTTTTAATTAT...TTC"),
                ("ponAbe2.chr6", "AGTCTTCATAAGTGGAAATATAAGTTTTAATTAT...TTC"),
            ],
        )
        self.assertEqual(len(alignments[2]), 5)
        self.assertEqual(alignments[2].get_alignment_length(), 127)
        self.check_alignment_columns(
            alignments[2], ["TTTTT", "GGGGG", "GGGGG", "GGGGG", "TTTTC", "CCCCC"]
        )
        self.assertEqual(len(alignments[47]), 6)
        self.assertEqual(alignments[47].get_alignment_length(), 46)
        self.check_alignment_columns(
            alignments[47], ["TTTTTT", "GGGGGG", "TTTTTT", "TTTTTT", "TGGGAT", "tTTTT-"]
        )
        self.check_summary(alignments[47], "DNA")
        self.check_reverse_write_read(alignments)

    def test_reading_alignments_mauve(self):
        path = "Mauve/simple.xmfa"
        alignments = self.check_iterator_for_loop_handle(path, "mauve", 5)
        self.check_iterator_for_loop_filename(path, "mauve", 5)
        self.check_iterator_next(path, "mauve", 5)
        self.check_iterator_next_and_list(path, "mauve", 5)
        self.check_iterator_next_for_loop(path, "mauve", 5)
        self.check_read_fails(path, "mauve")
        self.assertEqual(len(alignments[0]), 2)
        self.assertEqual(alignments[0].get_alignment_length(), 5670)
        self.check_alignment_rows(
            alignments[0],
            [
                ("1/0-5670", "ATATTAGGTTTTTACCTACCCAGGAAAAGCCAAC...AAT"),
                ("2/0-5670", "ATATTAGGTTTTTACCTACCCAGGAAAAGCCAAC...AAT"),
            ],
        )
        self.assertEqual(len(alignments[1]), 2)
        self.assertEqual(alignments[1].get_alignment_length(), 4420)
        self.check_alignment_rows(
            alignments[1],
            [
                ("1/5670-9940", "GAACATCAGCACCTGAGTTGCTAAAGTCATTTAG...CTC"),
                ("2/7140-11410", "GAACATCAGCACCTGAGTTGCTAAAGTCATTTAG...CTC"),
            ],
        )
        self.assertEqual(len(alignments[2]), 1)
        self.assertEqual(alignments[2].get_alignment_length(), 4970)
        self.check_alignment_rows(
            alignments[2],
            [("1/9940-14910", "TCTACCAACCACCACAGACATCAATCACTTCTGC...GAC")],
        )
        self.assertEqual(len(alignments[3]), 1)
        self.assertEqual(alignments[3].get_alignment_length(), 1470)
        self.assertEqual(len(alignments[4]), 1)
        self.assertEqual(alignments[4].get_alignment_length(), 1470)
        self.check_alignment_rows(
            alignments[4],
            [("2/11410-12880", "ATTCGCACATAAGAATGTACCTTGCTGTAATTTA...ATA")],
        )
        self.check_summary(alignments[4], "DNA")
        self.check_reverse_write_read(alignments)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
