# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for seq module."""

import array
import copy
import unittest
import warnings

from Bio import BiopythonWarning
from Bio import Seq
from Bio.Data.IUPACData import (
    ambiguous_dna_complement,
    ambiguous_rna_complement,
    ambiguous_dna_values,
    ambiguous_rna_values,
)
from Bio.Data.CodonTable import TranslationError, standard_dna_table
from Bio.Seq import MutableSeq

test_seqs = [
    Seq.Seq("TCAAAAGGATGCATCATG"),
    Seq.Seq("T"),
    Seq.Seq("ATGAAACTG"),
    Seq.Seq("ATGAARCTG"),
    Seq.Seq("AWGAARCKG"),  # Note no U or T
    Seq.Seq("".join(ambiguous_rna_values)),
    Seq.Seq("".join(ambiguous_dna_values)),
    Seq.Seq("AWGAARCKG"),
    Seq.Seq("AUGAAACUG"),
    Seq.Seq("ATGAAA-CTG"),
    Seq.Seq("ATGAAACTGWN"),
    Seq.Seq("AUGAAA==CUG"),
    Seq.Seq("AUGAAACUGWN"),
    Seq.Seq("AUGAAACTG"),  # U and T
    Seq.MutableSeq("ATGAAACTG"),
    Seq.MutableSeq("AUGaaaCUG"),
    Seq.Seq("ACTGTCGTCT"),
]
protein_seqs = [
    Seq.Seq("ATCGPK"),
    Seq.Seq("T.CGPK"),
    Seq.Seq("T-CGPK"),
    Seq.Seq("MEDG-KRXR*"),
    Seq.MutableSeq("ME-K-DRXR*XU"),
    Seq.Seq("MEDG-KRXR@"),
    Seq.Seq("ME-KR@"),
    Seq.Seq("MEDG.KRXR@"),
]


class TestSeq(unittest.TestCase):
    def setUp(self):
        self.s = Seq.Seq("TCAAAAGGATGCATCATG")

    def test_as_string(self):
        """Test converting Seq to string."""
        self.assertEqual("TCAAAAGGATGCATCATG", str(self.s))

    def test_construction_using_a_seq_object(self):
        """Test using a Seq object to initialize another Seq object."""
        with self.assertRaises(TypeError):
            Seq.Seq(self.s)

    def test_repr(self):
        """Test representation of Seq object."""
        self.assertEqual("Seq('TCAAAAGGATGCATCATG')", repr(self.s))

    def test_truncated_repr(self):
        seq = "TCAAAAGGATGCATCATGTCAAAAGGATGCATCATGTCAAAAGGATGCATCATGTCAAAAGGA"
        expected = "Seq('TCAAAAGGATGCATCATGTCAAAAGGATGCATCATGTCAAAAGGATGCATCATG...GGA')"
        self.assertEqual(expected, repr(Seq.Seq(seq)))

    def test_length(self):
        """Test len method on Seq object."""
        self.assertEqual(18, len(self.s))

    def test_first_nucleotide(self):
        """Test getting first nucleotide of Seq."""
        self.assertEqual("T", self.s[0])

    def test_last_nucleotide(self):
        """Test getting last nucleotide of Seq."""
        self.assertEqual("G", self.s[-1])

    def test_slicing(self):
        """Test slicing of Seq."""
        self.assertEqual("AA", str(self.s[3:5]))

    def test_reverse(self):
        """Test reverse using -1 stride."""
        self.assertEqual("GTACTACGTAGGAAAACT", self.s[::-1])

    def test_extract_third_nucleotide(self):
        """Test extracting every third nucleotide (slicing with stride 3)."""
        self.assertEqual("TAGTAA", str(self.s[0::3]))
        self.assertEqual("CAGGTT", str(self.s[1::3]))
        self.assertEqual("AAACCG", str(self.s[2::3]))

    def test_concatenation_of_seq(self):
        t = Seq.Seq("T")
        u = self.s + t
        self.assertEqual(str(self.s) + "T", str(u))
        self.assertEqual(self.s + Seq.Seq("T"), "TCAAAAGGATGCATCATGT")

    def test_ungap(self):
        self.assertEqual("ATCCCA", str(Seq.Seq("ATC-CCA").ungap("-")))

        with self.assertRaises(ValueError):
            Seq.Seq("ATC-CCA").ungap("--")

        with self.assertRaises(ValueError):
            Seq.Seq("ATC-CCA").ungap(gap=None)


class TestSeqStringMethods(unittest.TestCase):
    def setUp(self):
        self.s = Seq.Seq("TCAAAAGGATGCATCATG")
        self.dna = [
            Seq.Seq("ATCG"),
            Seq.Seq("gtca"),
            Seq.MutableSeq("GGTCA"),
            Seq.Seq("CTG-CA"),
        ]
        self.rna = [
            Seq.Seq("AUUUCG"),
            Seq.MutableSeq("AUUCG"),
            Seq.Seq("uCAg"),
            Seq.MutableSeq("UC-AG"),
            Seq.Seq("U.CAG"),
        ]
        self.nuc = [Seq.Seq("ATCG")]
        self.protein = [
            Seq.Seq("ATCGPK"),
            Seq.Seq("atcGPK"),
            Seq.Seq("T.CGPK"),
            Seq.Seq("T-CGPK"),
            Seq.Seq("MEDG-KRXR*"),
            Seq.MutableSeq("ME-K-DRXR*XU"),
            Seq.Seq("MEDG-KRXR@"),
            Seq.Seq("ME-KR@"),
            Seq.Seq("MEDG.KRXR@"),
        ]
        self.test_chars = ["-", Seq.Seq("-"), Seq.Seq("*"), "-X@"]

    def test_string_methods(self):
        for a in self.dna + self.rna + self.nuc + self.protein:
            if isinstance(a, Seq.Seq):
                self.assertEqual(str(a.strip()), str(a).strip())
                self.assertEqual(str(a.lstrip()), str(a).lstrip())
                self.assertEqual(str(a.rstrip()), str(a).rstrip())
                self.assertEqual(str(a.lower()), str(a).lower())
                self.assertEqual(str(a.upper()), str(a).upper())

    def test_hash(self):
        with warnings.catch_warnings(record=True):
            hash(self.s)

    def test_not_equal_comparsion(self):
        """Test __ne__ comparison method."""
        self.assertNotEqual(Seq.Seq("TCAAA"), Seq.Seq("TCAAAA"))

    def test_less_than_comparison(self):
        """Test __lt__ comparison method."""
        self.assertLess(self.s[:-1], self.s)

    def test_less_than_comparison_of_incompatible_types(self):
        """Test incompatible types __lt__ comparison method."""
        with self.assertRaises(TypeError):
            self.s < 1

    def test_less_than_or_equal_comparison(self):
        """Test __le__ comparison method."""
        self.assertLessEqual(self.s, self.s)

    def test_less_than_or_equal_comparison_of_incompatible_types(self):
        """Test incompatible types __le__ comparison method."""
        with self.assertRaises(TypeError):
            self.s <= 1

    def test_greater_than_comparison(self):
        """Test __gt__ comparison method."""
        self.assertGreater(self.s, self.s[:-1])

    def test_greater_than_comparison_of_incompatible_types(self):
        """Test incompatible types __gt__ comparison method."""
        with self.assertRaises(TypeError):
            self.s > 1

    def test_greater_than_or_equal_comparison(self):
        """Test __ge__ comparison method."""
        self.assertGreaterEqual(self.s, self.s)

    def test_greater_than_or_equal_comparison_of_incompatible_types(self):
        """Test incompatible types __ge__ comparison method."""
        with self.assertRaises(TypeError):
            self.s >= 1

    def test_add_method_using_wrong_object(self):
        with self.assertRaises(TypeError):
            self.s + {}

    def test_radd_method(self):
        self.assertEqual(
            "TCAAAAGGATGCATCATGTCAAAAGGATGCATCATG", str(self.s.__radd__(self.s))
        )

    def test_radd_method_using_wrong_object(self):
        with self.assertRaises(TypeError):
            self.s.__radd__({})

    def test_contains_method(self):
        self.assertIn("AAAA", self.s)

    def test_startswith(self):
        self.assertTrue(self.s.startswith("TCA"))
        self.assertTrue(self.s.startswith(("CAA", "CTA"), 1))

    def test_endswith(self):
        self.assertTrue(self.s.endswith("ATG"))
        self.assertTrue(self.s.endswith(("ATG", "CTA")))

    def test_append_nucleotides(self):
        self.test_chars.append(Seq.Seq("A"))
        self.assertEqual(5, len(self.test_chars))

    def test_append_proteins(self):
        self.test_chars.append(Seq.Seq("K"))
        self.test_chars.append(Seq.Seq("K-"))
        self.test_chars.append(Seq.Seq("K@"))

        self.assertEqual(7, len(self.test_chars))

    def test_stripping_characters(self):
        for a in self.dna + self.rna + self.nuc + self.protein:
            for char in self.test_chars:
                str_char = str(char)
                if isinstance(a, Seq.Seq):
                    self.assertEqual(str(a.strip(char)), str(a).strip(str_char))
                    self.assertEqual(str(a.lstrip(char)), str(a).lstrip(str_char))
                    self.assertEqual(str(a.rstrip(char)), str(a).rstrip(str_char))

    def test_finding_characters(self):
        for a in self.dna + self.rna + self.nuc + self.protein:
            for char in self.test_chars:
                str_char = str(char)
                if isinstance(a, Seq.Seq):
                    self.assertEqual(a.find(char), str(a).find(str_char))
                    self.assertEqual(a.find(char, 2, -2), str(a).find(str_char, 2, -2))
                    self.assertEqual(a.rfind(char), str(a).rfind(str_char))
                    self.assertEqual(
                        a.rfind(char, 2, -2), str(a).rfind(str_char, 2, -2)
                    )

    def test_counting_characters(self):
        for a in self.dna + self.rna + self.nuc + self.protein:
            for char in self.test_chars:
                str_char = str(char)
                if isinstance(a, Seq.Seq):
                    self.assertEqual(a.count(char), str(a).count(str_char))
                    self.assertEqual(
                        a.count(char, 2, -2), str(a).count(str_char, 2, -2)
                    )

    def test_splits(self):
        for a in self.dna + self.rna + self.nuc + self.protein:
            for char in self.test_chars:
                str_char = str(char)
                if isinstance(a, Seq.Seq):
                    self.assertEqual(
                        [str(x) for x in a.split(char)], str(a).split(str_char)
                    )
                    self.assertEqual(
                        [str(x) for x in a.rsplit(char)], str(a).rsplit(str_char)
                    )

                    for max_sep in [0, 1, 2, 999]:
                        self.assertEqual(
                            [str(x) for x in a.split(char, max_sep)],
                            str(a).split(str_char, max_sep),
                        )


class TestSeqAddition(unittest.TestCase):
    def setUp(self):
        self.dna = [
            Seq.Seq("ATCG"),
            Seq.Seq("gtca"),
            Seq.MutableSeq("GGTCA"),
            Seq.Seq("CTG-CA"),
            "TGGTCA",
        ]
        self.rna = [
            Seq.Seq("AUUUCG"),
            Seq.MutableSeq("AUUCG"),
            Seq.Seq("uCAg"),
            Seq.MutableSeq("UC-AG"),
            Seq.Seq("U.CAG"),
            "UGCAU",
        ]
        self.nuc = [
            Seq.Seq("ATCG"),
            "UUUTTTACG",
        ]
        self.protein = [
            Seq.Seq("ATCGPK"),
            Seq.Seq("atcGPK"),
            Seq.Seq("T.CGPK"),
            Seq.Seq("T-CGPK"),
            Seq.Seq("MEDG-KRXR*"),
            Seq.MutableSeq("ME-K-DRXR*XU"),
            "TEDDF",
        ]

    def test_addition_dna_rna_with_generic_nucleotides(self):
        for a in self.dna + self.rna:
            for b in self.nuc:
                c = a + b
                self.assertEqual(str(c), str(a) + str(b))

    def test_addition_dna_rna_with_generic_nucleotides_inplace(self):
        for a in self.dna + self.rna:
            for b in self.nuc:
                c = b + a
                b += a  # can't change 'a' as need value next iteration
                self.assertEqual(c, b)

    def test_addition_rna_with_rna(self):
        self.rna.pop(3)
        for a in self.rna:
            for b in self.rna:
                c = a + b
                self.assertEqual(str(c), str(a) + str(b))

    def test_addition_rna_with_rna_inplace(self):
        self.rna.pop(3)
        for a in self.rna:
            for b in self.rna:
                c = b + a
                b += a
                self.assertEqual(c, b)

    def test_addition_dna_with_dna(self):
        for a in self.dna:
            for b in self.dna:
                c = a + b
                self.assertEqual(str(c), str(a) + str(b))

    def test_addition_dna_with_dna_inplace(self):
        for a in self.dna:
            for b in self.dna:
                c = b + a
                b += a
                self.assertEqual(c, b)

    def test_addition_dna_with_rna(self):
        self.dna.pop(4)
        self.rna.pop(5)
        for a in self.dna:
            for b in self.rna:
                self.assertEqual(str(a) + str(b), a + b)
                self.assertEqual(str(b) + str(a), b + a)
                # Check in place works
                c = a
                c += b
                self.assertEqual(c, str(a) + str(b))
                c = b
                c += a
                self.assertEqual(c, str(b) + str(a))

    def test_addition_proteins(self):
        self.protein.pop(2)
        for a in self.protein:
            for b in self.protein:
                c = a + b
                self.assertEqual(str(c), str(a) + str(b))

    def test_addition_proteins_inplace(self):
        self.protein.pop(2)
        for a in self.protein:
            for b in self.protein:
                c = b + a
                b += a
                self.assertEqual(c, b)

    def test_adding_protein_with_nucleotides(self):
        for a in self.protein[0:5]:
            for b in self.dna[0:3] + self.rna[0:4]:
                self.assertEqual(str(a) + str(b), a + b)
                a += b

    def test_adding_generic_nucleotide_with_other_nucleotides(self):
        for a in self.nuc:
            for b in self.dna + self.rna + self.nuc:
                c = a + b
                self.assertEqual(str(c), str(a) + str(b))

    def test_adding_generic_nucleotide_with_other_nucleotides_inplace(self):
        for a in self.nuc:
            for b in self.dna + self.rna + self.nuc:
                c = b + a
                b += a
                self.assertEqual(c, b)


class TestSeqMultiplication(unittest.TestCase):
    def test_mul_method(self):
        """Test mul method; relies on addition method."""
        for seq in test_seqs + protein_seqs:
            self.assertEqual(seq * 3, seq + seq + seq)

    def test_mul_method_exceptions(self):
        """Test mul method exceptions."""
        for seq in test_seqs + protein_seqs:
            with self.assertRaises(TypeError):
                seq * 3.0
            with self.assertRaises(TypeError):
                seq * ""

    def test_rmul_method(self):
        """Test rmul method; relies on addition method."""
        for seq in test_seqs + protein_seqs:
            self.assertEqual(3 * seq, seq + seq + seq)

    def test_rmul_method_exceptions(self):
        """Test rmul method exceptions."""
        for seq in test_seqs + protein_seqs:
            with self.assertRaises(TypeError):
                3.0 * seq
            with self.assertRaises(TypeError):
                "" * seq

    def test_imul_method(self):
        """Test imul method; relies on addition and mull methods."""
        for seq in test_seqs + protein_seqs:
            original_seq = seq * 1  # make a copy
            seq *= 3
            self.assertEqual(seq, original_seq + original_seq + original_seq)

    def test_imul_method_exceptions(self):
        """Test imul method exceptions."""
        for seq in test_seqs + protein_seqs:
            with self.assertRaises(TypeError):
                seq *= 3.0
            with self.assertRaises(TypeError):
                seq *= ""


class TestMutableSeq(unittest.TestCase):
    def setUp(self):
        self.s = Seq.Seq("TCAAAAGGATGCATCATG")
        self.mutable_s = MutableSeq("TCAAAAGGATGCATCATG")

    def test_mutableseq_creation(self):
        """Test creating MutableSeqs in multiple ways."""
        mutable_s = MutableSeq("TCAAAAGGATGCATCATG")
        self.assertIsInstance(mutable_s, MutableSeq, "Creating MutableSeq")

        mutable_s = self.s.tomutable()
        self.assertIsInstance(mutable_s, MutableSeq, "Converting Seq to mutable")

        array_seq = MutableSeq(array.array("u", "TCAAAAGGATGCATCATG"))
        self.assertIsInstance(array_seq, MutableSeq, "Creating MutableSeq using array")

    def test_repr(self):
        self.assertEqual("MutableSeq('TCAAAAGGATGCATCATG')", repr(self.mutable_s))

    def test_truncated_repr(self):
        seq = "TCAAAAGGATGCATCATGTCAAAAGGATGCATCATGTCAAAAGGATGCATCATGTCAAAAGGA"
        expected = (
            "MutableSeq('TCAAAAGGATGCATCATGTCAAAAGGATGCATCATGTCAAAAGGATGCATCATG...GGA')"
        )
        self.assertEqual(expected, repr(MutableSeq(seq)))

    def test_equal_comparison(self):
        """Test __eq__ comparison method."""
        self.assertEqual(self.mutable_s, "TCAAAAGGATGCATCATG")

    def test_not_equal_comparison(self):
        """Test __ne__ comparison method."""
        self.assertNotEqual(self.mutable_s, "other thing")

    def test_less_than_comparison(self):
        """Test __lt__ comparison method."""
        self.assertLess(self.mutable_s[:-1], self.mutable_s)

    def test_less_than_comparison_of_incompatible_types(self):
        with self.assertRaises(TypeError):
            self.mutable_s < 1

    def test_less_than_comparison_with_str(self):
        self.assertLessEqual(self.mutable_s[:-1], "TCAAAAGGATGCATCATG")

    def test_less_than_or_equal_comparison(self):
        """Test __le__ comparison method."""
        self.assertLessEqual(self.mutable_s[:-1], self.mutable_s)

    def test_less_than_or_equal_comparison_of_incompatible_types(self):
        with self.assertRaises(TypeError):
            self.mutable_s <= 1

    def test_less_than_or_equal_comparison_with_str(self):
        self.assertLessEqual(self.mutable_s[:-1], "TCAAAAGGATGCATCATG")

    def test_greater_than_comparison(self):
        """Test __gt__ comparison method."""
        self.assertGreater(self.mutable_s, self.mutable_s[:-1])

    def test_greater_than_comparison_of_incompatible_types(self):
        with self.assertRaises(TypeError):
            self.mutable_s > 1

    def test_greater_than_comparison_with_str(self):
        self.assertGreater(self.mutable_s, "TCAAAAGGATGCATCAT")

    def test_greater_than_or_equal_comparison(self):
        """Test __ge__ comparison method."""
        self.assertGreaterEqual(self.mutable_s, self.mutable_s)

    def test_greater_than_or_equal_comparison_of_incompatible_types(self):
        with self.assertRaises(TypeError):
            self.mutable_s >= 1

    def test_greater_than_or_equal_comparison_with_str(self):
        self.assertGreaterEqual(self.mutable_s, "TCAAAAGGATGCATCATG")

    def test_add_method(self):
        """Test adding wrong type to MutableSeq."""
        with self.assertRaises(TypeError):
            self.mutable_s + 1234

    def test_radd_method(self):
        self.assertEqual(
            "TCAAAAGGATGCATCATGTCAAAAGGATGCATCATG",
            self.mutable_s.__radd__(self.mutable_s),
        )

    def test_radd_method_using_mutalbeseq_object(self):
        self.assertEqual(
            "UCAAAAGGATCAAAAGGATGCATCATG",
            self.mutable_s.__radd__(MutableSeq("UCAAAAGGA")),
        )

    def test_radd_method_using_seq_object(self):
        self.assertEqual(
            "TCAAAAGGATGCATCATGTCAAAAGGATGCATCATG", self.mutable_s.__radd__(self.s)
        )

    def test_radd_method_wrong_type(self):
        with self.assertRaises(TypeError):
            self.mutable_s.__radd__(1234)

    def test_as_string(self):
        self.assertEqual("TCAAAAGGATGCATCATG", str(self.mutable_s))

    def test_length(self):
        self.assertEqual(18, len(self.mutable_s))

    def test_converting_to_immutable(self):
        self.assertIsInstance(self.mutable_s.toseq(), Seq.Seq)

    def test_first_nucleotide(self):
        self.assertEqual("T", self.mutable_s[0])

    def test_setting_slices(self):
        self.assertEqual(
            MutableSeq("CAAA"), self.mutable_s[1:5], "Slice mutable seq",
        )

        self.mutable_s[1:3] = "GAT"
        self.assertEqual(
            MutableSeq("TGATAAAGGATGCATCATG"),
            self.mutable_s,
            "Set slice with string and adding extra nucleotide",
        )

        self.mutable_s[1:3] = self.mutable_s[5:7]
        self.assertEqual(
            MutableSeq("TAATAAAGGATGCATCATG"),
            self.mutable_s,
            "Set slice with MutableSeq",
        )

        self.mutable_s[1:3] = array.array("u", "GAT")
        self.assertEqual(
            MutableSeq("TGATTAAAGGATGCATCATG"), self.mutable_s, "Set slice with array",
        )

    def test_setting_item(self):
        self.mutable_s[3] = "G"
        self.assertEqual(MutableSeq("TCAGAAGGATGCATCATG"), self.mutable_s)

    def test_deleting_slice(self):
        del self.mutable_s[4:5]
        self.assertEqual(MutableSeq("TCAAAGGATGCATCATG"), self.mutable_s)

    def test_deleting_item(self):
        del self.mutable_s[3]
        self.assertEqual(MutableSeq("TCAAAGGATGCATCATG"), self.mutable_s)

    def test_appending(self):
        self.mutable_s.append("C")
        self.assertEqual(MutableSeq("TCAAAAGGATGCATCATGC"), self.mutable_s)

    def test_inserting(self):
        self.mutable_s.insert(4, "G")
        self.assertEqual(MutableSeq("TCAAGAAGGATGCATCATG"), self.mutable_s)

    def test_popping_last_item(self):
        self.assertEqual("G", self.mutable_s.pop())

    def test_remove_items(self):
        self.mutable_s.remove("G")
        self.assertEqual(
            MutableSeq("TCAAAAGATGCATCATG"), self.mutable_s, "Remove first G"
        )

        self.assertRaises(ValueError, self.mutable_s.remove, "Z")

    def test_count(self):
        self.assertEqual(7, self.mutable_s.count("A"))
        self.assertEqual(2, self.mutable_s.count("AA"))

    def test_index(self):
        self.assertEqual(2, self.mutable_s.index("A"))
        self.assertRaises(ValueError, self.mutable_s.index, "8888")

    def test_reverse(self):
        """Test using reverse method."""
        self.mutable_s.reverse()
        self.assertEqual(MutableSeq("GTACTACGTAGGAAAACT"), self.mutable_s)

    def test_reverse_with_stride(self):
        """Test reverse using -1 stride."""
        self.assertEqual(MutableSeq("GTACTACGTAGGAAAACT"), self.mutable_s[::-1])

    def test_complement(self):
        self.mutable_s.complement()
        self.assertEqual("AGTTTTCCTACGTAGTAC", str(self.mutable_s))

    def test_complement_rna(self):
        seq = Seq.MutableSeq("AUGaaaCUG")
        seq.complement()
        self.assertEqual("UACuuuGAC", str(seq))

    def test_complement_mixed_aphabets(self):
        seq = Seq.MutableSeq("AUGaaaCTG")
        with self.assertRaises(ValueError):
            seq.complement()

    def test_complement_rna_string(self):
        seq = Seq.MutableSeq("AUGaaaCUG")
        seq.complement()
        self.assertEqual("UACuuuGAC", str(seq))

    def test_complement_dna_string(self):
        seq = Seq.MutableSeq("ATGaaaCTG")
        seq.complement()
        self.assertEqual("TACtttGAC", str(seq))

    def test_reverse_complement(self):
        self.mutable_s.reverse_complement()
        self.assertEqual("CATGATGCATCCTTTTGA", str(self.mutable_s))

    def test_extend_method(self):
        self.mutable_s.extend("GAT")
        self.assertEqual(MutableSeq("TCAAAAGGATGCATCATGGAT"), self.mutable_s)

    def test_extend_with_mutable_seq(self):
        self.mutable_s.extend(MutableSeq("TTT"))
        self.assertEqual(MutableSeq("TCAAAAGGATGCATCATGTTT"), self.mutable_s)

    def test_delete_stride_slice(self):
        del self.mutable_s[4 : 6 - 1]
        self.assertEqual(MutableSeq("TCAAAGGATGCATCATG"), self.mutable_s)

    def test_extract_third_nucleotide(self):
        """Test extracting every third nucleotide (slicing with stride 3)."""
        self.assertEqual(MutableSeq("TAGTAA"), self.mutable_s[0::3])
        self.assertEqual(MutableSeq("CAGGTT"), self.mutable_s[1::3])
        self.assertEqual(MutableSeq("AAACCG"), self.mutable_s[2::3])

    def test_set_wobble_codon_to_n(self):
        """Test setting wobble codon to N (set slice with stride 3)."""
        self.mutable_s[2::3] = "N" * len(self.mutable_s[2::3])
        self.assertEqual(MutableSeq("TCNAANGGNTGNATNATN"), self.mutable_s)


class TestUnknownSeq(unittest.TestCase):
    def setUp(self):
        self.s = Seq.UnknownSeq(6)

    def test_construction(self):
        self.assertEqual("??????", str(Seq.UnknownSeq(6)))
        self.assertEqual("NNNNNN", str(Seq.UnknownSeq(6, character="N")))
        self.assertEqual("XXXXXX", str(Seq.UnknownSeq(6, character="X")))
        self.assertEqual("??????", str(Seq.UnknownSeq(6, character="?")))

        with self.assertRaises(ValueError):
            Seq.UnknownSeq(-10)

        with self.assertRaises(ValueError):
            Seq.UnknownSeq(6, character="??")

    def test_length(self):
        self.assertEqual(6, len(self.s))

    def test_repr(self):
        self.assertEqual("UnknownSeq(6, character='?')", repr(self.s))

    def test_add_method(self):
        seq1 = Seq.UnknownSeq(3, character="N")
        self.assertEqual("??????NNN", str(self.s + seq1))

        seq2 = Seq.UnknownSeq(3, character="N")
        self.assertEqual("NNNNNN", str(seq1 + seq2))

    def test_getitem_method(self):
        self.assertEqual("", self.s[-1:-1])
        self.assertEqual("?", self.s[1])
        self.assertEqual("?", self.s[5:])
        self.assertEqual("?", self.s[:1])
        self.assertEqual("??", self.s[1:3])
        self.assertEqual("???", self.s[1:6:2])
        self.assertEqual("????", self.s[1:-1])
        with self.assertRaises(ValueError):
            self.s[1:6:0]

    def test_count(self):
        self.assertEqual(6, self.s.count("?"))
        self.assertEqual(3, self.s.count("??"))
        self.assertEqual(0, Seq.UnknownSeq(6, character="N").count("?"))
        self.assertEqual(0, Seq.UnknownSeq(6, character="N").count("??"))
        self.assertEqual(4, Seq.UnknownSeq(6, character="?").count("?", start=2))
        self.assertEqual(2, Seq.UnknownSeq(6, character="?").count("??", start=2))

    def test_complement(self):
        self.s.complement()
        self.assertEqual("??????", str(self.s))

    def test_reverse_complement(self):
        self.s.reverse_complement()
        self.assertEqual("??????", str(self.s))

    def test_transcribe(self):
        self.assertEqual("??????", self.s.transcribe())

    def test_back_transcribe(self):
        self.assertEqual("??????", self.s.back_transcribe())

    def test_upper(self):
        seq = Seq.UnknownSeq(6, character="N")
        self.assertEqual("NNNNNN", str(seq.upper()))

    def test_lower(self):
        seq = Seq.UnknownSeq(6, character="N")
        self.assertEqual("nnnnnn", str(seq.lower()))

    def test_translation(self):
        self.assertEqual("XX", str(self.s.translate()))

    def test_ungap(self):
        seq = Seq.UnknownSeq(7, character="N")
        self.assertEqual("NNNNNNN", str(seq.ungap("-")))

        seq = Seq.UnknownSeq(20, character="-")
        self.assertEqual("", seq.ungap("-"))


class TestAmbiguousComplements(unittest.TestCase):
    def test_ambiguous_values(self):
        """Test that other tests do not introduce characters to our values."""
        self.assertNotIn("-", ambiguous_dna_values)
        self.assertNotIn("?", ambiguous_dna_values)


class TestComplement(unittest.TestCase):
    def test_complement_ambiguous_dna_values(self):
        for ambig_char, values in sorted(ambiguous_dna_values.items()):
            compl_values = str(Seq.Seq(values).complement())
            ambig_values = ambiguous_dna_values[ambiguous_dna_complement[ambig_char]]
            self.assertEqual(set(compl_values), set(ambig_values))

    def test_complement_ambiguous_rna_values(self):
        for ambig_char, values in sorted(ambiguous_rna_values.items()):
            # Will default to DNA if neither T nor U found...
            compl_values = str(Seq.Seq(values).complement().transcribe())
            ambig_values = ambiguous_rna_values[ambiguous_rna_complement[ambig_char]]
            self.assertEqual(set(compl_values), set(ambig_values))

    def test_complement_incompatible_letters(self):
        seq = Seq.Seq("CAGGTU")
        with self.assertRaises(ValueError):
            seq.complement()

    def test_complement_of_mixed_dna_rna(self):
        seq = "AUGAAACTG"  # U and T
        self.assertRaises(ValueError, Seq.complement, seq)

    def test_complement_of_rna(self):
        seq = "AUGAAACUG"
        self.assertEqual("UACUUUGAC", Seq.complement(seq))

    def test_complement_of_dna(self):
        seq = "ATGAAACTG"
        self.assertEqual("TACTTTGAC", Seq.complement(seq))


class TestReverseComplement(unittest.TestCase):
    def test_reverse_complement(self):
        test_seqs_copy = copy.copy(test_seqs)
        test_seqs_copy.pop(13)

        for nucleotide_seq in test_seqs_copy:
            if isinstance(nucleotide_seq, Seq.Seq):
                expected = Seq.reverse_complement(nucleotide_seq)
                self.assertEqual(
                    repr(expected), repr(nucleotide_seq.reverse_complement())
                )
                self.assertEqual(
                    repr(expected[::-1]), repr(nucleotide_seq.complement())
                )
                self.assertEqual(
                    str(nucleotide_seq.complement()),
                    str(Seq.reverse_complement(nucleotide_seq))[::-1],
                )
                self.assertEqual(
                    str(nucleotide_seq.reverse_complement()),
                    str(Seq.reverse_complement(nucleotide_seq)),
                )

    def test_reverse_complement_of_mixed_dna_rna(self):
        seq = "AUGAAACTG"  # U and T
        self.assertRaises(ValueError, Seq.reverse_complement, seq)

    def test_reverse_complement_of_rna(self):
        seq = "AUGAAACUG"
        self.assertEqual("CAGUUUCAU", Seq.reverse_complement(seq))

    def test_reverse_complement_of_dna(self):
        seq = "ATGAAACTG"
        self.assertEqual("CAGTTTCAT", Seq.reverse_complement(seq))


class TestDoubleReverseComplement(unittest.TestCase):
    def test_reverse_complements(self):
        """Test double reverse complement preserves the sequence."""
        sorted_amb_rna = sorted(ambiguous_rna_values)
        sorted_amb_dna = sorted(ambiguous_dna_values)
        for sequence in [
            Seq.Seq("".join(sorted_amb_rna)),
            Seq.Seq("".join(sorted_amb_dna)),
            Seq.Seq("".join(sorted_amb_rna).replace("X", "")),
            Seq.Seq("".join(sorted_amb_dna).replace("X", "")),
            Seq.Seq("AWGAARCKG"),
        ]:  # Note no U or T
            reversed_sequence = sequence.reverse_complement()
            self.assertEqual(str(sequence), str(reversed_sequence.reverse_complement()))


class TestTranscription(unittest.TestCase):
    def test_transcription_dna_into_rna(self):
        for nucleotide_seq in test_seqs:
            expected = Seq.transcribe(nucleotide_seq)
            self.assertEqual(
                str(nucleotide_seq).replace("t", "u").replace("T", "U"), str(expected),
            )

    def test_transcription_dna_string_into_rna(self):
        seq = "ATGAAACTG"
        self.assertEqual("AUGAAACUG", Seq.transcribe(seq))

    def test_seq_object_transcription_method(self):
        for nucleotide_seq in test_seqs:
            if isinstance(nucleotide_seq, Seq.Seq):
                self.assertEqual(
                    repr(Seq.transcribe(nucleotide_seq)),
                    repr(nucleotide_seq.transcribe()),
                )

    def test_back_transcribe_rna_into_dna(self):
        for nucleotide_seq in test_seqs:
            expected = Seq.back_transcribe(nucleotide_seq)
            self.assertEqual(
                str(nucleotide_seq).replace("u", "t").replace("U", "T"), str(expected),
            )

    def test_back_transcribe_rna_string_into_dna(self):
        seq = "AUGAAACUG"
        self.assertEqual("ATGAAACTG", Seq.back_transcribe(seq))

    def test_seq_object_back_transcription_method(self):
        for nucleotide_seq in test_seqs:
            if isinstance(nucleotide_seq, Seq.Seq):
                expected = Seq.back_transcribe(nucleotide_seq)
                self.assertEqual(repr(nucleotide_seq.back_transcribe()), repr(expected))


class TestTranslating(unittest.TestCase):
    def setUp(self):
        self.test_seqs = [
            Seq.Seq("TCAAAAGGATGCATCATG"),
            Seq.Seq("ATGAAACTG"),
            Seq.Seq("ATGAARCTG"),
            Seq.Seq("AWGAARCKG"),  # Note no U or T
            Seq.Seq("".join(ambiguous_rna_values)),
            Seq.Seq("".join(ambiguous_dna_values)),
            Seq.Seq("AUGAAACUG"),
            Seq.Seq("ATGAAACTGWN"),
            Seq.Seq("AUGAAACUGWN"),
            Seq.MutableSeq("ATGAAACTG"),
            Seq.MutableSeq("AUGaaaCUG"),
        ]

    def test_translation(self):
        for nucleotide_seq in self.test_seqs:
            nucleotide_seq = nucleotide_seq[: 3 * (len(nucleotide_seq) // 3)]
            if isinstance(nucleotide_seq, Seq.Seq) and "X" not in str(nucleotide_seq):
                expected = Seq.translate(nucleotide_seq)
                self.assertEqual(repr(expected), repr(nucleotide_seq.translate()))

    def test_gapped_seq_with_gap_char_given(self):
        seq = Seq.Seq("ATG---AAACTG")
        self.assertEqual("M-KL", seq.translate(gap="-"))
        self.assertRaises(TranslationError, seq.translate, gap="~")

        seq = Seq.Seq("GTG---GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
        self.assertEqual("V-AIVMGR*KGAR*", seq.translate(gap="-"))
        self.assertRaises(TranslationError, seq.translate, gap=None)

        seq = Seq.Seq("ATG~~~AAACTG")
        self.assertRaises(TranslationError, seq.translate, gap="-")

        seq = Seq.Seq("ATG---AAACTGTAG")
        self.assertEqual("M-KL*", seq.translate(gap="-"))
        self.assertEqual("M-KL@", seq.translate(gap="-", stop_symbol="@"))
        self.assertRaises(TranslationError, seq.translate, gap="~")

        seq = Seq.Seq("ATG~~~AAACTGTAG")
        self.assertRaises(TranslationError, seq.translate, gap="-")

    def test_gapped_seq_no_gap_char_given(self):
        seq = Seq.Seq("ATG---AAACTG")
        self.assertRaises(TranslationError, seq.translate, gap=None)

    def test_translation_wrong_type(self):
        """Test translation table cannot be CodonTable."""
        seq = Seq.Seq("ATCGTA")
        with self.assertRaises(ValueError):
            seq.translate(table=ambiguous_dna_complement)

    def test_translation_of_string(self):
        seq = "GTGGCCATTGTAATGGGCCGC"
        self.assertEqual("VAIVMGR", Seq.translate(seq))

    def test_translation_of_gapped_string_with_gap_char_given(self):
        seq = "GTG---GCCATTGTAATGGGCCGC"
        expected = "V-AIVMGR"
        self.assertEqual(expected, Seq.translate(seq, gap="-"))
        self.assertRaises(TypeError, Seq.translate, seq, gap=[])
        self.assertRaises(ValueError, Seq.translate, seq, gap="-*")

    def test_translation_of_gapped_string_no_gap_char_given(self):
        seq = "GTG---GCCATTGTAATGGGCCGC"
        self.assertRaises(TranslationError, Seq.translate, seq)

    def test_translation_to_stop(self):
        for nucleotide_seq in self.test_seqs:
            nucleotide_seq = nucleotide_seq[: 3 * (len(nucleotide_seq) // 3)]
            if isinstance(nucleotide_seq, Seq.Seq) and "X" not in str(nucleotide_seq):
                short = Seq.translate(nucleotide_seq, to_stop=True)
                self.assertEqual(
                    str(short), str(Seq.translate(nucleotide_seq).split("*")[0])
                )

        seq = "GTGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
        self.assertEqual("VAIVMGRWKGAR", Seq.translate(seq, table=2, to_stop=True))

    def test_translation_on_proteins(self):
        """Check translation fails on a protein."""
        for s in protein_seqs:
            with self.assertRaises(TranslationError):
                Seq.translate(s)

            if isinstance(s, Seq.Seq):
                with self.assertRaises(TranslationError):
                    s.translate()

    def test_translation_of_invalid_codon(self):
        for codon in ["TA?", "N-N", "AC_", "Ac_"]:
            with self.assertRaises(TranslationError):
                Seq.translate(codon)

    def test_translation_of_glutamine(self):
        for codon in ["SAR", "SAG", "SAA"]:
            self.assertEqual("Z", Seq.translate(codon))

    def test_translation_of_asparagine(self):
        for codon in ["RAY", "RAT", "RAC"]:
            self.assertEqual("B", Seq.translate(codon))

    def test_translation_of_leucine(self):
        for codon in ["WTA", "MTY", "MTT", "MTW", "MTM", "MTH", "MTA", "MTC", "HTA"]:
            self.assertEqual("J", Seq.translate(codon))

    def test_translation_with_bad_table_argument(self):
        table = {}
        with self.assertRaises(ValueError):
            Seq.translate("GTGGCCATTGTAATGGGCCGC", table=table)

    def test_translation_with_codon_table_as_table_argument(self):
        table = standard_dna_table
        self.assertEqual("VAIVMGR", Seq.translate("GTGGCCATTGTAATGGGCCGC", table=table))

    def test_translation_incomplete_codon(self):
        with self.assertWarns(BiopythonWarning):
            Seq.translate("GTGGCCATTGTAATGGGCCG")

    def test_translation_extra_stop_codon(self):
        seq = "GTGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGTAG"
        with self.assertRaises(TranslationError):
            Seq.translate(seq, table=2, cds=True)

    def test_translation_using_cds(self):
        seq = "GTGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
        self.assertEqual("MAIVMGRWKGAR", Seq.translate(seq, table=2, cds=True))

        seq = "GTGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCG"  # not multiple of three
        with self.assertRaises(TranslationError):
            Seq.translate(seq, table=2, cds=True)

        seq = "GTGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA"  # no stop codon
        with self.assertRaises(TranslationError):
            Seq.translate(seq, table=2, cds=True)

        seq = "GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"  # no start codon
        with self.assertRaises(TranslationError):
            Seq.translate(seq, table=2, cds=True)

    def test_translation_using_tables_with_ambiguous_stop_codons(self):
        """Check for error and warning messages.

        Here, 'ambiguous stop codons' means codons of unambiguous sequence
        but with a context sensitive encoding as STOP or an amino acid.
        Thus, these codons appear within the codon table in the forward
        table as well as in the list of stop codons.
        """
        seq = "ATGGGCTGA"
        with self.assertRaises(ValueError):
            Seq.translate(seq, table=28, to_stop=True)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            Seq.translate(seq, table=28)
            message = str(w[-1].message)
            self.assertTrue(message.startswith("This table contains"))
            self.assertTrue(message.endswith("be translated as amino acid."))


class TestStopCodons(unittest.TestCase):
    def setUp(self):
        self.misc_stops = "TAATAGTGAAGAAGG"

    def test_stops(self):
        for nucleotide_seq in [self.misc_stops, Seq.Seq(self.misc_stops)]:
            self.assertEqual("***RR", str(Seq.translate(nucleotide_seq)))
            self.assertEqual("***RR", str(Seq.translate(nucleotide_seq, table=1)))
            self.assertEqual("***RR", str(Seq.translate(nucleotide_seq, table="SGC0")))
            self.assertEqual("**W**", str(Seq.translate(nucleotide_seq, table=2)))
            self.assertEqual(
                "**WRR", str(Seq.translate(nucleotide_seq, table="Yeast Mitochondrial"))
            )
            self.assertEqual("**WSS", str(Seq.translate(nucleotide_seq, table=5)))
            self.assertEqual("**WSS", str(Seq.translate(nucleotide_seq, table=9)))
            self.assertEqual(
                "**CRR", str(Seq.translate(nucleotide_seq, table="Euplotid Nuclear"))
            )
            self.assertEqual("***RR", str(Seq.translate(nucleotide_seq, table=11)))
            self.assertEqual(
                "***RR", str(Seq.translate(nucleotide_seq, table="Bacterial"))
            )

    def test_translation_of_stops(self):
        self.assertEqual(Seq.translate("TAT"), "Y")
        self.assertEqual(Seq.translate("TAR"), "*")
        self.assertEqual(Seq.translate("TAN"), "X")
        self.assertEqual(Seq.translate("NNN"), "X")

        self.assertEqual(Seq.translate("TAt"), "Y")
        self.assertEqual(Seq.translate("TaR"), "*")
        self.assertEqual(Seq.translate("TaN"), "X")
        self.assertEqual(Seq.translate("nnN"), "X")

        self.assertEqual(Seq.translate("tat"), "Y")
        self.assertEqual(Seq.translate("tar"), "*")
        self.assertEqual(Seq.translate("tan"), "X")
        self.assertEqual(Seq.translate("nnn"), "X")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
