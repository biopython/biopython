# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for seq module."""

import array
import copy
import unittest
import warnings

try:
    import numpy
except ImportError:
    numpy = None

from Bio import BiopythonWarning, BiopythonDeprecationWarning
from Bio import Seq
from Bio.Data.IUPACData import (
    ambiguous_dna_complement,
    ambiguous_rna_complement,
    ambiguous_dna_values,
    ambiguous_rna_values,
)
from Bio.Data.CodonTable import TranslationError, standard_dna_table

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
        self.assertEqual("TCAAAAGGATGCATCATG", self.s)

    def test_seq_construction(self):
        """Test Seq object initialization."""
        sequence = bytes(self.s)
        s = Seq.Seq(sequence)
        self.assertIsInstance(s, Seq.Seq, "Creating MutableSeq using bytes")
        self.assertEqual(s, self.s)
        s = Seq.Seq(bytearray(sequence))
        self.assertIsInstance(s, Seq.Seq, "Creating MutableSeq using bytearray")
        self.assertEqual(s, self.s)
        s = Seq.Seq(sequence.decode("ASCII"))
        self.assertIsInstance(s, Seq.Seq, "Creating MutableSeq using str")
        self.assertEqual(s, self.s)
        s = Seq.Seq(self.s)
        self.assertIsInstance(s, Seq.Seq, "Creating MutableSeq using Seq")
        self.assertEqual(s, self.s)
        s = Seq.Seq(Seq.MutableSeq(sequence))
        self.assertIsInstance(s, Seq.Seq, "Creating MutableSeq using MutableSeq")
        self.assertEqual(s, self.s)
        self.assertRaises(
            UnicodeEncodeError, Seq.Seq, "ÄþÇÐ"
        )  # All are Latin-1 characters
        self.assertRaises(UnicodeEncodeError, Seq.Seq, "あいうえお")  # These are not

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
        self.assertEqual("AA", self.s[3:5])

    def test_reverse(self):
        """Test reverse using -1 stride."""
        self.assertEqual("GTACTACGTAGGAAAACT", self.s[::-1])

    def test_extract_third_nucleotide(self):
        """Test extracting every third nucleotide (slicing with stride 3)."""
        self.assertEqual("TAGTAA", self.s[0::3])
        self.assertEqual("CAGGTT", self.s[1::3])
        self.assertEqual("AAACCG", self.s[2::3])

    def test_concatenation_of_seq(self):
        t = Seq.Seq("T")
        u = self.s + t
        self.assertEqual(str(self.s) + "T", u)
        self.assertEqual(self.s + Seq.Seq("T"), "TCAAAAGGATGCATCATGT")

    def test_replace(self):
        self.assertEqual("ATCCCA", Seq.Seq("ATC-CCA").replace("-", ""))

    def test_cast_to_list(self):
        self.assertEqual(list("ATC"), list(Seq.Seq("ATC")))
        self.assertEqual(list("ATC"), list(Seq.MutableSeq("ATC")))
        self.assertEqual(list(""), list(Seq.MutableSeq("")))
        self.assertEqual(list(""), list(Seq.Seq("")))
        with self.assertRaises(Seq.UndefinedSequenceError):
            list(Seq.Seq(None, length=3))
        with self.assertRaises(Seq.UndefinedSequenceError):
            list(Seq.Seq({3: "ACGT"}, length=10))


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
            self.assertEqual(a.lower(), str(a).lower())
            self.assertEqual(a.upper(), str(a).upper())
            self.assertEqual(a.islower(), str(a).islower())
            self.assertEqual(a.isupper(), str(a).isupper())
            self.assertEqual(a.strip(), str(a).strip())
            self.assertEqual(a.lstrip(), str(a).lstrip())
            self.assertEqual(a.rstrip(), str(a).rstrip())

    def test_mutableseq_upper_lower(self):
        seq = Seq.MutableSeq("ACgt")
        lseq = seq.lower()
        self.assertEqual(lseq, "acgt")
        self.assertEqual(seq, "ACgt")
        self.assertTrue(lseq.islower())
        self.assertFalse(seq.islower())
        lseq = seq.lower(inplace=False)
        self.assertEqual(lseq, "acgt")
        self.assertEqual(seq, "ACgt")
        self.assertTrue(lseq.islower())
        self.assertFalse(seq.islower())
        lseq = seq.lower(inplace=True)
        self.assertEqual(lseq, "acgt")
        self.assertIs(lseq, seq)
        self.assertTrue(lseq.islower())
        self.assertTrue(lseq.islower())
        seq = Seq.MutableSeq("ACgt")
        useq = seq.upper()
        self.assertEqual(useq, "ACGT")
        self.assertEqual(seq, "ACgt")
        self.assertTrue(useq.isupper())
        self.assertFalse(seq.isupper())
        useq = seq.upper(inplace=False)
        self.assertEqual(useq, "ACGT")
        self.assertEqual(seq, "ACgt")
        self.assertTrue(useq.isupper())
        self.assertFalse(seq.isupper())
        useq = seq.upper(inplace=True)
        self.assertEqual(useq, "ACGT")
        self.assertIs(useq, seq)
        self.assertTrue(useq.isupper())
        self.assertTrue(seq.isupper())

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

    def test_radd_method_using_wrong_object(self):
        self.assertEqual(self.s.__radd__({}), NotImplemented)

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
                self.assertEqual(a.strip(char), str(a).strip(str_char))
                self.assertEqual(a.lstrip(char), str(a).lstrip(str_char))
                self.assertEqual(a.rstrip(char), str(a).rstrip(str_char))

    def test_finding_characters(self):
        for a in self.dna + self.rna + self.nuc + self.protein:
            for char in self.test_chars:
                str_char = str(char)
                self.assertEqual(a.find(char), str(a).find(str_char))
                self.assertEqual(a.find(char, 2, -2), str(a).find(str_char, 2, -2))
                self.assertEqual(a.rfind(char), str(a).rfind(str_char))
                self.assertEqual(a.rfind(char, 2, -2), str(a).rfind(str_char, 2, -2))

    def test_counting_characters(self):
        from Bio.SeqRecord import SeqRecord

        for a in self.dna + self.rna + self.nuc + self.protein:
            r = SeqRecord(a)
            for char in self.test_chars:
                str_char = str(char)
                n = str(a).count(str_char)
                self.assertEqual(a.count(char), n)
                self.assertEqual(r.count(char), n)
                n = str(a).count(str_char, 2, -2)
                self.assertEqual(a.count(char, 2, -2), n)
                self.assertEqual(r.count(char, 2, -2), n)

    def test_splits(self):
        for a in self.dna + self.rna + self.nuc + self.protein:
            for char in self.test_chars:
                str_char = str(char)
                self.assertEqual(a.split(char), str(a).split(str_char))
                self.assertEqual(a.rsplit(char), str(a).rsplit(str_char))

                for max_sep in [0, 1, 2, 999]:
                    self.assertEqual(
                        a.split(char, max_sep), str(a).split(str_char, max_sep)
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
        self.nuc = [Seq.Seq("ATCG"), "UUUTTTACG"]
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
                self.assertEqual(c, str(a) + str(b))

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
                self.assertEqual(c, str(a) + str(b))

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
                self.assertEqual(c, str(a) + str(b))

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
                self.assertEqual(c, str(a) + str(b))

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
                self.assertEqual(c, str(a) + str(b))

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
        if numpy is not None:
            factor = numpy.intc(3)  # numpy integer
            for seq in test_seqs + protein_seqs:
                self.assertEqual(seq * factor, seq + seq + seq)

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
        if numpy is not None:
            factor = numpy.intc(3)  # numpy integer
            for seq in test_seqs + protein_seqs:
                self.assertEqual(factor * seq, seq + seq + seq)

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
        if numpy is not None:
            factor = numpy.intc(3)  # numpy integer
            for seq in test_seqs + protein_seqs:
                original_seq = seq * 1  # make a copy
                seq *= factor
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
        sequence = b"TCAAAAGGATGCATCATG"
        self.s = Seq.Seq(sequence)
        self.mutable_s = Seq.MutableSeq(sequence)

    def test_mutableseq_construction(self):
        """Test MutableSeq object initialization."""
        sequence = bytes(self.s)
        mutable_s = Seq.MutableSeq(sequence)
        self.assertIsInstance(
            mutable_s, Seq.MutableSeq, "Initializing MutableSeq from bytes"
        )
        self.assertEqual(mutable_s, self.s)
        mutable_s = Seq.MutableSeq(bytearray(sequence))
        self.assertIsInstance(
            mutable_s, Seq.MutableSeq, "Initializing MutableSeq from bytearray"
        )
        self.assertEqual(mutable_s, self.s)
        mutable_s = Seq.MutableSeq(sequence.decode("ASCII"))
        self.assertIsInstance(
            mutable_s, Seq.MutableSeq, "Initializing MutableSeq from str"
        )
        self.assertEqual(mutable_s, self.s)
        mutable_s = Seq.MutableSeq(self.s)
        self.assertIsInstance(
            mutable_s, Seq.MutableSeq, "Initializing MutableSeq from Seq"
        )
        self.assertEqual(mutable_s, self.s)
        mutable_s = Seq.MutableSeq(Seq.MutableSeq(sequence))
        self.assertEqual(mutable_s, self.s)
        self.assertIsInstance(
            mutable_s, Seq.MutableSeq, "Initializing MutableSeq from MutableSeq"
        )
        # Deprecated:
        with self.assertWarns(BiopythonDeprecationWarning):
            mutable_s = Seq.MutableSeq(array.array("u", sequence.decode("ASCII")))
        self.assertIsInstance(
            mutable_s, Seq.MutableSeq, "Creating MutableSeq using array"
        )
        self.assertEqual(mutable_s, self.s)
        self.assertRaises(
            UnicodeEncodeError, Seq.MutableSeq, "ÄþÇÐ"
        )  # All are Latin-1 characters
        self.assertRaises(UnicodeEncodeError, Seq.MutableSeq, "あいうえお")  # These are not

    def test_repr(self):
        self.assertEqual("MutableSeq('TCAAAAGGATGCATCATG')", repr(self.mutable_s))

    def test_truncated_repr(self):
        seq = "TCAAAAGGATGCATCATGTCAAAAGGATGCATCATGTCAAAAGGATGCATCATGTCAAAAGGA"
        expected = (
            "MutableSeq('TCAAAAGGATGCATCATGTCAAAAGGATGCATCATGTCAAAAGGATGCATCATG...GGA')"
        )
        self.assertEqual(expected, repr(Seq.MutableSeq(seq)))

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

    def test_radd_method_wrong_type(self):
        self.assertEqual(self.mutable_s.__radd__(1234), NotImplemented)

    def test_contains_method(self):
        self.assertIn("AAAA", self.mutable_s)

    def test_startswith(self):
        self.assertTrue(self.mutable_s.startswith("TCA"))
        self.assertTrue(self.mutable_s.startswith(("CAA", "CTA"), 1))

    def test_endswith(self):
        self.assertTrue(self.mutable_s.endswith("ATG"))
        self.assertTrue(self.mutable_s.endswith(("ATG", "CTA")))

    def test_as_string(self):
        self.assertEqual("TCAAAAGGATGCATCATG", self.mutable_s)

    def test_length(self):
        self.assertEqual(18, len(self.mutable_s))

    def test_converting_to_immutable(self):
        self.assertIsInstance(Seq.Seq(self.mutable_s), Seq.Seq)

    def test_first_nucleotide(self):
        self.assertEqual("T", self.mutable_s[0])

    def test_setting_slices(self):
        self.assertEqual(
            Seq.MutableSeq("CAAA"), self.mutable_s[1:5], "Slice mutable seq"
        )

        self.mutable_s[1:3] = "GAT"
        self.assertEqual(
            Seq.MutableSeq("TGATAAAGGATGCATCATG"),
            self.mutable_s,
            "Set slice with string and adding extra nucleotide",
        )

        self.mutable_s[1:3] = self.mutable_s[5:7]
        self.assertEqual(
            Seq.MutableSeq("TAATAAAGGATGCATCATG"),
            self.mutable_s,
            "Set slice with MutableSeq",
        )
        if numpy is not None:
            one, three, five, seven = numpy.array([1, 3, 5, 7])  # numpy integers
            self.assertEqual(
                Seq.MutableSeq("AATA"), self.mutable_s[one:five], "Slice mutable seq"
            )

            self.mutable_s[one:three] = "GAT"
            self.assertEqual(
                Seq.MutableSeq("TGATTAAAGGATGCATCATG"),
                self.mutable_s,
                "Set slice with string and adding extra nucleotide",
            )

            self.mutable_s[one:three] = self.mutable_s[five:seven]
            self.assertEqual(
                Seq.MutableSeq("TAATTAAAGGATGCATCATG"),
                self.mutable_s,
                "Set slice with MutableSeq",
            )

    def test_setting_item(self):
        self.mutable_s[3] = "G"
        self.assertEqual(Seq.MutableSeq("TCAGAAGGATGCATCATG"), self.mutable_s)
        if numpy is not None:
            i = numpy.intc(3)
            self.mutable_s[i] = "X"
            self.assertEqual(Seq.MutableSeq("TCAXAAGGATGCATCATG"), self.mutable_s)

    def test_deleting_slice(self):
        del self.mutable_s[4:5]
        self.assertEqual(Seq.MutableSeq("TCAAAGGATGCATCATG"), self.mutable_s)

    def test_deleting_item(self):
        del self.mutable_s[3]
        self.assertEqual(Seq.MutableSeq("TCAAAGGATGCATCATG"), self.mutable_s)

    def test_appending(self):
        self.mutable_s.append("C")
        self.assertEqual(Seq.MutableSeq("TCAAAAGGATGCATCATGC"), self.mutable_s)

    def test_inserting(self):
        self.mutable_s.insert(4, "G")
        self.assertEqual(Seq.MutableSeq("TCAAGAAGGATGCATCATG"), self.mutable_s)

    def test_popping_last_item(self):
        self.assertEqual("G", self.mutable_s.pop())

    def test_remove_items(self):
        self.mutable_s.remove("G")
        self.assertEqual(
            Seq.MutableSeq("TCAAAAGATGCATCATG"), self.mutable_s, "Remove first G"
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
        self.assertEqual(Seq.MutableSeq("GTACTACGTAGGAAAACT"), self.mutable_s)

    def test_reverse_with_stride(self):
        """Test reverse using -1 stride."""
        self.assertEqual(Seq.MutableSeq("GTACTACGTAGGAAAACT"), self.mutable_s[::-1])

    def test_complement_old(self):
        # old approach
        with self.assertWarns(BiopythonDeprecationWarning):
            self.mutable_s.complement()
        self.assertEqual("AGTTTTCCTACGTAGTAC", self.mutable_s)

    def test_complement(self):
        # new approach
        self.mutable_s.complement(inplace=True)
        self.assertEqual("AGTTTTCCTACGTAGTAC", self.mutable_s)

    def test_complement_rna(self):
        m = self.mutable_s.complement_rna()
        self.assertEqual(self.mutable_s, "TCAAAAGGATGCATCATG")
        self.assertIsInstance(m, Seq.MutableSeq)
        self.assertEqual(m, "AGUUUUCCUACGUAGUAC")
        m = self.mutable_s.complement_rna(inplace=True)
        self.assertEqual(self.mutable_s, "AGUUUUCCUACGUAGUAC")
        self.assertIsInstance(m, Seq.MutableSeq)
        self.assertEqual(m, "AGUUUUCCUACGUAGUAC")

    def test_reverse_complement_rna(self):
        m = self.mutable_s.reverse_complement_rna()
        self.assertEqual(self.mutable_s, "TCAAAAGGATGCATCATG")
        self.assertIsInstance(m, Seq.MutableSeq)
        self.assertEqual(m, "CAUGAUGCAUCCUUUUGA")
        m = self.mutable_s.reverse_complement_rna(inplace=True)
        self.assertEqual(self.mutable_s, "CAUGAUGCAUCCUUUUGA")
        self.assertIsInstance(m, Seq.MutableSeq)
        self.assertEqual(m, "CAUGAUGCAUCCUUUUGA")

    def test_transcribe(self):
        r = self.mutable_s.transcribe()
        self.assertEqual(self.mutable_s, "TCAAAAGGATGCATCATG")
        self.assertIsInstance(r, Seq.MutableSeq)
        self.assertEqual(r, "UCAAAAGGAUGCAUCAUG")
        r = self.mutable_s.transcribe(inplace=True)
        self.assertEqual(self.mutable_s, "UCAAAAGGAUGCAUCAUG")
        self.assertIsInstance(r, Seq.MutableSeq)
        self.assertEqual(r, "UCAAAAGGAUGCAUCAUG")
        d = self.mutable_s.back_transcribe()
        self.assertEqual(self.mutable_s, "UCAAAAGGAUGCAUCAUG")
        self.assertIsInstance(d, Seq.MutableSeq)
        self.assertEqual(d, "TCAAAAGGATGCATCATG")
        d = self.mutable_s.back_transcribe(inplace=True)
        self.assertEqual(self.mutable_s, "TCAAAAGGATGCATCATG")
        self.assertIsInstance(d, Seq.MutableSeq)
        self.assertEqual(d, "TCAAAAGGATGCATCATG")

    def test_complement_mixed_aphabets(self):
        # new approach
        seq = Seq.MutableSeq("AUGaaaCTG")
        seq.complement_rna(inplace=True)
        self.assertEqual("UACuuuGAC", seq)
        # old approach
        seq = Seq.MutableSeq("AUGaaaCTG")
        with self.assertWarns(BiopythonDeprecationWarning):
            with self.assertRaises(ValueError):
                seq.complement()

    def test_complement_rna_string(self):
        # new approach
        seq = Seq.MutableSeq("AUGaaaCUG")
        seq.complement_rna(inplace=True)
        self.assertEqual("UACuuuGAC", seq)
        # old approach
        seq = Seq.MutableSeq("AUGaaaCUG")
        with self.assertWarns(BiopythonDeprecationWarning):
            seq.complement()
        self.assertEqual("UACuuuGAC", seq)

    def test_complement_dna_string(self):
        # new approach
        seq = Seq.MutableSeq("ATGaaaCTG")
        seq.complement(inplace=True)
        self.assertEqual("TACtttGAC", seq)
        # old approach
        seq = Seq.MutableSeq("ATGaaaCTG")
        with self.assertWarns(BiopythonDeprecationWarning):
            seq.complement()
        self.assertEqual("TACtttGAC", seq)

    def test_reverse_complement(self):
        # new approach
        self.mutable_s.reverse_complement(inplace=True)
        self.assertEqual("CATGATGCATCCTTTTGA", self.mutable_s)

    def test_reverse_complement_old(self):
        # old approach
        with self.assertWarns(BiopythonDeprecationWarning):
            self.mutable_s.reverse_complement()
        self.assertEqual("CATGATGCATCCTTTTGA", self.mutable_s)

    def test_extend_method(self):
        self.mutable_s.extend("GAT")
        self.assertEqual(Seq.MutableSeq("TCAAAAGGATGCATCATGGAT"), self.mutable_s)

    def test_extend_with_mutable_seq(self):
        self.mutable_s.extend(Seq.MutableSeq("TTT"))
        self.assertEqual(Seq.MutableSeq("TCAAAAGGATGCATCATGTTT"), self.mutable_s)

    def test_delete_stride_slice(self):
        del self.mutable_s[4 : 6 - 1]
        self.assertEqual(Seq.MutableSeq("TCAAAGGATGCATCATG"), self.mutable_s)

    def test_extract_third_nucleotide(self):
        """Test extracting every third nucleotide (slicing with stride 3)."""
        self.assertEqual(Seq.MutableSeq("TAGTAA"), self.mutable_s[0::3])
        self.assertEqual(Seq.MutableSeq("CAGGTT"), self.mutable_s[1::3])
        self.assertEqual(Seq.MutableSeq("AAACCG"), self.mutable_s[2::3])

    def test_set_wobble_codon_to_n(self):
        """Test setting wobble codon to N (set slice with stride 3)."""
        self.mutable_s[2::3] = "N" * len(self.mutable_s[2::3])
        self.assertEqual(Seq.MutableSeq("TCNAANGGNTGNATNATN"), self.mutable_s)
        if numpy is not None:
            start, step = numpy.array([2, 3])  # numpy integers
            self.mutable_s[start::step] = "X" * len(self.mutable_s[2::3])
            self.assertEqual(Seq.MutableSeq("TCXAAXGGXTGXATXATX"), self.mutable_s)


class TestUnknownSeq(unittest.TestCase):
    def setUp(self):
        warnings.simplefilter("ignore", BiopythonDeprecationWarning)
        self.s = Seq.UnknownSeq(6)
        self.u = Seq.Seq(None, length=6)

    def tearDown(self):
        warnings.simplefilter("default", BiopythonDeprecationWarning)

    def test_unknownseq_construction(self):
        self.assertEqual("??????", Seq.UnknownSeq(6))
        self.assertEqual("NNNNNN", Seq.UnknownSeq(6, character="N"))
        self.assertEqual("XXXXXX", Seq.UnknownSeq(6, character="X"))
        self.assertEqual("??????", Seq.UnknownSeq(6, character="?"))
        with self.assertRaises(ValueError):
            "??????" == self.u
        with self.assertRaises(ValueError):
            self.u == "??????"

        with self.assertRaises(ValueError):
            Seq.UnknownSeq(-10)

        with self.assertRaises(ValueError):
            Seq.Seq(None, length=-10)

        with self.assertRaises(ValueError):
            Seq.UnknownSeq(6, character="??")

    def test_length(self):
        self.assertEqual(6, len(self.s))
        self.assertEqual(6, len(self.u))

    def test_repr(self):
        self.assertEqual("UnknownSeq(6, character='?')", repr(self.s))
        self.assertEqual("Seq(None, length=6)", repr(self.u))

    def test_add_method(self):
        seq1 = Seq.UnknownSeq(3, character="N")
        self.assertEqual("??????NNN", self.s + seq1)

        seq2 = Seq.UnknownSeq(3, character="N")
        self.assertEqual("NNNNNN", seq1 + seq2)

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
        with self.assertRaises(ValueError):
            self.u[1:6:0]

    def test_count(self):
        self.assertEqual(6, self.s.count("?"))
        self.assertEqual(3, self.s.count("??"))
        self.assertEqual(0, Seq.UnknownSeq(6, character="N").count("?"))
        self.assertEqual(0, Seq.UnknownSeq(6, character="N").count("??"))
        self.assertEqual(4, Seq.UnknownSeq(6, character="?").count("?", start=2))
        self.assertEqual(2, Seq.UnknownSeq(6, character="?").count("??", start=2))
        self.assertRaises(ValueError, self.u.count, "?")

    def test_complement(self):
        self.s.complement()
        self.assertEqual("??????", self.s)
        t = self.u.complement()
        self.assertEqual(len(t), 6)
        self.assertRaises(ValueError, str, t)

    def test_reverse_complement(self):
        self.s.reverse_complement()
        self.assertEqual("??????", self.s)
        t = self.u.reverse_complement()
        self.assertEqual(len(t), 6)
        self.assertRaises(ValueError, str, t)

    def test_transcribe(self):
        self.assertEqual("??????", self.s.transcribe())
        t = self.u.transcribe()
        self.assertEqual(len(t), 6)
        self.assertRaises(ValueError, str, t)

    def test_back_transcribe(self):
        self.assertEqual("??????", self.s.back_transcribe())
        t = self.u.back_transcribe()
        self.assertEqual(len(t), 6)
        self.assertRaises(ValueError, str, t)

    def test_upper(self):
        seq = Seq.UnknownSeq(6, character="N")
        self.assertEqual("NNNNNN", seq.upper())
        self.assertEqual("Seq(None, length=6)", repr(self.u.upper()))

    def test_lower(self):
        seq = Seq.UnknownSeq(6, character="N")
        self.assertEqual("nnnnnn", seq.lower())
        self.assertEqual("Seq(None, length=6)", repr(self.u.lower()))

    def test_translation(self):
        self.assertEqual("XX", self.s.translate())
        t = self.u.translate()
        self.assertEqual(len(t), 2)
        self.assertRaises(ValueError, str, t)

    def test_ungap(self):
        seq = Seq.UnknownSeq(7, character="N")
        self.assertEqual("NNNNNNN", seq.ungap("-"))

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
            compl_values = Seq.Seq(values).complement()
            ambig_values = ambiguous_dna_values[ambiguous_dna_complement[ambig_char]]
            self.assertCountEqual(compl_values, ambig_values)

    def test_complement_ambiguous_rna_values(self):
        for ambig_char, values in sorted(ambiguous_rna_values.items()):
            # Will default to DNA if neither T nor U found...
            if "u" in values or "U" in values:
                compl_values = Seq.Seq(values).complement_rna().transcribe()
            else:
                compl_values = Seq.Seq(values).complement().transcribe()
            ambig_values = ambiguous_rna_values[ambiguous_rna_complement[ambig_char]]
            self.assertCountEqual(compl_values, ambig_values)

    def test_complement_incompatible_letters(self):
        seq = Seq.Seq("CAGGTU")
        # new approach
        dna = seq.complement(inplace=False)  # TODO: remove inplace=False
        self.assertEqual("GTCCAA", dna)
        rna = seq.complement_rna()
        self.assertEqual("GUCCAA", rna)
        # old approach
        with self.assertWarns(BiopythonDeprecationWarning):
            with self.assertRaises(ValueError):
                seq.complement()

    def test_complement_of_mixed_dna_rna(self):
        seq = "AUGAAACTG"  # U and T
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=BiopythonDeprecationWarning)
            self.assertRaises(ValueError, Seq.complement, seq)

    def test_complement_of_rna(self):
        seq = "AUGAAACUG"
        # new approach
        rna = Seq.complement_rna(seq)
        self.assertEqual("UACUUUGAC", rna)
        # old approach
        with self.assertWarns(BiopythonDeprecationWarning):
            rna = Seq.complement(seq)
        self.assertEqual("UACUUUGAC", rna)

    def test_complement_of_dna(self):
        seq = "ATGAAACTG"
        self.assertEqual("TACTTTGAC", Seq.complement(seq))

    def test_immutable(self):
        from Bio.SeqRecord import SeqRecord

        r = SeqRecord(Seq.Seq("ACGT"))
        with self.assertRaises(TypeError) as cm:
            Seq.complement(r, inplace=True)
        self.assertEqual(str(cm.exception), "SeqRecords are immutable")
        with self.assertRaises(TypeError) as cm:
            Seq.complement("ACGT", inplace=True)
        self.assertEqual(str(cm.exception), "strings are immutable")
        with self.assertRaises(TypeError) as cm:
            Seq.complement_rna(r, inplace=True)
        self.assertEqual(str(cm.exception), "SeqRecords are immutable")
        with self.assertRaises(TypeError) as cm:
            Seq.complement_rna("ACGT", inplace=True)
        self.assertEqual(str(cm.exception), "strings are immutable")


class TestReverseComplement(unittest.TestCase):
    def test_reverse_complement(self):
        test_seqs_copy = copy.copy(test_seqs)
        test_seqs_copy.pop(13)

        for nucleotide_seq in test_seqs_copy:
            if not isinstance(nucleotide_seq, Seq.Seq):
                continue
            if "u" in nucleotide_seq or "U" in nucleotide_seq:
                expected = Seq.reverse_complement_rna(nucleotide_seq)
                self.assertEqual(
                    repr(expected), repr(nucleotide_seq.reverse_complement_rna())
                )
                self.assertEqual(
                    repr(expected[::-1]), repr(nucleotide_seq.complement_rna())
                )
                self.assertEqual(
                    nucleotide_seq.complement_rna(),
                    Seq.reverse_complement_rna(nucleotide_seq)[::-1],
                )
                self.assertEqual(
                    nucleotide_seq.reverse_complement_rna(),
                    Seq.reverse_complement_rna(nucleotide_seq),
                )
            else:
                expected = Seq.reverse_complement(nucleotide_seq)
                self.assertEqual(
                    repr(expected), repr(nucleotide_seq.reverse_complement())
                )
                self.assertEqual(
                    repr(expected[::-1]), repr(nucleotide_seq.complement())
                )
                self.assertEqual(
                    nucleotide_seq.complement(),
                    Seq.reverse_complement(nucleotide_seq)[::-1],
                )
                self.assertEqual(
                    nucleotide_seq.reverse_complement(),
                    Seq.reverse_complement(nucleotide_seq),
                )

    def test_reverse_complement_of_mixed_dna_rna(self):
        seq = "AUGAAACTG"  # U and T
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=BiopythonDeprecationWarning)
            self.assertRaises(ValueError, Seq.reverse_complement, seq)

    def test_reverse_complement_of_rna(self):
        # old approach
        seq = "AUGAAACUG"
        with self.assertWarns(BiopythonDeprecationWarning):
            rna = Seq.reverse_complement(seq)
        self.assertEqual("CAGUUUCAU", rna)
        # new approach
        dna = Seq.reverse_complement(seq, inplace=False)  # TODO: remove inplace=False
        self.assertEqual("CAGTTTCAT", dna)

    def test_reverse_complement_of_dna(self):
        seq = "ATGAAACTG"
        self.assertEqual("CAGTTTCAT", Seq.reverse_complement(seq))

    def test_immutable(self):
        from Bio.SeqRecord import SeqRecord

        r = SeqRecord(Seq.Seq("ACGT"))
        with self.assertRaises(TypeError) as cm:
            Seq.reverse_complement(r, inplace=True)
        self.assertEqual(str(cm.exception), "SeqRecords are immutable")
        with self.assertRaises(TypeError) as cm:
            Seq.reverse_complement("ACGT", inplace=True)
        self.assertEqual(str(cm.exception), "strings are immutable")
        with self.assertRaises(TypeError) as cm:
            Seq.reverse_complement_rna(r, inplace=True)
        self.assertEqual(str(cm.exception), "SeqRecords are immutable")
        with self.assertRaises(TypeError) as cm:
            Seq.reverse_complement_rna("ACGT", inplace=True)
        self.assertEqual(str(cm.exception), "strings are immutable")


class TestDoubleReverseComplement(unittest.TestCase):
    def test_reverse_complements(self):
        """Test double reverse complement preserves the sequence."""
        sorted_amb_rna = sorted(ambiguous_rna_values)
        sorted_amb_dna = sorted(ambiguous_dna_values)
        for sequence in [
            Seq.Seq("".join(sorted_amb_dna)),
            Seq.Seq("".join(sorted_amb_dna).replace("X", "")),
            Seq.Seq("AWGAARCKG"),  # Note no U or T
        ]:
            reversed_sequence = sequence.reverse_complement()
            self.assertEqual(sequence, reversed_sequence.reverse_complement())
        for sequence in [
            Seq.Seq("".join(sorted_amb_rna)),
            Seq.Seq("".join(sorted_amb_rna).replace("X", "")),
            Seq.Seq("AWGAARCKG"),  # Note no U or T
        ]:
            reversed_sequence = sequence.reverse_complement_rna()
            self.assertEqual(sequence, reversed_sequence.reverse_complement_rna())


class TestTranscription(unittest.TestCase):
    def test_transcription_dna_into_rna(self):
        for nucleotide_seq in test_seqs:
            expected = Seq.transcribe(nucleotide_seq)
            self.assertEqual(
                str(nucleotide_seq).replace("t", "u").replace("T", "U"), expected
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
                str(nucleotide_seq).replace("u", "t").replace("U", "T"), expected
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
            if "X" not in nucleotide_seq:
                expected = Seq.translate(nucleotide_seq)
                self.assertEqual(expected, nucleotide_seq.translate())

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
            if "X" not in nucleotide_seq:
                short = Seq.translate(nucleotide_seq, to_stop=True)
                self.assertEqual(short, Seq.translate(nucleotide_seq).split("*")[0])

        seq = "GTGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
        self.assertEqual("VAIVMGRWKGAR", Seq.translate(seq, table=2, to_stop=True))

    def test_translation_on_proteins(self):
        """Check translation fails on a protein."""
        for s in protein_seqs:
            if len(s) % 3 != 0:
                with self.assertWarns(BiopythonWarning):
                    with self.assertRaises(TranslationError):
                        Seq.translate(s)

                with self.assertWarns(BiopythonWarning):
                    with self.assertRaises(TranslationError):
                        s.translate()
            else:
                with self.assertRaises(TranslationError):
                    Seq.translate(s)

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
        with self.assertRaises(ValueError) as cm:
            Seq.translate("GTGGCCATTGTAATGGGCCGC", table=table)
        self.assertEqual(str(cm.exception), "Bad table argument")
        table = b"0x"
        with self.assertRaises(TypeError) as cm:
            Seq.translate("GTGGCCATTGTAATGGGCCGC", table=table)
        self.assertEqual(str(cm.exception), "table argument must be integer or string")

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
            self.assertEqual("***RR", Seq.translate(nucleotide_seq))
            self.assertEqual("***RR", Seq.translate(nucleotide_seq, table=1))
            self.assertEqual("***RR", Seq.translate(nucleotide_seq, table="SGC0"))
            self.assertEqual("**W**", Seq.translate(nucleotide_seq, table=2))
            self.assertEqual(
                "**WRR", Seq.translate(nucleotide_seq, table="Yeast Mitochondrial")
            )
            self.assertEqual("**WSS", Seq.translate(nucleotide_seq, table=5))
            self.assertEqual("**WSS", Seq.translate(nucleotide_seq, table=9))
            self.assertEqual(
                "**CRR", Seq.translate(nucleotide_seq, table="Euplotid Nuclear")
            )
            self.assertEqual("***RR", Seq.translate(nucleotide_seq, table=11))
            self.assertEqual("***RR", Seq.translate(nucleotide_seq, table="Bacterial"))

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


class TestAttributes(unittest.TestCase):
    def test_seq(self):
        s = Seq.Seq("ACGT")
        with self.assertRaises(AttributeError):
            s.dog
        s.dog = "woof"
        self.assertIn("dog", dir(s))
        self.assertEqual(s.dog, "woof")
        del s.dog
        with self.assertRaises(AttributeError):
            s.dog
        self.assertNotIn("dog", dir(s))
        with self.assertRaises(AttributeError):
            s.cat
        s.dog = "woof"
        s.cat = "meow"
        self.assertIn("dog", dir(s))
        self.assertIn("cat", dir(s))
        self.assertEqual(s.dog, "woof")
        self.assertEqual(s.cat, "meow")
        del s.dog
        with self.assertRaises(AttributeError):
            s.dog
        self.assertNotIn("dog", dir(s))
        self.assertIn("cat", dir(s))
        self.assertEqual(s.cat, "meow")
        del s.cat
        with self.assertRaises(AttributeError):
            s.cat
        self.assertNotIn("cat", dir(s))
        s.dog = "woof"
        s.dog = "bark"
        self.assertIn("dog", dir(s))
        self.assertEqual(s.dog, "bark")
        del s.dog
        with self.assertRaises(AttributeError):
            s.dog
        self.assertNotIn("dog", dir(s))

    def test_mutable_seq(self):
        s = Seq.MutableSeq("ACGT")
        with self.assertRaises(AttributeError):
            s.dog
        s.dog = "woof"
        self.assertIn("dog", dir(s))
        self.assertEqual(s.dog, "woof")
        del s.dog
        with self.assertRaises(AttributeError):
            s.dog
        self.assertNotIn("dog", dir(s))
        with self.assertRaises(AttributeError):
            s.cat
        s.dog = "woof"
        s.cat = "meow"
        self.assertIn("dog", dir(s))
        self.assertIn("cat", dir(s))
        self.assertEqual(s.dog, "woof")
        self.assertEqual(s.cat, "meow")
        del s.dog
        with self.assertRaises(AttributeError):
            s.dog
        self.assertNotIn("dog", dir(s))
        self.assertIn("cat", dir(s))
        self.assertEqual(s.cat, "meow")
        del s.cat
        with self.assertRaises(AttributeError):
            s.cat
        self.assertNotIn("cat", dir(s))
        s.dog = "woof"
        s.dog = "bark"
        self.assertIn("dog", dir(s))
        self.assertEqual(s.dog, "bark")
        del s.dog
        with self.assertRaises(AttributeError):
            s.dog
        self.assertNotIn("dog", dir(s))


class TestSeqDefined(unittest.TestCase):
    def test_zero_length(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=BiopythonDeprecationWarning)
            zero_length_seqs = [
                Seq.Seq(""),
                Seq.Seq(None, length=0),
                Seq.Seq({}, length=0),
                Seq.UnknownSeq(length=0),
                Seq.MutableSeq(""),
            ]

        for seq in zero_length_seqs:
            self.assertTrue(seq.defined, msg=repr(seq))
            self.assertEqual(seq.defined_ranges, (), msg=repr(seq))

    def test_undefined(self):
        seq = Seq.Seq(None, length=1)
        self.assertFalse(seq.defined)
        self.assertEqual(seq.defined_ranges, ())
        seq = Seq.Seq({3: "ACGT"}, length=10)
        self.assertFalse(seq.defined)
        self.assertEqual(seq.defined_ranges, ((3, 7),))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=BiopythonDeprecationWarning)
            seq = Seq.UnknownSeq(length=1)
        self.assertFalse(seq.defined)
        self.assertEqual(seq.defined_ranges, ())

    def test_defined(self):
        seqs = [
            Seq.Seq("T"),
            Seq.Seq({0: "A"}, length=1),
            Seq.Seq({0: "A", 1: "C"}, length=2),
        ]

        for seq in seqs:
            self.assertTrue(seq.defined, msg=repr(seq))
            self.assertEqual(seq.defined_ranges, ((0, len(seq)),), msg=repr(seq))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
