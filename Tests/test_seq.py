# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from __future__ import print_function
import array
import copy
import sys
import warnings

# Remove unittest2 import after dropping support for Python 2
if sys.version_info[0] < 3:
    try:
        import unittest2 as unittest
    except ImportError:
        from Bio import MissingPythonDependencyError
        raise MissingPythonDependencyError("Under Python 2 this test needs the unittest2 library")
else:
    import unittest

from Bio import BiopythonWarning
from Bio import Alphabet
from Bio import Seq
from Bio.Alphabet import IUPAC, Gapped
from Bio.Data.IUPACData import (ambiguous_dna_complement,
                                ambiguous_rna_complement,
                                ambiguous_dna_values, ambiguous_rna_values)
from Bio.Data.CodonTable import TranslationError, standard_dna_table
from Bio.Seq import MutableSeq


if sys.version_info[0] == 3:
    array_indicator = "u"
else:
    array_indicator = "c"

test_seqs = [
    Seq.Seq("TCAAAAGGATGCATCATG", IUPAC.unambiguous_dna),
    Seq.Seq("T", IUPAC.ambiguous_dna),
    Seq.Seq("ATGAAACTG"),
    Seq.Seq("ATGAARCTG"),
    Seq.Seq("AWGAARCKG"),  # Note no U or T
    Seq.Seq("".join(ambiguous_rna_values)),
    Seq.Seq("".join(ambiguous_dna_values)),
    Seq.Seq("".join(ambiguous_rna_values), Alphabet.generic_rna),
    Seq.Seq("".join(ambiguous_dna_values), Alphabet.generic_dna),
    Seq.Seq("".join(ambiguous_rna_values), IUPAC.IUPACAmbiguousRNA()),
    Seq.Seq("".join(ambiguous_dna_values), IUPAC.IUPACAmbiguousDNA()),
    Seq.Seq("AWGAARCKG", Alphabet.generic_dna),
    Seq.Seq("AUGAAACUG", Alphabet.generic_rna),
    Seq.Seq("ATGAAACTG", IUPAC.unambiguous_dna),
    Seq.Seq("ATGAAA-CTG", Alphabet.Gapped(IUPAC.unambiguous_dna)),
    Seq.Seq("ATGAAACTGWN", IUPAC.ambiguous_dna),
    Seq.Seq("AUGAAACUG", Alphabet.generic_rna),
    Seq.Seq("AUGAAA==CUG", Alphabet.Gapped(Alphabet.generic_rna, "=")),
    Seq.Seq("AUGAAACUG", IUPAC.unambiguous_rna),
    Seq.Seq("AUGAAACUGWN", IUPAC.ambiguous_rna),
    Seq.Seq("ATGAAACTG", Alphabet.generic_nucleotide),
    Seq.Seq("AUGAAACTG", Alphabet.generic_nucleotide),  # U and T
    Seq.MutableSeq("ATGAAACTG", Alphabet.generic_dna),
    Seq.MutableSeq("AUGaaaCUG", IUPAC.unambiguous_rna),
    Seq.Seq("ACTGTCGTCT", Alphabet.generic_protein),
]
protein_seqs = [
    Seq.Seq("ATCGPK", IUPAC.protein),
    Seq.Seq("T.CGPK", Alphabet.Gapped(IUPAC.protein, ".")),
    Seq.Seq("T-CGPK", Alphabet.Gapped(IUPAC.protein, "-")),
    Seq.Seq("MEDG-KRXR*",
            Alphabet.Gapped(Alphabet.HasStopCodon(IUPAC.extended_protein, "*"),
                            "-")),
    Seq.MutableSeq("ME-K-DRXR*XU",
                   Alphabet.Gapped(Alphabet.HasStopCodon(
                       IUPAC.extended_protein, "*"), "-")),
    Seq.Seq("MEDG-KRXR@",
            Alphabet.HasStopCodon(Alphabet.Gapped(IUPAC.extended_protein, "-"),
                                  "@")),
    Seq.Seq("ME-KR@",
            Alphabet.HasStopCodon(Alphabet.Gapped(IUPAC.protein, "-"), "@")),
    Seq.Seq("MEDG.KRXR@",
            Alphabet.Gapped(Alphabet.HasStopCodon(IUPAC.extended_protein, "@"),
                            ".")),
]


class TestSeq(unittest.TestCase):
    def setUp(self):
        self.s = Seq.Seq("TCAAAAGGATGCATCATG", IUPAC.unambiguous_dna)

    def test_as_string(self):
        """Test converting Seq to string"""
        self.assertEqual("TCAAAAGGATGCATCATG", str(self.s))

    def test_construction_using_a_seq_object(self):
        """Test using a Seq object to initialize another Seq object"""
        with self.assertRaises(TypeError):
            Seq.Seq(self.s)

    def test_repr(self):
        """Test representation of Seq object"""
        self.assertEqual("Seq('TCAAAAGGATGCATCATG', IUPACUnambiguousDNA())",
                         repr(self.s))

    def test_truncated_repr(self):
        seq = "TCAAAAGGATGCATCATGTCAAAAGGATGCATCATGTCAAAAGGATGCATCATGTCAAAAGGA"
        expected = "Seq('TCAAAAGGATGCATCATGTCAAAAGGATGCATCATGTCAAAAGGATGC" + \
                   "ATCATG...GGA', IUPACAmbiguousDNA())"
        self.assertEqual(expected, repr(Seq.Seq(seq, IUPAC.ambiguous_dna)))

    def test_length(self):
        """Test len method on Seq object"""
        self.assertEqual(18, len(self.s))

    def test_first_nucleotide(self):
        """Test getting first nucleotide of Seq"""
        self.assertEqual("T", self.s[0])

    def test_last_nucleotide(self):
        """Test getting last nucleotide of Seq"""
        self.assertEqual("G", self.s[-1])

    def test_slicing(self):
        """Test slicing of Seq"""
        self.assertEqual("AA", str(self.s[3:5]))

    def test_reverse(self):
        """Test reverse using -1 stride"""
        self.assertEqual("GTACTACGTAGGAAAACT", self.s[::-1])

    def test_extract_third_nucleotide(self):
        """Test extracting every third nucleotide (slicing with stride 3)"""
        self.assertEqual("TAGTAA", str(self.s[0::3]))
        self.assertEqual("CAGGTT", str(self.s[1::3]))
        self.assertEqual("AAACCG", str(self.s[2::3]))

    def test_alphabet_letters(self):
        """Test nucleotides in DNA Seq"""
        self.assertEqual("GATC", self.s.alphabet.letters)

    def test_alphabet(self):
        """Test alphabet of derived Seq object"""
        t = Seq.Seq("T", IUPAC.unambiguous_dna)
        u = self.s + t
        self.assertEqual("IUPACUnambiguousDNA()", str(u.alphabet))

    def test_length_concatenated_unambiguous_seq(self):
        """Test length of concatenated Seq object with unambiguous DNA"""
        t = Seq.Seq("T", IUPAC.unambiguous_dna)
        u = self.s + t
        self.assertEqual(19, len(u))

    def test_concatenation_of_seq(self):
        t = Seq.Seq("T", IUPAC.unambiguous_dna)
        u = self.s + t
        self.assertEqual(str(self.s) + "T", str(u))

    def test_concatenation_error(self):
        """Test DNA Seq objects cannot be concatenated with Protein Seq
        objects"""
        with self.assertRaises(TypeError):
            self.s + Seq.Seq("T", IUPAC.protein)

    def test_concatenation_of_ambiguous_and_unambiguous_dna(self):
        """Test concatenated Seq object with ambiguous and unambiguous DNA
        returns ambiguous Seq"""
        t = Seq.Seq("T", IUPAC.ambiguous_dna)
        u = self.s + t
        self.assertEqual("IUPACAmbiguousDNA()", str(u.alphabet))

    def test_ungap(self):
        self.assertEqual("ATCCCA", str(Seq.Seq("ATC-CCA").ungap("-")))

        with self.assertRaises(ValueError):
            Seq.Seq("ATC-CCA").ungap("--")

        with self.assertRaises(ValueError):
            Seq.Seq("ATC-CCA").ungap()


class TestSeqStringMethods(unittest.TestCase):
    def setUp(self):
        self.s = Seq.Seq("TCAAAAGGATGCATCATG", IUPAC.unambiguous_dna)
        self.dna = [
            Seq.Seq("ATCG", IUPAC.ambiguous_dna),
            Seq.Seq("gtca", Alphabet.generic_dna),
            Seq.MutableSeq("GGTCA", Alphabet.generic_dna),
            Seq.Seq("CTG-CA", Alphabet.Gapped(IUPAC.unambiguous_dna, "-")),
        ]
        self.rna = [
            Seq.Seq("AUUUCG", IUPAC.ambiguous_rna),
            Seq.MutableSeq("AUUCG", IUPAC.ambiguous_rna),
            Seq.Seq("uCAg", Alphabet.generic_rna),
            Seq.MutableSeq("UC-AG",
                           Alphabet.Gapped(Alphabet.generic_rna, "-")),
            Seq.Seq("U.CAG", Alphabet.Gapped(Alphabet.generic_rna, ".")),
        ]
        self.nuc = [Seq.Seq("ATCG", Alphabet.generic_nucleotide)]
        self.protein = [
            Seq.Seq("ATCGPK", IUPAC.protein),
            Seq.Seq("atcGPK", Alphabet.generic_protein),
            Seq.Seq("T.CGPK", Alphabet.Gapped(IUPAC.protein, ".")),
            Seq.Seq("T-CGPK", Alphabet.Gapped(IUPAC.protein, "-")),
            Seq.Seq("MEDG-KRXR*",
                    Alphabet.Gapped(
                        Alphabet.HasStopCodon(IUPAC.extended_protein, "*"),
                        "-")),
            Seq.MutableSeq("ME-K-DRXR*XU",
                           Alphabet.Gapped(
                               Alphabet.HasStopCodon(IUPAC.extended_protein,
                                                     "*"), "-")),
            Seq.Seq("MEDG-KRXR@",
                    Alphabet.HasStopCodon(
                        Alphabet.Gapped(IUPAC.extended_protein, "-"), "@")),
            Seq.Seq("ME-KR@",
                    Alphabet.HasStopCodon(Alphabet.Gapped(IUPAC.protein, "-"),
                                          "@")),
            Seq.Seq("MEDG.KRXR@",
                    Alphabet.Gapped(Alphabet.HasStopCodon(
                        IUPAC.extended_protein, "@"), ".")),
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

    def test_equal_comparison_of_incompatible_alphabets(self):
        """Test __eq__ comparison method"""
        with warnings.catch_warnings(record=True):
            Seq.Seq("TCAAAA", IUPAC.ambiguous_dna) == \
                              Seq.Seq("TCAAAA", IUPAC.ambiguous_rna)

    def test_not_equal_comparsion(self):
        """Test __ne__ comparison method"""
        self.assertNotEqual(Seq.Seq("TCAAA", IUPAC.ambiguous_dna),
                            Seq.Seq("TCAAAA", IUPAC.ambiguous_dna))

    def test_less_than_comparison_of_incompatible_alphabets(self):
        """Test __lt__ comparison method"""
        seq1 = Seq.Seq("TCAAA", IUPAC.ambiguous_dna)
        seq2 = Seq.Seq("UCAAAA", IUPAC.ambiguous_rna)
        with self.assertWarns(BiopythonWarning):
            self.assertTrue(seq1 < seq2)

    def test_less_than_or_equal_comparison_of_incompatible_alphabets(self):
        """Test __lt__ comparison method"""
        seq1 = Seq.Seq("TCAAA", IUPAC.ambiguous_dna)
        seq2 = Seq.Seq("UCAAAA", IUPAC.ambiguous_rna)
        with self.assertWarns(BiopythonWarning):
            self.assertTrue(seq1 <= seq2)

    def test_add_method_using_wrong_object(self):
        with self.assertRaises(TypeError):
            self.s + dict()

    def test_radd_method(self):
        self.assertEqual("TCAAAAGGATGCATCATGTCAAAAGGATGCATCATG",
                         str(self.s.__radd__(self.s)))

    def test_radd_method_using_incompatible_alphabets(self):
        rna_seq = Seq.Seq("UCAAAA", IUPAC.ambiguous_rna)
        with self.assertRaises(TypeError):
            self.s.__radd__(rna_seq)

    def test_radd_method_using_wrong_object(self):
        with self.assertRaises(TypeError):
            self.s.__radd__(dict())

    def test_to_string_deprecated_method(self):
        with self.assertWarns(BiopythonWarning):
            self.s.tostring()

    def test_contains_method(self):
        self.assertIn("AAAA", self.s)

    def test_startswith(self):
        self.assertTrue(self.s.startswith("TCA"))
        self.assertTrue(self.s.startswith(("CAA", "CTA"), 1))

    def test_endswith(self):
        self.assertTrue(self.s.endswith("ATG"))
        self.assertTrue(self.s.endswith(("ATG", "CTA")))

    def test_append_nucleotides(self):
        self.test_chars.append(Seq.Seq("A", IUPAC.ambiguous_dna))
        self.test_chars.append(Seq.Seq("A", IUPAC.ambiguous_rna))
        self.test_chars.append(Seq.Seq("A", Alphabet.generic_nucleotide))

        self.assertEqual(7, len(self.test_chars))

    def test_append_proteins(self):
        self.test_chars.append(Seq.Seq("K", Alphabet.generic_protein))
        self.test_chars.append(Seq.Seq("K-",
                                       Alphabet.Gapped(
                                           Alphabet.generic_protein, "-")))
        self.test_chars.append(Seq.Seq("K@",
                                       Alphabet.Gapped(IUPAC.protein, "@")))

        self.assertEqual(7, len(self.test_chars))

    def test_exception_when_clashing_alphabets(self):
        """Test by setting up clashing alphabet sequences"""
        b = Seq.Seq("-", Alphabet.generic_nucleotide)
        self.assertRaises(TypeError, self.protein[0].strip, b)

        b = Seq.Seq("-", Alphabet.generic_protein)
        self.assertRaises(TypeError, self.dna[0].strip, b)

    def test_stripping_characters(self):
        for a in self.dna + self.rna + self.nuc + self.protein:
            for char in self.test_chars:
                str_char = str(char)
                if isinstance(a, Seq.Seq):
                    self.assertEqual(str(a.strip(char)),
                                     str(a).strip(str_char))
                    self.assertEqual(str(a.lstrip(char)),
                                     str(a).lstrip(str_char))
                    self.assertEqual(str(a.rstrip(char)),
                                     str(a).rstrip(str_char))

    def test_finding_characters(self):
        for a in self.dna + self.rna + self.nuc + self.protein:
            for char in self.test_chars:
                str_char = str(char)
                if isinstance(a, Seq.Seq):
                    self.assertEqual(a.find(char), str(a).find(str_char))
                    self.assertEqual(a.find(char, 2, -2),
                                     str(a).find(str_char, 2, -2))
                    self.assertEqual(a.rfind(char), str(a).rfind(str_char))
                    self.assertEqual(a.rfind(char, 2, -2),
                                     str(a).rfind(str_char, 2, -2))

    def test_counting_characters(self):
        for a in self.dna + self.rna + self.nuc + self.protein:
            for char in self.test_chars:
                str_char = str(char)
                if isinstance(a, Seq.Seq):
                    self.assertEqual(a.count(char), str(a).count(str_char))
                    self.assertEqual(a.count(char, 2, -2),
                                     str(a).count(str_char, 2, -2))

    def test_splits(self):
        for a in self.dna + self.rna + self.nuc + self.protein:
            for char in self.test_chars:
                str_char = str(char)
                if isinstance(a, Seq.Seq):
                    self.assertEqual([str(x) for x in a.split(char)],
                                     str(a).split(str_char))
                    self.assertEqual([str(x) for x in a.rsplit(char)],
                                     str(a).rsplit(str_char))

                    for max_sep in [0, 1, 2, 999]:
                        self.assertEqual(
                            [str(x) for x in a.split(char, max_sep)],
                            str(a).split(str_char, max_sep))


class TestSeqAddition(unittest.TestCase):
    def setUp(self):
        self.dna = [
            Seq.Seq("ATCG", IUPAC.ambiguous_dna),
            Seq.Seq("gtca", Alphabet.generic_dna),
            Seq.MutableSeq("GGTCA", Alphabet.generic_dna),
            Seq.Seq("CTG-CA", Alphabet.Gapped(IUPAC.unambiguous_dna, "-")),
            "TGGTCA",
        ]
        self.rna = [
            Seq.Seq("AUUUCG", IUPAC.ambiguous_rna),
            Seq.MutableSeq("AUUCG", IUPAC.ambiguous_rna),
            Seq.Seq("uCAg", Alphabet.generic_rna),
            Seq.MutableSeq("UC-AG",
                           Alphabet.Gapped(Alphabet.generic_rna, "-")),
            Seq.Seq("U.CAG",
                    Alphabet.Gapped(Alphabet.generic_rna, ".")),
            "UGCAU",
        ]
        self.nuc = [
            Seq.Seq("ATCG", Alphabet.generic_nucleotide),
            "UUUTTTACG",
        ]
        self.protein = [
            Seq.Seq("ATCGPK", IUPAC.protein),
            Seq.Seq("atcGPK", Alphabet.generic_protein),
            Seq.Seq("T.CGPK", Alphabet.Gapped(IUPAC.protein, ".")),
            Seq.Seq("T-CGPK", Alphabet.Gapped(IUPAC.protein, "-")),
            Seq.Seq("MEDG-KRXR*",
                    Alphabet.Gapped(Alphabet.HasStopCodon(
                        IUPAC.extended_protein, "*"), "-")),
            Seq.MutableSeq("ME-K-DRXR*XU",
                           Alphabet.Gapped(Alphabet.HasStopCodon(
                               IUPAC.extended_protein, "*"), "-")),
            "TEDDF",
        ]

    def test_addition_dna_rna_with_generic_nucleotides(self):
        for a in self.dna + self.rna:
            for b in self.nuc:
                c = a + b
                self.assertEqual(str(c), str(a) + str(b))

    def test_addition_rna_with_rna(self):
        self.rna.pop(3)
        for a in self.rna:
            for b in self.rna:
                c = a + b
                self.assertEqual(str(c), str(a) + str(b))

    def test_exception_when_added_rna_has_more_than_one_gap_type(self):
        """Test resulting sequence has gap types '-' and '.'"""
        with self.assertRaises(ValueError):
            self.rna[3] + self.rna[4]

    def test_addition_dna_with_dna(self):
        for a in self.dna:
            for b in self.dna:
                c = a + b
                self.assertEqual(str(c), str(a) + str(b))

    def test_addition_dna_with_rna(self):
        self.dna.pop(4)
        self.rna.pop(5)
        for a in self.dna:
            for b in self.rna:
                with self.assertRaises(TypeError):
                    a + b
                with self.assertRaises(TypeError):
                    b + a

    def test_addition_proteins(self):
        self.protein.pop(2)
        for a in self.protein:
            for b in self.protein:
                c = a + b
                self.assertEqual(str(c), str(a) + str(b))

    def test_exception_when_added_protein_has_more_than_one_gap_type(self):
        """Test resulting protein has gap types '-' and '.'"""
        a = Seq.Seq("T.CGPK", Alphabet.Gapped(IUPAC.protein, "."))
        b = Seq.Seq("T-CGPK", Alphabet.Gapped(IUPAC.protein, "-"))
        with self.assertRaises(ValueError):
            a + b

    def test_exception_when_added_protein_has_several_stop_codon_types(self):
        """Test resulting protein has stop codon types '*' and '@'"""
        a = Seq.Seq("MEDG-KRXR@", Alphabet.HasStopCodon(
            Alphabet.Gapped(IUPAC.extended_protein, "-"), "@"))
        b = Seq.Seq("MEDG-KRXR*", Alphabet.Gapped(
            Alphabet.HasStopCodon(IUPAC.extended_protein, "*"), "-"))
        with self.assertRaises(ValueError):
            a + b

    def test_exception_when_adding_protein_with_nucleotides(self):
        for a in self.protein[0:5]:
            for b in self.dna[0:3] + self.rna[0:4]:
                with self.assertRaises(TypeError):
                    a + b

    def test_adding_generic_nucleotide_with_other_nucleotides(self):
        for a in self.nuc:
            for b in self.dna + self.rna + self.nuc:
                c = a + b
                self.assertEqual(str(c), str(a) + str(b))


class TestMutableSeq(unittest.TestCase):
    def setUp(self):
        self.s = Seq.Seq("TCAAAAGGATGCATCATG", IUPAC.unambiguous_dna)
        self.mutable_s = MutableSeq("TCAAAAGGATGCATCATG", IUPAC.ambiguous_dna)

    def test_mutableseq_creation(self):
        """Test creating MutableSeqs in multiple ways"""
        mutable_s = MutableSeq("TCAAAAGGATGCATCATG", IUPAC.ambiguous_dna)
        self.assertIsInstance(mutable_s, MutableSeq, "Creating MutableSeq")

        mutable_s = self.s.tomutable()
        self.assertIsInstance(mutable_s, MutableSeq,
                              "Converting Seq to mutable")

        array_seq = MutableSeq(array.array(array_indicator,
                                           "TCAAAAGGATGCATCATG"),
                               IUPAC.ambiguous_dna)
        self.assertIsInstance(array_seq, MutableSeq,
                              "Creating MutableSeq using array")

    def test_repr(self):
        self.assertEqual(
            "MutableSeq('TCAAAAGGATGCATCATG', IUPACAmbiguousDNA())",
            repr(self.mutable_s))

    def test_truncated_repr(self):
        seq = "TCAAAAGGATGCATCATGTCAAAAGGATGCATCATGTCAAAAGGATGCATCATGTCAAAAGGA"
        expected = "MutableSeq('TCAAAAGGATGCATCATGTCAAAAGGATGCATCATGTCAAA" + \
                   "AGGATGCATCATG...GGA', IUPACAmbiguousDNA())"
        self.assertEqual(expected, repr(MutableSeq(seq, IUPAC.ambiguous_dna)))

    def test_equal_comparison(self):
        """Test __eq__ comparison method"""
        self.assertEqual(self.mutable_s, "TCAAAAGGATGCATCATG")

    def test_equal_comparison_of_incompatible_alphabets(self):
        with self.assertWarns(BiopythonWarning):
            self.mutable_s == MutableSeq('UCAAAAGGA', IUPAC.ambiguous_rna)

    def test_not_equal_comparison(self):
        """Test __ne__ comparison method"""
        self.assertNotEqual(self.mutable_s, "other thing")

    def test_less_than_comparison(self):
        """Test __lt__ comparison method"""
        self.assertTrue(self.mutable_s[:-1] < self.mutable_s)

    def test_less_than_comparison_of_incompatible_alphabets(self):
        with self.assertWarns(BiopythonWarning):
            self.mutable_s[:-1] < MutableSeq("UCAAAAGGAUGCAUCAUG",
                                             IUPAC.ambiguous_rna)

    def test_less_than_comparison_without_alphabet(self):
        self.assertTrue(self.mutable_s[:-1] < "TCAAAAGGATGCATCATG")

    def test_less_than_or_equal_comparison(self):
        """Test __le__ comparison method"""
        self.assertTrue(self.mutable_s[:-1] <= self.mutable_s)

    def test_less_than_or_equal_comparison_of_incompatible_alphabets(self):
        with self.assertWarns(BiopythonWarning):
            self.mutable_s[:-1] <= MutableSeq("UCAAAAGGAUGCAUCAUG",
                                              IUPAC.ambiguous_rna)

    def test_less_than_or_equal_comparison_without_alphabet(self):
        self.assertTrue(self.mutable_s[:-1] <= "TCAAAAGGATGCATCATG")

    def test_add_method(self):
        """Test adding wrong type to MutableSeq"""
        with self.assertRaises(TypeError):
            self.mutable_s + 1234

    def test_radd_method(self):
        self.assertEqual("TCAAAAGGATGCATCATGTCAAAAGGATGCATCATG",
                         self.mutable_s.__radd__(self.mutable_s))

    def test_radd_method_incompatible_alphabets(self):
        with self.assertRaises(TypeError):
            self.mutable_s.__radd__(MutableSeq("UCAAAAGGA",
                                               IUPAC.ambiguous_rna))

    def test_radd_method_using_seq_object(self):
        self.assertEqual("TCAAAAGGATGCATCATGTCAAAAGGATGCATCATG",
                         self.mutable_s.__radd__(self.s))

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
        self.assertEqual('T', self.mutable_s[0])

    def test_setting_slices(self):
        self.assertEqual(MutableSeq('CAAA', IUPAC.ambiguous_dna),
                         self.mutable_s[1:5], "Slice mutable seq")

        self.mutable_s[1:3] = "GAT"
        self.assertEqual(MutableSeq("TGATAAAGGATGCATCATG",
                                    IUPAC.ambiguous_dna),
                         self.mutable_s,
                         "Set slice with string and adding extra nucleotide")

        self.mutable_s[1:3] = self.mutable_s[5:7]
        self.assertEqual(MutableSeq("TAATAAAGGATGCATCATG",
                                    IUPAC.ambiguous_dna),
                         self.mutable_s, "Set slice with MutableSeq")

        self.mutable_s[1:3] = array.array(array_indicator, "GAT")
        self.assertEqual(MutableSeq("TGATTAAAGGATGCATCATG",
                                    IUPAC.ambiguous_dna),
                         self.mutable_s, "Set slice with array")

    def test_setting_item(self):
        self.mutable_s[3] = "G"
        self.assertEqual(MutableSeq("TCAGAAGGATGCATCATG", IUPAC.ambiguous_dna),
                         self.mutable_s)

    def test_deleting_slice(self):
        del self.mutable_s[4:5]
        self.assertEqual(MutableSeq("TCAAAGGATGCATCATG", IUPAC.ambiguous_dna),
                         self.mutable_s)

    def test_deleting_item(self):
        del self.mutable_s[3]
        self.assertEqual(MutableSeq("TCAAAGGATGCATCATG", IUPAC.ambiguous_dna),
                         self.mutable_s)

    def test_appending(self):
        self.mutable_s.append("C")
        self.assertEqual(MutableSeq("TCAAAAGGATGCATCATGC",
                                    IUPAC.ambiguous_dna),
                         self.mutable_s)

    def test_inserting(self):
        self.mutable_s.insert(4, "G")
        self.assertEqual(MutableSeq("TCAAGAAGGATGCATCATG",
                                    IUPAC.ambiguous_dna),
                         self.mutable_s)

    def test_popping_last_item(self):
        self.assertEqual("G", self.mutable_s.pop())

    def test_remove_items(self):
        self.mutable_s.remove("G")
        self.assertEqual(MutableSeq("TCAAAAGATGCATCATG", IUPAC.ambiguous_dna),
                         self.mutable_s, "Remove first G")

        self.assertRaises(ValueError, self.mutable_s.remove, 'Z')

    def test_count(self):
        self.assertEqual(7, self.mutable_s.count("A"))
        self.assertEqual(2, self.mutable_s.count("AA"))

    def test_index(self):
        self.assertEqual(2, self.mutable_s.index("A"))
        self.assertRaises(ValueError, self.mutable_s.index, "8888")

    def test_reverse(self):
        """Test using reverse method"""
        self.mutable_s.reverse()
        self.assertEqual(MutableSeq("GTACTACGTAGGAAAACT", IUPAC.ambiguous_dna),
                         self.mutable_s)

    def test_reverse_with_stride(self):
        """Test reverse using -1 stride"""
        self.assertEqual(MutableSeq("GTACTACGTAGGAAAACT", IUPAC.ambiguous_dna),
                         self.mutable_s[::-1])

    def test_complement(self):
        self.mutable_s.complement()
        self.assertEqual(str("AGTTTTCCTACGTAGTAC"), str(self.mutable_s))

    def test_complement_rna(self):
        seq = Seq.MutableSeq("AUGaaaCUG", IUPAC.unambiguous_rna)
        seq.complement()
        self.assertEqual(str("UACuuuGAC"), str(seq))

    def test_complement_mixed_aphabets(self):
        seq = Seq.MutableSeq("AUGaaaCTG")
        with self.assertRaises(ValueError):
            seq.complement()

    def test_complement_rna_string(self):
        seq = Seq.MutableSeq("AUGaaaCUG")
        seq.complement()
        self.assertEqual('UACuuuGAC', str(seq))

    def test_complement_dna_string(self):
        seq = Seq.MutableSeq("ATGaaaCTG")
        seq.complement()
        self.assertEqual('TACtttGAC', str(seq))

    def test_reverse_complement(self):
        self.mutable_s.reverse_complement()
        self.assertEqual("CATGATGCATCCTTTTGA", str(self.mutable_s))

    def test_reverse_complement_of_protein(self):
        seq = Seq.MutableSeq("ACTGTCGTCT", Alphabet.generic_protein)
        with self.assertRaises(ValueError):
            seq.reverse_complement()

    def test_to_string_method(self):
        """This method is currently deprecated, probably will need to remove
        this test soon"""
        with self.assertWarns(BiopythonWarning):
            self.mutable_s.tostring()

    def test_extend_method(self):
        self.mutable_s.extend("GAT")
        self.assertEqual(MutableSeq("TCAAAAGGATGCATCATGGAT",
                                    IUPAC.ambiguous_dna),
                         self.mutable_s)

    def test_extend_with_mutable_seq(self):
        self.mutable_s.extend(MutableSeq("TTT", IUPAC.ambiguous_dna))
        self.assertEqual(MutableSeq("TCAAAAGGATGCATCATGTTT",
                                    IUPAC.ambiguous_dna),
                         self.mutable_s)

    def test_delete_stride_slice(self):
        del self.mutable_s[4:6 - 1]
        self.assertEqual(MutableSeq("TCAAAGGATGCATCATG", IUPAC.ambiguous_dna),
                         self.mutable_s)

    def test_extract_third_nucleotide(self):
        """Test extracting every third nucleotide (slicing with stride 3)"""
        self.assertEqual(MutableSeq("TAGTAA", IUPAC.ambiguous_dna),
                         self.mutable_s[0::3])
        self.assertEqual(MutableSeq("CAGGTT", IUPAC.ambiguous_dna),
                         self.mutable_s[1::3])
        self.assertEqual(MutableSeq("AAACCG", IUPAC.ambiguous_dna),
                         self.mutable_s[2::3])

    def test_set_wobble_codon_to_n(self):
        """Test setting wobble codon to N (set slice with stride 3)"""
        self.mutable_s[2::3] = "N" * len(self.mutable_s[2::3])
        self.assertEqual(MutableSeq("TCNAANGGNTGNATNATN", IUPAC.ambiguous_dna),
                         self.mutable_s)


class TestUnknownSeq(unittest.TestCase):
    def setUp(self):
        self.s = Seq.UnknownSeq(6)

    def test_construction(self):
        self.assertEqual("??????", str(Seq.UnknownSeq(6)))
        self.assertEqual("NNNNNN",
                         str(Seq.UnknownSeq(6, Alphabet.generic_dna)))
        self.assertEqual("XXXXXX",
                         str(Seq.UnknownSeq(6, Alphabet.generic_protein)))
        self.assertEqual("??????", str(Seq.UnknownSeq(6, character="?")))

        with self.assertRaises(ValueError):
            Seq.UnknownSeq(-10)

        with self.assertRaises(ValueError):
            Seq.UnknownSeq(6, character='??')

    def test_length(self):
        self.assertEqual(6, len(self.s))

    def test_repr(self):
        self.assertEqual(
            "UnknownSeq(6, alphabet = Alphabet(), character = '?')",
            repr(self.s))

    def test_add_method(self):
        seq1 = Seq.UnknownSeq(3, Alphabet.generic_dna)
        self.assertEqual("??????NNN", str(self.s + seq1))

        seq2 = Seq.UnknownSeq(3, Alphabet.generic_dna)
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
        self.assertEqual(4,
                         Seq.UnknownSeq(6, character="?").count("?", start=2))
        self.assertEqual(2,
                         Seq.UnknownSeq(6, character="?").count("??", start=2))

    def test_complement(self):
        self.s.complement()
        self.assertEqual(str("??????"), str(self.s))

    def test_complement_of_protein(self):
        """Test reverse complement shouldn't work on a protein!"""
        seq = Seq.UnknownSeq(6, Alphabet.generic_protein)
        with self.assertRaises(ValueError):
            seq.complement()

    def test_reverse_complement(self):
        self.s.reverse_complement()
        self.assertEqual("??????", str(self.s))

    def test_reverse_complement_of_protein(self):
        seq = Seq.UnknownSeq(6, Alphabet.generic_protein)
        self.assertRaises(ValueError, seq.reverse_complement)

    def test_transcribe(self):
        self.assertEqual("??????", self.s.transcribe())

    def test_back_transcribe(self):
        self.assertEqual("??????", self.s.back_transcribe())

    def test_upper(self):
        seq = Seq.UnknownSeq(6, Alphabet.generic_dna)
        self.assertEqual("NNNNNN", str(seq.upper()))

    def test_lower(self):
        seq = Seq.UnknownSeq(6, Alphabet.generic_dna)
        self.assertEqual("nnnnnn", str(seq.lower()))

    def test_translation(self):
        self.assertEqual("XX", str(self.s.translate()))

    def test_translation_of_proteins(self):
        seq = Seq.UnknownSeq(6, IUPAC.protein)
        self.assertRaises(ValueError, seq.translate)

    def test_ungap(self):
        seq = Seq.UnknownSeq(7,
                             alphabet=Alphabet.Gapped(Alphabet.DNAAlphabet(),
                                                      "-"))
        self.assertEqual("NNNNNNN", str(seq.ungap("-")))

        seq = Seq.UnknownSeq(20,
                             alphabet=Alphabet.Gapped(Alphabet.DNAAlphabet(),
                                                      "-"), character='-')
        self.assertEqual("", seq.ungap("-"))


class TestAmbiguousComplements(unittest.TestCase):
    def test_ambiguous_values(self):
        """Test that other tests do not introduce characters to our values"""
        self.assertFalse("-" in ambiguous_dna_values)
        self.assertFalse("?" in ambiguous_dna_values)


class TestComplement(unittest.TestCase):
    def test_complement_ambiguous_dna_values(self):
        for ambig_char, values in sorted(ambiguous_dna_values.items()):
            compl_values = str(
                Seq.Seq(values, alphabet=IUPAC.ambiguous_dna).complement())
            ambig_values = (
                ambiguous_dna_values[ambiguous_dna_complement[ambig_char]])
            self.assertEqual(set(compl_values), set(ambig_values))

    def test_complement_ambiguous_rna_values(self):
        for ambig_char, values in sorted(ambiguous_rna_values.items()):
            compl_values = str(
                Seq.Seq(values, alphabet=IUPAC.ambiguous_rna).complement())
            ambig_values = (
                ambiguous_rna_values[ambiguous_rna_complement[ambig_char]])
            self.assertEqual(set(compl_values), set(ambig_values))

    def test_complement_incompatible_alphabets(self):
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

    def test_complement_on_proteins(self):
        """Test complement shouldn't work on a protein!"""
        for s in protein_seqs:
            with self.assertRaises(ValueError):
                Seq.complement(s)

            with self.assertRaises(ValueError):
                s.complement()


class TestReverseComplement(unittest.TestCase):
    def test_reverse_complement(self):
        test_seqs_copy = copy.copy(test_seqs)
        test_seqs_copy.pop(21)

        for nucleotide_seq in test_seqs_copy:
            if not isinstance(nucleotide_seq.alphabet,
                              Alphabet.ProteinAlphabet) and \
                              isinstance(nucleotide_seq, Seq.Seq):
                expected = Seq.reverse_complement(nucleotide_seq)
                self.assertEqual(
                    repr(expected), repr(nucleotide_seq.reverse_complement()))
                self.assertEqual(
                    repr(expected[::-1]), repr(nucleotide_seq.complement()))
                self.assertEqual(
                    str(nucleotide_seq.complement()),
                    str(Seq.reverse_complement(nucleotide_seq))[::-1])
                self.assertEqual(str(nucleotide_seq.reverse_complement()),
                                 str(Seq.reverse_complement(nucleotide_seq)))

    def test_reverse_complement_of_mixed_dna_rna(self):
        seq = "AUGAAACTG"  # U and T
        self.assertRaises(ValueError, Seq.reverse_complement, seq)

    def test_reverse_complement_of_rna(self):
        seq = "AUGAAACUG"
        self.assertEqual("CAGUUUCAU", Seq.reverse_complement(seq))

    def test_reverse_complement_of_dna(self):
        seq = "ATGAAACTG"
        self.assertEqual("CAGTTTCAT", Seq.reverse_complement(seq))

    def test_reverse_complement_on_proteins(self):
        """Test reverse complement shouldn't work on a protein!"""
        for s in protein_seqs:
            with self.assertRaises(ValueError):
                Seq.reverse_complement(s)

            with self.assertRaises(ValueError):
                s.reverse_complement()


class TestDoubleReverseComplement(unittest.TestCase):
    def test_reverse_complements(self):
        """Test double reverse complement preserves the sequence"""
        sorted_amb_rna = sorted(ambiguous_rna_values)
        sorted_amb_dna = sorted(ambiguous_dna_values)
        for sequence in [Seq.Seq("".join(sorted_amb_rna)),
                         Seq.Seq("".join(sorted_amb_dna)),
                         Seq.Seq("".join(sorted_amb_rna),
                                 Alphabet.generic_rna),
                         Seq.Seq("".join(sorted_amb_dna),
                                 Alphabet.generic_dna),
                         Seq.Seq("".join(sorted_amb_rna).replace("X", ""),
                                 IUPAC.IUPACAmbiguousRNA()),
                         Seq.Seq("".join(sorted_amb_dna).replace("X", ""),
                                 IUPAC.IUPACAmbiguousDNA()),
                         Seq.Seq("AWGAARCKG")]:  # Note no U or T
            reversed_sequence = sequence.reverse_complement()
            self.assertEqual(str(sequence),
                             str(reversed_sequence.reverse_complement()))


class TestSequenceAlphabets(unittest.TestCase):
    def test_sequence_alphabets(self):
        """Sanity test on the test sequence alphabets (see also enhancement
        bug 2597)"""
        for nucleotide_seq in test_seqs:
            if "U" in str(nucleotide_seq).upper():
                self.assertNotIsInstance(nucleotide_seq.alphabet,
                                         Alphabet.DNAAlphabet)
            if "T" in str(nucleotide_seq).upper():
                self.assertNotIsInstance(nucleotide_seq.alphabet,
                                         Alphabet.RNAAlphabet)


class TestTranscription(unittest.TestCase):
    def test_transcription_dna_into_rna(self):
        for nucleotide_seq in test_seqs:
            if isinstance(nucleotide_seq.alphabet, Alphabet.DNAAlphabet):
                expected = Seq.transcribe(nucleotide_seq)
                self.assertEqual(
                    str(nucleotide_seq).replace("t", "u").replace("T", "U"),
                    str(expected))

    def test_transcription_dna_string_into_rna(self):
        seq = "ATGAAACTG"
        self.assertEqual("AUGAAACUG", Seq.transcribe(seq))

    def test_seq_object_transcription_method(self):
        for nucleotide_seq in test_seqs:
            if isinstance(nucleotide_seq.alphabet, Alphabet.DNAAlphabet) and \
                    isinstance(nucleotide_seq, Seq.Seq):
                self.assertEqual(repr(Seq.transcribe(nucleotide_seq)),
                                 repr(nucleotide_seq.transcribe()))

    def test_transcription_of_rna(self):
        """Test transcription shouldn't work on RNA!"""
        seq = Seq.Seq("AUGAAACUG", IUPAC.ambiguous_rna)
        with self.assertRaises(ValueError):
            seq.transcribe()

    def test_transcription_of_proteins(self):
        """Test transcription shouldn't work on a protein!"""
        for s in protein_seqs:
            with self.assertRaises(ValueError):
                Seq.transcribe(s)

            if isinstance(s, Seq.Seq):
                with self.assertRaises(ValueError):
                    s.transcribe()

    def test_back_transcribe_rna_into_dna(self):
        for nucleotide_seq in test_seqs:
            if isinstance(nucleotide_seq.alphabet, Alphabet.RNAAlphabet):
                expected = Seq.back_transcribe(nucleotide_seq)
                self.assertEqual(
                    str(nucleotide_seq).replace("u", "t").replace("U", "T"),
                    str(expected))

    def test_back_transcribe_rna_string_into_dna(self):
        seq = "AUGAAACUG"
        self.assertEqual("ATGAAACTG", Seq.back_transcribe(seq))

    def test_seq_object_back_transcription_method(self):
        for nucleotide_seq in test_seqs:
            if isinstance(nucleotide_seq.alphabet, Alphabet.RNAAlphabet) and \
                    isinstance(nucleotide_seq, Seq.Seq):
                expected = Seq.back_transcribe(nucleotide_seq)
                self.assertEqual(repr(nucleotide_seq.back_transcribe()),
                                 repr(expected))

    def test_back_transcription_of_proteins(self):
        """Test back-transcription shouldn't work on a protein!"""
        for s in protein_seqs:
            with self.assertRaises(ValueError):
                Seq.back_transcribe(s)

            if isinstance(s, Seq.Seq):
                with self.assertRaises(ValueError):
                    s.back_transcribe()

    def test_back_transcription_of_dna(self):
        """Test back-transcription shouldn't work on DNA!"""
        seq = Seq.Seq("ATGAAACTG", IUPAC.ambiguous_dna)
        with self.assertRaises(ValueError):
            seq.back_transcribe()


class TestTranslating(unittest.TestCase):
    def setUp(self):
        self.test_seqs = [
            Seq.Seq("TCAAAAGGATGCATCATG", IUPAC.unambiguous_dna),
            Seq.Seq("ATGAAACTG"),
            Seq.Seq("ATGAARCTG"),
            Seq.Seq("AWGAARCKG"),  # Note no U or T
            Seq.Seq("".join(ambiguous_rna_values)),
            Seq.Seq("".join(ambiguous_dna_values)),
            Seq.Seq("".join(ambiguous_rna_values), Alphabet.generic_rna),
            Seq.Seq("".join(ambiguous_dna_values), Alphabet.generic_dna),
            Seq.Seq("".join(ambiguous_rna_values), IUPAC.IUPACAmbiguousRNA()),
            Seq.Seq("".join(ambiguous_dna_values), IUPAC.IUPACAmbiguousDNA()),
            Seq.Seq("AWGAARCKG", Alphabet.generic_dna),
            Seq.Seq("AUGAAACUG", Alphabet.generic_rna),
            Seq.Seq("ATGAAACTG", IUPAC.unambiguous_dna),
            Seq.Seq("ATGAAACTGWN", IUPAC.ambiguous_dna),
            Seq.Seq("AUGAAACUG", Alphabet.generic_rna),
            Seq.Seq("AUGAAACUG", IUPAC.unambiguous_rna),
            Seq.Seq("AUGAAACUGWN", IUPAC.ambiguous_rna),
            Seq.Seq("ATGAAACTG", Alphabet.generic_nucleotide),
            Seq.MutableSeq("ATGAAACTG", Alphabet.generic_dna),
            Seq.MutableSeq("AUGaaaCUG", IUPAC.unambiguous_rna),
        ]

    def test_translation(self):
        for nucleotide_seq in self.test_seqs:
            nucleotide_seq = nucleotide_seq[:3 * (len(nucleotide_seq) // 3)]
            if isinstance(nucleotide_seq, Seq.Seq) and \
               'X' not in str(nucleotide_seq):
                expected = Seq.translate(nucleotide_seq)
                self.assertEqual(repr(expected),
                                 repr(nucleotide_seq.translate()))

    def test_alphabets_of_translated_seqs(self):

        def triple_pad(s):
            """Add N to ensure length is a multiple of three (whole codons)."""
            while len(s) % 3:
                s += "N"
            return s

        self.assertEqual("IUPACProtein()",
                         repr(self.test_seqs[0].translate().alphabet))
        self.assertEqual("ExtendedIUPACProtein()",
                         repr(self.test_seqs[1].translate().alphabet))
        self.assertEqual("ExtendedIUPACProtein()",
                         repr(self.test_seqs[2].translate().alphabet))
        self.assertEqual("ExtendedIUPACProtein()",
                         repr(self.test_seqs[3].translate().alphabet))
        self.assertEqual("ExtendedIUPACProtein()",
                         repr(self.test_seqs[10].translate().alphabet))
        self.assertEqual("ExtendedIUPACProtein()",
                         repr(self.test_seqs[11].translate().alphabet))
        self.assertEqual("IUPACProtein()",
                         repr(self.test_seqs[12].translate().alphabet))
        self.assertEqual(
            "ExtendedIUPACProtein()",
            repr(triple_pad(self.test_seqs[13]).translate().alphabet))
        self.assertEqual("ExtendedIUPACProtein()",
                         repr(self.test_seqs[14].translate().alphabet))
        self.assertEqual("IUPACProtein()",
                         repr(self.test_seqs[15].translate().alphabet))
        self.assertEqual(
            "ExtendedIUPACProtein()",
            repr(triple_pad(self.test_seqs[16]).translate().alphabet))
        self.assertEqual(
            "ExtendedIUPACProtein()",
            repr(triple_pad(self.test_seqs[17]).translate().alphabet))

    def test_gapped_seq_with_gap_char_given(self):
        seq = Seq.Seq("ATG---AAACTG")
        self.assertEqual("M-KL", seq.translate(gap="-"))
        self.assertRaises(TranslationError, seq.translate, gap="~")

    def test_gapped_seq_with_stop_codon_and_gap_char_given(self):
        seq = Seq.Seq("GTG---GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
        self.assertEqual("V-AIVMGR*KGAR*", seq.translate(gap="-"))
        self.assertRaises(TranslationError, seq.translate)

    def test_gapped_seq_with_gap_char_given_and_inferred_from_alphabet(self):
        seq = Seq.Seq("ATG---AAACTG", Gapped(IUPAC.unambiguous_dna))
        self.assertEqual("M-KL", seq.translate(gap="-"))
        self.assertRaises(ValueError, seq.translate, gap="~")

        seq = Seq.Seq("ATG~~~AAACTG", Gapped(IUPAC.unambiguous_dna))
        self.assertRaises(ValueError, seq.translate, gap="~")
        self.assertRaises(TranslationError, seq.translate, gap="-")

    def test_gapped_seq_with_gap_char_given_and_inferred_from_alphabet2(self):
        """Test using stop codon in sequence"""
        seq = Seq.Seq("ATG---AAACTGTAG", Gapped(IUPAC.unambiguous_dna))
        self.assertEqual("M-KL*", seq.translate(gap="-"))
        self.assertRaises(ValueError, seq.translate, gap="~")

        seq = Seq.Seq("ATG---AAACTGTAG", Gapped(IUPAC.unambiguous_dna))
        self.assertEqual("M-KL@", seq.translate(gap="-", stop_symbol="@"))
        self.assertRaises(ValueError, seq.translate, gap="~")

        seq = Seq.Seq("ATG~~~AAACTGTAG", Gapped(IUPAC.unambiguous_dna))
        self.assertRaises(ValueError, seq.translate, gap="~")
        self.assertRaises(TranslationError, seq.translate, gap="-")

    def test_gapped_seq_no_gap_char_given(self):
        seq = Seq.Seq("ATG---AAACTG")
        self.assertRaises(TranslationError, seq.translate)

    def test_gapped_seq_no_gap_char_given_and_inferred_from_alphabet(self):
        seq = Seq.Seq("ATG---AAACTG", Gapped(IUPAC.unambiguous_dna))
        self.assertEqual("M-KL", seq.translate())

        seq = Seq.Seq("ATG~~~AAACTG", Gapped(IUPAC.unambiguous_dna))
        self.assertRaises(TranslationError, seq.translate)

        seq = Seq.Seq("ATG~~~AAACTG", Gapped(IUPAC.unambiguous_dna, "~"))
        self.assertEqual("M~KL", seq.translate())

    def test_alphabet_of_translated_gapped_seq(self):
        seq = Seq.Seq("ATG---AAACTG", Gapped(IUPAC.unambiguous_dna))
        self.assertEqual("Gapped(ExtendedIUPACProtein(), '-')",
                         repr(seq.translate().alphabet))

        seq = Seq.Seq("ATG---AAACTG", Gapped(IUPAC.unambiguous_dna, "-"))
        self.assertEqual("Gapped(ExtendedIUPACProtein(), '-')",
                         repr(seq.translate().alphabet))

        seq = Seq.Seq("ATG~~~AAACTG", Gapped(IUPAC.unambiguous_dna, "~"))
        self.assertEqual("Gapped(ExtendedIUPACProtein(), '~')",
                         repr(seq.translate().alphabet))

        seq = Seq.Seq("ATG---AAACTG")
        self.assertEqual("Gapped(ExtendedIUPACProtein(), '-')",
                         repr(seq.translate(gap="-").alphabet))

        seq = Seq.Seq("ATG~~~AAACTG")
        self.assertEqual("Gapped(ExtendedIUPACProtein(), '~')",
                         repr(seq.translate(gap="~").alphabet))

        seq = Seq.Seq("ATG~~~AAACTGTAG")
        self.assertEqual(
            "HasStopCodon(Gapped(ExtendedIUPACProtein(), '~'), '*')",
            repr(seq.translate(gap="~").alphabet))

        seq = Seq.Seq("ATG---AAACTGTGA")
        self.assertEqual(
            "HasStopCodon(Gapped(ExtendedIUPACProtein(), '-'), '*')",
            repr(seq.translate(gap="-").alphabet))

        seq = Seq.Seq("ATG---AAACTGTGA")
        self.assertEqual(
            "HasStopCodon(Gapped(ExtendedIUPACProtein(), '-'), '@')",
            repr(seq.translate(gap="-", stop_symbol="@").alphabet))

    def test_translation_wrong_type(self):
        """Test translation table cannot be CodonTable"""
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
            nucleotide_seq = nucleotide_seq[:3 * (len(nucleotide_seq) // 3)]
            if isinstance(nucleotide_seq, Seq.Seq) and \
               'X' not in str(nucleotide_seq):
                short = Seq.translate(nucleotide_seq, to_stop=True)
                self.assertEqual(
                    str(short),
                    str(Seq.translate(nucleotide_seq).split('*')[0]))

        seq = "GTGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
        self.assertEqual("VAIVMGRWKGAR", Seq.translate(seq, table=2,
                                                       to_stop=True))

    def test_translation_on_proteins(self):
        """Test translation shouldn't work on a protein!"""
        for s in protein_seqs:
            with self.assertRaises(ValueError):
                Seq.translate(s)

            if isinstance(s, Seq.Seq):
                with self.assertRaises(ValueError):
                    s.translate()

    def test_translation_of_invalid_codon(self):
        for codon in ["TA?", "N-N", "AC_", "Ac_"]:
            with self.assertRaises(TranslationError):
                Seq.translate(codon)

    def test_translation_of_glutamine(self):
        for codon in ['SAR', 'SAG', 'SAA']:
            self.assertEqual('Z', Seq.translate(codon))

    def test_translation_of_asparagine(self):
        for codon in ['RAY', 'RAT', 'RAC']:
            self.assertEqual('B', Seq.translate(codon))

    def test_translation_of_leucine(self):
        for codon in ['WTA', 'MTY', 'MTT', 'MTW', 'MTM', 'MTH', 'MTA', 'MTC',
                      'HTA']:
            self.assertEqual('J', Seq.translate(codon))

    def test_translation_with_bad_table_argument(self):
        table = dict()
        with self.assertRaises(ValueError):
            Seq.translate("GTGGCCATTGTAATGGGCCGC", table=table)

    def test_translation_with_codon_table_as_table_argument(self):
        table = standard_dna_table
        self.assertEqual("VAIVMGR", Seq.translate("GTGGCCATTGTAATGGGCCGC",
                                                  table=table))

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


class TestStopCodons(unittest.TestCase):
    def setUp(self):
        self.misc_stops = "TAATAGTGAAGAAGG"

    def test_stops(self):
        for nucleotide_seq in [self.misc_stops, Seq.Seq(self.misc_stops),
                               Seq.Seq(self.misc_stops,
                                       Alphabet.generic_nucleotide),
                               Seq.Seq(self.misc_stops,
                                       Alphabet.DNAAlphabet()),
                               Seq.Seq(self.misc_stops,
                                       IUPAC.unambiguous_dna)]:
            self.assertEqual("***RR", str(Seq.translate(nucleotide_seq)))
            self.assertEqual("***RR", str(Seq.translate(nucleotide_seq,
                                                        table=1)))
            self.assertEqual("***RR", str(Seq.translate(nucleotide_seq,
                                                        table="SGC0")))
            self.assertEqual("**W**", str(Seq.translate(nucleotide_seq,
                                                        table=2)))
            self.assertEqual("**WRR", str(Seq.translate(nucleotide_seq,
                                          table='Yeast Mitochondrial')))
            self.assertEqual("**WSS", str(Seq.translate(nucleotide_seq,
                                                        table=5)))
            self.assertEqual("**WSS", str(Seq.translate(nucleotide_seq,
                                                        table=9)))
            self.assertEqual("**CRR", str(Seq.translate(nucleotide_seq,
                                          table='Euplotid Nuclear')))
            self.assertEqual("***RR", str(Seq.translate(nucleotide_seq,
                                                        table=11)))
            self.assertEqual("***RR", str(Seq.translate(nucleotide_seq,
                                                        table='Bacterial')))

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
