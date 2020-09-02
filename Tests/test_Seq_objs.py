# Copyright 2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Unittests for the Seq objects."""

import warnings
import unittest

from Bio import BiopythonWarning
from Bio import SeqIO
from Bio.Data.IUPACData import ambiguous_dna_values, ambiguous_rna_values
from Bio.Seq import Seq, UnknownSeq, MutableSeq, translate
from Bio.Data.CodonTable import TranslationError, CodonTable


# This is just the standard table with less stop codons
# (replaced with coding for O as an artificial example)

# Turn black code style off
# fmt: off

special_table = CodonTable(forward_table={
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "O",
    "TGT": "C", "TGC": "C", "TGA": "O", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"},
    start_codons=["TAA", "TAG", "TGA"],
    stop_codons=["TAG"],
)

Chilodonella_uncinata_table = CodonTable(forward_table={
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y",             "TAG": "Q",  # noqa: E241
    "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"},
    start_codons=["ATG"],
    stop_codons=["TAA"],
)

# Turn black code style on
# fmt: on


class StringMethodTests(unittest.TestCase):
    _examples = [
        # These are length 9, a multiple of 3 for translation tests:
        Seq("ACGTGGGGT"),
        Seq("ACGUGGGGU"),
        Seq("GG"),
        Seq("A"),
        UnknownSeq(1),
        UnknownSeq(1, character="n"),
        UnknownSeq(1, character="N"),
        UnknownSeq(12, character="N"),
        UnknownSeq(12, character="X"),
        UnknownSeq(12),
    ]
    for seq in _examples[:]:
        if isinstance(seq, Seq):
            _examples.append(seq.tomutable())
    _start_end_values = [0, 1, 2, 1000, -1, -2, -999, None]

    def _test_method(self, method_name, pre_comp_function=None, start_end=False):
        """Check this method matches the plain string's method."""
        if pre_comp_function is None:
            # Define a no-op function:

            def pre_comp_function(x):
                return x

        self.assertIsInstance(method_name, str)
        for example1 in self._examples:
            if not hasattr(example1, method_name):
                # e.g. MutableSeq does not support find
                continue
            str1 = str(example1)

            for example2 in self._examples:
                if not hasattr(example2, method_name):
                    # e.g. MutableSeq does not support find
                    continue
                if (
                    method_name in ("index", "rindex")
                    and isinstance(example1, MutableSeq)
                    and len(example2) > 1
                ):
                    # MutableSeq index only supports single entries
                    continue
                str2 = str(example2)

                try:
                    i = pre_comp_function(getattr(example1, method_name)(str2))
                except ValueError:
                    i = ValueError
                try:
                    j = pre_comp_function(getattr(str1, method_name)(str2))
                except ValueError:
                    j = ValueError
                if i != j:
                    raise ValueError(
                        "%r.%s(%r) = %r, not %r" % (example1, method_name, str2, i, j)
                    )

                try:
                    i = pre_comp_function(getattr(example1, method_name)(example2))
                except ValueError:
                    i = ValueError
                try:
                    j = pre_comp_function(getattr(str1, method_name)(str2))
                except ValueError:
                    j = ValueError
                if i != j:
                    raise ValueError(
                        "%r.%s(%r) = %r, not %r"
                        % (example1, method_name, example2, i, j)
                    )

                if start_end:
                    if isinstance(example1, MutableSeq):
                        # Does not support start/end arguments
                        continue
                    for start in self._start_end_values:
                        try:
                            i = pre_comp_function(
                                getattr(example1, method_name)(str2, start)
                            )
                        except ValueError:
                            i = ValueError
                        try:
                            j = pre_comp_function(
                                getattr(str1, method_name)(str2, start)
                            )
                        except ValueError:
                            j = ValueError
                        if i != j:
                            raise ValueError(
                                "%r.%s(%r, %i) = %r, not %r"
                                % (example1, method_name, str2, start, i, j)
                            )

                        for end in self._start_end_values:
                            try:
                                i = pre_comp_function(
                                    getattr(example1, method_name)(str2, start, end)
                                )
                            except ValueError:
                                i = ValueError
                            try:
                                j = pre_comp_function(
                                    getattr(str1, method_name)(str2, start, end)
                                )
                            except ValueError:
                                j = ValueError
                            if i != j:
                                raise ValueError(
                                    "%r.%s(%r, %i, %i) = %r, not %r"
                                    % (example1, method_name, str2, start, end, i, j,)
                                )

    def test_str_count(self):
        """Check matches the python string count method."""
        self._test_method("count", start_end=True)

    def test_str_count_overlap_GG(self):
        """Check our count_overlap method using GG."""
        # Testing with self._examples
        expected = [
            3,
            3,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        ]
        expected *= 2  # MutableSeq() Tests

        assert len(self._examples) == len(expected)

        for seq, exp in zip(self._examples, expected):
            # Using search term GG as a string
            self.assertEqual(seq.count_overlap("GG"), exp)
            self.assertEqual(seq.count_overlap("G" * 5), 0)
            # Using search term GG as a Seq
            self.assertEqual(seq.count_overlap(Seq("GG")), exp)
            self.assertEqual(seq.count_overlap(Seq("G" * 5)), 0)

    def test_count_overlap_start_end_GG(self):
        """Check our count_overlap method using GG with variable ends and starts."""
        # Testing Seq() and MutableSeq() with variable start and end arguments
        start_end_exp = [
            (1, 7, 3),
            (3, None, 3),
            (3, 6, 2),
            (4, 6, 1),
            (4, -1, 2),
            (-5, None, 2),
            (-5, 7, 2),
            (7, -5, 0),
            (-100, None, 3),
            (None, 100, 3),
            (-100, 1000, 3),
        ]

        testing_seq = "GTAGGGGAG"

        for start, end, exp in start_end_exp:
            self.assertEqual(Seq(testing_seq).count_overlap("GG", start, end), exp)
            self.assertEqual(
                MutableSeq(testing_seq).count_overlap("GG", start, end), exp
            )

        # Testing Seq() and MutableSeq() with a more heterogeneous sequenece
        self.assertEqual(Seq("GGGTGGTAGGG").count_overlap("GG"), 5)
        self.assertEqual(MutableSeq("GGGTGGTAGGG").count_overlap("GG"), 5)
        self.assertEqual(Seq("GGGTGGTAGGG").count_overlap("GG", 2, 8), 1)
        self.assertEqual(MutableSeq("GGGTGGTAGGG").count_overlap("GG", 2, 8), 1)
        self.assertEqual(Seq("GGGTGGTAGGG").count_overlap("GG", -11, 6), 3)
        self.assertEqual(MutableSeq("GGGTGGTAGGG").count_overlap("GG", -11, 6), 3)
        self.assertEqual(Seq("GGGTGGTAGGG").count_overlap("GG", 7, 2), 0)
        self.assertEqual(MutableSeq("GGGTGGTAGGG").count_overlap("GG", 7, 2), 0)
        self.assertEqual(Seq("GGGTGGTAGGG").count_overlap("GG", -2, -10), 0)

        # Testing UnknownSeq() with variable start and end arguments
        char_start_end_exp = [
            ("N", 1, 7, 0),
            ("N", 1, 7, 0),
            ("N", -4, None, 0),
            ("N", -4, None, 0),
            ("X", 1, 7, 0),
        ]

        for char, start, end, exp in char_start_end_exp:
            self.assertEqual(
                UnknownSeq(12, character=char).count_overlap("GG", start, end), exp
            )
        self.assertEqual(UnknownSeq(12, character="X").count_overlap("GG", 1, 7), 0)

        # Testing UnknownSeq() with some more cases including unusual edge cases
        substr_start_end_exp = [
            ("G", 100, 105, 0),
            ("G", -1, 4, 0),
            ("G", 4, -1, 0),
            ("G", -8, -2, 0),
            ("G", -2, -8, 0),
            ("G", 8, 2, 0),
            ("G", 2, 8, 0),
            ("GG", 8, 2, 0),
            ("GG", 2, 8, 0),
            ("GG", -5, -1, 0),
            ("GG", 1, 5, 0),
            ("GGG", None, None, 0),
            ("GGGGGGGGG", None, None, 0),
            ("GGG", 1, 2, 0),
        ]

        for substr, start, end, exp in substr_start_end_exp:
            self.assertEqual(
                UnknownSeq(7, character="N").count_overlap(substr, start, end), exp
            )
        self.assertEqual(UnknownSeq(7, character="N").count_overlap("GG", 1), 0)

    def test_str_count_overlap_NN(self):
        """Check our count_overlap method using NN."""
        # Testing with self._examples
        expected = [
            0,
            0,
            0,
            0,  # Seq() Tests
            0,
            0,
            0,
            11,
            0,
            0,
        ]  # UnknownSeq() Tests
        expected *= 2  # MutableSeq() Tests

        assert len(self._examples) == len(expected)

        for seq, exp in zip(self._examples, expected):
            # Using search term NN as a string
            self.assertEqual(seq.count_overlap("NN"), exp)
            self.assertEqual(seq.count_overlap("N" * 13), 0)
            # Using search term NN as a Seq
            self.assertEqual(seq.count_overlap(Seq("NN")), exp)
            self.assertEqual(seq.count_overlap(Seq("N" * 13)), 0)

    def test_count_overlap_start_end_NN(self):
        """Check our count_overlap method using NN with variable ends and starts."""
        # Testing Seq() and MutableSeq() with variable start and end arguments
        start_end_exp = [
            (1, 7, 0),
            (3, None, 0),
            (3, 6, 0),
            (4, 6, 0),
            (4, -1, 0),
            (-5, None, 0),
            (-5, 7, 0),
            (7, -5, 0),
            (-100, None, 0),
            (None, 100, 0),
            (-100, 1000, 0),
        ]

        testing_seq = "GTAGGGGAG"

        for start, end, exp in start_end_exp:
            self.assertEqual(Seq(testing_seq).count_overlap("NN", start, end), exp)
            self.assertEqual(
                MutableSeq(testing_seq).count_overlap("NN", start, end), exp
            )

        # Testing Seq() and MutableSeq() with a more heterogeneous sequenece
        self.assertEqual(Seq("GGGTGGTAGGG").count_overlap("NN"), 0)
        self.assertEqual(MutableSeq("GGGTGGTAGGG").count_overlap("NN"), 0)
        self.assertEqual(Seq("GGGTGGTAGGG").count_overlap("NN", 2, 8), 0)
        self.assertEqual(MutableSeq("GGGTGGTAGGG").count_overlap("NN", 2, 8), 0)
        self.assertEqual(Seq("GGGTGGTAGGG").count_overlap("NN", -11, 6), 0)
        self.assertEqual(MutableSeq("GGGTGGTAGGG").count_overlap("NN", -11, 6), 0)
        self.assertEqual(Seq("GGGTGGTAGGG").count_overlap("NN", 7, 2), 0)
        self.assertEqual(MutableSeq("GGGTGGTAGGG").count_overlap("NN", 7, 2), 0)
        self.assertEqual(Seq("GGGTGGTAGGG").count_overlap("NN", -10, -2), 0)

        # Testing UnknownSeq() with variable start and end arguments
        char_start_end_exp = [
            ("N", 1, 7, 5),
            ("N", 1, 7, 5),
            ("N", -4, None, 3),
            ("N", -4, None, 3),
            ("X", 1, 7, 0),
        ]

        for char, start, end, exp in char_start_end_exp:
            self.assertEqual(
                UnknownSeq(12, character=char).count_overlap("NN", start, end), exp
            )
        self.assertEqual(UnknownSeq(12, character="X").count_overlap("NN", 1, 7), 0)

        # Testing UnknownSeq() with some more cases including unusual edge cases
        substr_start_end_exp = [
            ("N", 100, 105, 0),
            ("N", -1, 4, 0),
            ("N", 4, -1, 2),
            ("N", -8, -2, 5),
            ("N", -2, -8, 0),
            ("N", 8, 2, 0),
            ("N", 2, 8, 5),
            ("NN", 8, 2, 0),
            ("NN", 2, 8, 4),
            ("NN", -5, -1, 3),
            ("NN", 1, 5, 3),
            ("NNN", None, None, 5),
            ("NNNNNNNNN", None, None, 0),
            ("NNN", 1, 2, 0),
        ]

        for substr, start, end, exp in substr_start_end_exp:
            self.assertEqual(
                UnknownSeq(7, character="N").count_overlap(substr, start, end), exp
            )
        self.assertEqual(UnknownSeq(7, character="N").count_overlap("NN", 1), 5)

    def test_str_find(self):
        """Check matches the python string find method."""
        self._test_method("find", start_end=True)

    def test_str_rfind(self):
        """Check matches the python string rfind method."""
        self._test_method("rfind", start_end=True)

    def test_str_index(self):
        """Check matches the python string index method."""
        self._test_method("index", start_end=True)

    def test_str_rindex(self):
        """Check matches the python string rindex method."""
        self._test_method("rindex", start_end=True)

    def test_str_startswith(self):
        """Check matches the python string startswith method."""
        self._test_method("startswith", start_end=True)
        self.assertTrue("ABCDE".startswith(("ABE", "OBE", "ABC")))

        # Now check with a tuple of sub sequences
        for example1 in self._examples:
            if not hasattr(example1, "startswith"):
                # e.g. MutableSeq does not support this
                continue
            subs = tuple(
                example1[start : start + 2] for start in range(0, len(example1) - 2, 3)
            )
            subs_str = tuple(str(s) for s in subs)

            self.assertEqual(
                str(example1).startswith(subs_str), example1.startswith(subs)
            )
            self.assertEqual(
                str(example1).startswith(subs_str), example1.startswith(subs_str)
            )  # strings!
            self.assertEqual(
                str(example1).startswith(subs_str, 3), example1.startswith(subs, 3)
            )
            self.assertEqual(
                str(example1).startswith(subs_str, 2, 6),
                example1.startswith(subs, 2, 6),
            )

    def test_str_endswith(self):
        """Check matches the python string endswith method."""
        self._test_method("endswith", start_end=True)
        self.assertTrue("ABCDE".endswith(("ABE", "OBE", "CDE")))

        # Now check with a tuple of sub sequences
        for example1 in self._examples:
            if not hasattr(example1, "endswith"):
                # e.g. MutableSeq does not support this
                continue
            subs = tuple(
                example1[start : start + 2] for start in range(0, len(example1) - 2, 3)
            )
            subs_str = tuple(str(s) for s in subs)

            self.assertEqual(str(example1).endswith(subs_str), example1.endswith(subs))
            self.assertEqual(
                str(example1).startswith(subs_str), example1.startswith(subs_str)
            )  # strings!
            self.assertEqual(
                str(example1).endswith(subs_str, 3), example1.endswith(subs, 3)
            )
            self.assertEqual(
                str(example1).endswith(subs_str, 2, 6), example1.endswith(subs, 2, 6)
            )

    def test_str_strip(self):
        """Check matches the python string strip method."""
        self._test_method("strip", pre_comp_function=str)

    def test_str_rstrip(self):
        """Check matches the python string rstrip method."""
        self._test_method("rstrip", pre_comp_function=str)

    def test_str_split(self):
        """Check matches the python string rstrip method."""
        # Calling split should return a list of Seq-like objects, we'll
        # just apply str() to each of them so it matches the string method
        self._test_method(
            "split", pre_comp_function=lambda x: [str(y) for y in x]  # noqa: E731
        )

    def test_str_rsplit(self):
        """Check matches the python string rstrip method."""
        # Calling rsplit should return a list of Seq-like objects, we'll
        # just apply str() to each of them so it matches the string method
        self._test_method(
            "rsplit", pre_comp_function=lambda x: [str(y) for y in x]  # noqa: E731
        )

    def test_str_lsplit(self):
        """Check matches the python string rstrip method."""
        # Calling lsplit should return a list of Seq-like objects, we'll
        # just apply str() to each of them so it matches the string method
        self._test_method(
            "lsplit", pre_comp_function=lambda x: [str(y) for y in x]  # noqa: E731
        )

    def test_str_length(self):
        """Check matches the python string __len__ method."""
        for example1 in self._examples:
            str1 = str(example1)
            self.assertEqual(len(example1), len(str1))

    def test_str_upper(self):
        """Check matches the python string upper method."""
        for example1 in self._examples:
            if isinstance(example1, MutableSeq):
                continue
            str1 = str(example1)
            self.assertEqual(str(example1.upper()), str1.upper())

    def test_str_lower(self):
        """Check matches the python string lower method."""
        for example1 in self._examples:
            if isinstance(example1, MutableSeq):
                continue
            str1 = str(example1)
            self.assertEqual(str(example1.lower()), str1.lower())

    def test_str_encode(self):
        """Check matches the python string encode method."""
        for example1 in self._examples:
            if isinstance(example1, MutableSeq):
                continue
            str1 = str(example1)
            self.assertEqual(example1.encode("ascii"), str1.encode("ascii"))
            self.assertEqual(example1.encode(), str1.encode())

    def test_str_hash(self):
        for example1 in self._examples:
            if isinstance(example1, MutableSeq):
                continue
            with warnings.catch_warnings():
                # Silence change in behaviour warning
                warnings.simplefilter("ignore", BiopythonWarning)
                self.assertEqual(
                    hash(str(example1)),
                    hash(example1),
                    "Hash mismatch, %r for %r vs %r for %r"
                    % (hash(str(example1)), id(example1), hash(example1), example1),
                )

    def test_str_comparison(self):
        for example1 in self._examples:
            for example2 in self._examples:
                with warnings.catch_warnings():
                    self.assertEqual(
                        str(example1) == str(example2),
                        example1 == example2,
                        "Checking %r == %r" % (example1, example2),
                    )
                    self.assertEqual(
                        str(example1) != str(example2),
                        example1 != example2,
                        "Checking %r != %r" % (example1, example2),
                    )
                    self.assertEqual(
                        str(example1) < str(example2),
                        example1 < example2,
                        "Checking %r < %r" % (example1, example2),
                    )
                    self.assertEqual(
                        str(example1) <= str(example2),
                        example1 <= example2,
                        "Checking %r <= %r" % (example1, example2),
                    )
                    self.assertEqual(
                        str(example1) > str(example2),
                        example1 > example2,
                        "Checking %r > %r" % (example1, example2),
                    )
                    self.assertEqual(
                        str(example1) >= str(example2),
                        example1 >= example2,
                        "Checking %r >= %r" % (example1, example2),
                    )

    def test_str_getitem(self):
        """Check slicing and indexing works like a string."""
        for example1 in self._examples:
            str1 = str(example1)
            for i in self._start_end_values:
                if i is not None and abs(i) < len(example1):
                    self.assertEqual(str(example1[i]), str1[i])
                self.assertEqual(str(example1[:i]), str1[:i])
                self.assertEqual(str(example1[i:]), str1[i:])
                for j in self._start_end_values:
                    self.assertEqual(str(example1[i:j]), str1[i:j])
                    for step in range(-3, 4):
                        if step == 0:
                            try:
                                print(example1[i:j:step])
                                self._assert(False)  # Should fail!
                            except ValueError:
                                pass
                        else:
                            self.assertEqual(str(example1[i:j:step]), str1[i:j:step])

    def test_tomutable(self):
        """Check obj.tomutable() method."""
        for example1 in self._examples:
            if isinstance(example1, MutableSeq):
                continue
            mut = example1.tomutable()
            self.assertIsInstance(mut, MutableSeq)
            self.assertEqual(str(mut), str(example1))

    def test_toseq(self):
        """Check obj.toseq() method."""
        for example1 in self._examples:
            try:
                seq = example1.toseq()
            except AttributeError:
                self.assertIsInstance(example1, Seq)
                continue
            self.assertIsInstance(seq, Seq)
            self.assertEqual(str(seq), str(example1))

    def test_the_complement(self):
        """Check obj.complement() method."""
        mapping = ""
        for example1 in self._examples:
            if isinstance(example1, MutableSeq):
                continue
            try:
                comp = example1.complement()
            except ValueError as e:
                self.assertEqual(str(e), "Proteins do not have complements!")
                continue
            str1 = str(example1)
            if "U" in str1 or "u" in str1:
                mapping = str.maketrans("ACGUacgu", "UGCAugca")
            else:
                # Default to DNA, e.g. complement("A") -> "T" not "U"
                mapping = str.maketrans("ACGTacgt", "TGCAtgca")
            self.assertEqual(str1.translate(mapping), str(comp))

    def test_the_reverse_complement(self):
        """Check obj.reverse_complement() method."""
        mapping = ""
        for example1 in self._examples:
            if isinstance(example1, MutableSeq):
                continue
            try:
                comp = example1.reverse_complement()
            except ValueError as e:
                self.assertEqual(str(e), "Proteins do not have complements!")
                continue
            str1 = str(example1)
            if "U" in str1 or "u" in str1:
                mapping = str.maketrans("ACGUacgu", "UGCAugca")
            else:
                # Defaults to DNA, so reverse_complement("A") --> "T" not "U"
                mapping = str.maketrans("ACGTacgt", "TGCAtgca")
            self.assertEqual(str1.translate(mapping)[::-1], str(comp))

    def test_the_transcription(self):
        """Check obj.transcribe() method."""
        mapping = ""
        for example1 in self._examples:
            if isinstance(example1, MutableSeq):
                continue
            try:
                tran = example1.transcribe()
            except ValueError as e:
                if str(e) == "Proteins cannot be transcribed!":
                    continue
                if str(e) == "RNA cannot be transcribed!":
                    continue
                raise
            str1 = str(example1)
            if len(str1) % 3 != 0:
                # TODO - Check for or silence the expected warning?
                continue
            self.assertEqual(str1.replace("T", "U").replace("t", "u"), str(tran))

    def test_the_back_transcription(self):
        """Check obj.back_transcribe() method."""
        mapping = ""
        for example1 in self._examples:
            if isinstance(example1, MutableSeq):
                continue
            try:
                tran = example1.back_transcribe()
            except ValueError as e:
                if str(e) == "Proteins cannot be back transcribed!":
                    continue
                if str(e) == "DNA cannot be back transcribed!":
                    continue
                raise
            str1 = str(example1)
            self.assertEqual(str1.replace("U", "T").replace("u", "t"), str(tran))

    def test_the_translate(self):
        """Check obj.translate() method."""
        mapping = ""
        for example1 in self._examples:
            if isinstance(example1, MutableSeq):
                continue
            if len(example1) % 3 != 0:
                # TODO - Check for or silence the expected warning?
                continue
            try:
                tran = example1.translate()
            except ValueError as e:
                if str(e) == "Proteins cannot be translated!":
                    continue
                raise
            # Try with positional vs named argument:
            self.assertEqual(example1.translate(11), example1.translate(table=11))

            # TODO - check the actual translation, and all the optional args

    def test_the_translation_of_stops(self):
        """Check obj.translate() method with stop codons."""
        misc_stops = "TAATAGTGAAGAAGG"
        nuc = Seq(misc_stops)
        self.assertEqual("***RR", str(nuc.translate()))
        self.assertEqual("***RR", str(nuc.translate(1)))
        self.assertEqual("***RR", str(nuc.translate("SGC0")))
        self.assertEqual("**W**", str(nuc.translate(table=2)))
        self.assertEqual("**WRR", str(nuc.translate(table="Yeast Mitochondrial")))
        self.assertEqual("**WSS", str(nuc.translate(table=5)))
        self.assertEqual("**WSS", str(nuc.translate(table=9)))
        self.assertEqual("**CRR", str(nuc.translate(table="Euplotid Nuclear")))
        self.assertEqual("***RR", str(nuc.translate(table=11)))
        self.assertEqual("***RR", str(nuc.translate(table="11")))
        self.assertEqual("***RR", str(nuc.translate(table="Bacterial")))
        self.assertEqual("**GRR", str(nuc.translate(table=25)))
        self.assertEqual("", str(nuc.translate(to_stop=True)))
        self.assertEqual("O*ORR", str(nuc.translate(table=special_table)))
        self.assertEqual("*QWRR", str(nuc.translate(table=Chilodonella_uncinata_table)))
        # These test the Bio.Seq.translate() function - move these?:
        self.assertEqual(
            "*QWRR", translate(str(nuc), table=Chilodonella_uncinata_table)
        )
        self.assertEqual("O*ORR", translate(str(nuc), table=special_table))
        self.assertEqual("", translate(str(nuc), to_stop=True))
        self.assertEqual("***RR", translate(str(nuc), table="Bacterial"))
        self.assertEqual("***RR", translate(str(nuc), table="11"))
        self.assertEqual("***RR", translate(str(nuc), table=11))
        self.assertEqual("**W**", translate(str(nuc), table=2))
        self.assertEqual(str(Seq("TAT").translate()), "Y")
        self.assertEqual(str(Seq("TAR").translate()), "*")
        self.assertEqual(str(Seq("TAN").translate()), "X")
        self.assertEqual(str(Seq("NNN").translate()), "X")
        self.assertEqual(str(Seq("TAt").translate()), "Y")
        self.assertEqual(str(Seq("TaR").translate()), "*")
        self.assertEqual(str(Seq("TaN").translate()), "X")
        self.assertEqual(str(Seq("nnN").translate()), "X")
        self.assertEqual(str(Seq("tat").translate()), "Y")
        self.assertEqual(str(Seq("tar").translate()), "*")
        self.assertEqual(str(Seq("tan").translate()), "X")
        self.assertEqual(str(Seq("nnn").translate()), "X")

    def test_the_translation_of_invalid_codons(self):
        """Check obj.translate() method with invalid codons."""
        for codon in ["TA?", "N-N", "AC_", "Ac_"]:
            nuc = Seq(codon)
            try:
                nuc.translate()
                self.fail("Translating %s should fail" % codon)
            except TranslationError:
                pass

    def test_the_translation_of_ambig_codons(self):
        """Check obj.translate() method with ambiguous codons."""
        for ambig_values in [ambiguous_dna_values, ambiguous_rna_values]:
            ambig = set(ambig_values.keys())
            ambig.remove("X")
            for c1 in ambig:
                for c2 in ambig:
                    for c3 in ambig:
                        values = {
                            str(Seq(a + b + c).translate())
                            for a in ambig_values[c1]
                            for b in ambig_values[c2]
                            for c in ambig_values[c3]
                        }
                        t = str(Seq(c1 + c2 + c3).translate())
                        if t == "*":
                            self.assertEqual(values, set("*"))
                        elif t == "X":
                            self.assertGreater(
                                len(values),
                                1,
                                "translate('%s') = '%s' not '%s'"
                                % (c1 + c2 + c3, t, ",".join(values)),
                            )
                        elif t == "Z":
                            self.assertEqual(values, set("EQ"))
                        elif t == "B":
                            self.assertEqual(values, set("DN"))
                        elif t == "J":
                            self.assertEqual(values, set("LI"))
                        else:
                            self.assertEqual(values, set(t))
                        # TODO - Use the Bio.Data.IUPACData module for the
                        # ambiguous protein mappings?

    def test_init_typeerror(self):
        """Check Seq __init__ gives TypeError exceptions."""
        # Only expect it to take strings - not Seq objects!
        self.assertRaises(TypeError, Seq, 1066)
        self.assertRaises(TypeError, Seq, Seq("ACGT"))

    def test_MutableSeq_init_typeerror(self):
        """Check MutableSeq __init__ gives TypeError exceptions."""
        self.assertRaises(TypeError, MutableSeq, (Seq("A")))
        self.assertRaises(TypeError, MutableSeq, (UnknownSeq(1)))
        self.assertRaises(TypeError, MutableSeq, 1)
        self.assertRaises(TypeError, MutableSeq, 1.0)

    def test_join_Seq_TypeError(self):
        """Checks that a TypeError is thrown for all non-iterable types."""
        # No iterable types which contain non-accepted types either.

        spacer = Seq("NNNNN")
        self.assertRaises(TypeError, spacer.join, 5)
        self.assertRaises(TypeError, spacer.join, ["ATG", "ATG", 5, "ATG"])

    def test_join_UnknownSeq_TypeError_iter(self):
        """Checks that a TypeError is thrown for all non-iterable types."""
        # No iterable types which contain non-accepted types either.

        spacer = UnknownSeq(5, character="-")
        self.assertRaises(TypeError, spacer.join, 5)
        self.assertRaises(TypeError, spacer.join, ["ATG", "ATG", 5, "ATG"])

    def test_join_MutableSeq_TypeError_iter(self):
        """Checks that a TypeError is thrown for all non-iterable types."""
        # No iterable types which contain non-accepted types either.

        spacer = MutableSeq("MMMMM")
        self.assertRaises(TypeError, spacer.join, 5)
        self.assertRaises(TypeError, spacer.join, ["ATG", "ATG", 5, "ATG"])

    def test_join_Seq(self):
        """Checks if Seq join correctly concatenates sequence with the spacer."""
        spacer = Seq("NNNNN")
        self.assertEqual(
            "N" * 15, spacer.join([Seq("NNNNN"), Seq("NNNNN")]),
        )

        spacer1 = Seq("")
        spacers = [spacer1, Seq("NNNNN"), Seq("GGG")]
        example_strings = ["ATG", "ATG", "ATG", "ATG"]
        example_strings_seqs = ["ATG", "ATG", Seq("ATG"), "ATG"]

        # strings with empty spacer
        str_concatenated = spacer1.join(example_strings)

        self.assertEqual(str(str_concatenated), "".join(example_strings))

        for spacer in spacers:
            seq_concatenated = spacer.join(example_strings_seqs)
            self.assertEqual(str(seq_concatenated), str(spacer).join(example_strings))
            # Now try single sequence arguments, should join the letters
            for target in example_strings + example_strings_seqs:
                self.assertEqual(
                    str(spacer).join(str(target)), str(spacer.join(target))
                )

    def test_join_UnknownSeq(self):
        """Checks if UnknownSeq join correctly concatenates sequence with the spacer."""
        spacer1 = UnknownSeq(5, character="-")
        spacer2 = UnknownSeq(0, character="-")
        spacers = [spacer1, spacer2]

        self.assertEqual(
            "-" * 15,
            spacer1.join([UnknownSeq(5, character="-"), UnknownSeq(5, character="-")]),
        )
        self.assertEqual(
            "N" * 5 + "-" * 10,
            spacer1.join([Seq("NNNNN"), UnknownSeq(5, character="-")]),
        )

        example_strings = ["ATG", "ATG", "ATG", "ATG"]
        example_strings_seqs = ["ATG", "ATG", Seq("ATG"), "ATG"]

        # strings with empty spacer
        str_concatenated = spacer2.join(example_strings)

        self.assertEqual(str(str_concatenated), "".join(example_strings))

        for spacer in spacers:
            seq_concatenated = spacer.join(example_strings_seqs)
            self.assertEqual(str(seq_concatenated), str(spacer).join(example_strings))
            # Now try single sequence arguments, should join the letters
            for target in example_strings + example_strings_seqs:
                self.assertEqual(
                    str(spacer).join(str(target)), str(spacer.join(target))
                )

    def test_join_MutableSeq_mixed(self):
        """Check MutableSeq objects can be joined."""
        spacer = MutableSeq("NNNNN")
        self.assertEqual(
            "N" * 15, spacer.join([MutableSeq("NNNNN"), MutableSeq("NNNNN")]),
        )
        self.assertRaises(
            TypeError, spacer.join([Seq("NNNNN"), MutableSeq("NNNNN")]),
        )

    def test_join_Seq_with_file(self):
        """Checks if Seq join correctly concatenates sequence from a file with the spacer."""
        filename = "Fasta/f003"
        seqlist = [record.seq for record in SeqIO.parse(filename, "fasta")]
        seqlist_as_strings = [str(_) for _ in seqlist]

        spacer = Seq("NNNNN")
        spacer1 = Seq("")
        # seq objects with spacer
        seq_concatenated = spacer.join(seqlist)
        # seq objects with empty spacer
        seq_concatenated1 = spacer1.join(seqlist)

        ref_data = ref_data1 = ""
        ref_data = str(spacer).join(seqlist_as_strings)
        ref_data1 = str(spacer1).join(seqlist_as_strings)

        self.assertEqual(str(seq_concatenated), ref_data)
        self.assertEqual(str(seq_concatenated1), ref_data1)
        with self.assertRaises(TypeError):
            spacer.join(SeqIO.parse(filename, "fasta"))

    def test_join_UnknownSeq_with_file(self):
        """Checks if UnknownSeq join correctly concatenates sequence from a file with the spacer."""
        filename = "Fasta/f003"
        seqlist = [record.seq for record in SeqIO.parse(filename, "fasta")]
        seqlist_as_strings = [str(_) for _ in seqlist]

        spacer = UnknownSeq(0, character="-")
        spacer1 = UnknownSeq(5, character="-")
        # seq objects with spacer
        seq_concatenated = spacer.join(seqlist)
        # seq objects with empty spacer
        seq_concatenated1 = spacer1.join(seqlist)

        ref_data = ref_data1 = ""
        ref_data = str(spacer).join(seqlist_as_strings)
        ref_data1 = str(spacer1).join(seqlist_as_strings)

        self.assertEqual(str(seq_concatenated), ref_data)
        self.assertEqual(str(seq_concatenated1), ref_data1)
        with self.assertRaises(TypeError):
            spacer.join(SeqIO.parse(filename, "fasta"))

    def test_join_MutableSeq(self):
        """Checks if MutableSeq join correctly concatenates sequence with the spacer."""
        # Only expect it to take Seq objects and/or strings in an iterable!

        spacer1 = MutableSeq("")
        spacers = [
            spacer1,
            MutableSeq("NNNNN"),
            MutableSeq("GGG"),
        ]
        example_strings = ["ATG", "ATG", "ATG", "ATG"]
        example_strings_seqs = ["ATG", "ATG", Seq("ATG"), "ATG"]

        # strings with empty spacer
        str_concatenated = spacer1.join(example_strings)

        self.assertEqual(str(str_concatenated), "".join(example_strings))

        for spacer in spacers:
            seq_concatenated = spacer.join(example_strings_seqs)
            self.assertEqual(str(seq_concatenated), str(spacer).join(example_strings))

    def test_join_MutableSeq_with_file(self):
        """Checks if MutableSeq join correctly concatenates sequence from a file with the spacer."""
        filename = "Fasta/f003"
        seqlist = [record.seq for record in SeqIO.parse(filename, "fasta")]
        seqlist_as_strings = [str(_) for _ in seqlist]

        spacer = MutableSeq("NNNNN")
        spacer1 = MutableSeq("")
        # seq objects with spacer
        seq_concatenated = spacer.join(seqlist)
        # seq objects with empty spacer
        seq_concatenated1 = spacer1.join(seqlist)

        ref_data = ref_data1 = ""
        ref_data = str(spacer).join(seqlist_as_strings)
        ref_data1 = str(spacer1).join(seqlist_as_strings)

        self.assertEqual(str(seq_concatenated), ref_data)
        self.assertEqual(str(seq_concatenated1), ref_data1)
        with self.assertRaises(TypeError):
            spacer.join(SeqIO.parse(filename, "fasta"))

    # TODO - Addition...


class FileBasedTests(unittest.TestCase):
    """Test Seq objects created from files by SeqIO."""

    def test_unknown_seq_ungap(self):
        """Test ungap() works properly on UnknownSeq instances."""
        rec = SeqIO.read("GenBank/NT_019265.gb", "genbank")
        self.assertIsInstance(rec.seq, UnknownSeq)

        ungapped_seq = rec.features[1].extract(rec.seq).ungap("-")
        self.assertIsInstance(ungapped_seq, UnknownSeq)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
