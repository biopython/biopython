# Copyright 2009-2017 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""SeqFeature related tests for Seq objects from Bio.SeqIO.

Initially this takes matched tests of GenBank and FASTA files from the NCBI
and confirms they are consistent using our different parsers.
"""
import unittest

from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from Bio.SeqFeature import WithinPosition, BeforePosition, AfterPosition, OneOfPosition


class SeqRecordCreation(unittest.TestCase):
    """Test basic creation of SeqRecords."""

    def test_annotations(self):
        """Pass in annotations to SeqRecords."""
        rec = Seq("ACGT", id="Test", name="Test", description="Test")
        self.assertEqual(rec.annotations, {})
        rec = Seq(
            "ACGT",
            id="Test",
            name="Test",
            description="Test",
            annotations={"test": ["a test"]},
        )
        self.assertEqual(rec.annotations["test"], ["a test"])

    def test_letter_annotations(self):
        """Pass in letter annotations to SeqRecords."""
        rec = Seq("ACGT", id="Test", name="Test", description="Test")
        self.assertEqual(rec.annotations, {})
        rec = Seq("ACGT",
            id="Test",
            name="Test",
            description="Test",
            letter_annotations={"test": (1, 2, 3, 4)},
        )
        self.assertEqual(rec.letter_annotations["test"], (1, 2, 3, 4))
        # Now try modifying it to a bad value...
        try:
            rec.letter_annotations["bad"] = "abc"
            self.fail("Adding a bad letter_annotation should fail!")
        except (TypeError, ValueError) as e:
            pass
        # Now try setting it afterwards to a bad value...
        rec = Seq("ACGT", id="Test", name="Test", description="Test")
        try:
            rec.letter_annotations = {"test": [1, 2, 3]}
            self.fail("Changing to bad letter_annotations should fail!")
        except (TypeError, ValueError) as e:
            pass
        # Now try setting it at creation time to a bad value...
        try:
            rec = Seq("ACGT",
                id="Test",
                name="Test",
                description="Test",
                letter_annotations={"test": [1, 2, 3]},
            )
            self.fail("Wrong length letter_annotations should fail!")
        except (TypeError, ValueError) as e:
            pass

    def test_valid_id(self):
        with self.assertRaises(TypeError):
            Seq("ACGT", id={})

    def test_valid_name(self):
        with self.assertRaises(TypeError):
            Seq("ACGT", name={})

    def test_valid_description(self):
        with self.assertRaises(TypeError):
            Seq("ACGT", description={})

    def test_valid_dbxrefs(self):
        with self.assertRaises(TypeError):
            Seq("ACGT", dbxrefs={})

    def test_valid_annotations(self):
        with self.assertRaises(TypeError):
            Seq("ACGT", annotations=[])

    def test_valid_features(self):
        with self.assertRaises(TypeError):
            Seq("ACGT", features={})


class SeqRecordMethods(unittest.TestCase):
    """Test SeqRecord methods."""

    def setUp(self):
        f0 = SeqFeature(
            FeatureLocation(0, 26),
            type="source",
            qualifiers={"mol_type": ["fake protein"]},
        )
        f1 = SeqFeature(FeatureLocation(0, ExactPosition(10)))
        f2 = SeqFeature(
            FeatureLocation(WithinPosition(12, left=12, right=15), BeforePosition(22))
        )
        f3 = SeqFeature(
            FeatureLocation(
                AfterPosition(16),
                OneOfPosition(26, [ExactPosition(25), AfterPosition(26)]),
            )
        )
        self.record = Seq("ABCDEFGHIJKLMNOPQRSTUVWZYX",
            id="TestID",
            name="TestName",
            description="TestDescr",
            dbxrefs=["TestXRef"],
            annotations={"k": "v"},
            letter_annotations={"fake": "X" * 26},
            features=[f0, f1, f2, f3],
        )

    def test_iter(self):
        for amino in self.record:
            self.assertEqual("A", amino)
            break

    def test_contains(self):
        self.assertIn(Seq("ABC"), self.record)

    def test_str(self):
        expected = """
ABCDEFGHIJKLMNOPQRSTUVWZYX
ID: TestID
Name: TestName
Description: TestDescr
Database cross-references: TestXRef
Number of features: 4
/k=v
Per letter annotation for: fake"""
        self.assertEqual(expected.lstrip(), format(self.record, "debug"))

    def test_repr(self):
        expected = (
            "Seq('ABCDEFGHIJKLMNOPQRSTUVWZYX', "
            "id='TestID', name='TestName', description='TestDescr')"
        )
        self.assertEqual(expected, repr(self.record))

    def test_format(self):
        expected = ">TestID TestDescr\nABCDEFGHIJKLMNOPQRSTUVWZYX\n"
        self.assertEqual(expected, format(self.record, "fasta"))

    def test_format_str(self):
        expected = ">TestID TestDescr\nABCDEFGHIJKLMNOPQRSTUVWZYX\n"
        self.assertEqual(expected, f"{self.record:fasta}")

    def test_format_str_binary(self):
        with self.assertRaisesRegex(
            ValueError, "Binary format sff cannot be used with Seq format method"
        ):
            f"{self.record:sff}"

    def test_format_spaces(self):
        rec = Seq("ABCDEFGHIJKLMNOPQRSTUVWZYX",
            id="TestID",
            name="TestName",
            description="TestDescr",
        )
        rec.description = "TestDescr     with5spaces"
        expected = ">TestID TestDescr     with5spaces\nABCDEFGHIJKLMNOPQRSTUVWZYX\n"
        self.assertEqual(expected, format(rec, "fasta"))

    def test_upper(self):
        self.assertEqual("ABCDEFGHIJKLMNOPQRSTUVWZYX", self.record.lower().upper())

    def test_lower(self):
        self.assertEqual("abcdefghijklmnopqrstuvwzyx", self.record.lower())

    def test_slicing(self):
        self.assertEqual("B", self.record[1])
        self.assertEqual("BC", self.record[1:3])
        with self.assertRaises(TypeError):
            c = self.record["a"].seq

    def test_slice_variants(self):
        """Simple slices using different start/end values."""
        for start in list(range(-30, 30)) + [None]:
            for end in list(range(-30, 30)) + [None]:
                if start is None and end is None:
                    continue
                seq = self.record[start:end]
                seq_str = str(self.record)[start:end]
                self.assertEqual(seq_str, str(seq))
                self.assertEqual("X" * len(seq_str), seq.letter_annotations["fake"])

    def test_slice_simple(self):
        """Simple slice."""
        rec = self.record
        self.assertEqual(len(rec), 26)
        left = rec[:10]
        self.assertEqual(str(left), str(rec)[:10])
        right = rec[-10:]
        self.assertEqual(str(right), str(rec)[-10:])
        mid = rec[12:22]
        self.assertEqual(str(mid), str(rec)[12:22])
        for sub in [left, right, mid]:
            self.assertEqual(len(sub), 10)
            self.assertEqual(sub.id, "TestID")
            self.assertEqual(sub.name, "TestName")
            self.assertEqual(sub.description, "TestDescr")
            self.assertEqual(sub.letter_annotations, {"fake": "X" * 10})
            self.assertEqual(sub.dbxrefs, [])  # May change this...
            self.assertEqual(sub.annotations, {})  # May change this...
            self.assertEqual(len(sub.features), 1)
            # By construction, each feature matches the full sliced region:
            self.assertEqual(str(sub.features[0].extract(sub)), str(sub))
            self.assertEqual(sub.features[0].extract(str(sub)), str(sub))

    def test_slice_zero(self):
        """Zero slice."""
        rec = self.record
        self.assertEqual(len(rec), 26)
        self.assertEqual(len(rec[2:-2]), 22)
        self.assertEqual(len(rec[5:2]), 0)
        self.assertEqual(len(rec[5:2][2:-2]), 0)

    def test_add_simple(self):
        """Simple addition."""
        rec = self.record + self.record
        self.assertEqual(len(rec), 52)
        self.assertEqual(rec.id, "TestID")
        self.assertEqual(rec.name, "TestName")
        self.assertEqual(rec.description, "TestDescr")
        self.assertEqual(rec.dbxrefs, ["TestXRef"])
        self.assertEqual(rec.annotations, {"k": "v"})
        self.assertEqual(rec.letter_annotations, {"fake": "X" * 52})
        self.assertEqual(len(rec.features), 2 * len(self.record.features))

    def test_add_seq(self):
        """Simple addition of string."""
        rec = self.record + "BIO"  # will use Seq's __add__ method
        self.assertEqual(len(rec), 26 + 3)
        self.assertEqual(str(rec), str(self.record) + "BIO")
        self.assertEqual(rec.id, "TestID")
        self.assertEqual(rec.name, "TestName")
        self.assertEqual(rec.description, "TestDescr")
        self.assertEqual(rec.dbxrefs, ["TestXRef"])
        self.assertEqual(rec.annotations, {"k": "v"})
        self.assertEqual(rec.letter_annotations, {})
        self.assertEqual(len(rec.features), len(self.record.features))
        self.assertEqual(rec.features[0].type, "source")
        self.assertEqual(rec.features[0].location.nofuzzy_start, 0)
        self.assertEqual(rec.features[0].location.nofuzzy_end, 26)  # not +3

    def test_add_seqrecord(self):
        """Simple left addition of SeqRecord from genbank file."""
        other = SeqIO.read("GenBank/dbsource_wrap.gb", "gb")
        other.dbxrefs = ["dummy"]
        rec = self.record + other
        self.assertEqual(len(rec), len(self.record) + len(other))
        self.assertEqual(str(rec), str(self.record) + str(other))
        self.assertEqual(rec.id, "")
        self.assertEqual(rec.name, "")
        self.assertEqual(rec.description, "")
        self.assertEqual(rec.dbxrefs, ["TestXRef", "dummy"])
        self.assertEqual(len(rec.annotations), 0)
        self.assertEqual(len(rec.letter_annotations), 0)
        self.assertEqual(
            len(rec.features), len(self.record.features) + len(other.features)
        )
        self.assertEqual(rec.features[0].type, "source")
        self.assertEqual(rec.features[0].location.nofuzzy_start, 0)
        self.assertEqual(
            rec.features[0].location.nofuzzy_end, len(self.record)
        )  # not +3
        i = len(self.record.features)
        self.assertEqual(rec.features[i].type, "source")
        self.assertEqual(rec.features[i].location.nofuzzy_start, len(self.record))
        self.assertEqual(rec.features[i].location.nofuzzy_end, len(rec))

    def test_add_seq_left(self):
        """Simple left addition of Seq or string."""
        rec = "BIO" + self.record  # will use Seq's __radd__ method
        self.assertEqual(len(rec), 26 + 3)
        self.assertEqual(str(rec), "BIO" + str(self.record))
        self.assertEqual(rec.id, "TestID")
        self.assertEqual(rec.name, "TestName")
        self.assertEqual(rec.description, "TestDescr")
        self.assertEqual(rec.dbxrefs, ["TestXRef"])
        self.assertEqual(rec.annotations, {"k": "v"})
        self.assertEqual(rec.letter_annotations, {})
        self.assertEqual(len(rec.features), len(self.record.features))
        self.assertEqual(rec.features[0].type, "source")
        self.assertEqual(rec.features[0].location.nofuzzy_start, 3)
        self.assertEqual(rec.features[0].location.nofuzzy_end, 26 + 3)

    def test_slice_add_simple(self):
        """Simple slice and add."""
        for cut in range(27):
            rec = self.record[:cut] + self.record[cut:]
            self.assertEqual(str(rec), str(self.record))
            self.assertEqual(len(rec), 26)
            self.assertEqual(rec.id, "TestID")
            self.assertEqual(rec.name, "TestName")
            self.assertEqual(rec.description, "TestDescr")
            self.assertEqual(rec.dbxrefs, [])  # May change this...
            self.assertEqual(rec.annotations, {})  # May change this...
            self.assertEqual(rec.letter_annotations, {"fake": "X" * 26})
            self.assertLessEqual(len(rec.features), len(self.record.features))

    def test_slice_add_shift(self):
        """Simple slice and add to shift."""
        for cut in range(27):
            rec = self.record[cut:] + self.record[:cut]
            self.assertEqual(
                str(rec), str(self.record[cut:] + self.record[:cut])
            )
            self.assertEqual(len(rec), 26)
            self.assertEqual(rec.id, "TestID")
            self.assertEqual(rec.name, "TestName")
            self.assertEqual(rec.description, "TestDescr")
            self.assertEqual(rec.dbxrefs, [])  # May change this...
            self.assertEqual(rec.annotations, {})  # May change this...
            self.assertEqual(rec.letter_annotations, {"fake": "X" * 26})
            self.assertLessEqual(len(rec.features), len(self.record.features))


class SeqRecordMethodsMore(unittest.TestCase):
    """Test SeqRecord methods cont."""

    # This class does not have a setUp defining self.record

    def test_reverse_complement_seq(self):
        s = Seq("ACTG",
            id="TestID",
            name="TestName",
            description="TestDescription",
            dbxrefs=["TestDbxrefs"],
            features=[SeqFeature(FeatureLocation(0, 3), type="Site")],
            annotations={"organism": "bombyx"},
            letter_annotations={"test": "abcd"},
        )
        rc = s.reverse_complement(
            id=True,
            name=True,
            description=True,
            dbxrefs=True,
            features=True,
            annotations=True,
            letter_annotations=True,
        )

        self.assertEqual("CAGT", rc)
        self.assertEqual("TestID", rc.id)
        self.assertEqual("TestID", s.reverse_complement(id="TestID").id)

        self.assertEqual("TestName", rc.name)
        self.assertEqual("TestName", s.reverse_complement(name="TestName").name)

        self.assertEqual("TestDescription", rc.description)
        self.assertEqual(
            "TestDescription",
            s.reverse_complement(description="TestDescription").description,
        )

        self.assertEqual(["TestDbxrefs"], rc.dbxrefs)
        self.assertEqual(
            ["TestDbxrefs"], s.reverse_complement(dbxrefs=["TestDbxrefs"]).dbxrefs
        )

        self.assertEqual(
            "[SeqFeature(FeatureLocation(ExactPosition(1), ExactPosition(4)), type='Site')]",
            repr(rc.features),
        )
        rc2 = s.reverse_complement(
            features=[SeqFeature(FeatureLocation(1, 4), type="Site")]
        )
        self.assertEqual(
            "[SeqFeature(FeatureLocation(ExactPosition(1), ExactPosition(4)), type='Site')]",
            repr(rc2.features),
        )

        self.assertEqual({"organism": "bombyx"}, rc.annotations)
        self.assertEqual(
            {"organism": "bombyx"},
            s.reverse_complement(annotations={"organism": "bombyx"}).annotations,
        )

        self.assertEqual({"test": "dcba"}, rc.letter_annotations)
        self.assertEqual(
            {"test": "abcd"},
            s.reverse_complement(
                letter_annotations={"test": "abcd"}
            ).letter_annotations,
        )

    def test_reverse_complement_mutable_seq(self):
        s = MutableSeq("ACTG")
        self.assertEqual("CAGT", str(s.reverse_complement(inplace=False)))

    def test_translate(self):
        s = Seq("ATGGTGTAA",
            id="TestID",
            name="TestName",
            description="TestDescription",
            dbxrefs=["TestDbxrefs"],
            features=[SeqFeature(FeatureLocation(0, 3), type="Site")],
            annotations={"organism": "bombyx"},
            letter_annotations={"test": "abcdefghi"},
        )

        t = s.translate()
        self.assertEqual(t, "MV*")
        self.assertEqual(t.id, "")
        self.assertEqual(t.name, "")
        self.assertEqual(t.description, "")
        self.assertFalse(t.dbxrefs)
        self.assertFalse(t.features)
        self.assertEqual(t.annotations, {"molecule_type": "protein"})
        self.assertFalse(t.letter_annotations)

        t = s.translate(
            cds=True,
            id=True,
            name=True,
            description=True,
            dbxrefs=True,
            annotations=True,
        )
        self.assertEqual(t, "MV")
        self.assertEqual(t.id, "TestID")
        self.assertEqual(t.name, "TestName")
        self.assertEqual(t.description, "TestDescription")
        self.assertEqual(t.dbxrefs, ["TestDbxrefs"])
        self.assertFalse(t.features)
        self.assertEqual(
            t.annotations, {"organism": "bombyx", "molecule_type": "protein"}
        )
        self.assertFalse(t.letter_annotations)


class TestTranslation(unittest.TestCase):
    def setUp(self):
        self.s = Seq("ATGGTGTAA",
            id="TestID",
            name="TestName",
            description="TestDescription",
            dbxrefs=["TestDbxrefs"],
            features=[SeqFeature(FeatureLocation(0, 3), type="Site")],
            annotations={"organism": "bombyx"},
            letter_annotations={"test": "abcdefghi"},
        )

    def test_defaults(self):
        t = self.s.translate()
        self.assertEqual(t, "MV*")
        self.assertEqual(t.id, "")
        self.assertEqual(t.name, "")
        self.assertEqual(t.description, "")
        self.assertFalse(t.dbxrefs)
        self.assertFalse(t.features)
        self.assertEqual(t.annotations, {"molecule_type": "protein"})
        self.assertFalse(t.letter_annotations)

    def test_preserve(self):
        t = self.s.translate(
            cds=True,
            id=True,
            name=True,
            description=True,
            dbxrefs=True,
            annotations=True,
        )
        self.assertEqual(t, "MV")
        self.assertEqual(t.id, "TestID")
        self.assertEqual(t.name, "TestName")
        self.assertEqual(t.description, "TestDescription")
        self.assertEqual(t.dbxrefs, ["TestDbxrefs"])
        self.assertFalse(t.features)
        self.assertEqual(
            t.annotations, {"organism": "bombyx", "molecule_type": "protein"}
        )
        self.assertFalse(t.letter_annotations)

        # Should not preserve these
        self.assertRaises(TypeError, self.s.translate, features=True)
        self.assertRaises(TypeError, self.s.translate, letter_annotations=True)

    def test_new_annot(self):
        t = self.s.translate(
            1,
            to_stop=True,
            gap="-",
            id="Foo",
            name="Bar",
            description="Baz",
            dbxrefs=["Nope"],
            features=[SeqFeature(FeatureLocation(0, 3), type="Site")],
            annotations={"a": "team"},
            letter_annotations={"aa": ("Met", "Val")},
        )
        self.assertEqual(t, "MV")
        self.assertEqual(t.id, "Foo")
        self.assertEqual(t.name, "Bar")
        self.assertEqual(t.description, "Baz")
        self.assertEqual(t.dbxrefs, ["Nope"])
        self.assertEqual(len(t.features), 1)
        self.assertEqual(t.annotations, {"a": "team", "molecule_type": "protein"})
        self.assertEqual(t.letter_annotations, {"aa": ("Met", "Val")})


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
