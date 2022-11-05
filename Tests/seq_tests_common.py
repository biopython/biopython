# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Common code for SeqRecord object tests."""

import unittest

from Bio.Seq import UnknownSeq, UndefinedSequenceError
from Bio.SeqUtils.CheckSum import seguid
from Bio.SeqFeature import ExactPosition, UnknownPosition
from Bio.SeqFeature import SimpleLocation, CompoundLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from test_SeqIO import SeqIOTestBaseClass


class SeqRecordTestBaseClass(unittest.TestCase):
    def compare_reference(self, r1, r2):
        """Compare two Reference objects.

        Note r2 is assumed to be a BioSQL DBSeqRecord, due to limitations
        of the BioSQL table structure.
        """
        self.assertEqual(r1.title, r2.title)
        self.assertEqual(r1.authors, r2.authors)
        self.assertEqual(r1.journal, r2.journal)
        self.assertEqual(r1.medline_id, r2.medline_id)
        if r1.pubmed_id and r2.pubmed_id:
            self.assertEqual(r1.pubmed_id, r2.pubmed_id)
            # Looking at BioSQL/BioSeq.py function _retrieve_reference
            # it seems that it will get either the MEDLINE or PUBMED,
            # but not both.  I *think* the current schema does not allow
            # us to store both... must confirm this.

        # TODO - assert r1.comment == r2.comment
        # Looking at the tables, I *think* the current schema does not
        # allow us to store a reference comment.  Must confirm this.
        if r2.comment:
            self.assertEqual(r1.comment, r2.comment)

        # TODO - assert r1.consrtm == r2.consrtm
        # Looking at the tables, I *think* the current schema does not
        # allow us to store a consortium.
        if r2.consrtm:
            self.assertEqual(r1.consrtm, r2.consrtm)

        if len(r1.location) == 0:
            self.assertEqual(len(r2.location), 0)
        else:
            # BioSQL can only store ONE location!
            # TODO - Check BioPerl with a GenBank file with multiple ref locations
            self.assertIsInstance(r1.location[0], SimpleLocation)
            self.assertIsInstance(r2.location[0], SimpleLocation)
            self.assertEqual(r1.location[0].start, r2.location[0].start)
            self.assertEqual(r1.location[0].end, r2.location[0].end)

    def compare_feature(self, old_f, new_f):
        """Compare two SeqFeature objects."""
        self.assertIsInstance(old_f, SeqFeature)
        self.assertIsInstance(new_f, SeqFeature)
        self.assertEqual(old_f.type, new_f.type)
        self.assertEqual(old_f.strand, new_f.strand)
        self.assertEqual(old_f.ref, new_f.ref)
        self.assertEqual(old_f.ref_db, new_f.ref_db)
        # TODO - BioSQL does not store/retrieve feature's id (Bug 2526)
        if new_f.id != "<unknown id>":
            self.assertEqual(old_f.id, new_f.id)

        # TODO - Work out how the location_qualifier_value table should
        # be used, given BioPerl seems to ignore it (Bug 2766)
        # assert old_f.location_operator == new_f.location_operator, \
        #        "%s -> %s" % (old_f.location_operator, new_f.location_operator)

        # We dont store fuzzy locations:
        if not (
            isinstance(old_f.location.start, UnknownPosition)
            and isinstance(new_f.location.start, UnknownPosition)
        ):
            self.assertEqual(old_f.location.start, new_f.location.start)
        if not (
            isinstance(old_f.location.end, UnknownPosition)
            and isinstance(new_f.location.end, UnknownPosition)
        ):
            self.assertEqual(old_f.location.end, new_f.location.end)

        if isinstance(old_f.location, CompoundLocation):
            self.assertIsInstance(new_f.location, CompoundLocation)
        else:
            self.assertNotIsInstance(new_f.location, CompoundLocation)
        if isinstance(old_f.location, CompoundLocation):
            self.assertEqual(len(old_f.location.parts), len(new_f.location.parts))
            for old_l, new_l in zip(old_f.location.parts, new_f.location.parts):
                self.assertEqual(old_l.start, new_l.start)
                self.assertEqual(old_l.end, new_l.end)
                self.assertEqual(old_l.strand, new_l.strand)
                self.assertEqual(old_l.ref, new_l.ref)
                self.assertEqual(old_l.ref_db, new_l.ref_db)

        self.assertEqual(len(old_f.location.parts), len(new_f.location.parts))
        for old_sub, new_sub in zip(old_f.location.parts, new_f.location.parts):
            # These are SimpleLocation objects
            # Note UnknownPosition != UnknownPosition (just like NaN != NaN)
            if isinstance(old_sub.start, UnknownPosition):
                self.assertIsInstance(new_sub.start, UnknownPosition)
            else:
                self.assertEqual(old_sub.start, new_sub.start)
            if isinstance(old_sub.end, UnknownPosition):
                self.assertIsInstance(new_sub.end, UnknownPosition)
            else:
                self.assertEqual(old_sub.end, new_sub.end)
            self.assertEqual(old_sub.strand, new_sub.strand)

        self.assertCountEqual(old_f.qualifiers, new_f.qualifiers)
        for key in old_f.qualifiers:
            if isinstance(old_f.qualifiers[key], str):
                if isinstance(new_f.qualifiers[key], str):
                    self.assertEqual(old_f.qualifiers[key], new_f.qualifiers[key])
                elif isinstance(new_f.qualifiers[key], list):
                    # Maybe a string turning into a list of strings?
                    self.assertEqual([old_f.qualifiers[key]], new_f.qualifiers[key])
                else:
                    self.fail(f"Problem with feature's '{key}' qualifier")
            else:
                # Should both be lists of strings...
                self.assertEqual(old_f.qualifiers[key], new_f.qualifiers[key])

    def compare_features(self, old_list, new_list):
        """Compare two lists of SeqFeature objects."""
        self.assertIsInstance(old_list, list)
        self.assertIsInstance(new_list, list)
        self.assertEqual(len(old_list), len(new_list))
        for old_f, new_f in zip(old_list, new_list):
            self.compare_feature(old_f, new_f)

    def compare_sequence(self, old, new):
        """Compare two Seq objects."""
        self.assertEqual(len(old), len(new))

        if isinstance(old, UnknownSeq):
            self.assertIsInstance(new, UnknownSeq)
        else:
            self.assertNotIsInstance(new, UnknownSeq)

        self.assertEqual(len(old), len(new))
        try:
            bytes(old)
        except UndefinedSequenceError:
            self.assertRaises(UndefinedSequenceError, bytes, new)
            return

        self.assertEqual(old, new)

        ln = len(old)
        s = str(old)

        # Don't check every single element; for long sequences
        # this takes far far far too long to run!
        # Test both positive and negative indices
        if ln < 50:
            indices = list(range(-ln, ln))
        else:
            # A selection of end cases, and the mid point
            indices = [-ln, -ln + 1, -(ln // 2), -1, 0, 1, ln // 2, ln - 2, ln - 1]

        # Test element access,
        for i in indices:
            expected = s[i]
            self.assertEqual(expected, old[i])
            self.assertEqual(expected, new[i])

        # Test slices
        indices.append(ln)  # check copes with overflows
        indices.append(ln + 1000)  # check copes with overflows
        for i in indices:
            for j in indices:
                expected = s[i:j]
                self.assertEqual(expected, old[i:j])
                self.assertEqual(expected, new[i:j])
                # Slicing with step of 1 should make no difference.
                # Slicing with step 3 might be useful for codons.
                for step in [1, 3]:
                    expected = s[i:j:step]
                    self.assertEqual(expected, old[i:j:step])
                    self.assertEqual(expected, new[i:j:step])

            # Check automatic end points
            expected = s[i:]
            self.assertEqual(expected, old[i:])
            self.assertEqual(expected, new[i:])

            expected = s[:i]
            self.assertEqual(expected, old[:i])
            self.assertEqual(expected, new[:i])

        # Check "copy" splice
        self.assertEqual(s, old[:])
        self.assertEqual(s, new[:])

    def compare_record(self, old, new):
        """Compare two SeqRecord or DBSeqRecord objects."""
        self.assertIsInstance(old, SeqRecord)
        self.assertIsInstance(new, SeqRecord)
        # Sequence:
        self.compare_sequence(old.seq, new.seq)
        # Basics:
        self.assertEqual(old.id, new.id)
        self.assertEqual(old.name, new.name)
        self.assertEqual(old.description, new.description)
        self.assertEqual(old.dbxrefs, new.dbxrefs)
        # Features:
        self.compare_features(old.features, new.features)

        # Annotation:
        # We are expecting to see some "extra" annotations appearing,
        # such as 'cross_references', 'dates', 'data_file_division',
        # 'ncbi_taxon' and 'gi'.
        # TODO - address these, see Bug 2681?
        new_keys = set(new.annotations).difference(old.annotations)
        new_keys = new_keys.difference(
            ["cross_references", "date", "data_file_division", "ncbi_taxid", "gi"]
        )
        self.assertEqual(
            len(new_keys),
            0,
            msg=f"Unexpected new annotation keys: {', '.join(new_keys)}",
        )
        missing_keys = set(old.annotations).difference(new.annotations)
        missing_keys = missing_keys.difference(
            ["gene_name", "ncbi_taxid", "structured_comment"]  # Can't store chimeras
        )
        self.assertEqual(
            len(missing_keys),
            0,
            msg=f"Unexpectedly missing annotation keys: {', '.join(missing_keys)}",
        )

        # In the short term, just compare any shared keys:
        for key in set(old.annotations).intersection(new.annotations):
            if key == "references":
                self.assertEqual(len(old.annotations[key]), len(new.annotations[key]))
                for old_r, new_r in zip(old.annotations[key], new.annotations[key]):
                    self.compare_reference(old_r, new_r)
            elif key == "comment":
                # Turn them both into containing strings for comparison - due to
                # line wrapping in GenBank etc we don't really expect the white
                # space to be 100% the same.
                if isinstance(old.annotations[key], list):
                    old_comment = " ".join(old.annotations[key])
                else:
                    old_comment = old.annotations[key]
                if isinstance(new.annotations[key], list):
                    new_comment = " ".join(new.annotations[key])
                else:
                    new_comment = new.annotations[key]
                old_comment = old_comment.replace("\n", " ").replace("  ", " ")
                new_comment = new_comment.replace("\n", " ").replace("  ", " ")
                self.assertEqual(
                    old_comment,
                    new_comment,
                    msg="Comment annotation changed by load/retrieve",
                )
            elif key in ["taxonomy", "organism", "source"]:
                # If there is a taxon id recorded, these fields get overwritten
                # by data from the taxon/taxon_name tables.  There is no
                # guarantee that they will be identical after a load/retrieve.
                self.assertTrue(
                    isinstance(new.annotations[key], str)
                    or isinstance(new.annotations[key], list)
                )
            elif isinstance(old.annotations[key], type(new.annotations[key])):
                self.assertEqual(
                    old.annotations[key],
                    new.annotations[key],
                    msg=f"Annotation '{key}' changed by load/retrieve",
                )
            elif isinstance(old.annotations[key], str) and isinstance(
                new.annotations[key], list
            ):
                # Any annotation which is a single string gets turned into
                # a list containing one string by BioSQL at the moment.
                self.assertEqual(
                    [old.annotations[key]],
                    new.annotations[key],
                    msg=f"Annotation '{key}' changed by load/retrieve",
                )
            elif isinstance(old.annotations[key], list) and isinstance(
                new.annotations[key], str
            ):
                self.assertEqual(
                    old.annotations[key],
                    [new.annotations[key]],
                    msg=f"Annotation '{key}' changed by load/retrieve",
                )

    def compare_records(self, old_list, new_list):
        """Compare two lists of SeqRecord objects."""
        self.assertEqual(len(old_list), len(new_list))
        for old_r, new_r in zip(old_list, new_list):
            self.compare_record(old_r, new_r)
