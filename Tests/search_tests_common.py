# Copyright 2012 by Wibowo Arindrarto. All rights reserved.
# Revisions Copyright 2012-2015 by Peter Cock. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Common code for SearchIO tests."""


import os
import gzip
import unittest

try:
    import sqlite3
except ImportError:
    sqlite3 = None

from Bio import SearchIO
from Bio.SeqRecord import SeqRecord


class SearchTestBaseClass(unittest.TestCase):
    def compare_attrs(self, obj_a, obj_b, attrs):
        """Compare attribute values of two objects."""
        for attr in attrs:
            # don't check for contained items, they are handled separately
            if attr.startswith("_items"):
                continue
            # get attribute values from each objects
            val_a = getattr(obj_a, attr)
            val_b = getattr(obj_b, attr)

            # special case for HSP and HSPFragment {hit,query}
            # since they are seqrecords, we compare the strings only
            # comparing using compare_record is too slow
            if attr in ("_hit", "_query") and (val_a is not None and val_b is not None):
                # compare seq directly if it's a contiguous hsp
                if isinstance(val_a, SeqRecord) and isinstance(val_b, SeqRecord):
                    msg = f"Comparing attribute {attr}"
                    self.assertEqual(str(val_a.seq), str(val_b.seq), msg=msg)
                elif isinstance(val_a, list) and isinstance(val_b, list):
                    for seq_a, seq_b in zip(val_a, val_b):
                        msg = f"Comparing attribute {attr}"
                        self.assertEqual(str(seq_a.seq), str(seq_b.seq), msg=msg)
            else:
                self.assertIsInstance(val_b, type(val_a))
                msg = f"Comparing attribute {attr}"
                self.assertEqual(val_a, val_b)

    def compare_search_obj(self, obj_a, obj_b):
        """Compare attribute values of two QueryResult objects."""
        # check that both qresults contain the same instance attributes
        self.assertEqual(_num_difference(obj_a, obj_b), 0)

        # compare qresult attributes
        # if the above assertion pass, doesn't matter if we use a or be here
        self.compare_attrs(obj_a, obj_b, list(obj_a.__dict__))

        # compare objects recursively if it's not an HSPFragment
        if not isinstance(obj_a, SearchIO.HSPFragment):
            # check the number of hits contained
            msg = f"comparing {obj_a!r} vs {obj_b!r}"
            self.assertEqual(len(obj_a), len(obj_b), msg=msg)

            for item_a, item_b in zip(obj_a, obj_b):
                self.compare_search_obj(item_a, item_b)


class CheckRaw(unittest.TestCase):
    """Base class for testing index's get_raw method."""

    fmt = None  # define this in subclasses!

    def check_raw(self, filename, id, raw, **kwargs):
        """Index filename using keyword arguments, check get_raw(id)==raw."""
        idx = SearchIO.index(filename, self.fmt, **kwargs)
        raw = raw.encode()
        # Anticipate cases where the raw string and/or file uses different
        # newline characters ~ we set everything to \n.
        new = idx.get_raw(id)
        self.assertIsInstance(new, bytes, f"Didn't get bytes from {self.fmt} get_raw")
        self.assertEqual(raw.replace(b"\r\n", b"\n"), new.replace(b"\r\n", b"\n"))
        idx.close()

        # Now again, but using SQLite backend
        if sqlite3:
            idx = SearchIO.index_db(":memory:", filename, self.fmt, **kwargs)
            new = idx.get_raw(id)
            self.assertIsInstance(
                new, bytes, f"Didn't get bytes from {self.fmt} get_raw"
            )
            self.assertEqual(raw.replace(b"\r\n", b"\n"), new.replace(b"\r\n", b"\n"))
            idx.close()

        if os.path.isfile(filename + ".bgz"):
            # Do the tests again with the BGZF compressed file
            print(f"[BONUS {filename}.bgz]")
            self.check_raw(filename + ".bgz", id, raw, **kwargs)


class CheckIndex(SearchTestBaseClass):
    """Base class for testing indexing."""

    def check_index(self, filename, format, **kwargs):
        if filename.endswith(".bgz"):
            with gzip.open(filename) as handle:
                parsed = list(SearchIO.parse(handle, format, **kwargs))
        else:
            parsed = list(SearchIO.parse(filename, format, **kwargs))
        # compare values by index
        indexed = SearchIO.index(filename, format, **kwargs)
        self.assertEqual(
            len(parsed),
            len(indexed),
            "Should be %i records in %s, index says %i"
            % (len(parsed), filename, len(indexed)),
        )
        # compare values by index_db, only if sqlite3 is present
        if sqlite3 is not None:
            db_indexed = SearchIO.index_db(":memory:", [filename], format, **kwargs)
            self.assertEqual(
                len(parsed),
                len(db_indexed),
                "Should be %i records in %s, index_db says %i"
                % (len(parsed), filename, len(db_indexed)),
            )

        for qres in parsed:
            idx_qres = indexed[qres.id]
            # parsed and indexed qresult are different objects!
            self.assertNotEqual(id(qres), id(idx_qres))
            # but they should have the same attribute values
            self.compare_search_obj(qres, idx_qres)
            # sqlite3 comparison, only if it's present
            if sqlite3 is not None:
                dbidx_qres = db_indexed[qres.id]
                self.assertNotEqual(id(qres), id(dbidx_qres))
                self.compare_search_obj(qres, dbidx_qres)

        indexed.close()
        if sqlite3 is not None:
            db_indexed.close()
            db_indexed._con.close()

        if os.path.isfile(filename + ".bgz"):
            # Do the tests again with the BGZF compressed file
            print(f"[BONUS {filename}.bgz]")
            self.check_index(filename + ".bgz", format, **kwargs)


def _num_difference(obj_a, obj_b):
    """Return the number of instance attributes present only in one object."""
    attrs_a = set(obj_a.__dict__)
    attrs_b = set(obj_b.__dict__)
    diff = attrs_a.symmetric_difference(attrs_b)
    privates = len([x for x in diff if x.startswith("_")])
    return len(diff) - privates
