# Copyright 2012 by Wibowo Arindrarto. All rights reserved.
# Revisions Copyright 2012 by Peter Cock. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from __future__ import print_function

import os
import gzip
import unittest
try:
    import sqlite3
except ImportError:
    sqlite3 = None

from Bio import SearchIO
from Bio._py3k import _as_bytes
from Bio.SeqRecord import SeqRecord


class CheckRaw(unittest.TestCase):

    """Base class for testing index's get_raw method."""

    def check_raw(self, filename, id, raw, **kwargs):
        """Index filename using **kwargs, check get_raw(id)==raw."""
        idx = SearchIO.index(filename, self.fmt, **kwargs)
        raw = _as_bytes(raw)
        # Anticipate cases where the raw string and/or file uses different
        # newline characters ~ we set everything to \n.
        self.assertEqual(raw.replace(b'\r\n', b'\n'),
                idx.get_raw(id).replace(b'\r\n', b'\n'))
        idx.close()

        #Now again, but using SQLite backend
        if sqlite3:
            idx = SearchIO.index_db(":memory:", filename, self.fmt, **kwargs)
            self.assertEqual(raw.replace(b'\r\n', b'\n'),
                    idx.get_raw(id).replace(b'\r\n', b'\n'))
            idx.close()

        if os.path.isfile(filename + ".bgz"):
            #Do the tests again with the BGZF compressed file
            print("[BONUS %s.bgz]" % filename)
            self.check_raw(filename + ".bgz", id, raw, **kwargs)


class CheckIndex(unittest.TestCase):

    """Base class for testing indexing."""

    def check_index(self, filename, format, **kwargs):
        # check if Python3 installation has sqlite3
        try:
            import sqlite3
        except ImportError:
            sqlite3 = None

        if filename.endswith(".bgz"):
            handle = gzip.open(filename)
            parsed = list(SearchIO.parse(handle, format, **kwargs))
            handle.close()
        else:
            parsed = list(SearchIO.parse(filename, format, **kwargs))
        # compare values by index
        indexed = SearchIO.index(filename, format, **kwargs)
        self.assertEqual(len(parsed), len(indexed),
                         "Should be %i records in %s, index says %i" \
                         % (len(parsed), filename, len(indexed)))
        # compare values by index_db, only if sqlite3 is present
        if sqlite3 is not None:
            db_indexed = SearchIO.index_db(':memory:', [filename], format, **kwargs)
            self.assertEqual(len(parsed), len(db_indexed),
                             "Should be %i records in %s, index_db says %i" \
                             % (len(parsed), filename, len(db_indexed)))

        for qres in parsed:
            idx_qres = indexed[qres.id]
            # parsed and indexed qresult are different objects!
            self.assertNotEqual(id(qres), id(idx_qres))
            # but they should have the same attribute values
            self.assertTrue(compare_search_obj(qres, idx_qres))
            # sqlite3 comparison, only if it's present
            if sqlite3 is not None:
                dbidx_qres = db_indexed[qres.id]
                self.assertNotEqual(id(qres), id(dbidx_qres))
                self.assertTrue(compare_search_obj(qres, dbidx_qres))

        indexed.close()
        if sqlite3 is not None:
            db_indexed.close()
            db_indexed._con.close()

        if os.path.isfile(filename + ".bgz"):
            #Do the tests again with the BGZF compressed file
            print("[BONUS %s.bgz]" % filename)
            self.check_index(filename + ".bgz", format, **kwargs)

def _num_difference(obj_a, obj_b):
    """Returns the number of instance attributes presence only in one object."""
    attrs_a = set(obj_a.__dict__.keys())
    attrs_b = set(obj_b.__dict__.keys())
    diff = attrs_a.symmetric_difference(attrs_b)
    privates = len([x for x in diff if x.startswith('_')])
    return len(diff) - privates


def compare_search_obj(obj_a, obj_b):
    """Compares attribute values of two QueryResult objects."""

    # check that both qresults contain the same instance attributes
    assert _num_difference(obj_a, obj_b) == 0

    # compare qresult attributes
    # if the above assertion pass, doesn't matter if we use a or be here
    compare_attrs(obj_a, obj_b, obj_a.__dict__.keys())

    # compare objects recursively if it's not an HSPFragment
    if not isinstance(obj_a, SearchIO.HSPFragment):
        # check the number of hits contained
        assert len(obj_a) == len(obj_b), "length: %r vs %r" % (len(obj_a),
                len(obj_b), obj_a, obj_b)
        for item_a, item_b in zip(obj_a, obj_b):
            assert compare_search_obj(item_a, item_b)

    return True


def compare_attrs(obj_a, obj_b, attrs):
    """Compares attribute values of two objects."""
    for attr in attrs:
        # don't check for contained items, they are handled separately
        if attr.startswith('_items'):
            continue
        # get attribute values from each objects
        val_a = getattr(obj_a, attr)
        val_b = getattr(obj_b, attr)

        # special case for HSP and HSPFragment {hit,query}
        # since they are seqrecords, we compare the strings only
        # comparing using compare_record is too slow
        if attr in ('_hit', '_query') and (val_a is not None and val_b is
                not None):
            # compare seq directly if it's a contiguous hsp
            if isinstance(val_a, SeqRecord) and isinstance(val_b, SeqRecord):
                assert str(val_a.seq) == str(val_b.seq), \
                        "%s: %r vs %r" % (attr, val_a, val_b)
            elif isinstance(val_a, list) and isinstance(val_b, list):
                for seq_a, seq_b in zip(val_a, val_b):
                    assert str(seq_a.seq) == str(seq_b.seq), \
                            "%s: %r vs %r" % (attr, seq_a, seq_b)
        # if it's a dictionary, compare values and keys
        elif isinstance(val_a, dict):
            assert isinstance(val_b, dict)
            keys_a = sorted(val_a.keys())
            values_a = sorted(val_a.values())
            keys_b = sorted(val_b.keys())
            values_b = sorted(val_b.values())
            assert keys_a == keys_b, "%s: %r vs %r" % (attr, keys_a, keys_b)
            assert values_a == values_b, "%s: %r vs %r" % (attr, values_a,
                    values_b)
        # if it's an alphabet, check the class names as alphabets are instances
        elif attr == '_alphabet':
            alph_a = val_a.__class__.__name__
            alph_b = val_b.__class__.__name__
            assert alph_a == alph_b, "%s: %r vs %r" % (attr, alph_a, alph_b)
        else:
            assert val_a == val_b, "%s: %r vs %r" % (attr, val_a, val_b)

    return True
