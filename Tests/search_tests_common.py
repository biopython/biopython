# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio.SearchIO import HSPFragment
from Bio.SeqRecord import SeqRecord


def _num_difference(obj_a, obj_b):
    """Returns the number of instance attributes presence only in one object."""
    attrs_a = obj_a.__dict__.keys()
    attrs_b = obj_b.__dict__.keys()
    diff = set(attrs_a).symmetric_difference(set(attrs_b))
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
    if not isinstance(obj_a, HSPFragment):
        # check the number of hits contained
        assert len(obj_a) == len(obj_b), "length: %r vs %r" % (len(obj_a), \
                len(obj_b), obj_a, obj_b)
        for item_a, item_b in zip(obj_a, obj_b):
            assert compare_search_obj(item_a, item_b)

    return True


def compare_attrs(obj_a, obj_b, attrs):
    """Compares attribute values of two objects."""
    for attr in attrs:
        # don't check for private attributes
        if attr.startswith('_'):
            # except for _query and _hit
            if attr not in ('_query', '_hit'):
                continue
        # get attribute values from each objects
        val_a = getattr(obj_a, attr)
        val_b = getattr(obj_b, attr)

        # special case for HSP and HSPFragment {hit,query}
        # since they are seqrecords, we compare the strings only
        # comparing using compare_record is too slow
        if attr in ('_hit', '_query') and (val_a is not None and val_b is \
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
            keys_a, values_a = val_a.keys(), val_a.values()
            keys_b, values_b = val_b.keys(), val_b.values()
            # sort all values and keys
            [x.sort() for x in (keys_a, values_a, keys_b, values_b)]
            assert keys_a == keys_b
            assert values_a == values_b
        else:
            assert val_a == val_b, "%s: %r vs %r" % (attr, val_a, val_b)

    return True
