# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio.SearchIO import _ITERATOR_MAP


SUPPORTED = _ITERATOR_MAP.keys() + ['mock'] # mock format for other tests


def compare_qresult(qres_a, qres_b, fmt=None):
    """Compares attribute values of two QueryResult objects."""
    # check that file format is supported
    assert fmt in SUPPORTED, "File format %r not supported" % fmt

    # check the number of hits contained
    assert len(qres_a) == len(qres_b), "length: %r vs %r" % (len(qres_a), \
            len(qres_b))

    # set common attributes
    attrs = ['id', 'program', 'target', 'version', 'hit_keys', ]
    # and add format-specific attributes
    if fmt == 'blast-xml':
        from Bio.SearchIO.BlastIO.blastxml import _ELEM_QRESULT_OPT, _ELEM_META
        attrs += [x[0] for x in _ELEM_QRESULT_OPT.values()] + \
                [x[0] for x in _ELEM_META.values()] + ['desc', 'seq_len']
    elif fmt == 'blast-tab':
        from Bio.SearchIO.BlastIO.blasttab import _COLUMN_QRESULT
        attrs += _COLUMN_QRESULT.values() + ['desc', 'rid', 'fields']
    elif fmt == 'hmmer-text':
        attrs += ['acc', 'desc', 'seq_len']

    # compare qresult attributes
    compare_attrs(qres_a, qres_b, attrs)

    # compare hit objects within
    for hit_a, hit_b in zip(qres_a, qres_b):
        assert compare_hit(hit_a, hit_b, fmt)

    return True


def compare_hit(hit_a, hit_b, fmt=None):
    """Compares attribute values of two Hit objects."""
    # check the number of hsps contained
    assert len(hit_a) == len(hit_b), "length: %r vs %r" % (len(hit_a), \
            len(hit_b))

    # set common attributes
    attrs = ['id', 'query_id', ]
    # and add format-specific attributes
    if fmt == 'blast-xml':
        from Bio.SearchIO.BlastIO.blastxml import _ELEM_HIT
        attrs += [x[0] for x in _ELEM_HIT.values()]
    elif fmt == 'blast-tab':
        from Bio.SearchIO.BlastIO.blasttab import _COLUMN_HIT
        attrs += _COLUMN_HIT.values()
    elif fmt == 'hmmer-text':
        attrs += ['evalue', 'bitscore', 'bias', 'domain_exp_num', \
                'domain_obs_num', 'desc', 'is_in_inclusion', ]

    # compare hit attributes
    compare_attrs(hit_a, hit_b, attrs)

    # compare hit objects within
    for hsp_a, hsp_b in zip(hit_a, hit_b):
        assert compare_hsp(hsp_a, hsp_b, fmt)

    return True


def compare_hsp(hsp_a, hsp_b, fmt=None):
    """Compares attribute values of two HSP objects."""

    # set common attributes
    # not comparing alignment attribute, since it's basically comparing
    # hit and query seqs
    attrs = ['hit_id', 'query_id', 'hit', 'query', 'hit_strand', \
            'hit_from', 'hit_to', 'hit_span', 'query_strand', 'query_from', \
            'query_to', 'query_span', 'alignment_annotation', ]

    if fmt == 'blast-xml':
        from Bio.SearchIO.BlastIO.blastxml import _ELEM_HSP
        attrs += [x[0] for x in _ELEM_HSP.values()]
    elif fmt == 'blast-tab':
        from Bio.SearchIO.BlastIO.blasttab import _COLUMN_HSP
        attrs += _COLUMN_HSP.values()
    elif fmt == 'hmmer-text':
        attrs += ['domain_index', 'is_in_inclusion', 'bitscore', 'bias', \
                'evalue', 'evalue_cond', 'hit_endtype', 'query_endtype', \
                'env_from', 'env_to', 'env_endtype', 'acc_avg']

    # compare hit attributes
    compare_attrs(hsp_a, hsp_b, attrs)

    return True


def compare_attrs(obj_a, obj_b, attrs):
    """Compares attribute values of two objects."""
    for attr in attrs:
        try:
            # get attribute values from each objects
            val_a = getattr(obj_a, attr)
            val_b = getattr(obj_b, attr)
        except AttributeError:
            # make sure if the attribute doesn't exist, it doesn't exist
            # in both objects
            assert not hasattr(obj_a, attr)
            assert not hasattr(obj_b, attr)
            # set proxy values for comparison below, which will
            # always return True
            val_a, val_b = None, None

        # special case for HSP.{hit,query}
        # since they are seqrecords, we compare the strings only
        # comparing using compare_record is too slow
        if attr in ['hit', 'query'] and (val_a is not None and val_b is \
                not None):
            assert val_a.seq.tostring() == val_b.seq.tostring(), \
                    "%s: %r vs %r" % (attr, val_a, val_b)
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
